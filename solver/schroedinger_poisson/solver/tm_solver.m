%
% monacoQC: An object-oriented Matlab-based device engineering tool for
% quantum cascade structures.
%
% Copyright (C) 2025, Computational Photonics Group, Technical University of
% Munich
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

classdef tm_solver < handle
    % A solver class calculating eigenstates using the transfer matrix
    % method (TMM).

    properties
        sim_const % sim_constants-object: Simulation constants describing the TMM system.
        boundary_cond % boundary_conditions-object: Boundary conditions for the TMM system.
        error = 1e-5 % scalar: Error parameter for convergence of SP loop.
        max_iteration = 20 % scalar: Max iteration of SP loop.
        num_E = 2000 % scalar: Number of energy steps for eigenstates search.
    end

    methods
        function obj = tm_solver(d, s)
            % Constructs an object of type tm_solver.
            %
            % Syntax:
            %   obj = tm_solver(d, s)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.

            % Simulation constants.
            obj.sim_const = sim_constants(d, s);
            % Boundary conditions.
            obj.boundary_cond = boundary_condition.create_instance(d, ...
                s, obj.sim_const);
        end

        function [eig_system, cond_profile] = solve(obj, ...
                func_carr_distr, eig_system_old)
            % Solves the coupled system of Schrödinger and Poisson
            % equations for a given setup and boundary conditions.
            %
            % Syntax:
            %   [eig_system, cond_profile] = solve(obj)
            %   [eig_system, cond_profile] = solve(obj, func_carr_distr)
            %   [eig_system, cond_profile] = solve(obj, func_carr_distr, eig_system_old)
            %
            % Input Arguments:
            %   func_carr_distr (carrier_distribution-object): Contains
            %     subband occupations which are required as input for
            %     the Poisson solver. A uniform disitribution is assumed
            %     by default.
            %   eig_system_old (eigenstates-object): Contains eigenstates
            %     which are used to initialize the electrostatic
            %     potential originating from the mobile charge carriers.
            %
            % Output Arguments:
            %   eig_system (eigenstates-object): Contains the converged
            %     wavefunctions, eigenenergies and energy dependent
            %     effective masses.
            %   cond_profile: Contains the biased conduction band
            %     potential profile including the effects of ionized
            %     space charges and mobile charge carriers.

            % Do inputs of eigenstates and/or function handle for
            % carrier distribution, with an eigenstate instance as input
            % argument, exist?
            if ~exist("func_carr_distr")
                func_carr_distr = @(eig_st)carrier_distribution. ...
                    generate_uniform(obj.sim_const.num_wavefct);
            end
            if ~exist('eig_system_old', 'var')
                eig_system_old = [];
            end
            % Calculate vector with charge density.
            rho = poisson_solver.get_rho(obj.sim_const, ...
                func_carr_distr, eig_system_old);

            total_time = tic;
            fprintf('Solving Schrödinger-Poisson equation ...\n\n')

            % SP loop
            for iter = 1:obj.max_iteration
                % Solve Poisson equation for additional electrostatic
                % potential.
                dV = poisson_solver.calc_potential(obj.sim_const, rho);
                % Calculate total potential in J.
                V = obj.sim_const.vec_V_0 + dV;
                % Generate object of conduction band profile.
                cond_profile = conduction_band(obj.sim_const.vec_z_tm, V);

                % Calculate eigenenergies and wavefunctions.
                [psi, E, m_E_eff, E_bound_CBO] = obj.calc_wavefcts(V);

                % Check the number of found wavefunctions.
                if (length(E_bound_CBO) < obj.sim_const.num_wavefct)
                    error(['Could not find the', ...
                        ' appropriate number of wavefunctions!']);
                end

                % Returns object of system eigenstates.
                eig_system = obj.gen_eig_syst(psi, E, m_E_eff, ...
                    E_bound_CBO, cond_profile);

                rho_old = rho;
                % Calculate updated charge density from new eigenstates.
                rho = poisson_solver.get_rho(obj.sim_const, ...
                    func_carr_distr, eig_system);
                % Termination condition (relative change in sheet density).
                err = sum((rho - rho_old).^2) / sum(rho.^2);

                disp('('+string(iter)+'/'+string(obj.max_iteration)+ ...
                    ')'+' residual: '+string(err));
                if err < obj.error
                    fprintf(['\n', 'Schrödinger-Poisson solver reached', ...
                        ' desired residual of ', num2str(obj.error), ...
                        '. Elapsed time: ', num2str(toc(total_time)), ...
                        's\n\n'])
                    break
                end
                if iter == obj.max_iteration
                    fprintf(['\n', 'Schrödinger-Poisson solver did not ', ...
                        'reach desired residual of ', num2str(obj.error), ...
                        '. Elapsed time: ', num2str(toc(total_time)), ...
                        's\n\n'])
                end
            end
        end

        function [psi, E, m_E_eff, E_bound_CBO] = ...
                calc_wavefcts(obj, V)
            % Calculates eigenstates for bound states using transfer
            % matrix method and a root finding algorithm (fzero).
            %
            % Syntax:
            %   [psi, E, m_E_eff, E_bound_CBO] = calc_wavefcts(obj, V)
            %
            % Input Arguments:
            %   V (vector): Potential profile.
            %
            % Output Arguments:
            %   psi (matrix): Matrix containing the wavefunction of the
            %     states in one period as columns.
            %   E (vector): Eigenergies of the states in one period [J].
            %   m_E_eff (vector): In-plane effective mass of the states
            %     in one period [-].
            %   E_bound_CBO: Eigenenergies with respect to conduction
            %     band edge of the states in one period [J].

            % Energy interval (Emin, Emax)
            Emin = V(obj.boundary_cond.ind_E_min);
            Emax = Emin + obj.boundary_cond.Eperiod;

            % Energy spacing DE
            DE = (Emax - Emin) / obj.num_E;
            opt = optimset('TolX', obj.boundary_cond.dE);

            % Number of scans depends on boundary condition.
            % For instance, for QCDs it is necessary to divide the
            % simulation domain into two scans.
            for n_scans = 1:obj.boundary_cond.num_scans
                % Function to find wavefunctions.
                fun_tmm = @(x) obj.calc_tmm_result(V, x, n_scans);
                % Energy interval [Ea Eb] with fun_tmm evaluations a, b.
                b = 0;
                Eb = Emin;

                % Solutions for eigenfunctions and eigenenergies.
                psi_tm = [];
                E_tm = [];
                % Calculation of wavefunctions and energy eigenvalues
                % using fzero and TMM (method transfermat).
                for E_i = Emin:DE:Emax
                    a = b;
                    Ea = Eb;
                    % Start value for root finding algorithm.
                    [b] = fun_tmm(E_i);
                    % Start interval [Ea E_i].
                    Eb = E_i;
                    if (a * b < 0)
                        % Find root eigenenergy in interval [Ea Eb].
                        E_tm(end+1) = fzero(fun_tmm, [Ea, Eb], opt);

                        % Calculate corresponding wavefunction psi.
                        tm_coeff_AB = ...
                            obj.transfermat(V, E_tm(end), n_scans);
                        % Extend wavefunction over whole simulation domain.
                        psi_tm(:, end+1) = ...
                            obj.extend_wavefct(E_tm(end), V, ...
                            tm_coeff_AB, n_scans);
                    end
                end

                % Eigenenergies with respect to conduction band edge
                % in J.
                [E_bound_CBO_tm] = obj.calc_E_bound_CBO(E_tm, psi_tm, ...
                    V, n_scans);

                psi_scan{n_scans} = psi_tm(obj.sim_const.vec_ind_tm, :);
                % Eigenenergies in J.
                E_scan{n_scans} = E_tm;
                E_bound_CBO_scan{n_scans} = E_bound_CBO_tm;
            end
            % Concatenate scans in one vector.
            E_bound_CBO = [E_bound_CBO_scan{:}];
            psi = [psi_scan{:}];
            E = [E_scan{:}];
            % Sort eigenenergies ascending relative to the
            % conduction band edge. Eq. 6 in Jirauschek and Kubis 2014
            % (https://aip.scitation.org/doi/abs/10.1063/1.4863665)
            [E_bound_CBO, ind] = sort(E_bound_CBO);
            psi = psi(:, ind);
            E = E(ind);
            % Delete unbound wavefunctions.
            E(isnan(E_bound_CBO)) = [];
            psi(:, isnan(E_bound_CBO)) = [];
            E_bound_CBO(isnan(E_bound_CBO)) = [];
            % Calculate vector with effective masses.
            m_E_eff = obj.get_meff_E(V, E, psi);
        end

        function psi = extend_wavefct(obj, Eh, V, ...
                tm_coeff_AB, n_scans)
            % Constructs wavefunction from the transfer matrix
            % coefficients. The wavefunctions extend over all considered
            % QCL periods.
            %
            % Syntax:
            %   psi = extend_wavefct(obj, Eh, V, tm_coeff_AB, n_scans)
            %
            % Input Arguments:
            %   Eh (scalar): Eigenenergy of a specific state.
            %   V (vector): Potential profile [J].
            %   tm_coeff_AB (matrix): 2-by-N matrix containing the
            %     transfer matrix coefficiences of complete device.
            %   n_scans (scalar): Current scan number.
            %
            % Output Arguments:
            %   psi (vector): Wavefunction extending over whole device.

            meff_np = obj.sim_const.get_meff_np(Eh, V);
            vec_z = obj.sim_const.vec_z_tm;
            ind_z_L = obj.boundary_cond.ind_z_L(n_scans);
            ind_z_R = obj.boundary_cond.ind_z_R(n_scans);
            % Wavenumber left boundary.
            KL = sqrt(2*meff_np(ind_z_L).* ...
                (V(ind_z_L) - Eh)) / phys_const.hbar;
            % Wavenumber right boundary.
            KR = sqrt(2*meff_np(ind_z_R).* ...
                (V(ind_z_R) - Eh)) / phys_const.hbar;
            % Calculated wavefunctions within the simulation period
            % based on the TMM amplitudes A and B.
            psi(ind_z_L:ind_z_R) = real(tm_coeff_AB(1, ind_z_L:ind_z_R) ...
                +tm_coeff_AB(2, ind_z_L:ind_z_R));

            % Extend wavefunctions to left and right simulation domain
            % beyond the TMM bounds.
            psi(1:(ind_z_L - 1)) = ...
                exp(KL*(obj.sim_const.z_unit * ...
                vec_z ...
                (1:(ind_z_L - 1)) - ...
                obj.sim_const.z_unit * ...
                vec_z(ind_z_L)));
            psi((ind_z_R + 1): ...
                size(vec_z, 2)) = psi(ind_z_R) * ...
                exp(-KR*(obj.sim_const.z_unit * ...
                vec_z((ind_z_R + 1):end) - ...
                obj.sim_const.z_unit * vec_z(ind_z_R)));
            % Real valued wavefunctions.
            psi = real(psi);
            % Normalize wavefunctions,
            % integral of the probability density has to be one.
            psi = psi / sqrt(trapz(vec_z, ...
                psi.^2));
        end

        function meff_E = get_meff_E(obj, V, E, psi)
            % Calculating in-plane effective mass according to Eq. 9 in
            % Jirauschek and Kubis 2014,
            % https://doi.org/10.1063/1.4863665.
            %
            % Syntax:
            %   meff_E = get_meff_E(obj, V, E, psi)
            %
            % Input Arguments:
            %   V (vector): Potential profile [J].
            %   E (vector): Eigenergies of the states in one period [J].
            %   psi (matrix): Wavefunctions of the states in one period.
            %
            % Output Arguments:
            %   meff_E (vector): Energy dependent effective mass of each
            %     state in one period [-].

            me_eff = ...
                obj.sim_const.vec_meff(obj.sim_const.vec_ind_tm) ...
                / phys_const.me;
            % Vector including nonparabolicity parameters.
            vec_parab = ...
                obj.sim_const.vec_parab(obj.sim_const.vec_ind_tm);
            % Calculated averaged in-plane effective mass
            % for the respective subband.
            vec_meff_E = me_eff .* ...
                (1 + vec_parab .* (E' - ...
                V(obj.sim_const.vec_ind_tm))) .* psi'.^2;
            for i = 1:size(vec_meff_E, 1)
                meff_E(i) = trapz(obj.sim_const.vec_z, ...
                    vec_meff_E(i, :));
            end
        end

        function tm_coeff_AB = transfermat(obj, V, E, n_scans)
            % Calculates the coefficiences of the transfer matrix method
            % according to chapter II D.1 in Jirauschek and Kubis 2014,
            % https://doi.org/10.1063/1.4863665.
            %
            % Syntax:
            %   tm_coeff_AB = transfermat(obj, V, E, n_scans)
            %
            % Input Arguments:
            %   V (vector): Potential profile [J].
            %   E (vector): Eigenenergies [J].
            %   n_scans (scalar): Current scan index.
            %
            % Output Arguments:
            %   tm_coeff_AB (matrix): 2-by-N matrix containing the
            %     transfer matrix coefficiences of complete device.

            % Indices of simulation domain, defined in boundary condition.
            indices = obj.boundary_cond.get_indices(n_scans);
            % Increment for simulation domain.
            incr_z = obj.boundary_cond.incr_z;

            vec_z = obj.sim_const.z_unit ...
                * obj.sim_const.vec_z_tm;
            % Calculate effective mass with respect to energy,
            % accounting for nonparabolicity effects.
            meff_np = obj.sim_const.get_meff_np(E, V);

            % Matrix dimension n.
            n = length(V);

            tm_coeff_AB = zeros(2, n);
            % Boundary conditions at start position.
            tm_coeff_AB(:, indices(1)) = obj.boundary_cond.get_bc_start();
            % Wavevector k.
            k = sqrt(2*meff_np.*(E - V)) / phys_const.hbar;
            % Delta z vector.
            D = diff(vec_z(indices));
            D(end+1) = D(end);

            if (sum(imag(k) < 0))
                error(['Sum over the imaginary part of wavevector k', ...
                    ' is smaller than zero!']);
            end
            % Normalization factor kn for wavenumber.
            kn = sqrt(2*phys_const.me*phys_const.e0) / phys_const.hbar;

            % Renorm for stepwise constant matrix methods.
            k = k / kn;
            D = kn * D;

            % Transfer matrix 2: step centered
            % See Eq. 17 in Jirauschek and Kubis 2014.
            T2 = [1, 0; 0, 1];
            for i = 1:length(indices)
                m = indices(i);

                k1 = k(m);
                k2 = k(m + incr_z);
                b1 = k(m) / meff_np(m);
                b2 = k(m + incr_z) / ...
                    meff_np(m + incr_z);
                T0 = [0.5 * (1 + b1 / b2) * exp(1i*(k1 + k2)*D(i)/2), ...
                    0.5 * (1 - b1 / b2) * exp(1i*(k2 - k1)*D(i)/2); ...
                    0.5 * (1 - b1 / b2) * exp(1i*(k1 - k2)*D(i)/2), ...
                    0.5 * (1 + b1 / b2) * exp(-1i*(k1 + k2)*D(i)/2)];
                T2 = T0 * T2;
                tm_coeff_AB(:, m+incr_z) = T2 * tm_coeff_AB(:, indices(1));
            end
        end

        function tm_coeff = calc_tmm_result(obj, V, E, n_scans)
            % Helper function which returns the transfer matrix
            % coefficient for a specific state at the boundary.
            %
            % Syntax:
            %   tm_coeff = calc_tmm_result(obj, V, E, n_scans)
            %
            % Input Argument:
            %   V (vector): Potential profile [J].
            %   E (scalar): Eigenenergy of a specific state [J].
            %   n_scans (scalar): Scan index.
            %
            % Output Arguments:
            %   tm_coeff (scalar): Transfer matrix coefficient.

            tm_coeff_AB = transfermat(obj, V, E, n_scans);
            index_bc_end = ...
                obj.boundary_cond.get_index_bc_end(size(tm_coeff_AB), ...
                n_scans);
            tm_coeff = real(tm_coeff_AB(index_bc_end));
        end

        function E_bound = calc_E_bound_CBO(obj, E, psi_tm, V, n_scans)
            % Calculates eigenenergy with respect to the conduction
            % band potential to find the most strongly bound levels.
            %
            % Syntax:
            %   E_bound = calc_E_bound_CBO(obj, E, psi_tm, V, n_scans)
            %
            % Input Arguments:
            %   E (vector): Eigenenergy of the states on one period
            %     [J].
            %   psi_tm (matrix): Matrix containing the wavefunctions of
            %     the states in one period as columns.
            %   V (vector): Potential profile [J].
            %   n_scans (scalar): Scan index.
            %
            % Output Arguments:
            %   E_bound (vector): Eigenenergies of the states in one
            %     period with respect to the conduction band edge.

            E_bound = zeros(1, length(E));
            % Check, whether solution eigenstate is located within
            % the simulation interval [z_L z_R].
            z_L = ...
                obj.boundary_cond.get_z_L(obj.sim_const.vec_z_tm, n_scans);
            z_R = ...
                obj.boundary_cond.get_z_R(obj.sim_const.vec_z_tm, n_scans);
            for i = 1:length(E)
                E_bound(i) = E(i) - trapz(obj.sim_const.vec_z_tm, ...
                    psi_tm(:, i).^2.*V');
                % If 50 % of the probability density is outside of the
                % TMM simulation domain, E_bound is set to NaN to
                % characterize the wavefunction solution as improper.
                if (trapz(obj.sim_const.vec_z_tm, ...
                        psi_tm(:, i).^2.*z_L'.*z_R') < 0.5)
                    E_bound(i) = NaN;
                end
            end
        end

        function eig_system = gen_eig_syst(obj, psi, E, m_E_eff, ...
                E_bound_CBO, cond_profile)
            % Constructs an eigenstates object containing the solution
            % to the Schrödinger-Poisson equations.
            %
            % Syntax:
            %   eig_system = gen_eig_syst(obj, psi, E, m_E_eff, E_bound_CBO, cond_profile)
            %
            % Input Arguments:
            %   psi (matrix): Matrix containing the wavefunctions of
            %     the states in one period as columns.
            %   E (vector): Eigenenergies of the states in period [J].
            %   m_E_eff (vector): Effective masses of the states in one
            %     period [-].
            %   E_bound_CBO: Eigenenergies of the states in one period
            %     with respect to conduction band edge [J].
            %   cond_profile (vector): Conduction band profile [J].

            % Save matrix A containing information about the found
            % eigenenergies, to be written into test.dat.
            A_test = [E', E_bound_CBO'] / phys_const.e0;
            % Interpolate wavefunctions over 4 periods.
            [psi_eig, E_eig, E_bound_CBO_eig, m_E_eff_eig] = ...
                interp_wfs(obj, psi, E, E_bound_CBO, m_E_eff);
            % Select boundary condition to generate system hamiltonian.
            if (strcmp('tb', obj.boundary_cond.name))
                % Use tight-binding base transform.
                Vt = tb_base_transform.get_V_t ...
                    (obj.boundary_cond.ind_z_L, ...
                    obj.boundary_cond.ind_z_R, cond_profile.Vh);
                % Find tight-binding period.
                ind_period_tb = ...
                    tb_base_transform.find_index_tb_period( ...
                    obj.boundary_cond.ind_z_L, ...
                    obj.boundary_cond.ind_z_R, obj.sim_const.vec_z, ...
                    psi_eig, obj.sim_const.num_periods_wf);
                % Calculate hamiltonian with respect to
                % the tight-binding potential.
                hamiltonian_eig = tb_base_transform.calc_hamiltonian( ...
                    obj.sim_const.vec_z, psi_eig, E_eig, ind_period_tb, ...
                    Vt(obj.sim_const.vec_ind_tm), ...
                    cond_profile.Vh(obj.sim_const.vec_ind_tm), ...
                    obj.sim_const.num_periods_wf);
            elseif (strcmp('ez', obj.boundary_cond.name))
                % Use ez-base transform.
                % Define start object of period eigenstates.
                eig_start = eigenstates(E_eig, psi_eig, ...
                    obj.sim_const.vec_z, m_E_eff_eig, ...
                    E_bound_CBO/phys_const.e0);
                % Finds multiplets marker for ez-transformation.
                marker = ...
                    ez_base_transform.find_multiplets(eig_start, ...
                    obj.boundary_cond.e_multiplet, ...
                    obj.boundary_cond.d_multiplet);
                % Transform period eigenstates into ez-states.
                [hamiltonian_eig, psi_eig, m_E_eff_eig, ~] = ...
                    ez_base_transform.transform(eig_start, marker);
            else
                hamiltonian_eig = diag(E_eig);
            end

            % Generate eigenstate object for the whole system.
            eig_system = eigenstates(hamiltonian_eig, psi_eig, ...
                obj.sim_const.vec_z, m_E_eff_eig, ...
                E_bound_CBO(1:obj.sim_const.num_wavefct)/phys_const.e0);
            % Save matrix A in eigenstate object.
            eig_system.A_test = A_test;
        end

        function [psi_interp, E_interp, E_bound_CBO_eig, ...
                m_E_eff_interp] = interp_wfs(obj, psi, E, ...
                E_bound_CBO, m_E_eff)
            % Interpolates wavefunctions and copies the wavefunctions
            % from a single period to all other periods of the device.
            %
            % Syntax:
            %   [psi_interp, E_interp, E_bound_CBO_eig, m_E_eff_interp] = interp_wfs(obj, psi, E, E_bound_CBO, m_E_eff)
            %
            % Input Arguments:
            %   psi (matrix): Matrix containing the wavefunctions of
            %     the states in one period as columns.
            %   E (vector): Eigenenergies of the states in one period [J].
            %   E_bound_CBO: Eigenenergies of the states in one period
            %     with respect to conduction band edge [J].
            %   m_E_eff (vector): Effective masses of the states in one
            %     period [-].
            %
            % Output Arguments:
            %   psi_interp (matrix): Matrix containing the wavefunctions
            %     in all periods as columns.
            %   E_interp (vector): Eigenenergies of the states in all
            %     periods [J].
            %   E_bound_CBO_eig (vector): Eigenenergies with respect to
            %     conduction band edge of the states in all periods [J].
            %   m_E_eff_interp (vector): Effective masses of the states
            %     in all periods [-].

            % Sort eigenstates ascending with respect to the eigenenergies.
            E = E(1:obj.sim_const.num_wavefct);
            [E, ind] = sort(E);
            psi = psi(:, ind);
            m_E_eff = m_E_eff(ind);
            E_bound_CBO_eig = E_bound_CBO(ind);
            % Interpolate wavefunctions into the most likely adjacent
            % periods.
            % The wavefunctions will be filled to the left until the stop
            % criterion is reached. Here, we stop if the probability
            % density value is below the margin:
            % 1 - boundary_cond.delta_psi.
            % The missing number of wavefunctions will be interpolated
            % to periods on the right side.
            psi_interp_cell = cell(length(E), 1);
            E_interp_cell = cell(length(E), 1);
            zv_interp_R = ...
                obj.sim_const.vec_z_tm(obj.boundary_cond.ind_interp);
            for i = 1:length(E)
                psi_interp_cell{i} = psi(:, i);
                E_interp_cell{i} = E(i);
                num_period = 1;
                while (num_period < obj.sim_const.num_periods_wf)
                    while (trapz(obj.sim_const.vec_z, ...
                            obj.sim_const.vec_z'.* ...
                            psi_interp_cell{i}(:, 1).^2) < zv_interp_R)
                        psi_R = interp1(obj.sim_const.vec_z, ...
                            psi_interp_cell{i}(:, end), ...
                            obj.sim_const.vec_z-1* ...
                            obj.sim_const.l_period, 'v5cubic');
                        psi_R(isnan(psi_R)) = 0;
                        if (trapz(obj.sim_const.vec_z, psi_R.^2) >= 1 - ...
                                obj.boundary_cond.delta_psi)
                            psi_interp_cell{i} = [psi_interp_cell{i}, ...
                                psi_R'];
                            E_interp_cell{i} = [E_interp_cell{i}, ...
                                E_interp_cell{i}(end) + ...
                                obj.sim_const.E_period];
                            num_period = num_period + 1;
                        else
                            break
                        end
                    end
                    psi_L = interp1(obj.sim_const.vec_z, ...
                        psi_interp_cell{i}(:, 1), ...
                        obj.sim_const.vec_z+1*obj.sim_const.l_period, ...
                        'v5cubic');
                    psi_L(isnan(psi_L)) = 0;
                    psi_interp_cell{i} = [psi_L', psi_interp_cell{i}];
                    E_interp_cell{i} = [E_interp_cell{i}(1) ...
                        - obj.sim_const.E_period, E_interp_cell{i}];
                    num_period = num_period + 1;
                end
            end

            % Bring interpolated eigenstates into matrix form, which serves
            % then as input for the generation of a class object
            % eigenstates.
            E_interp = zeros(obj.sim_const.num_periods_wf* ...
                obj.sim_const.num_wavefct, 1);
            psi_interp = zeros(length(obj.sim_const.vec_z), ...
                obj.sim_const.num_periods_wf*obj.sim_const.num_wavefct);
            m_E_eff_interp = ...
                zeros(obj.sim_const.num_periods_wf* ...
                obj.sim_const.num_wavefct, 1);
            for ni = 1:obj.sim_const.num_wavefct
                for n_period = 0:obj.sim_const.num_periods_wf - 1
                    psi_interp(:, ...
                        ni+n_period*obj.sim_const.num_wavefct) = ...
                        psi_interp_cell{ni}(:, 1 + n_period);
                    E_interp(ni+n_period*obj.sim_const.num_wavefct) = ...
                        E_interp_cell{ni}(1 + n_period) / phys_const.e0;
                    m_E_eff_interp(ni+n_period* ...
                        obj.sim_const.num_wavefct) = ...
                        m_E_eff(ni);
                end
            end
        end
    end
end
