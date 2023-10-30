classdef tm_solver < solver
    %tm_solver A solver class calculating eigenstates
    % using the transfer matrix method.
    properties
        
        sim_const % Simulation constants describing the TMM system.
        boundary_cond % Boundary conditions for the TMM system.
        error = 1e-5; % error parameter for convergence of SP loop.
        num_E = 2000; % Number of energy steps.
        
    end
    
    methods
        % Constructs tm_solver.
        function obj = tm_solver(d, s)
            obj = obj@solver('transfer matrix method');
            % Simulation constants.
            obj.sim_const = sim_constants(d, s);
            % Boundary conditions.
            obj.boundary_cond = boundary_condition.create_instance(d, ...
                s, obj.sim_const);
        end
        % Solves the Schrödinger equation for the given setup and boundary
        % conditions, and returns an eigenstates and a conduction band
        % profile object.
        function [eig_system, cond_profile] = solve(obj, ...
                func_carr_distr, eig_system_old)
            
            % Do inputs of eigenstates and/or function handle for
            % carrier distribution, with an eigenstate instance as input
            % argument, exist? Available carrier_distribution frontends:
            % emc_carrier_distribution, eqdist_carrier_distribution or
            % thermal_carrier_distribution.
            if ~exist("func_carr_distr")
                func_carr_distr = @(eig_st) ...
                    eqdist_carrier_distribution. ...
                    generate(obj.sim_const.num_wavefct);
            end
            if ~exist('eig_system_old', 'var')
                eig_system_old = [];
            end
            % Calculate vector with charge density.
            rho = poisson_solver.get_rho(obj.sim_const, ...
                func_carr_distr, eig_system_old);
            
            % Setting up convergence parameter.
            err_new = 2 * obj.error;
            err_old = [];
            err_count = 1;
            
            %% SP loop
            while (err_new > obj.error)
                % Solve Poisson equation for additional electrostatic
                % potential.
                dV = poisson_solver.calc_potential(obj.sim_const, rho);
                % Calculate total potential in J.
                V = obj.sim_const.vec_V_0 + dV;
                % Generate object of conduction band profile.
                cond_profile = conduction_band(obj.sim_const.vec_z_tm, V);
                
                % Calculate eigenenergies and wavefunctions.
                [psi, E, m_E_eff, E_bound_V] = obj.calc_wavefcts(V);
                
                % Check the number of found wavefunctions.
                if (length(E_bound_V) < obj.sim_const.num_wavefct)
                    error(['Could not find the', ...
                        ' appropriate number of wavefunctions!']);
                end
                
                % Returns object of system eigenstates.
                eig_system = obj.gen_eig_syst(psi, E, m_E_eff, ...
                    E_bound_V, cond_profile);
                % Save output parameter of first SP-run.
                if isempty(err_old)
                    cond_init = cond_profile;
                    eig_init = eig_system;
                end
                
                rho_old = rho;
                % Calculate updated charge density from new eigenstates.
                rho = poisson_solver.get_rho(obj.sim_const, ...
                    func_carr_distr, eig_system);
                % Termination condition (relative change in sheet density).
                err_new = sum((rho - rho_old).^2) / sum(rho.^2);
                disp(num2str(err_new));
                
                % Compare error of every second run to prevent infinite
                % loop. Return saved output parameters.
                if mod(err_count-1, 2)
                    if((abs(err_old-err_new) / err_old) < 0.01 || ...
                        err_count > 20)
                        warning(['<strong>RuntimeError:', ...
                            '</strong> Infinite loop. Return initial', ...
                            ' values.'])
                        eig_system = eig_init;
                        cond_profile = cond_init;
                        break
                    end
                end
                
                % Save every error.
                err_old = err_new;
                err_count = err_count + 1;
            end
            
            % Reset minimum energy for root finding algorithm.
            obj.sim_const.Emin = [];
        end
        
        function [psi, E, m_E_eff, E_bound_V] = ...
                calc_wavefcts(obj, V)
            % Calculates eigenstates using transfer
            % matrix method and fzero (root finding algorithm).
            
            % Set/ get energy interval (Emin, Emax) to always find the same
            % wavefunctions.
            if isempty(obj.sim_const.Emin)
                Emin = V(obj.boundary_cond.ind_E_min);
                obj.sim_const.Emin = Emin;
            else
                Emin = obj.sim_const.Emin;
            end
            Emax = Emin + obj.boundary_cond.Eperiod;
            
            % Energy spacing DE
            DE = (Emax - Emin) / obj.num_E;
            opt = optimset('TolX', obj.boundary_cond.dE);
            % Solutions for eigenfunctions and eigenenergies.
            % TM function to be solved.
            disp('Now calculating wavefunctions.');
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
                E_V_tm = ...
                    obj.check_wavefct_interval(E_tm, psi_tm, n_scans);
                psi_scan{n_scans} = psi_tm(obj.sim_const.vec_ind_tm, :);
                % Eigenenergies in J.
                E_scan{n_scans} = E_tm;
                E_V_scan{n_scans} = E_V_tm;
            end
            % Concatenate scans in one vector.
            E_V = [E_V_scan{:}];
            psi = [psi_scan{:}];
            E = [E_scan{:}];
            E_bound_V = obj.calc_E_bound_potential(E_V, psi, V);
            % Sort eigenenergies ascending relative to the
            % conduction band edge. Eq. 6 in Jirauschek and Kubis 2014
            % (https://aip.scitation.org/doi/abs/10.1063/1.4863665)
            [E_bound_V, ind] = sort(E_bound_V);
            psi = psi(:, ind);
            E = E(ind);
            % Delete unbound wavefunctions.
            E(isnan(E_bound_V)) = [];
            psi(:, isnan(E_bound_V)) = [];
            E_bound_V(isnan(E_bound_V)) = [];
            % Calculate vector with effective masses.
            m_E_eff = obj.get_meff_E(V, E, psi);
            
        end
        
        function psi = extend_wavefct(obj, Eh, V, ...
                tm_coeff_AB, n_scans)
            % Extend wavefunction to the considered stacks of QCL periods.
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
            % Calculating in-plane effective mass,
            % see Eq.9 in Jirauschek and Kubis 2014
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
            % Transfer matrix method,
            % see Chapter D.1 in Jirauschek and Kubis 2014.
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
        
        function [tm_coeff] = calc_tmm_result(obj, V, E, n_scans)
            
            % TMM helper function to return boundary solution value of
            % considered wavefunction.
            
            tm_coeff_AB = transfermat(obj, V, E, n_scans);
            index_bc_end = ...
                obj.boundary_cond.get_index_bc_end(size(tm_coeff_AB), ...
                n_scans);
            tm_coeff = real(tm_coeff_AB(index_bc_end));
            
        end
        
        function E_bound = ...
                calc_E_bound_potential(obj, E, psi_tm, V)
            % Relative energy with respect to potential V, measure to find
            % the most strongly bound levels.
            E_bound = zeros(1, length(E));
            % Reduced system with already applied boundaries.
            z = obj.sim_const.vec_z_tm(obj.sim_const.vec_ind_tm);
            for i = 1:length(E)
                E_bound(i) = E(i) - trapz(z, ...
                    psi_tm(:, i).^2.*V(obj.sim_const.vec_ind_tm)');
            end
        end
        
        function E_check = ...
                check_wavefct_interval(obj, E, psi_tm, n_scans)
            z = obj.sim_const.vec_z_tm;
            % Check, whether solution eigenstate is located within
            % the simulation interval [z_L z_R].
            z_L = ...
                obj.boundary_cond.get_z_L(z, n_scans);
            z_R = ...
                obj.boundary_cond.get_z_R(z, n_scans);
            for i = 1:length(E)
                % If 50 % of the probability density is outside of the
                % TMM simulation domain, E is set to NaN to
                % characterize the wavefunction solution as improper.
                if (trapz(z, ...
                        psi_tm(:, i).^2.*z_L'.*z_R') < 0.5)
                    E(i) = NaN;
                end
            end
            E_check = E;
        end
        
        function eig_system = gen_eig_syst(obj, psi, E, m_E_eff, ...
                E_bound_V, cond_profile)
            % Get solution object of system eigenstates.
            
            % Save matrix A containing information about the found
            % eigenenergies, to be written into test.dat.
            A_test = [E', E_bound_V'] / phys_const.e0;
            
            % Select boundary condition to generate system hamiltonian.
            if (strcmp('tb', obj.boundary_cond.name) || ...
                    strcmp('tb-ez', obj.boundary_cond.name))
                % Take the most strongly bound states.
                E = E(1:obj.sim_const.num_wavefct);
                psi = psi(:, 1:obj.sim_const.num_wavefct);
                m_E_eff = m_E_eff(1:obj.sim_const.num_wavefct);
                E_bound_V = E_bound_V(1:obj.sim_const.num_wavefct);
                % Interpolate wavefunctions over 4 periods.
                [psi_eig, E_eig, E_bound_CBO_eig, m_E_eff_eig] = ...
                    interp_wfs(obj, psi, E, E_bound_V, m_E_eff);
                % Use tight-binding base transform.
                Vt = tb_base_transform.get_V_t ...
                    (obj.boundary_cond.ind_z_L, ...
                    obj.boundary_cond.ind_z_R, cond_profile.Vh);
                % Find tight-binding period.
                % Left and right position of tight-binding period.
                z_ind_L = ...
                    obj.sim_const.vec_z_tm(obj.boundary_cond.ind_z_L);
                z_ind_R = ...
                    obj.sim_const.vec_z_tm(obj.boundary_cond.ind_z_R);
                % Return index of tight-binding period.
                ind_period_tb = ...
                    tb_base_transform.find_index_tb_period(z_ind_L, ...
                    z_ind_R, obj.sim_const.vec_z, ...
                    psi_eig, obj.sim_const.num_periods_wf);
                % Calculate hamiltonian with respect to
                % the tight-binding potential.
                hamiltonian_eig = tb_base_transform.calc_hamiltonian( ...
                    obj.sim_const.vec_z, psi_eig, E_eig, ind_period_tb, ...
                    Vt(obj.sim_const.vec_ind_tm), ...
                    cond_profile.Vh(obj.sim_const.vec_ind_tm), ...
                    obj.sim_const.num_periods_wf);
                if (strcmp('tb-ez', obj.boundary_cond.name))
                    % Use ez-base transform after tight-binding.
                    % Define start object of period eigenstates.
                    eig_start = eigenstates(hamiltonian_eig, psi_eig, ...
                        obj.sim_const.vec_z, m_E_eff_eig, ...
                        E_bound_CBO_eig/phys_const.e0);
                    % Finds multiplets marker for ez-transformation.
                    marker = ...
                        ez_base_transform.find_multiplets(eig_start, ...
                        obj.boundary_cond.e_multiplet, ...
                        obj.boundary_cond.d_multiplet);
                    % Transform period eigenstates into ez-states.
                    [hamiltonian_eig, psi_eig, m_E_eff_eig, ~] = ...
                        ez_base_transform.transform(eig_start, marker);
                end
            elseif (strcmp('ext-ez', obj.boundary_cond.name))
                % Interpolate wavefunctions over 4 periods.
                [psi_eig, E_eig, E_bound_CBO_eig, m_E_eff_eig] = ...
                    interp_wfs(obj, psi, E, E_bound_V, m_E_eff);
                % Use ez-base transform.
                % Define start object of period eigenstates.
                eig_start = eigenstates(E_eig, psi_eig, ...
                    obj.sim_const.vec_z, m_E_eff_eig, ...
                    E_bound_CBO_eig/phys_const.e0);
                % Finds multiplets marker for ez-transformation.
                marker = ...
                    ez_base_transform.find_multiplets(eig_start, ...
                    obj.boundary_cond.e_multiplet, ...
                    obj.boundary_cond.d_multiplet);
                % Transform period eigenstates into ez-states.
                [hamiltonian_eig, psi_eig, m_E_eff_eig, ~] = ...
                    ez_base_transform.transform(eig_start, marker);
                % Calculate E_bound_V of the whole system.
                E_bound_CBO_eig = obj.calc_E_bound_potential( ...
                    diag(hamiltonian_eig)', psi_eig, ...
                    cond_profile.Vh/ ...
                    phys_const.e0);
                % Reduce the system to the most strongly bound states.
                [~, ind] = sort(E_bound_CBO_eig);
                ind = ind((obj.sim_const.num_wavefct * ...
                    obj.sim_const.num_periods_wf + 1):end);
                hamiltonian_eig(ind, :) = [];
                hamiltonian_eig(:, ind) = [];
                psi_eig(:, ind) = [];
                m_E_eff_eig(ind) = [];
                E_bound_CBO_eig(ind) = [];
            else
                % Take the most strongly bound states.
                E = E(1:obj.sim_const.num_wavefct);
                psi = psi(:, 1:obj.sim_const.num_wavefct);
                m_E_eff = m_E_eff(1:obj.sim_const.num_wavefct);
                E_bound_V = E_bound_V(1:obj.sim_const.num_wavefct);
                % Interpolate wavefunctions over 4 periods.
                [psi_eig, E_eig, E_bound_CBO_eig, m_E_eff_eig] = ...
                    interp_wfs(obj, psi, E, E_bound_V, m_E_eff);
                hamiltonian_eig = diag(E_eig);
            end
            
            % Generate eigenstate object for the whole system.
            eig_system = eigenstates(hamiltonian_eig, psi_eig, ...
                obj.sim_const.vec_z, m_E_eff_eig, ...
                E_bound_CBO_eig/phys_const.e0);
            % Save matrix A in eigenstate object.
            eig_system.A_test = A_test;
        end
        
        function [psi_interp, E_interp, E_bound_CBO_interp, ...
                m_E_eff_interp] = interp_wfs(obj, psi, E, ...
                E_bound_V, m_E_eff)
            % Interpolate wavefunctions
            % over the number of considered periods.
            
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
                length(E), 1);
            psi_interp = zeros(length(obj.sim_const.vec_z), ...
                obj.sim_const.num_periods_wf*length(E));
            m_E_eff_interp = ...
                zeros(obj.sim_const.num_periods_wf* ...
                length(E), 1);
            E_bound_CBO_interp = ...
                zeros(obj.sim_const.num_periods_wf* ...
                length(E), 1);
            for ni = 1:length(E)
                for n_period = 0:obj.sim_const.num_periods_wf - 1
                    psi_interp(:, ...
                        ni+n_period*length(E)) = ...
                        psi_interp_cell{ni}(:, 1 + n_period);
                    E_interp(ni+n_period*length(E)) = ...
                        E_interp_cell{ni}(1 + n_period) / phys_const.e0;
                    m_E_eff_interp(ni+n_period* ...
                        length(E)) = ...
                        m_E_eff(ni);
                    E_bound_CBO_interp(ni+n_period* ...
                        length(E)) = ...
                        E_bound_V(ni);
                end
            end
        end
    end
end
