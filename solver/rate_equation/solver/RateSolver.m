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

classdef RateSolver < handle
    % Solves the rate equations self-consistently for the two central
    % QCL-periods and assuming periodic boundary conditions.

    properties
        scattering_rates % cell-array: Array containing scattering mechanism objects.
        tau_inv % matrix: Transition rates of the two central periods.
        tau_inv_k % 3-d array: K-resolved transition rates of the two central periods.
        pure_dephasing % matrix: Pure dephasing rates of the two central periods.
        pure_dephasing_k % 3-d array: K-resolved pure dephasing rates of the two central periods.
        electron_density % vector: Sheet densities of each subband [1/m^2].
        num_states % scalar: Number of considered subbands.
    end

    properties (SetAccess = private)
        optics_flag % logical: Flag enabling simulation with optical field.
    end

    methods (Access = public)
        function obj = RateSolver(scattering_rates)
            % Constructs an object of type RateSolver.
            %
            % Syntax:
            %   obj = RateSolver(scattering_rates)
            %
            % Input Arguments:
            %   scattering_rates (cell-array): Cell array containing the
            %     scattering rate objects to be considered in the
            %     simulation.

            obj.optics_flag = 0; %< bla
            obj.set_scattering_rates(scattering_rates);
            obj.num_states = obj.scattering_rates{1}.num_states;
            obj.electron_density = scattering_rates{1}.electron_density;

            obj.tau_inv = zeros(obj.num_states);
            obj.tau_inv_k = zeros(size(obj.scattering_rates{1}.W));
            obj.pure_dephasing = zeros(obj.num_states);
            obj.pure_dephasing_k = zeros(size( ...
                obj.scattering_rates{1}.pure_dephasing_k));
        end

        function carr_dist = solve(obj, options)
            % Determines subband populations in steady-state by iteratively
            % solving the rate equations.
            %
            % Syntax:
            %   carr_dist = solve(obj)
            %   carr_dist = solve(__, name, value)
            %
            % Name Value Arguments:
            %   ns_init: Vector containing the normalized, initial subband
            %     populations of one QCL-period. A uniform distribution is
            %     assumed by default.
            %   nmax (15 (default) | scalar): Maximum number of iterations.
            %   relTol (0.01 (default) | scalar): Desired relative error
            %     between two subsequent iterations to terminate the loop.
            %   parallel (false (default) | logical): Flag
            %     enabling parallel computation.
            %
            % Output Arguments:
            %   carr_dist: Returns object of type carrier_distribution.

            arguments
                obj
                options.ns_init
                options.nmax = 15
                options.relTol = 1e-2
                options.parallel = false
            end

            keys = fieldnames(options);
            values = struct2cell(options);
            options = {};
            for i = 1:length(keys)
                options{2*i-1} = keys{i};
                options{2*i} = values{i};
            end

            if obj.optics_flag
                carr_dist = solve_with_optics(obj, options{:});
            else
                carr_dist = solve_no_optics(obj, options{:});
            end
        end

        function carr_dist = solve_no_optics(obj, options)
            % Determines subband populations in steady-state without
            % inclusion of the optical cavity field.
            %
            % Syntax:
            %   carr_dist = solve(obj)
            %   carr_dist = solve(__, name, value)
            %
            % Name Value Arguments:
            %   ns_init (vector): Normalized initial subband populations in
            %     one period.
            %   nmax (scalar): Maximum number of iterations.
            %   relTol (scalar): Desired relative error between two
            %     subsequent iterations to terminate the loop.
            %   parallel (logical): Flag enabling parallel computation.
            %
            % Output Arguments:
            %   carr_dist: Returns object of type carrier_distribution.

            arguments
                obj
                options.ns_init
                options.nmax
                options.relTol
                options.parallel
            end

            max_iteration = options.nmax;
            parallel = options.parallel;
            relTol = options.relTol;
            if any(fieldnames(options) == "ns_init")
                occ = reshape(options.ns_init, [], 1);
                obj.electron_density = repmat(occ, 4, 1) .* ...
                    sum(obj.electron_density) / 4;
            end

            fprintf('Solving rate equations ...\n\n')
            time_rates = tic;

            ns_old = obj.electron_density;

            % self-consistent loop
            for iter = 1:max_iteration

                % calculate transition rates
                obj.calc_scattering_rates(ns_old, parallel, iter);

                % calculate new sheet densities
                ns_new = calc_sheet_dens(obj);
                ns_new = 0.5 * (ns_new + ns_old);

                % determine new residuals
                err = norm(ns_old-ns_new) / norm(abs(ns_old));

                fprintf(['(', num2str(iter), '/', num2str(max_iteration), ...
                    ') residual: ', num2str(err), '\n\n'])
                if (err < relTol)
                    fprintf(['Rate solver reached desired residual ', ...
                        'of ', num2str(relTol), '. Elapsed time: ', ...
                        num2str(toc(time_rates)), 's\n\n'])
                    break
                end
                if iter == max_iteration
                    fprintf(['Rate solver did not reach desired ', ...
                        'residual of ', num2str(relTol), '. Elapsed ', ...
                        'time: ', num2str(toc(time_rates)), 's\n\n'])
                end
                ns_old = ns_new;
            end

            obj.electron_density = ns_new;

            % calculate electron distributions and occupations for
            % a single period
            [dist, Ekin] = obj.calc_distribution();
            occupations = obj.electron_density(1:obj.num_states/2) / ...
                sum(obj.electron_density(1:obj.num_states/2));
            % k-resolved sheet density
            nS_k = obj.calc_sheet_dens_k(dist);
            % return carrier distribution object
            carr_dist = carrier_distribution(dist, ...
                occupations, Ekin/phys_const.e0, nS_k);
        end

        function carr_dist = solve_with_optics(obj, options)
            % Determines subband populations in steady-state with inclusion
            % of the optical cavity field.
            %
            % Syntax:
            %   carr_dist = solve(obj)
            %   carr_dist = solve(__, name, value)
            %
            % Name Value Arguments:
            %   ns_init (vector): Normalized initial subband populations in
            %     one period.
            %   nmax (scalar): Maximum number of iterations.
            %   relTol (scalar): Desired relative error between two
            %     subsequent iterations to terminate the loop.
            %   parallel (logical): Flag enabling parallel computation.
            %
            % Output Arguments:
            %   carr_dist: Returns object of type carrier_distribution.

            arguments
                obj
                options.ns_init
                options.nmax
                options.relTol
                options.parallel
            end

            parallel = options.parallel;
            max_iterations = options.nmax;
            error = options.relTol;

            keys = fieldnames(options);
            values = struct2cell(options);
            options = {};
            for i = 1:length(keys)
                options{2*i-1} = keys{i};
                options{2*i} = values{i};
            end

            % In the first iteration calculate transition rates and
            % occupations without considering the optical field.
            obj.solve_no_optics(options{:});
            ns_old = obj.electron_density;
            I = [];

            % After obtaining a converged solution for the transition rates
            % and occupations, include the optical field in the simulation.
            fprintf('Coupling rate equations with optical field ...\n\n')
            for iter = 1:max_iterations
                % Calculate steady-state solution of the optical field
                % and the optical transition rates, which also results
                % in new sheet densities.
                [ns, ~] = calc_optical_field(obj, ns_old, I);
                if ns == ns_old
                    ns_new = ns;
                    fprintf('\n\n')
                    break
                end
                % Re-calculate the non-radiative transition rates for the
                % new sheet densities.
                obj.electron_density = ns;
                obj.calc_scattering_rates(ns, parallel, iter);

                % Add optical transition rates to total transition rates
                obj.scattering_rates{end}.update(ns);
                rates = obj.tau_inv + obj.scattering_rates{end}.calculate();
                % Calculate new sheet densities
                ns_new = calc_sheet_dens(obj, ns(1:obj.num_states), rates);
                ns_new = 0.5 * (ns_new + ns_old);

                % Check convergence
                err = norm(ns_old-ns_new) / norm(abs(ns_old));
                fprintf(['(', num2str(iter), '/', num2str( ...
                    max_iterations), ') residual: ', num2str(err), ...
                    '\n\n'])
                if (err < error)
                    fprintf(['Rate-Optics solver reached desired ', ...
                        'residual of ', num2str(error), '.\n\n'])
                    break
                end
                if iter == max_iterations
                    fprintf(['Rate-Optics solver did not ', ...
                        'reach desired residual of ', num2str(error), ...
                        '.\n\n'])
                end
                ns_old = ns_new;
            end

            obj.electron_density = ns_new;

            % Calculate electron distributions and occupations for
            % a single period
            [dist, Ekin] = obj.calc_distribution();
            occupations = obj.electron_density(1:obj.num_states/2) / ...
                sum(obj.electron_density(1:obj.num_states/2));
            % k-resolved sheet density
            nS_k = obj.calc_sheet_dens_k(dist);
            % return carrier distribution object
            carr_dist = carrier_distribution(dist, ...
                occupations, Ekin/phys_const.e0, nS_k);
        end

        function calc_scattering_rates(obj, ns_old, parallel, iteration)
            % Calculates the scattering rates and pure dephasing rates
            % between subbands of the two central periods.
            %
            % Syntax:
            %   calc_scattering_rates(obj, ns_old, parallel, iteration)
            %
            % Input Arguments:
            %   ns_old (vector): Sheet densities from previous iteration.
            %   parallel (logical): Flag enabling parallel computation.
            %   iteration (scalar): Current iteration index.

            % reset rates
            obj.tau_inv = zeros(size(obj.tau_inv));
            obj.tau_inv_k = zeros(size(obj.tau_inv_k));
            obj.pure_dephasing = zeros(size(obj.pure_dephasing));
            obj.pure_dephasing_k = zeros(size(obj.pure_dephasing_k));

            % compute scattering rates for each scattering mechanism
            % in obj.scattering_rates
            for i = 1:length(obj.scattering_rates) - obj.optics_flag
                tclock = tic;
                % determine if k-resolved transition rates have to
                % be recalculated
                re_calc = obj.re_calculate_rates(iteration, ...
                    obj.scattering_rates{i}.Name, ...
                    obj.scattering_rates{1}.screening_model);
                % update dynamical properties of all
                % scattering rate objects
                if obj.scattering_rates{i}.Name == "tunneling"
                    % k-averaged
                    inv_lifetimes = sum(obj.total_module_rate( ...
                        obj.tau_inv), 2);
                    obj.scattering_rates{i}.update(ns_old, ...
                        inv_lifetimes, obj.pure_dephasing);
                    % k-resolved
                    inv_lifetimes = reshape(sum(obj. ...
                        total_module_rate(obj.tau_inv_k), 2), ...
                        size(obj.tau_inv, 1), []);
                    obj.scattering_rates{i}.update(ns_old, ...
                        inv_lifetimes, obj.pure_dephasing_k);
                else
                    % set new sheet densities and quasi-Fermi levels
                    obj.scattering_rates{i}.update(ns_old);
                    % calculate pure dephasing rates
                    if parallel
                        obj.pure_dephasing = obj.pure_dephasing + ...
                            obj.scattering_rates{i}. ...
                            calulate_pure_dephasing_rate_parallel(re_calc);
                    else
                        obj.pure_dephasing = obj.pure_dephasing + ...
                            obj.scattering_rates{i}. ...
                            calulate_pure_dephasing_rate(re_calc);
                    end
                    obj.pure_dephasing_k = obj.pure_dephasing_k + ...
                        obj.scattering_rates{i}.pure_dephasing_k;
                end
                % calculate scattering rate
                if parallel
                    obj.scattering_rates{i}.calculate_parallel(re_calc);
                else
                    obj.scattering_rates{i}.calculate(re_calc);
                end
                obj.tau_inv = obj.tau_inv + ...
                    obj.scattering_rates{i}.tau_inv;
                obj.tau_inv_k = obj.tau_inv_k + ...
                    obj.scattering_rates{i}.W;
                fprintf(['%s: ', num2str(toc(tclock)), 's\n'], ...
                    obj.scattering_rates{i}.Name)
            end
        end

        function [ns, I] = calc_optical_field(obj, ns_init, I_init)
            % Calculates the optical field and optical transition rates.
            %
            % Syntax:
            %   calc_scattering_rates(obj, ns_old, parallel, iteration)
            %
            % Input Arguments:
            %   ns_init (vector): Normalized initial subband populations in
            %     one period.
            %   I_init (vector): Initial intensities of the cavity modes.
            %
            % Output Arguments:
            %   ns (vector): New sheet densities.
            %   I (vector): New intensities of the cavity modes.

            tstart = tic;
            if isempty(I_init) || nargin == 2
                I = 300 + (500 - 300) .* rand( ...
                    size(obj.scattering_rates{end}.freq));
            else
                I = I_init;
            end

            % set dephasing rates
            lvl_broad = sum(obj.total_module_rate(obj.tau_inv), 2);
            obj.scattering_rates{end}.set_dephasing_rates( ...
                lvl_broad, obj.pure_dephasing);
            lvl_broad = sum(obj.total_module_rate(obj.tau_inv_k), 2);
            obj.scattering_rates{end}.set_dephasing_rates( ...
                lvl_broad, obj.pure_dephasing_k);

            % time evolution of the optical field
            optics_iteration = 3000;
            ns = ns_init;
            log_string = [];
            for i = 1:optics_iteration
                % Calculate new  mode intensities
                obj.scattering_rates{end}.update(ns);
                [dist, ~] = obj.calc_distribution(ns);
                ns_k = obj.calc_sheet_dens_k(dist, ns);
                [I_new, gain] = obj.scattering_rates{end}. ...
                    calc_intensity(I, ns_k);
                if max(gain) < 0
                    log_string = '(Gain < 0)';
                    break
                end

                % Add optical transition rates to the total scattering
                % rates and calculate new sheet densities. Note that the
                % non-radiative transition rates are kept constant in this
                % loop to speed up convergence.
                rates = obj.tau_inv + ...
                    obj.scattering_rates{end}.calculate();
                ns_new = calc_sheet_dens(obj, ns(1:obj.num_states), rates);

                % Check convergence of the optical field
                err = norm(abs(I_new-I)) ./ norm(I);
                I = I_new;
                ns = ns_new;
                if err < 1e-3 || i == optics_iteration
                    log_string = '(Sub-iterations: ' + string(i) + ...
                        '/' + string(optics_iteration) + ')';
                    break
                end
            end
            fprintf(['%s: ', num2str(toc(tstart)), 's %s\n'], ...
                obj.scattering_rates{end}.Name, log_string)
        end

        function update_eigenstates(obj, eigen)
            % Sets new eigenenergies, effective masses and wavefunctions
            % for every scattering mechanism object.
            %
            % Syntax:
            %   update_eigenstates(obj, eigen)
            %
            % Input Arguments:
            %   eigen (eigenstates-object): Eigenstates object containing
            %     the updated energies, wavefunctions and effective masses
            %     from the Schroedinger-Poisson solver.

            for i = 1:length(obj.scattering_rates)
                obj.scattering_rates{i}.set_eigenstates(eigen);
            end
        end

        function deph_rate = get_dephasing_rates(obj, carr_dist, ...
                scat_rates_mech)
            % Returns object of type dephasing_rates for postprocessing.
            %
            % Syntax:
            %   deph_rate = get_dephasing_rates(obj, carr_dist)
            %   deph_rate = get_dephasing_rates(obj, carr_dist, scat_rates_mech)
            %
            % Input Arguments:
            %   carr_dist (carrier_distribution-object): Information about
            %     the carrier distribution.
            %   scat_rates_mech (vector): Vector containing names of
            %     scattering mechanisms which should contribute to dephasing.
            %     Valid names are ``interface roughness scattering``,
            %     ``impurity scattering``, ``LO phonon absorption``,
            %     ``LO phonon emission`` and ``acoustic phonon scattering``.
            %
            % Output Arguments:
            %   deph_rate (dephasing_rates-object): Returns object of type
            %     dephasing_rates.

            if nargin < 3
                mech_indices = 1:length(obj.scattering_rates) - ...
                    obj.optics_flag;
            else
                mech_indices = zeros(1, length(scat_rates_mech));
                c = 1;
                for ind = 1:length(obj.scattering_rates)
                    name = obj.scattering_rates{ind}.Name;
                    if ismember(name, scat_rates_mech)
                        mech_indices(c) = ind;
                        c = c + 1;
                    end
                end
            end

            num_k = length(obj.scattering_rates{1}.k); % number of k states
            num_E = length(obj.scattering_rates{1}.E); % number of subbands

            keys_deph = {'interface roughness', 'impurity', ...
                'optical absorption', 'optical emission', ...
                'acoustic phonon', 'total'};
            % initialize array with pure dephasing rates for each mechanism
            pd_mech = {};
            for j = 1:length(keys_deph)
                pd_mech{j} = zeros(num_E, num_E, num_k);
            end
            pure_dephasing_mech = containers.Map(keys_deph, pd_mech);

            keys_rates = {'interface roughness scattering', ...
                'impurity scattering', 'LO phonon absorption', ...
                'LO phonon emission', 'acoustic phonon scattering'};

            % map names of scattering rates to names of dephasing rates
            map_rates_deph = containers.Map(keys_rates, keys_deph(1:end-1));

            % level broadening for states in one module (=period)
            lvl_broad_module = zeros(obj.num_states/2, num_k);

            pd_total = zeros(num_E, num_E, num_k);
            % collect k-resolved pure-dephasing rates and transition rates
            % from scattering rate objects
            n_wfs = obj.num_states / 2;
            for i = mech_indices
                key = obj.scattering_rates{i}.Name;
                pd = zeros(num_E, num_E, num_k);
                if isKey(map_rates_deph, key)
                    % get dephasing rates as block matrices
                    r_middle = obj.scattering_rates{i}.pure_dephasing_k( ...
                        1:n_wfs, 1:n_wfs, :);
                    r_middle = (r_middle + obj.scattering_rates{i}. ...
                        pure_dephasing_k(n_wfs+1:2*n_wfs, ...
                        n_wfs+1:2*n_wfs, :)) / 2;
                    r_left = obj.scattering_rates{i}.pure_dephasing_k( ...
                        n_wfs+1:2*n_wfs, 1:n_wfs, :);
                    r_right = obj.scattering_rates{i}.pure_dephasing_k( ...
                        1:n_wfs, n_wfs+1:2*n_wfs, :);
                    % fill pure dephasing rate matrix with block matrices
                    for j = 1:4
                        pd((j - 1)*n_wfs+1:j*n_wfs, ...
                            (j - 1)*n_wfs+1:j*n_wfs, :) = r_middle;
                        if j ~= 1
                            pd((j - 1)*n_wfs+1:j*n_wfs, ...
                                (j - 2)*n_wfs+1:(j - 1)*n_wfs, :) = r_left;
                        end
                        if j ~= 4
                            pd((j - 1)*n_wfs+1:j*n_wfs, ...
                                j*n_wfs+1:(j + 1)*n_wfs, :) = r_right;
                        end
                    end
                    % pure dephasing rates divided into the different
                    % scattering mechanisms
                    pure_dephasing_mech(map_rates_deph(key)) = pd;
                    % pure dephasing rates summed over all the different
                    % scattering mechanisms
                    pd_total = pd_total + pd;
                end
                % k resolved transistion rates
                W_ikj = obj.scattering_rates{i}.W;
                m = obj.num_states / 2;
                W_ikj = W_ikj(m+1:end, 1:m, :) + W_ikj(1:m, m+1:end, :) + ...
                    0.5 * (W_ikj(1:m, 1:m, :) + W_ikj(m+1:end, m+1:end, :));
                lvl_broad_module = lvl_broad_module + ...
                    reshape(sum(W_ikj, 2), m, num_k);
            end
            pure_dephasing_mech('total') = pd_total;

            num_periods = num_E / (obj.num_states / 2);
            lvl_broadening = repmat(lvl_broad_module, num_periods, 1);

            deph_rate = dephasing_rates(pure_dephasing_mech, ...
                lvl_broadening, carr_dist, keys_deph);
        end

        function scat_rates = get_scattering_rates(obj)
            % Returns object of class scattering_rates for postprocessing.
            %
            % Syntax:
            %   scat_rates = get_scattering_rates(obj)
            %
            % Output Arguments:
            %   scat_rates (scattering_rates-object): Returns object of
            %     type scattering_rates.

            keys_scat = {'electron-electron', ...
                'optical absorption', 'optical emission', 'tunneling', ...
                'alloy disorder', 'acoustic phonon', 'impurity', ...
                'interface roughness'};

            keys_rates = {'electron-electron scattering', ...
                'LO phonon absorption', 'LO phonon emission', ...
                'tunneling', 'alloy scattering', ...
                'acoustic phonon scattering', 'impurity scattering', ...
                'interface roughness scattering'};

            % map names of scattering rates to names of dephasing rates
            map_rates = containers.Map(keys_rates, keys_scat);

            n_wfs = obj.num_states / 2;
            nE = length(obj.scattering_rates{1}.E);
            scat = {};
            for j = 1:length(keys_scat)
                scat{j} = zeros(nE, nE);
            end
            scattering = containers.Map(keys_scat, scat);

            for i = 1:length(obj.scattering_rates)
                key = obj.scattering_rates{i}.Name;
                if isKey(map_rates, key)
                    rate = zeros(nE, nE);

                    % get scattering rates as block matrices
                    r_middle = obj.scattering_rates{i}.tau_inv( ...
                        1:n_wfs, 1:n_wfs);
                    r_middle = (r_middle + ...
                        obj.scattering_rates{i}.tau_inv( ...
                        n_wfs+1:2*n_wfs, n_wfs+1:2*n_wfs)) / 2;
                    r_left = obj.scattering_rates{i}.tau_inv( ...
                        n_wfs+1:2*n_wfs, 1:n_wfs);
                    r_right = obj.scattering_rates{i}.tau_inv( ...
                        1:n_wfs, n_wfs+1:2*n_wfs);

                    % fill the rate matrix for the full structure
                    % (4 periods) with the block matrices.
                    for j = 1:4
                        rate((j - 1)*n_wfs+1:j*n_wfs, ...
                            (j - 1)*n_wfs+1:j*n_wfs) = r_middle;
                        if j ~= 1
                            rate((j - 1)*n_wfs+1:j*n_wfs, ...
                                (j - 2)*n_wfs+1:(j - 1)*n_wfs) = r_left;
                        end
                        if j ~= 4
                            rate((j - 1)*n_wfs+1:j*n_wfs, ...
                                j*n_wfs+1:(j + 1)*n_wfs) = r_right;
                        end
                    end
                    scattering(map_rates(key)) = rate;
                end
            end
            scat_rates = scattering_rates(scattering);
        end

        function field = get_optical_field(obj)
            % Returns object of class optical_field for postprocessing.
            %
            % Syntax:
            %   field = get_optical_field(obj)
            %
            % Output Arguments:
            %   field (optical_field-object): Returns object of type
            %     optical_field.

            if class(obj.scattering_rates{end}) == "OpticalRate"
                intensity = obj.scattering_rates{end}.intensity;
                freq = obj.scattering_rates{end}.freq;
                field = optical_field(intensity, freq);
            else
                error("Optical rate not found!")
            end
        end

        function J = calc_current_density(obj)
            % Calculates current density across interface between two
            % QCL-periods.
            %
            % Syntax:
            %   J = calc_current_density(obj)
            %
            % Output Arguments:
            %   J (scalar): Current density [A/m^2].

            n = size(obj.tau_inv, 1) / 2;
            rate_left = obj.tau_inv(n+1:2*n, 1:n);
            rate_right = obj.tau_inv(1:n, n+1:2*n);
            n2D = obj.electron_density(1:n);
            J = -sum(sum((rate_right - rate_left).*n2D)) .* phys_const.e0;
        end
    end

    methods (Access = private)
        function ns_all = calc_sheet_dens(obj, ns_old, tau_inv)
            % Calculates new sheet densities.
            %
            % Syntax:
            %   ns_all = calc_sheet_dens(obj)
            %   ns_all = calc_sheet_dens(obj, ns_old, tau_inv)
            %
            % Input Arguments:
            %   ns_old (vector): Sheet densities for the two central
            %     periods [1/m^2].
            %   tau_inv (matrix): Matrix containing the transition rates
            %     for the two central periods.
            %
            % Output Arguments:
            %   ns_all (vector): Sheet densities for all four periods
            %     [1/m^2].

            if nargin < 2
                tau_inv = obj.tau_inv;
                ns_old = obj.scattering_rates{1}.electron_density( ...
                    obj.num_states/2+1:3*obj.num_states/2);
            end

            % sheet densities for all modules
            ns_all = zeros(length(obj.scattering_rates{1}.E), 1);
            % sheet densities for a single module
            ns_module = zeros(obj.num_states/2, 1);

            % transition rates
            tau_inv_tot = obj.total_module_rate(tau_inv);

            % calculate new occupation for each level
            for i = 1:(obj.num_states / 2)
                ns_module(i) = sum(tau_inv_tot(:, i) ...
                    .*ns_old) / sum(tau_inv_tot(i, :));
            end
            % conserve total electron density
            ns_module = ns_module / sum(ns_module) * sum(ns_old) / 2;
            % assume same occupations in each period
            num_periods = length(obj.scattering_rates{1}.E) ...
                / length(ns_module);
            for i = 1:num_periods
                ns_all((i - 1)*length(ns_module)+ ...
                    1:i*length(ns_module)) = ns_module;
            end
        end

        function [dist, E_kin] = calc_distribution(obj, n_sheet)
            % Returns carrier distribution in k-space for each subband
            % according to Fermi-Dirac statistics.
            %
            % Syntax:
            %   [dist, E_kin] = calc_distribution(obj)
            %   [dist, E_kin] = calc_distribution(obj, n_sheet)
            %
            % Input Arguments:
            %   n_sheet (vector): Sheet densities [1/m^2]. If n_sheet is
            %     not provided as input, then the values of the class
            %     property electron_density are used instead.
            %
            % Output Arguments:
            %   dist (matrix): Matrix, where each column contains the
            %     carrier distribution of a specific subband in k-space.
            %   E_kin (matrix): Matrix containing the kinetic energy
            %     vectors E_kin=(hbar*k)^2/(2m*) for each subband as
            %     columns.

            % extract relevant quantities from first scattering rate object
            rate_obj = obj.scattering_rates{1};
            k = rate_obj.k;
            E = rate_obj.E;
            T_e = rate_obj.T_e;
            mEff = rate_obj.mEff;
            num = obj.num_states / 2;

            % calculate new quasi Fermi-levels
            if nargin < 2
                n_sheet = obj.electron_density;
            end
            E_F = rate_obj.calc_fermi_energy(E(1:num), mEff(1:num), ...
                n_sheet(1:num), T_e(1:num));

            % calculate distribution
            dist = zeros(length(k), num);
            E_kin = zeros(length(k), num);
            for i = 1:num
                E_kin(:, i) = phys_const.hbar^2 .* ...
                    reshape(k.^2, [], 1) ./ (2 * mEff(i));
                dist(:, i) = 1 ./ (exp((E(i) + E_kin(:, i) - E_F(i)) ...
                    ./(phys_const.kB * T_e(i))) + 1);
            end
        end

        function nS_k = calc_sheet_dens_k(obj, dist, n_sheet)
            % Calculates k-resolved sheet densities for all subbands in a
            % single QCL-period.
            %
            % Syntax:
            %   nS_k = calc_sheet_dens_k(obj, dist)
            %   nS_k = calc_sheet_dens_k(obj, dist, ns)
            %
            % Input Arguments:
            %   dist (matrix): Matrix, where each column contains the
            %     carrier distribution of a specific subband in k-space.
            %   n_sheet (vector): Sheet densities [1/m^2]. If n_sheet is
            %     not provided as input, then the values of the class
            %     property electron_density are used instead.
            %
            % Output Arguments:
            %   nS_k (matrix): Matrix containing the k-resolved sheet
            %     densities [1/m^2] for each subband as columns.


            if nargin < 3
                n_sheet = obj.electron_density;
            end
            num_wfs = obj.num_states / 2;
            nS_k = dist ./ sum(dist) .* reshape(n_sheet(1:num_wfs), 1, []);
        end

        function set_scattering_rates(obj, scattering_rates)
            % Sets correctly ordered cell-array of scattering rate
            % mechanism objects. The elements will be ordered such that
            % tunneling- and optial_rate-object (if present) are moved to
            % the end of the cell-array.
            %
            % Syntax:
            %   set_scattering_rates(obj, scattering_rates)
            %
            % Input Arguments:
            %   scattering_rates (cell-array): Array containing scattering
            %     mechanism objects.

            index = reshape(1:length(scattering_rates), ...
                size(scattering_rates));
            index_tunnel = 0;
            index_optfield = 0;

            % get index of tunneling and optical_field rate
            for i = 1:length(scattering_rates)
                if scattering_rates{i}.Name == "tunneling"
                    index_tunnel = i;
                    index(index == i) = [];
                end
                if scattering_rates{i}.Name == "optical field"
                    index_optfield = i;
                    index(index == i) = [];
                    obj.optics_flag = 1;
                end
            end

            % move tunneling and optical field rate object to the end
            % of the index array
            if index_tunnel
                index(end+1) = index_tunnel;
            end
            if index_optfield
                index(end+1) = index_optfield;
            end

            obj.scattering_rates = scattering_rates(index);
        end
    end

    methods (Static)
        function rate_tot = total_module_rate(rate)
            % Calculates the total transition rates for two QCL-periods
            % assuming periodic boundary conditions for the transition
            % rates.
            %
            % Syntax:
            %   rate_tot = total_module_rate(rate)
            %
            % Input Arguments:
            %   rate (matrix | 3-d array): Transition rates (either
            %     k-averaged or k-resolved) for two QCL-periods [1/s].
            %
            % Output Arguments:
            %   rate_tot (matrix | 3-d array): Transition rates for the two
            %     central periods with periodic boundary conditions applied.

            middle = size(rate, 1) / 2;

            if length(size(rate)) == 2 % k-independent rates
                % remove diagonal elements
                rate = rate - diag(diag(rate));
                % intermodule scattering rates
                r_left = rate(middle+1:end, 1:middle);
                r_right = rate(1:middle, middle+1:end);
                % intramodule scattering rates
                r_middle = 0.5 * (rate(1:middle, 1:middle) + ...
                    rate(middle+1:end, middle+1:end));
            elseif length(size(rate)) == 3 % k-dependent rates
                % intermodule scattering rates
                r_left = rate(middle+1:end, 1:middle, :);
                r_right = rate(1:middle, middle+1:end, :);
                % intramodule scattering rates
                r_middle = 0.5 * (rate(1:middle, 1:middle, :) + ...
                    rate(middle+1:end, middle+1:end, :));
            end

            rate_tot = [r_middle, r_left + r_right; ...
                r_left + r_right, r_middle];
        end

        function cal = re_calculate_rates(iteration, mechanism, ...
                screening_model)
            % Evaluates if scattering rates have to  be recalculated
            % depending on the screening model and the current iteration
            % index.
            %
            % Syntax:
            %   cal = re_calculate_rates(iteration, mechanism, screening_model)
            %
            % Input Arguments:
            %   iteration (scalar): Current iteration index.
            %   mechanism (char): Name of scattering mechanism.
            %   screening_model (char): Name of the screening model.
            %
            % Output Arguments:
            %   cal (logical): Flag specifying if scattering rates have to
            %     be recalculated.

            if iteration ~= 1
                if any(strcmp(mechanism, {'impurity scattering', ...
                        'electron-electron scattering', ['LO phonon ', ...
                        'emission'], 'LO phonon absorption'}))
                    if any(strcmp(screening_model, ["debye", ...
                            "thomas-fermi"]))
                        cal = 0;
                    else
                        cal = 1;
                    end
                else
                    cal = 0;
                end
            else
                cal = 1;
            end
        end
    end
end
