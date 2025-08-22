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

classdef setup_mb < handle
    % Writes input files for the dynamic maxwell-bloch solvers
    % mbsolve (full-wave) and mb_rwa (rotating-wave approximation).

    properties (Constant)
        % Container map of the different cavity shapes.
        shape = containers.Map ...
            ({'Fabry-Perot', 'ring', 'ring_right_travel'}, [0, 1, 2]);
    end

    properties (SetAccess = private)
        gain % sim_gain-object: Contains information about optical gain.
        sc % scattering_rates-object: Contains scattering rates.
        OAB = 3 % scalar: Order of Adams-Bashforth.
        Int = 1e4 % scalar: Seed intensity [W/m^2].
        num_reg = 1 % scalar: Number of simulated regions.
        r_trip = 100 % scalar: Number of roundtrips.
        N_z_min = 0 % scalar: Starting grid point for E-field plot.
        dN_z = 2 % scalar: Grid point steps for E-field plot.
        dN_t = 4 % scalar: Time grid point steps for E-field plot.
        bc = 0 % scalar: Indicator for cavity shape.
        dens_carrier % scalar: Space-averaged electron density.
        GVD = 0 % scalar: Group velocity dispersion.
        sat_abs_intens = 0 % scalar: Saturable absorber saturation intensity.
        sat_abs_rt = 0 % scalar: Saturable absorber (loss, recovery time).
        T_bloch = 0 % scalar: Electron temp. for Bloch gain - PRL. 127, 093902.
        rel_permeab = 1 % scalar: Relative permeability.
        i_wf = [] % vector: Period of interest (reduced system).
        t % scalar: Simulation time for mbsolve
        num_wfs % scalar: Number of wavefunctions of reduced system.
        l_wg % scalar: Waveguide length.
        rel_perm % scalar: Relative permittivity.
        n_eff % scalar: Refractive index.
        a_power % scalar: Power absorption coefficient.
        overlap_factor % scalar: Overlap factor.
        D = 0 % scalar: Diffusion coefficient.
        r10 % scalar: Reflection coefficient from region 1 to 0.
        r12 % scalar: Reflection coefficient from region 1 to 2.
        t10 % scalar: Transmission coefficient from region 1 to 0.
        t12 % scalar: Transmission coefficient from region 1 to 2.
        l_period % scalar: Period length.
        A_act % scalar: Cross-sectional area of the active region.
        freq_range % vector: Min and max frequency.
        f_c % scalar: Center frequency in Hz.
        N_z % scalar: Number of grid points.
        N_z_max % scalar: Last grid point for E-field plot.
        N_t % scalar: Number of time grid points for E-field plot.
        N_t_min % Starting time grid point for E-field plot.
        scat_rates_2lvl % matrix: Scattering rates for 2 lvl representation
        deph_rate_2lvl % scalar: Dephasing rate for 2 lvl representation.
    end

    properties
        dir % char | string: Folder where input file should be stored.
        name % char | string: Name of input file.
        mb_input_data % input_data-object: Stationary simulation results.
        injection % vector:
    end

    methods
        function obj = setup_mb(dir, s, gain, sc)
            % Constructs an object of type setup_mb.
            %
            % Syntax:
            %   obj = setup_mb(dir, s, gain, sc)
            %
            % Input Arguments:
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %   gain (sim_gain-object): Contains description of the optical
            %     gain.
            %   sc (scattering_rates): Contains information about the
            %     scattering rates.

            obj.name = s.name;
            obj.dir = dir;
            % Set gain with object of class emc_sim_gain.
            obj.gain = gain;
            % Set scattering rates with object of class
            % emc_scattering_rates
            obj.sc = sc;
            % Set num_wfs from the property of gain.eigen.
            obj.num_wfs = gain.eigen.num_wfs;
            % Set relative permittivity.
            obj.rel_perm = gain.device.rel_permittivity;
            % Set the refractive index.
            obj.n_eff = gain.device.n_eff;
            % Set frequency range
            obj.freq_range = [s.fmin, s.fmax];
            % Calculate the center frequency with default settings.
            obj.calc_fc(s.fmin, s.fmax, 'total');
            % Set the waveguide length.
            obj.l_wg = gain.device.waveguide.l_waveguide;
            % Set the power absorption coefficient.
            obj.a_power = gain.device.waveguide.a_power;
            % Set the overlap factor.
            obj.overlap_factor = gain.device.waveguide.overlap_factor;
            % Set the space-averaged electron density.
            obj.dens_carrier = gain.device.dens_carrier;
            % Set reflection and transmission coefficients.
            obj.set_facet_refl_coeff( ...
                gain.device.waveguide.refl_coeff_left, ...
                gain.device.waveguide.refl_coeff_right);
            % Set the period length.
            obj.l_period = gain.device.l_period;
            % Set default injection current.
            obj.injection = zeros(1, obj.num_wfs);
            % Set cross-sectional area of the active region.
            obj.set_A_act(gain.device.waveguide.A_act);
            % Set default simulation time for mbsolve
            obj.t = 2 * gain.device.waveguide.l_waveguide ...
                * gain.device.n_eff / phys_const.c0;
        end

        function set_OAB(obj, oab)
            % Sets the Order of Adams-Bashforth for mb (manuel).
            obj.OAB = oab;
        end

        function set_Int(obj, int)
            % Sets seed intensity for mb (manuel).
            obj.Int = int;
        end

        function set_num_reg(obj, numreg)
            % Sets number of simulation regions for mb (manuel).
            obj.num_reg = numreg;
        end

        function set_r_trip(obj, num_trip)
            % Sets the number of round trips for mb (manuel).
            obj.r_trip = num_trip;
            % Change in time-grid points for mb.
            obj.generate_grid();
            % Change simulation time for mbsolve.
            obj.t = 2 * obj.gain.device.waveguide.l_waveguide ...
                * obj.gain.device.n_eff / phys_const.c0 * obj.r_trip;
        end

        function set_N_z_min(obj, z_min)
            % Sets the starting grid point for mb (manuel).
            obj.N_z_min = z_min;
        end

        function set_dN_z(obj, delta)
            % Sets the grid point steps for mb (manuel).
            obj.dN_z = delta;
        end

        function set_dN_t(obj, delta)
            % Sets the time-grid point steps for mb (manuel).
            obj.dN_t = delta;
        end

        function set_shape_cavity(obj, type)
            % Sets a specific shape for mb (manuel).
            obj.bc = obj.shape(type);
            % Change in time-grid.
            obj.generate_grid();
        end

        function set_dens_carrier(obj, n)
            % Sets the space-averaged electron density for mb (manual).
            obj.dens_carrier = n;
        end

        function set_GVD(obj, gvd)
            % Sets group velocity dispersion for mb (manual).
            obj.GVD = gvd;
        end

        function set_rel_permeab(obj, mu)
            % Sets dispersion for mbsolve (manual).
            obj.rel_permeab = mu;
        end

        function set_t_sim(obj, t)
            % Sets simulation time for mbsolve (manual).
            obj.t = t;
        end

        function set_i_wf(obj, wf)
            % Sets period of interest of the reduced system (manual).
            obj.i_wf = wf;
        end

        function set_num_wfs(obj, num)
            % Sets number of wavefunctions (manual).
            obj.num_wfs = num;
        end

        function set_l_wg(obj, l)
            % Sets waveguide length (manual).
            obj.l_wg = l;
        end

        function set_rel_perm(obj, perm)
            % Sets relative permittivity for mbsolve (manual).
            obj.rel_perm = perm;
        end

        function set_n_eff(obj, neff)
            % Sets refractive index (manual).
            obj.n_eff = neff;
        end

        function set_a_power(obj, a)
            % Sets power absorption coefficient (manual).
            obj.a_power = a;
        end

        function set_overlap_factor(obj, Gamma)
            % Sets overlap factor (manual).
            obj.overlap_factor = Gamma;
        end

        function set_D(obj, D)
            % Sets diffusion coefficient for mb (manual).
            obj.D = D;
        end

        function set_facet_refl_coeff(obj, r10, r12)
            % Sets reflection/ transmission coefficient.
            t_10 = sqrt(1-r10^2);
            t_12 = sqrt(1-r12^2);
            if isreal(r10)
                obj.r10 = r10;
            else
                % Consider imaginary values.
                obj.r10 = [real(r10), imag(r10)];
            end
            if isreal(t_10)
                obj.t10 = t_10;
            else
                % Consider imaginary values.
                obj.t10 = [real(t_10), imag(t_10)];
            end
            if isreal(r12)
                obj.r12 = r12;
            else
                % Consider imaginary values.
                obj.r12 = [real(r12), imag(r12)];
            end
            if isreal(t_12)
                obj.t12 = t_12;
            else
                % Consider imaginary values.
                obj.t12 = [real(t_12), imag(t_12)];
            end
        end

        function set_l_period(obj, l)
            % Sets period lenght for mb (manual).
            obj.l_period = l;
        end

        function set_A_act(obj, A)
            % Sets a cross-sectional area of the active region for
            % spontaneous emission for mb (manuel).
            if A == 0
                obj.A_act = A;
            elseif A > 0
                % Automatically change the seed intensity when spontaneous
                % emission is set.
                obj.A_act = A;
                obj.Int = 0;
            end
        end

        function set_Nz(obj, Nz)
            % Sets the grid points (manual).
            obj.N_z = Nz;
            obj.generate_grid();
        end

        function set_freq_range(obj, fmin, fmax)
            % Sets minimum and maximum frequency.
            obj.freq_range = [fmin, fmax];
        end

        function set_input_data(obj)
            % Sets object of class input_data as property.
            obj.mb_input_data = input_data(obj.name, obj.gain.device, ...
                obj.i_wf, obj.f_c, obj.gain.eigen, obj.gain.c_dist, ...
                obj.gain.dep, obj.sc);
        end

        function set_rates_2lvl(obj, fmin, fmax, ...
                I, t_recovery, r_pure_deph)
            % Calculates scattering rates and dephasing rates for a two
            % level representation of the full system by fitting a
            % Lorentzian lineshape function to the total gain spectrum.
            %
            % Syntax:
            %   set_rates_2lvl(obj)
            %   set_rates_2lvl(obj, fmin, fmax)
            %   set_rates_2lvl(obj, fmin, fmax, I)
            %   set_rates_2lvl(obj, fmin, fmax, I, t_recovery)
            %   set_rates_2lvl(obj, fmin, fmax, I, t_recovery, r_pure_deph)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of the gain spectrum.
            %   fmax (scalar): Maximum frequency of the gain spectrum.
            %   I ([] (default) | scalar): Optical intensity for which the
            %     gain recovery time is calculated [W/m^2].
            %   t_recovery ([] (default) | scalar): Gain recovery time
            %     used for initializing fit-parameter for the dephasing rate.
            %   r_pure_deph ([] (default) | scalar): Pure dephasing rate
            %     used for initializing fit-parameter for the dephasing rate.

            if nargin < 6
                r_pure_deph = [];
                if nargin < 5
                    t_recovery = [];
                    if nargin < 4
                        I = [];
                        if nargin < 3
                            fmin = [];
                            fmax = [];
                        end
                    end
                end
            end
            % Default frequency range.
            if isempty(fmin) || isempty(fmax)
                fmin = obj.freq_range(1);
                fmax = obj.freq_range(2);
            end

            % Calc dephasing and transition rates.
            [scat, deph] = obj.mb_input_data.calc_rates_two_level( ...
                fmin, fmax, I, t_recovery, r_pure_deph);
            obj.scat_rates_2lvl = scat;
            obj.deph_rate_2lvl = deph;
        end


        function [ind_wf, varargout] = find_i_wf(obj, ...
                num_injectors, conduction_band, cut_off)
            % Determines the wavefunction indices of the reduced system by
            % finding the injector states. The injector states are
            % determined based on their contribution to the current
            % density. If not provided as input, the number of injectors is
            % determined automatically based on the cut_off value.
            %
            % Syntax:
            %   [ind_wf] = find_i_wf(obj)
            %   [ind_wf, num_injectors] = find_i_wf(obj)
            %   [ind_wf, num_injectors] = find_i_wf(__, name, value)
            %
            % Name Value Arguments:
            %   num_injectors (scalar): Fixed number of injector states.
            %   conduction_band (conduction_band-object): Contains
            %     conduction band profile of the QCL device. If provided it
            %     is provided as input, the wavefunction and conduction
            %     band profile of the reduced system are automatically
            %     plotted for inspection.
            %   cut_off (0.1 (default) | scalar): Specifies the minimum,
            %     relative contribution of a state to the current density
            %     in order to be classified as an injector state.
            %
            % Output Arguments:
            %   ind_wf (scalar):
            %   num_injectors (scalar): Number of automatically detected
            %     injectors.

            arguments
                obj
                num_injectors = [] % number of injectors
                conduction_band = [] % conduction band object
                cut_off = 0.1 % cutoff for automatic injector detection
            end

            ind_wf = 1:obj.num_wfs;
            % Determine relative contribution of each state to the
            % current density
            r_left = obj.sc.get_scattering_matrix( ...
                "left", obj.num_wfs+ind_wf);
            r_right = obj.sc.get_scattering_matrix( ...
                "right", obj.num_wfs+ind_wf);
            occ = obj.gain.c_dist.get_occupation(obj.num_wfs+ind_wf);
            occ = reshape(occ, [], 1) * obj.gain.device.dens_sheet;
            Ji = sum(r_left.*occ-r_right.*occ, 2) * phys_const.e0 / 1000;
            J_rel = Ji / max(Ji);
            % Find injector rates. If the number of injector levels is not
            % explicitly specified by the user, the number of injectors is
            % determined automatically based on the respective current
            % contributions.
            if ~isempty(num_injectors)
                [~, sortedInds] = sort(J_rel, 'descend');
                inj_ind = sortedInds(1:num_injectors);
            else
                arrayfun(@(ind)fprintf('='), 1:30);
                fprintf('\n')
                fprintf('%5s%4s%-21s\n', 'State', ' ', 'Current contribution')
                arrayfun(@(ind)fprintf('='), 1:30);
                fprintf('\n')
                for i = 1:obj.num_wfs
                    fprintf('%5s%4s%8.2g kV/cm (%.0f%%)\n', ...
                        num2str(i+obj.num_wfs), ' ', Ji(i), J_rel(i)*100);
                end
                arrayfun(@(ind)fprintf('='), 1:30);
                fprintf('\n')
                inj_ind = ind_wf(J_rel >= cut_off);
                num_injectors = length(inj_ind);
                fprintf('Number of injectors detected: %d\n\n', num_injectors);
            end
            varargout{1} = num_injectors;
            % create array with indices for one period and replace the
            % injector states with indices from the next period.
            ind_inj_left = inj_ind + obj.num_wfs;
            for i = 1:length(inj_ind)
                ind_wf(inj_ind(i) == ind_wf) = ind_inj_left(i);
            end
            % plot wavefunctions of reduced system
            if ~isempty(conduction_band)
                figure
                ax = subplot(1, 1, 1);
                hold on
                plot(ax, -conduction_band.zv/10, conduction_band.Vh/phys_const.e0, ...
                    "Color", [0, 0, 0], "LineWidth", 1, ...
                    "DisplayName", "CB")
                psis2 = obj.gain.eigen.get_psi_squared();
                for i = 1:2 * obj.num_wfs
                    if any(i == ind_wf)
                        if any(i == ind_inj_left)
                            plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                                psis2(:, i), "Color", ...
                                [162, 173, 0]./255, "LineWidth", 2, ...
                                "DisplayName", [num2str(i), ' (injector)'])
                        else
                            plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                                psis2(:, i), "Color", ...
                                [0, 101, 189]./255, "LineWidth", 2, ...
                                "DisplayName", num2str(i))
                        end
                    else
                        plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                            psis2(:, i), "Color", ...
                            [217, 218, 219]./255, "LineWidth", 1, ...
                            "DisplayName", num2str(i))
                    end
                end
                legend();
                xend = (-obj.gain.eigen.z_wf(1) - ...
                    obj.gain.device.l_period) / 10;
                xstart = (-obj.gain.eigen.z_wf(1) - ...
                    3 * obj.gain.device.l_period - ...
                    obj.gain.device.layers{1}.length) / 10;
                xlim([xstart, xend]) % units: nm
            end
        end

        function plot_wavefunctions(obj, cond_band)
            % Plots the wavefunctions of the reduced system, highlighting
            % upper/lower laser levels as well as injector states.
            %
            % Syntax:
            %   plot_wavefunctions(obj, cond_band)
            %
            % Input Arguments:
            %   cond_band (conduction_band-object): Biased conduction band
            %     profile.

            figure
            ax = subplot(1, 1, 1);
            hold on
            plot(ax, -cond_band.zv/10, cond_band.Vh/phys_const.e0, ...
                "Color", [0, 0, 0], "LineWidth", 1, ...
                "DisplayName", "CB")
            psis2 = obj.gain.eigen.get_psi_squared();
            ind_wfs = obj.mb_input_data.ind_wfs;

            % Extract upper and lower laser states
            ul = [];
            ll = [];
            for i = 1:length(obj.mb_input_data.pairs_dipole)
                ul = [ul, obj.mb_input_data.pairs_dipole{1, i}(1)];
                ll = [ll, obj.mb_input_data.pairs_dipole{1, i}(2)];
            end
            ul = unique(ul);
            ll = unique(ll);

            % Extract injectors
            periods = ceil(ind_wfs/length(ind_wfs));
            periods = periods - (min(periods) - 1);
            inj = ind_wfs(periods == 2);

            % Plot wavefunctions
            for i = obj.num_wfs + 1:3 * obj.num_wfs
                if any(i == ind_wfs)
                    if any(i == ul)
                        plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                            psis2(:, i), "Color", ...
                            [196, 7, 27]./255, "LineWidth", 2, ...
                            "DisplayName", [num2str(i), ' (upper LL)'])
                    elseif any(i == ll)
                        plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                            psis2(:, i), "Color", ...
                            [227, 114, 34]./255, "LineWidth", 2, ...
                            "DisplayName", [num2str(i), ' (lower LL)'])
                    elseif any(i == inj)
                        plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                            psis2(:, i), "Color", ...
                            [162, 173, 0]./255, "LineWidth", 2, ...
                            "DisplayName", [num2str(i), ' (injector)'])
                    else
                        plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                            psis2(:, i), "Color", ...
                            [0, 101, 189]./255, "LineWidth", 2, ...
                            "DisplayName", num2str(i))
                    end
                else
                    plot(ax, -(obj.gain.eigen.z_wf / 10), ...
                        psis2(:, i), "Color", ...
                        [217, 218, 219]./255, "LineWidth", 1, ...
                        "DisplayName", num2str(i))
                end
            end
            legend()
            xend = (-obj.gain.eigen.z_wf(1) - ...
                2 * obj.gain.device.l_period) / 10;
            xstart = (-obj.gain.eigen.z_wf(1) - ...
                4 * obj.gain.device.l_period - ...
                obj.gain.device.layers{1}.length) / 10;
            xlim([xstart, xend]) % units: nm
        end

        function generate_mbsolve(obj)
            % Generates the python input file for the mbsolve simulation
            % tool.
            %
            % Syntax:
            %   generate_mbsolve(obj)

            if isempty(obj.mb_input_data)
                obj.set_input_data();
            end
            % Open and generate mbsolve python script.
            name_py = append(obj.name, '.py');
            obj.mkdir_if_not_exist(obj.dir);
            fileID = fopen(fullfile(obj.dir, name_py), 'w');
            fprintf(fileID, ['import mbsolve.lib as mb\n', ...
                'import mbsolve.solvercpu\n', ...
                'import mbsolve.writerhdf5\n', ...
                '\nimport math\nimport time\n\n']);

            % Write Hamiltonian description.
            fprintf(fileID, '# Hamiltonian\n');
            energies = diag(obj.mb_input_data.hamiltonian/phys_const.e0);
            string_e = obj.generate_string(energies, 'energies = [ ', ...
                ' * mb.E0', '%.4f');
            fprintf(fileID, append(string_e, ' ]\n'));
            ind = triu(true(size(obj.mb_input_data.hamiltonian)), 1);
            off_diagonales = obj.mb_input_data.hamiltonian(ind) / ...
                phys_const.e0;
            string_off = obj.generate_string(off_diagonales, ...
                'off_diagonales = [ ', ' * mb.E0', '%.4f');
            fprintf(fileID, append(string_off, ' ]\n'));
            fprintf(fileID, ...
                'H = mb.qm_operator(energies, off_diagonales)\n');
            fprintf(fileID, '\n');

            % Write dipole moment operator description.
            fprintf(fileID, '# dipole moment operator\n');
            ind = triu(true(size(obj.mb_input_data.dipole_matrix)), 1);
            dipoles = obj.mb_input_data.dipole_matrix(ind) / phys_const.e0;
            string_dipoles = obj.generate_string(dipoles, ...
                'off_dipoles = [ ', ' * mb.E0', '%.4e');
            fprintf(fileID, append(string_dipoles, ' ]\n'));
            string_main_dipoles = ...
                obj.generate_string(zeros( ...
                size(obj.mb_input_data.dipole_matrix, 1), 1), ... .
            'diag_dipoles = [ ', ' * mb.E0', '%.4e');
            fprintf(fileID, append(string_main_dipoles, ' ]\n'));
            fprintf(fileID, ...
                'u = mb.qm_operator(diag_dipoles, off_dipoles)\n');
            fprintf(fileID, '\n');

            % Write relaxation superoperator description.
            fprintf(fileID, '# relaxation superoperator\n');
            fprintf(fileID, '# scattering rate matrix R\n');
            scat = obj.mb_input_data.transition_rates;
            scat = scat - diag(diag(scat));
            for i = 1:length(scat) - 1
                if (i == 1)
                    string_scat = obj.generate_string(scat(:, i), ...
                        'rates = [ [ ', '', '%.4e');
                else
                    string_scat = obj.generate_string(scat(:, i), ...
                        '          [ ', '', '%.4e');
                end
                fprintf(fileID, append(string_scat, ' ],\n'));
            end
            string_scat = obj.generate_string(scat(:, length(scat)), ...
                '          [ ', '', '%.4e');
            fprintf(fileID, append(string_scat, ' ] ]\n'));
            fprintf(fileID, '\n');

            % Write pure dephasing description.
            fprintf(fileID, '# pure dephasing rates\n');
            % Column-major ordered pure dephasing rates vector.
            ind = triu(true(size( ...
                obj.mb_input_data.pure_dephasing_rates)), 1);
            pure_deph = obj.mb_input_data.pure_dephasing_rates(ind);
            string_pure_deph = obj.generate_string(pure_deph, ...
                'pure_deph = [ ', '', '%.4e');
            fprintf(fileID, append(string_pure_deph, ' ]\n'));
            fprintf(fileID, ...
                append('relax_sop', ...
                ' = mb.qm_lindblad_relaxation(rates, pure_deph)\n'));
            fprintf(fileID, '\n');

            % Initialize density matrix diagonal elements (occupations).
            fprintf(fileID, '# initial density matrix \n');
            occ = zeros(length(obj.mb_input_data.ind_wfs));
            for i = 1:length(obj.mb_input_data.ind_wfs)
                occ(i) = obj.mb_input_data.carr_dist.get_occupation( ...
                    obj.mb_input_data.ind_wfs(i));
            end
            string_occ = obj.generate_string(occ, ...
                'rho_init = mb.qm_operator([ ', '', '%.4f');
            fprintf(fileID, append(string_occ, '])\n'));
            fprintf(fileID, '\n');

            % Write quantum mechanical description.
            fprintf(fileID, '# quantum mechanical description\n');
            % Length of one period in m
            doping_dens = obj.gain.device.dens_carrier();
            string_qm = append('qm = mb.qm_description(', ...
                num2str(doping_dens, '%.4e'), ', H, u, relax_sop)');
            fprintf(fileID, append(string_qm, '\n'));
            fprintf(fileID, ['mat_ar = mb.material("AR_', ...
                obj.name(1:5), '", qm, ', ...
                num2str(obj.gain.device.rel_permittivity), ', ', ...
                num2str(obj.gain.device.waveguide.overlap_factor), ', ', ...
                num2str(obj.gain.device.waveguide.a_field), ', ', ...
                num2str(1), ')\n']);
            fprintf(fileID, 'mb.material.add_to_library(mat_ar)\n');
            fprintf(fileID, '\n');
            fprintf(fileID, ['dev = mb.device("', ...
                num2str(length(obj.mb_input_data.ind_wfs)), 'lvl")\n']);
            fprintf(fileID, ['dev.add_region(mb.region("Active region"', ...
                ', mat_ar, 0.0, ', ...
                num2str(obj.gain.device.waveguide.l_waveguide), '))\n']);
            fprintf(fileID, '\n');
            fprintf(fileID, '# Scenario\n');
            fprintf(fileID, 'ic_d = mb.ic_density_const(rho_init)\n');
            fprintf(fileID, 'ic_e = mb.ic_field_const(0.0)\n');
            dn_z = phys_const.c0 / ...
                (obj.n_eff * obj.f_c) / 20; % Grid spacing
            fprintf(fileID, ['sce = mb.scenario("', obj.name, ...
                '", ', num2str(floor(obj.l_wg/dn_z)), ', ', ...
                num2str(obj.t), ', ', 'rho_init)\n\n']);

            % Write ouput options for electric and magnetic field
            fprintf(fileID, ['sampling_interval = ', ...
                num2str(0.1*obj.t), '\n']);
            fprintf(fileID, ['sce.add_record(mb.record("e", ', ...
                'sampling_interval, -1))\n']);
            fprintf(fileID, ['sce.add_record(mb.record("h", ', ...
                'sampling_interval, -1))\n']);
            fprintf(fileID, ['sce.add_record(mb.record("e1", ', ...
                '0, ', num2str(obj.l_wg), '))\n']);
            fprintf(fileID, ['sce.add_record(mb.record("h1", ', ...
                '0, ', num2str(obj.l_wg), '))\n']);
            % Write ouput options for density matrix elements
            for i = 1:length(obj.mb_input_data.ind_wfs)
                fprintf(fileID, ['sce.add_record(mb.record("d', ...
                    num2str(i), num2str(i), '", mb.record.density, ', ...
                    num2str(i), ', ', num2str(i), ', ', ...
                    'sampling_interval, -1))\n']);
            end
            fprintf(fileID, '\n');

            fprintf(fileID, '# run solver\n');
            fprintf(fileID, ['sol = mb.solver.create_instance(', ...
                '"cpu-fdtd-', num2str(length(obj.mb_input_data.ind_wfs)), ...
                'lvl-reg-cayley", dev, sce)\n']);
            fprintf(fileID, ['print(', '''Solver ''', ...
                ' + sol.get_name() + ', ''' started''', ')\n']);
            fprintf(fileID, 'tic = time.time()\n');
            fprintf(fileID, 'sol.run()\n');
            fprintf(fileID, 'toc = time.time()\n');
            fprintf(fileID, ['print(', '''Solver ''', ...
                ' + sol.get_name() + ', ''' finished in ''', ...
                ' + str(toc - tic) + ', ''' sec'')\n\n']);

            fprintf(fileID, '# write results\n');
            fprintf(fileID, 'wri = mb.writer.create_instance("hdf5")\n');
            fprintf(fileID, ['outfile = dev.get_name() + "_" + ', ...
                'sce.get_name() + "." + wri.get_extension()\n']);
            fprintf(fileID, 'results = sol.get_results()\n');
            fprintf(fileID, ['wri.write(outfile, ', ...
                'sol.get_results(), dev, sce)\n']);
            fclose(fileID);
        end

        function generate_mbsolve_2lvl(obj)
            % Generates the python input file for the mbsolve simulation
            % tool for a two-level representation of the full system.
            %
            % Syntax:
            %   generate_mbsolve_2lvl(obj)

            if isempty(obj.mb_input_data)
                obj.set_input_data();
            end

            % Determine indices of two level system
            i_ind = obj.mb_input_data.pairs_dipole{1, 1}(1); % full system
            f_ind = obj.mb_input_data.pairs_dipole{1, 1}(2); % full system
            ni = find(obj.mb_input_data.ind_wfs == i_ind); % reduced system
            nf = find(obj.mb_input_data.ind_wfs == f_ind); % reduced system

            % Open and generate mbsolve python script.
            name_py = append(obj.name, '_2lvl.py');
            obj.mkdir_if_not_exist(obj.dir);
            fileID = fopen(fullfile(obj.dir, name_py), 'w');
            fprintf(fileID, ['import mbsolve.lib as mb\n', ...
                'import mbsolve.solvercpu\n', ...
                'import mbsolve.writerhdf5\n', ...
                '\nimport math\nimport time\n\n']);

            % Write Hamiltonian description.
            fprintf(fileID, '# Hamiltonian\n');
            energies = [obj.mb_input_data.get_trans_freq(i_ind, f_ind) ...
                * phys_const.h / phys_const.e0, 0];
            string_e = obj.generate_string(energies, 'energies = [ ', ...
                ' * mb.E0', '%.4f');
            fprintf(fileID, append(string_e, ' ]\n'));
            off_diagonales = 0;
            string_off = obj.generate_string(off_diagonales, ...
                'off_diagonales = [ ', ' * mb.E0', '%.4f');
            fprintf(fileID, append(string_off, ' ]\n'));
            fprintf(fileID, ...
                'H = mb.qm_operator(energies, off_diagonales)\n');
            fprintf(fileID, '\n');

            % Write dipole moment operator description.
            fprintf(fileID, '# dipole moment operator\n');
            dipoles = obj.mb_input_data.dipole_matrix(ni, nf) ./ ...
                phys_const.e0;
            string_dipoles = obj.generate_string(dipoles, ...
                'off_dipoles = [ ', ' * mb.E0', '%.4e');
            fprintf(fileID, append(string_dipoles, ' ]\n'));
            string_main_dipoles = ...
                obj.generate_string(zeros(2, 1), ...
                'diag_dipoles = [ ', ' * mb.E0', '%.4e');
            fprintf(fileID, append(string_main_dipoles, ' ]\n'));
            fprintf(fileID, ...
                'u = mb.qm_operator(diag_dipoles, off_dipoles)\n');
            fprintf(fileID, '\n');

            % Write relaxation superoperator description.
            fprintf(fileID, '# relaxation superoperator\n');
            fprintf(fileID, '# scattering rate matrix R\n');
            string_scat_row_1 = obj.generate_string( ...
                obj.scat_rates_2lvl(:, 1), 'rates = [ [ ', '', '%.4e');
            fprintf(fileID, append(string_scat_row_1, ' ],\n'));
            string_scat_row_2 = obj.generate_string( ...
                obj.scat_rates_2lvl(:, 2), '          [ ', '', '%.4e');
            fprintf(fileID, append(string_scat_row_2, ' ] ]\n'));
            fprintf(fileID, '\n');

            % Write pure dephasing description.
            fprintf(fileID, '# pure dephasing rates\n');
            pure_deph = obj.deph_rate_2lvl - ...
                0.5 * sum(sum(obj.scat_rates_2lvl));
            string_pure_deph = obj.generate_string(pure_deph, ...
                'pure_deph = [ ', '', '%.4e');
            fprintf(fileID, append(string_pure_deph, ' ]\n'));
            fprintf(fileID, ...
                append('relax_sop', ...
                ' = mb.qm_lindblad_relaxation(rates, pure_deph)\n'));
            fprintf(fileID, '\n');

            % Initialize density matrix diagonal elements (occupations).
            fprintf(fileID, '# initial density matrix \n');
            rate21 = obj.scat_rates_2lvl(2, 1);
            rate12 = obj.scat_rates_2lvl(1, 2);
            occ = [rate21, rate12] / (rate21 + rate12);
            string_occ = obj.generate_string(occ, ...
                'rho_init = mb.qm_operator([ ', '', '%.4f');
            fprintf(fileID, append(string_occ, '])\n'));
            fprintf(fileID, '\n');

            % Write quantum mechanical description.
            fprintf(fileID, '# quantum mechanical description\n');
            doping_dens = obj.gain.device.dens_carrier();
            string_qm = append('qm = mb.qm_description(', ...
                num2str(doping_dens, '%.4e'), ', H, u, relax_sop)');
            fprintf(fileID, append(string_qm, '\n'));
            fprintf(fileID, ['mat_ar = mb.material("AR_', ...
                obj.name(1:5), '", qm, ', ...
                num2str(obj.gain.device.rel_permittivity), ', ', ...
                num2str(obj.gain.device.waveguide.overlap_factor), ', ', ...
                num2str(obj.gain.device.waveguide.a_field), ', ', ...
                num2str(1), ')\n']);
            fprintf(fileID, 'mb.material.add_to_library(mat_ar)\n');
            fprintf(fileID, '\n');
            fprintf(fileID, ['dev = mb.device("', ...
                num2str(2), 'lvl")\n']);
            fprintf(fileID, ['dev.add_region(mb.region("Active region"', ...
                ', mat_ar, 0.0, ', ...
                num2str(obj.gain.device.waveguide.l_waveguide), '))\n']);
            fprintf(fileID, '\n');
            fprintf(fileID, '# Scenario\n');
            fprintf(fileID, 'ic_d = mb.ic_density_const(rho_init)\n');
            fprintf(fileID, 'ic_e = mb.ic_field_const(0.0)\n');
            dn_z = phys_const.c0 / ...
                (obj.n_eff * obj.f_c) / 20; % Grid spacing
            fprintf(fileID, ['sce = mb.scenario("', obj.name, ...
                '", ', num2str(floor(obj.l_wg/dn_z)), ', ', ...
                num2str(obj.t), ', ', 'rho_init)\n\n']);

            % Write output options for electric and magnetic field
            fprintf(fileID, ['sampling_interval = ', ...
                num2str(0.1*obj.t), '\n']);
            fprintf(fileID, ['sce.add_record(mb.record("e", ', ...
                'sampling_interval, -1))\n']);
            fprintf(fileID, ['sce.add_record(mb.record("h", ', ...
                'sampling_interval, -1))\n']);
            fprintf(fileID, ['sce.add_record(mb.record("e1", ', ...
                '0, ', num2str(obj.l_wg), '))\n']);
            fprintf(fileID, ['sce.add_record(mb.record("h1", ', ...
                '0, ', num2str(obj.l_wg), '))\n']);
            % Write ouput options for density matrix elements
            for i = 1:2
                fprintf(fileID, ['sce.add_record(mb.record("d', ...
                    num2str(i), num2str(i), '", mb.record.density, ', ...
                    num2str(i), ', ', num2str(i), ', ', ...
                    'sampling_interval, -1))\n']);
            end
            fprintf(fileID, '\n');

            fprintf(fileID, '# run solver\n');
            fprintf(fileID, ['sol = mb.solver.create_instance(', ...
                '"cpu-fdtd-', num2str(2), 'lvl-reg-cayley", ', ...
                'dev, sce)\n']);
            fprintf(fileID, ['print(', '''Solver ''', ...
                ' + sol.get_name() + ', ''' started''', ')\n']);
            fprintf(fileID, 'tic = time.time()\n');
            fprintf(fileID, 'sol.run()\n');
            fprintf(fileID, 'toc = time.time()\n');
            fprintf(fileID, ['print(', '''Solver ''', ...
                ' + sol.get_name() + ', ''' finished in ''', ...
                ' + str(toc - tic) + ', ''' sec'')\n\n']);

            fprintf(fileID, '# write results\n');
            fprintf(fileID, 'wri = mb.writer.create_instance("hdf5")\n');
            fprintf(fileID, ['outfile = dev.get_name() + "_" + ', ...
                'sce.get_name() + "." + wri.get_extension()\n']);
            fprintf(fileID, 'results = sol.get_results()\n');
            fprintf(fileID, ['wri.write(outfile, ', ...
                'sol.get_results(), dev, sce)\n']);
            fclose(fileID);
        end
    end

    methods (Access = private)
        function calc_fc(obj, fmin, fmax, id)
            % Calculates and sets the center frequency of the laser
            % emission.
            %
            % Syntax:
            %   calc_fc(obj, fmin, fmax, id)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum considered frequency [Hz].
            %   fmax (scalar): Maximum considered frequency [Hz].
            %   id (char): Broadening mechanism, either pure dephasing +
            %     lifetime broadening (``id=total``) or only lifetime
            %     broadening (``id=ltbroad``).

            % Frequency vector between fmin and fmax of length 100.
            f = (fmin + (fmax - fmin) * (1:100) / 100);
            % Gain vector between fmin and fmax of length 100.
            g = obj.gain.get_gain(fmin, fmax, id);
            % Get index of max gain.
            [~, imax] = max(g);
            % Calculate the FWHM in the given frequency range.
            FWHM = full_width_half_max.calc(g) * (f(2) - f(1)) / 1e12;
            % Check for gain peak.
            if isnan(FWHM)
                disp('No gain in this frequency spectrum.');
                return
            else
                % Set center frequency at max gain index.
                obj.f_c = f(imax);
            end
            % Calculate the grid points.
            obj.generate_grid();
        end

        function generate_grid(obj)
            % Sets the properties for the spatial and temporal grid for the
            % mb-tool based on the center frequency, the cavity shape and
            % the number of round trips.
            %
            % Syntax:
            %   generate_grid(obj)

            if isempty(obj.N_z)
                % Default grid point number
                obj.N_z = 1000;
            end
            % Max. grid point for the E-field plot.
            obj.N_z_max = obj.N_z - 1;
            if obj.bc == 0
                % Time-grid points for one round trip * number of round
                % trips (Fabry-Perot).
                obj.N_t = obj.N_z * 2 * obj.r_trip;
                % Min. time-grid point for the E-field plot. Only consider
                % last round trip.
                obj.N_t_min = obj.N_t - obj.N_z * 2;
            else
                % Time-grid points for one round trip * number of round
                % trips (ring).
                obj.N_t = obj.N_z * obj.r_trip;
                % Min. time-grid point for the E-field plot. Only consider
                % last round trip.
                obj.N_t_min = obj.N_t - obj.N_z;
            end
        end
    end

    methods (Static, Access = private)
        function str_data = generate_string(data, str_init, str_add, formatSpec)
            % Generates string from input data with predefined line length.
            %
            % Syntax:
            %    str_data = generate_string(data, str_init, str_add, formatSpec)
            %
            % Input Arguments:
            %   data (vector | scalar): Vector
            %   str_init (char): String which is printed before the data
            %     array.
            %   str_add (char): String added after each element of the data
            %     array.
            %   formatSpec (char): Format specifiers.
            %
            % Output Arguments:
            %   str_data (char): Formatted string of the data.

            str_data = str_init;
            l_line = length(str_data);
            for i = 1:length(data)
                if (data(i) == 0)
                    str_i = '0';
                else
                    str_i = append(num2str(data(i), formatSpec), str_add);
                end
                if i < length(data)
                    str_i = append(str_i, ', ');
                end
                l_line = l_line + length(str_i);
                if (l_line > 78)
                    str_data = append(str_data, '\n');
                    l_line = length(str_i);
                end
                str_data = append(str_data, str_i);
            end
        end

        function mkdir_if_not_exist(directory)
            % Generates directory if it is not existing.
            %
            % Syntax:
            %   mkdir_if_not_exist(directory)
            %
            % Input Arguments:
            %   directory (char | string): Absolute or relative path to the
            %     directory to be generated.

            if ~exist(directory, 'dir')
                mkdir(directory)
            end
        end
    end
end
