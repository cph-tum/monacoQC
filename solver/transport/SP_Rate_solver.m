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

classdef SP_Rate_solver < Stationary_Carrier_Transport
    % Solves carrier transport self-consistently based on rate-equations.

    properties
        options_rates % struct: Options for rate equation solver.
        scat_mechanisms % cell-array: Array with names of scattering mechanisms.
    end

    methods
        function obj = SP_Rate_solver(d, s, mechanisms)
            % Constructs an object of type SP_Rate_solver.
            %
            % Syntax:
            %   obj = SP_Rate_solver(d, s)
            %   obj = SP_Rate_solver(d, s, mechanisms)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %   mechanisms (char cell-array): Array containing names of the
            %     scattering mechanisms to be taken into account in the
            %     simulation. Valid names are ``lo phonon absorption``,
            %     ``lo phonon emission``, ``alloy disorder``, ``impurity``,
            %     ``acoustic phonon``, ``interface roughness``, ``electron
            %     electron``, ``tunneling`` and ``optical field``. By
            %     default all mechanisms are included except for
            %     ``optical field``.

            obj@Stationary_Carrier_Transport(d, s);

            if nargin == 3
                obj.scat_mechanisms = mechanisms;
            else
                % Set default scattering mechanisms
                if obj.scenario.basis_sp == "tb"
                    % For tight-binding states include tunneling rates
                    obj.scat_mechanisms = {'lo phonon absorption', ...
                        'lo phonon emission', 'alloy disorder', ...
                        'acoustic phonon', 'impurity', ...
                        'interface roughness', 'electron electron', ...
                        'tunneling'};
                else
                    obj.scat_mechanisms = {'lo phonon absorption', ...
                        'lo phonon emission', 'alloy disorder', ...
                        'acoustic phonon', 'impurity', ...
                        'interface roughness', 'electron electron'};
                end
            end

            % Default options for rate solver
            Te = obj.scenario.T * ones(1, obj.scenario.num_wavefct*4);
            obj.options_rates = struct( ...
                "num_k", 50, ...
                "num_theta", 30, ...
                "parallel", false, ...
                "num_interp", 200, ...
                "max_iteration", 10, ...
                "relTol", 1e-2, ...
                "T_e", Te, ...
                "screening_model", "static-lindhard");
        end

        function [cond, eigen, curr_dens, dist, deph, sc, gain, ...
                rateSolver, m, varargout] = solve(obj)
            % Simulates carrier transport by solving rate equation and
            % Schrödinger-Poisson equation self-consistently.
            %
            % Syntax:
            %   [cond, eigen, J_new, dist, deph, sc, gain, rateSolver, m] = solve(obj)
            %   [cond, eigen, J_new, dist, deph, sc, gain, rateSolver, m, cavity_field] = solve(obj)
            %
            % Output Arguments:
            %   cond (conduction_band-object): Contains information about
            %     the biased conduction band profile including space charge
            %     effects.
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.
            %   dist (carrier_distribution-object): Contains information
            %     about the subband occupations/ distributions.
            %   deph (dephasing_rates): Contains information about the
            %     dephasing rates (lifetime +  pure dephasing).
            %   sc (scattering_rates-object): Contains information about
            %     the scattering rates.
            %   gain (sim_gain-object): Contains information about the gain
            %     spectrum.
            %   rateSolver (RateSolver-object): Contains information about
            %     the rate equations.
            %   curr_dens (current_density-object): Contains information
            %     about the current density.
            %   m (scalar): Current iteration number when loop is
            %     terminated.
            %   cavity_field (optical_field-object): Contains information
            %     about the cavity field.

            % self consistent loop
            fprintf(['Start solving coupled Schrödinger-Poisson-Rate ', ...
                'equations ...\n\n'])
            total_time = tic;
            for m = 1:obj.max_iteration
                if m == 1
                    % Solve Schrödinger-Poisson equation
                    SP_solver = tm_solver(obj.device, obj.scenario);
                    [eigen, cond] = SP_solver.solve();
                    cond_band = conduction_band(cond.zv, ...
                        SP_solver.sim_const.vec_V_0);
                    % Initialize rate rate solver
                    rates = obj.init_rates(eigen, cond_band);
                    rateSolver = RateSolver(rates);
                else
                    % Solve Schrödinger-Poisson equation
                    [eigen, cond] = SP_solver.solve(@(eig)dist, eigen);
                    % Update rate solver
                    rateSolver.update_eigenstates(eigen);
                end
                % Solve rate equations
                dist = rateSolver.solve( ...
                    "parallel", obj.options_rates.parallel, ...
                    "relTol", obj.options_rates.relTol, ...
                    "nmax", obj.options_rates.max_iteration);

                % Postprocessing
                J_new = rateSolver.calc_current_density();
                deph = rateSolver.get_dephasing_rates(dist);
                gain = sim_gain(obj.device, eigen, deph, dist);
                sc = rateSolver.get_scattering_rates();
                curr_dens = current_density(0, J_new);
                % Check if optical field is included
                if rateSolver.scattering_rates{end}.Name == "optical field" ...
                        && nargout > 9
                    cavity_field = rateSolver.get_optical_field();
                    cavity_field.set_waveguide_properties(obj.device);
                    varargout = {cavity_field};
                end

                if m > 1
                    err = abs((J_new - J_old)/J_old);
                    fprintf(['Schrödinger-Poisson-Rate solver has ', ...
                        'been passed ', num2str(m), ' times. ', ...
                        'Residual: ', num2str(err), '.\n\n'])
                    if err < obj.relTol
                        fprintf(['Schrödinger-Poisson-Rate solver ', ...
                            'reached desired residual of ', num2str( ...
                            obj.relTol), '. Elapsed time: ', num2str( ...
                            toc(total_time)), 's\n\n\n'])
                        break
                    end
                    if m == obj.max_iteration
                        fprintf(['Schrödinger-Poisson-Rate solver did ', ...
                            'not reach desired residual of ', num2str( ...
                            obj.relTol), '. Elapsed time: ', num2str( ...
                            toc(total_time)), 's\n\n'])
                        break
                    end
                end
                fprintf('\n')
                J_old = J_new;
            end
        end

        function obj = set_rate_eq_options(obj, options)
            % Sets options for rate equation solver.
            %
            % Syntax:
            %   obj = set_rate_eq_options(obj, options)
            %
            % Name Value Arguments:
            %   num_k (scalar): Number of k-values.
            %   num_theta (scalar): Number of values for angle theta.
            %   parallel (logical): Flag enabling parallel computation of
            %     scattering rates.
            %   Te (vector): Electron temperatures of each subband in one
            %     period.
            %   num_interp: Number of points for the pre-calculation of the
            %     formfactor.
            %   max_iteration: Maximum number of iterations.
            %   relTol: Desired relative error between two subsequent
            %     iterations to terminate the loop.
            %   screening_model (char): Name of screening model. Valid
            %     names are ``static-lindhard``, ``debye``,
            %     ``thomas-fermi`` and ``modified-single-subband``.

            arguments
                obj
                options.num_k
                options.num_theta
                options.parallel
                options.num_interp
                options.max_iteration
                options.relTol
                options.T_e
                options.screening_model
            end

            keys = fieldnames(options);
            values = struct2cell(options);

            for i = 1:length(keys)
                obj.options_rates.(keys{i}) = values{i};
            end
        end

        function rates = init_rates(obj, eigen, cond_band)
            % Creates cell-array of scattering rate objects.
            %
            % Syntax:
            %   rates = init_rates(obj, eigen, cond_band)
            %
            % Input Arguments:
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses
            %   cond_band (conduction_band-object): Contains information
            %     about the conduction band profile.

            % get options for rate-equations
            rateOptions = {"num_k", obj.options_rates.num_k, ...
                "Te", obj.options_rates.T_e, "screening_model", ...
                obj.options_rates.screening_model};
            eeRateOptions = {"num_k", obj.options_rates.num_k, ...
                "num_theta", obj.options_rates.num_theta, ...
                "Te", obj.options_rates.T_e, ...
                "num_interp", obj.options_rates.num_interp, ...
                "screening_model", obj.options_rates.screening_model};

            % create cell-array with scattering mechanisms
            rates = {length(obj.scat_mechanisms), 1};
            for i = 1:length(obj.scat_mechanisms)
                name = obj.scat_mechanisms{i};
                switch name
                    case 'lo phonon absorption'
                        rates{i} = LOPhononRate(eigen, obj.device, ...
                            obj.scenario, 'absorption', rateOptions);
                    case 'lo phonon emission'
                        rates{i} = LOPhononRate(eigen, obj.device, ...
                            obj.scenario, 'emission', rateOptions);
                    case 'alloy disorder'
                        rates{i} = AlloyRate(eigen, obj.device, ...
                            obj.scenario, rateOptions);
                    case 'acoustic phonon'
                        rates{i} = AcousticPhononRate(eigen, ...
                            obj.device, obj.scenario, rateOptions);
                    case 'impurity'
                        rates{i} = ImpurityRate(eigen, ...
                            obj.device, obj.scenario, rateOptions);
                    case 'interface roughness'
                        rates{i} = InterfaceRoughnessRate(eigen, ...
                            obj.device, obj.scenario, rateOptions);
                    case 'electron electron'
                        rates{i} = ElectronElectronRate(eigen, ...
                            obj.device, obj.scenario, eeRateOptions);
                    case 'tunneling'
                        rates{i} = TunnelingRate(eigen, obj.device, ...
                            obj.scenario, cond_band, rateOptions);
                    case 'optical field'
                        rates{i} = OpticalRate(eigen, obj.device, ...
                            obj.scenario, rateOptions);
                    otherwise
                        error(string(name)+' is not a valid name '+ ...
                            'for a scattering mechanism! Valid names '+ ...
                            'are ''lo phonon absorption'', '+ ...
                            '''lo phonon emission'', '+ ...
                            '''alloy disorder'', '+ ...
                            '''acoustic phonon'', '+ ...
                            '''interface roughness'', '+ ...
                            '''electron electron'', '+ ...
                            '''tunneling'', '+ ...
                            '''optical field''.')
                end
            end
        end
    end
end
