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

classdef boundary_condition_tb < boundary_condition
    % Defines the boundary condition for tight-binding eigenstates.

    methods
        function obj = boundary_condition_tb(d, s, sim_c)
            % Creates an object of type boundary_condition_tb.
            %
            % Syntax:
            %   obj = boundary_condition_tb(d, s, sim_c)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %   sim_c (sim_constants-object): Contains simulation
            %    constants.

            obj = obj@boundary_condition(s.basis_sp);

            obj.incr_z = 1;
            obj.delta_psi = 0.001;
            obj.dE = 1e-12 * phys_const.e0;

            % Gets index number of one period.
            ind_plus_period = sim_c.get_vec_index(d.num_layers_period);
            % Set index for interpolation.
            obj.ind_interp = 1 + (d.num_periods - 1) * ind_plus_period;

            % Shifting by one period, accounting for double values at
            % interfaces.
            middle_structure = floor(d.num_periods/2);
            obj.ind_z_L = 1 + sim_c.get_vec_index(middle_structure ...
                *d.num_layers_period);
            % Setting end point of period.
            obj.ind_z_R = sim_c.get_vec_index((middle_structure + 1) ...
                *d.num_layers_period+1);
            % Set indices for wavefunction bounding box.
            % Left index at middle of first barrier in period.
            obj.ind_z_L_bound = sim_c.get_vec_index_mid( ...
                middle_structure*d.num_layers_period+1);
            % Right index at middle of first barrier of next period.
            obj.ind_z_R_bound = sim_c.get_vec_index_mid( ...
                (middle_structure + 1)*d.num_layers_period+1);
            % Point of lowest potential directly next to widest barrier.
            obj.ind_E_min = 1 + sim_c.get_vec_index(middle_structure ...
                *d.num_layers_period+1);

            % Define Energy period in J.
            obj.Eperiod = 1.2 * max(diff(sim_c.vec_V_0));
        end

        function indices = get_indices(obj, n_scan)
            % Gets index vector of considered localization interval
            % (tight-binding period).
            %
            % Input Arguments:
            %   n_scan (scalar): Scan number.
            %
            % Output Arguments:
            %   indices (vector): Index vector.

            indices = obj.ind_z_L:obj.incr_z:obj.ind_z_R;
        end

        function index_end = get_index_bc_end(obj, size_m, n_scan)
            % Gets index of end position of localization interval.
            %
            % Syntax:
            %   index_end = get_index_bc_end(obj, size_m, n_scan)
            %
            % Input Arguments:
            %   size_m (vector): Vector containing the size of the array to
            %     be indexed.
            %   n_scan (scalar): Scan number.
            %
            % Output Arguments:
            %   index_end (scalar): End index of boundary condition interval.

            index_end = sub2ind(size_m, 2, obj.ind_z_R);
        end
    end

    methods (Static)
        function bc_start = get_bc_start()
            % Gets transfer matrix coefficiences at the start position of
            % the localization interval.

            bc_start = [0; 1];
        end
    end
end
