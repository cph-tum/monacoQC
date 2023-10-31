%
% monacoQC: An object-oriented Matlab-based device engineering tool for
% quantum cascade structures.
%
% Copyright (C) 2023, Computational Photonics Group, Technical University of
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

classdef boundary_condition_ext < boundary_condition
    %boundary_condition_ext Boundary condition class for extended
    % eigenstates.
    properties
        
    end
    
    methods
        function obj = boundary_condition_ext(d, s, sim_c)
            % Creates an object with properties describing the simulation
            % boundary conditions for extended state solutions.
            obj = obj@boundary_condition(s.basis_sp);
            obj.delta_psi = 0.005;
            % The increment is negative here, assuring to apply the
            % boundary condition for the root-finding algorithm at the
            % approriate boundary position. Otherwise not all wavefunctions
            % will be found.
            obj.incr_z = -1;
            obj.dE = 1e-12 * phys_const.e0;
            % Set index for interpolation.
            obj.ind_interp = 1 + sim_c.get_vec_index( ...
                (d.num_periods - 1)*d.num_layers_period);
            % Point of lowest potential directly next to thickest barrier.
            obj.ind_E_min = 1 + sim_c.get_vec_index( ...
                ceil(d.num_periods/2)*d.num_layers_period+1);
            % Set start point of scan.
            obj.ind_z_L = 2;
            % Set end point of scan.
            obj.ind_z_R = length(sim_c.vec_z_tm);
            % Set indices for wavefunction bounding box.
            obj.ind_z_L_bound = sim_c.get_vec_index_mid(1);
            obj.ind_z_R_bound = sim_c.get_vec_index_mid( ...
                d.num_periods*d.num_layers_period+1);
            % Define Energy period in J.
            obj.Eperiod = sim_c.E_period;
        end
        
        function indices = get_indices(obj, n_scans)
            indices = obj.ind_z_R:obj.incr_z:obj.ind_z_L;
        end
        function index_end = get_index_bc_end(obj, size_m, n_scans)
            index_end = sub2ind(size_m, 1, obj.ind_z_L);
        end
    end
    methods (Static)
        function bc_start = get_bc_start()
            bc_start = [1; 0];
        end
    end
end