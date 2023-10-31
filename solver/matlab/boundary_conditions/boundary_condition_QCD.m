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

classdef boundary_condition_QCD < boundary_condition
    %boundary_condition_QCD Boundary condition class for QCD eigenstates.
    
    properties
        
    end
    methods
        function obj = boundary_condition_QCD(d, s, sim_c)
            % Creates an object with properties describing the simulation
            % boundary conditions for extended QCD state solutions.
            obj = obj@boundary_condition(s.basis_sp);
            obj.incr_z = 1;
            obj.delta_psi = 0.02;
            obj.dE = 1e-16 * phys_const.e0;
            % Set index for interpolation.
            obj.ind_interp = 1 + sim_c.get_vec_index( ...
                (d.num_periods - 1)*d.num_layers_period);
            % Set number of scans.
            obj.num_scans = 2;
            % Set number of layers overlaping the period end.
            num_overlap_layer = 1;
            ind_barr = d.get_ind_barrier();
            % Get indices of barrier layers in first period.
            i_b = find(ind_barr < d.num_layers_period);
            % Find reference index of a barrier
            % in the middle of one period.
            i_b_ref = i_b(ceil(length(i_b)/2));
            ind_b_ref = ind_barr(i_b_ref);
            % Get indices for left scan.
            % Left index at end barrier minus overlap.
            ind_b1_L = ind_barr(i_b(end-num_overlap_layer));
            % Right index at reference barrier plus overlap.
            ind_b1_R = ind_barr(i_b_ref+num_overlap_layer);
            % Get indices for right scan.
            % Left index at ref barrier minus overlap.
            ind_b2_L = ind_barr(i_b_ref-num_overlap_layer);
            % Right index at first barrier plus overlap.
            ind_b2_R = ind_barr(i_b(1+num_overlap_layer));
            % Set start and end points of scan.
            % Left indices at beginning of barrier.
            obj.ind_z_L(1) = 1 + sim_c.get_vec_index(ind_b1_L-1);
            obj.ind_z_L(2) = 1 + sim_c.get_vec_index( ...
                d.num_layers_period+ind_b2_L-1);
            % Right indices at ending of barrier.
            obj.ind_z_R(1) = sim_c.get_vec_index( ...
                d.num_layers_period+ind_b1_R);
            obj.ind_z_R(2) = sim_c.get_vec_index( ...
                2*d.num_layers_period+ind_b2_R);
            % Set indices for wavefunction bounding box.
            % Left index at middle of first barrier in period.
            obj.ind_z_L_bound(1) = ...
                sim_c.get_vec_index_mid(d.num_layers_period+1);
            % Right index at middle of reference barrier.
            obj.ind_z_R_bound(1) = sim_c.get_vec_index_mid( ...
                d.num_layers_period+ind_b_ref);
            % Left index at middle of reference barrier.
            obj.ind_z_L_bound(2) = obj.ind_z_R_bound(1);
            % Right index at middle of first barrier of next period.
            obj.ind_z_R_bound(2) = sim_c.get_vec_index_mid( ...
                2*d.num_layers_period+1);
            % Point of lowest potential directly next to thickest barrier.
            obj.ind_E_min = 1 + sim_c.get_vec_index( ...
                d.num_layers_period+1);
            % Define Energy period in J.
            obj.Eperiod = max(diff(sim_c.vec_V_0));
        end
        
        function indices = get_indices(obj, n_scans)
            indices = obj.ind_z_L(n_scans):obj.incr_z:obj.ind_z_R(n_scans);
        end
        function index_end = get_index_bc_end(obj, size_m, n_scans)
            index_end = sub2ind(size_m, 2, obj.ind_z_R(n_scans));
        end
    end
    methods (Static)
        function bc_start = get_bc_start()
            bc_start = [0; 1];
        end
    end
end