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

classdef boundary_condition < handle
    %boundary_condition Base class for simulation boundary conditions
    % for the tm_solver class to calculate the system eigenstates.
    properties (SetAccess = protected)
        name = '' % Boundary condition name.
        ind_z_L % Index start of structure.
        ind_z_R % Index end of structure.
        ind_interp % Index of z-position for interpolation.
        ind_E_min % Lowest energy value of potential.
        incr_z % Increment size of z-vector.
        Eperiod % Energy of bias over one period.
        delta_psi % Cut-off parameter for wavefunction copying to
        % the leftmost period before interpolation.
        dE % Termination tolerance dE for root of eigenenergies.
        num_scans = 1; % Number of energy scan intervals.
        ind_z_L_bound % Bounding box index left.
        ind_z_R_bound % Bounding box index right.
    end
    methods
        % Constructs boundary_condition.
        function obj = boundary_condition(name)
            obj.name = name;
            obj.register_bc();
        end
        % Set Eperiod.
        function set_Eperiod(obj, E_p)
            obj.Eperiod = E_p;
        end
        % Get left vector for localization interval.
        function z_L = get_z_L(obj, vec_z, n_scans)
            z_L = vec_z >= vec_z(obj.ind_z_L_bound(n_scans));
        end
        % Get right vector for localization interval.
        function z_R = get_z_R(obj, vec_z, n_scans)
            z_R = vec_z <= vec_z(obj.ind_z_R_bound(n_scans));
        end
    end
    methods (Static)
        function m_bootstrap = register_bc()
            % Determines boundary condition regarding there key value.
            % Condition type (tb=tightbinding, ext=extended, ez=EZ states,
            % QCD= Quantum cascade detector states).
            keySet = {'tb', 'ext', 'ext-ez', 'tb-ez', 'QCD'};
            b_tb = @(d, s, sim_c, key) boundary_condition_tb(d, s, sim_c);
            b_ext = @(d, s, sim_c, key) boundary_condition_ext(d, ...
                s, sim_c);
            b_ext_ez = @(d, s, sim_c, key) ...
                boundary_condition_ext_ez(d, s, sim_c);
            b_tb_ez = @(d, s, sim_c, key) ...
                boundary_condition_tb_ez(d, s, sim_c);
            b_QCD = @(d, s, sim_c, key) boundary_condition_QCD(d, ...
                s, sim_c);
            valueSet = {b_tb, b_ext, b_ext_ez, b_tb_ez, b_QCD};
            m_bootstrap = containers.Map(keySet, valueSet);
        end
        
        function bc = create_instance(d, s, sim_c)
            % Creates instance of object boundary condition.
            bc_bootstraps = boundary_condition.register_bc();
            bc_reg = bc_bootstraps(s.basis_sp);
            bc = bc_reg(d, s, sim_c);
        end
    end
    methods (Abstract)
        % Gets indices for TMM simulation domain.
        indices = get_indices(obj, n_scans)
        % Gets index of boundary condition at endpoint.
        index_end = get_index_bc_end(obj, size_m, n_sans)
        % Gets boundary condition at startpoint.
        bc_start = get_bc_start()
    end
end
