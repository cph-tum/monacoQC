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

classdef boundary_condition < handle
    % Base class for boundary conditions for the tm_solver class to
    % calculate the eigenstates of the system.

    properties (SetAccess = protected)
        name = '' % char| string: Name of boundary condition.
        ind_z_L % scalar: Start index of local interval.
        ind_z_R % scalar: End index of local intercal.
        ind_interp % scalar: Index of z-position for interpolation.
        ind_E_min % scalar: Lowest energy value of potential.
        incr_z % scalar: Increment size of z-vector.
        Eperiod % scalar: Energy of bias over one period.
        delta_psi % scalar: Cut-off parameter for wavefunction copying to the leftmost period before interpolation.
        dE % scalar: Termination tolerance dE for root of eigenenergies.
        num_scans = 1 % scalar: Number of energy scan intervals.
        ind_z_L_bound % scalar: Bounding box index left.
        ind_z_R_bound % scalar: Bounding box index right.
    end

    methods
        function obj = boundary_condition(name)
            % Constructs an object of type boundary_condition.
            %
            % Syntax:
            %   obj = boundary_condition(name)
            %
            % Input Arguments:
            %   name (char | string): Name of boundary condition. Valid
            %     names are ``tb``, ``ext``, ``ez`` and ``QCD``.

            obj.name = name;
            obj.register_bc();
        end

        function set_Eperiod(obj, E_p)
            % Sets Eperiod.
            obj.Eperiod = E_p;
        end

        function z_L = get_z_L(obj, vec_z, n_scans)
            % Gets left vector for localization interval.
            z_L = vec_z >= vec_z(obj.ind_z_L_bound(n_scans));
        end

        function z_R = get_z_R(obj, vec_z, n_scans)
            % Get right vector for localization interval.
            z_R = vec_z <= vec_z(obj.ind_z_R_bound(n_scans));
        end
    end

    methods (Static)
        function m_bootstrap = register_bc()
            % Determines boundary condition regarding there key value.
            % Condition type (tb=tightbinding, ext=extended, ez=EZ states,
            % QCD= Quantum cascade detector states).
            %
            % Syntax:
            %   m_bootstrap = register_bc()

            keySet = {'tb', 'ext', 'ez', 'QCD'};
            b_tb = @(d, s, sim_c, key) boundary_condition_tb(d, s, sim_c);
            b_ext = @(d, s, sim_c, key) boundary_condition_ext(d, ...
                s, sim_c);
            b_ez = @(d, s, sim_c, key) boundary_condition_ez(d, s, sim_c);
            b_QCD = @(d, s, sim_c, key) boundary_condition_QCD(d, ...
                s, sim_c);
            valueSet = {b_tb, b_ext, b_ez, b_QCD};
            m_bootstrap = containers.Map(keySet, valueSet);
        end

        function bc = create_instance(d, s, sim_c)
            % Creates an instance of object boundary condition.
            %
            % Syntax:
            %   bc = create_instance(d, s, sim_c)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %   sim_c (sim_constants): Contains simulation constants.

            bc_bootstraps = boundary_condition.register_bc();
            bc_reg = bc_bootstraps(s.basis_sp);
            bc = bc_reg(d, s, sim_c);
        end
    end

    methods (Abstract)
        indices = get_indices(obj, n_scans) % Gets indices for TMM
        % simulation domain.
        index_end = get_index_bc_end(obj, size_m, n_sans) % Gets index of
        % boundary condition at endpoint.
        bc_start = get_bc_start() % Gets boundary condition at startpoint.
    end
end
