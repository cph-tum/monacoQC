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

classdef sim_constants < handle
    % Class containing the simulation constants for the tm_solver.

    properties
        vec_num_indices % vector: Contains number of TM indices (grid points) per layer.
        vec_alpha % vector: Nonparabolicity parameter alpha at each grid point.
        vec_parab % vector: Nonparabolicity from material at each grid point.
        vec_meff % vector: Effective mass at each grid point [kg].
        vec_V_0 % vector: Biased conduction band profile [J].
        vec_z_tm % vector: Position vector with double values at interfaces [Angstrom].
        vec_z % vector: Reduced position vector [Angstrom].
        vec_ind_tm % vector: Grid indices with double values at interfaces.
        z_unit = 1e-10 % scalar: Factor to change z-grid units from Angstrom to m.
        dz_poisson % scalar: Uniform grid spacing for Poisson solver [Angstrom].
        vec_z_poisson % vector: Uniform grid used in the Poisson solver [Angstrom].
        l_period % scalar: Length of one period [Angstrom].
        E_period % scalar: Energy period, bias drop over one QCL period [J].
        vec_rhod % vector: Donor charge density [1/m^3].
        rel_permittivity % scalar: Effective relative permittivity [-].
        num_periods_wf % scalar: Number of periods to fill with wavefunctions.
        num_wavefct % scalar: Number of wavefunctions per period.
    end

    methods
        function obj = sim_constants(d, s)
            % Constructs an object of type sim_constants.
            %
            % Syntax:
            %   obj = sim_constants(d, s)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.

            % Grid spacing.
            dz = s.dz_sp;
            obj = obj.set_vec_z(d.int_pos, dz);
            [obj.vec_z, obj.vec_ind_tm, ~] ...
                = unique(obj.vec_z_tm);
            obj.dz_poisson = min(diff(obj.vec_z));
            obj.vec_z_poisson = ...
                obj.vec_z(1):obj.dz_poisson:obj.vec_z(end);
            obj.vec_alpha = obj.set_vec_sim_const(d.get_vec_alpha);
            obj.vec_parab = obj.set_vec_sim_const(d.get_vec_parab);
            obj.vec_meff = obj.set_vec_sim_const(d.get_vec_meff);
            obj.vec_V_0 = obj.set_vec_sim_const(d.get_vec_Ec) ...
                +obj.get_vec_bias(d, s);
            obj.l_period = d.l_period;
            obj.E_period = d.get_Eperiod(s.V) * phys_const.e0;
            obj = obj.set_rhod(d);
            obj.rel_permittivity = d.rel_permittivity;
            obj.num_periods_wf = d.num_periods - 1;
            obj.num_wavefct = s.num_wavefct;
        end

        function obj = set_vec_z(obj, int_pos, dz)
            % Sets vector containing z-coordinates.
            %
            % Syntax:
            %   obj = set_vec_z(obj, int_pos, dz)
            %
            % Input Arguments:
            %   int_pos (vector): Position of interface positions.
            %   dz (scalar): Grid spacing in z-direction.

            z = [];
            vec_n_ind = zeros(1, length(int_pos)-1);
            for k = 1:length(int_pos) - 1
                z_k = int_pos(k):dz:int_pos(k+1);
                if (z_k(end) ~= int_pos(k+1))
                    z_k = [z_k, int_pos(k+1)];
                end
                vec_n_ind(k) = length(z_k);
                z = [z, z_k];
            end

            obj.vec_z_tm = z;
            obj.vec_num_indices = vec_n_ind;
        end

        function vec_sim_const = ...
                set_vec_sim_const(obj, sim_const_layers)
            % Creates vector of specific simulation constant at each grid
            % point based on the material properties of the individual
            % layers.
            %
            % Syntax:
            %   vec_sim_const = set_vec_sim_const(obj, sim_const_layers)
            %
            % Input Arguments:
            %   sim_const_layers (vector): Simulation constant for each layer.
            %
            % Output Arguments:
            %    vec_sim_const (vector): Simulation constant at each grid point.

            % Fill vector for specific simulation constant.
            vec_sim_const = [];
            for i = 1:length(obj.vec_num_indices)
                vec_sim_const = [vec_sim_const, ...
                    sim_const_layers(i) .* ...
                    ones(1, obj.vec_num_indices(i))];
            end
        end

        function meff_np = get_meff_np(obj, E, Vt)
            % Gets effective mass vector energy resolved, taking into
            % account nonparabolicity effects.
            %
            % Syntax:
            %   meff_np = get_meff_np(obj, E, Vt)
            %
            % Input Arguments:
            %   E (scalar): Eigenenergy [J].
            %   Vt (vector): Potential profile [J].
            %
            % Output Arguments:
            %   meff_np (vector): Energy dependent effective masses [kg].

            x = (E - Vt) ./ obj.vec_alpha;
            meff_np = obj.vec_meff .* ((1 + x) .* (x >= 0) + ...
                (1 - sqrt(1-4*x)) / 2 ./ x .* (x < 0));
        end

        function vec_bias = get_vec_bias(obj, d, s)
            % Creates vector with potential profile caused by applied bias.
            %
            % Syntax:
            %   vec_bias = get_vec_bias(obj, d, s)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %
            % Output Arguments:
            %   vec_bias (vector): Potential profile of applied bias [J].

            % Potential at beginning and end of structure in V.
            phiniz = 0;
            phi_fin = s.V * 1e5 * 1e-10 * d.tot_length;
            % Bias in J.
            vec_bias = phys_const.e0 * ((obj.vec_z_tm ...
                -obj.vec_z_tm(1)) * (phi_fin - phiniz) ...
                / (obj.vec_z_tm(end) - obj.vec_z_tm(1)));

        end

        function obj = set_rhod(obj, d)
            % Sets vector of charge carrier density of ionized donors.
            %
            % Syntax:
            %   obj = set_rhod(obj, d)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.

            % zvu: Position vector without double entries at interfaces
            zvu = obj.vec_z_poisson;
            zvud = zvu' + max(-zvu);
            rhod = zvud * 0;
            % Calculate vector with space charges of positive donor atoms.
            for k = 1:size(d.doping, 1)
                vdop = ((zvud >= d.doping(k, 1)) .* ...
                    (zvud <= d.doping(k, 2)) ...
                    +(zvud > d.doping(k, 1)) ...
                    .* (zvud < d.doping(k, 2))) / 2;
                % Count doping density at interface gridpoints only 50% to
                % get correct donor density.
                rhod = rhod + vdop .* d.doping(k, 3) .* ...
                    d.doping(k, 4) * 1e6;
            end
            % Donor density in m^(-3).
            obj.vec_rhod = rhod;
        end

        function index_vec = get_vec_index(obj, index_dev_layer)
            % Gets index of simulation vector, which corresponds to the
            % position where a specific layer starts.
            %
            % Syntax:
            %   index_vec = get_vec_index(obj, index_dev_layer)
            %
            % Input Arguments:
            %   index_dev_layer (scalar): Index of layer.
            %
            % Output Arguments:
            %   index_vec (scalar): Index of simulation vector corresponing
            %     to the start position of the layer.

            index_vec = sum(obj.vec_num_indices(1:index_dev_layer));
        end

        function index_vec = get_vec_index_mid(obj, index_dev_layer)
            % Gets index of simulation vector, which corresponds to the
            % position where a specific layer starts by taking into account
            % the last layer length only half.
            %
            % Syntax:
            %   index_vec = get_vec_index_mid(obj, index_dev_layer)
            %
            % Input Arguments:
            %   index_dev_layer (scalar): Index of layer.
            %
            % Output Arguments:
            %   index_vec (scalar): Index of simulation vector corresponing
            %     to the start position of the layer.

            index_vec = sum(obj.vec_num_indices(1:index_dev_layer-1));
            index_vec = index_vec ...
                +ceil(obj.vec_num_indices(index_dev_layer)/2);
        end
    end
end
