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

classdef sim_constants < handle
    %sim_constants Class containing the simulation constants for the
    % tm_solver.
    properties
        vec_num_indices % Vector containing number of TM indices per layer.
        vec_alpha % Simulation vector of nonparabolicity parameter alpha.
        vec_parab % Simulation vector of nonparabolicity from material.
        vec_meff % Simulation vector of effective mass.
        vec_V_0 % Simulation vector of biased conduction band profile.
        vec_z_tm % Simulation vector of position vector with double values
        % at interfaces.
        vec_z % Simulation vector of reduced position vector.
        vec_ind_tm % Simulation vector of indices with double values
        % at interfaces
        z_unit = 1e-10; % Factor to change z-grid units from Angstrom to m.
        dz_poisson % Uniform grid spacing Poisson solver.
        vec_z_poisson % Uniform grid used in the Poisson solver.
        l_period; % Length of one period.
        E_period; % Energy period, bias drop over a QCL period.
        vec_rhod; % Vector with donor charge density.
        rel_permittivity; % Relative permittivity.
        num_periods_wf; % Number of periods to fill with wavefunctions.
        num_wavefct; % Number of wavefunctions per period.
        Emin = []; % Minimum energy for root finding algorithm.
    end
    
    methods
        function obj = sim_constants(d, s)
            % Constructs object with simulation constants for the
            % transfer matrix method solver.
            % Grid spacing.
            dz = s.dz_sp;
            obj = obj.set_vec_z(d.int_pos, dz);
            [obj.vec_z, obj.vec_ind_tm, ~] ...
                = unique(obj.vec_z_tm);
            obj.dz_poisson = min(diff(obj.vec_z));
            obj.vec_z_poisson = ...
                obj.vec_z(1):obj.dz_poisson:obj.vec_z(end);
            obj.vec_alpha = obj.set_vec_sim_const(d.get_vec_alpha(s.T));
            obj.vec_parab = obj.set_vec_sim_const(d.get_vec_parab(s.T));
            obj.vec_meff = obj.set_vec_sim_const(d.get_vec_meff(s.T));
            obj.vec_V_0 = obj.set_vec_sim_const(d.get_vec_Ec(s.T)) ...
                +obj.get_vec_bias(d, s);
            obj.l_period = d.l_period;
            obj.E_period = d.get_Eperiod(s.V) * phys_const.e0;
            obj = obj.set_rhod(d);
            obj.rel_permittivity = d.rel_permittivity;
            obj.num_periods_wf = d.num_periods - 1;
            obj.num_wavefct = s.num_wavefct;
            
        end
        % Set vector containing z-coordinate.
        function obj = set_vec_z(obj, int_pos, dz)
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
            %
        end
        % Set vector for corresponding simulation constant.
        function vec_sim_const = ...
                set_vec_sim_const(obj, sim_const_layers)
            % Fill vector for specific simulation constant.
            vec_sim_const = [];
            for i = 1:length(obj.vec_num_indices)
                vec_sim_const = [vec_sim_const, ...
                    sim_const_layers(i) .* ...
                    ones(1, obj.vec_num_indices(i))];
            end
        end
        % Gets effective mass vector energy resolved, taking into account
        % nonparabolicity effects.
        function meff_np = get_meff_np(obj, E, Vt)
            x = (E - Vt) ./ obj.vec_alpha;
            meff_np = obj.vec_meff .* ((1 + x) .* (x >= 0) + ...
                (1 - sqrt(1-4*x)) / 2 ./ x .* (x < 0));
        end
        % Gets vector with applied bias.
        function vec_bias = get_vec_bias(obj, d, s)
            % Potential at beginning and end of structure in V.
            phiniz = 0;
            phi_fin = s.V * 1e5 * 1e-10 * d.tot_length;
            % Bias in J.
            vec_bias = phys_const.e0 * ((obj.vec_z_tm ...
                -obj.vec_z_tm(1)) * (phi_fin - phiniz) ...
                / (obj.vec_z_tm(end) - obj.vec_z_tm(1)));
            
        end
        % Set rho vector of positive charged donors.
        function obj = set_rhod(obj, d)
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
            % Donor density in cm^(-3).
            obj.vec_rhod = rhod;
        end
        %
        % Get index of simulation vector for a specifc layer index in the
        % device.
        function index_vec = get_vec_index(obj, index_dev_layer)
            index_vec = sum(obj.vec_num_indices(1:index_dev_layer));
        end
        % Get index of simulation vector for a specifc layer index in the
        % device by taking into account the last layer length only half.
        function index_vec = get_vec_index_mid(obj, index_dev_layer)
            index_vec = sum(obj.vec_num_indices(1:index_dev_layer-1));
            index_vec = index_vec ...
                +ceil(obj.vec_num_indices(index_dev_layer)/2);
        end
    end
end