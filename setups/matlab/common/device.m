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

classdef device < handle
    %device Represents an active region design.
    % Maintains internally a list of layers that represent the active
    % region.
    %
    properties (SetAccess = protected)
        layers = {} % Internal list of layers
        num_periods = 0 % Number of periods
        substrate % Substrate
        material_system % Material system
        orientation = '001'; % Interface orientation
        waveguide % Waveguide description
    end
    %
    properties (Dependent)
        doping = [] % Doping matrix
        tot_length % Total device length in Angstrom
        tot_length_well % Total well layer length in Angstrom
        tot_length_barrier % Total barrier layer length in Angstrom
        int_pos % Well-Barrier Interface position in Angstrom
        dens_sheet % Total sheet density (in m^-2) in a period
        dens_carrier % Total charge carrier density (in m^-3)
        l_dop % length of doped layers in Angstrom
        rel_permittivity % Relative permittivity
        l_period % Length of one period in Angstrom
        CBO % Conduction band offset in eV
        num_layers_period % Number of layers per period.
    end
    %
    properties
        n_eff % Effective refractive index
    end
    %
    methods
        % Constructs device.
        function obj = device(period, num_periods)
            % number of periods
            obj.num_periods = num_periods;
            obj.layers = obj.set_layers(period, num_periods);
        end
        
        % Determines material system.
        function sys = set_material_system(obj, period)
            % available materials
            materials = {'GaAs', 'AlGaAs', 'Si', 'SiGe', 'AlAs', ...
                'InAs', 'InGaAs', 'InAlAs'};
            
            % material histogram
            map = containers.Map(materials, zeros(size(materials)));
            for i = 1:length(period)
                count = map(period{i}.material.name);
                map(period{i}.material.name) = count + 1;
            end
            % determine matarial system
            if (map('GaAs') + map('AlGaAs') == length(period))
                sys = GaAs_AlGaAs('GaAs/AlGaAs');
                for i = 1:length(period)
                    if (strcmp(period{i}.material.name, materials{1}))
                        period{i}.layer_type = 'w';
                    elseif (strcmp(period{i}.material.name, materials{2}))
                        period{i}.layer_type = 'b';
                    end
                end
                obj.substrate = GaAs();
            elseif (map('InGaAs') + map('InAlAs') == length(period))
                sys = InGaAs_InAlAs('InGaAs/InAlAs');
                for i = 1:length(period)
                    if (strcmp(period{i}.material.name, materials{7}))
                        period{i}.layer_type = 'w';
                    elseif (strcmp(period{i}.material.name, materials{8}))
                        period{i}.layer_type = 'b';
                    end
                end
                obj.substrate = InP();
            elseif (map('GaAs') + map('AlAs') == length(period))
                error('Not implemented yet!');
            elseif (map('Si') + map('SiGe') == length(period))
                error('Not implemented yet!');
            else
                error('Could not determine material system!');
            end
        end
        % Determines total length of device in Angstrom.
        function l = get.tot_length(obj)
            l = 0;
            for i = 1:length(obj.layers)
                l = l + obj.layers{i}.length;
            end
        end
        % Determines total length of well layers in Angstrom.
        function l = get.tot_length_well(obj)
            l = 0;
            ind_w = obj.get_ind_well;
            for i = 1:length(ind_w)
                l = l + obj.layers{ind_w(i)}.length;
            end
        end
        % Determines total length of barrier layers in Angstrom.
        function l = get.tot_length_barrier(obj)
            l = 0;
            ind_b = obj.get_ind_barrier;
            for i = 1:length(ind_b)
                l = l + obj.layers{ind_b(i)}.length;
            end
        end
        %
        % Determines layer interface position in Angstrom.
        function r = get.int_pos(obj)
            r = zeros(1, length(obj.layers)+1);
            r(1) = -obj.tot_length;
            r(length(obj.layers)+1) = 0;
            for i = 2:length(obj.layers)
                r(i) = r(i-1) + obj.layers{i-1}.length;
            end
        end
        %
        % Determines total sheet density in m^-2.
        function n2D = get.dens_sheet(obj)
            n2D = 0;
            for i = 1:length(obj.layers) - 1
                if (obj.layers{i}.doping ~= 0)
                    n2D = n2D + obj.layers{i}.length * 1e-8 ...
                        * obj.layers{i}.doping;
                end
            end
            n2D = n2D / obj.num_periods * 1e4;
        end
        %
        % Determines total charge carrier density in m^-3
        function n3D = get.dens_carrier(obj)
            n3D = obj.dens_sheet / obj.l_period / 1e-10;
        end
        % Calculates total length of doped layers in Angstrom.
        function l = get.l_dop(obj)
            l = 0;
            for i = 1:length(obj.layers) - 1
                if (obj.layers{i}.doping ~= 0)
                    l = l + obj.layers{i}.length;
                end
            end
        end
        %
        % Calculates length of one period in Angstrom.
        function l = get.l_period(obj)
            l = 0;
            for i = 1:obj.num_layers_period
                l = l + obj.layers{i}.length;
            end
        end
        %
        % Prepares doping matrix.
        function dop = get.doping(obj)
            position = 0;
            dop = [];
            for i = 1:length(obj.layers) - 1
                if (obj.layers{i}.doping > 0)
                    if (size(dop, 1) > 0 && obj.layers{i}.doping ...
                            == obj.layers{i - 1}.doping && ...
                            strcmp(obj.layers{i}.doping_type, ...
                            obj.layers{i - 1}.doping_type))
                        dop(end, 2) = dop(end, 2) ...
                            +obj.layers{i}.length;
                    else
                        doping_entry = zeros(1, 4);
                        %
                        % start position of doped layer in Angstrom
                        doping_entry(1) = position;
                        %
                        % end position of doped layer in Angstrom
                        doping_entry(2) = position + obj.layers{i}.length;
                        % doping concentration in 1/cm^3
                        doping_entry(3) = obj.layers{i}.doping;
                        % doping type. for electron donors this is set to 1
                        if (strcmp(obj.layers{i}.doping_type, 'n'))
                            doping_entry(4) = 1;
                        else
                            error('Unknown doping type!');
                        end
                        dop = [dop; doping_entry];
                    end
                end
                %
                % update layer start position
                position = position + obj.layers{i}.length;
            end
        end
        %
        function eps = get.rel_permittivity(obj)
            % Gets static relative permittivity of device.
            ind_w = obj.get_ind_well;
            eps_mat = zeros(1, length(ind_w));
            for i = 1:length(ind_w)
                eps_mat(i) = (obj.layers{ind_w(i)}.material.get_eps() * ...
                    obj.layers{ind_w(i)}.length);
            end
            l_well = obj.tot_length_well;
            eps = sum(eps_mat) / l_well;
        end
        %
        % Set effective refractive index.
        function set.n_eff(obj, neff)
            obj.n_eff = neff;
        end
        %
        % Gets effective refractive index.
        function neff = get.n_eff(obj)
            if (isempty(obj.n_eff))
                neff = sqrt(obj.rel_permittivity);
            else
                neff = obj.n_eff;
            end
        end
        %
        % Gets conduction band offset in eV.
        function cbo = get_CBO(obj, T)
            % TODO Monte Carlo requires one value for the CBO per device.
            % we assume here that this is the difference between the first
            % and second material.
            % 1) revise whether this makes sense, in particular for
            %    devices with more than two different materials.
            % 2) find a better way to determine representative barrier and
            %    well material.
            %
            b = zeros(1, length(obj.layers));
            for i = 2:length(obj.layers)
                b(i) = obj.layers{i-1}.material.get_Ec(obj.substrate, ...
                    T, obj.orientation) - obj.layers{i}.material.get_Ec ...
                    (obj.substrate, T, obj.orientation);
            end
            b_count = gt(b, 0);
            cbo = max(nonzeros(unique(b.*b_count)));
        end
        %
        % Calculates total mismatch h and total critical thickness h_crit
        % in one period. For abs(h)> abs(hcrit) the structure
        % is estimated to relaxe.
        function h2 = get_strain(obj, T, N)
            % Lattice constant of the substrate
            asub = obj.substrate.get_lattice_constant(T);
            h = 0;
            h_crit = 0;
            
            for i = 1:(length(obj.layers) - 1) / 5 * 2
                h = h + ...
                    obj.layers{1, i}.material.get_eps_II( ...
                    obj.substrate(), T) * obj.layers{1, i}.length;
                if (nargin < 3)
                    h_crit = h_crit + ...
                        (obj.layers{1, i}.material ...
                        .get_lattice_constant(T) - asub) / asub * ...
                        obj.layers{1, i} ...
                        .material.get_hcrit(obj.substrate, T, ...
                        obj.orientation);
                else
                    h_crit = h_crit + ...
                        (obj.layers{1, i}.material ...
                        .get_lattice_constant(T) - ...
                        asub) / asub * obj.layers{1, i} ...
                        .material.get_hcrit(obj.substrate, T, ...
                        obj.orientation, N);
                end
            end
            h2.h = h / 2;
            h2.h_crit = h_crit / 2 / ((length(obj.layers) - 1) / 5);
        end
        %
        % Does an active strain compensation of a specific
        % device period, e.g. generated by Bayesian optimization
        function new_x = strain_comp(obj, T, CBO_target)
            cell_mat = obj.find_mat();
            if (size(cell_mat, 2) ~= 2)
                error(['Active strain compensation is implemented ', ...
                    'for', ' devices comprising two materials!']);
            end
            for i = 1:size(cell_mat, 2)
                if (cell_mat{1, i}.n_comp == 2)
                    lb(i) = 0;
                    ub(i) = 0;
                    x0(i) = 0;
                elseif (cell_mat{1, i}.n_comp == 3)
                    lb(i) = 0;
                    ub(i) = 1;
                    x0(i) = cell_mat{1, i}.get_conc();
                else
                    error(['Materials with %d components', ...
                        ' cannot be considered', 'for', ...
                        'active strain compensation!'], ...
                        cell_mat{1, i}.n_comp);
                end
                
            end
            if (nargin < 3)
                CBO_target = obj.get_CBO(T);
            end
            fun = @(x) obj.f_strain(T, x, CBO_target);
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            new_x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
        end
        %
        % Find different materials of device
        function cell_mat = find_mat(obj)
            cell_mat = {};
            for i = 1:size(obj.layers, 2)
                if (i == 1)
                    cell_mat{1} = obj.layers{i}.material;
                else
                    flag = 0;
                    for j = 1:size(cell_mat, 2)
                        flag = flag + eq(cell_mat{j}, ...
                            obj.layers{i}.material);
                    end
                    if (~flag)
                        cell_mat{1+end} = obj.layers{i}.material;
                    end
                end
            end
        end
        %
        function obj = set_orientation(obj, indizes)
            if (ischar(indizes))
                obj.orientation = indizes;
            else
                error('Indizes must be a character array, e.g. 001.');
            end
        end
        %
        % width should be multiple of lattice constant to
        % reduce interface roughness(Wolf 2017)
        function nwidth = get_layer_length_reduced_IR(obj, T)
            z = (length(obj.layers) - 1) / 5;
            nwidth = zeros(1, z);
            % Get cell with used materials.
            mat_d = obj.find_mat();
            for i = 1:z
                for j = 1:size(mat_d, 2)
                    if (eq(obj.layers{i}.material, mat_d{1, j}))
                        % Get lattice constant of material j.
                        a_mat = mat_d{1, j}.get_lattice_constant(T);
                        % Get ratio of layer length and lattice constant.
                        ratio_len_a = round(obj.layers{i}.length/a_mat);
                        % Calculate optimal layer length
                        % for IR scattering reduction.
                        nwidth(i) = ratio_len_a * a_mat;
                    end
                end
            end
        end
        % Get vector containing conduction band offsets for device layers.
        function V_c = get_vec_Ec(obj, T)
            
            % find minimum conduction band edge
            ec_ref = min(cellfun(@(x) (x.material.get_Ec(obj.substrate, ...
                T, obj.orientation)), obj.layers));
            % Generate vector with size of layers.
            V_c = zeros(1, length(obj.layers));
            for i = 1:length(obj.layers)
                V_c(i) = obj.layers{i}.material.get_Ec(obj.substrate, ...
                    T, obj.orientation) - ec_ref;
            end
            % Returns CBO vector in J.
            V_c = phys_const.e0 * V_c;
        end
        % Get vector containing nonparabolicity coefficient
        % for device layers.
        function parab_xy = get_vec_parab(obj, T)
            % Generate vector with size of layers.
            parab_xy = zeros(1, length(obj.layers));
            for i = 1:length(obj.layers)
                parab_xy(i) = ...
                    obj.layers{i}.material.get_parab(obj.substrate, ...
                    T, obj.orientation);
            end
            % Returns nonparabolicity vector in J^-1.
            parab_xy = parab_xy / phys_const.e0;
        end
        % Get vector containing effective mass for device layers.
        function meff = get_vec_meff(obj, T)
            % Generate vector with size of layers.
            meff = zeros(1, length(obj.layers));
            for i = 1:length(obj.layers)
                meff(i) = obj.layers{i}.material.get_meffII( ...
                    obj.substrate, T, obj.orientation);
            end
            % Returns effective mass in kg.
            meff = meff * phys_const.me;
            
        end
        % Get vector containing nonparabolicity parameter alpha
        % for device layers.
        function vec_alpha = get_vec_alpha(obj, T)
            % Generate vector with size of layers.
            vec_alpha = zeros(1, length(obj.layers));
            for i = 1:length(obj.layers)
                vec_alpha(i) = obj.layers{i}.material.get_1_alpha( ...
                    obj.substrate, T, obj.orientation);
            end
            % Returns alpah in J.
            vec_alpha = phys_const.e0 * vec_alpha;
        end
        %
        function Eperiod = get_Eperiod(obj, V)
            % Calculates applied energy over one period
            phi_fin = V * 1e5 * 1e-10 * obj.tot_length;
            Eperiod = abs(phi_fin/min((-1)*obj.tot_length)*obj.l_period);
            
        end
        
        function ind_b = get_ind_barrier(obj)
            % Gets indices of layers with barrier materials.
            type = obj.get_type_layer;
            ind_b = find([type{:}] == 'b');
        end
        %
        function ind_w = get_ind_well(obj)
            % Gets indices of layers with well materials.
            type = obj.get_type_layer;
            ind_w = find([type{:}] == 'w');
        end
        %
        function type_l = get_type_layer(obj)
            % Find type layers.
            type_l = cell(1, length(obj.layers));
            for k = 1:length(obj.layers)
                type_l{k} = obj.layers{k}.layer_type;
            end
        end
        %
        function num_l_period = get.num_layers_period(obj)
            % Get number of layers per period.
            num_l_period = (length(obj.layers) - 1) / obj.num_periods;
        end
        % Set waveguide description
        function obj = set_waveguide(obj, wg)
            obj.waveguide = wg;
        end
    end
    
    methods (Access = private)
        % Calculates strain for specific material combination (x)
        % weighted with the CBO of the given device object.
        function res = f_strain(obj, T, x, CBO_target)
            mat_d = obj.find_mat();
            mat_test = obj.copy_mat();
            Lmat = zeros(size(mat_d, 2), 1);
            nlayers = (length(obj.layers) - 1) / 5;
            % calculates total length of different materials in
            % one period
            for i = 1:nlayers
                for j = 1:length(mat_d)
                    if (eq(obj.layers{i}.material, mat_d{1, j}))
                        Lmat(j) = Lmat(j) + obj.layers{1, i}.length;
                    end
                end
            end
            for i = 1:size(mat_test, 2)
                mat_test{1, i}.set_conc(x(i));
                epsII(i) = mat_test{1, i}.get_eps_II(obj.substrate, T);
                a(i) = mat_test{1, i}.get_lattice_constant(T);
                G(i) = mat_test{1, i}.get_G(T, obj.orientation);
                Ec_test(i) = mat_test{1, i}.get_Ec(obj.substrate, ...
                    T, obj.orientation);
            end
            delta_CBO = abs(abs(Ec_test(1)-Ec_test(2))/CBO_target-1);
            res = abs(epsII.*G.*a*Lmat) ...
                * (1 - exp((-1)*(delta_CBO)^2));
        end
        % Copy materials
        function copy_mat = copy_mat(obj)
            mat = obj.find_mat();
            copy_mat = cell(size(mat));
            for i = 1:size(mat, 2)
                copy_mat{1, i} = copy(mat{1, i});
            end
        end
        % Allocates layers.
        function l = set_layers(obj, period, num_periods)
            % Determine material system
            obj.material_system = obj.set_material_system(period);
            % Find thickest barrier.
            lt = zeros(length(period), 1);
            for k = 1:length(period)
                lt(k) = period{k}.length;
                type_l{k} = period{k}.layer_type;
            end
            index_layer_period_b = find([type_l{:}] == 'b');
            % Returns index of thickest barrier.
            [~, I_b] = max(lt(index_layer_period_b));
            I = index_layer_period_b(I_b);
            % Change period order, start with thickest barrier.
            period_ordered = {period{I:end}, period{1:I-1}};
            period_ordered = reshape(period_ordered, 1, ...
                length(period_ordered));
            % Generate layer cell with 5 periods plus one barrier at the
            % end.
            l = {};
            for i = 1:num_periods
                l = [l, period_ordered];
            end
            l = [l, l(1)];
        end
    end
    
end
