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

classdef ternary < material & dynamicprops

    properties (SetAccess = public)
        temp
        conc{mustBeNonnegative, mustBeLessThanOrEqual(conc, 1)}
    end
    properties (SetAccess = private)
        comp1
        comp2
    end
    properties (Dependent)
        name_III
        name_V
    end

    methods
        function obj = ternary(name, concentration, temperature)
            arguments
                name{mustBeNonzeroLengthText}
                concentration = 0.5
                temperature = 300
            end

            filename = char(name+".yaml");
            if isempty(which(filename))
                error("Material file '%s' not found.", filename)
            end
            material_data = yaml.ReadYaml(filename);
            if ~strcmp(material_data.composition, "ternary")
                error("Material '%s' is not a ternary.", name)
            end

            obj.name = char(material_data.name);
            obj.conc = concentration;
            obj.temp = temperature;
            obj.comp1 = binary(material_data.component_1, obj.temp);
            obj.comp2 = binary(material_data.component_2, obj.temp);

            components = [obj.name_III, obj.name_V];
            if length(components) ~= 3
                error("Binary components '%s' and '%s' do not form a ternary.", ...
                    obj.comp1.name, obj.comp2.name)
            end
            if ~strcmp(obj.name, cell2mat(components))
                error("Ternary material name '%s' not consistent with binary components '%s' and '%s' ---> '%s'", ...
                    obj.name, obj.comp1.name, obj.comp2.name, cell2mat(components))
            end

            for prop_cell = material_data.parameters
                prop = prop_cell{1};
                p = addprop(obj, prop.name);
                p.SetAccess = "private";
                if strcmp(prop.name, "E_lo")
                    % phonon energy is treated as a special case (for now)
                    p.GetMethod = @(obj)ternary.lo_phonon_dep(obj, prop);
                elseif ~isfield(prop, "conc_dep")
                    obj.(prop.name) = prop.value;
                elseif strcmp(prop.conc_dep.model, "bowing")
                    p.GetMethod = @(obj)ternary.bowing_dep(obj, prop);
                elseif strcmp(prop.conc_dep.model, "lattice")
                    p.GetMethod = @(obj)ternary.lattice_dep(obj, prop);
                elseif strcmp(prop.conc_dep.model, "linear")
                    p.GetMethod = @(obj)ternary.linear_dep(obj, prop);
                elseif strcmp(prop.conc_dep.model, "quadratic")
                    p.GetMethod = @(obj)ternary.quadratic_dep(obj, prop);
                elseif strcmp(prop.conc_dep.model, "pwl")
                    p.GetMethod = @(obj)ternary.pwl_dep(obj, prop);
                else
                    error("%s: model not implemented.", prop.name)
                end
            end
        end

        function str = name_conc(obj, symbol)
            if nargin ~= 2
                % use concentration
                c1 = sprintf("%.2f", obj.conc);
                c2 = sprintf("%.2f", 1-obj.conc);
            else
                s = char(symbol);
                if isscalar(s) && isletter(s)
                    % use symbol
                    c1 = sprintf("%s", symbol);
                    c2 = sprintf("1-%s", symbol);
                else
                    error("Symbol must be a single letter!")
                end
            end
            g3 = obj.name_III;
            g5 = obj.name_V;
            if length(g3) == 2
                str = sprintf("%s_{%s}%s_{%s}%s", g3{1}, c1, g3{2}, c2, g5{1});
            else
                str = sprintf("%s%s_{%s}%s_{%s}", g3{1}, g5{1}, c1, g5{2}, c2);
            end
        end
    end

    % getter methods for dependent properties
    methods
        function out = get.name_III(obj)
            out = unique({obj.comp1.name_III, obj.comp2.name_III}, "stable");
        end

        function out = get.name_V(obj)
            out = unique({obj.comp1.name_V, obj.comp2.name_V}, "stable");
        end
    end

    % setter methods for public properties
    methods
        function set.temp(obj, value)
            obj.temp = value;
            % set temperature of binary components
            obj.comp1.temp = value;
            obj.comp2.temp = value;
        end
    end

    % static methods for parameter models
    methods (Static)
        function out = bowing_dep(obj, prop)
            if ~isfield(prop.conc_dep, "param")
                bow = 0;
            elseif isscalar(prop.conc_dep.param)
                bow = prop.conc_dep.param;
            else
                a = cell2mat(prop.conc_dep.param);
                bow = a(1) + a(2) * obj.conc;
            end
            c = [obj.conc, -obj.conc * (1 - obj.conc), 1 - obj.conc];
            p = [obj.comp1.(prop.name); bow; obj.comp2.(prop.name)];
            out = c * p;
        end

        function out = lattice_dep(obj, prop)
            c = [obj.conc, 1 - obj.conc];
            p = [obj.comp1.(prop.name) * obj.comp1.a_lattice; ...
                obj.comp2.(prop.name) * obj.comp2.a_lattice];
            out = c * p / obj.a_lattice;
        end

        function out = linear_dep(obj, prop)
            a0 = prop.value;
            a1 = prop.conc_dep.param;
            out = a0 + a1 * obj.conc;
        end

        function out = quadratic_dep(obj, prop)
            a0 = prop.value;
            a = cell2mat(prop.conc_dep.param);
            out = a0 + a(1) * obj.conc + a(2) * obj.conc^2;
        end

        function out = pwl_dep(obj, prop)
            x = cell2mat(prop.conc_dep.concs);
            y = cell2mat(prop.conc_dep.values);
            out = interp1(x, y, obj.conc);
        end

        function out = lo_phonon_dep(obj, prop)
            out = cell2mat(prop.value);
            if isfield(prop, "conc_dep")
                params = cell2mat(prop.conc_dep.param);
                if strcmp(prop.conc_dep.model, "linear")
                    out = out + params * obj.conc;
                elseif strcmp(prop.conc_dep.model, "quadratic")
                    out = out + (params * [obj.conc; obj.conc^2])';
                end
            end
        end
    end
end
