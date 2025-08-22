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

classdef binary < material & dynamicprops

    properties (SetAccess = public)
        temp
    end
    properties (Dependent)
        name_III
        name_V
    end

    methods
        function obj = binary(name, temperature)
            arguments
                name{mustBeNonzeroLengthText}
                temperature = 300
            end

            filename = char(name+".yaml");
            if isempty(which(filename))
                error("Material file '%s' not found.", filename)
            end
            material_data = yaml.ReadYaml(filename);
            if ~strcmp(material_data.composition, "binary")
                error("Material '%s' is not a binary.", name)
            end

            obj.name = char(material_data.name);
            obj.temp = temperature;

            for prop_cell = material_data.parameters
                prop = prop_cell{1};
                p = addprop(obj, prop.name);
                p.SetAccess = "private";
                if ~isfield(prop, "temp_dep")
                    obj.(prop.name) = prop.value;
                elseif strcmp(prop.temp_dep.model, "linear")
                    p.GetMethod = @(obj)binary.linear_dep(obj, prop);
                elseif strcmp(prop.temp_dep.model, "Varshni")
                    p.GetMethod = @(obj)binary.varshni_dep(obj, prop);
                else
                    error("%s: model not implemented.", prop.name);
                end
            end
        end
    end

    % getter methods for dependent properties
    methods
        function out = get.name_III(obj)
            caps = find(isstrprop(obj.name, "upper"));
            out = obj.name(caps(1):caps(2)-1);
        end

        function out = get.name_V(obj)
            caps = find(isstrprop(obj.name, "upper"));
            out = obj.name(caps(2):end);
        end
    end

    % static methods for parameter models
    methods (Static)
        function out = linear_dep(obj, prop)
            out = prop.value + prop.temp_dep.param * obj.temp;
        end

        function out = varshni_dep(obj, prop)
            % Source: Vurgaftman et al. 2001, Eq. (2.13)
            % https://doi.org/10.1063/1.1368156
            out = prop.value - prop.temp_dep.param{1} * obj.temp^2 / ...
                (obj.temp + prop.temp_dep.param{2});
        end
    end
end
