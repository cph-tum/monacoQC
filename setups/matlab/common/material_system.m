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

classdef material_system < handle
    %material_system TODO add description.
    %
    properties (SetAccess = protected)
        name = ''; % Material system name
        Gamma = 0; % Interface roughness corr. length in m
        Delta = 0; % Interface roughness height in m
        flag_material_system = 0;
    end
    %
    methods
        % Constructs material system.
        function obj = material_system(name)
            obj.name = name;
        end
        %
        % Gets interface roughness corr. length.
        function gamma = get.Gamma(obj)
            gamma = obj.Gamma;
        end
        %
        % Gets interface roughness height.
        function delta = get.Delta(obj)
            delta = obj.Delta;
        end
    end
end
