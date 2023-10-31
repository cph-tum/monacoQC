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

classdef InGaAs_InAlAs < material_system
    %InGaAs_InAlAs Represents material system InGaAs/InAlAs.
    %
    properties (SetAccess = protected)
    end
    %
    methods
        % Constructs material system.
        function obj = InGaAs_InAlAs(name)
            obj = obj@material_system(name);
            % Interface roughness parameter
            obj.Gamma = 1e-8;
            obj.Delta = 1e-18 / obj.Gamma;
            % Material system flag
            obj.flag_material_system = 1;
        end
    end
end
