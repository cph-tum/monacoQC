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

classdef GaAs_AlGaAs < material_system
    %GaAs_AlGaAs Represents material system GaAs/AlGaAs.
    %
    properties (SetAccess = protected)
    end
    %
    methods
        % Constructs material system.
        function obj = GaAs_AlGaAs(name)
            obj = obj@material_system(name);
            % Interface roughness parameter
            obj.Gamma = 1e-8;
            obj.Delta = 1.2e-10;
            % Material system flag
            obj.flag_material_system = 2;
        end
    end
end
