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

classdef InAs < binary
    %InAs Represents InAs material.
    %
    properties
    end
    %
    methods
        % Constructs InAs.
        function obj = InAs()
            name = 'InAs';
            % intialize base class and properties at temperature t
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table III
            obj.c_parab = 0;
            obj.param = parameter(-0.892, 0.390, 15.1, 0.417, 0.833, ...
                0.453, 0.396, 6.0583-2.74e-5*300, ...
                1.00, -5.08, -2.90, 21.5, 0, -1.8);
            % Source elastic constants C11, C12, C44: Ioffe
            obj.param_T = parameter(0, 0, 0, [0.276e-3, 94], ...
                0, 0, 0, 2.74e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
