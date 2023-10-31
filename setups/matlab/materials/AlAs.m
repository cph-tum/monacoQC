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

classdef AlAs < binary
    %AlAs Represents AlAs material.
    properties
    end
    %
    methods
        % Constructs AlAs.
        function obj = AlAs()
            name = 'AlAs';
            % intialize base class and properties
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table II
            obj.c_parab = 0;
            obj.param = parameter(0.99, 0.280, 10.06, 3.099, 1.250, ...
                0.534, 0.542, 5.6611-2.90e-5*300, ...
                2.47, -5.64, -0.48, 21.1, 0, -2.3);
            obj.param_T = parameter(0, 0, 0, [0.885e-3, 530], ...
                0, 0, 0, 2.9e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
