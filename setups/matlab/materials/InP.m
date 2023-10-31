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

classdef InP < binary
    %InP Represents InP material.
    %
    properties
    end
    %
    methods
        % Constructs InP.
        function obj = InP()
            name = 'InP';
            % intialize base class and properties at temperature t
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table VI
            % Source: aftershoq
            obj.param = parameter(-0.2854, 0.108, 12.5, 1.4236, 1.011, ...
                0.561, 0.456, 5.8697-2.79e-5*300, ...
                1.27, -5.04, -1.31, 20.7, 0, -2);
            obj.param_T = parameter(0, 0, 0, [0.363e-3, 162], ...
                0, 0, 0, 2.79e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
