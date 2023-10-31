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

classdef AlP < binary
    %InAs Represents InAs material.
    %
    properties
    end
    %
    methods
        % Constructs InP.
        function obj = AlP()
            name = 'AlP';
            % intialize base class and properties at temperature t
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table VI
            % Source: aftershoq
            obj.param = parameter(1.171, 0.07, 9.8, 3.63, 1.330, ...
                0.63, 0.615, 5.4672-2.92e-5*300, ...
                3, -5.7, -0.65, 17.7, 0, -1.5);
            obj.param_T = parameter(0, 0, 0, [0.5771e-3, 372], ...
                0, 0, 0, 2.92e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
