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

classdef InSb < binary
    % InSb Represents InSb material.
    %
    properties
    end
    %
    methods
        % Constructs InSb.
        function obj = InSb()
            name = 'InSb';
            % intialize base class and properties
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table VIII
            % http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/mechanic.html
            obj.param = parameter(-0.484, 0.81, 16.8, 0.235, 0.6847, ...
                0.3735, 0.3111, 6.4794-3.48e-5*300, ...
                0.36, -6.94, -0.23, 23.3, 0, -2);
            obj.param_T = parameter(0, 0, 0, [0.32e-3, 170], ...
                0, 0, 0, 3.48e-5, 0, 0, 0, 0, 0, 0);
        end
    end
    
end
