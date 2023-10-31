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

classdef GaAs < binary
    %GaAs Represents GaAs material.
    %
    properties
    end
    %
    methods
        % Constructs GaAs.
        function obj = GaAs()
            name = 'GaAs';
            % intialize base class and properties
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table I
            obj.c_parab = 3.2;
            % Source: Ekenberg et al. 1989, Hendorfer et al. 1993
            % https://doi.org/10.1103/PhysRevB.48.2328
            % http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/mechanic.html
            obj.param = parameter(0, 0.341, 12.9, 1.519, 1.217, ...
                0.546, 0.616, 5.65325-3.88e-5*300, ...
                1.16, -7.17, -1.94, 28.8, 0, -2.0);
            % Source elastic constants C11, C12, C44: ioffe
            obj.param_T = parameter(0, 0, 0, [0.5405e-3, 204], ...
                -1.44e-4, -6.4e-6, -7e-5, 3.88e-5, 0, 0, 0, 0, 0, 0);
        end
    end
    
end
