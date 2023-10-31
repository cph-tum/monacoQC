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

classdef (Abstract) solver < handle
    %solver Abstract solver class for solving the Schrödinger equation.
    properties (SetAccess = protected)
        name = '' % Solver name
    end
    
    methods
        % Constructs solver.
        function obj = solver(name)
            obj.name = name;
        end
    end
    methods (Abstract)
        [eig_system, cond_profile] = solve(obj, func_carr_distr, ...
            eig_system_old);
    end
end
