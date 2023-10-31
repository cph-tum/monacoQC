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

classdef boundary_condition_ext_ez < boundary_condition_ext
    %boundary_condition_ez Boundary condition class for ez eigenstates.
    
    properties
        e_multiplet % Maximum energy spacing of multiplets.
        d_multiplet % Minimum dipole moment of multiplets.
    end
    
    methods
        function obj = boundary_condition_ext_ez(d, s, sim_c)
            % Creates an object with properties describing the simulation
            % boundary conditions for ez state solutions.
            obj = obj@boundary_condition_ext(d, s, sim_c);
            obj.e_multiplet = s.e_multiplet;
            obj.d_multiplet = s.d_multiplet;
        end
         
    end
    methods (Static)
        function bc_start = get_bc_start()
            bc_start = [1; 0];
        end
    end
end