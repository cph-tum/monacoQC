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

classdef parameter < handle
    % bowing_factor includes all used bowing factors C.
    properties (SetAccess = private)
        name = '' % Material name
        conc = 0 % Concentration
        
        Ec % Conduction band (eV)
        Dso % Split off energy gap (eV)
        eps % Relative permittivity
        Eg % Energy gap (eV)
        c11 % Elastic constant (TPa)
        c12 % Elastic constant (TPa)
        c44 % Elastic constant (TPa)
        a_lattice % Lattice constant(A)
        a_v % Valence band deformation potential (eV)
        a_c % Conduction band deformaton potential (eV)
        F % Kane parameter()
        Ep % Matrix element (eV)
        m_eff % Effective mass (1/m0)
        b % Shear deformation potential (eV)
    end
    
    methods
        function obj = parameter(ec, dso, eps, eg, c11, c12, ...
                c44, a, av, ac, f, ep, meff, b)
            obj.Ec = ec;
            obj.Dso = dso;
            obj.eps = eps;
            obj.Eg = eg;
            obj.c11 = c11;
            obj.c12 = c12;
            obj.c44 = c44;
            obj.a_lattice = a;
            obj.a_v = av;
            obj.a_c = ac;
            obj.F = f;
            obj.Ep = ep;
            obj.m_eff = meff;
            obj.b = b;
        end
    end
    
end
