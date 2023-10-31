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

classdef binary < material
    %binary Contains all required properties of a certain binary material.
    %
    properties (SetAccess = protected)
        param;
        param_T;
    end
    %
    methods
        % Constructs binary material.
        function obj = binary(name)
            obj = obj@material(name);
            obj.n_comp = 2;
        end
        % Gets name of the material
        function name_mat = get_name(obj)
            name_mat = obj.name;
        end
        % Calculates the band gap energy without strain.
        function eg = get_Eg0(obj, T)
            % Source: Vurgaftman et al. 2001, Eq. (2.13)
            % https://doi.org/10.1063/1.1368156
            alpha = obj.param_T.Eg(1);
            beta = obj.param_T.Eg(2);
            Eg_0 = obj.param.Eg;
            eg = Eg_0 - (alpha * T^2) / (T + beta);
        end
        % Gets Kane parameter F
        function F = get_F(obj)
            F = obj.param.F;
        end
        % Gets the spin-orbit splitting in eV
        function dso = get_Dso(obj)
            dso = obj.param.Dso;
        end
        % Gets relative permittivity
        function eps = get_eps(obj)
            eps = obj.param.eps;
        end
        % Gets conduction band offset with respect to GaAs (eV).
        function ec = get_Ec0(obj)
            ec = obj.param.Ec;
        end
        % Gets elastic constant C11.
        function c11 = get_C11(obj, T)
            c11 = obj.calc_temp_depend(obj.param.c11, ...
                obj.param_T.c11, T);
        end
        % Gets elastic constant C12.
        function c12 = get_C12(obj, T)
            c12 = obj.calc_temp_depend(obj.param.c12, ...
                obj.param_T.c12, T);
        end
        % Gets conduction band deformation potential
        function ac = get_ac(obj, T)
            ac = obj.calc_temp_depend(obj.param.a_c, ...
                obj.param_T.a_c, T);
        end
        % Gets valence band deformation potential
        function av = get_av(obj)
            av = obj.param.a_v;
        end
        % Gets shear deformation potential
        function b = get_b(obj)
            b = obj.param.b;
        end
        % Gets lattice constant in dependence of temperature T (Angstrom).
        function a_latt = get_lattice_constant(obj, T)
            a_latt = obj.calc_temp_depend(obj.param.a_lattice, ...
                obj.param_T.a_lattice, T);
        end
        % Gets momentum matrix element in energy units (eV)
        function Ep = get_Ep(obj)
            Ep = obj.param.Ep;
        end
        function disp_param(obj, substrate, T, orientation)
            disp('Parameter: ');
            disp(obj.name(1:end));
            print_mat_par(obj, substrate, T, orientation);
        end
    end
    %
    methods (Static)
        function param = calc_temp_depend(x1, x2, T)
            param = x1 + x2 * T;
        end
    end
end
