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

classdef ternary < material
    %ternary Contains all required properties of a certain ternary
    %        material A_(x)B_(1-x)C.
    %
    properties (SetAccess = protected)
        b1 % Binary BC
        b2 % Binary AC
        C % Bowing factors
    end
    %
    methods
        % Constructs ternary material.
        function obj = ternary(name)
            obj = obj@material(name);
            obj.n_comp = 3;
        end
        %
        % Gets conduction band offset with respect to GaAs (eV).
        function ec = get_Ec0(obj)
            i1 = obj.b1.get_Ec0;
            i2 = obj.b2.get_Ec0;
            ec = obj.calc_param(obj.conc, obj.C.Ec, i1, i2);
        end
        % Calculates the effective mass.
        function me = interp_meff(obj, T, orientation)
            i1 = obj.b1.get_meffL(obj.b1, T, orientation);
            i2 = obj.b2.get_meffL(obj.b2, T, orientation);
            me = obj.calc_param(obj.conc, obj.C.m_eff, i1, i2);
        end
        % Gets the spin-orbit splitting in eV.
        function dso = get_Dso(obj)
            i1 = obj.b1.get_Dso;
            i2 = obj.b2.get_Dso;
            dso = obj.calc_param(obj.conc, obj.C.Dso, i1, i2);
        end
        %
        % Gets relative permittivity.
        function eps = get_eps(obj)
            i1 = obj.b1.get_eps;
            i2 = obj.b2.get_eps;
            eps = obj.calc_param(obj.conc, obj.C.eps, i1, i2);
        end
        %
        % Calculates the band gap energy without strain.
        function eg = get_Eg0(obj, T)
            i1 = obj.b1.get_Eg0(T);
            i2 = obj.b2.get_Eg0(T);
            eg = obj.calc_param(obj.conc, obj.C.Eg, i1, i2);
        end
        % Calculates elastic constant c11.
        function c11 = get_C11(obj, T)
            i1 = obj.b1.get_C11(T) * obj.b1.get_lattice_constant(T);
            i2 = obj.b2.get_C11(T) * obj.b2.get_lattice_constant(T);
            a = obj.get_lattice_constant(T);
            c11 = obj.calc_param(obj.conc, 0, i1, i2);
            c11 = c11 / a;
        end
        % Calculates elastic constant c12.
        function c12 = get_C12(obj, T)
            i1 = obj.b1.get_C12(T) * obj.b1.get_lattice_constant(T);
            i2 = obj.b2.get_C12(T) * obj.b2.get_lattice_constant(T);
            a = obj.get_lattice_constant(T);
            c12 = obj.calc_param(obj.conc, 0, i1, i2);
            c12 = c12 / a;
        end
        
        % Calculates deformation potential of conduction band.
        function ac = get_ac(obj, T)
            i1 = obj.b1.get_ac(T);
            i2 = obj.b2.get_ac(T);
            ac = obj.calc_param(obj.conc, obj.C.a_c, i1, i2);
        end
        % Calculates deformation potential of valence band.
        function av = get_av(obj)
            i1 = obj.b1.get_av;
            i2 = obj.b2.get_av;
            av = obj.calc_param(obj.conc, obj.C.a_v, i1, i2);
        end
        function b = get_b(obj)
            i1 = obj.b1.get_b;
            i2 = obj.b2.get_b;
            b = obj.calc_param(obj.conc, 0, i1, i2);
        end
        %
        function F = get_F(obj)
            i1 = obj.b1.get_F;
            i2 = obj.b2.get_F;
            F = obj.calc_param(obj.conc, obj.C.F, i1, i2);
        end
        function Ep = get_Ep(obj)
            i1 = obj.b1.get_Ep;
            i2 = obj.b2.get_Ep;
            Ep = obj.calc_param(obj.conc, obj.C.Ep, i1, i2);
        end
        function a_latt = get_lattice_constant(obj, T)
            a1 = obj.b1.get_lattice_constant(T);
            a2 = obj.b2.get_lattice_constant(T);
            a_latt = obj.calc_param(obj.conc, 0, a1, a2);
        end
        function disp_param(obj, substrate, T, orientation)
            disp('Parameter: ');
            disp(strcat(obj.name(1:2), num2str(obj.get_conc), ...
                obj.name(3:4), num2str(1-obj.get_conc), obj.name(5:6)));
            print_mat_par(obj, substrate, T, orientation);
        end
        function name_mat = get_name(obj)
            n_III = {obj.b2.name(1:2), obj.b1.name(1:2)};
            n_V = {obj.b2.name(3:end), obj.b1.name(3:end)};
            if (isequal(n_III{1}, n_III{2}))
                name_mat = strcat(n_III{1}, n_V{1}, '_{', num2str(obj.conc), ...
                    '}', n_V{2}, '_{', num2str(1-obj.conc), '}');
            elseif (isequal(n_V{1}, n_V{2}))
                name_mat = strcat(n_III{1}, '_{', num2str(obj.conc), '}', ...
                    n_III{2}, '_{', num2str(1-obj.conc), '}', n_V{1});
            else
                error('Cannot generate proper material name!');
            end
        end
        
    end
    
    %
    methods (Static)
        % Calculates ternary material parameter including bowing factors.
        function param = calc_param(x, C, input1, input2)
            param = [input1, C, input2];
            % Source: Vurgaftman et al. 2001, Eq. (4.1)
            X = [1 - x; (-x * (1 - x)); x];
            param = param * X;
        end
    end
end
