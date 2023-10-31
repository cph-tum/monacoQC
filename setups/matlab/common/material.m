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

classdef material < matlab.mixin.Copyable
    %material Contains all required properties of a certain material.
    %
    properties (SetAccess = protected)
        name = '' % Material name
        conc = 0 % Concentration
        n_comp; % Number material components
        strain = true; % Parameter calculations including
        % strain dependencies (default: true)
        c_parab = 0 % Constant for nonparabolicity calculations
        % c_parab expresses the mass increasing factor of
        % in-plane to perpendicular mass (c_parab = (2*alpha+beta)/alpha)
        % Hendorfer et al. 1993, Equation 13
        % https://doi.org/10.1103/PhysRevB.48.232
    end
    %
    methods
        % Constructs material.
        function obj = material(name)
            obj.name = name;
        end
        % In-plane strain tensor component eps_II
        function eps_II = get_eps_II(obj, substrate, T)
            if (obj.strain == true)
                % relaxed lattice constant
                ar = obj.get_lattice_constant(T);
                % in plane lattice constant for fully strained material
                aII = substrate.get_lattice_constant(T);
                % epsilon in-plane
                eps_II = aII / ar - 1;
            else
                eps_II = 0;
            end
        end
        % In-growth strain tensor component eps_L
        function eps_L = get_eps_L(obj, substrate, T, orientation)
            % van der Walle 1989,
            % https://journals.aps.org/prb/abstract ...
            % /10.1103/PhysRevB.39.1871
            if (obj.strain == true)
                % Constant D depending on elastic constants
                if (strcmp(orientation, '001'))
                    D = 2 * obj.get_C12(T) / obj.get_C11(T);
                else
                    error(['Strain calculations ', 'for', ...
                        ' this interface orientation are not', ...
                        ' implemented, please use orientation 001.']);
                end
                % epsilon in-plane
                epsII = obj.get_eps_II(substrate, T);
                % relaxed lattice constant
                ar = obj.get_lattice_constant(T);
                % lattice constant perpendicular to plane
                aL = ar * (1 - D * epsII);
                % epsilon perpendicular
                eps_L = aL / ar - 1;
            else
                eps_L = 0;
            end
        end
        % Poisson ratio \nu
        function nu = get_nu(obj, T, orientation)
            % Hoke 2001, https://aip.scitation.org/doi/10.1063/1.1425954
            % Ioffe, accessed 04.2021
            if (strcmp(orientation, '001'))
                nu = obj.get_C12(T) / (obj.get_C12(T) + obj.get_C11(T));
            else
                error(['Strain calculations ', 'for', ...
                    ' this interface orientation are not ', ...
                    'implemented,  please use orientation 001.']);
            end
        end
        % Shear modulus G
        function G = get_G(obj, T, orientation)
            % van der Walle 1989,
            % https://journals.aps.org/prb/abstract ...
            % /10.1103/PhysRevB.39.1871
            % Constant D depending on elastic constants
            if (strcmp(orientation, '001'))
                D = 2 * obj.get_C12(T) / obj.get_C11(T);
            else
                error(['Strain calculations ', 'for', ...
                    ' this interface orientation are not ', ...
                    'implemented, please use orientation 001.']);
            end
            %
            G = 2 * (obj.get_C11(T) + 2 * obj.get_C12(T)) * (1 - D / 2);
        end
        % Calculates energy gap in eV
        function eg = get_Eg(obj, substrate, T, orientation)
            % Energy gap
            eg = obj.get_Eg0(T);
            % In-plane strain tensor component
            eII = obj.get_eps_II(substrate, T);
            % In-growth strain tensor component
            eL = obj.get_eps_L(substrate, T, orientation);
            % Conduction band shift
            Pec = obj.get_ac(T) * (2 * eII + eL);
            % Valence band shift component Pev
            Pev = obj.get_av() * (2 * eII + eL);
            eg = eg + Pec + Pev;
        end
        % Calculates conduction band edge energy Ec in eV
        function ec = get_Ec(obj, substrate, T, orientation)
            % Conduction band edge energy without strain
            ec = obj.get_Ec0();
            % epsilon in-plane
            epsII = obj.get_eps_II(substrate, T);
            % epsilon perpendicular
            epsL = obj.get_eps_L(substrate, T, orientation);
            % Fractional volume change due to strain
            dOm = 2 * epsII + epsL;
            % Conduction band shift
            ac = obj.get_ac(T);
            % Model solid theory in
            % van der Walle 1989,
            % https://journals.aps.org/prb/abstract ...
            % /10.1103/PhysRevB.39.1871
            ec = ec + ac * dOm;
        end
        % In-growth effective mass
        function me = get_meffL(obj, substrate, T, orientation, D)
            if (nargin ~= 5)
                D = 2 * obj.get_F;
            end
            me = [1, 0] * obj.calc_meff(substrate, T, orientation, D);
        end
        % In-plane effective mass
        function me = get_meffII(obj, substrate, T, orientation, D)
            if (nargin ~= 5)
                D = 2 * obj.get_F;
            end
            me = [0, 1] * obj.calc_meff(substrate, T, orientation, D);
        end
        % Gets material concentration in alloy.
        function c = get_conc(obj)
            c = obj.conc;
        end
        % Sets material concentration in alloy.
        function obj = set_conc(obj, value)
            obj.conc = value;
        end
        % Gets parabolicity factor in (eV)^-1.
        function par = get_parab(obj, substrate, T, orientation)
            par = calc_parab(obj, substrate, T, orientation);
        end
        % Set variable strain: Calculation of parameters with or without
        % strain
        function obj = set_strain(obj, value)
            obj.strain = value;
        end
        % Get variable strain: Calculation of parameters with or without
        % strain
        function s = get_strain(obj)
            s = obj.strain;
        end
        % Get critical layer thickness (Angstrom)
        function hcrit = get_hcrit(obj, substrate, T, orientation, N, id)
            % Returns the critical thickness
            % in Angstrom in iterative procedure
            if (nargin < 6)
                id = 'people';
                if (nargin < 5)
                    N = 10; % N: number of self-consistent iterations
                end
            end
            hcrit = obj.calc_hcrit(substrate, T, orientation, N, id);
        end
        %
        % Gets nonparabolicity parameter alpha (Eg + Dso/3)^-1.
        % Nelson 1987, https://doi.org/10.1103/PhysRevB.35.7770.
        % Sirtori 1994, https://doi.org/10.1103/PhysRevB.50.8663.
        % Eq. 7 of https://aip.scitation.org/doi/full/10.1063/1.4863665
        function a_p = get_1_alpha(obj, substrate, T, orientation)
            a_p = obj.get_Eg(substrate, T, orientation) + obj.get_Dso / 3;
        end
    end
    
    %
    methods (Access = protected)
        % Print material parameter
        function print_mat_par(obj, substrate, T, orientation)
            disp(['Ec = ', num2str(obj.get_Ec(substrate, T, ...
                orientation))]);
            disp(['Eg = ', num2str(obj.get_Eg(T))]);
            disp(['m_effL = ', num2str(obj.get_meffL(substrate, T, ...
                orientation))]);
            disp(['m_effII = ', num2str(obj.get_meffII(substrate, T, ...
                orientation))]);
            disp(['Parab = ', num2str(obj.get_parab(substrate, T, ...
                orientation))]);
            disp(['eps = ', num2str(obj.get_eps())]);
            disp(['lattice constant = ', ...
                num2str(obj.get_lattice_constant(T))]);
        end
        % Calculates the parabolicity factor.
        function par = calc_parab(obj, substrate, T, orientation)
            Eg = get_1_alpha(obj, substrate, T, orientation);
            par = 1 / Eg;
            par = par * obj.c_parab;
        end
        %
        function meff = calc_meff(obj, substrate, T, orientation, D)
            % Calculates in plane and in growth direction effective mass.
            % Sugawara 1993: https://doi.org/10.1103/PhysRevB.48.8102
            % Default value D = -6 for InGaAs.
            % Energy gap
            Eg = obj.get_Eg0(T);
            % In-plane strain tensor component
            eII = obj.get_eps_II(substrate, T);
            % In-growth strain tensor component
            eL = obj.get_eps_L(substrate, T, orientation);
            % Spin-orbit splitting
            dso = obj.get_Dso;
            % Conduction band shift
            Pec = obj.get_ac(T) * (2 * eII + eL);
            % Valence band shift component Pev
            Pev = obj.get_av() * (2 * eII + eL);
            % Valence band shift component Qev
            Qev = (-1) * obj.get_b * (eII - eL);
            % Calculation of strain parameters
            A = dso + Qev;
            %
            B = sqrt(dso^2+2*Qev*dso+9*Qev^2);
            %
            if (eII == 0)
                alpha = 1;
                beta = 0;
            else
                C = sqrt(2*B*(B - A));
                alpha = 2 * sqrt(2) * abs(Qev) / C;
                beta = (A - B) * abs(Qev) / (C * Qev);
            end
            % Energy-eigenstates of heavy-hole (HH), light hole (LH)
            % and split-off (SO) bands.
            Ehh = (-1) * Pev - Qev;
            Elh = (-1) * Pev + 0.5 * (Qev - dso + B);
            Eso = (-1) * Pev + 0.5 * (Qev - dso - B);
            % Momentum matrix element in energy units (eV)
            Ep = obj.get_Ep();
            % Calculation of inverse mass componentes
            % In-growth direction
            invmefL = (1 + D) + (Ep / 3) ...
                * ((sqrt(2) * alpha - beta)^2 / (Eg + Pec - Elh) ...
                +(sqrt(2) * beta + alpha)^2 / (Eg + Pec - Eso));
            % In-plane direction
            invmefII = (1 + D) + 1 / 2 * Ep * 1 / 3 * (3 ...
                / (Eg + Pec - Ehh) + (alpha - sqrt(2) * beta)^2 ...
                / (Eg + Pec - Elh) + (beta +sqrt(2) * alpha)^2 ...
                / (Eg + Pec - Eso));
            % Effective mass
            meff = [1 / invmefL; 1 / invmefII];
        end
        function hc = calc_hcrit(obj, substrate, T, orientation, N, id)
            % Calculates the critical thickness
            % in Angstrom in iterative procedure
            if (strcmp(obj.name, substrate.name))
                hc = inf;
                return
            end
            % Use model from People and Bean 1986, because of better
            % agreement with experimental results
            % P. J. Orders and B. F. Usher:
            % https://doi.org/10.1063/1.98004
            % (Werner Prost 1997, Springer)
            hc = obj.calc_hc(substrate, T, orientation, id);
            for i = 1:N
                hc = obj.calc_hc(substrate, T, orientation, id, hc);
            end
            
        end
        function hc = calc_hc(obj, substrate, T, orientation, id, hc)
            % Calculation of critical layer thickness (Angstrom)
            % People and Bean 1986: https://doi.org/10.1063/1.96206,
            % https://doi.org/10.1063/1.97637
            % Matthews and Blakeslee 1974
            % https://doi.org/10.1016/S0022-0248(74)80055-2
            % Lattice constants
            a = obj.get_lattice_constant(T);
            a_sub = substrate.get_lattice_constant(T);
            % Burgers vector
            b = a / sqrt(2);
            % Lattice mismatch
            f = abs(a-a_sub) / a_sub;
            % Poisson ratio
            nu = obj.get_nu(T, orientation);
            %
            if (nargin < 6)
                hc = 10;
            end
            % Critical thickness hc
            if (strcmp('people', id))
                hc = (1 - nu) * b^2 / f^2 / (1 + nu) ...
                    / (16 * pi * sqrt(2) * a) * log(hc/b);
            else
                hc = (1 - nu / 4) * b / (1 + nu) ...
                    / (4 * pi * f) * (log(hc/b) + 1);
            end
        end
    end
    
    % Calculation of material specific parameters
    methods (Abstract)
        % Gets name of the material
        name_mat = get_name(obj)
        % Gets band gap (Gamma valley) in eV.
        eg = get_Eg0(obj, T);
        % Gets relative permittivity.
        eps = get_eps(obj);
        % Gets the spin-orbit splitting in eV.
        dso = get_Dso(obj);
        % Gets conduction band offset with respect to GaAs (eV).
        ec = get_Ec0(obj);
        % Gets lattice constant in dependence of temperature T (Angstrom).
        a_latt = get_lattice_constant(obj, T);
        % Gets elastic constants.
        c11 = get_C11(obj, T);
        c12 = get_C12(obj, T);
        % Gets conduction band deformation potential
        ac = get_ac(obj, T);
        % Gets valence band deformation potential
        av = get_av(obj);
        % Gets shear deformation potential
        b = get_b(obj);
        % Gets Kane parameter F
        F = get_F(obj);
        % Gets momentum matrix element in energy units (eV)
        Ep = get_Ep(obj);
        % Displays material parameters
        disp_param(obj, substrate, T, orientation);
    end
end
