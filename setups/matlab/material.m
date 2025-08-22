%
% monacoQC: An object-oriented Matlab-based device engineering tool for
% quantum cascade structures.
%
% Copyright (C) 2025, Computational Photonics Group, Technical University of
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

    properties (Abstract, SetAccess = public)
        temp{mustBeNonnegative}
    end
    properties (Abstract, Dependent)
        name_III
        name_V
    end
    properties (SetAccess = public)
        strain = true
    end
    properties (SetAccess = protected)
        name
    end

    % setter methods for public properties
    methods
        function set.strain(obj, value)
            if islogical(value)
                obj.strain = value;
            else
                error("Strain value must be boolean!")
            end
        end
    end

    methods
        % In-plane strain tensor component eps_II
        function eps_II = get_eps_II(obj, substrate, ~)
            if (obj.strain == true)
                % epsilon in-plane
                eps_II = substrate.a_lattice / obj.a_lattice - 1;
            else
                eps_II = 0;
            end
        end

        % In-growth strain tensor component eps_L
        function eps_L = get_eps_L(obj, substrate, orientation)
            % van der Walle 1989,
            % https://journals.aps.org/prb/abstract/10.1103/PhysRevB.39.1871
            if (obj.strain == true)
                % Constant D depending on elastic constants
                if (strcmp(orientation, '001'))
                    D = 2 * obj.C12 / obj.C11;
                else
                    error(['Strain calculations ', 'for', ...
                        ' this interface orientation are not', ...
                        ' implemented, please use orientation 001.']);
                end
                % epsilon in-plane
                epsII = obj.get_eps_II(substrate);
                % relaxed lattice constant
                ar = obj.a_lattice;
                % lattice constant perpendicular to plane
                aL = ar * (1 - D * epsII);
                % epsilon perpendicular
                eps_L = aL / ar - 1;
            else
                eps_L = 0;
            end
        end

        % Poisson ratio \nu
        function nu = get_nu(obj, ~, orientation)
            % Hoke 2001, https://aip.scitation.org/doi/10.1063/1.1425954
            % Ioffe, accessed 04.2021
            if (strcmp(orientation, '001'))
                nu = obj.C12 / (obj.C12 + obj.C11);
            else
                error(['Strain calculations ', 'for', ...
                    ' this interface orientation are not ', ...
                    'implemented,  please use orientation 001.']);
            end
        end

        % Shear modulus G
        function G = get_G(obj, ~, orientation)
            % van der Walle 1989,
            % https://journals.aps.org/prb/abstract/10.1103/PhysRevB.39.1871
            % Constant D depending on elastic constants
            if (strcmp(orientation, '001'))
                D = 2 * obj.C12 / obj.C11;
            else
                error(['Strain calculations ', 'for', ...
                    ' this interface orientation are not ', ...
                    'implemented, please use orientation 001.']);
            end

            G = 2 * (obj.C11 + 2 * obj.C12) * (1 - D / 2);
        end

        % Calculates energy gap in eV
        function eg = get_Eg(obj, substrate, orientation)
            % Energy gap
            eg = obj.Eg;
            % In-plane strain tensor component
            eII = obj.get_eps_II(substrate);
            % In-growth strain tensor component
            eL = obj.get_eps_L(substrate, orientation);
            % Conduction band shift
            Pec = obj.a_cb * (2 * eII + eL);
            % Valence band shift component Pev
            Pev = obj.a_vb * (2 * eII + eL);
            eg = eg + Pec + Pev;
        end

        % Calculates conduction band edge energy Ec in eV
        function ec = get_Ec(obj, substrate, orientation)
            % Conduction band edge energy without strain
            ec = obj.Ec;
            % epsilon in-plane
            epsII = obj.get_eps_II(substrate);
            % epsilon perpendicular
            epsL = obj.get_eps_L(substrate, orientation);
            % Fractional volume change due to strain
            dOm = 2 * epsII + epsL;
            % Conduction band shift
            ac = obj.a_cb;
            % Model solid theory in van der Walle 1989,
            % https://journals.aps.org/prb/abstract ...
            % /10.1103/PhysRevB.39.1871
            ec = ec + ac * dOm;
        end

        % In-growth effective mass
        function me = get_meffL(obj, substrate, orientation, D)
            if (nargin ~= 4)
                D = 2 * obj.F;
            end
            me = [1, 0] * obj.calc_meff(substrate, orientation, D);
        end

        % In-plane effective mass
        function me = get_meffII(obj, substrate, orientation, D)
            if (nargin ~= 4)
                D = 2 * obj.F;
            end
            me = [0, 1] * obj.calc_meff(substrate, orientation, D);
        end

        % Gets parabolicity factor in (eV)^-1.
        function par = get_parab(obj, substrate, orientation)
            par = calc_parab(obj, substrate, orientation);
        end

        % Get critical layer thickness (Angstrom)
        function hcrit = get_hcrit(obj, substrate, orientation, N, id)
            % Returns the critical thickness
            % in Angstrom in iterative procedure
            if (nargin < 5)
                id = 'people';
                if (nargin < 4)
                    N = 10; % N: number of self-consistent iterations
                end
            end
            hcrit = obj.calc_hcrit(substrate, orientation, N, id);
        end

        % Gets nonparabolicity parameter alpha (Eg + Dso/3)^-1.
        % Nelson 1987, https://doi.org/10.1103/PhysRevB.35.7770.
        % Sirtori 1994, https://doi.org/10.1103/PhysRevB.50.8663.
        % Eq. 7 of https://aip.scitation.org/doi/full/10.1063/1.4863665
        function a_p = get_1_alpha(obj, substrate, orientation)
            a_p = obj.get_Eg(substrate, orientation) + obj.delta_so / 3;
        end
    end

    methods (Access = protected)
        % Calculates the parabolicity factor.
        function par = calc_parab(obj, substrate, orientation)
            par = obj.c_parab / get_1_alpha(obj, substrate, orientation);
        end

        function meff = calc_meff(obj, substrate, orientation, D)
            % Calculates in plane and in growth direction effective mass.
            % Sugawara 1993: https://doi.org/10.1103/PhysRevB.48.8102
            % Default value D = -6 for InGaAs.
            % Energy gap
            Eg = obj.Eg;
            % In-plane strain tensor component
            eII = obj.get_eps_II(substrate);
            % In-growth strain tensor component
            eL = obj.get_eps_L(substrate, orientation);
            % Spin-orbit splitting
            dso = obj.delta_so;
            % Conduction band shift
            Pec = obj.a_cb * (2 * eII + eL);
            % Valence band shift component Pev
            Pev = obj.a_vb * (2 * eII + eL);
            % Valence band shift component Qev
            Qev = (-1) * obj.b * (eII - eL);
            % Calculation of strain parameters
            A = dso + Qev;

            B = sqrt(dso^2+2*Qev*dso+9*Qev^2);

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
            Ep = obj.Ep;
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

        function hc = calc_hcrit(obj, substrate, orientation, N, id)
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
            hc = obj.calc_hc(substrate, orientation, id);
            for i = 1:N
                hc = obj.calc_hc(substrate, orientation, id, hc);
            end
        end

        function hc = calc_hc(obj, substrate, orientation, id, hc)
            % Calculation of critical layer thickness (Angstrom)
            % People and Bean 1986: https://doi.org/10.1063/1.96206,
            % https://doi.org/10.1063/1.97637
            % Matthews and Blakeslee 1974
            % https://doi.org/10.1016/S0022-0248(74)80055-2
            % Lattice constants
            a = obj.a_lattice;
            a_sub = substrate.a_lattice;
            % Burgers vector
            b = a / sqrt(2);
            % Lattice mismatch
            f = abs(a-a_sub) / a_sub;
            % Poisson ratio
            nu = obj.get_nu(substrate, orientation);

            if (nargin < 5)
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
end
