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

classdef AcousticPhononRate < FermiGoldenRule
    % Calculates intersubband transition rates due to scattering of
    % electrons by acoustic phonons using Fermi's golden rule.

    properties
        density % vector: Density of each layer [kg/m^3].
        v_speed % vector: Speed of sound for each layer [m/s].
        Xi % vector: Deformation potential for each layer [J].
    end

    methods
        function obj = AcousticPhononRate(eigen, device, scenario, options)
            % Constructs an object of type AcousticPhononRate.
            %
            % Syntax:
            %   obj = AcousticPhononRate(eigen, device, scenario)
            %   obj = AcousticPhononRate(eigen, device, scenario, options)
            %
            % Input Arguments:
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.
            %   device (device-object): Contains information about the
            %     structure/ geometry and materials of the QCL.
            %   scenario (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %   options (cell-array): Array containing name-value pairs for
            %     changing values of default properties. Valid names are
            %     ``num_k``, ``Te`` and ``screening_model``.

            if nargin < 4
                options = {};
            else
                validStrings = ["num_k", "Te", "screening_model"];
                for l = 1:2:length(options)
                    validatestring(options{l}, validStrings);
                end
            end

            obj = obj@FermiGoldenRule(eigen, device, scenario, options);
            obj.Name = 'acoustic phonon scattering';

            % get physical properties of each layer
            for j = 1:length(device.layers)
                obj.density(j) = device.layers{j}.material.rho * 1e3;
                obj.v_speed(j) = device.layers{j}.material.v_sound;
                obj.Xi(j) = device.layers{j}.material.Xi * phys_const.e0;
            end
        end

        function W_ikj = calc_state_rate(obj, i, j, k)
            % Calculates transition rate from state ik to subband j
            % according to https://doi.org/10.1063/1.4863665.
            %
            % Syntax:
            %   W_ikj = calc_state_rate(obj, i, j, k)
            %
            % Input Arguments:
            %   i (scalar): Inital subband.
            %   j (scalar): Final subband.
            %   k (scalar): Wavevector of the inital state.
            %
            % Output Arguments:
            %   W_ikj (scalar): Transition rate.

            % energy of electron in subband i with in-plane wave vector k
            E_ik = obj.E(i) + (phys_const.hbar * k)^2 / (2 * obj.mEff(i));

            if obj.E(j) < E_ik
                % get properties between interfaces
                rho = zeros(1, length(obj.z));
                vs = zeros(1, length(obj.z));
                defPot = zeros(1, length(obj.z));
                PosI = obj.PosInterf(2:end);
                c = 1;
                for n = 1:length(obj.z)
                    if c < length(PosI) && obj.z(n) >= PosI(c)
                        c = c + 1;
                    end
                    rho(n) = obj.density(c);
                    vs(n) = obj.v_speed(c);
                    defPot(n) = obj.Xi(c);
                end

                f = defPot.^2 ./ (rho .* vs.^2) .* ...
                    abs((obj.psi(:, i) .* obj.psi(:, j)))'.^2;

                W_ikj = phys_const.kB * obj.T / phys_const.hbar^3 * ...
                    obj.mEff(j) * trapz(obj.z*1e-10, f);

            else
                % no scattering for E_j > E_ik
                W_ikj = 0;
            end
        end
    end
end
