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

classdef AlloyRate < FermiGoldenRule
    % Calculates intersubband transition rates due to scattering of
    % electrons in ternary alloys using Fermi's golden rule.

    properties
        layer_conc % vector: Material concentration x in ternary alloys.
        lattice_const % vector: Lattice constant of each layer [m].
        dV_alloy % vector: Alloy scattering potential for each layer [J].
    end

    methods
        function obj = AlloyRate(eigen, device, scenario, options)
            % Constructs an object of type AlloyRate.
            %
            % Syntax:
            %   obj = AlloyRate(eigen, device, scenario)
            %   obj = AlloyRate(eigen, device, scenario, options)
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
            obj.Name = 'alloy scattering';

            % get physical properties of each layer
            for i = 1:length(device.layers)
                dV = device.layers{i}.material.delta_V;
                if dV == 0 && ~isa(device.layers{i}.material, 'binary')
                    warning(['Alloy potential is zero ', ...
                        'for certain layer. Assuming ', ...
                        'material parameter is not set.']);
                end
                if isa(device.layers{i}.material, 'ternary')
                    obj.layer_conc(i) = device.layers{i}.material.conc;
                else
                    obj.layer_conc(i) = 0;
                end

                % convert from A to m
                obj.lattice_const(i) = device.layers{i}.material.a_lattice * 1e-10;
                obj.dV_alloy(i) = dV * phys_const.e0;
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
            %   k (scalar): Wavevector of inital state.
            %
            % Output Arguments:
            %   W_ikj (scalar): Transition rate.

            % energy of electron in subband i with in-plane wave vector k
            E_ik = obj.E(i) + (phys_const.hbar * k)^2 / (2 * obj.mEff(i));

            if obj.E(j) < E_ik
                % Get lattice const, concentration, alloy potential between
                % interfaces
                latticeConst = zeros(1, length(obj.z));
                co = zeros(1, length(obj.z));
                dV = zeros(1, length(obj.z));
                PosI = obj.PosInterf(2:end);

                c = 1;
                for n = 1:length(obj.z)
                    if c < length(PosI) && obj.z(n) >= PosI(c)
                        c = c + 1;
                    end
                    latticeConst(n) = obj.lattice_const(c);
                    co(n) = obj.layer_conc(c);
                    dV(n) = obj.dV_alloy(c);
                end

                Omega0 = latticeConst.^3 / 4; % volume of crystal unit cell
                f = obj.mEff(j) / phys_const.hbar^3 * Omega0 .* ...
                    co .* (1 - co) .* dV.^2 .* (abs(obj.psi(:, i) ...
                    .*obj.psi(:, j)))'.^2;
                W_ikj = trapz(obj.z*1e-10, f);

            else
                % no scattering for E_j > E_ik
                W_ikj = 0;
            end
        end
    end
end
