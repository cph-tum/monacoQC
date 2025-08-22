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

classdef LOPhononRate < FermiGoldenRule
    % Calculates intersubband transition rates due to scattering of
    % electrons by longitudinal optical phonons using Fermi's golden rule.

    properties
        epsr_stat % vector: Static relative permitivity at each grid point.
        epsr_inf % vector: High freq rel permitivity at each grid point.
        E_phLO % vector: LO Phonon energy for each layer [J].
        layer_conc % vector: Material concentration for each layer.
        em_ab % scalar: ``em_ab=1`` (emission), ``em_ab=-1`` (absorption).
        psi2_ft % 3-d array: Fourier transform of product of wavefunctions.
    end

    methods
        function obj = LOPhononRate(eigen, device, ...
                scenario, em_ab, options)
            % Constructs an object of type LOPhononRate.
            %
            % Syntax:
            %   obj = LOPhononRate(eigen, device, scenario, em_ab)
            %   obj = LOPhononRate(eigen, device, scenario, em_ab, options)
            %
            % Input:
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.
            %   device (device-object): Contains information about the
            %     structure/ geometry and materials of the QCL.
            %   scenario (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %   em_ab (scalar): Set absorption or emission mechanism.
            %     Valid choices: ``absorption`` or ``emission``.
            %   options (cell-array): Array containing name-value pairs for
            %     changing values of default properties. Valid names are
            %     ``num_k``, ``Te`` and ``screening_model``.

            if nargin < 5
                options = {};
            else
                validStrings = ["num_k", "Te", "screening_model"];
                for l = 1:2:length(options)
                    validatestring(options{l}, validStrings);
                end
            end

            obj = obj@FermiGoldenRule(eigen, device, scenario, options);

            if strcmp(em_ab, 'emission')
                obj.Name = 'LO phonon emission';
                obj.em_ab = 1;
            elseif strcmp(em_ab, 'absorption')
                obj.Name = 'LO phonon absorption';
                obj.em_ab = -1;
            else
                error("Rate needs to be classified as"+ ...
                    "'emission' or 'absorption'");
            end

            % get physical properties of each layer
            for i = 1:length(device.layers)
                epsr_stat(i) = device.layers{i}.material.eps_r;
                epsr_inf(i) = device.layers{i}.material.eps_r_inf;
                PhLO = device.layers{i}.material.E_lo;
                obj.E_phLO(:, i) = PhLO * phys_const.e0 * 1e-3;
                if isa(device.layers{i}.material, 'ternary')
                    obj.layer_conc(i) = device.layers{i}.material.conc;
                else
                    obj.layer_conc(i) = 0;
                end
            end

            % get physical properties for each grid point
            obj.epsr_stat = zeros(1, length(obj.z));
            obj.epsr_inf = zeros(1, length(obj.z));
            PosI = obj.PosInterf(2:end);
            c = 1;
            for n = 1:length(obj.z)
                if c < length(PosI) && obj.z(n) >= PosI(c)
                    c = c + 1;
                end
                obj.epsr_stat(n) = epsr_stat(c);
                obj.epsr_inf(n) = epsr_inf(c);
            end

            % calculate screened wave vector squared qs2
            obj.qs2 = obj.calc_qs2();

            % calculate FFT of product of wavefunctions psi2_ft
            obj.set_psi2_ft();
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

            % grid points for integration
            x_i = linspace(0, pi, 100);

            % get phonon energy of first well
            co = [obj.layer_conc(2), 1 - obj.layer_conc(2)];
            Eph = obj.E_phLO(:, 2);
            if Eph(1) == Eph(2)
                Eph = max(Eph);
                co = 1;
            end

            % Bose-Einstein distribution for the the
            % phonon occupation number
            Nph = 1 ./ (exp(Eph/(phys_const.kB * obj.T)) - 1);
            % Energy of state |ik>
            E_ik = obj.E(i) + (phys_const.hbar * k)^2 / (2 * obj.mEff(i));

            % Final wavevector squared
            kPrime2 = k^2 * obj.mEff(j) / obj.mEff(i) + ...
                2 * obj.mEff(j) * (obj.E(i) - obj.E(j) - ...
                1 * obj.em_ab * Eph) / phys_const.hbar^2;
            W_ikj = 0;

            for m = 1:length(Eph)
                if obj.E(j) <= (E_ik - obj.em_ab * ...
                        Eph(m)) && kPrime2(m) >= 0
                    % shift in-plane phonon wave vector by screening
                    % wave vector
                    q = sqrt(k^2+kPrime2(m)-2*k*sqrt( ...
                        kPrime2(m))*cos(x_i)+obj.qs2);

                    r = obj.mEff(j) * phys_const.e0^2 / ...
                        (4 * pi^2 * phys_const.hbar^2 * phys_const.eps0);

                    % multiply with concentration co to consider both Eph
                    % (in ternary alloy)
                    r = r * co(m) * Eph(m) / phys_const.hbar * ...
                        (Nph(m) + 0.5 + 0.5 * obj.em_ab);

                    zPos = obj.z * 1e-10;
                    delta_z = abs(zPos(end)-zPos(end-1));

                    J = calc_formfactor_fft(obj, i, j, q, delta_z) * ...
                        mean((1 ./ obj.epsr_inf - 1 ./ obj.epsr_stat)) * ...
                        pi ./ q;

                    W_ikj = W_ikj + r * trapz(x_i, J);
                end
            end
        end

        function integral = calc_formfactor_fft(obj, i, j, q, delta_z)
            % Calculates formfactor using Fourier transform method
            % described in https://doi.org/10.1063/5.0041392.
            %
            % Syntax:
            %   integral = calc_form_factor_fft(obj, i, j, q, delta_z)
            %
            % Input Arguments:
            %   i (scalar): Initial subband.
            %   j (scalar): Final subband.
            %   q (vector): Exchanged wavevector [1/m].
            %   delta_z (scalar): Grid spacing [m].
            %
            % Output Arguments:
            %   integral (vector): Formfactor for the transition from
            %     subband i to subband j in dependence on q.

            integral = zeros(size(q));
            psi_prod_ft = reshape(obj.psi2_ft(i, j, :), 1, []);

            N = length(psi_prod_ft);
            N_f = (N - 1) / 2;
            c = -N_f:N_f;

            prod_psi = psi_prod_ft .* flip(psi_prod_ft);
            two_pi_c2 = (2 * pi * c).^2;

            for m = 1:length(q)
                denominator1 = (q(m) * N * delta_z)^2 + two_pi_c2;
                nominator1 = 2 * q(m) * (1 - obj.qs2 / (2 * q(m).^2));
                denominator2 = ((q(m) * N * delta_z)^2 + two_pi_c2).^2;
                nominator2 = -obj.qs2 / q(m) * ((q(m) * ...
                    N * delta_z)^2 - two_pi_c2);

                integral(m) = N * delta_z^3 * real( ...
                    sum(prod_psi.*(nominator1 ./ ...
                    denominator1 + nominator2 ./ denominator2)));
            end
            integral(q == 0) = 0;
        end

        function integral = calc_formfactor(obj, i, j, q)
            % Calculates the formfactor using normal integration method
            % to solve the double integral.
            %
            % Syntax:
            %   integral = calc_form_factor_fft(obj, i, j, q)
            %
            % Input Arguments:
            %   i (scalar): Initial subband.
            %   j (scalar): Final subband.
            %   q (vector): Exchanged wavevector [1/m].
            %
            % Output Arguments:
            %   integral (vector): Formfactor J_ij(q) for the transition
            %     from subband i to subband j in dependence on q.


            integral = zeros(size(q));
            zPos = obj.z * 1e-10;
            psi1 = obj.psi(:, i) .* conj(obj.psi(:, j));
            psi2 = obj.psi(:, j) .* conj(obj.psi(:, i));

            for n = 1:length(q)
                % constants
                c1 = obj.qs2 ./ (2 * q(n)^2);
                c2 = obj.qs2 ./ (2 * q(n));

                % z' > z
                I = (1 - c1 - c2 .* zPos) .* exp(-q(n)*zPos) .* ...
                    cumtrapz(zPos, psi1.*exp(q(n)*zPos));
                I = I + exp(-q(n)*zPos) .* c2 .* ...
                    cumtrapz(zPos, zPos.*psi1.*exp(q(n)*zPos));
                % z' < z
                I = I + (1 - c1 + c2 .* zPos) .* exp(q(n)*zPos) .* ...
                    flip(cumtrapz(zPos, flip(psi1.*exp(-q(n)*zPos))));
                I = I - exp(q(n)*zPos) .* c2 .* flip( ...
                    cumtrapz(zPos, flip(zPos.*psi1.*exp(-q(n)*zPos))));
                % Integration over dz'
                integral(n) = trapz(zPos, I.*psi2);
            end
            integral(q == 0) = 0;
        end

        function E = get_energy_pauli_blocking(obj)
            % Calculates the energy of the absorbed or emitted phonon.
            %
            % Syntax:
            %   E = get_energy_pauli_blocking(obj)
            %
            % Output Arguments:
            %   E (vector): Returns the phonon energy for the quantum well
            %     and barrier material. ``E<0`` corresponds to the emission
            %     and ``E>0`` to the absorption of a phonon. For  ternary
            %     alloys, E is linearly interpolated according to the mole
            %     fraction x.

            co = [obj.layer_conc(2), 1 - obj.layer_conc(2)];
            Eph = obj.E_phLO(:, 2);
            if Eph(1) == Eph(2)
                % binary alloys
                E = Eph(1);
            else
                % ternary alloys: use linear
                % interpolation of phonon energy
                E = Eph(1) * co(1) + Eph(2) * co(2);
            end
            E = -obj.em_ab .* E;
        end

        function set_eigenstates(obj, eigen)
            % Sets new eigenenergies, effective masses and wavefunctions.
            %
            % Syntax:
            %   set_eigenstates(obj, eigen)
            %
            % Input Arguments:
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.

            set_eigenstates@ScatteringRate(obj, eigen);
            obj.set_psi2_ft();
        end

        function update(obj, ns)
            % Updates the sheet densities, quasi Fermi-levels and screening
            % wavevector.
            %
            % Syntax:
            %   update(obj, ns)
            %
            % Input:
            %   ns (vector): Sheet density for four QCL-periods [1/m^2].

            update@FermiGoldenRule(obj, ns);
            % Calculate new screened wave vector qs2
            obj.qs2 = obj.calc_qs2();
        end

        function psi2_ft = set_psi2_ft(obj)
            % Calculates Fourier transform of product of wavefunctions
            % for fast calculation of the form factor as described in
            % https://doi.org/10.1063/5.0041392.
            %
            % Syntax:
            %   psi2_ft = set_psi2_ft(obj)
            %
            % Output:
            %   psi2_ft (3-d array): Returns the Fourier transforms of
            %     products of two wavefunctions in dependence on the
            %     initial subband (1st dimension), final subband (2nd
            %     dimension) and k-space variable of the Fourier transform
            %     (3rd dimension).

            % collect products of all wavefunctions
            psi2 = NaN(size(obj.psi, 2), ...
                size(obj.psi, 2), size(obj.psi, 1));

            for i = 1:size(obj.psi, 2)
                for j = 1:size(obj.psi, 2)
                    psi2(i, j, :) = obj.psi(:, i) .* obj.psi(:, j);
                end
            end

            % need uneven number of sampling points for special FFT form
            N = size(psi2, 3);
            if rem(N, 2) == 0
                psi2 = psi2(:, :, 1:end-1);
                N = size(psi2, 3);
            end

            % transform matlab convention of FFT to convention of
            % https://doi.org/10.1063/5.0041392
            N_f = (N - 1) / 2;
            m = 1:N;
            phase_factor = fftshift(exp(1i*2*pi*N_f*(m - 1)/N)');
            psi2_ft = fftshift(fft(psi2, [], 3), 3);
            for i = 1:size(obj.psi, 2)
                for j = 1:size(obj.psi, 2)
                    psi2_ft(i, j, :) = reshape(psi2_ft(i, j, :), ...
                        [], 1) .* phase_factor;
                end
            end
            obj.psi2_ft = psi2_ft;
        end
    end
end
