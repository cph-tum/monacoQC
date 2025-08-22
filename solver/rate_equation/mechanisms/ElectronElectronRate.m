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

classdef ElectronElectronRate < FermiGoldenRule
    % Calculates intersubband transition rates due to collision of two
    % electrons using Fermi's golden rule.

    properties
        epsr_stat % vector: Static relative permitivity for each layer.
        psi2_ft % 3-d array: Fourier transform of products of the wavefunctions.
        formf % 5-d array: Pre calculated form factors.
        qval % vector: q values used for pre-calculating the form factor.
        theta % vector: Grid for phase space integration.
        W_ikfk_jg % 6-d array: K-resolved transition rates.
    end

    methods
        function obj = ElectronElectronRate(eigen, ...
                device, scenario, options)
            % Constructs an object of type ElectronElectronRate.
            %
            % Syntax:
            %   obj = ElectronElectronRate(eigen, device, scenario)
            %   obj = ElectronElectronRate(eigen, device, scenario, options)
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
            %     ``num_k``, ``Te``, ``num_theta``, ``screening_model`` and
            %     ``num_interp`` (number of interpolation points
            %     for the form factor).

            if nargin < 4
                options = {};
            else
                validStrings = ["num_k", "Te", "num_theta", ...
                    "num_interp", "screening_model"];
                for l = 1:2:length(options)
                    validatestring(options{l}, validStrings);
                end
            end

            obj = obj@FermiGoldenRule(eigen, device, scenario, options);
            obj.Name = 'electron-electron scattering';

            % get physical properties of each layer
            for i = 1:length(device.layers)
                obj.epsr_stat(i) = device.layers{i}.material.eps_r;
            end

            % calculate screened wave vector qs2
            obj.qs2 = obj.calc_qs2();
            % calculate FFT of product of wavefunctions psi2_ft
            obj.set_psi2_ft();

            % default values
            nwf = obj.num_states; % number of wavefunctions
            n_alpha_theta = 10; % points for phase space integration
            nq_ff = 200; % points for pre-calulating the form factor

            obj.formf = zeros(nq_ff, nwf, nwf, nwf, nwf);
            obj.qval = zeros(1, nq_ff);
            obj.theta = linspace(0, pi, n_alpha_theta);

            % overwrite default values by user defined options
            for i = 1:2:length(options)
                key = options{i};
                value = options{i+1};
                if key == "num_theta"
                    obj.theta = linspace(0, pi, value);
                elseif key == "num_interp"
                    obj.qval = zeros(1, value);
                    obj.formf = zeros(value, nwf, nwf, nwf, nwf);
                end
            end
            obj.W_ikfk_jg = zeros(length(obj.k), ...
                length(obj.k), nwf, nwf, nwf, nwf);
        end

        function tau_inv = calculate(obj, calc_W)
            % Calculates electron-electron scattering rates between all
            % subbands.
            %
            % Syntax:
            %   tau_inv = calculate(obj, calc_W)
            %
            % Input Arguments:
            %   calc_W (logical): Flag determining if k-resolved transition
            %     rates have to be re-calculated.
            %
            % Output Arguments:
            %   tau_inv (matrix): Matrix containing the transition rates
            %     between the subbands of the two central periods.

            if nargin == 1
                calc_W = 1;
            end
            nwf = obj.num_states;
            if calc_W
                obj.set_formf
                obj.W_ikfk_jg = zeros(length(obj.k), ...
                    length(obj.k), nwf, nwf, nwf, nwf);
            end
            obj.W = zeros(nwf, nwf, length(obj.k));
            tau_inv = zeros(nwf, nwf);

            % loop over all states and calculate transition rates
            for ii = 1:nwf
                for jj = 1:nwf
                    if ii ~= jj
                        W = zeros(1, length(obj.k));
                        for ff = 1:nwf
                            for gg = 1:nwf
                                % the transition rate W_ijji does not
                                % result in a net change in  the
                                % carrier density
                                if ~(ii == gg && jj == ff)
                                    if calc_W
                                        % calculate k-resolved rates
                                        W_ikfk_jg = obj.calc_state_rate( ...
                                            ii, ff, jj, gg);
                                        if ff == ii && gg == jj
                                            W_ikfk_jg = 2 * W_ikfk_jg;
                                        end
                                        obj.W_ikfk_jg(:, :, ii, ...
                                            ff, jj, gg) = W_ikfk_jg;
                                    end
                                    f = obj.indices_i_j(ff);
                                    i = obj.indices_i_j(ii);
                                    W_kk = reshape(obj.W_ikfk_jg( ...
                                        :, :, ii, ff, jj, gg), length( ...
                                        obj.k), length(obj.k));
                                    fermi_f = FermiGoldenRule. ...
                                        fermi_dirac(obj.E(f)+ ...
                                        phys_const.hbar^2* ...
                                        obj.k.^2./(2 * obj.mEff(i)), ...
                                        obj.E_F(f), obj.T_e(f));
                                    W_ikf_jg = trapz(obj.k, W_kk.* ...
                                        reshape(fermi_f.*obj.k, [], 1));
                                    % W_iijj contributes twice to the
                                    % total scattering rates, because
                                    % two electrons are transferred
                                    % between the subbands
                                    W = W + W_ikf_jg;
                                end
                            end
                        end
                        obj.W(ii, jj, :) = W; % update
                        % calculate k-averaged transition rates
                        i = obj.indices_i_j(ii);
                        j = obj.indices_i_j(jj);
                        tau_inv(ii, jj) = obj.fermi_average(obj.E_F(i), ...
                            obj.E_F(j), obj.E(i), obj.mEff(i), ...
                            reshape(obj.W(ii, jj, :), size(obj.k)), ...
                            obj.T_e(i), obj.T_e(j));
                    end
                end
            end
            obj.tau_inv = tau_inv;
        end

        function tau_inv = calculate_parallel(obj, calc_W)
            % Calculates electron-electron scattering rates between all
            % subbands using parallelization technique.
            %
            % Syntax:
            %   tau_inv = calculate_parallel(obj, calc_W)
            %
            % Input Arguments:
            %   calc_W (logical): Flag determining if k-resolved transition
            %     rates have to be re-calculated.
            %
            % Output Arguments:
            %   tau_inv (matrix): Matrix containing the transition rates
            %     between the subbands of the two central periods.

            if nargin < 2
                calc_W = 1;
            end
            nwf = obj.num_states;
            nk = length(obj.k);
            if calc_W
                obj.set_formf
                obj.W_ikfk_jg = zeros(nk, nk, nwf, nwf, nwf, nwf);
            end

            % calculate k-resolved scattering rates
            sz = [nwf, nwf, nwf, nwf];
            maxind = prod(sz);
            W_ikfk_jg = zeros(nk, nk, maxind);
            if calc_W
                parfor ind = 1:maxind
                    [ii, ff, jj, gg] = ind2sub(sz, ind);
                    if ii ~= jj && ~(ii == gg && jj == ff)
                        W_ikfk_jg(:, :, ind) = ...
                            obj.calc_state_rate(ii, ff, jj, gg);
                        if ff == ii && gg == jj
                            W_ikfk_jg(:, :, ind) = ...
                                2 * W_ikfk_jg(:, :, ind);
                        end
                    end
                end
                for kf = 1:nk
                    for ki = 1:nk
                        obj.W_ikfk_jg(kf, ki, :, :, :, :) = ...
                            reshape(W_ikfk_jg(kf, ki, :), sz);
                    end
                end
            end

            % calculate k-averaged scattering rates
            obj.W = zeros(nwf, nwf, nk);
            tau_inv = zeros(nwf, nwf);
            for ii = 1:nwf
                for jj = 1:nwf
                    W = zeros(1, length(obj.k));
                    for ff = 1:nwf
                        for gg = 1:nwf
                            % the transition rate W_ijji does not
                            % result in a net change of carrier density
                            if ~(ii == gg && jj == ff) && ~(ii == jj)
                                i = obj.indices_i_j(ii);
                                j = obj.indices_i_j(jj);
                                f = obj.indices_i_j(ff);
                                % average over kf
                                fermi_f = FermiGoldenRule.fermi_dirac( ...
                                    obj.E(f)+phys_const.hbar^2* ...
                                    obj.k.^2./(2 * obj.mEff(i)), ...
                                    obj.E_F(f), obj.T_e(f));
                                W_ikf_jg = trapz(obj.k, reshape( ...
                                    obj.W_ikfk_jg(:, :, ii, ff, jj, gg), ...
                                    nk, nk).*reshape(fermi_f.* ...
                                    obj.k, [], 1));
                                W = W + W_ikf_jg;
                                W_ifjg = obj.fermi_average( ...
                                    obj.E_F(i), obj.E_F(j), ...
                                    obj.E(i), obj.mEff(i), ...
                                    W_ikf_jg, obj.T_e(i), obj.T_e(j));
                                tau_inv(ii, jj) = tau_inv(ii, jj) ...
                                    +W_ifjg;
                            end
                        end
                    end
                    obj.W(ii, jj, :) = W;
                end
            end
            obj.tau_inv = tau_inv;

        end

        function tau_inv = calculate_fast(obj)
            % Simplified formla for electron-electron scattering rate
            % according to "P. Hyldgaard and J. Wilkins, Phys. Rev. B
            % 53(11), 6889–6892 (1996)".
            %
            % Note:
            %   This method is currently not used!
            %
            % Syntax:
            %   tau_inv = calculate_fast(obj)
            %
            % Output Arguments:
            %   tau_inv (matrix): Matrix containing the transition rates
            %     between the subbands of the two central periods.

            nwf = obj.num_states;

            qmin = sqrt(obj.qs2);
            [Emin, Emax] = bounds(obj.E(obj.indices_i_j));
            dkmax2 = 8 * max(obj.mEff) ./ ...
                phys_const.hbar^2 * (Emax - Emin);
            kmax2 = max(obj.k).^2;
            qmax = sqrt(0.25*(2 * kmax2 + dkmax2 + 2 * ...
                sqrt(kmax2.^2+dkmax2.*kmax2))+obj.qs2);

            nk = length(obj.qval);
            obj.qval = linspace(qmin, qmax, nk);
            delta_z = abs(obj.z(2)-obj.z(1)) * 1e-10;

            for i = 1:nwf
                for j = 1:nwf
                    for f = 1:nwf
                        for g = 1:nwf
                            q = obj.qval;
                            if ~((f > g) || (i > j) || (i > f) || (i > g))
                                formfactor = obj.calc_form_factor_fft( ...
                                    i, f, j, g, q, delta_z);
                                obj.formf(:, i, f, j, g) = formfactor;
                                obj.formf(:, j, f, i, g) = formfactor;
                                obj.formf(:, i, g, j, f) = formfactor;
                                obj.formf(:, f, i, g, j) = formfactor;
                                obj.formf(:, g, j, f, i) = formfactor;
                                obj.formf(:, g, i, f, j) = formfactor;
                                obj.formf(:, j, g, i, f) = formfactor;
                                obj.formf(:, f, j, g, i) = formfactor;
                            end
                        end
                    end
                end
            end

            Ry = 13.605693122994 * phys_const.e0; % Rydberg energy (J)
            Ip = 0.785; % integrated phase space constant
            W = zeros(obj.num_states, obj.num_states, ...
                obj.num_states, obj.num_states);
            tau_inv = zeros(obj.num_states);

            for ii = 1:nwf
                i = obj.indices_i_j(ii);
                for jj = 1:nwf
                    j = obj.indices_i_j(jj);
                    if ii ~= jj
                        dE = abs(obj.E(i)-obj.E(j)); % energy transfer
                        q = sqrt(2*obj.mEff(i)*dE/ ...
                            phys_const.hbar.^2); % momentum transfer
                        mu = obj.E_F(i);
                        if dE > 0
                            factor = Ry * mu / ...
                                (pi^2 * dE * phys_const.hbar);
                            tau_inv(ii, jj) = factor * (2 * pi * ...
                                obj.calc_form_factor_fft(i, i, j, j, ...
                                q, delta_z)).^2 * Ip;
                        end
                    end
                end
            end
            obj.tau_inv = tau_inv;
        end

        function W_ikfk_jg = calc_state_rate(obj, ii, ff, jj, gg)
            %CALC_STATE_RATE Calculates the transition rate for the process
            % in which the first electron scatters from state ik to subband
            % j and that second electron scatters from state fk to subband
            % g according to https://doi.org/10.1063/1.4863665 (Note that
            % in the literature the indeces f and j are often interchanged.)
            %
            % Syntax:
            %   W_ikfk_jg = calc_state_rate(obj, ii, ff, jj, gg)
            %
            % Input Arguments:
            %   ii (scalar): Inital subband index of first electron.
            %   ff (scalar): Inital subband index of second electron.
            %   jj (scalar): Final subband index of first electron.
            %   gg (scalar): Final subband index of second electron.
            %
            % Output Arguments:
            %   W_ikfk_jg (scalar): Transition rate from states ik, fk to
            %     subbands j, g.

            i = obj.indices_i_j(ii);
            f = obj.indices_i_j(ff);
            j = obj.indices_i_j(jj);
            g = obj.indices_i_j(gg);

            nk = length(obj.k);
            W_ikfk_jg = zeros(nk);

            % The prefactor 4 is introduced, since the integration
            % over theta and alpha is performed over the range [0, pi].
            factor = 4 * obj.mEff(j) * phys_const.e0^4 / ...
                (4 * pi * phys_const.hbar^3 * ...
                (4 * pi * obj.epsr_stat(2) * phys_const.eps0)^2);

            % interpolate pre-calculated values of the form factor
            x = obj.qval;
            y = obj.formf(:, ii, ff, jj, gg);
            for ki = 1:nk
                q = obj.calc_q_values(i, f, j, g, ki, 1:ki);
                value = interp1(x, y, q, 'spline');
                if any(isnan(value))
                    disp("NaN value")
                end
                % integration over alpha and theta
                step_theta = abs(obj.theta(2)-obj.theta(1));
                W_ikfk_jg(1:ki, ki) = trapz(step_theta, ...
                    trapz(obj.theta, factor.*value.^2./ ...
                    (q + sqrt(obj.qs2)).^2));
                W_ikfk_jg(ki, 1:ki-1) = W_ikfk_jg(1:ki-1, ki)';
            end
        end

        function q_values = calc_q_values(obj, i, f, j, g, ki, kf)
            %CALC_Q_VALUES Calculates all q values which are needed for
            % the calculation of a specifc rate W_ikfjg for a fixed number
            % of alpha and theta values.
            %
            % Syntax:
            %   q_values = calc_q_values(obj, i, f, j, g, ki, kf)
            %
            % Input Arguments:
            %   i (scalar): Initial subband (first electron).
            %   f (scalar): Initial subband (second electron).
            %   j (scalar): Final subband (first electron).
            %   g (scalar): Final subband (second electron).
            %   ki (scalar): Wavevector of inital state (first electron).
            %   kf (scalar | vector): Wavevector of inital state (second
            %    electron).
            %
            % Output Arguments:
            %   q_values (3-d array): Returns array of q-values in
            %     dependence on angle theta (1st dimension), angle alpha
            %     (2nd dimension) und wavevector kf (3rd dimension).

            nth = length(obj.theta);
            mEff = obj.mEff(i); % approximation: m* = m_i
            delta_k02 = 4 * mEff / (phys_const.hbar^2) * ...
                (obj.E(i) + obj.E(f) - obj.E(j) - obj.E(g));
            q_values = zeros(nth, nth, length(kf));

            for m = 1:nth % loop over theta
                for n = 1:nth % loop over alpha
                    for o = kf
                        q2 = obj.calc_q2(obj.k(ki), obj.k(o), ...
                            delta_k02, obj.theta(n), obj.theta(m));
                        % only real and positive q2 are valid
                        if (q2 > 0 && isreal(q2))
                            q_values(m, n, o) = sqrt(q2);
                        end
                    end
                end
            end
        end

        function set_formf(obj)
            %SET_FORMF Inititalize form factor for range of homogeneously
            % spaced q values.
            %
            % Syntax:
            %   set_formf(obj)

            tinit = tic;

            nwf = obj.num_states;
            nk = length(obj.qval);
            obj.formf = zeros(length(obj.qval), nwf, nwf, nwf, nwf);

            % calculate upper and lower bound of q values
            [Emin, Emax] = bounds(obj.E(obj.indices_i_j));
            dkmax2 = 8 * max(obj.mEff) ./ ...
                phys_const.hbar^2 * (Emax - Emin);
            kmax2 = max(obj.k).^2;
            qmax = sqrt(0.25*(2 * kmax2 + dkmax2 + ...
                2 * sqrt(kmax2.^2+dkmax2.*kmax2)));

            obj.qval = linspace(0, qmax, nk);
            delta_z = abs(obj.z(2)-obj.z(1)) * 1e-10;

            q = obj.qval;

            for i = 1:nwf
                for j = 1:nwf
                    for f = 1:nwf
                        for g = 1:nwf
                            if ~((f > g) || (i > j) || (i > f) || (i > g))
                                formfactor = obj.calc_form_factor_fft( ...
                                    i, f, j, g, q, delta_z);
                                obj.formf(:, i, f, j, g) = formfactor;

                                obj.formf(:, j, f, i, g) = formfactor;
                                obj.formf(:, i, g, j, f) = formfactor;

                                obj.formf(:, f, i, g, j) = formfactor;
                                obj.formf(:, g, j, f, i) = formfactor;

                                obj.formf(:, g, i, f, j) = formfactor;
                                obj.formf(:, j, g, i, f) = formfactor;
                                obj.formf(:, f, j, g, i) = formfactor;
                            end
                        end
                    end
                end
            end
            fprintf(['initialize form factor: ', num2str( ...
                toc(tinit)), 's\n'])
        end

        function result = calc_form_factor_fft(obj, i, f, j, g, q, delta_z)
            %CALC_FORM_FACTOR_FFT Calculates formfactor using Fourier
            % transform method described in
            % https://doi.org/10.1063/5.0041392.
            %
            % Syntax:
            %   result = calc_form_factor_fft(obj, i, f, j, g, q, delta_z)
            %
            % Input Arguments:
            %   i (scalar): Initial subband (first electron).
            %   f (scalar): Initial subband (second electron).
            %   j (scalar): Final subband (first electron).
            %   g (scalar): Final subband (second electron).
            %   q (matrix): Exchanged wavevector in dependence on the
            %     scattering angles theta and alpha.
            %   delta_z (scalar): Grid spacing [m].
            %
            % Output Arguments:
            %   result (matrix): Formfactor for the transitions from
            %     subbands if -> jg in dependence on the scattering angles
            %     theta and alpha.

            psi_prod1_ft = reshape(obj.psi2_ft(i, j, :), 1, []);
            psi_prod2_ft = reshape(obj.psi2_ft(f, g, :), 1, []);

            N = length(psi_prod1_ft);
            N_f = (N - 1) / 2;
            c = -N_f:N_f;

            nominator = psi_prod1_ft .* flip(psi_prod2_ft);
            two_pi_c2 = (2 * pi * c).^2;

            integral = zeros(size(q));

            for m = 1:size(q, 1)
                for n = 1:size(q, 2)
                    denominator = (q(m, n) * N * delta_z)^2 + two_pi_c2;
                    integral(m, n) = 2 * N * q(m, n) * ...
                        delta_z^3 * sum(nominator./denominator);
                end
            end
            result = abs(integral);
            if any(q < 1e2)
                % Since FT method diverges for small q, use conventional
                % integration scheme instead
                q_small = q(q < 1e2);
                integral_small = abs(calc_form_factor( ...
                    obj, i, f, j, g, q_small, delta_z));
                result(q < 1e2) = integral_small;
            end
        end

        function result = calc_form_factor(obj, i, f, j, g, q, delta_z)
            %CALC_FORM_FACTOR Calculates form factor A_ifjg(q) using normal
            % integration method (without FFT). Used for small q values,
            % since FFT-method might diverge.
            %
            % Syntax:
            %   result = calc_form_factor(obj, i, f, j, g, q, delta_z)
            %
            % Input Arguments:
            %   i (scalar): Initial subband (first electron).
            %   f (scalar): Initial subband (second electron).
            %   j (scalar): Final subband (first electron).
            %   g (scalar): Final subband (second electron).
            %   q (matrix): Exchanged wavevector in dependence on the
            %     scattering angles theta and alpha.
            %   delta_z (scalar): Grid spacing [m].
            %
            % Output Arguments:
            %   result (matrix): Formfactor for the transitions from
            %     subbands if to subbands jg in dependence on the scattering
            %     angles theta and alpha.

            zPos = obj.z * 1e-10;
            result = zeros(length(q), 1);

            exp_qZpos = exp((q' .* zPos));
            psi1 = obj.psi(:, i) .* conj(obj.psi(:, j));
            psi2 = obj.psi(:, f) .* conj(obj.psi(:, g));

            for n = 1:length(q)
                % Integration over dz
                I = psi2 ./ exp_qZpos(:, n) .* ...
                    cumtrapz(psi1.*exp_qZpos(:, n)) * delta_z;
                I = I + psi2 .* exp_qZpos(:, n) .* ...
                    flip(delta_z*cumtrapz(flip(psi1./exp_qZpos(:, n))));
                % Integration over dz'
                result(n) = trapz(I) * delta_z;
            end
        end

        function update(obj, ns)
            %UPDATE Updates the sheet densities, quasi Fermi-levels and
            % screening wavevector.
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

        function set_eigenstates(obj, eigen)
            %SET_EIGENSTATES Sets new eigenenergies, effective masses and
            % wavefunctions.
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

        function psi2_ft = set_psi2_ft(obj)
            % Calculate Fourier transform of product of wavefunctions
            % for fast calculation of the form factor as described in
            % https://doi.org/10.1063/5.0041392.
            %
            % Syntax:
            %   psi2_ft = set_psi2_ft(obj)
            %
            % Output Arguments:
            %   psi2_ft (3-d array): Returns the Fourier transforms of
            %     products of two wavefunctions in dependence on the initial
            %     subband (1st dimension), final subband (2nd dimension) and
            %     k-space variable of the Fourier transform (3rd dimension).

            % product of two wavefunctions
            psi2 = zeros(size(obj.psi, 2), ...
                size(obj.psi, 2), size(obj.psi, 1));
            for i = 1:size(obj.psi, 2)
                for j = 1:size(obj.psi, 2)
                    psi2(i, j, :) = obj.psi(:, i) .* obj.psi(:, j);
                end
            end

            N = size(psi2, 3); % number of grid points in z-direction
            if rem(N, 2) == 0
                psi2 = psi2(:, :, 1:end-1);
                N = size(psi2, 3);
            end
            N_f = (N - 1) / 2;
            m = 1:N;

            % phase factor could be ignored, but will be kept for correct
            % transformation of FT conventions
            phase_factor = fftshift(exp(1i*2*pi*N_f*(m - 1)/N)');
            psi2_ft = fftshift(fft(psi2, [], 3), 3);

            for i = 1:size(obj.psi, 2)
                for j = 1:size(obj.psi, 2)
                    psi2_ft(i, j, :) = reshape( ...
                        psi2_ft(i, j, :), [], 1) .* phase_factor;
                end
            end
            obj.psi2_ft = psi2_ft;
        end
    end

    methods (Static)

        function q2 = calc_q2(k_i, k_f, delta_k02, alpha, theta)
            %CALC_Q2 Calculates the exchanged wavevector squared q2, which
            % accounts for energy and momentum transfer in the scattering
            % process.
            %
            % Syntax:
            %   q2 = calc_q2(k_i, k_f, delta_k02, alpha, theta)
            %
            % Input Arguments:
            %   k_i (scalar): Wavevector of inital state.
            %   k_f (scalar): Wavevector of final state.
            %   delta_k02 (scalar): Difference in wavevectors.
            %   alpha (scalar): Scattering angle.
            %   theta (scalar): Scattering angle.
            %
            % Output Arguments:
            %   q2: Exchanged wavevector squared.

            k_if2 = k_i^2 + k_f.^2 - 2 * k_i .* k_f .* cos(alpha);
            k_jg2 = k_if2 + delta_k02;

            k_jg2(k_jg2 <= 0) = 0;
            k_if2(k_if2 <= 0) = 0;
            q2 = (k_if2 + k_jg2 - 2 * sqrt(k_if2.*k_jg2) .* ...
                cos(theta)) / 4;
        end
    end
end
