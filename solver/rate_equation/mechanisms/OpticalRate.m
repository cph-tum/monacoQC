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

classdef OpticalRate < FermiGoldenRule
    % Calculates intersubband transition rates due to the scattering of
    % electrons by photons using Fermi's golden rule. This scattering
    % mechanism accounts for stimulated emission and absorption.

    properties
        intensity % vector: Optical intensity of the cavity modes [W/m^2].
        gain % vector: Optical gain of the cavity modes [1/m].
        deph_rates %matrix: Dephasing rates [1/s].
        deph_rates_k % 3-d array: Dephasing rates k-resolved [1/s].
        freq % vector: Frequencies of cavity modes [1/s].
        tstep % Scalar: Time step for field evolution [s].
        Gamma % scalar | vector: Confinement factor (overlap factor).
        a_loss % scalar: Total power loss coefficient [1/m].
        n_eff % scalar: Effective refractive index of the device.
        d_ij % matrix: Transition dipole moments [Cm].
    end

    methods
        function obj = OpticalRate(eigen, device, scenario, options)
            % Constructs an object of type OpticalRate.
            %
            % Syntax:
            %   obj = OpticalRate(eigen, device, scenario)
            %   obj = OpticalRate(eigen, device, scenario, options)
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
            obj.Name = 'optical field';

            % set properties of the cavity
            obj.Gamma = device.waveguide.overlap_factor;
            obj.a_loss = device.waveguide.a_power;
            L = device.waveguide.l_waveguide;
            obj.n_eff = device.n_eff;

            obj.deph_rates = zeros(obj.num_states, obj.num_states);
            obj.deph_rates_k = zeros(obj.num_states, ...
                obj.num_states, length(obj.k));
            obj.d_ij = zeros(obj.num_states);
            for ii = 1:obj.num_states
                for jj = 1:obj.num_states
                    i = obj.indices_i_j(ii);
                    j = obj.indices_i_j(jj);
                    obj.d_ij(ii, jj) = eigen.get_dipole_element(i, j);
                end
            end

            obj.tstep = 5e-12;
            df = phys_const.c0 / (2 * L * obj.n_eff);
            obj.freq = scenario.fmin:df:scenario.fmax;
            if isempty(obj.freq)
                error("Frequency bounds fmin and fmax are not specified. "+ ...
                    "Please set the corresponding properties in "+ ...
                    "scenario accordingly!")
            end
            obj.intensity = zeros(1, length(obj.freq));
        end

        function tau_inv = calculate(obj)
            % Calculates transition rates between all subbands.
            %
            % Syntax:
            %   tau_inv = calculate(obj)
            %
            % Output Arguments:
            %   tau_inv (matrix): Matrix containing the transition rates
            %     between the subbands of the two central periods.

            obj.W = zeros(obj.num_states, obj.num_states, length(obj.k));
            tau_inv = zeros(obj.num_states, obj.num_states);

            for ii = 1:obj.num_states
                for jj = 1:obj.num_states
                    if ii ~= jj
                        i = obj.indices_i_j(ii);
                        j = obj.indices_i_j(jj);
                        for l = 1:length(obj.k)
                            obj.W(ii, jj, l) = calc_state_rate( ...
                                obj, i, j, obj.k(l));
                        end
                        % Calculates subband transition rate i->j from
                        % state transition rates by averaging over the
                        % subbands
                        tau_inv(ii, jj) = obj.fermi_average(obj.E_F(i), ...
                            obj.E_F(j), obj.E(i), obj.mEff(i), ...
                            reshape(obj.W(ii, jj, :), size(obj.k)), ...
                            obj.T_e(i), obj.T_e(j));
                    end
                end
            end
            obj.tau_inv = tau_inv;
        end

        function tau_inv = calc_state_rate_avg(obj, i, j)
            % Calculates transition rate from subband i to subband j
            % using k-averaged dephasing rates according to
            % https://doi.org/10.1063/1.3284523.
            %
            % Syntax:
            %   W_ikj = calc_state_rate_avg(obj, i, j, k)
            %
            % Input Arguments:
            %   i (scalar): Inital subband.
            %   j (scalar): Final subband.
            %   k (scalar): Wavevector of the inital state.
            %
            % Output Arguments:
            %   W_ikj (scalar): Transition rate.

            ii = i - min(obj.indices_i_j) + 1;
            jj = j - min(obj.indices_i_j) + 1;

            omega = 2 * pi * obj.freq; % mode frequency

            gamma_ij = obj.deph_rates(ii, jj); % dephasing rate
            omega_ij = (obj.E(i) - obj.E(j)) ./ ...
                phys_const.hbar; % transition frequency
            lorentz = 1 ./ (1 + (omega - abs(omega_ij)).^2 ./ ...
                gamma_ij.^2) ./ gamma_ij; % lineshape

            tau_inv = 1 ./ (phys_const.hbar.^2 * phys_const.eps0 .* ...
                obj.n_eff .* phys_const.c0) .* obj.d_ij(ii, jj).^2 .* ...
                sum(obj.intensity.*lorentz);
        end

        function W_ikj = calc_state_rate(obj, i, j, k)
            % Calculates transition rate from state ik to subband j using
            % k-dependent dephasing rates according to
            % https://doi.org/10.1063/1.3284523.
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

            ii = i - min(obj.indices_i_j) + 1;
            jj = j - min(obj.indices_i_j) + 1;
            kk = find(obj.k == k);

            % mode frequencies
            omega = 2 * pi * obj.freq;
            % dephasing rate
            gamma_ij = obj.deph_rates_k(ii, jj, kk);
            % transition frequency
            Eik = obj.E(i) + (k * phys_const.hbar).^2 ./ (2 * obj.mEff(i));
            Ejk = obj.E(j) + (k * phys_const.hbar).^2 ./ (2 * obj.mEff(j));
            omega_ij = (Eik - Ejk) / phys_const.hbar;
            % lorentzian lineshape
            lorentz = 1 ./ (1 + (omega - abs(omega_ij)).^2 ./ ...
                gamma_ij.^2) ./ gamma_ij;

            W_ikj = 1 ./ (phys_const.hbar.^2 * phys_const.eps0 .* ...
                obj.n_eff .* phys_const.c0) .* obj.d_ij(ii, jj).^2 .* ...
                sum(obj.intensity.*lorentz);
        end

        function g = calc_gain(obj, nS)
            % Calculates gain spectrum without considering the k-dependance
            % of the dephasing rates and sheet densities.
            %
            % Syntax:
            %   g = calc_gain(obj, nS)
            %
            % Input Arguments:
            %   nS (vector): Sheet densities of the two central QCL-periods
            %     [1/m^2].
            %
            % Output Arguments:
            %   g (vector): Optical gain of the cavity modes [1/m].

            omega = 2 * pi * obj.freq;
            g = zeros(size(omega));

            for ii = 1:obj.num_states
                for jj = 1:obj.num_states
                    % ensure that gain by intramodule scattering is not
                    % calculated twice
                    if ii < obj.num_states / 2 + 1 || ...
                            jj < obj.num_states / 2 + 1
                        i = obj.indices_i_j(ii);
                        j = obj.indices_i_j(jj);
                        % transition rates
                        gamma_ij = obj.deph_rates(ii, jj);
                        % lorentzian lineshape
                        omega_ij = (obj.E(i) - ...
                            obj.E(j)) ./ phys_const.hbar;
                        lorentz = 1 ./ (1 + (omega - ...
                            abs(omega_ij)).^2 ./ gamma_ij.^2) ./ gamma_ij;
                        % gain
                        g_ij = sign(omega_ij) * omega ./ ...
                            (phys_const.hbar * phys_const.eps0 ...
                            * phys_const.c0 * obj.n_eff * ...
                            obj.l_period) .* nS(i) .* ...
                            obj.d_ij(ii, jj).^2 .* lorentz;
                        g = g + g_ij;
                    end
                end
            end
        end

        function g = calc_gain_k_resolved(obj, nS)
            % Calculates gain spectrum using k-dependent dephasing rates
            % and sheet densities.
            %
            % Syntax:
            %   g = calc_gain_k_resolved(obj, nS)
            %
            % Input Arguments:
            %   nS (matrix): Sheet densities in dependence on wavevector k
            %     (1st dimension) for all states in the two central
            %     QCL-periods (2nd dimension) [1/m^2].
            %
            % Output Arguments:
            %   g (vector): Optical gain of the cavity modes [1/m].

            omega = 2 * pi * obj.freq;
            g = zeros(size(omega));

            % k-resolved sheet density
            num_wfs = obj.num_states / 2;

            for ii = 1:obj.num_states
                for jj = 1:obj.num_states
                    % ensure that gain by intramodule scattering is not
                    % calculated twice
                    if ii < obj.num_states / 2 + 1 || ...
                            jj < obj.num_states / 2 + 1
                        for kk = 1:length(obj.k)
                            i = obj.indices_i_j(ii);
                            j = obj.indices_i_j(jj);
                            % dephasing rates
                            gamma_ijk = obj.deph_rates_k(ii, jj, kk);
                            % transition frequency
                            Eik = obj.E(i) + obj.k(kk).^2 * ...
                                phys_const.hbar.^2 / (2 * obj.mEff(i));
                            Ejk = obj.E(j) + obj.k(kk).^2 * ...
                                phys_const.hbar.^2 / (2 * obj.mEff(j));
                            omega_ijk = (Eik - Ejk) / phys_const.hbar;
                            % lorentzian lineshape
                            lorentz = 1 ./ (1 + (omega - ...
                                abs(omega_ijk)).^2 ./ ...
                                gamma_ijk.^2) ./ gamma_ijk;
                            % gain
                            g_ijk = sign(omega_ijk) * omega ./ ...
                                (phys_const.hbar * phys_const.eps0 ...
                                * phys_const.c0 * obj.n_eff * ...
                                obj.l_period) .* nS(kk, ...
                                mod(i-1, num_wfs)+1) .* ...
                                obj.d_ij(ii, jj).^2 .* lorentz;
                            g = g + g_ijk;
                        end
                    end
                end
            end
        end

        function [I, gain] = calc_intensity(obj, I, ns)
            % Calculates intensity and gain for a new time step and updates
            % the respective properties. Depending on the dimensions of ns,
            % k-dependent or k-averaged dephasing rates are used for the
            % calculation. If k-dependent dephasing rates should be used,
            % the first dimension of ns must have the same size as the
            % discretized k-vector.
            %
            % Syntax:
            %   [I, gain] = calc_intensity(obj, I, ns)
            %
            % Input Arguements:
            %   I (vector): Intensities of the cavity modes from previous
            %     iteration.
            %   ns (vector | matrix): Sheet densities of the two central
            %     QCL-periods.
            %
            % Output Arguments:
            %   I (vector): New mode intensities.
            %   gain (vector): New gain spectrum.

            if size(ns, 1) == length(obj.k)
                gain = obj.calc_gain_k_resolved(ns);
            else
                gain = obj.calc_gain(ns);
            end
            I = I .* exp((obj.Gamma .* gain - ...
                obj.a_loss).*phys_const.c0.*obj.tstep./obj.n_eff);
            obj.intensity = I;
            obj.gain = gain;
        end

        function update(obj, ns, lvl_broadening, pure_dephasing)
            % Updates electron density, quasi Fermi-levels and optionally
            % also the dephasing rates.
            %
            % Syntax:
            %   update(obj, ns)
            %   update(obj, ns, lvl_broadening, pure_dephasing)
            %
            % Input Arguments:
            %   ns (vector): Sheet densities of all four QCL-periods.
            %   lvl_broadening (vector | matrix): Lifetime broadening of
            %     each subband in the two central QCL-periods (1st
            %     dimension). The k-dependence can be specified
            %     as well (2nd dimension).
            %   pure_dephasing (matrix | 3-d array): Pure dephasing
            %     rates for all subbands in the two central QCL-eriods (1st/
            %     2nd dimension). The k-dependence can be specified
            %     as well (3rd dimension).

            obj.E_F = FermiGoldenRule.calc_fermi_energy( ...
                obj.E, obj.mEff, ns, obj.T_e);
            obj.electron_density = ns;

            if nargin > 2
                obj.set_dephasing_rates(lvl_broadening, pure_dephasing);
            end
        end

        function set_dephasing_rates(obj, lvl_broadening, pure_dephasing)
            % Sets new dephasing rates based on level broadening and pure
            % dephasing rates.
            %
            % Syntax:
            %   set_dephasing_rates(obj, lvl_broadening, pure_dephasing)
            %
            % Input Arguments:
            %   lvl_broadening (vector | matrix): Lifetime broadening of
            %     each subband in the two central QCL-periods.
            %   pure_dephasing (matrix | 3-d array): Pure dephasing
            %     rates for all subbands in the two central QCL-periods.

            for i = 1:obj.num_states
                for j = 1:obj.num_states
                    if length(size(pure_dephasing)) == 2
                        % k-independent
                        obj.deph_rates(i, j) = 0.5 * (lvl_broadening(i) ...
                            +lvl_broadening(j)) + pure_dephasing(i, j);
                    elseif length(size(pure_dephasing)) == 3
                        % k-dependent
                        obj.deph_rates_k(i, j, :) = 0.5 * ...
                            (lvl_broadening(i, :) + lvl_broadening(j, :)) ...
                            +reshape(pure_dephasing(i, j, :), ...
                            1, length(obj.k));
                    end
                end
            end
        end
    end
end
