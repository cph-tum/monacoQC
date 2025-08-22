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

classdef ImpurityRate < FermiGoldenRule
    % Calculates intersubband transition rates due to scattering of
    % electrons by ionized impurities based on Fermi's golden rule.

    properties
        dop % vector: Doping concentration at each grid point [1/m^3].
        eps % vector: Static permittivity eps_0*eps_r at each grid point.
    end

    methods
        function obj = ImpurityRate(eigen, device, scenario, options)
            % Constructs an object of type ImpurityRate.
            %
            % Syntax:
            %   obj = ImpurityRate(eigen, device, scenario)
            %   obj = ImpurityRate(eigen, device, scenario, options)
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
            obj.Name = 'impurity scattering';

            Ndop = [];
            epsr_stat = [];
            % get physical properties of each layer
            for j = 1:length(device.layers)
                Ndop(j) = device.layers{j}.doping * 1e6;
                epsr_stat(j) = device.layers{j}.material.eps_r;
            end

            % get properties at each grid point
            obj.dop = zeros(1, length(obj.z));
            obj.eps = zeros(1, length(obj.z));
            PosI = obj.PosInterf(2:end);
            c = 1;
            for n = 1:length(obj.z)
                if c < length(PosI) && obj.z(n) >= PosI(c)
                    c = c + 1;
                end
                obj.dop(n) = Ndop(c);
                obj.eps(n) = phys_const.eps0 * epsr_stat(c);
            end

            % Shift position vector in order to avoid numerical
            % issues when calculating transition rates
            obj.z = (obj.z(end) - obj.z(1)) / 2 + obj.z;

            % calculate screened wave vector qs2
            obj.qs2 = obj.calc_qs2();
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
            theta = linspace(0, pi, 100);

            kPrime2 = k^2 * obj.mEff(j) / obj.mEff(i) + 2 * obj.mEff(j) ...
                * (obj.E(i) - obj.E(j)) / phys_const.hbar^2;
            q = sqrt(max(k^2+kPrime2-2*k*sqrt(kPrime2)* ...
                cos(theta), 0));

            % exclude imaginary q values
            q = q(q >= 0);
            theta = theta(q >= 0);

            psi2 = obj.psi(:, i) .* conj(obj.psi(:, j)) ./ obj.eps';
            pre_fact = obj.mEff(j) * phys_const.e0^4 / ...
                (4 * pi * phys_const.hbar^3);
            zPos = obj.z * 1e-10;

            if kPrime2 >= 0
                F_ij = zeros(1, length(theta));
                for n = 1:length(q)
                    % Integration over dz
                    I = exp(-q(n)*zPos) .* cumtrapz(zPos, ...
                        psi2.*exp(q(n)*zPos)) + exp(q(n)*zPos) .* ...
                        flip(cumtrapz(zPos, flip(psi2.*exp(-q(n)*zPos))));
                    % Integration over dz'
                    F_ij(n) = trapz(zPos, obj.dop'.*I.^2);
                end
                % Integration over theta
                W_ikj = pre_fact * trapz(theta, ...
                    F_ij./(q + sqrt(obj.qs2)).^2);
            else
                % no scattering for kPrime2 < 0
                W_ikj = 0;
            end
            if isnan(W_ikj) || isinf(abs(W_ikj))
                error('NaN or Inf value!')
            end
        end

        function puredep = calulate_pure_dephasing_rate(obj, calc_rate)
            % Calculates pure dephasing rates for all subbands.
            %
            % Syntax:
            %   puredep = calulate_pure_dephasing_rate(obj)
            %   puredep = calulate_pure_dephasing_rate(obj, calc_rate)
            %
            % Input Arguments:
            %   calc_rate (logical): Flag specifying if k-dependent pure
            %     dephasing rates have to be recalculated.
            %
            % Output Arguments:
            %   puredep (matrix): Matrix containing the pure dephasing
            %     rates for the two central periods.

            if nargin < 2
                calc_rate = 1;
            end

            if calc_rate
                % Call superclass method to initialize arrays
                calulate_pure_dephasing_rate@FermiGoldenRule( ...
                    obj, calc_rate);
            end

            for ii = 1:obj.num_states
                for jj = ii:obj.num_states
                    if ii ~= jj
                        i = obj.indices_i_j(ii);
                        j = obj.indices_i_j(jj);

                        % calculate k-resolved pure dephasing rates
                        if calc_rate
                            for l = 1:length(obj.k)
                                g_ikj = calc_puredep(obj, i, j, obj.k(l));
                                obj.pure_dephasing_k(ii, jj, l) = g_ikj;
                                obj.pure_dephasing_k(jj, ii, l) = g_ikj;
                            end
                        end

                        % calculate k-averaged pure dephasing rates
                        if mod(ii, obj.num_states/2) == ...
                                mod(jj, obj.num_states/2)
                            % Check if subband populations of |i> and |j>
                            % are equal, for which the avg pure dephasing
                            % rate is 0
                            obj.pure_dephasing(ii, jj) = ...
                                obj.fermi_average_no_blocking( ...
                                obj.E_F(i), obj.E(i), obj.mEff(i), ...
                                reshape(obj.pure_dephasing_k(ii, jj, :), ...
                                size(obj.k)), obj.T_e(i));
                            obj.pure_dephasing(jj, ii) = ...
                                obj.fermi_average_no_blocking( ...
                                obj.E_F(j), obj.E(j), obj.mEff(j), ...
                                reshape(obj.pure_dephasing_k(jj, ii, :), ...
                                size(obj.k)), obj.T_e(j));
                        else
                            obj.pure_dephasing(ii, jj) = ...
                                obj.pop_inv_average(obj.E_F(i), ...
                                obj.E_F(j), obj.E(i), obj.E(j), ...
                                obj.mEff(i), obj.mEff(j), reshape( ...
                                obj.pure_dephasing_k(ii, jj, :), ...
                                size(obj.k)), obj.T_e(i), obj.T_e(j));
                            obj.pure_dephasing(jj, ii) = ...
                                obj. pop_inv_average(obj.E_F(j), ...
                                obj.E_F(i), obj.E(j), obj.E(i), ...
                                obj.mEff(j), obj.mEff(i), reshape( ...
                                obj.pure_dephasing_k(jj, ii, :), ...
                                size(obj.k)), obj.T_e(j), obj.T_e(i));
                        end
                    end
                end
            end
            puredep = obj.pure_dephasing;
        end

        function puredep = calulate_pure_dephasing_rate_parallel( ...
                obj, calc_rate)
            % Calculates pure dephasing rates for all subbands on parallel
            % workes.
            %
            % Syntax:
            %   puredep = calulate_pure_dephasing_rate(obj)
            %   puredep = calulate_pure_dephasing_rate(obj, calc_rate)
            %
            % Input Arguments:
            %   calc_rate (logical): Flag specifying if k-dependent pure
            %     dephasing rates have to be recalculated.
            %
            % Output Arguments:
            %   puredep (matrix): Matrix containing the pure dephasing
            %     rates for the two central periods.

            if nargin < 2
                calc_rate = 1;
            end

            nwf = obj.num_states;

            % calculate k-resolved pure dephasing rates
            if calc_rate
                calulate_pure_dephasing_rate_parallel@FermiGoldenRule( ...
                    obj, calc_rate); % initialize arrays
                sz = [nwf, nwf, length(obj.k)];
                maxind = prod(sz);
                g_ikj = zeros(maxind, 1);
                parfor ind = 1:maxind
                    [ii, jj, l] = ind2sub(sz, ind);
                    i = obj.indices_i_j(ii);
                    j = obj.indices_i_j(jj);
                    % use symmetry of pure dephasing rates
                    if i > j
                        g_ikj(ind) = calc_puredep(obj, i, j, obj.k(l));
                    end
                end
                obj.pure_dephasing_k = reshape(g_ikj, sz);
            end

            % calculate k-averaged pure dephasing rates
            for ii = 1:obj.num_states
                for jj = ii:obj.num_states
                    obj.pure_dephasing_k(ii, jj, :) = ...
                        obj.pure_dephasing_k(jj, ii, :);
                    % Check if subband populations of |i> and |j> are
                    % equal, for which the avg pure dephasing rate is 0
                    i = obj.indices_i_j(ii);
                    j = obj.indices_i_j(jj);
                    if mod(ii, obj.num_states/2) == ...
                            mod(jj, obj.num_states/2)
                        obj.pure_dephasing(ii, jj) = ...
                            obj.fermi_average_no_blocking(obj.E_F(i), ...
                            obj.E(i), obj.mEff(i), reshape( ...
                            obj.pure_dephasing_k(ii, jj, :), ...
                            size(obj.k)), obj.T_e(i));
                        obj.pure_dephasing(jj, ii) = ...
                            obj.fermi_average_no_blocking(obj.E_F(j), ...
                            obj.E(j), obj.mEff(j), ...
                            reshape(obj.pure_dephasing_k(jj, ii, :), ...
                            size(obj.k)), obj.T_e(j));
                    else
                        % Average over the inversion between ik and jk
                        obj.pure_dephasing(ii, jj) = ...
                            obj.pop_inv_average(obj.E_F(i), ...
                            obj.E_F(j), obj.E(i), obj.E(j), ...
                            obj.mEff(i), obj.mEff(j), reshape( ...
                            obj.pure_dephasing_k(ii, jj, :), ...
                            size(obj.k)), obj.T_e(i), obj.T_e(j));
                        obj.pure_dephasing(jj, ii) = ...
                            obj. pop_inv_average(obj.E_F(j), ...
                            obj.E_F(i), obj.E(j), obj.E(i), ...
                            obj.mEff(j), obj.mEff(i), reshape( ...
                            obj.pure_dephasing_k(jj, ii, :), ...
                            size(obj.k)), obj.T_e(j), obj.T_e(i));
                    end
                end
            end
            puredep = obj.pure_dephasing;
        end

        function puredep = calc_puredep(obj, i, j, k)
            % Calculates k-resolved pure dephasing rates according to
            % https://doi.org/10.1063/1.5005618.
            %
            % Syntax:
            %   puredep = calc_puredep(obj, i, j, k)
            %
            % Input Arguments:
            %   i (scalar): Inital subband.
            %   j (scalar): Final subband.
            %   k (scalar): Wavevector of inital state.
            %
            % Output Arguments:
            %   puredep (scalar): Rure dephasing rate.

            % grid points for integration
            theta = linspace(0, pi, 100);

            q2 = 2 * k.^2 * (1 - cos(theta));
            q = sqrt(max(q2, 0));
            pi2 = obj.psi(:, i).^2;
            pf2 = obj.psi(:, j).^2;

            % Calculate rate
            I = zeros(1, length(q));
            zPos = obj.z * 1e-10;

            for n = 1:length(q)
                % Calculate integrals over dz
                Ii = exp(-q(n)*zPos) .* cumtrapz(zPos, pi2.*exp(q(n)*zPos));
                If = exp(-q(n)*zPos) .* cumtrapz(zPos, pf2.*exp(q(n)*zPos));
                Ii = Ii + exp(q(n)*zPos) .* flip(cumtrapz(zPos, ...
                    flip(pi2.*exp(-q(n)*zPos))));
                If = If + exp(q(n)*zPos) .* flip(cumtrapz(zPos, ...
                    flip(pf2.*exp(-q(n)*zPos))));
                f = obj.mEff(i) * Ii.^2 - (obj.mEff(i) + obj.mEff(j)) * ...
                    Ii .* If + obj.mEff(j) * If.^2;
                % Integrate over dz'
                I(n) = trapz(zPos, obj.dop.*f');
            end

            % Integrate over Theta
            I2 = trapz(theta, I./(q + sqrt(obj.qs2)).^2);
            eps = phys_const.eps0 * obj.epsr_const;
            puredep = phys_const.e0^4 / ...
                (8 * pi * eps^2 * phys_const.hbar^3) * I2;
        end

        function update(obj, ns)
            % Updates the sheet densities, quasi Fermi-levels and screening
            % wavevector.
            %
            % Syntax:
            %   update(obj, ns)
            %
            % Input Arguments:
            %   ns (vector): Sheet densities for all four QCL-periods [1/m^2].

            update@FermiGoldenRule(obj, ns);
            % Calculate new screened wave vector qs2
            obj.qs2 = obj.calc_qs2();
        end
    end
end
