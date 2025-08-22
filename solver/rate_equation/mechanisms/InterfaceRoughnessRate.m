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

classdef InterfaceRoughnessRate < FermiGoldenRule
    % Calculates intersubband transition rates due to scattering of
    % electrons at material interfaces using Fermi's golden rule.

    properties
        Lambda % scalar: Correlation length [m].
        Delta % scalar: Standard deviation of interface position [m].
        CBoffset % scalar: Conduction band offset [J].
    end

    methods
        function obj = InterfaceRoughnessRate(eigen, device, scenario, ...
                options)
            % Construct an object of type InterfaceRoughnessRate.
            %
            % Syntax:
            %   obj = InterfaceRoughnessRate(eigen, device, scenario)
            %   obj = InterfaceRoughnessRate(eigen, device, scenario, options)
            %
            % Input:
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

            obj.Name = 'interface roughness scattering';
            obj.Lambda = device.material_system.interface_roughness.gamma;
            obj.Delta = device.material_system.interface_roughness.delta;
            obj.CBoffset = device.get_CBO * phys_const.e0;
        end

        function W_ikj = calc_state_rate(obj, i, j, k)
            %CALC_STATE_RATE calculates transition rate from state ik to
            % subband j according to https://doi.org/10.1063/1.4863665.
            %
            % Syntax:
            %   W_ikj = calc_state_rate(obj, i, j, k)
            %
            % Input:
            %   i (scalar): Inital subband.
            %   j (scalar): Final subband.
            %   k (scalar): Wavevector of the inital state.
            %
            % Output:
            %   W_ikj (scalar): Transition rate.

            % grid points for integration
            theta = linspace(0, pi, 100);

            q02 = 2 * obj.mEff(j) * (obj.E(i) - ...
                obj.E(j)) / phys_const.hbar^2;
            kPrime2 = k^2 * obj.mEff(j) / obj.mEff(i) + q02;

            if kPrime2 > 0
                q2 = k^2 + kPrime2 - 2 * k * sqrt(kPrime2) .* cos(theta);
                psi2 = abs(obj.psi(:, i).*conj(obj.psi(:, j))).^2;

                % interpolating to get values at interfaces
                sum_psi2 = sum(interp1(obj.z, psi2, ...
                    obj.PosInterf(2:end-1), 'spline'));

                integrand = exp(-obj.Lambda^2*q2/4);
                integrand(q2 < 0) = 0; % no scattering for q^2 < 0

                W_ikj = obj.mEff(j) * (obj.CBoffset * obj.Delta * ...
                    obj.Lambda)^2 * sum_psi2 * trapz(theta, ...
                    integrand) / phys_const.hbar^3;
            else
                % no scattering for kPrime2 < 0
                W_ikj = 0;
            end
        end

        function puredep = calulate_pure_dephasing_rate(obj, calc_rate)
            %CALCULATE_PURE_DEPHASING_RATES Calculates pure dephasing rates
            % for all subbands.
            %
            % Syntax:
            %   puredep = calulate_pure_dephasing_rate(obj)
            %   puredep = calulate_pure_dephasing_rate(obj, calc_rate)
            %
            % Input:
            %   calc_rate (logical): Flag specifying if k-dependent pure
            %     dephasing rates have to be recalculated.
            %
            % Output:
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
                            % Check if subband populations of
                            % |i> and |j> are equal
                            obj.pure_dephasing(ii, jj) = ...
                                obj.fermi_average_no_blocking(obj.E_F(i), ...
                                obj.E(i), obj.mEff(i), reshape( ...
                                obj.pure_dephasing_k(ii, jj, :), ...
                                size(obj.k)), obj.T_e(i));
                            obj.pure_dephasing(jj, ii) = ...
                                obj.fermi_average_no_blocking(obj.E_F(j), ...
                                obj.E(j), obj.mEff(j), reshape( ...
                                obj.pure_dephasing_k(jj, ii, :), ...
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

        function puredep = calulate_pure_dephasing_rate_parallel(obj, ...
                calc_rate)
            %CALCULATE_PURE_DEPHASING_PARALLEL Calculates pure dephasing
            % rates for all subbands using parallelization technique.
            %
            % Syntax:
            %   puredep = calulate_pure_dephasing_rate(obj)
            %   puredep = calulate_pure_dephasing_rate(obj, calc_rate)
            %
            % Input:
            %   calc_rate (logical): Flag specifying if k-dependent pure
            %     dephasing rates have to be recalculated.
            %
            % Output:
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
                            obj.E(j), obj.mEff(j), reshape( ...
                            obj.pure_dephasing_k(jj, ii, :), ...
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
            puredep = obj.pure_dephasing;
        end

        function gamma_ijk = calc_puredep(obj, i, j, k)
            %CALC_PUREDEP Calculates k-resolved pure dephasing rates
            % according to https://doi.org/10.1063/1.5005618.
            %
            % Syntax:
            %   puredep = calc_puredep(obj, i, j, k)
            %
            % Input:
            %   i (scalar): Inital subband.
            %   j (scalar): Final subband.
            %   k (scalar): Wavevector of inital state.
            %
            % Output:
            %   puredep (scalar): Rure dephasing rate.

            gamma_ijk = (obj.CBoffset * obj.Delta ...
                * obj.Lambda)^2 / (2 * phys_const.hbar^3);
            f = obj.mEff(i) * obj.psi(:, i).^4 - (obj.mEff(i) + ...
                obj.mEff(j)) .* obj.psi(:, i).^2 .* ...
                obj.psi(:, j).^2 + obj.mEff(j) * obj.psi(:, j).^4;
            x = 0.5 * obj.Lambda^2 * k^2;
            I = pi * exp(-x) * besseli(0, x);
            gamma_ijk = gamma_ijk * sum(interp1(obj.z, f, ...
                obj.PosInterf(1:end-1))) * I;
        end
    end
end
