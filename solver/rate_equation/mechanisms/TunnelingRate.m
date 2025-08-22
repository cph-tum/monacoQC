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

classdef TunnelingRate < FermiGoldenRule
    % Calculates transition rates due to incoherent tunneling across the
    % injection barrier between two adjacent periods. Note that this
    % scattering mechanism requires wavefunctions obtained from a tight-
    % binding description.

    properties
        deph_rates % matrix: Dephasing rates k-averaged [1/s].
        deph_rates_k % 3-d array: Dephasing rates k-resolved [1/s].
        V_scatter % matrix: Scattering potential for each period.
        index_period % vector: Period number for each subband.
        k_resolved % logical: Flag specifying if k-resolved dephasing rates are used.
    end

    methods (Access = public)
        function obj = TunnelingRate(eigen, ...
                device, scenario, cond, options)
            % Constructs an object of type TunnelingRate.
            %
            % Syntax:
            %   obj = TunnelingRate(eigen, device, scenario, cond)
            %   obj = TunnelingRate(eigen, device, scenario, cond, options)
            %
            % Input Arguments:
            %   eigen (eigenstates-object): Contains information about
            %      eigenenergies, wavefunctions and effective masses.
            %   device (device-object): Contains information about the
            %      structure/ geometry and materials of the QCL.
            %   scenario (scenario-object): Contains information about the
            %      specific scenario considered for the simulation.
            %   cond (conduction_band-object): Contains information about
            %      the conduction band profile. The influence of the
            %      electrostatic potential arising from the solution of
            %      Poisson's equation must not be included.
            %   options (cell-array): Array containing name-value pairs for
            %      changing values of default properties. Valid names are
            %      ``num_k``, ``Te`` and ``screening_model``.

            if nargin < 5
                options = {};
            else
                validStrings = ["num_k", "Te", "screening_model"];
                for l = 1:2:length(options)
                    validatestring(options{l}, validStrings);
                end
            end

            obj@FermiGoldenRule(eigen, device, scenario, options);

            obj.Name = "tunneling";
            obj.deph_rates = zeros(obj.num_states, obj.num_states);
            obj.deph_rates_k = zeros(obj.num_states, ...
                obj.num_states, length(obj.k));

            obj.index_period = zeros(obj.num_states, 1);
            obj.W = zeros(obj.num_states, obj.num_states, length(obj.k));
            obj.V_scatter = zeros(length(obj.z), device.num_periods);
            obj.k_resolved = false;

            % check if correct boundary conditions for Schrödinger
            % equation were chosen.
            if ~(scenario.basis_sp == "tb")
                error("TunnelingRate requires tight-binding states.")
            end

            % remove double values from conduction band vector
            [~, ind1] = unique(cond.zv, "first");
            [~, ind2] = unique(cond.zv, "last");
            for i = 1:length(ind1)
                if cond.Vh(ind1(i)) < cond.Vh(ind2(i))
                    ind1(i) = ind2(i);
                end
            end

            % unbiased conduction band vector
            V_unbiased = (cond.Vh(ind1) - cond.zv(ind1) .* ...
                scenario.V * 1e-5 * phys_const.e0);
            % potential of barrier material
            V_barrier = max(V_unbiased) * ones(1, length(obj.z));
            % get index of center of first tunneling barrier
            offset = round((device.tot_length - device.l_period * ...
                device.num_periods)/(2 * scenario.dz_sp));
            % number of points in each period
            num = device.l_period / scenario.dz_sp;

            % for each period compute scattering potential as
            % the difference between conduction band profile and
            % tight-binding potential profile
            num_states_periods = length(eigen.E) / scenario.num_wavefct;
            for i = 1:num_states_periods
                % create scattering potential for each period
                deltaV = V_unbiased - V_barrier;
                % indeces of period boundaries
                index_first = offset + num * (i - 1) + 1;
                index_last = offset + num * i;
                deltaV(index_first:index_last) = zeros(1, num);
                % save scattering potential
                obj.V_scatter(:, i) = reshape(deltaV, [], 1);
                % assign period number to each subband
                start = (i - 1) * eigen.num_wfs + 1;
                stop = i * eigen.num_wfs;
                obj.index_period(start:stop) = i + 1;
            end
        end

        function tau_inv = calculate(obj, calc_W)
            % Calculates transition rates between all subbands. Depending
            % on the property flag k_resolved, either k-dependend or
            % k-averaged dephasing rates are used for the calculation.
            %
            % Note:
            %   The tunneling rates are always re-calculated independent of
            %   the value of the input variable calc_W.
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


            if obj.k_resolved
                tau_inv = obj.calculate_k_resolved(calc_W);
            else
                tau_inv = obj.calculate_avg();
                obj.set_W(); % k-resolved rates
            end
        end

        function tau_inv = calculate_parallel(obj, calc_W)
            % Calculates transition rates between all subbands.
            %
            % Note:
            %   The computation of the tunneling rate is not parallelized,
            %   because it does not result in any significant speedup of
            %   the overall calculation.
            %
            % Note:
            %   The tunneling rates are always re-calculated independent of
            %   the value of the input variable calc_W.
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

            tau_inv = obj.calculate(calc_W);
        end

        function tau_inv = calculate_avg(obj)
            % Calculates k-independent tunneling rates between all subbands
            % spanning the tunneling barrier (i.e. only transitions between
            % subbands residing in different periods are considered).
            %
            % Syntax:
            %   tau_inv = calculate_avg(obj)
            %
            % Input Arguments:
            %   None
            %
            % Output Arguments:
            %   tau_inv (matrix): Matrix containing the transition rates
            %     between the subbands of the two central periods.

            tau_inv = zeros(obj.num_states);
            matrix_elem = zeros(obj.num_states);

            for ii = 1:obj.num_states
                for jj = 1:obj.num_states
                    i = obj.indices_i_j(ii);
                    j = obj.indices_i_j(jj);
                    % exclude transitions between states
                    % from the same period
                    if obj.index_period(j) ~= obj.index_period(i)
                        % asymmetric matrix elements
                        M_ij = obj.calc_matrix_element(i, j);
                        M_ij = M_ij .* obj.calc_matrix_element(j, i);
                        if M_ij < 0
                            M_ij = abs(M_ij);
                        end
                        % dephasing rates
                        gamma_ij = obj.deph_rates(ii, jj);
                        % transition frequency
                        omega_ij = (obj.E(i) - obj.E(j)) / phys_const.hbar;
                        % save calculations
                        matrix_elem(ii, jj) = M_ij;
                        tau_inv(ii, jj) = 2 * M_ij ./ ...
                            phys_const.hbar^2 .* gamma_ij ./ ...
                            (omega_ij.^2 + gamma_ij.^2);
                    end
                end
            end
            obj.tau_inv = tau_inv;
        end

        function tau_inv = calculate_k_resolved(obj, calc_W)
            % Calculates k-dependent tunneling rates between all subbands
            % spanning the tunneling barrier (i.e. only transitions between
            % subbands residing in different periods are considered).
            %
            % Note:
            %   The tunneling rates are always recalculated independent of
            %   value of the input variable calc_W.
            %
            % Syntax:
            %   tau_inv = calculate_k_resolved(obj, obj_W)
            %
            % Input Arguments:
            %   calc_W (logical): Flag determining if k-resolved transition
            %     rates have to be re-calculated.
            %
            % Output Arguments:
            %   tau_inv (matrix): Matrix containing the transition rates
            %     between the subbands of the two central periods.

            if calc_W
                obj.W = zeros(obj.num_states, ...
                    obj.num_states, length(obj.k));
            end
            tau_inv = zeros(obj.num_states, obj.num_states);

            for ii = 1:obj.num_states
                for jj = 1:obj.num_states
                    if ii ~= jj
                        i = obj.indices_i_j(ii);
                        j = obj.indices_i_j(jj);
                        % calculate k-resolved transition rates
                        if calc_W
                            for l = 1:length(obj.k)
                                obj.W(ii, jj, l) = calc_state_rate( ...
                                    obj, i, j, obj.k(l));
                            end
                        end
                        % calculate k-averaged transition rates
                        tau_inv(ii, jj) = obj.fermi_average( ...
                            obj.E_F(i), obj.E_F(j), obj.E(i), ...
                            obj.mEff(i), reshape(obj.W(ii, jj, :), ...
                            size(obj.k)), obj.T_e(i), obj.T_e(j));
                    end
                end
            end
            obj.tau_inv = tau_inv;
        end

        function W_ikj = calc_state_rate(obj, i, j, k)
            % Calculates transition rate from state ik to subband j
            % according to https://doi.org/10.1063/1.5005618.
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

            if obj.index_period(j) ~= obj.index_period(i)
                % get indices
                kk = find(obj.k == k);
                ii = i - min(obj.indices_i_j) + 1;
                jj = j - min(obj.indices_i_j) + 1;
                % asymmetric matrix elements
                M_ij = obj.calc_matrix_element(i, j);
                M_ij = M_ij .* obj.calc_matrix_element(j, i);
                if M_ij < 0
                    M_ij = abs(M_ij);
                end
                % dephasing rates
                gamma_ij = obj.deph_rates_k(ii, jj, kk);
                % transition frequency
                w_ik = obj.E(i) / phys_const.hbar + ...
                    phys_const.hbar .* k.^2 ./ obj.mEff(i);
                w_jk = obj.E(j) / phys_const.hbar + ...
                    phys_const.hbar .* k.^2 ./ obj.mEff(j);
                omega_ij = w_ik - w_jk;
                W_ikj = 2 * M_ij ./ phys_const.hbar^2 .* ...
                    gamma_ij ./ (omega_ij.^2 + gamma_ij.^2);
            else
                W_ikj = 0;
            end
        end

        function set_W(obj)
            % Sets k-resolved transition rates.
            %
            % Syntax:
            %   set_W(obj)

            obj.W = zeros(obj.num_states, obj.num_states, length(obj.k));
            for ii = 1:obj.num_states
                for jj = 1:obj.num_states
                    if ii ~= jj
                        i = obj.indices_i_j(ii);
                        j = obj.indices_i_j(jj);
                        for l = 1:length(obj.k)
                            obj.W(ii, jj, l) = calc_state_rate( ...
                                obj, i, j, obj.k(l));
                        end
                    end
                end
            end
        end

        function update(obj, ns, inv_lifetime, pure_dephasing)
            % Updates electron density, quasi Fermi-levels and dephasing
            % rates.
            %
            % Syntax:
            %   update(obj, ns, lvl_broadening, pure_dephasing)
            %
            % Input Arguments:
            %   ns (vector): Sheet densities of all four QCL-periods.
            %   inv_lifetime (vector | matrix): Lifetime broadening of
            %     each subband in the two central QCL-periods (1st
            %     dimension). The k-dependence can be specified
            %     as well (2nd dimension).
            %   pure_dephasing (matrix | 3-d array): Pure dephasing
            %     rates for all subbands in the two central QCL-eriods (1st/
            %     2nd dimension). The k-dependence can be specified
            %     as well (3rd dimension).

            obj.E_F = FermiGoldenRule.calc_fermi_energy(obj.E, ...
                obj.mEff, ns, obj.T_e);
            obj.electron_density = ns;
            obj.set_dephasing_rates(inv_lifetime, pure_dephasing);
        end
    end

    methods (Access = private)
        function M_ij = calc_matrix_element(obj, i, j)
            % Calculates the transition matrix element for a transition
            % from subband i to subband j.
            %
            % Syntax:
            %   M_ij = calc_matrix_element(obj, i, j)
            %
            % Input Arguments:
            %   i (scalar): Inital subband.
            %   j (scalar): Final subband.
            %
            % Output Arguments:
            %   M_ij (scalar): Transition matrix element <i|V_scatter|j>.

            % get index of period where state i resides
            period_i = obj.index_period(i);
            % evaluate matrix elements
            M_ij = trapz(obj.z*1e-10, conj(obj.psi(:, i)).* ...
                obj.V_scatter(:, period_i).*obj.psi(:, j));
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
                        % k-independent dephasing rates
                        obj.deph_rates(i, j) = 0.5 * (lvl_broadening(i) ...
                            +lvl_broadening(j)) + pure_dephasing(i, j);
                    elseif length(size(pure_dephasing)) == 3
                        % k-dependent dephasing rates
                        obj.deph_rates_k(i, j, :) = 0.5 * ...
                            (lvl_broadening(i, :) + lvl_broadening(j, :)) ...
                            +reshape(pure_dephasing(i, j, :), 1, []);
                    end
                end
            end
        end
    end
end
