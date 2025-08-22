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

classdef (Abstract) FermiGoldenRule < ScatteringRate
    % Base class for implementing scattering rate mechanisms based on
    % Fermi's golden rule.

    properties
        W % 3-d array: Transisition rates from state ik to subband j.
        E_F % vector: Quasi fermi-levels for each subband [J].
        electron_density % vector: 2D electron density [1/m^2].
        pure_dephasing_k % 3-d array: Pure dephasing rates k-resolved.
        pure_dephasing % matrix: Pure dephasing rates k-averaged.
        PosInterf % vector: Position of heterojunction interfaces.
        epsr_const % scalar: Relative permitivity of the qw-material.
        qs2 % scalar: Screening wavevector squared [1/m^2].
        l_period % scalar: Length of one QCL-period [m].
        screening_model = "static-lindhard" % char: Name of screening model.
    end

    methods
        function obj = FermiGoldenRule(eigen, device, scenario, options)
            % Constructs an object of type FermiGoldenRule.
            %
            % Syntax:
            %   obj = FermiGoldenRule(eigen, device, scenario)
            %   obj = FermiGoldenRule(eigen, device, scenario, options)
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

            obj@ScatteringRate(eigen, device, scenario);
            % overwrite default values by user defined options
            if nargin > 3
                obj.set_options(options);
            end

            obj.l_period = device.l_period * 1e-10;
            obj.PosInterf = device.int_pos;
            obj.epsr_const = device.layers{2}.material.eps_r;
            obj.W = zeros(obj.num_states, obj.num_states, length(obj.k));
            obj.qs2 = 0;

            % Calculate Fermi energy
            nE = length(obj.E);
            obj.electron_density = device.dens_sheet / ...
                eigen.num_wfs * ones(nE, 1);
            obj.E_F = FermiGoldenRule.calc_fermi_energy( ...
                obj.E, obj.mEff, obj.electron_density, obj.T_e);

            % Pure dephasing rates
            obj.pure_dephasing_k = zeros(obj.num_states, ...
                obj.num_states, length(obj.k));
            obj.pure_dephasing = zeros(obj.num_states, obj.num_states);
        end

        function inv_tau_ij = fermi_average(obj, EF_i, EF_j, E_i, ...
                mEffi, v_k, Te_i, Te_j)
            % Averages vector v_k over Fermi-Dirac distribution including
            % final state blocking. If EF_i-E_i << kT, the Boltzmann
            % distribution is used for averaging.
            %
            % Syntax:
            %   inv_tau_ij = fermi_average(obj, EF_i, EF_j, E_i, mEffi, v_k, Te_i, Te_j)
            %
            % Input Arguments:
            %   EF_i (scalar): Quasi Fermi-level of inital subband.
            %   EF_j (scalar): Quasi Fermi-level of final subband.
            %   E_i (scalar): Energy level of inital subband.
            %   mEffi (scalar): In-plane effective mass of subband i.
            %   v_k (scalar): Transistion rates from state ik to subband j.
            %   Te_i (scalar): Electron temperature of subband i.
            %   Te_j (scalar): Electron temperature of subband j.
            %
            % Output Arguments:
            %   inv_tau_ij (scalar): K-averaged transition rate from
            %     subband i to j.

            % energy of state ik
            kT = phys_const.kB * Te_i;
            Ekin = phys_const.hbar^2 * obj.k.^2 / (2 * mEffi);
            E_ik = E_i + Ekin;

            if exp((EF_i - E_i)/kT) < 1e-10
                % Boltzmann distribution
                fkin = phys_const.hbar^2 / (mEffi * kT) * exp(-Ekin/kT);
            else
                % Fermi-Dirac distribution
                ns_i = mEffi * kT / (pi * phys_const.hbar^2) * ...
                    ((EF_i - E_i) / kT + log(1+exp((E_i - EF_i)/kT)));
                fkin = FermiGoldenRule.fermi_dirac(E_ik, EF_i, Te_i) / ...
                    (pi * ns_i);
            end

            fkin2 = fkin .* (1 - FermiGoldenRule.fermi_dirac(E_ik+ ...
                get_energy_pauli_blocking(obj), EF_j, Te_j));
            inv_tau_ij = trapz(obj.k, fkin2.*v_k.*obj.k);
        end

        function inv_tau_ij = pop_inv_average(obj, EF_i, EF_j, E_i, E_j, ...
                mEff_i, mEff_j, v_k, Te_i, Te_j)
            % Averages vector v_k over k-resolved population inversion
            % between subband i and j, where the distribution in each
            % subband is given by the Fermi-Dirac distribution. This
            % averaging scheme is needed e.g. for calculating pure
            % dephasing rates.
            %
            % Syntax:
            %   inv_tau_ij = pop_inv_average(obj, EF_i, EF_j, E_i, mEffi, mEff_j, v_k, Te_i, Te_j)
            %
            % Input Arguments:
            %   EF_i (scalar): Quasi Fermi-level of inital subband.
            %   EF_j (scalar): Quasi Fermi-level of final subband.
            %   E_i (scalar): Energy level of inital subband.
            %   E_j (scalar): Energy level of final subband.
            %   mEffi (scalar): In-plane effective mass of subband i.
            %   mEffj (scalar): In-plane effective mass of subband j.
            %   v_k (scalar): Transistion rates from state ik to subband j.
            %   Te_i (scalar): Electron temperature of subband i.
            %   Te_j (scalar): Electron temperature of subband j.
            %
            % Output Arguments:
            %   inv_tau_ij (scalar): K-averaged transition rate from
            %     subband i to j.

            % energy of state ik
            Ekin_i = phys_const.hbar^2 * obj.k.^2 ./ (2 * mEff_i);
            E_ik = E_i + Ekin_i;

            % energy of state jk
            Ekin_j = phys_const.hbar^2 * obj.k.^2 ./ (2 * mEff_j);
            E_jk = E_j + Ekin_j;

            fkin = FermiGoldenRule.fermi_dirac(E_ik, EF_i, Te_i);
            fkin = fkin - FermiGoldenRule.fermi_dirac(E_jk, EF_j, Te_j);

            norm = trapz(obj.k, abs(fkin).*obj.k);

            inv_tau_ij = trapz(obj.k, abs(fkin).*v_k.*obj.k) / norm;
        end

        function inv_tau_ij = fermi_average_no_blocking(obj, EF_i, E_i, ...
                mEffi, v_k, Te_i)
            % Averages vector v_k over Fermi-Dirac distribution neglecting
            % final state blocking. If EF_i-E_i << kT, the Boltzmann
            % distribution is used for averaging.
            %
            % Syntax:
            %   inv_tau_ij = fermi_average(obj, EF_i, EF_j, E_i, mEffi, v_k, Te_i, Te_j)
            %
            % Input Arguments:
            %   EF_i (scalar): Quasi Fermi-level of inital subband.
            %   EF_j (scalar): Quasi Fermi-level of final subband.
            %   E_i (scalar): Energy level of inital subband.
            %   mEffi (scalar): In-plane effective mass of subband i.
            %   v_k (scalar): Transistion rates from state ik to subband j.
            %   Te_i (scalar): Electron temperature of subband i.
            %   Te_j (scalar): Electron temperature of subband j.
            %
            % Output Arguments:
            %   inv_tau_ij (scalar): K-averaged transition rate from
            %     subband i to j.

            % energy of state ik
            kT = phys_const.kB * Te_i;
            Ekin = phys_const.hbar^2 * obj.k.^2 / (2 * mEffi);
            E_ik = E_i + Ekin;

            if exp((EF_i - E_i)/kT) < 1e-10
                % Boltzmann distribution
                fkin = phys_const.hbar^2 / (mEffi * kT) * exp(-Ekin/kT);
            else
                % Fermi-Dirac distribution
                ns_i = mEffi * kT / (pi * phys_const.hbar^2) * ...
                    ((EF_i - E_i) / kT + log(1+exp((E_i - EF_i)/kT)));
                fkin = FermiGoldenRule.fermi_dirac(E_ik, EF_i, Te_i) / ...
                    (pi * ns_i);
            end

            inv_tau_ij = trapz(obj.k, fkin.*v_k.*obj.k);
        end

        function tau_inv = calculate(obj, calc_W)
            % Calculates transition rates between all considered subbands.
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
                obj.W = zeros(obj.num_states, obj.num_states, length(obj.k));
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
                                obj.W(ii, jj, l) = obj.calc_state_rate( ...
                                    i, j, obj.k(l));
                            end
                        end

                        % calculate k-averaged transition rates
                        tau_inv(ii, jj) = obj.fermi_average(obj.E_F(i), ...
                            obj.E_F(j), obj.E(i), obj.mEff(i), ...
                            reshape(obj.W(ii, jj, :), size(obj.k)), ...
                            obj.T_e(i), obj.T_e(j));
                    end
                end
            end
            obj.tau_inv = tau_inv;
        end

        function puredep = calulate_pure_dephasing_rate(obj, calc_rate)
            % Base method for the calculation of pure dephasing rates.
            % Initializes pure dephasing rates and k-resolved pure
            % dephasing rates properties.
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

            obj.pure_dephasing_k = zeros(size(obj.pure_dephasing_k));
            obj.pure_dephasing = zeros(size(obj.pure_dephasing));

            puredep = obj.pure_dephasing;
        end

        function puredep = calulate_pure_dephasing_rate_parallel(obj, ...
                calc_rate)
            % Base method for the calculation of pure dephasing rates
            % using parallelization technique. Initializs pure dephasing
            % rates and k-resolved pure dephasing rates properties.
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

            obj.pure_dephasing_k = zeros(size(obj.pure_dephasing_k));
            obj.pure_dephasing = zeros(size(obj.pure_dephasing));

            puredep = obj.pure_dephasing;
        end

        function gamma_ijk = calc_puredep(obj, i, j, k)
            % Base method for the calculation of k-resolved pure dephasing
            % rate from state ik to subband j. Returns zero if not
            % implemented for a certain subclass.
            %
            % Syntax:
            %   puredep = calc_puredep(obj, i, j, k)
            %
            % Input Arguments:
            %   i (scalar): Inital subband.
            %   j (scalar): Final subband.
            %   k (scalar): Wavevector of the inital state.
            %
            % Output Arguments:
            %   puredep (scalar): K-dependent pure dephasing rate.

            gamma_ijk = 0;
        end

        function tau_inv = calculate_parallel(obj, calc_W)
            % Calculates transition rates between all considered
            % subbands on parallel workers.
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

            % calculate k-resolved scattering rates
            if calc_W
                obj.W = zeros(nwf, nwf, length(obj.k));
                indexSpaceIJK = [nwf, nwf, length(obj.k)];
                numLoopsIJK = prod(indexSpaceIJK);
                W_ikj = zeros(1, numLoopsIJK);
                parfor idx = 1:numLoopsIJK
                    [ii, jj, l] = ind2sub(indexSpaceIJK, idx);
                    i = obj.indices_i_j(ii);
                    j = obj.indices_i_j(jj);
                    if i ~= j
                        W_ikj(idx) = calc_state_rate(obj, i, j, obj.k(l));
                    end
                end
                obj.W = reshape(W_ikj, indexSpaceIJK);
            end

            % calculate k-averaged scattering rates
            indexSpaceIJ = [nwf, nwf];
            numLoopsIJ = prod(indexSpaceIJ);
            tau_inv = zeros(numLoopsIJ, 1);
            parfor idx = 1:numLoopsIJ
                [ii, jj] = ind2sub(indexSpaceIJ, idx);
                i = obj.indices_i_j(ii);
                j = obj.indices_i_j(jj);
                if i ~= j
                    tau_inv(idx) = obj.fermi_average(obj.E_F(i), ...
                        obj.E_F(j), obj.E(i), obj.mEff(i), ...
                        reshape(obj.W(ii, jj, :), size(obj.k)), ...
                        obj.T_e(i), obj.T_e(j));
                end
            end
            tau_inv = reshape(tau_inv, indexSpaceIJ);
            obj.tau_inv = tau_inv;
        end

        function obj = update(obj, ns)
            %UPDATE Updates the sheet densities and quasi Fermi-levels.
            %
            % Syntax:
            %   update(obj, ns)
            %
            % Input Arguments:
            %   ns (vector): Sheet density for all four QCL-periods [1/m^2].

            obj.E_F = FermiGoldenRule.calc_fermi_energy( ...
                obj.E, obj.mEff, ns, obj.T_e);
            obj.electron_density = ns;
        end

        function E = get_energy_pauli_blocking(obj)
            %GET_ENERGY_PAULI_BLOCKING Returns energy used for calculation
            % with final state blocking due Pauli's exclusion principle.
            %
            % Syntax:
            %   E = get_energy_pauli_blocking(obj)
            %
            % Output Arguments:
            %   E (scalar): Exchanged energy due to the scattering process.
            %     Set to zero for elastic processes.

            E = 0;
        end

        function qs2 = calc_qs2(obj)
            % Calculates the screening wavevector based on the selected
            % screening model.
            %
            % Syntax:
            %   qs2 = calc_qs2(obj)
            %
            % Output Arguments:
            %   qs2 (scalar): Screening wavevector squared.

            switch obj.screening_model
                case "static-lindhard"
                    qs2 = obj.calc_qs2_static_lindhard();
                case "modified-single-subband"
                    qs2 = obj.calc_qs2_modified_single_subband();
                case "debye"
                    qs2 = obj.calc_qs2_debye();
                case "thomas-fermi"
                    qs2 = obj.calc_qs2_thomas_fermi();
            end
        end

        function set.screening_model(obj, name)
            % Set screening model.
            if ismember(name, ["static-lindhard", "debye", ...
                    "modified-single-subband", "thomas-fermi"])
                obj.screening_model = name;
            else
                error(['Selected screening model is invalid. ', ...
                    'Choose one of the following: ', ...
                    '''static-lindhard''', ', ', ...
                    '''modified-single-subband''', ', ', ...
                    '''debye''', ', ', ...
                    '''thomas-fermi''', '.'])
            end
        end

        function qs2 = calc_qs2_static_lindhard(obj)
            % Calculates screening wave vector squared qs^2 based on
            % static lindhard formula (eq. 7 in
            % https://doi.org/10.1063/1.4940192).
            %
            % Syntax:
            %   qs2 = calc_qs2_static_lindhard(obj)
            %
            % Output Arguments:
            %   qs2 (scalar): Screening wavevector squared.

            qs2 = 0;
            for ii = length(obj.indices_i_j) / 2
                i = obj.indices_i_j(ii);
                qs2 = qs2 + obj.mEff(i) * FermiGoldenRule.fermi_dirac( ...
                    obj.E(i), obj.E_F(i), obj.T_e(i));
            end
            eps = obj.epsr_const * phys_const.eps0;
            qs2 = qs2 * phys_const.e0^2 / ...
                (pi * phys_const.hbar^2 * eps * obj.l_period);
        end

        function qs2 = calc_qs2_modified_single_subband(obj)
            % Calculate screening wave vector squared qs^2 based on the
            % modified single subband model (eq. 95 in
            % https://doi.org/10.1063/1.4863665).
            %
            % Syntax:
            %   qs2 = calc_qs2_modified_single_subband(obj)
            %
            % Output Arguments:
            %   qs2 (scalar): Screening wavevector squared.

            qs = 0;
            for ii = length(obj.indices_i_j) / 2
                i = obj.indices_i_j(ii);
                qs = qs + obj.mEff(i) * FermiGoldenRule.fermi_dirac( ...
                    obj.E(i), obj.E_F(i), obj.T_e(i));
            end
            eps = obj.epsr_const * phys_const.eps0;
            qs2 = (qs * phys_const.e0^2 / ...
                (2 * eps * pi * phys_const.hbar^2))^2;
        end

        function qs2 = calc_qs2_debye(obj)
            % Calculates screening wave vector squared qs^2 based on the
            % Debye model (eq. 87 in https://doi.org/10.1063/1.4863665).
            %
            % Syntax:
            %   qs2 = calc_qs2_debye(obj)
            %
            % Output Arguments:
            %   qs2 (scalar): Screening wavevector squared.

            dens_e_3d = sum(obj.electron_density) / ...
                (4 * obj.l_period); % electron density in 1/m^3
            eps = obj.epsr_const * phys_const.eps0;
            EF = (3 * pi^2 * dens_e_3d * phys_const.hbar^3)^(2 / 3) / ...
                (2 * mean(obj.mEff));

            if phys_const.kB * obj.T < 2 * EF / 3
                % Use Thomas-Fermi model for small temperatures
                qs2 = 3 * phys_const.e0^2 * dens_e_3d / ...
                    (2 * eps * EF);
            else
                qs2 = phys_const.e0^2 * dens_e_3d / ...
                    (eps * phys_const.kB * obj.T);
            end
        end

        function qs2 = calc_qs2_thomas_fermi(obj)
            % Calculates screening wave vector squared qs^2 based on the
            % Thomas-Fermi model (eq. 5 in
            % https://doi.org/10.1063/1.4940192).
            %
            % Syntax:
            %   qs2 = calc_qs2_thomas_fermi(obj)
            %
            % Output Arguments:
            %   qs2 (scalar): Screening wavevector squared.

            dens_e_3d = sum(obj.electron_density) / ...
                (4 * obj.l_period); % electron density in 1/m^3
            eps = obj.epsr_const * phys_const.eps0;
            EF = (3 * pi^2 * dens_e_3d * phys_const.hbar^3)^(2 / 3) / ...
                (2 * mean(obj.mEff));
            qs2 = 3 * phys_const.e0^2 * dens_e_3d / (2 * eps * EF);
        end

        function set_options(obj, options)
            % Changes default values of some selected properties.
            %
            % Syntax:
            %   set_options(obj, options)
            %
            % Input Arguments:
            %   options (cell-array): Array containing name-value pairs for
            %     changing default properties. Valid names are ``num_k``,
            %     ``Te`` and ``screening_model``.

            set_options@ScatteringRate(obj, options);
            for i = 1:2:length(options)
                key = options{i};
                value = options{i+1};
                if key == "screening_model"
                    obj.screening_model = value;
                end
            end
        end
    end

    methods (Static)
        function E_F = calc_fermi_energy(E, mEff, ns, T_e)
            % Calclates quasi Fermi-level for each subband.
            %
            % Syntax:
            %   E_F = calc_fermi_energy(E, mEff, ns, T_e)
            %
            % Input Arguments:
            %   E (vector): Energies of subbands [J].
            %   mEff (vector): Effective masses of subbands [kg].
            %   ns (vector): Sheet densities of subbands [1/m^2].
            %   T_e (vector): Eelectron temperatures of subbands [K].
            %
            % Output Arguments:
            %   E_F (vector): Quasi Fermi-levels of all subbands [J].

            nE = length(E);
            E_F = zeros(nE, 1);
            for i = 1:nE
                kT = phys_const.kB * T_e(i);
                E_F(i) = kT * log(exp(ns(i)*pi*phys_const.hbar^2/ ...
                    (mEff(i) * kT))-1) + E(i);
            end
        end
        function f_ik = fermi_dirac(E, E_F, T)
            % Returns Fermi-Dirac distribution.
            %
            % Syntax:
            %   f_ik = fermi_dirac(E, E_F, T)
            %
            % Input Arguments:
            %   E (scalar | vector): Energy [J].
            %   E_F (scalar | vector): Quasi Fermi-level [J].
            %   T (scalar | vector): Temperature in [K].
            %
            % Output Arguments:
            %   f_ik (scalar | vector): Function value of Fermi-Dirac
            %     distribution.

            f_ik = 1 ./ (exp((E - E_F)./(phys_const.kB .* T)) + 1);
        end
    end

    methods (Abstract)
        W_ikj = calc_state_rate(obj, i, j, k); % transition rate ik->j
    end
end
