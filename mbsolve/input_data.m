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

classdef input_data < handle
    % Contains collection of simulation results from the stationary carrier
    % transport solver, which is provided as input for the dynamic
    % maxwell-bloch solver.

    properties (SetAccess = private)
        project_name = "" % string: Project name.
        ind_wfs = [] % vector: Wavefunction indices of considered period.
        pairs_dipole = [] % cell-array: Dipole pairs for the simulations.
        pairs_tunneling = [] % cell-array: Tunnling pairs for the simulations.
        carr_dist % carrier_distribution-object: Contains occupations/distributions.
        deph % dephasing_rates-object: Contains dephasing rates.
        r_sc % scattering_rates-object: Contains scattering rates.
        device % device-object: Contains information about the QCL device.
        f_c % scalar: Center frequency of the laser emission spectrum.
    end

    properties (Dependent)
        dipole_matrix % matrix: Dipole matrix of the reduced system.
        hamiltonian % matrix: Hamiltonian of the reduced system.
        pure_dephasing_rates % matrix: Pure dephasing rates of the reduced system.
        dephasing_rates % matrix: Dephasing rates of the reduced system.
        transition_rates % matrix: Transition rates of the reduced system.
        scattering_rates_left % matrix: Scattering rates to the left period.
        scattering_rates_right % matrix: Scattering rates to the right period.
    end

    properties
        eigen_states % eignestates-object: Eigenstates of the specific QCL.
        lim_chi2 = 1e-10; % scalar: Limit in searching method for 2nd order nonlinear susceptibility.
    end

    methods
        function obj = input_data(name, d, i_wf, f_c, ...
                eigen, carr_dist, deph, sc)
            % Constructs an object of type input_data.
            %
            % Syntax:
            %   obj = input_data(name, d, i_wf, f_c, eigen, carr_dist, deph, sc)
            %
            % Input Arguments:
            %   name (scalar): Name of the project.
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   i_wf (vector): Wavefunction indices of considered period.
            %   f_c (scalar): Center frequency of laser emission.
            %   eigen (eigenstates_struct): Contains information about
            %     eigenenergies, wavefunctions and effective masses.
            %   carr_dist (carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.
            %   deph (dephasing_rates-object): Contains information about
            %     pure dephasing rates and lifetime broadening values.
            %   sc (scattering_rates-object): Contains information about
            %     the scattering rates.

            % Project name.
            obj.project_name = name;
            % Properties describing the quantum system.
            obj.eigen_states = eigen;
            obj.carr_dist = carr_dist;
            obj.deph = deph;
            obj.r_sc = sc;

            % Sort indices of the reduced system period vector.
            if isempty(i_wf)
                error(['Index array i_wf of reduced system ', ...
                    'must not be empty!'])
            end
            E_wf = obj.eigen_states.E(i_wf);
            i_wf = obj.sort_vector(i_wf, E_wf);
            obj.ind_wfs = obj.check_period(i_wf);
            % Set device property.
            obj.device = d;
            % Set center frequency.
            obj.f_c = f_c;
        end

        function add_dipole_pair(obj, ind_i, ind_j)
            % Add a dipole pair with state indices ind_i, ind_j.
            %
            % Syntax:
            %   add_dipole_pair(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_i (scalar): Index of initial state.
            %   ind_j (scalar): Index of final state.

            if (rem(ind_i, 1) ~= 0 || rem(ind_j, 1) ~= 0)
                error('Indices have to be of type integer!');
            end
            p_d = obj.pairs_dipole;
            % Check, whether indices are included in the predefined period.
            dipole_pair = obj.check_index([ind_i, ind_j], obj.ind_wfs);
            % Sort dipole pair with respect to the eigenenergy of the
            % corresponding eigenstates.
            E_pair = obj.eigen_states.E(dipole_pair);
            p_d{end+1} = obj.sort_vector(dipole_pair, E_pair);
            % Get unique dipole pair array.
            [obj.pairs_dipole, ~] = obj.get_unique_pairs(p_d);
            % Sort dipole pairs.
            obj.sort_dipole_pairs()
        end

        function add_tunnel_pair(obj, ind_i, ind_j)
            % Add a tunneling pair with indices ind_i, ind_j.
            %
            % Syntax:
            %   add_tunnel_pair(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_i (scalar): Wavefunction of initial state.
            %   ind_j (scalar): Wavefunction of final state.

            if (rem(ind_i, 1) ~= 0 || rem(ind_j, 1) ~= 0)
                error('Indices have to be of type integer!');
            end
            % Check, whether tight binding states are considered.
            if (isempty(obj.eigen_states.hamiltonian- ...
                    diag(diag(obj.eigen_states.hamiltonian))))
                warning(['Here extended states are considered,', ...
                    ' no tunneling included!']);
            end
            % Check, whether indices are included in the predefined period.
            tunnel_pair = obj.check_index([ind_i, ind_j], obj.ind_wfs);
            % Add tunneling pair to tunneling transition vector.
            p_t = obj.pairs_tunneling;
            p_t{end+1} = tunnel_pair;
            % Get unique tunneling pair array.
            [obj.pairs_tunneling, ~] = obj.get_unique_pairs(p_t);
        end

        function sort_ind_wfs(obj, sorting)
            % Sorts indices of the reduced system according to selected
            % sorting strategy.
            %
            % Syntax:
            %   sort_ind_wfs(obj, sorting)
            %
            % Input Arguments:
            %   sorting (char): Sorting strategy, sorts the states indices
            %     either by decending energy values (``sorting=descending``)
            %     or by period (``sorting=period``).

            i_wf = reshape(obj.ind_wfs, [], length(obj.ind_wfs));
            switch sorting
                case "descending"
                    % Indices are sorted with descending eigenenergies.
                    E_wf = obj.eigen_states.E(i_wf);
                    i_wf = obj.sort_vector(i_wf, E_wf);
                    ind_wf = obj.check_period(i_wf);
                case "period"
                    % Indices are sorted for each period seperately,
                    % starting with the period of the injector states.
                    i_wf = obj.check_period(i_wf);
                    period = floor((i_wf - 1)/length(i_wf));
                    p = unique(period);
                    p = sort(p, 'descend');
                    ind_wf = [];
                    for i = p
                        ind = i_wf(period == i);
                        E_wf = obj.eigen_states.E(ind);
                        ind_wf = [ind_wf, obj.sort_vector(ind, E_wf)];
                    end
                otherwise
                    error("Invalid value for input argument ''sorting''. "+ ...
                        "Valid values are ''descending'', ''period''.")
            end
            obj.ind_wfs = ind_wf;
        end

        function u = get.dipole_matrix(obj)
            % Gets the dipole matrix for the mb simulation
            % with relevant dipole pairs.
            u = zeros(4*obj.eigen_states.num_wfs, ...
                4*obj.eigen_states.num_wfs);
            for i = 1:length(obj.pairs_dipole)
                ind_i = obj.pairs_dipole{1, i}(1);
                ind_j = obj.pairs_dipole{1, i}(2);
                u(ind_i, ind_j) = ...
                    obj.eigen_states.get_dipole_element(ind_i, ind_j);
                u(ind_j, ind_i) = u(ind_i, ind_j);
            end
            u = u(obj.ind_wfs, obj.ind_wfs);
        end

        function u = get_dipole_mb(obj)
            % Get dipole moment matrix and add dipole elements between
            % injector states and lower laser levels for mb.

            u = obj.dipole_matrix;
            d = obj.eigen_states.dipoles(obj.ind_wfs, obj.ind_wfs);
            for i = 1:length(obj.pairs_tunneling)
                ni = find(obj.ind_wfs == ...
                    obj.pairs_tunneling{1, i}(1));
                nj = find(obj.ind_wfs == ...
                    obj.pairs_tunneling{1, i}(2));
                tun_pair = [ni, nj];
                % find lower laser level
                for j = 1:length(obj.pairs_dipole)
                    nUL = find(obj.ind_wfs == ...
                        obj.pairs_dipole{1, j}(1));
                    nLL = find(obj.ind_wfs == ...
                        obj.pairs_dipole{1, j}(2));
                    nInj = tun_pair(tun_pair ~= nUL);
                    if length(nInj) == 1
                        dij = d(nInj, nLL);
                        u(nInj, nLL) = dij;
                        u(nLL, nInj) = dij;
                    end
                end
            end
        end

        function H = get.hamiltonian(obj)
            % Gets the Hamiltonian for the mb simulations
            % including the relevant anticrossing energies for the given
            % tunneling pairs.
            % Add eigenenergies to Hamiltonian.
            E = obj.eigen_states.E(obj.ind_wfs);
            % Shift energy values to reduce time discretization points.
            E_shift = (max(E) + min(E)) / 2;
            H = diag(E-E_shift);
            num_wf = 4 * obj.eigen_states.num_wfs;
            Ea = zeros(num_wf, num_wf);
            % Add anticrossing energies to Hamiltonian.
            for i = 1:length(obj.pairs_tunneling)
                ind_i = obj.pairs_tunneling{1, i}(1);
                ind_j = obj.pairs_tunneling{1, i}(2);
                Ea(ind_i, ind_j) = ...
                    obj.eigen_states.get_anticrossing_energy(ind_i, ind_j);
                Ea(ind_j, ind_i) = Ea(ind_i, ind_j);
            end
            H = (H + Ea(obj.ind_wfs, obj.ind_wfs)) * phys_const.e0;
        end

        function mat_pure_deph = get.pure_dephasing_rates(obj)
            % Gets the vector of pure dephasing rates
            % for the reduced system.
            n_wf = obj.eigen_states.num_wfs;
            mat_pure_deph = zeros(n_wf, n_wf);
            % Add dephasing rates of the considered coherences.
            % Tight binding theory
            for i = 1:length(obj.pairs_tunneling)
                ni = obj.pairs_tunneling{1, i}(1);
                nj = obj.pairs_tunneling{1, i}(2);
                ind_i = obj.ind_wfs == ni;
                ind_j = obj.ind_wfs == nj;
                mat_pure_deph(ind_i, ind_j) = ...
                    obj.deph.get_pure_dephasing_rate(ni, nj);
                mat_pure_deph(ind_j, ind_i) = mat_pure_deph(ind_i, ind_j);
            end
            % Dipole matrix elements
            for i = 1:length(obj.pairs_dipole)
                ni = obj.pairs_dipole{1, i}(1);
                nj = obj.pairs_dipole{1, i}(2);
                ind_i = obj.ind_wfs == ni;
                ind_j = obj.ind_wfs == nj;
                mat_pure_deph(ind_i, ind_j) = ...
                    obj.deph.get_pure_dephasing_rate(ni, nj);
                mat_pure_deph(ind_j, ind_i) = mat_pure_deph(ind_i, ind_j);
            end
        end

        function mat_deph = get.dephasing_rates(obj)
            % Gets the matrix with dephasing rates for the reduced system.
            n_wf = obj.eigen_states.num_wfs;
            mat_deph = zeros(n_wf, n_wf);
            % Add dephasing rates of the considered coherences.
            % Tight binding theory
            for i = 1:length(obj.pairs_tunneling)
                ni = obj.pairs_tunneling{1, i}(1);
                nj = obj.pairs_tunneling{1, i}(2);
                ind_i = obj.ind_wfs == ni;
                ind_j = obj.ind_wfs == nj;
                mat_deph(ind_i, ind_j) = ...
                    obj.deph.get_dephasing_rate(ni, nj);
                mat_deph(ind_j, ind_i) = mat_deph(ind_i, ind_j);
            end
            % Dipole matrix elements
            for i = 1:length(obj.pairs_dipole)
                ni = obj.pairs_dipole{1, i}(1);
                nj = obj.pairs_dipole{1, i}(2);
                ind_i = obj.ind_wfs == ni;
                ind_j = obj.ind_wfs == nj;
                mat_deph(ind_i, ind_j) = ...
                    obj.deph.get_dephasing_rate(ni, nj);
                mat_deph(ind_j, ind_i) = mat_deph(ind_i, ind_j);
            end
        end

        function deph_red = get_dephasing_rates_mb(obj)
            % Get dephasing rate matrix and add elements between
            % injector states and lower laser levels for mb.
            %
            % Syntax:
            %   deph_red = get_dephasing_rates_mb(obj)
            %
            % Output Arguments:
            %   deph_red (matrix): Dephasing rates for reduced system.

            deph_red = obj.dephasing_rates;
            deph_tot = obj.deph.get_dephasing_rate(obj.ind_wfs, obj.ind_wfs);
            for i = 1:length(obj.pairs_tunneling)
                ni = find(obj.ind_wfs == ...
                    obj.pairs_tunneling{1, i}(1));
                nj = find(obj.ind_wfs == ...
                    obj.pairs_tunneling{1, i}(2));
                tun_pair = [ni, nj];
                % find lower laser level
                for j = 1:length(obj.pairs_dipole)
                    nUL = find(obj.ind_wfs == ...
                        obj.pairs_dipole{1, j}(1));
                    nLL = find(obj.ind_wfs == ...
                        obj.pairs_dipole{1, j}(2));
                    nInj = tun_pair(tun_pair ~= nUL);
                    if length(nInj) == 1
                        gamma_ij = deph_tot(nInj, nLL);
                        deph_red(nInj, nLL) = gamma_ij;
                        deph_red(nLL, nInj) = gamma_ij;
                    end
                end
            end
        end

        function m_trans_rates = get.transition_rates(obj)
            % Gets the matrix of transition rates for the reduced system.
            m_trans_rates = obj.r_sc.get_scattering_matrix('middle', ...
                obj.ind_wfs);
            m_rates_tun = ...
                obj.r_sc.get_scattering_matrix_id('tunneling', ...
                'middle', obj.ind_wfs);
            % Remove tunneling from rate,
            % if tunneling is considered in Hamiltonian.
            for i = 1:length(obj.pairs_tunneling)
                ind_i = (obj.ind_wfs == obj.pairs_tunneling{1, i}(1));
                ind_j = (obj.ind_wfs == obj.pairs_tunneling{1, i}(2));
                m_trans_rates(ind_i, ind_j) = ...
                    m_trans_rates(ind_i, ind_j) ...
                    -m_rates_tun(ind_i, ind_j);
                m_trans_rates(ind_j, ind_i) = ...
                    m_trans_rates(ind_j, ind_i) ...
                    -m_rates_tun(ind_j, ind_i);
                m_trans_rates(ind_i, ind_i) = ...
                    m_trans_rates(ind_i, ind_i) ...
                    +m_rates_tun(ind_i, ind_j);
                m_trans_rates(ind_j, ind_j) = ...
                    m_trans_rates(ind_j, ind_j) ...
                    +m_rates_tun(ind_j, ind_i);
            end
            m_trans_rates = m_trans_rates + obj.scattering_rates_left ...
                +obj.scattering_rates_right;
        end

        function m_scatt_rates_l = get.scattering_rates_left(obj)
            % Gets the matrix of scattering rates to the left period of the
            % reduced system.
            m_scatt_rates_l = obj.r_sc.get_scattering_matrix('left', ...
                obj.ind_wfs);
        end

        function m_scatt_rates_r = get.scattering_rates_right(obj)
            % Gets the matrix of scattering rates to the right period
            %of the reduced system.
            m_scatt_rates_r = obj.r_sc.get_scattering_matrix('right', ...
                obj.ind_wfs);
        end

        function [pairs_opt, num_opt] = find_optical_trans(obj, f_opt, num_opt)
            % Find the states of the optical transitions for a given
            % frequency by calculating the gain contributions of the
            % possible transition pairs for the give optical frequency.
            %
            % Syntax:
            %   [pairs_opt, num_opt] = find_optical_trans(obj, f_opt)
            %   [pairs_opt, num_opt] = find_optical_trans(obj, f_opt, num_opt)
            %
            % Input Arguments:
            %   f_opt (scalar): Optical frequency [Hz].
            %   num_opt (scalar): Number of optical transitions. If the
            %     number of optical transitions is not specified all
            %     transitions for which the gain @f_opt is greater than 0
            %     are considered.
            %
            % Output Arguments:
            %   pairs_opt (cell-array): Array containing the pairs of
            %     states of all considered optical transitions.
            %   num_opt (scalar): Number of optical transitions which were
            %     found (if not provided as input).

            % Calculate gain of all transitions in the given period of the
            % reduced system.
            gain = zeros(4*length(obj.ind_wfs), 4*length(obj.ind_wfs));
            for i = 1:length(obj.ind_wfs)
                for j = (i + 1):length(obj.ind_wfs)
                    ind_i = obj.ind_wfs(i);
                    ind_j = obj.ind_wfs(j);
                    occ_i = ...
                        obj.carr_dist.get_occupation(ind_i);
                    occ_j = ...
                        obj.carr_dist.get_occupation(ind_j);
                    if (occ_i > occ_j)
                        % Calculate gain for the given level pair and
                        % optical frequency.
                        gain(ind_i, ind_j) = ...
                            obj.get_gain(ind_i, ind_j, f_opt);
                    end
                end
            end
            % Determine number of possible optical transitions, if not
            % specified as input argument.
            if (nargin < 3)
                num_opt = nnz(gain);
            end
            % Get pairs of indices for the given number of
            % optical transitions at the investigated optical frequency.
            [sortedGain, sortedInds] = sort(gain(:), 'descend');
            top_nt = sortedInds(1:num_opt);
            [ind_top_i, ind_top_j] = ind2sub(size(gain), top_nt);
            pairs_opt = {};
            for i = 1:num_opt
                if sortedGain(i) == 0
                    warning("Gain peak is smaller than 0!")
                    pairs_opt{end+1} = [];
                else
                    pairs_opt{end+1} = [ind_top_i(i), ind_top_j(i)];
                end
            end
        end

        function sort_dipole_pairs(obj)
            % Sort dipole pairs according to their respective gain
            % contribution, i.e. strongest optical transition comes first.
            %
            % Syntax:
            %   sort_dipole_pairs(obj)

            dipoles = obj.pairs_dipole;
            g = zeros(length(dipoles), 1);
            for i = 1:length(dipoles)
                ind_i = dipoles{i}(1);
                ind_f = dipoles{i}(2);
                g(i) = obj.get_gain(ind_i, ind_f, obj.f_c);
            end
            [~, idx] = sort(g, 'descend');
            obj.pairs_dipole = dipoles(idx);
        end

        function [chi2, pair_DFG_ind, triplet_DFG_ind] = ...
                find_triplet_DFG(obj, num_triplets, f_IR_1, f_IR_2, ...
                f_THz, carr_dens)
            % Find level triplets for DFG generation by calculating the
            % second-order nonlinear susceptibility for the given
            % optical frequencies.
            %
            % Syntax:
            %   [chi2, pair_DFG_ind, triplet_DFG_ind] = find_triplet_DFG(obj, num_triplets, f_IR_1, f_IR_2, f_THz, carr_dens)
            %
            % Input Arguments:
            %   num_triplets (scalar): Number of triples to search for.
            %   f_IR_1 (scalar): Frequency of first IR transition [Hz].
            %   f_IR_2 (scalar): Frequency of second IR transition [Hz].
            %   f_THz (scalar): Generated THz frequency in DFG process [Hz].
            %   carr_dens (scalar): Charge carrier density in one QCL
            %     period [1/m^3].
            %
            % Output Arguments:
            %   chi2 (scalar | vector): Second order nonlinear
            %     susceptibility.
            %   pair_DFG_ind (cell-array):
            %   triplet_DFG_ind (cell-array):

            % Calculate averaged mid-IR frequencies.
            f_IR = (f_IR_1 + f_IR_2) / 2;
            % Get pairs for optical transition.
            chi2 = [];
            pair_DFG = {};
            triplet_DFG = {};
            [pairs_opt, num_opt] = obj.find_optical_trans(f_IR);
            for i = 1:num_opt
                for j = i + 1:num_opt
                    [chi2_new, trans_new] = ...
                        obj.get_suscept_2( ...
                        pairs_opt{i}, ...
                        pairs_opt{j}, f_IR_1, f_IR_2, ...
                        f_THz, carr_dens, 0);
                    if (chi2_new ~= 0)
                        chi2(end+1) = chi2_new;
                        pair_DFG{end+1} = trans_new;
                        triplet_DFG{end+1} = ...
                            {trans_new, pairs_opt{i}, pairs_opt{j}};
                    end
                end
            end
            if (length(chi2) < num_triplets)
                error(append('Cannot find the ', ...
                    num2str(num_triplets), ...
                    ' wanted DFG triplets.'))
            end
            %             [pair_DFG, idx_pair] = ...
            %                 obj.get_unique_pairs(pair_DFG);
            %             chi2 = chi2(idx_pair);
            %             triplet_DFG = triplet_DFG(idx_pair);
            [~, ind_chi2] = sort(abs(chi2), 'descend');
            chi2 = chi2(ind_chi2(1:num_triplets));
            %             pair_DFG(ind_chi2(num_triplets+1:end)) = [];
            pair_DFG_ind = pair_DFG(ind_chi2(1:num_triplets));
            triplet_DFG_ind = triplet_DFG(ind_chi2(1:num_triplets));
            %             triplet_DFG(ind_chi2(num_triplets+1:end)) = [];
        end

        function pairs_tun = find_tunneling_trans(obj, num_t)
            % Find tunneling transitions by calculating the tunneling rates.
            %
            % Syntax:
            %   pairs_tun = find_tunneling_trans(obj, num_t)
            %
            % Input Arguments:
            %   num_t (scalar): Number of tunneling transions.
            %
            % Output Arguments:
            %    pairs_tun (cell-array): Array containing the pairs of
            %     states of the tunneling transitions.

            % Calculate tunneling rates of all transitions
            % in the given period of the reduced system.
            rates_tun = ...
                zeros(4*length(obj.ind_wfs), 4*length(obj.ind_wfs));
            for i = 1:length(obj.ind_wfs)
                for j = i:length(obj.ind_wfs)
                    ind_i = obj.ind_wfs(i);
                    ind_j = obj.ind_wfs(j);
                    resonance_energy_ij = ...
                        obj.eigen_states.get_resonance_energy(ind_i, ...
                        ind_j);
                    if (resonance_energy_ij ~= 0)
                        rates_tun(ind_i, ind_j) = ...
                            obj.get_tunneling_rate(ind_i, ind_j);
                    end
                end
            end
            % Get pairs of indices for the given number of
            % tunneling transitions.
            [~, sortedInds] = sort(rates_tun(:), 'descend');
            top_nt = sortedInds(1:num_t);
            [ind_top_i, ind_top_j] = ind2sub(size(rates_tun), top_nt);
            pairs_tun = {};
            for i = 1:num_t
                pairs_tun{end+1} = [ind_top_i(i), ind_top_j(i)];
            end
        end

        function r_tun = get_tunneling_rate(obj, ind_i, ind_j)
            % Gets tunneling rate for a given level pair.
            %
            % Syntax:
            %   r_tun = get_tunneling_rate(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_i (scalar): Initial state index.
            %   ind_j (scalar): Final state index.
            %
            % Output Arguments:
            %   r_tun (scalar): Tunneling rate.

            % Resonance frequency.
            omega_res_ij = ...
                obj.eigen_states.get_resonance_freq(ind_i, ...
                ind_j);
            % Dephasing rate.
            deph_rate_ij = ...
                obj.deph.get_dephasing_rate(ind_i, ind_j);
            % Rabi frequency.
            omega_rabi_ij = obj.eigen_states.get_rabi_freq(ind_i, ind_j);
            % Tunneling rate.
            r_tun = ...
                2 * omega_rabi_ij^2 * deph_rate_ij ...
                / (deph_rate_ij^2 + omega_res_ij^2);
        end

        function [chi2, trans_THz] = ...
                get_suscept_2(obj, trans_1, trans_2, ...
                f_IR_1, f_IR_2, f_THz, carr_dens, sw_warning)
            % Gets second order nonlinear susceptibility chi_2 for a given
            % triplet of states.
            %
            % Syntax:
            %   [chi2, trans_THz] = get_suscept_2(obj, trans_1, trans_2, f_IR_1, f_IR_2, f_THz, carr_dens)
            %   [chi2, trans_THz] = get_suscept_2(obj, trans_1, trans_2, f_IR_1, f_IR_2, f_THz, carr_dens, sw_warning)
            %
            % Input Arguments:
            %   trans_1 (vector): Contains the indices of the two states of
            %     the first transition.
            %   trans_2 (vector): Contains the indices of the two states of
            %     the second transition.
            %   f_IR_1 (scalar): Frequency of first IR transition [Hz].
            %   f_IR_2 (scalar): Frequency of second IR transition [Hz].
            %   f_THz (scalar): Generated THz frequency in DFG process [Hz].
            %   carr_dens (scalar): Charge carrier density in one QCL
            %     period [1/m^3].
            %   sw_warning (logical): Specifies if warning should be
            %     displayed when selected transitions don't form a DFG
            %     triplet.
            %
            % Output Arguments:
            %   chi2 (scalar): Second order nonlinear susceptibility for
            %     the selected DFG process.
            %   trans_THz (vector): Contains the indices of the two states
            %     of the THz transition.

            if (nargin < 8)
                sw_warning = true;
            end
            % Find THz transition.
            if (trans_1(1) == trans_2(1))
                E_trans(1) = obj.eigen_states.get_E_i(trans_1(2));
                E_trans(2) = obj.eigen_states.get_E_i(trans_2(2));
                trans_THz = ...
                    obj.sort_vector([trans_1(2), trans_2(2)], E_trans);
                trans_IR_1 = [trans_1(1), trans_THz(2)];
                trans_IR_2 = [trans_1(1), trans_THz(1)];
                pref_chi = 1;
            elseif (trans_1(2) == trans_2(2))
                E_trans(1) = obj.eigen_states.get_E_i(trans_1(1));
                E_trans(2) = obj.eigen_states.get_E_i(trans_2(1));
                trans_THz = ...
                    obj.sort_vector([trans_1(1), trans_2(1)], E_trans);
                trans_IR_1 = [trans_THz(1), trans_1(2)];
                trans_IR_2 = [trans_THz(2), trans_1(2)];
                pref_chi = -1;
            else
                if (sw_warning == true)
                    warning(append('The input optical transition ', ...
                        'does not', ' form', ' a DFG triplet.'))
                end
                chi2 = 0;
                trans_THz = [];
                return
            end
            % Calculate population inversion between upper state i and
            % lower state j.
            inv_IR_1 = carr_dens ...
                * (obj.carr_dist.get_occupation(trans_IR_1(1)) ...
                -obj.carr_dist.get_occupation(trans_IR_1(2)));
            inv_IR_2 = carr_dens ...
                * (obj.carr_dist.get_occupation(trans_IR_2(1)) ...
                -obj.carr_dist.get_occupation(trans_IR_2(2)));
            % Dipole moments.
            d_IR_1 = obj.eigen_states.get_dipole_element(trans_IR_1(1), ...
                trans_IR_1(2));
            d_IR_2 = obj.eigen_states.get_dipole_element(trans_IR_2(1), ...
                trans_IR_2(2));
            d_THz = obj.eigen_states.get_dipole_element(trans_THz(1), ...
                trans_THz(2));
            % Frequencies.
            omega_IR_1 = (obj.eigen_states.get_E_i(trans_IR_1(1)) ...
                -obj.eigen_states.get_E_i(trans_IR_1(2))) ...
                * phys_const.e0 / phys_const.hbar;
            omega_IR_2 = (obj.eigen_states.get_E_i(trans_IR_2(1)) ...
                -obj.eigen_states.get_E_i(trans_IR_2(2))) ...
                * phys_const.e0 / phys_const.hbar;
            omega_THz = (obj.eigen_states.get_E_i(trans_THz(1)) ...
                -obj.eigen_states.get_E_i(trans_THz(2))) ...
                * phys_const.e0 / phys_const.hbar;
            % Lifetime broadening.
            deph_IR_1 = ...
                obj.deph.get_dephasing_rate(trans_IR_1(1), trans_IR_1(2));
            deph_IR_2 = ...
                obj.deph.get_dephasing_rate(trans_IR_2(1), trans_IR_2(2));
            deph_THz = ...
                obj.deph.get_dephasing_rate(trans_THz(1), trans_THz(2));
            % Calculate non-linear susceptibility with Eq. 5,
            % Fujita 2018 et al., https://doi.org/10.1515/nanoph-2018-0093.
            chi2 = ...
                pref_chi / (phys_const.hbar^2 * phys_const.eps0) ...
                * (d_IR_1 * d_IR_2 * d_THz) ...
                / (2 * pi * f_THz - omega_THz + 1i * deph_THz) ...
                * (inv_IR_1 / (2 * pi * f_IR_1 - omega_IR_1 ...
                +1i * deph_IR_1) + inv_IR_2 / (omega_IR_2 ...
                -2 * pi * f_IR_2 + 1i * deph_IR_2));
        end

        function m_scat_tun = get_tot_scattering(obj)
            % Gets total scattering matrix for the reduced period including
            % the tunneling rates.
            %
            % Syntax:
            %   m_scat_tun = get_tot_scattering(obj)
            %
            % Output Arguments:
            %   m_scat_tun (matrix): Scattering rates.

            m_trans_rates = obj.transition_rates;
            % Correct diagonal elements to have diagonal element as sum of
            % all in and outscattering rates. (without tunneling rates)
            m_trans_rates = m_trans_rates - diag(diag(m_trans_rates));
            m_sum_offdiag = diag(-sum(m_trans_rates, 2));
            m_trans_rates = m_trans_rates + m_sum_offdiag;

            % Adds tunneling rates through off-diagonal
            % elements of the hamiltonian.
            m_scat_tun = m_trans_rates;
            for i = 1:obj.eigen_states.num_wfs
                for j = 1:obj.eigen_states.num_wfs
                    if obj.hamiltonian(i, j) ~= 0 && j ~= i
                        m_scat_tun(i, j) = m_scat_tun(i, j) + ...
                            2 * obj.hamiltonian(i, j)^2 / ...
                            phys_const.hbar^2 ...
                            * obj.dephasing_rates(i, j) ...
                            / ((obj.hamiltonian(i, i) ...
                            -obj.hamiltonian(j, j))^2 ...
                            / phys_const.hbar^2 ...
                            +obj.dephasing_rates(i, j)^2);
                    end
                end
            end

            % Diagonal element as sum of all in and outscattering rates.
            m_scat_tun = m_scat_tun - diag(diag(m_scat_tun));
            m_sum_offdiag = diag(-sum(m_scat_tun, 2));
            m_scat_tun = m_scat_tun + m_sum_offdiag;
        end

        function m_tot_scat_int = get_tot_scattering_int(obj, I)
            % Gets total scattering matrix for the reduced period including
            % optical rates with respect to the given optical intensity I.
            %
            % Syntax:
            %   m_tot_scat_int = get_tot_scattering_int(obj, I)
            %
            % Output Arguments:
            %   m_tot_scat_int (matrix): Scattering rates.

            % Gets total scattering matrix including tunneling.
            m_tot_scat_int = obj.get_tot_scattering;

            for i = 1:obj.eigen_states.num_wfs
                for j = 1:obj.eigen_states.num_wfs
                    if (obj.dipole_matrix(i, j) ~= 0)
                        m_tot_scat_int(i, j) = ...
                            m_tot_scat_int(i, j) ...
                            +obj.dipole_matrix(i, j)^2 / ...
                            phys_const.hbar^2 ...
                            * obj.dephasing_rates(i, j) / ...
                            (obj.dephasing_rates(i, j)^2 ...
                            +(2 * pi * obj.f_c ...
                            -abs(obj.hamiltonian(i, i) ...
                            -obj.hamiltonian(j, j)) ...
                            / phys_const.hbar)^2) ...
                            / phys_const.eps0 / phys_const.c0 / ...
                            obj.device.n_eff * I;
                        %.* (1-exp(-hbar*(wmax-abs(Hamilton(i,i,s)-
                        % Hamilton(j,j,s))/hbar)/(kB*T_electron)))
                    end
                end
            end

            % Diagonal element as sum of all in and outscattering rates.
            m_tot_scat_int = m_tot_scat_int ...
                -diag(diag(m_tot_scat_int));
            m_sum_offdiag = diag(-sum(m_tot_scat_int, 2));
            m_tot_scat_int = m_tot_scat_int ...
                +m_sum_offdiag;
        end

        function occ = calc_occupations(obj, I)
            % Calculate level occupations including optical transitions
            % based on the optical intensity.
            %
            % Syntax:
            %   occ = calc_occupations(obj, I)
            %
            % Input Arguments:
            %   I (scalar): Optical intensity [W/m^2].
            %
            % Output Arguments:
            %   occ (vector): Normalized occupations in each subband [1/m^2].

            % Total scattering matrix.
            m_tot_scat_int = obj.get_tot_scattering_int(I);
            % Normalize matrix. Change of scattering units to 1/ps reduces
            % singularity of matrix.
            R0 = m_tot_scat_int' / 1e12;
            % In- and outscattering of levels is balanced.
            r0 = (obj.eigen_states.num_wfs) * 0;
            % Sum of occuaptions equals 1.
            R0(obj.eigen_states.num_wfs, :) = 1;
            r0(obj.eigen_states.num_wfs) = 1;
            % Solve equation system and calc level occupations.
            occ = (R0 \ r0')';
        end

        function [rates_tl, gam2] = calc_rates_two_level(obj, fmin, fmax, ...
                I, tau_g, pure_dep)
            % Calculate scattering rates and pure dephasin rate for a
            % two-level-representation of the full system. The two level
            % representation is obtained by fitting a Lorentzian lineshape
            % function to the total gain spectrum.
            %
            % Syntax:
            %   [rates_tl, gam2] = calc_rates_two_level(obj, fmin, fmax, I, tau_g, pure_dep)
            %   [rates_tl, gam2] = calc_rates_two_level(obj, fmin, fmax, [], [], [])
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of the gain spectrum.
            %   fmax (scalar): Maximum frequency of the gain spectrum.
            %   I (1e8 (Default) | scalar): Optical intensity for which the
            %     gain recovery time is calculated [W/m^2].
            %   tau_g (scalar): Gain recovery time used for initializing
            %     fit-parameter for the dephasing rate.
            %   pure_dep (scalar): Pure dephasing rate used for
            %     initializing fit-parameter for the dephasing rate.
            %
            % Output Arguments:
            %   rates_tl (matrix): 2-by-2 matrix containing the transition
            %     rates of the 2-level system.
            %   gam2 (scalar): Fitted dephasing rate of 2-level system.

            % States of lasing transition
            i_ind = obj.pairs_dipole{1, 1}(1); % index in full system
            f_ind = obj.pairs_dipole{1, 1}(2); % index in full system

            if isempty(I)
                I = 1e8; % light intensity
            end

            if isempty(tau_g)
                % calculate gain recovery time
                tau_g = obj.calc_gain_recovery_two_level(I);
            end
            gam1 = 1 / tau_g;

            if isempty(pure_dep)
                % calulate pure dephasing rate between the two levels
                pure_dep = obj.deph.get_pure_dephasing_rate(i_ind, f_ind);
            end

            % calculate gain
            gain = sim_gain(obj.device, obj.eigen_states, obj.deph, obj.carr_dist);
            g = gain.get_gain(fmin, fmax, 'total');
            f = (fmin + (fmax - fmin) * (1:100) / 100);
            omega = 2 * pi * f;
            [max_g, max_ind] = max(g);
            w21_start = 2 * pi * f(max_ind);

            % sample gain profile around the gain maximum
            num_samples = 5;
            for i = 1:length(g)
                if g(i) >= 0.3 * max_g
                    break
                end
            end
            ind_gain = i + round(2*(max_ind - i)/(num_samples - 1)) ...
                * (0:(num_samples - 1));
            xdata = omega(ind_gain);
            ydata = g(ind_gain);

            % fit lorentzian lineshape function to gain profile to extract
            % a single dephasing rate for the two level system
            gam2_start = gam1 / 2 + pure_dep;
            x0 = [max_g, gam2_start, w21_start];
            lorentz = @(x, xdata) x(1) * x(2)^2 ./ (x(2)^2 + (xdata - x(3)).^2);
            options = optimoptions('lsqcurvefit', ...
                'OptimalityTolerance', 1e-40, 'FunctionTolerance', 1e-40, ...
                'MaxFunctionEvaluations', 2e4, 'MaxIterations', 2e4);
            fit_result = lsqcurvefit(lorentz, x0, xdata, ...
                ydata, [], [], options);
            gam2 = fit_result(2);

            % plot
            y = fit_result(1) * fit_result(2)^2 ./ ...
                (fit_result(2)^2 + (omega - fit_result(3)).^2);
            figure;
            plot(f/1e12, g/100, 'LineWidth', 1.5, ...
                'DisplayName', 'Full System')
            hold on
            plot(omega/2/pi/1e12, y/100, 'LineWidth', 1.5, ...
                'DisplayName', 'Two Level')
            plot(xdata/2/pi/1e12, ydata/100, '*', ...
                'LineWidth', 1.5, 'DisplayName', 'Fitting Points')
            legend;
            xlabel('Frequency (THz)');
            ylabel('Gain (cm^{-1})');

            % Calculate transition rates between the upper laser level
            % (level 1) and lower laser level (level 2)
            Gamma = obj.device.waveguide.overlap_factor;
            d21 = obj.eigen_states.get_dipole_element(i_ind, f_ind);
            w21 = 2 * pi * obj.get_trans_freq(i_ind, f_ind);
            g_pref = Gamma * obj.device.dens_carrier * w21 * d21^2 / ...
                phys_const.eps0 / obj.device.n_eff / phys_const.c0 / ...
                phys_const.hbar / gam2;
            w_pump = max_g / g_pref;
            r12 = 0.5 * gam1 * (1 - w_pump);
            r21 = gam1 - r12;
            rates_tl = [0, r12; r21, 0];
        end

        function [A, tau_g] = plot_gain_recovery(obj, I, t)
            % Plots gain recovery for a given optical intensity I
            % and time t.
            %
            % Syntax:
            %   [A, tau_g] = plot_gain_recovery(obj, I)
            %   [A, tau_g] = plot_gain_recovery(obj, I, t)
            %
            % Input Arguments:
            %   I (scalar): Optical intensity [W/m^2].
            %   t (vector): Time vector [s]. The default time vector ranges
            %     from 0 to 30 ps in steps of 0.03 ps.
            %
            % Output Arguments:
            %   A (vector): Amplitudes of gain recovery times.
            %   tau_g (vector): Gain recovery times.

            if nargin > 2
                t_plot = t;
            else
                t_plot = (0:0.03:30) * 1e-12;
            end
            % Indices of laser levels in reduced system.
            [ULL_ind, LLL_ind] = find(triu(obj.dipole_matrix) ~= 0);

            tau_g = zeros(size(t_plot));
            % Get scattering rates without intensity dependency.
            m_scat_tun = obj.get_tot_scattering();
            % Get occupationions with intensity dependency.
            occ_I = obj.calc_occupations(I);

            for i = 1:length(ULL_ind)
                for t = 1:length(tau_g)
                    % Change in occupation over time.
                    occ_dt = expm(m_scat_tun'*t_plot(t)) ...
                        * occ_I';
                    % Nubering upper level.
                    nu = ULL_ind(i);
                    % Numbering lower level.
                    nl = LLL_ind(i);
                    % Calculate transition frequency of every optical
                    % transition.
                    w_ij(i) = (obj.hamiltonian(nu, nu) ...
                        -obj.hamiltonian(nl, nl)) / phys_const.hbar;
                    % If laser levels are reversed, change numbering and
                    % sign of the transition frequency.
                    if (w_ij(i) < 0)
                        [nu, nl] = deal(nl, nu);
                        w_ij(i) = -w_ij(i);
                    end
                    tau_g(t) = tau_g(t) + ...
                        obj.device.waveguide.overlap_factor ...
                        * 2 * pi * obj.f_c ...
                        / phys_const.hbar / phys_const.eps0 ...
                        / phys_const.c0 / obj.device.n_eff ...
                        * obj.device.dens_carrier ...
                        * obj.dipole_matrix(nu, nl)^2 ...
                        * (occ_dt(nu) - occ_dt(nl)) ...
                        .* obj.dephasing_rates(nu, nl) ...
                        ./ (obj.dephasing_rates(nu, nl)^2 ...
                        +(2 * pi * obj.f_c - w_ij(i)).^2);
                end
                % Calculate gain recovery time according to ref. 2019 paper
                % eq. 134
                fun = @(x, xdata)(x(1) * exp(-x(2)*xdata) + (1 - x(1)) ...
                    * exp(-x(3)*xdata));
                x = lsqcurvefit(fun, [1, 1, 0.1], t_plot/1e-12, ...
                    (tau_g(end) - tau_g)/(tau_g(end) - tau_g(1)));
                disp(['Amplitudes: ', num2str([x(1), 1 - x(1)]), ...
                    ', times: ', num2str([1 / x(2), 1 / x(3)]), ' ps']);
                figure;
                plot(t_plot*1e12, tau_g/max(tau_g), '-b');
                title('Gain-recovery');
                xlabel('time in ps')
            end

            % return gain recovery times and corresponding amplitudes
            A = [x(1), 1 - x(1)];
            tau_g = [1 / x(2), 1 / x(3)] * 1e-12;
        end

        function t_recovery = calc_gain_recovery_two_level(obj, I, t)
            % Calculates the gain recovery for a two level representation
            % of the full systems for a given optical intensity I and time
            % t.
            %
            % Syntax:
            %   t_recovery = calc_gain_recovery_two_level(obj, I)
            %   t_recovery = calc_gain_recovery_two_level(obj, I, t)
            %
            % Input Arguments:
            %   I (scalar): Optical intensity [W/m^2].
            %   t (vector): Time vector [s].
            %
            % Output Arguments:
            %   t_recovery (scalar): Gain recovery time [s].

            if nargin > 2
                tdata = t;
            else
                tdata = (0:0.03:30) * 1e-12;
            end
            % Indices of laser levels in reduced system.
            [ULL_ind, LLL_ind] = find(triu(obj.dipole_matrix) ~= 0);

            tau_g = zeros(size(tdata));
            % Get scattering rates without intensity dependency.
            m_scat_tun = obj.get_tot_scattering();
            % Get occupationions with intensity dependency.
            occ_I = obj.calc_occupations(I);

            for i = 1:length(ULL_ind)
                for t = 1:length(tau_g)
                    % Change in occupation over time.
                    occ_dt = expm(m_scat_tun'*tdata(t)) ...
                        * occ_I';
                    % Nubering of upper and lower level.
                    nULL = ULL_ind(i);
                    nLLL = LLL_ind(i);
                    % Calculate transition frequency of every optical
                    % transition.
                    w_ij = (obj.hamiltonian(nULL, nULL) ...
                        -obj.hamiltonian(nLLL, nLLL)) / phys_const.hbar;
                    % If laser levels are reversed, change numbering and
                    % sign of the transition frequency.
                    if (w_ij < 0)
                        [nULL, nLLL] = deal(nLLL, nULL);
                        w_ij = -w_ij;
                    end
                    tau_g(t) = tau_g(t) + ...
                        obj.device.waveguide.overlap_factor ...
                        * 2 * pi * obj.f_c ...
                        / phys_const.hbar / phys_const.eps0 ...
                        / phys_const.c0 / obj.device.n_eff ...
                        * obj.device.dens_carrier ...
                        * obj.dipole_matrix(nULL, nLLL)^2 ...
                        * (occ_dt(nULL) - occ_dt(nLLL)) ...
                        .* obj.dephasing_rates(nULL, nLLL) ...
                        ./ (obj.dephasing_rates(nULL, nLLL)^2 ...
                        +(2 * pi * obj.f_c - w_ij).^2);
                end
            end

            % sample time vector
            num_samples = 5;
            for i = 1:length(tau_g)
                if (tau_g(i) - tau_g(1)) / (max(tau_g) - tau_g(1)) >= 0.999
                    break
                end
            end
            index = 1 + round(i/num_samples) * (0:(num_samples - 1));
            xdata = tdata(index) / 1e-12;
            ydata = (tau_g(end) - tau_g(index)) / (tau_g(end) - tau_g(1));

            % Calculate gain recovery time
            fun = @(x, xvar)(exp(-x*xvar));
            options = optimoptions("lsqcurvefit", 'TolX', 1e-20, 'TolFun', 1e-20);
            x = lsqcurvefit(fun, 1, xdata, ydata, [], [], options);
            t_recovery = 1e-12 * 1 / x;
            disp(['Gain recovery time: ', num2str(1/x), ' ps']);

            % plot fit
            figure; hold on
            plot(tdata*1e12, tau_g/max(tau_g), '-b', ...
                'DisplayName', 'Multi level', 'LineWidth', 1.2);
            plot(tdata*1e12, (max(tau_g) - (max(tau_g) - tau_g(1)) * ...
                fun(x, tdata*1e12))/max(tau_g), ...
                'DisplayName', 'Fit', 'LineWidth', 1.2);
            plot(xdata, tau_g(index)/max(tau_g), '.', ...
                'DisplayName', 'Sample points', 'MarkerSize', 15);
            title('Gain-recovery');
            xlabel('time in ps');
            legend();
        end

        function freq = get_trans_freq(obj, ind_i, ind_j)
            % Get transition frequency between level i and j.
            %
            % Syntax:
            %   freq = get_trans_freq(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_i (scalar): Initial state index.
            %   ind_j (scalar): Final state index.
            %
            % Output Arguments:
            %   freq (scalar): Transition frequency.

            freq = (obj.eigen_states.get_E_i(ind_i) ...
                -obj.eigen_states.get_E_i(ind_j)) ...
                * phys_const.e0 / phys_const.h;
        end

        function g_ij = get_gain(obj, ind_i, ind_j, f, I)
            % Get gain for the transition from state i to j.
            %
            % Syntax:
            %   g_ij = get_gain(obj, ind_i, ind_j, f)
            %   g_ij = get_gain(obj, ind_i, ind_j, f, I)
            %
            % Input Arguments:
            %   ind_i (scalar): Initial state index (must be element of
            %     ind_wfs).
            %   ind_j (scalar): Initial state index (must be element of
            %     ind_wfs).
            %   f (scalar | vector): Frequency values for which the gain
            %     should be calculated.
            %   I (scalar): Optical intensity (If provided, the optical
            %      transitions are included for calculation of the
            %      occupations).
            %
            % Output Arguments:
            %   g_ij (scalar | vector): Gain spectrum [1/m].

            if ~max(ind_i == obj.ind_wfs) && ...
                    ~max(ind_j == obj.ind_wfs)
                str = append ...
                    ('Please choose one of the following wavefunctions:' ...
                    , newline, num2str(obj.ind_wfs, 2));
                error(str);
                return;
            end
            % Set the occupations based on the optical intensity of the
            % reduced system.
            if nargin > 4
                [~, i] = max(obj.ind_wfs == ind_i);
                [~, j] = max(obj.ind_wfs == ind_j);
                occ = obj.calc_occupations(I);
                occ_i = occ(i);
                occ_j = occ(j);
            else
                % Set occupation for the given indices i and j without
                % intensity.
                occ = obj.carr_dist.occupation;
                occ_i = ...
                    occ(mod(ind_i-1, obj.eigen_states.num_wfs)+1);
                occ_j = ...
                    occ(mod(ind_j-1, obj.eigen_states.num_wfs)+1);
            end
            % Set the transition frequency for the given indices.
            omega_ij = 2 * pi * obj.get_trans_freq(ind_i, ind_j);
            % Set the dipole element for the given indices.
            d_ij = ...
                obj.eigen_states.get_dipole_element(ind_i, ind_j);
            % Set the dephasing rates for the given indices.
            gamma_ij = ...
                obj.deph.get_dephasing_rate(ind_i, ind_j);
            % Calculate the lorentzian for the given indices.
            L_ij = gamma_ij ...
                ./ (gamma_ij^2 + (omega_ij - 2 * pi * f).^2);
            % Prefactor gain.
            pref_gain = obj.device.waveguide.overlap_factor / ...
                phys_const.hbar / phys_const.eps0 / phys_const.c0 / ...
                obj.device.n_eff * obj.device.dens_carrier;
            % Calculate the gain for the given indicees.
            g_ij = pref_gain * ...
                omega_ij * (occ_i - occ_j) * d_ij^2 * L_ij;
            if isempty(g_ij)
                disp("")
            end
        end

        function [f_spec, gain] = calc_gain(obj, ind, freq, I)
            % Calculate gain spectrum for the given optical transitions.
            %
            % Syntax:
            %   [f_spec, gain] = calc_gain(obj)
            %   [f_spec, gain] = calc_gain(obj, ind)
            %   [f_spec, gain] = calc_gain(obj, ind, freq, I)
            %
            % Input Arguments:
            %   ind (cell-array): Array containing the pairs of
            %     states of the considered optical transitions. (If not
            %     provided as input the optical transitions as specified by
            %     the property pairs_dipole is used).
            %   freq (vector): Frequency values of gain spectrum.
            %   I (scalar): Optical intensity (If provided, the optical
            %      transitions are included for calculation of the
            %      occupations).
            %
            % Output Arguments:
            %   f_spec (vector): Frequencies of gain spectrum [Hz].
            %   gain (vector): Gain spectrum [1/m].

            % Find trans_opt from the dipole pairs.
            trans_opt = obj.pairs_dipole;
            % Calculate and set the default frequency range f_spec
            % for the given trans_opt and transition frequencies.
            f = zeros(length(trans_opt), 1);
            for i = 1:length(trans_opt)
                % Get position of the original system in the reduced
                % system.
                [~, ind_i] = max(obj.ind_wfs == trans_opt{i}(1));
                [~, ind_j] = max(obj.ind_wfs == trans_opt{i}(2));
                % Calculate transition frequency for the given level pair.
                f(i) = abs(obj.hamiltonian(ind_i, ind_i) ...
                    -obj.hamiltonian(ind_j, ind_j)) / phys_const.h;
            end
            f_spec = (linspace(0.7*min(f), 1.3*max(f), 1000));
            % If passed as argument set trans_opt and f_spec.
            if nargin > 1
                trans_opt = ind;
                if nargin > 2
                    f_spec = freq;
                end
            end
            gain = zeros(size(f_spec));
            % For the considered number of optical transitions.
            for i = 1:length(trans_opt)
                % Only consider transitions within the reduced system.
                if ~max(trans_opt{i}(1) == obj.ind_wfs) && ...
                        ~max(trans_opt{i}(2) == obj.ind_wfs)
                    str = append ...
                        ('Please choose one of the following ', ...
                        'wavefunctions:', newline, num2str(obj.ind_wfs, 2));
                    error(str);
                    return;
                end
                % Get the gain of all considered optical transitions.
                if nargin > 3
                    % With occupations based on the optical intensity.
                    gain = gain + obj.get_gain(trans_opt{i}(1), ...
                        trans_opt{i}(2), f_spec, I);
                else
                    % Without occupations based on the optical intensity.
                    gain = gain + obj.get_gain(trans_opt{1, i}(1), ...
                        trans_opt{1, i}(2), f_spec);
                end
            end
        end

        function plot_gain(obj, ind, freq, I)
            % Plots gain over frequency.
            %
            % Syntax:
            %   plot_gain(obj)
            %   plot_gain(obj, ind, freq, I)
            %
            % Input Arguments:
            %   ind (cell-array): Array containing the pairs of
            %     states of the considered optical transitions. If not
            %     provided as input the optical transitions as specified by
            %     the property pairs_dipole is used).
            %   freq (vector): Frequencies of gain spectrum.
            %   I (scalar): Optical intensity (If provided, the optical
            %      transitions are included for calculation of the
            %      occupations).

            if (nargin > 1)
                [f_spec, gain] = obj.calc_gain(ind, freq, I);
            else
                [f_spec, gain] = obj.calc_gain();
            end
            % Plot the gain of every considered transition.
            plot(f_spec/1e12, gain/100, 'LineWidth', 1.5);
            xlabel('Frequency in THz');
            ylabel('Gain in cm^{-1}');
        end

        function [f_spec, gvd] = calc_gvd(obj, ind, freq, I)
            % Calculate group velocity dispersion (GVD) for selected
            % optical transitions.
            %
            % Syntax:
            %   [f_spec, gvd] = calc_gvd(obj, ind, freq, I)
            %
            % Input Arguments:
            %   ind (cell-array): Array containing the pairs of
            %     states of the considered optical transitions. (If not
            %     provided as input the optical transitions as specified by
            %     the property pairs_dipole is used).
            %   freq (vector): Frequencies of gain spectrum.
            %   I (scalar): Optical intensity [W/m^2] (If provided, the
            %     optical transitions are included for calculation of the
            %     occupations).
            %
            % Output Arguments:
            %   f_spec (vector): Frequency vector [Hz].
            %   gvd (vector): Group velocity dispersion.

            % Set default trans_opt.
            trans_opt = obj.pairs_dipole;
            f = zeros(length(trans_opt), 1);
            % Calculate and set the default frequency range f_spec with
            % the default trans_opt and the hamiltonian of the reduced
            % system.
            for i = 1:length(trans_opt)
                % Get position of the original system in the reduced
                % system.
                [~, ind_i] = max(obj.ind_wfs == trans_opt{i}(1));
                [~, ind_j] = max(obj.ind_wfs == trans_opt{i}(2));
                f(i) = abs(obj.hamiltonian(ind_i, ind_i) ...
                    -obj.hamiltonian(ind_j, ind_j)) / phys_const.h;
            end
            f_spec = (linspace(0.7*min(f), 1.3*max(f), 1000));
            % If passed as argument set trans_opt and f_spec.
            if nargin > 1
                trans_opt = ind;
                if nargin > 2
                    f_spec = freq;
                end
            end
            gvd = zeros(size(f_spec));
            % For the considered number of optical transitions.
            for k = 1:length(trans_opt)
                % Only consider transitions within the reduced system.
                if ~max(trans_opt{k}(1) == obj.ind_wfs) && ...
                        ~max(trans_opt{k}(2) == obj.ind_wfs)
                    str = append ...
                        ('Please choose one of the following ', ...
                        'wavefunctions:', newline, ...
                        num2str(obj.ind_wfs, 2));
                    error(str);
                    return;
                end
                % Set the transition frequency.
                omega_opt = 2 * pi * obj.get_trans_freq( ...
                    trans_opt{k}(1), trans_opt{k}(2));
                % Set the dephasing rate.
                gamma_opt = obj.deph.get_dephasing_rate( ...
                    trans_opt{k}(1), trans_opt{k}(2));
                % Set the dipole.
                d_opt = obj.eigen_states.get_dipole_element( ...
                    trans_opt{k}(1), trans_opt{k}(2));
                % Set the occupations based on the optical intensity of the
                % reduced.
                if nargin > 3
                    [~, i] = max(obj.ind_wfs == trans_opt{k}(1));
                    [~, j] = max(obj.ind_wfs == trans_opt{k}(2));
                    occ = obj.calc_occupations(I);
                    occ_i = occ(i);
                    occ_j = occ(j);
                else
                    % Set occupation of the given indicees i and j without
                    % intensity.
                    occ = obj.carr_dist.occupation;
                    occ_i = ...
                        occ(mod(trans_opt{k}(1)-1, ...
                        obj.eigen_states.num_wfs)+1);
                    occ_j = ...
                        occ(mod(trans_opt{k}(2)-1, ...
                        obj.eigen_states.num_wfs)+1);
                end
                % Calculate the delta omega.
                delta_omega = 2 * pi * f_spec - omega_opt;

                alpha_gain = obj.device.waveguide.overlap_factor ...
                    * 2 * pi * f_spec / phys_const.hbar ...
                    / phys_const.eps0 / phys_const.c0 ...
                    / obj.device.n_eff * obj.device.dens_carrier ...
                    * d_opt^2 * ((occ_i - occ_j) .* gamma_opt);
                % Calculate the group velocity dispersion (GVD)
                % for the considered number of optical transitions
                % analytical.
                gvd = gvd + (alpha_gain) * 0.5 / gamma_opt ...
                    .* (-6 * delta_omega ./ (gamma_opt^2 ...
                    +delta_omega.^2).^2 ...
                    +8 * delta_omega.^3 ./ ...
                    (gamma_opt^2 + delta_omega.^2).^3);
            end
        end

        function plot_gvd(obj, ind, freq, I)
            % Plots the group velocity dispersion over frequency.
            %
            % Syntax:
            %   plot_gvd(obj)
            %   plot_gvd(obj, ind, freq, I)
            %
            % Input Arguments:
            %   ind (cell-array): Array containing the pairs of
            %     states of the considered optical transitions. (If not
            %     provided as input the optical transitions as specified by
            %     the property pairs_dipole is used).
            %   freq (vector): Frequencies of gain spectrum.
            %   I (scalar): Optical intensity [W/m^2].

            if (nargin > 1)
                [f_spec, gvd] = obj.calc_gvd(ind, freq, I);
            else
                [f_spec, gvd] = obj.calc_gvd();
            end
            % Plot the GVD of every considered transition.
            plot(f_spec/1e12, gvd/1e-27, 'LineWidth', 1.5);
            xlabel('Frequency in Hz');
            ylabel('GVD in fs^{2}/m');
        end

        function j = calc_current(obj, I)
            % Calculates the current density based on the scattering rates
            % and occupations of the reduced system. The influence of the
            % optical field on the current density can be either inlcuded
            % or excluded.
            %
            % Syntax:
            %   j = calc_current(obj)
            %   j = calc_current(obj, I)
            %
            % Input Arguments:
            %   I (scalar): Optical intensity [W/m^2].
            %
            % Output Arguments:
            %   j (scalar): Current density [kA/cm^2].

            if nargin > 1
                occ = obj.calc_occupations(I);
            else
                occ = obj.carr_dist.occupation;
            end
            % Calculate current through the considered period.
            j = 0;
            n2D = obj.device.dens_sheet;
            for ind_i = 1:length(obj.ind_wfs)
                for ind_j = 1:length(obj.ind_wfs)
                    j = j + phys_const.e0 * n2D ...
                        * (obj.scattering_rates_left(ind_i, ind_j) ...
                        -obj.scattering_rates_right(ind_i, ind_j)) ...
                        * occ(mod(obj.ind_wfs(ind_i)-1, ...
                        length(obj.ind_wfs))+1);
                end
            end
            % Get current density in kA/cm^2.
            j = j * 1e-7;
        end
    end

    methods (Static, Access = private)
        function vector_s = sort_vector(vector, E_vector)
            % Sorts the elements in a vector, where the index with the
            % higher eigenenergy is brought to the first position.
            %
            % Syntax:
            %   vector_s = sort_vector(vector, E_vector)
            %
            % Input Arguments:
            %   vector (vector): Vector whose elements should be sorted.
            %   E_vector (vector): Energy vector for sorting.
            %
            % Output Arguments:
            %   vector_s (vector): Sorted vector.

            [~, ind_E] = sort(E_vector, 'descend');
            vector_s = vector(ind_E);
        end

        function [p_unique, idx] = get_unique_pairs(p_cell)
            % Removes repeated elements in the given cell array.
            %
            % Syntax:
            %   [p_unique, idx] = get_unique_pairs(p_cell)
            %
            % Input Arguments:
            %   p_cell (cell-array): Input cell-array.
            %
            % Output Arguments:
            %   p_unique (cell-array): Cell-array with no repetitions.
            %   idx (array): Indices of first occurrence of the elements in
            %     the cell-array.

            B = cellfun(@(x) num2str(x(:)'), p_cell, ...
                'UniformOutput', false);
            [~, idx] = unique(B);
            idx = sort(idx); % Keep correct order of elements
            p_unique = p_cell(idx);
        end

        function pair_p = check_index(pair, period)
            % Check, whether the indices of the added pair are included
            % in the period of the reduced system.
            %
            % Syntax:
            %   [pair_p] = check_index(pair, period)
            %
            % Input Arguments:
            %   pair (vector): Contains state indices of the pair.
            %   period (vector): State indices of the period of the reduced
            %     system.
            %
            % Output Arguments:
            %   pair_p (vector): Returns the state indices of the pair.

            pair_p = 0 * pair;
            n_wf = length(period);
            period_1 = mod(period-1, n_wf) + 1;
            for i = 1:length(pair)
                nx_1 = mod(pair(i)-1, n_wf) + 1;
                ind = period_1 == nx_1;
                pair_p(i) = period(ind);
            end
            if (pair ~= pair_p)
                error(['Given states are not in de predefined ', ...
                    'period ', 'for', ...
                    ' the reduced system of the mbsolve simulation!']);
            end
        end

        function period_corr = check_period(period)
            % Check, whether the period of the reduced system does not
            % spread over more than two spatial periods and adjust the
            % indices to those of the spatial periods 2 and 3.
            %
            % Syntax:
            %   period_corr = check_period(period)
            %
            % Input Arguments:
            %   period (vector): State indices of the period of the reduced
            %     system.
            %
            % Output Arguments:
            %   period_corr (vector): State indices of the period of the
            %     reduced system shifted to the central periods 2 and 3.

            n_wf = length(period);
            if (sum(period <= n_wf))
                period_corr = period + n_wf;
            elseif (period > 2 * n_wf)
                period_corr = period - n_wf;
            else
                period_corr = period;
            end
            if (~ismember(period_corr, n_wf+1:3*n_wf))
                error(['Period of the reduced system spreads over', ...
                    'more than two spatial periods!']);
            end
        end
    end
end
