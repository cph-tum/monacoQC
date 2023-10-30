classdef mbsolve_sim < handle
    % mbsolve_sim Backend as input writer for the mbsolve tool.
    %
    properties (SetAccess = private)
        project_name = "" % Project name.
        ind_wfs = [] % Vector with wavefunctions of considered period.
        pairs_dipole = [] % Dipole pairs for the mbsolve simulation.
        pairs_tunneling = [] % Tunnling pairs for the mbsolve simulation.
        eigen_states % Eigenstates of the specific QCL device.
        carr_dist % Carrier distribution.
        deph % Dephasing rates.
        r_sc % Scattering rates.
        device % Device.
        f_c % Center frequency.
    end
    %
    properties (Dependent)
        dipole_matrix % Dipole matrix of the reduced system.
        hamiltonian % Hamiltonian of the reduced system.
        pure_dephasing_rates % Pure dephasing rates of the reduced system.
        dephasing_rates % Dephasing rates of the reduced system.
        transition_rates % Transition rates of the reduced system.
        scattering_rates_left % Scattering rates to the left period.
        scattering_rates_right % Scattering rates to the right period.
    end
    %
    properties
        lim_chi2 = 1e-10; % Limit for second order nonlinear susceptibility
        % in the implemented searching method.
    end
    %
    methods
        function obj = mbsolve_sim(name, d, i_wf, f_c, ...
                eigen, carr_dist, deph, sc)
            % Constructs mbsolve_sim.
            % Project name.
            obj.project_name = name;
            % Properties describing the quantum system.
            obj.eigen_states = eigen;
            obj.carr_dist = carr_dist;
            obj.deph = deph;
            obj.r_sc = sc;
            
            % Sort indices of the reduced system period vector.
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
        end
        function add_tunnel_pair(obj, ind_i, ind_j)
            % Add a tunneling pair with indices ind_i, ind_j.
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
        
        function u = get.dipole_matrix(obj)
            % Gets the dipole matrix for the mbsolve simulation
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
        
        function H = get.hamiltonian(obj)
            % Gets the Hamiltonian for the mbsolve simulation
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
        
        function [pairs_opt, num_opt] = ...
                find_optical_trans(obj, f_opt, num_opt)
            % Find optical transitions for a given frequency by calculating
            % the gain contributions of the possible transition pairs
            % for the give optical frequency.
            
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
            [~, sortedInds] = sort(gain(:), 'descend');
            top_nt = sortedInds(1:num_opt);
            [ind_top_i, ind_top_j] = ind2sub(size(gain), top_nt);
            pairs_opt = {};
            for i = 1:num_opt
                pairs_opt{end+1} = [ind_top_i(i), ind_top_j(i)];
            end
        end
        
        function [chi2, pair_DFG_ind, triplet_DFG_ind] = ...
                find_triplet_DFG(obj, num_triplets, f_IR_1, f_IR_2, ...
                f_THz, carr_dens)
            % Find level triplets for DFG generation by calculating the
            % second-order nonlinear susceptibility  for the given
            % optical frequencies.
            
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
            % Find tunneling transitions
            % by calculating the tunneling rates.
            
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
            if (nargin < 8)
                sw_warning = true;
            end
            % Find Thz transition.
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
            % Gets total scattering matrix including tunneling.
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
            % Gets total scattering matrix including optical rates with respect
            % to the given optical intensity I.
            
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
        
        function plot_gain_recovery(obj, I, t)
            % Plots gain recovery for a given optical intensity I
            % and time t.
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
            % Get occupationions wiht intensity dependency.
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
        end
        
        function freq = get_trans_freq(obj, ind_i, ind_j)
            % Get transition frequency between level i and j.
            freq = (obj.eigen_states.get_E_i(ind_i) ...
                -obj.eigen_states.get_E_i(ind_j)) ...
                * phys_const.e0 / phys_const.h;
        end
        
        function g_ij = get_gain(obj, ind_i, ind_j, f, I)
            % Get gain for the transition i to j
            % (must be element of ind_wfs).
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
        end
        
        function [f_spec, gain] = calc_gain(obj, ind, freq, I)
            % Calculate gain spectrum for the given optical transitions
            % trans_opt and frequency spectrum f_spec.
            % Find trans_opt from the set dipole pairs.
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
            % Plot gain over frequency.
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
            % Calculate group velocity dispersion (GVD) for the given
            % optical transitions trans_opt and frequency spectrum f_spec.
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
            % Plot GVD over frequency.
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
            % Get the current through the reduced system.
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
        
        function generate(obj, path)
            % Open and generate mbsolve python script.
            name_py = append(obj.project_name, '.py');
            fileID = fopen(fullfile(path, name_py), 'w');
            fprintf(fileID, ['import mbsolve.lib as mb\n', ...
                'import mbsolve.solvercpu\n', ...
                'import mbsolve.writerhdf5\n', ...
                '\nimport math\nimport time\n\n']);
            % Write Hamiltonian description.
            fprintf(fileID, '# Hamiltonian\n');
            formatSpec = '%.4e';
            energies = diag(obj.hamiltonian/phys_const.e0);
            string_e = obj.generate_string(energies, 'energies = [ ', ...
                ' * mb.E0', '%.4f');
            fprintf(fileID, append(string_e, ' ]\n'));
            ind = triu(true(size(obj.hamiltonian)), 1);
            off_diagonales = obj.hamiltonian(ind) / phys_const.e0;
            string_off = obj.generate_string(off_diagonales, ...
                'off_diagonales = [ ', ' * mb.E0', '%.4f');
            fprintf(fileID, append(string_off, ' ]\n'));
            fprintf(fileID, ...
                'H = mb.qm_operator(energies, off_diagonales)\n');
            fprintf(fileID, '\n');
            % Write dipole moment operator description.
            fprintf(fileID, '# dipole moment operator\n');
            ind = triu(true(size(obj.dipole_matrix)), 1);
            dipoles = obj.dipole_matrix(ind) / phys_const.e0;
            string_dipoles = obj.generate_string(dipoles, ...
                'off_dipoles = [ ', ' * mb.E0', '%.4e');
            fprintf(fileID, append(string_dipoles, ' ]\n'));
            string_main_dipoles = ...
                obj.generate_string(zeros( ...
                size(obj.dipole_matrix, 1), 1), ... .
            'diag_dipoles = [ ', ' * mb.E0', '%.4e');
            fprintf(fileID, append(string_main_dipoles, ' ]\n'));
            fprintf(fileID, ...
                'u = mb.qm_operator(diag_dipoles, off_dipoles)\n');
            fprintf(fileID, '\n');
            % Write relaxation superoperator description.
            fprintf(fileID, '# relaxation superoperator\n');
            fprintf(fileID, '# scattering rate matrix R\n');
            scat = obj.transition_rates;
            scat = scat - diag(diag(scat));
            for i = 1:length(scat) - 1
                if (i == 1)
                    string_scat = obj.generate_string(scat(:, i), ...
                        'rates = [ [ ', '', '%.4e');
                else
                    string_scat = obj.generate_string(scat(:, i), ...
                        '          [ ', '', '%.4e');
                end
                fprintf(fileID, append(string_scat, ' ],\n'));
            end
            string_scat = obj.generate_string(scat(:, length(scat)), ...
                '          [ ', '', '%.4e');
            fprintf(fileID, append(string_scat, ' ] ]\n'));
            fprintf(fileID, '\n');
            % Write pure dephasing description.
            fprintf(fileID, '# pure dephasing rates\n');
            % Column-major ordered pure dephasing rates vector.
            ind = triu(true(size(obj.pure_dephasing_rates)), 1);
            pure_deph = obj.pure_dephasing_rates(ind);
            string_pure_deph = obj.generate_string(pure_deph, ...
                'pure_deph = [ ', '', '%.4e');
            fprintf(fileID, append(string_pure_deph, ' ]\n'));
            fprintf(fileID, ...
                append('relax_sop', ...
                ' = mb.qm_lindblad_relaxation(rates, pure_deph)\n'));
            fprintf(fileID, '\n');
            % Initialize density matrix diagonal elements (occupations).
            fprintf(fileID, '# initial density matrix \n');
            occ = zeros(length(obj.ind_wfs));
            for i = 1:length(obj.ind_wfs)
                occ(i) = obj.carr_dist.get_occupation( ...
                    obj.ind_wfs(i));
            end
            string_occ = obj.generate_string(occ, ...
                'rho_init = mb.qm_operator([ ', '', '%.4f');
            fprintf(fileID, append(string_occ, '])\n'));
            fprintf(fileID, '\n');
            % Write quantum mechanical description.
            fprintf(fileID, '# quantum mechanical description\n');
            % Length of one period in m
            doping_dens = obj.device.dens_carrier();
            string_qm = append('qm = mb.qm_description(', ...
                num2str(doping_dens, '%.4e'), ', H, u, relax_sop)');
            fprintf(fileID, append(string_qm, '\n'));
            fprintf(fileID, ['mat_ar = mb.material("AR_', ...
                obj.project_name(1:5), ...
                '", qm,', num2str(obj.device.rel_permittivity), ',', ...
                num2str(obj.device.waveguide.overlap_factor), ',', ...
                num2str(obj.device.waveguide.a_field), ',', num2str(1), ...
                ')\n']);
            fprintf(fileID, 'mb.material.add_to_library(mat_ar)\n');
            fprintf(fileID, '\n');
            fprintf(fileID, ['dev = mb.device("', ...
                num2str(length(obj.ind_wfs)), 'lvl")\n']);
            fprintf(fileID, ['dev.add_region(mb.region("Active region"', ...
                ', mat_ar, 0.0, ', ...
                num2str(obj.device.waveguide.l_waveguide), '))\n']);
            fprintf(fileID, '\n');
            fprintf(fileID, '# Scenario\n');
            fprintf(fileID, 'ic_d = mb.ic_density_const(rho_init)\n');
            fprintf(fileID, 'ic_e = mb.ic_field_const(0.0)\n');
            fprintf(fileID, ['sce = mb.scenario("', obj.project_name, ...
                '", N_z,t_end,rho_init)\n']);
            fclose(fileID);
        end
    end
    
    methods (Static, Access = private)
        function vector_s = sort_vector(vector, E_vector)
            % Sort the given pair in the chosen order, where we bring
            % the index with the higher eigenenergy to the first position.
            [~, ind_E] = sort(E_vector, 'descend');
            vector_s = vector(ind_E);
        end
        function [p_unique, idx] = get_unique_pairs(p_cell)
            % Gets unique pairs in the given cell.
            B = cellfun(@(x) num2str(x(:)'), p_cell, ...
                'UniformOutput', false);
            [~, idx] = unique(B);
            p_unique = p_cell(idx);
        end
        function [pair_p] = check_index(pair, period)
            % Check, whether the indices of the added pair are included
            % in the period of the reduced system.
            pair_p = 0 * pair;
            n_wf = length(period);
            period_1 = mod(period-1, n_wf) + 1;
            for i = 1:length(pair)
                nx_1 = mod(pair(i)-1, n_wf) + 1;
                ind = period_1 == nx_1;
                pair_p(i) = period(ind);
            end
            if (pair ~= pair_p)
                warning(['Given states are not in de predefined ', ...
                    'period ', 'for', ...
                    ' the reduced system of the mbsolve simulation!', ...
                    ' New pair: ', num2str(pair_p), '.']);
            end
        end
        function period_corr = check_period(period)
            % Check, whether the period of the reduced system does not
            % spread over more than two spatial periods and adjust the
            % indices to those of the spatial periods 2,3, gives better
            % statistics in the EMC tool.
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
        function str_data = ...
                generate_string(data, str_init, str_add, formatSpec)
            % Generate string from input data with predefined line length.
            str_data = str_init;
            l_line = length(str_data);
            for i = 1:length(data)
                if (data(i) == 0)
                    str_i = '0';
                else
                    str_i = append(num2str(data(i), formatSpec), str_add);
                end
                if i < length(data)
                    str_i = append(str_i, ', ');
                end
                l_line = l_line + length(str_i);
                if (l_line > 78)
                    str_data = append(str_data, '\n');
                    l_line = length(str_i);
                end
                str_data = append(str_data, str_i);
            end
        end
    end
    
end
