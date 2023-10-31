classdef dephasing_rates < handle
    % dephasing_rates Contains the dephasing rates, consisting of the
    % lifetime broadening and pure dephasing rates.
    properties (SetAccess = private)
        pure_dephasing % Map including pure dephasing rates divided
        pure_dephasing_k % into different mechanisms (k-resolved).
        lvl_broadening % Matrix including the level broadening
        lvl_broadening_k % (k-resolved).
    end
    methods
        function obj = ...
                dephasing_rates(pure_deph_mech, ...
                lvl_broad, carr_dist, key_deph)
            % Constructs dephasing_rates.
            
            % Level broadening k-resolved.
            obj.lvl_broadening_k = lvl_broad; % Level broadening.
            obj.set_lvl_broadening_avr(carr_dist);
            % Pure dephasing rates divided into different mechanisms
            % k-resolved.
            obj.pure_dephasing_k = pure_deph_mech;
            obj.set_pure_dephasing_rate_avr(carr_dist, key_deph);
        end
        % Return the lifetime broadening between state i and j in k-space.
        function lt_b = get_lt_broadening(obj, ind_i, ind_j, ind_k)
            % Returns lifetime broadening for the given pair of
            % states (k-resolved).
            if nargin > 3
                lt_b = (obj.lvl_broadening_k(ind_i, ind_k) ...
                    +obj.lvl_broadening_k(ind_j, ind_k)) / 2;
            else
                lt_b = (obj.lvl_broadening(ind_i) ...
                    +obj.lvl_broadening(ind_j)) / 2;
            end
        end
        % Return dephasing rate between state i and j (k-resolved).
        function gamma_ijk = ...
                get_dephasing_rate(obj, ind_i, ind_j, ind_k)
            % Return dephasing rate dependent on the inputarguments.
            % Default pure dephasinge averaged + livetime broadening
            % averaged.
            puredephasing = ...
                obj.pure_dephasing('total');
            puredephasing = puredephasing(ind_i, ind_j);
            lt_broad = obj.get_lt_broadening(ind_i, ind_j);
            % Return default k-resolved.
            if nargin > 3
                puredephasing = ...
                    obj.pure_dephasing_k('total');
                puredephasing = puredephasing(ind_i, ind_j, ind_k);
                lt_broad = obj.get_lt_broadening(ind_i, ind_j, ind_k);
            end
            % Returns dephasing rate for the given pair of states
            % (k-resolved).
            gamma_ijk = lt_broad + puredephasing;
        end
        % Return the total pure dephasing rate between state i and j
        % (k-resolved).
        function gamma_p_ijk = ...
                get_pure_dephasing_rate(obj, ind_i, ind_j, id, ind_k)
            % Returns pure dephasing rate for the given pair of states
            % k-resolved.
            gpm = obj.pure_dephasing('total');
            gamma_p_ijk = gpm(ind_i, ind_j);
            if nargin > 3 && ischar(id)
                gpm = obj.pure_dephasing(id);
                gamma_p_ijk = gpm(ind_i, ind_j);
                if nargin > 4
                    gpm = ...
                        obj.pure_dephasing_k(id);
                    gamma_p_ijk = gpm(ind_i, ind_j, ind_k);
                end
            end
            % Returns pure dephasing rate of mechanisms(char)
            % for the given pair of states k-resolved.
            % char: {'interface roughness', 'impurity',
            %           'optical absorption', 'optical emission',
            %           'acoustic phonon', 'total'}
        end
        %%%%%%%%%%%%%%%% Dephasing rates averaged over FE %%%%%%%%%%%%%%%%%
        % Return level broadening of state i averaged over the carrier
        % distribution.
        function set_lvl_broadening_avr(obj, c_dist)
            dist_carr = transpose(c_dist.distribution);
            num_wfs = length(obj.lvl_broadening_k(:, 1)) / 4;
            % Returns level broadening averaged over the carrier
            % distribution.
            for ind_i = 1:length(obj.lvl_broadening_k(:, 1))
                ni = mod(ind_i-1, num_wfs) + 1;
                obj.lvl_broadening(ind_i, 1) = trapz(c_dist.E_kin, ...
                    (obj.lvl_broadening_k(ind_i, :)).*dist_carr(ni, :)) / ...
                    trapz(c_dist.E_kin, dist_carr(ni, :));
            end
        end
        % Set the pure dephasing rate between averaged over the
        % carrier distribution.
        function set_pure_dephasing_rate_avr(obj, c_dist, key)
            % Averageing factor.
            dist_carr = transpose(c_dist.distribution);
            num_wfs = length(obj.lvl_broadening_k(:, 1)) / 4;
            obj.pure_dephasing = containers.Map();
            for k = 1:length(key)
                for ind_i = 1:length(obj.lvl_broadening_k(:, 1))
                    for ind_j = 1:length(obj.lvl_broadening_k(:, 1))
                        ni = mod(ind_i-1, num_wfs) + 1;
                        nj = mod(ind_j-1, num_wfs) + 1;
                        % Prevent dividing through zero by taking the
                        % distribution of state i.
                        if (ni == nj)
                            dif_dis = dist_carr(ni, :);
                        else
                            dif_dis = abs(dist_carr(ni, :) ...
                                -dist_carr(nj, :));
                        end
                        gp = obj.pure_dephasing_k(key{k});
                        gp = reshape(gp, ...
                            length(obj.lvl_broadening_k(:, 1)), ...
                            length(obj.lvl_broadening_k(:, 1)), ...
                            length(dist_carr(1, :)));
                        % Reshape to match dif_dis for multiplication.
                        gamh = reshape(gp(ind_i, ind_j, :), ...
                            1, length(gp(1, 1, :)));
                        % Returns averaged pure dephasing rate for the
                        % given pair
                        % of states.
                        pure_deph(ind_i, ind_j) = ...
                            sum((gamh .* dif_dis)) ./ sum(dif_dis);
                    end
                end
                obj.pure_dephasing(key{k}) = pure_deph;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%% Post processing %%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_dephasing_k(obj, ind_i, ind_j, id, c_dist)
            switch id
                case 'total'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_dephasing_rate(ind_i, ind_j, i);
                    end
                    g_ij = obj.get_dephasing_rate(ind_i, ind_j);
                case 'ltbroad'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = obj.get_lt_broadening(ind_i, ind_j, i);
                    end
                    g_ij = obj.get_lt_broadening(ind_i, ind_j);
                otherwise
                    error(['Error: Choose one of the following', ...
                        ' options {''total'', ''ltbroad''}.'])
            end
            hold on;
            set(gca, 'FontSize', 12);
            plot(c_dist.E_kin, deph/1e12/2/pi);
            grid on;
            title(id, 'FontSize', 17)
            x = xlabel('E_{kin} [eV]');
            x.FontSize = 17;
            x.Position(2) = x.Position(2) - 0.001;
            y = ylabel('\gamma_{ij} [ps^{-1}]');
            y.FontSize = 17;
            y.Position(1) = y.Position(1) - 0.001;
            hold off;
            disp(['Average dephasing rate: ', ...
                num2str(g_ij/1e12), 'e+12'])
        end
        
        function plot_pure_dephasing_k(obj, ind_i, ind_j, id, c_dist)
            switch id
                case 'puredeph'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_pure_dephasing_rate ...
                            (ind_i, ind_j, 'total', i);
                    end
                    g_ij = obj.get_pure_dephasing_rate ...
                        (ind_i, ind_j, 'total');
                case 'interface roughness'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_pure_dephasing_rate ...
                            (ind_i, ind_j, id, i);
                    end
                    g_ij = obj.get_pure_dephasing_rate ...
                        (ind_i, ind_j, id);
                case 'impurity'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_pure_dephasing_rate ...
                            (ind_i, ind_j, id, i);
                    end
                    g_ij = obj.get_pure_dephasing_rate ...
                        (ind_i, ind_j, id);
                case 'optical absorption'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_pure_dephasing_rate ...
                            (ind_i, ind_j, id, i);
                    end
                    g_ij = obj.get_pure_dephasing_rate ...
                        (ind_i, ind_j, id);
                case 'optical emission'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_pure_dephasing_rate ...
                            (ind_i, ind_j, id, i);
                    end
                    g_ij = obj.get_pure_dephasing_rate ...
                        (ind_i, ind_j, id);
                case 'acoustic phonon'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_pure_dephasing_rate ...
                            (ind_i, ind_j, id, i);
                    end
                    g_ij = obj.get_pure_dephasing_rate ...
                        (ind_i, ind_j, id);
                otherwise
                    error(['Error: Choose one of the following ', ...
                        ['options {''puredeph'', ', ...
                        '''interface roughness''', ...
                        ', ''impurity'', ''optical absorption'''], ...
                        ',', ' ''optical emission''', ...
                        ', ''acoustic phonon''}.'])
            end
            hold on;
            set(gca, 'FontSize', 12);
            plot(c_dist.E_kin, deph/1e12/2/pi);
            grid on;
            title(id, 'FontSize', 17)
            x = xlabel('E_{kin} [eV]');
            x.FontSize = 17;
            x.Position(2) = x.Position(2) - 0.001;
            y = ylabel('\gamma_{ij} [ps^{-1}]');
            y.FontSize = 17;
            y.Position(1) = y.Position(1) - 0.001;
            hold off;
            disp(['Average dephasing rate: ', ...
                num2str(g_ij/1e12), 'e+12'])
        end
    end
end
