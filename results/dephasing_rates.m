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

classdef dephasing_rates < handle
    % Contains the dephasing rates including lifetime broadening and the
    % pure dephasing rates.

    properties (SetAccess = private)
        pure_dephasing % containers.Map: Pure dephasing rates devided into different scattering mechanisms.
        pure_dephasing_k % containers.Map: K-dependent pure dephasing rates devided into different scattering mechanisms
        lvl_broadening % vector: Level broadening values for each subband.
        lvl_broadening_k % matrix: Level broadening k-resolved for each subband.
    end

    properties (Constant)
        keys_deph = {'interface roughness', 'impurity', ...
            'optical absorption', 'optical emission', 'acoustic phonon', ...
            'total'};
        hdf5_ds_names = {'interface_roughness', 'impurity', ...
            'lo_phonon_absorption', 'lo_phonon_emission', 'acoustic_phonon', ...
            'total'};
    end

    methods
        function obj = ...
                dephasing_rates(pure_deph_mech, ...
                lvl_broad, carr_dist, key_deph)
            % Constructs an object of type dephasing_rates.
            %
            % Syntax:
            %   obj = dephasing_rates(pure_deph_mech, lvl_broad, carr_dist, key_deph)
            %
            % Input Arguments:
            %   pure_deph_mech (containers.Map): Contains the k-resolved
            %     pure dephasing rates (3-d array) for all four QCL-periods
            %     devided into the different scattering mechanisms.
            %   lvl_broad (matrix): K-resolved inverse lifetimes for all
            %     subbands.
            %   carr_dist (carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.
            %   key_deph (cell-array): Array containing the names of the
            %     considered dephasing mechanisms, which correspond to the
            %     keys of the containers.Map variable pure_deph_mech.

            % Level broadening k-resolved.
            obj.lvl_broadening_k = lvl_broad; % Level broadening.
            obj.set_lvl_broadening_avr(carr_dist);
            % Pure dephasing rates divided into different mechanisms
            % k-resolved.
            obj.pure_dephasing_k = pure_deph_mech;
            obj.set_pure_dephasing_rate_avr(carr_dist, key_deph);
        end

        function lt_b = get_lt_broadening(obj, ind_i, ind_j, ind_k)
            % Returns the lifetime broadening for a given pair of states
            % either k-dependent or k-averaged.
            %
            % Syntax:
            %   lt_b = get_lt_broadening(obj, ind_i, ind_j)
            %   lt_b = get_lt_broadening(obj, ind_i, ind_j, ind_k)
            %
            % Input Arguments:
            %   ind_i (scalar | vector): Initial subband index.
            %   ind_j (scalar | vector): Final subband index.
            %   ind_k (scalar | vector): Index of k-vector. If ind_k is not
            %     provided as input, the k-averaged value(s) are returned.
            %
            % Output Arguments:
            %   lt_b (scalar | matrix): Lifetime broadening in 1/s.

            if length(ind_i) > 1 || length(ind_j) > 1
                if nargin > 3
                    lt_b = (reshape(obj.lvl_broadening_k(ind_i, ind_k), ...
                        length(ind_i), 1, length(ind_k)) ...
                        +reshape(obj.lvl_broadening_k(ind_j, ind_k), ...
                        1, length(ind_j), length(ind_k))) / 2;
                else
                    lt_b = (reshape(obj.lvl_broadening(ind_i), [], 1) ...
                        +reshape(obj.lvl_broadening(ind_j), 1, [])) / 2;
                end
            else
                if nargin > 3
                    lt_b = (obj.lvl_broadening_k(ind_i, ind_k) ...
                        +obj.lvl_broadening_k(ind_j, ind_k)) / 2;
                else
                    lt_b = (obj.lvl_broadening(ind_i) ...
                        +obj.lvl_broadening(ind_j)) / 2;
                end
            end
        end

        function gamma_ijk = ...
                get_dephasing_rate(obj, ind_i, ind_j, ind_k)
            % Returns total dephasing rates, i.e. pure dephasing rates +
            % lifetime broadenings.
            %
            % Syntax:
            %   lt_b = get_lt_broadening(obj, ind_i, ind_j)
            %   lt_b = get_lt_broadening(obj, ind_i, ind_j, ind_k)
            %
            % Input Arguments:
            %   ind_i (scalar | vector): Initial subband index.
            %   ind_j (scalar | vector): Final subband index.
            %   ind_k (scalar | vector): Index of k-vector. If ind_k is not
            %     provided as input, the k-averaged value(s) are returned.
            %
            % Output Arguments:
            %   lt_b (scalar | matrix): Dephasing rate(s) in 1/s.

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

        function obj = set_lvl_broadening(obj, lvl_broad)
            % Set lvl broadening.
            obj.lvl_broadening = lvl_broad;
        end

        function gamma_p_ijk = ...
                get_pure_dephasing_rate(obj, ind_i, ind_j, id, ind_k)
            % Returns pure dephasing rate for the given pair of states
            % k-resolved.
            %
            % Syntax:
            %   gamma_p_ijk = get_pure_dephasing_rate(obj, ind_i, ind_j)
            %   gamma_p_ijk = get_pure_dephasing_rate(obj, ind_i, ind_j, id)
            %   gamma_p_ijk = get_pure_dephasing_rate(obj, ind_i, ind_j, id, ind_k)
            %
            % Input Arguments:
            %   ind_i (scalar | vector): Initial subband index.
            %   ind_j (scalar | vector): Final subband index.
            %   id ('total' (default) | char): Mechanism for which the pure
            %     dephasing rates should be returned. Valid names are
            %     ``interface roughness``, ``impurity``,
            %     ``optical absorption``, ``optical emission``,
            %     ``acoustic phonon`` and ``total``.
            %   ind_k (scalar | vector): Index of k-vector. If ind_k is not
            %     provided as input, the k-averaged value(s) are returned.

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

        function set_lvl_broadening_avr(obj, c_dist)
            % Sets lvl_broadening property by averaging k-dependent level
            % broadening values for each subband over the corresponding
            % carrier distribution.
            %
            % Syntax:
            %   set_lvl_broadening_avr(obj, c_dist)
            %
            % Input Arguments:
            %   c_dist (carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.

            dist_carr = transpose(c_dist.distribution);
            num_wfs = length(obj.lvl_broadening_k(:, 1)) / 4;
            if size(c_dist.E_kin, 2) == 1
                Ekin = repmat(c_dist.E_kin, 1, num_wfs);
            else
                Ekin = c_dist.E_kin;
            end
            % Returns level broadening averaged over the carrier
            % distribution.
            for ind_i = 1:length(obj.lvl_broadening_k(:, 1))
                ni = mod(ind_i-1, num_wfs) + 1;
                obj.lvl_broadening(ind_i, 1) = trapz(Ekin(:, ni), ...
                    (obj.lvl_broadening_k(ind_i, :)).*dist_carr(ni, :)) / ...
                    trapz(Ekin(:, ni), dist_carr(ni, :));
            end
        end

        function set_pure_dephasing_rate_avr(obj, c_dist, key)
            % Sets pure_dephasing property by averaging the k-dependent
            % pure dephasing rates between the subbands over the carrier
            % distributions.
            %
            % Syntax:
            %   set_pure_dephasing_rate_avr(obj, c_dist, key)
            %
            % Input Arguments:
            %   c_dist (carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.
            %   key (cell-array): Names of the scattering mechanisms for
            %     which the averaged pure dephasing rates should be
            %     calculated.

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

        function plot_dephasing_k(obj, ind_i, ind_j, id, c_dist)
            % Plots the k-dependent dephasing rate over kinetik energy for
            % a specific pair of subbands by either including or excluding
            % the the contribution of the pure dephasing rates.
            %
            % Syntax:
            %   plot_dephasing_k(obj, ind_i, ind_j, id, c_dist)
            %
            % Input Arguments:
            %   ind_i (scalar): Initial subband index.
            %   ind_j (scalar): Final subband index.
            %   id (char): Flag for including (``id = total``) or
            %     excluding (``id = ltbroad``) the pure dephasing rates.
            %   c_dist(carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.

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
            num_states = length(c_dist.occupation);
            plot(c_dist.E_kin(:, mod(ind_i-1, num_states)+1), deph/1e12/2/pi);
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
            % Plots the k-dependent pure dephasing rate over kinetik energy
            % for a specific pair of subbands and specific scattering
            % mechanism.
            %
            % Syntax:
            %   plot_dephasing_k(obj, ind_i, ind_j, id, c_dist)
            %
            % Input Arguments:
            %   ind_i (scalar): Initial subband index.
            %   ind_j (scalar): Final subband index.
            %   id (char): Mechanism for which the pure
            %     dephasing rates should be plottet. Valid names are
            %     ``interface roughness``, ``impurity``,
            %     ``lo phonon absorption``, ``lo phonon emission``,
            %     ``acoustic phonon`` and ``puredeph`` (= total).
            %   c_dist(carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.

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
                case 'lo phonon absorption'
                    for i = 1:length(obj.lvl_broadening_k(1, :))
                        deph(1, i) = ...
                            obj.get_pure_dephasing_rate ...
                            (ind_i, ind_j, id, i);
                    end
                    g_ij = obj.get_pure_dephasing_rate ...
                        (ind_i, ind_j, id);
                case 'lo phonon emission'
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
                        ', ''impurity'', ''lo phonon absorption'''], ...
                        ',', ' ''lo phonon emission''', ...
                        ', ''acoustic phonon''}.'])
            end
            hold on;
            set(gca, 'FontSize', 12);
            num_states = length(c_dist.occupation);
            plot(c_dist.E_kin(:, mod(ind_i-1, num_states)+1), deph/1e12/2/pi);
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

        function to_hdf5(obj, filename)
            % Saves dephasing_rates object in hdf5 format.
            %
            % Syntax:
            %   to_hdf5(obj, filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            write_to_hdf5(filename, "/level_broadening", obj.lvl_broadening_k');
            for i = 1:length(obj.keys_deph)
                puredep = obj.pure_dephasing_k(obj.keys_deph{i});
                write_to_hdf5(filename, ...
                    "/pure_dephasing/"+string(obj.hdf5_ds_names{i}), ...
                    permute(puredep, [3, 2, 1]));
            end
        end
    end

    methods (Static)
        function obj = from_hdf5(filename, carr_dist)
            % Constructs dephasing_rates object from hdf5 file.
            %
            % Syntax:
            %   obj = from_hdf5(filename)
            %   obj = from_hdf5(filename, carr_dist)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.
            %   carr_dist (carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.

            if nargin < 2
                carr_dist = carrier_distribution.from_hdf5(filename);
            end

            lvl_broadening_k = h5read(filename, "/level_broadening");
            pure_dephasing_mech = containers.Map();
            for i = 1:length(dephasing_rates.keys_deph)
                puredep = h5read(filename, ...
                    "/pure_dephasing/"+string(dephasing_rates.hdf5_ds_names{i}));
                pure_dephasing_mech(dephasing_rates.keys_deph{i}) = permute(puredep, [3, 2, 1]);
            end
            obj = dephasing_rates(pure_dephasing_mech, lvl_broadening_k', ...
                carr_dist, dephasing_rates.keys_deph);
        end
    end
end
