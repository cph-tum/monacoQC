classdef scattering_rates < handle
    %scattering_rates Contains all required scattering mechanisms and the
    % correspinding rates for the reasonable description of a QCL device.
    properties (SetAccess = private)
        scattering; % Map including all scattering rates.
    end
    properties (Constant)
        dir_scat = containers.Map({'left', 'middle', 'right'}, [-1, 0, 1]);
        % Different scattering mechanisms.
    end
    
    methods
        function obj = scattering_rates(scat)
            % Constructs scattering_rates.
            obj.scattering = scat;
        end
        function R_scat = ...
                get_scattering_matrix(obj, direction, ind_wfs)
            % Returns total scattering matrix.
            direction = validatestring(direction, keys(obj.dir_scat), ...
                'get_total_scattering_matrix', 'direction');
            % Find scattering matrix for the given scattering direction.
            ind_wfs_f = ind_wfs + length(ind_wfs) ...
                * obj.dir_scat(direction);
            R_scat = values(obj.scattering); % Get scattering matrices.
            R_scat = cat(3, R_scat{:});
            % Sum over all given scattering mechanisms.
            R_scat = sum(R_scat, 3);
            n_wf = length(ind_wfs);
            % Average over the rates in the two middle spatial periods.
            R_scat((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) = ...
                (R_scat((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) + ...
                R_scat((2 * n_wf + 1):(3 * n_wf), ...
                (n_wf + 1):(4 * n_wf))) / 2;
            % Copy values to second middle spatial period.
            R_scat((2 * n_wf + 1):(3 * n_wf), (n_wf + 1):(4 * n_wf)) = ...
                R_scat((n_wf + 1):(2 * n_wf), 1:(3 * n_wf));
            % Set intraband scattering rates to zero.
            R_scat = R_scat - diag(diag(R_scat));
            % Calculate outscattering rates of levels.
            R_scat = R_scat - diag(sum(R_scat, 2));
            % Return scattering matrix only for the relevant levels.
            R_scat = R_scat(ind_wfs, ind_wfs_f);
        end
        function R_id = ...
                get_scattering_matrix_id(obj, id, direction, ind_wfs)
            % Returns scattering matrix of the specified mechanism.
            % Find the scattering matrix for the given scattering
            % mechanism.
            id = validatestring(id, keys(obj.scattering), ...
                'get_scattering_matrix', 'id scattering mechanism');
            % Find scattering matrix for the given scattering direction.
            direction = validatestring(direction, keys(obj.dir_scat), ...
                'get_total_scattering_matrix', 'direction');
            ind_wfs_f = ind_wfs + length(ind_wfs) ...
                * obj.dir_scat(direction);
            R_id = obj.scattering(id);
            % Average over the rates in the two middle spatial periods.
            n_wf = length(ind_wfs);
            R_id((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) = ...
                (R_id((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) + ...
                R_id((2 * n_wf + 1):(3 * n_wf), ...
                (n_wf + 1):(4 * n_wf))) / 2;
            % Copy values to second middle spatial period.
            R_id((2 * n_wf + 1):(3 * n_wf), (n_wf + 1):(4 * n_wf)) = ...
                R_id((n_wf + 1):(2 * n_wf), 1:(3 * n_wf));
            % Return scattering matrix only for the relevant levels.
            R_id = R_id(ind_wfs, ind_wfs_f);
        end
        function gamma_scat = ...
                get_level_broadening(obj, ind)
            % Get level broadening for levels index ind.
            R_scat = values(obj.scattering); % Get scattering matrices.
            R_scat = cat(3, R_scat{:});
            % Sum over all given scattering mechanisms.
            R_scat = sum(R_scat, 3);
            n_wf = size(R_scat, 2) / 4;
            % Average over the rates in the two middle spatial periods.
            R_scat((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) = ...
                (R_scat((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) + ...
                R_scat((2 * n_wf + 1):(3 * n_wf), ...
                (n_wf + 1):(4 * n_wf))) / 2;
            % Copy values to second middle spatial period.
            R_scat((2 * n_wf + 1):(3 * n_wf), (n_wf + 1):(4 * n_wf)) = ...
                R_scat((n_wf + 1):(2 * n_wf), 1:(3 * n_wf));
            % Set intraband scattering rates to zero.
            R_scat = R_scat - diag(diag(R_scat));
            % Calculate outscattering rates of levels.
            gamma_scat = sum(R_scat(ind, :), 2) / 2 / pi;
        end
        
        function gamma_scat_id = ...
                get_level_broadening_id(obj, id, ind)
            % Returns scattering matrix of the specified mechanism.
            % Find the scattering matrix for the given scattering
            % mechanism.
            id = validatestring(id, keys(obj.scattering), ...
                'get_scattering_matrix', 'id scattering mechanism');
            R_id = obj.scattering(id);
            % Average over the rates in the two middle spatial periods.
            n_wf = size(R_id, 2) / 4;
            R_id((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) = ...
                (R_id((n_wf + 1):(2 * n_wf), 1:(3 * n_wf)) + ...
                R_id((2 * n_wf + 1):(3 * n_wf), ...
                (n_wf + 1):(4 * n_wf))) / 2;
            % Copy values to second middle spatial period.
            R_id((2 * n_wf + 1):(3 * n_wf), (n_wf + 1):(4 * n_wf)) = ...
                R_id((n_wf + 1):(2 * n_wf), 1:(3 * n_wf));
            % Set intraband scattering rates to zero.
            R_id = R_id - diag(diag(R_id));
            gamma_scat_id = sum(R_id(ind, :), 2) / 2 / pi;
        end
        
        function gamma_lt = ...
                get_lifetime_broadening(obj, ind_i, ind_j)
            gamma_lt = (obj.get_level_broadening(ind_i) + ...
                obj.get_level_broadening(ind_j)) / 2;
        end
        function gamma_lt_id = ...
                get_lifetime_broadening_id(obj, id, ind_i, ind_j)
            gamma_lt_id = (obj.get_level_broadening_id(id, ind_i) + ...
                obj.get_level_broadening_id(id, ind_j)) / 2;
        end
    end
    
end
