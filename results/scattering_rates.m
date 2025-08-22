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

classdef scattering_rates < handle
    % Contains scattering rates for all considered scattering
    % mechanisms for the QCL.

    properties
        scattering % containers.Map: Scattering rates divided into different scattering mechanisms.
    end

    properties (Constant)
        % containers.Map: possible scattering directions.
        dir_scat = containers.Map({'left', 'middle', 'right'}, [-1, 0, 1]);
        % cell-array: Name of scattering rate mechanisms.
        keys_scat = { ...
            'electron-electron', 'optical absorption', 'optical emission', ...
            'tunneling', 'alloy disorder', 'acoustic phonon', 'impurity', ... .
        'interface roughness'};
        hdf5_ds_names = {; ...
            'electron_electron', 'lo_phonon_absorption', 'lo_phonon_emission', ...
            'tunneling', 'alloy', 'acoustic_phonon', 'impurity', ...
            'interface_roughness'};
    end

    methods
        function obj = scattering_rates(scat)
            % Constructs an object of type scattering_rates.
            %
            % Syntax:
            %   obj = scattering_rates(scat)
            %
            % Input Arguments:
            %   scat (containers.Map): Scattering rates for 4 QCL periods
            %     divided into the different scattering mechanisms.

            obj.scattering = scat;
        end

        function R_scat = ...
                get_scattering_matrix(obj, direction, ind_wfs)
            % Returns total scattering rate matrix between the specified
            % states.
            %
            % Syntax:
            %   R_scat = get_scattering_matrix(obj, direction, ind_wfs)
            %
            % Input Arguments:
            %   direction (char): Specifies the scattering direction. Valid
            %     choices are ``middle`` (scattering within a period),
            %     ``left`` (scattering to the left period) and ``right``
            %     (scattering to the right period).
            %   ind_wfs (vector): Indices of states for which the
            %      scattering rates should be returned.
            %
            % Output Argument:
            %   R_scat (matrix): Total scattering rate matrix.

            direction = validatestring(direction, keys(obj.dir_scat), ...
                'get_total_scattering_matrix', 'direction');
            % Find scattering matrix for the given scattering direction.
            ind_wfs_f = ind_wfs + length(ind_wfs) ...
                * obj.dir_scat(direction);
            R_scat = values(obj.scattering); % Get scattering matrices.
            R_scat = cat(3, R_scat{:});
            % Sum over all given scattering mechanisms.
            R_scat = sum(R_scat, 3);
            n_wf = size(R_scat, 1) / 4;
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
            % Returns scattering rate matrix between specified states
            % for a specific scattering mechanism.
            %
            % Syntax:
            %   R_scat = get_scattering_matrix(obj, direction, ind_wfs)
            %
            % Input Arguments:
            %   id (char): Name of scattering rate mechanism. Valid names
            %     are ``electron-electron``, ``optical absorption``,
            %     ``optical emission``, ``tunneling``, ``alloy disorder``,
            %     ``acoustic phonon``, ``impurity`` and
            %     ``interface roughness``.
            %   direction (char): Specifies the scattering direction. Valid
            %     choices are ``middle`` (scattering within a period),
            %     ``left`` (scattering to the left period) and ``right``
            %     (scattering to the right period).
            %   ind_wfs (vector): Indices of states for which the
            %      scattering rates should be returned.
            %
            % Output Arguments:
            %   R_id (matrix): Scattering rate matrix.

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
            n_wf = size(R_id, 1) / 4;
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

        function gamma_scat = get_level_broadening(obj, ind)
            % Returns level broadening of specified states.
            %
            % Syntax:
            %   gamma_scat = get_level_broadening(obj, ind)
            %
            % Input Arguments:
            %   ind (scalar | vector): Indices of the states.
            %
            % Output Arguments:
            %   gamma_scat (scalar | vector): Inverse lifetime divided by
            %     2*pi.

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

        function gamma_scat_id = get_level_broadening_id(obj, id, ind)
            % Returns level broadening of specified states for a
            % specific scattering mechanism.
            %
            % Syntax:
            %   gamma_scat_id = get_level_broadening_id(obj, id, ind)
            %
            % Input Arguments:
            %   id (char): Name of scattering mechanism. Valid names are
            %     ``electron-electron``, ``optical absorption``,
            %     ``optical emission``, ``tunneling``, ``alloy disorder``,
            %     ``acoustic phonon``, ``impurity`` and
            %     ``interface roughness``.
            %   ind (scalar | vector): Indices of the states.
            %
            % Output Arguments:
            %   gamma_scat_id (scalar | vector): Inverse lifetime divided
            %     by 2*pi.

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

        function gamma_lt = get_lifetime_broadening(obj, ind_i, ind_j)
            % Returns lifetime broadening value for a specific transition
            % between two states.
            %
            % Syntax:
            %   gamma_lt = get_lifetime_broadening(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_i (scalar): Index of initial state.
            %   ind_j (scalar): Index of final state.
            %
            % Output arguments:
            %   gamma_lt (scalar): Lifetime broadening value.

            gamma_lt = (obj.get_level_broadening(ind_i) + ...
                obj.get_level_broadening(ind_j)) / 2;
        end

        function gamma_lt_id = ...
                get_lifetime_broadening_id(obj, id, ind_i, ind_j)
            % Returns lifetime broadening value for a specific transition
            % between two states for a specific scattering mechanism.
            %
            % Syntax:
            %   gamma_lt_id = get_lifetime_broadening_id(obj, id, ind_i, ind_j)
            %
            % Input Arguments:
            %   id (char): Name of scattering mechanism. Valid names are
            %     ``electron-electron``, ``optical absorption``,
            %     ``optical emission``, ``tunneling``, ``alloy disorder``,
            %     ``acoustic phonon``, ``impurity`` and
            %     ``interface roughness``.
            %   ind_i (scalar): Index of initial state.
            %   ind_j (scalar): Index of final state.
            %
            % Output arguments:
            %   gamma_lt_id (scalar): Lifetime broadening value.

            gamma_lt_id = (obj.get_level_broadening_id(id, ind_i) + ...
                obj.get_level_broadening_id(id, ind_j)) / 2;
        end

        function names = get_scattering_names(obj)
            % Returns the names of the scattering mechanisms.
            %
            % Syntax:
            %   names = get_scattering_names(obj)
            %
            % Output Arguments:
            %   names (cell-array): Array containing the names of the
            %     scattering mechanisms.

            names = keys(obj.scattering);
        end

        function add_scattering_rate(obj, name, scattering_matrix)
            % Adds the scattering rates of an additional scattering
            % mechanism. If the mechanism already exists, the new values
            % replace the old ones.
            %
            % Syntax:
            %   add_scattering_rate(obj, name, scattering_matrix)
            %
            % Input Arguments:
            %   name (char): Name of the new scattering mechanism.
            %   scattering_matrix (matrix): Scattering rates for the 4 QCL
            %      periods.

            obj.scattering(name) = scattering_matrix;
        end

        function to_hdf5(obj, filename)
            % Saves scattering_rates object in hdf5 format.
            %
            % Syntax:
            %   to_hdf5(obj, filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            for i = 1:length(obj.keys_scat)
                write_to_hdf5(filename, ...
                    "/scattering_matrix/"+string(obj.hdf5_ds_names{i}), ...
                    obj.scattering(obj.keys_scat{i})');
            end
        end
    end

    methods (Static)
        function obj = from_hdf5(filename)
            % Constructs scattering_rates object from hdf5 file.
            %
            % Syntax:
            %   obj = from_hdf5(filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            scattering = containers.Map();
            for i = 1:length(scattering_rates.keys_scat)
                mechanism = string(scattering_rates.hdf5_ds_names{i});
                rates = h5read(filename, "/scattering_matrix/"+mechanism);
                scattering(scattering_rates.keys_scat{i}) = rates';
            end
            obj = scattering_rates(scattering);
        end
    end
end
