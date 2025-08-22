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

classdef conduction_band < handle
    % Describes the conduction band profile including space charge effects.

    properties (SetAccess = private)
        zv % vector: Position vector [Angstrom].
        Vh % vector: Conduction band potential profile [J].
    end

    methods
        function obj = conduction_band(zv, Vh)
            % Constructs an object of type conduction_band.
            %
            % Syntax:
            %   obj = conduction_band(zv, Vh)
            %
            % Input Arguments:
            %   zv (vector): Position vector of the potential [Angstrom].
            %   Vh (vector): Conduction band potential profile [J].

            obj.zv = zv;
            obj.Vh = Vh;
        end

        function plot_profile(obj)
            % Plots conduction band profile with space charge effects
            % and probability densities.
            %
            % Syntax:
            %   plot_profile(obj)

            plot(-(obj.zv / 10), obj.Vh/phys_const.e0, 'Color', ...
                [0, 0, 0], 'linestyle', '--', 'linewidth', 2);
            xlabel('z/nm'); %'Position x in nm'); %
            ylabel('E/eV'); %'Energy in eV'); %
        end

        function write_cond_profile(obj, filename)
            % Write conduction band potential into filename.csv file.
            %
            % Syntax:
            %   write_cond_profile(obj, filename)
            %
            % Input Arguments:
            %   filenanme (char): Specifies absolute or relativ path of
            %     csv-file (without extension).

            T_cond_profile = array2table([(-1) * obj.zv' / 10, obj.Vh' ./ phys_const.e0]);
            header = {'z', 'V'};
            T_cond_profile.Properties.VariableNames(1:size(header, 2)) ...
                = header;
            writetable(T_cond_profile, filename);
        end

        function to_hdf5(obj, filename)
            % Saves conduction_band object in hdf5 format.
            %
            % Syntax:
            %   to_hdf5(obj, filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            z_pos = -obj.zv(end:-1:1);
            cb_pot = obj.Vh(end:-1:1);
            [uniq_val, uniq_idx] = unique(z_pos, 'stable');
            dupe_idx = setdiff(1:numel(z_pos), uniq_idx);

            if_pos = z_pos(dupe_idx);
            if_off = cb_pot(dupe_idx) - cb_pot(dupe_idx-1);

            write_to_hdf5(filename, "/interface_position", if_pos);
            write_to_hdf5(filename, "/interface_cb_offset", if_off)
            write_to_hdf5(filename, "/cb_potential", cb_pot(uniq_idx));
        end
    end

    methods (Static)
        function obj = from_hdf5(filename)
            % Constructs conduction_band object from hdf5 file.
            %
            % Syntax:
            %   obj = from_hdf5(filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            z_wf = h5read(filename, "/z_position");
            if_pos = h5read(filename, "/interface_position");
            if_off = h5read(filename, "/interface_cb_offset");
            cb = h5read(filename, "/cb_potential");
            % check if dimensions match
            if length(z_wf) == length(cb)
                z = z_wf;
            else
                z = linspace(0, z_wf(end), length(cb))';
            end
            % indices of duplicate z values
            if_idx = find(ismember(z, if_pos));
            z_concat = cat(1, z, if_pos);
            V_concat = cat(1, cb, cb(if_idx)+if_off);
            [z_cb, sort_idx] = sort(z_concat);
            V_cb = V_concat(sort_idx);
            % create object
            obj = conduction_band(-z_cb(end:-1:1)', V_cb(end:-1:1)');
        end
    end
end
