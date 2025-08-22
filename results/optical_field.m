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

classdef optical_field < handle
    % Contains the characterisitcs of the optical cavity field.

    properties
        intensity % vector: Intensities of the cavity modes inside the cavity [W/m^2].
        freq % vector: Frequencies of the cavity modes [Hz].
        gamma = 1 % scalar: Overlap factor [-].
        R = 0 % scalar: Power reflectance [-].
        A = 0 % scalar: Cross section area [m^2].
    end

    methods
        function obj = optical_field(intensity, freq)
            % Constructs an object of type optical_field.
            %
            % Syntax:
            %   obj = optical_field(intensity, freq)
            %
            % Input Arguments:
            %   intensity (vector): Intensities of the cavity modes inside
            %    the cavity [W/m^2].
            %   freq (vector): Frequencies of the cavity modes [Hz].

            obj.intensity = reshape(intensity, length(intensity), 1);
            obj.freq = reshape(freq, size(obj.intensity));
        end

        function set_waveguide_properties(obj, device)
            % Set waveguide properties from device object.
            %
            % Syntax:
            %   set_waveguide_properties(obj, device)
            %
            % Input Arguments:
            %   device (device-object): Contains information about the
            %     structure, geometry and materials of the QCL-device.

            obj.gamma = device.waveguide.overlap_factor;
            obj.R = device.waveguide.refl_left * ...
                device.waveguide.refl_right;
            obj.A = device.waveguide.A_act;
        end

        function power = get_power_spec(obj, id)
            % Returns the optical power for each cavity mode.
            %
            % Syntax:
            %   power = get_power_spec(obj, id)
            %
            % Input Arguments:
            %   id (char): Flag specifying if the power spectrum inside
            %     or outside the cavity should be returned. Valid choices
            %     for id are ``cavity`` and ``outcoupled``.
            %
            % Output arguments:
            %   power (vector): Optical power spectrum.

            arguments
                obj
                id{mustBeMember(id, {'cavity', 'outcoupled'})}
            end
            I = obj.get_intensity_spec(id);
            power = I .* obj.A;
        end

        function intensity = get_intensity_spec(obj, id)
            % Returns the optical intensity for each cavity mode.
            %
            % Syntax:
            %   intensity = get_intensity_spec(obj, id)
            %
            % Input Arguments:
            %   id (char): Flag specifying if the intensity spectrum inside
            %     or outside the cavity should be returned. Valid choices
            %     for id are ``cavity`` and ``outcoupled``.
            %
            % Output arguments:
            %   intensity (vector): Optical intensity spectrum.

            arguments
                obj
                id{mustBeMember(id, {'cavity', 'outcoupled'})}
            end
            I = obj.intensity ./ reshape(obj.gamma, [], 1);
            if id == "cavity"
                intensity = I;
            elseif id == "outcoupled"
                intensity = I .* (1 - obj.R) ./ (1 + obj.R);
            end
        end

        function power = get_power_total(obj, id)
            % Returns the total optical power of all cavity modes.
            %
            % Syntax:
            %   power = get_power_total(obj, id)
            %
            % Input Arguments:
            %   id (char): Flag specifying if the power inside or outside
            %     the cavity should be returned. Valid choices for id are
            %     ``cavity`` and ``outcoupled``.
            %
            % Output arguments:
            %   power (scalar): Total optical power.

            arguments
                obj
                id{mustBeMember(id, {'cavity', 'outcoupled'})}
            end
            power = sum(obj.get_power_spec(id));
        end

        function power = get_intensity_total(obj, id)
            % Returns the total optical intensity of all cavity modes.
            %
            % Syntax:
            %   power = get_intensity_total(obj, id)
            %
            % Input Arguments:
            %   id (char): Flag specifying if the power inside or outside
            %     the cavity should be returned. Valid choices for id are
            %     ``cavity`` and ``outcoupled``.
            %
            % Output arguments:
            %   intensity (scalar): Total optical intensity.

            arguments
                obj
                id{mustBeMember(id, {'cavity', 'outcoupled'})}
            end
            power = sum(obj.get_intensity_spec(id));
        end

        function plot_power_spec(obj, id)
            % Plots the optical power of the different modes over the
            % frequency.
            %
            % Syntax:
            %   plot_power_spec(obj, id)
            %
            % Input Arguments:
            %   id (char): Flag specifying if the power inside or outside
            %     the cavity should be returned. Valid choices for id are
            %     ``cavity`` and ``outcoupled``.

            power = obj.get_power_spec(id);
            bar(obj.freq/1e12, power*1e3, 0.8, 'EdgeColor', 'w')
            xlabel("Frequency (THz)")
            ylabel("Power (mW)")
            title("Power Spectrum")
        end

        function plot_intensity_spec(obj, id)
            % Plots the intensity of the different modes over the frequency.
            %
            % Syntax:
            %   plot_intensity_spec(obj, id)
            %
            % Input Arguments:
            %   id (char): Flag specifying if the power inside or outside
            %     the cavity should be returned. Valid choices for id are
            %     ``cavity`` and ``outcoupled``.

            I = obj.get_intensity_spec(id);
            bar(obj.freq/1e12, I/1e4, 0.8, 'EdgeColor', 'w')
            xlabel("Frequency (THz)")
            ylabel("Intensity (W/cm^2)")
            title("Intensity Spectrum")
        end

        function to_hdf5(obj, filename)
            % Saves optical_field object in hdf5 format.
            %
            % Syntax:
            %   to_hdf5(obj, filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            write_to_hdf5(filename, "/mode_intensity", obj.intensity);
            write_to_hdf5(filename, "/mode_frequency", obj.freq);
        end
    end

    methods (Static)
        function obj = from_hdf5(filename)
            % Constructs optical_field object from hdf5 file.
            %
            % Syntax:
            %   obj = from_hdf5(filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            I = h5read(filename, "/mode_intensity");
            f = h5read(filename, "/mode_frequency");
            obj = optical_field(I, f);
        end
    end
end
