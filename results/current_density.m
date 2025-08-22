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

classdef current_density < handle
    % Describes the current density in the QCL over time.

    properties
        J; % vector: Current density as a function of time [A/m^2].
        t; % vector: Time vector [s].
        window; % scalar: Time window (number of points) in which the current density is plotted.
        avring_time; % scalar: Start time to calculate the average current density [s].
    end

    methods
        function obj = current_density(t, J)
            % Constructs an object of type current_density.
            %
            % Syntax:
            %   obj = current_density(t, J)
            %
            % Input Arguments:
            %   t (vector): Time vector [s].
            %   J (vector): Current density over time [A/m^2].

            obj.t = t; % Time
            obj.J = J; % Current density
            obj.set_avring_time(9e-12); % Default start time
            obj.set_window(100) % Default window
        end

        function set_avring_time(obj, avring_time)
            obj.avring_time = avring_time;
        end

        function set_window(obj, window)
            obj.window = window;
        end

        function average_curr_dens = moving_average(obj)
            % Calculates the temporal average of the current density for
            % t > avring_time. This method can be used to calculates the
            % steady-state current density.
            %
            % Syntax:
            %   average_curr_dens = moving_average(obj)
            %
            % Output Arguments:
            %   average_curr_dens (scalar): Averaged current density [kA/cm^2].

            % Sum up all densities starting at t > avring_time and average
            % it over the number of added densities.
            average_curr_dens = sum((obj.t > obj.avring_time) ...
                .*obj.J) / ...
                sum((obj.t > obj.avring_time)) / 1e7;
            % Divide by 1e7 to get the average current density in [kA/cm^2]
        end

        function plot_curr_dens(obj)
            % Plots current density over a pre-defined time window.
            %
            % Syntax:
            %   plot_curr_dens(obj)

            I = smooth(obj.J, obj.window);
            p = (1:(length(obj.t) - obj.window)) + obj.window / 2;
            plot(obj.t(p)*1e12, I(p));
            xlabel('t/ps');
            ylabel('J/(A/m^2)');
        end

        function to_hdf5(obj, filename)
            % Saves current_density object in hdf5 format.
            %
            % Syntax:
            %   to_hdf5(obj, filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            write_to_hdf5(filename, "/current_density", obj.J);
            write_to_hdf5(filename, "/time", obj.t);
        end
    end

    methods (Static)
        function obj = from_hdf5(filename)
            % Constructs current_density object from hdf5 file.
            %
            % Syntax:
            %   obj = from_hdf5(filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            J = h5read(filename, "/current_density");
            t = h5read(filename, "/time");
            obj = current_density(t, J);
        end

        function J = curr2(sc, carr_dist, d)
            % Static method for calculating the current density from one
            % period to the next one.
            %
            % Syntax:
            %   J = curr2(sc, carr_dist, d)
            %
            % Input Arguments:
            %   sc (scattering_rates-object): Contains information about
            %     the scattering rates between the individual subbands.
            %   carr_dist (carrier_distribution-object): Contains
            %     information about the carrier distributions within each
            %     subband.
            %   d (device-object): Contains information about the structure,
            %     geometry and materials of the QCL.
            %
            % Output Arguments:
            %   J (scalar): Current density [kA/m^2].

            % Sheet density.
            n2D = d.dens_sheet;
            % Mean value of the occupations of the two middle periods.
            occ = carr_dist.occupation;
            occ = reshape(occ, length(occ), []); % line i contains occupation of ith level
            % Number of wavefunctions in one period.
            num_wfs = length(occ);
            % Index of wavefunctions in the observed period.
            ind_wfs = (num_wfs + 1):(2 * num_wfs);
            % Scattering rates form the observed period to the right
            % (e.g. from 1 to 2).
            scat_right = sc.get_scattering_matrix('right', ind_wfs); % entry ij contains rate i->j
            % Scattering rates form the observed period to the left
            % (e.g. from 2 to 1).
            scat_left = sc.get_scattering_matrix('left', ind_wfs); % entry ij contains rate i->j

            % Carriers moving to the next period.
            right = scat_right .* occ;
            % Carriers moving to the previous period.
            left = scat_left .* occ;

            % J %in kA/m^2 if n2D is in 1/m^2.
            % right: +lvl, left: -lvl.
            J = -sum(sum(right-left)) * n2D * phys_const.e0 / 1000;
        end
    end
end
