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

classdef laser_waveguide < handle
    % Describes the properties of a laser waveguide.

    properties (SetAccess = protected)
        l_waveguide = 0 % scalar: Length of waveguide [m].
        w_waveguide = 0 % scalar: Width of waveguide [m].
        h_waveguide = 0 % Height of waveguide [m].
        refl_coeff_left = 1 % scalar: Reflectivity coefficient r at left facet.
        refl_coeff_right = 1 % scalar: Reflectivity coefficient r at right facet.
        overlap_factor = 1 % scalar: Overlap factor.
        a_waveguide = [] % scalar: Waveguide loss coefficient in [m^-1].
    end

    properties (Dependent)
        a_power % scalar: Power loss coefficient [m^-1].
        a_field % scalar: Field loss coefficient [m^-1].
        a_mirror % scalar: Mirror loss coefficient [m^-1].
        A_act % Cross-sectional area of the active region [m^2].
        refl_left % Power reflectance R at left facet.
        refl_right % Power reflectance R at right facet.
    end

    methods
        function obj = laser_waveguide(l_wg, loss_wg, h_wg, ...
                overlap_f, refl_coeff_left, refl_coeff_right, w_wg)
            % Constructs an object of type laser_waveguide.
            %
            % Syntax:
            %   obj = laser_waveguide(l_wg, loss_wg, h_wg)
            %   obj = laser_waveguide(l_wg, loss_wg, h_wg, overlap_f)
            %   obj = laser_waveguide(l_wg, loss_wg, h_wg, overlap_f, refl_coeff_left)
            %   obj = laser_waveguide(l_wg, loss_wg, h_wg, overlap_f, refl_coeff_left, refl_coeff_right)
            %   obj = laser_waveguide(l_wg, loss_wg, h_wg, overlap_f, refl_coeff_left, refl_coeff_right, w_wg)
            %
            % Input Arguments:
            %   l_wg (scalar): Lenght of waveguide.
            %   loss_wg (scalar): Waveguide loss coefficient.
            %   h_wg (scalar): Height of waveguide.
            %   overlap_f (scalar): Overlap factor gamma.
            %   refl_coeff_left (scalar): Field reflectivity coefficient of
            %     the left facet.
            %   refl_coeff_right (scalar): Field reflectivity coefficient
            %     of the right facet.
            %   w_wg (scalar): Width of waveguide.

            % Optical device length.
            obj.l_waveguide = l_wg;
            % Waveguide power loss coefficient.
            obj.a_waveguide = loss_wg;
            % Optical device height.
            obj.h_waveguide = h_wg;
            if (nargin > 3)
                % Overlap factor.
                obj.overlap_factor = overlap_f;
                if (nargin > 4)
                    % Field reflectivity at left facet.
                    obj.refl_coeff_left = refl_coeff_left;
                    if (nargin > 5)
                        % Field reflectivity at right facet.
                        obj.refl_coeff_right = refl_coeff_right;
                        if (nargin == 7)
                            obj.w_waveguide = w_wg;
                        end
                    end
                end
            end
        end

        function a_m = get.a_mirror(obj)
            % Gets mirror power loss coefficient.

            % Facet reflectance.
            R = obj.refl_left * obj.refl_right;
            % Mirror power loss coefficient.
            a_m = (-1) * log(R) / obj.l_waveguide;
        end

        function a = get.a_power(obj)
            % Gets power loss coefficient.
            a = obj.a_mirror + obj.a_waveguide;
        end

        function a_f = get.a_field(obj)
            % Gets field loss coefficient.

            % Facet reflectivity.
            R = obj.refl_coeff_left * obj.refl_coeff_right;
            % Mirror field loss coefficient.
            a_f = (-1) * log(R) / obj.l_waveguide;
        end

        function A = get.A_act(obj)
            % Gets cross-sectional area of the active region.
            A = obj.w_waveguide * obj.h_waveguide;
        end

        function r_left = get.refl_left(obj)
            % Gets reflectance at left facet.
            r_left = abs(obj.refl_coeff_left)^2;
        end

        function r_right = get.refl_right(obj)
            % Gets reflectance at right facet
            r_right = abs(obj.refl_coeff_right)^2;
        end
    end
end
