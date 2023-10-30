classdef waveguide < handle
    %cavity Represents an optical cavity.
    %
    properties (SetAccess = protected)
        l_waveguide = 0; % Length of waveguide in m
        w_waveguide = 0; % Width of waveguide in m
        h_waveguide = 0; % Height of waveguide in m
        refl_left = 1; % Facet reflectance at left facet
        refl_right = 1; % Facet reflectance at right facet
        overlap_factor = 1; % Overlap factor
        a_waveguide = []; % Waveguide loss coefficient in m^-1
    end
    properties (Dependent)
        a_power % Power loss coefficient in m^-1
        a_field % Field loss coefficient in m^-1
        a_mirror % Mirror loss coefficient in m^-1
        A_act % Cross-sectional area of the active region.
    end
    %
    methods
        % Constructs cavity.
        function obj = waveguide(l_wg, loss_wg, ...
                h_wg, overlap_f, refl_left, refl_right, w_wg)
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
                    % Facet reflectance left.
                    obj.refl_left = refl_left;
                    if (nargin > 5)
                        % Facet reflectance right.
                        obj.refl_right = refl_right;
                        if (nargin == 7)
                            obj.w_waveguide = w_wg;
                        end
                    end
                end
            end
        end
        %
        % Get mirror power loss coefficient.
        function a_m = get.a_mirror(obj)
            % Facet refelctance.
            R = obj.refl_left * obj.refl_right;
            % Mirror power loss coefficient.
            a_m = (-1) * log(R) / obj.l_waveguide;
        end
        % Get power loss coefficient.
        function a = get.a_power(obj)
            a = obj.a_mirror + obj.a_waveguide;
        end
        % Get field loss coefficient.
        function a_f = get.a_field(obj)
            a_f = obj.a_power / 2;
        end
        % Get cross-sectional area of the active region.
        function A = get.A_act(obj)
            A = obj.w_waveguide * obj.h_waveguide;
        end
    end
    %
end
