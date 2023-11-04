%
% monacoQC: An object-oriented Matlab-based device engineering tool for
% quantum cascade structures.
%
% Copyright (C) 2023, Computational Photonics Group, Technical University of
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

classdef phys_const
    %phys_const Defines all necessary physical constants in SI units
    % used in the monacoQC project.
    %
    properties
    end
    %
    methods (Static)
        function E0 = e0()
            % Returns elementary charge.
            E0 = 1.60217646e-19;
        end
        function EPS0 = eps0()
            % Returns vacuum permittivity.
            EPS0 = 8.854187817e-12;
        end
        function MU0 = mu0()
            % Returns vacuum permeability.
            MU0 = pi * 4e-7;
        end
        function C0 = c0()
            % Returns vacuum speed of light.
            C0 = 1 / sqrt(phys_const.mu0*phys_const.eps0);
        end
        function KB = kB()
            % Returns Boltzmann constant.
            KB = 1.38064852e-23;
        end
        function HBAR = hbar()
            % Returns reduced Planck's constant.
            HBAR = 1.05457266e-34;
        end
        function H = h()
            % Returns Planck's constant.
            H = 2 * pi * phys_const.hbar;
        end
        function ME = me()
            % Returns electron mass.
            ME = 9.10938188e-31;
        end
    end
end
