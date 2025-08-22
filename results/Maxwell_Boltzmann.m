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

classdef Maxwell_Boltzmann
    % Describes the Maxwell_Boltzmann distribution.

    properties
    end

    methods (Static)
        function fE = calc(A, T, delta_E)
            % Returns Maxwell-Boltzmann distribution.
            %
            % Syntax:
            %   fE = calc(A, T, delta_E)
            %
            % Input Arguments:
            %   A (scalar): Amplitude A.
            %   T (scalar): Electron temperature [K].
            %   delta_E (scalar): Energy difference [eV].
            %
            % Output Arguments:
            %   fE (scalar): Evaluated Maxwell-Boltzmann distribution.

            fE = A .* exp((-1)*(phys_const.e0 * delta_E) ...
                ./(phys_const.kB * T));
        end
    end
end
