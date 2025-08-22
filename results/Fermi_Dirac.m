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

classdef Fermi_Dirac
    % Describes the Fermi Dirac distribution.

    methods (Static)
        function fE = calc(A, T, mu, E)
            if nargin < 4
                A = 1;
            end
            % Returns Fermi-Dirac distribution.
            %
            % Syntax:
            %   fE = calc(A, T, mu, E)
            %
            % Input Arguments:
            %   A (scalar): Amplitude A.
            %   T (scalar): Electron temperature [K].
            %   mu (scalar): Chemical potential [eV].
            %   E (scalar): Energy [eV].
            %
            % Output Arguments:
            %   fE (scalar): Evaluated Fermi-Dirac distribution.

            fE = A * (exp(phys_const.e0*(E - mu)./(phys_const.kB * T)) ...
                +1).^(-1);
        end
    end
end
