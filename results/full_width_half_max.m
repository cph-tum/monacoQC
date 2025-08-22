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

classdef full_width_half_max
    % Calculates the full width at half maximum.

    methods (Static)
        function fwhm = calc(peak)
            % Returns the full width at half maximum for a given peak.
            %
            % Syntax:
            %   fwhm = calc(peak)
            %
            % Input Arguments:
            %   peak (vector): Vector including one peak.
            %
            % Output Arguments:
            %   fwhm (scalar): Full-width at half maximum.

            [pmax, imax] = max(peak);
            % Initial index left of peak.
            il = imax;
            % Initial index right of peak.
            ir = imax;
            % Iterate index left until first point below half maximum.
            while (peak(il) > 0.5 * pmax)
                il = il - 1;
                if (il == 0)
                    break;
                end
            end
            % Iterate index right until first point below half maximum.
            while (peak(ir) > 0.5 * pmax)
                ir = ir + 1;
                if (ir > length(peak))
                    break;
                end
            end
            if ((il > 0) && (il < length(peak)) && ...
                    (ir <= length(peak)) && (ir > 1))
                % Difference between first point above and below half
                % maximum.
                ml = abs(peak(il+1)) - abs(peak(il));
                % Difference between first point below and above half
                % maximum.
                mr = abs(peak(ir)) - abs(peak(ir-1));
                % Consider linear interpolation between first point below
                % and above half maximum. Determine the fractional
                % dependency between half maximum (pmax * 0.5) and first
                % point below (abs(peak(ir))) by dividing through the
                % difference between first point above and below half
                % maximum (/mr | /ml).
                % This fraction can be added to the corresponding index due
                % to linearity.
                fwhm = ((pmax * 0.5 - abs(peak(ir))) / mr + ir) - ...
                    ((pmax * 0.5 - abs(peak(il))) / ml + il);
            else
                fwhm = NaN;
            end
        end
    end
end
