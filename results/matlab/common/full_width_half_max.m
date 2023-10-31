classdef full_width_half_max
    % Full_width_half_max calculates the full width at half maximum.
    properties
    end
    %
    methods (Static)
        function fwhm = calc(peak)
            % Returns the full width at half maximum for a given peak.
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
