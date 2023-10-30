classdef AlGaSb < ternary
    %AlGaSb represents AlGaSb material.
    %
    properties
    end
    %
    methods
        % Constructs AlGaSb.
        function obj = AlGaSb(conc)
            if (nargin == 0)
                disp('Generated a AlGaSb object with x = 0.5 of AlSb.');
                conc = 0.5;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_Sb must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class GaSb', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class AlSb', ...
                        ' to generate an object!']);
                end
            else
                name = 'AlGaSb';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = GaSb();
            obj.b2 = AlSb();
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XX
            obj.C = parameter(-0.044+1.22*obj.conc, 0.3, 0, ...
                -0.044+1.22*obj.conc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        end
    end
end
