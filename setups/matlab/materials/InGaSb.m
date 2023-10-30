classdef InGaSb < ternary
    %InGaSb represents InGaSb material.
    %
    properties
    end
    %
    methods
        % Constructs InGaSb.
        function obj = InGaSb(conc)
            if (nargin == 0)
                disp('Generated a InGaSb object with x = 0.5 of InSb.');
                conc = 0.5;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_In must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class GaSb', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class InSb', ...
                        ' to generate an object!']);
                end
            else
                name = 'InGaSb';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = GaSb();
            obj.b2 = InSb();
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XVIII
            % Source: Ioffe.ru
            obj.C = parameter(0.415, 0.1, 0, 0.415, 0, 0, 0, 0, 0, ...
                2.61, -6.84, 0, 0.0091, 0);
        end
    end
end
