classdef InAlP < ternary
    %InAlP represents InAlP material.
    %
    properties
    end
    %
    methods
        % Constructs InAlP.
        function obj = InAlP(conc)
            if (nargin == 0)
                disp('Generated a InAlP object with x = 0.48 of InP.');
                conc = 0.48;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_In must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class AlP', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class InP', ...
                        ' to generate an object!']);
                end
            else
                name = 'InAlP';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = AlP();
            obj.b2 = InP();
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XVI
            obj.C = parameter(-0.48, -0.19, 0, -0.48, 0, 0, 0, 0, 0, ...
                0, 0, 0, 0.22, 0);
        end
    end
end
