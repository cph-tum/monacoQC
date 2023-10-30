classdef InGaP < ternary
    %InGaP represents InGaP material.
    %
    properties
    end
    %
    methods
        % Constructs InGaP.
        function obj = InGaP(conc)
            if (nargin == 0)
                disp('Generated a InGaP object with x = 0.49 of InP.');
                conc = 0.49;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_In must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class GaP', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class InP', ...
                        ' to generate an object!']);
                end
            else
                name = 'InGaP';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = GaP();
            obj.b2 = InP();
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XV
            % Source: Ioffe.ru
            % Source: aftershoq
            obj.C = parameter(0.65, 0, 0, 0.65, 0, 0, 0, 0, 0, ...
                2.61, 0.78, 0, 0.051, 0);
        end
    end
end
