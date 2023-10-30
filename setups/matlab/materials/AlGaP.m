classdef AlGaP < ternary
    %AlGaP represents AlGaP material.
    %
    properties
    end
    %
    methods
        % Constructs AlGaP.
        function obj = AlGaP(conc)
            if (nargin == 0)
                disp('Generated a AlGaP object with x = 0.5 of AlP.');
                conc = 0.5;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_Al must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class GaP', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class AlP', ...
                        ' to generate an object!']);
                end
            else
                name = 'AlGaP';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = GaP();
            obj.b2 = AlP();
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XVII
            obj.C = parameter(0, 0, 0, 0, 0, 0, ...
                0, 0, 0, 0, 0, 0, 0, 0);
        end
    end
end
