classdef AlGaAs < ternary
    %AlGaAs represents AlGaAs material.
    %
    properties
    end
    %
    methods
        % Constructs AlGaAs.
        function obj = AlGaAs(conc)
            if (nargin == 0)
                disp('Generated a AlGaAs object with x = 0.15 of AlAs.');
                conc = 0.15;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_Al must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class GaAs', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class AlAs', ...
                        ' to generate an object!']);
                end
            else
                name = 'AlGaAs';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = GaAs();
            obj.b2 = AlAs();
            obj.c_parab = 3.2;
            % Source: Ekenberg et al. 1989, Hendorfer et al. 1993
            % https://doi.org/10.1103/PhysRevB.48.2328
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XII
            obj.C = parameter(0, 0, 0, -0.127+1.31*obj.conc, 0, 0, ...
                0, 0, 0, 0, 0, 0, 0, 0);
        end
        % Calculates conduction band offset (special case for AlGaAs).
        % Ec from Yi et al, PRB 81, 235325 (2010).
        function ec = get_Ec0(obj)
            if (obj.conc < 0.42)
                ec = 0.831 * obj.conc;
            else
                ec = -0.115 + 1.105 * obj.conc;
            end
        end
    end
end
