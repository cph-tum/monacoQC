classdef (Abstract) bisection_method < handle
    %bisection_method Abstract method class for bisection method.
    methods (Static)
        function Eh = solve(FunFcnIn, x, options)
            while (x(2) - x(1) > options.TolX)
                Eh = 0.5 * (x(1) + x(2));
                bh = FunFcnIn(Eh);
                b = FunFcnIn(x(2));
                if (bh * b < 0)
                    x(1) = Eh;
                else
                    b = bh;
                    x(2) = Eh;
                end
            end
        end
    end
end