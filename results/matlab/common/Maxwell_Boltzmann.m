classdef Maxwell_Boltzmann
    %Maxwell_Boltzmann Describes the Maxwell_Boltzmann distribution.
    %
    properties
    end
    %
    methods (Static)
        function fE = calc(A, T, delta_E)
            % Returns Maxwell-Boltzmann distribution.
            % Maxwell-Boltzmann distribution for given amplitude A,
            % energie difference delta_E and electron temperature T.
            fE = A .* exp((-1)*(phys_const.e0 * delta_E) ...
                ./(phys_const.kB * T));
        end
        
    end
end
