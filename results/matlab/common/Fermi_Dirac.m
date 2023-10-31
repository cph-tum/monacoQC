classdef Fermi_Dirac
    %Fermi_Dirac Describes the Fermi Dirac distribution.
    
    methods (Static)
        function fE = calc(A, T, mu, E)
            if nargin < 4
                A = 1;
            end
            % Returns Fermi-Dirac distribution for given energy E,
            % chemical potential mu and electron temperature T.
            % A represents Amplitude.
            fE = A * (exp(phys_const.e0*(E - mu)./(phys_const.kB * T)) ...
                +1).^(-1);
        end
        
    end
end