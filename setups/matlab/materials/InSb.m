classdef InSb < binary
    % InSb Represents InSb material.
    %
    properties
    end
    %
    methods
        % Constructs InSb.
        function obj = InSb()
            name = 'InSb';
            % intialize base class and properties
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table VIII
            % http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/mechanic.html
            obj.param = parameter(-0.484, 0.81, 16.8, 0.235, 0.6847, ...
                0.3735, 0.3111, 6.4794-3.48e-5*300, ...
                0.36, -6.94, -0.23, 23.3, 0, -2);
            obj.param_T = parameter(0, 0, 0, [0.32e-3, 170], ...
                0, 0, 0, 3.48e-5, 0, 0, 0, 0, 0, 0);
        end
    end
    
end
