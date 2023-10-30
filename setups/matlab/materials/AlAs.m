classdef AlAs < binary
    %AlAs Represents AlAs material.
    properties
    end
    %
    methods
        % Constructs AlAs.
        function obj = AlAs()
            name = 'AlAs';
            % intialize base class and properties
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table II
            obj.c_parab = 0;
            obj.param = parameter(0.99, 0.280, 10.06, 3.099, 1.250, ...
                0.534, 0.542, 5.6611-2.90e-5*300, ...
                2.47, -5.64, -0.48, 21.1, 0, -2.3);
            obj.param_T = parameter(0, 0, 0, [0.885e-3, 530], ...
                0, 0, 0, 2.9e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
