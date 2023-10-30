classdef InAs < binary
    %InAs Represents InAs material.
    %
    properties
    end
    %
    methods
        % Constructs InAs.
        function obj = InAs()
            name = 'InAs';
            % intialize base class and properties at temperature t
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table III
            obj.c_parab = 0;
            obj.param = parameter(-0.892, 0.390, 15.1, 0.417, 0.833, ...
                0.453, 0.396, 6.0583-2.74e-5*300, ...
                1.00, -5.08, -2.90, 21.5, 0, -1.8);
            % Source elastic constants C11, C12, C44: Ioffe
            obj.param_T = parameter(0, 0, 0, [0.276e-3, 94], ...
                0, 0, 0, 2.74e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
