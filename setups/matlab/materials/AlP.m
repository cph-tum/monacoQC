classdef AlP < binary
    %InAs Represents InAs material.
    %
    properties
    end
    %
    methods
        % Constructs InP.
        function obj = AlP()
            name = 'AlP';
            % intialize base class and properties at temperature t
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table VI
            % Source: aftershoq
            obj.param = parameter(1.171, 0.07, 9.8, 3.63, 1.330, ...
                0.63, 0.615, 5.4672-2.92e-5*300, ...
                3, -5.7, -0.65, 17.7, 0, -1.5);
            obj.param_T = parameter(0, 0, 0, [0.5771e-3, 372], ...
                0, 0, 0, 2.92e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
