classdef GaP < binary
    %GaP Represents GaP material.
    %
    properties
    end
    %
    methods
        % Constructs GaP.
        function obj = GaP()
            name = 'GaP';
            % intialize base class and properties at temperature t
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table IV
            obj.param = parameter(0.361, 0.08, 11.1, 2.35, 1.405, ...
                0.6203, 0.7033, 5.4505-2.92e-5*300, ...
                1.7, -8.2, -2.04, 31.4, 0, -1.6);
            obj.param_T = parameter(0, 0, 0, [0.5771e-3, 372], ...
                0, 0, 0, 2.92e-5, 0, 0, 0, 0, 0, 0);
        end
    end
end
