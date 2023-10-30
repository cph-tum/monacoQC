classdef GaSb < binary
    %GaSb Represents GaSb material.
    %
    properties
    end
    %
    methods
        % Constructs GaSb.
        function obj = GaSb()
            name = 'GaSb';
            % intialize base class and properties
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table VII
            % http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/mechanic.html
            obj.param = parameter(0.063, 0.76, 15.7, 0.812, 0.8842, ...
                0.4026, 0.4322, 6.0959-4.72e-5*300, ...
                0.8, -7.5, -1.63, 27, 0, -2.0);
            obj.param_T = parameter(0, 0, 0, [0.417e-3, 140], ...
                0, 0, 0, 4.72e-5, 0, 0, 0, 0, 0, 0);
        end
    end
    
end
