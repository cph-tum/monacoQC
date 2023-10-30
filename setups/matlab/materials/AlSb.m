classdef AlSb < binary
    %AlSb Represents AlSb material.
    %
    properties
    end
    %
    methods
        % Constructs AlSb.
        function obj = AlSb()
            name = 'AlSb';
            % intialize base class and properties
            obj = obj@binary(name);
            % Material parameter
            % Source: Vurgaftman et al. 2001, cf. Table VIII
            % http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/mechanic.html
            obj.param = parameter(1.257, 0.676, 12.04, 2.386, 0.8769, ...
                0.4341, 0.4076, 6.1355-2.6e-5*300, ...
                1.4, -4.5, -0.56, 18.7, 0, -1.35);
            obj.param_T = parameter(0, 0, 0, [0.42e-3, 140], ...
                0, 0, 0, 2.6e-5, 0, 0, 0, 0, 0, 0);
        end
    end
    
end
