classdef material_system < handle
    %material_system TODO add description.
    %
    properties (SetAccess = protected)
        name = ''; % Material system name
        Gamma = 0; % Interface roughness corr. length in m
        Delta = 0; % Interface roughness height in m
        flag_material_system = 0;
    end
    %
    methods
        % Constructs material system.
        function obj = material_system(name)
            obj.name = name;
        end
        %
        % Gets interface roughness corr. length.
        function gamma = get.Gamma(obj)
            gamma = obj.Gamma;
        end
        %
        % Gets interface roughness height.
        function delta = get.Delta(obj)
            delta = obj.Delta;
        end
    end
end
