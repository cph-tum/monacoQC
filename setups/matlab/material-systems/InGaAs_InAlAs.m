classdef InGaAs_InAlAs < material_system
    %InGaAs_InAlAs Represents material system InGaAs/InAlAs.
    %
    properties (SetAccess = protected)
    end
    %
    methods
        % Constructs material system.
        function obj = InGaAs_InAlAs(name)
            obj = obj@material_system(name);
            % Interface roughness parameter
            obj.Gamma = 1e-8;
            obj.Delta = 1e-18 / obj.Gamma;
            % Material system flag
            obj.flag_material_system = 1;
        end
    end
end
