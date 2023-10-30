classdef GaAs_AlGaAs < material_system
    %GaAs_AlGaAs Represents material system GaAs/AlGaAs.
    %
    properties (SetAccess = protected)
    end
    %
    methods
        % Constructs material system.
        function obj = GaAs_AlGaAs(name)
            obj = obj@material_system(name);
            % Interface roughness parameter
            obj.Gamma = 1e-8;
            obj.Delta = 1.2e-10;
            % Material system flag
            obj.flag_material_system = 2;
        end
    end
end
