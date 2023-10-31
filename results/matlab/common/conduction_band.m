classdef conduction_band < handle
    %CONDUCTION_BAND Describes the conduction band profile including space
    % charge effects.
    
    properties (SetAccess = private)
        zv % Position vector in Angstrom.
        Vh % Conduction band potential profile in J.
    end
    
    methods
        function obj = conduction_band(zv, Vh)
            % Position vector of the potential.
            obj.zv = zv;
            % Conduction band potential profile.
            obj.Vh = Vh;
        end
        % Plot conduction band profile with space charge effects
        % and probability densities.
        function plot_profile(obj)
            plot(-(obj.zv / 10), obj.Vh/phys_const.e0, 'Color', ...
                [0, 0, 0], 'linestyle', '--', 'linewidth', 2);
            xlabel('z/nm'); %'Position x in nm'); %
            ylabel('E/eV'); %'Energy in eV'); %
        end
        
    end
end
