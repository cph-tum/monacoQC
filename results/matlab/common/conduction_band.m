%
% monacoQC: An object-oriented Matlab-based device engineering tool for
% quantum cascade structures.
%
% Copyright (C) 2023, Computational Photonics Group, Technical University of
% Munich
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

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
