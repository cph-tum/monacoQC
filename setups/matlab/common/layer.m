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

classdef layer < matlab.mixin.Copyable
    %layer Represents a layer in the active region.
    %
    properties
        material = ''; % Material identifier string
        length = 0.0; % Length of layer in Angstrom
        doping = 0.0; % Doping concentration in 1/cm^3
        doping_type = 'n'; % Doping type (n or p)
        orientation = '001'; % Interface orientation
        layer_type; % Layer type: barrier or well
    end
    %
    methods
        % Constructs layer.
        function obj = layer(material, length, doping)
            
            if (nargin < 2)
                error('Too few arguments provided in layer().');
            elseif (nargin == 2)
                doping = 0.0;
            end
            %
            % reasonability check
            if (length < 0)
                error('Parameter length must be positive!');
            end
            if (doping < 0)
                error('Parameter doping must be positive!');
            end
            %
            % intialize properties
            obj.material = material;
            obj.length = length;
            obj.doping = doping;
        end
        
        function obj = set.orientation(obj, indizes)
            if (ischar(indizes))
                obj.orientation = indizes;
            else
                error('Indizes must be a character array, e.g. 001.')
            end
        end
        function orientation = get_orientation(obj)
            orientation = obj.orientation;
        end
        %
    end
end
