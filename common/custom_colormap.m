%
% monacoQC: An object-oriented Matlab-based device engineering tool for
% quantum cascade structures.
%
% Copyright (C) 2025, Computational Photonics Group, Technical University of
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

classdef custom_colormap
    % Defines the colormap for plot generation in Monaco.

    methods (Static)
        function B = colormap_old(num_wfs)
            % Returns old colormap for given number of states.
            %
            % Syntax:
            %   B = colormap_old(num_wfs)
            %
            % Input Arguments:
            %   num_wfs (scalar): Number of different colors included in
            %     the colormap. The maximum number is 12.
            %
            % Output Arguments:
            %   B (matrix): RGB color values for the selected number of
            %     states.

            colormap_old = [0, 0, 1; 0, 0.5, 0; 1, 0, 0; 0, 0.75, 0.75; ...
                0.75, 0, 0.75; 0.75, 0.75, 0; 0.25, 0.25, 0.25; ...
                0, 1, 0; 0.5, 0.5, 0.5; 0.5, 0.5, 0; 1, 0.5, 0; 1, 0.7, 0];
            B = colormap_old(1:num_wfs, :);
        end

        function B = colormap_TUM(num_wfs)
            % Returns colormap with TUM colors for given number of states.
            %
            % Syntax:
            %   B = colormap_TUM(num_wfs)
            %
            % Input Arguments:
            %   num_wfs (scalar): Number of different colors to be included
            %     in the colormap. The maximum number is 13.
            %
            % Output Arguments:
            %   B (matrix): RGB color values for the selected number of
            %     states.

            colormap_TUM = [0.4118, 0.0314, 0.3529; 0, 0.2000, 0.3490; ...
                0, 0.3216, 0.5765; 0, 0.3961, 0.7412; ...
                0.3922, 0.6275, 0.7843; 0.5961, 0.7765, 0.9176; ...
                0.6118, 0.6157, 0.6235; 0.8510, 0.8549, 0.8588; ...
                0.8549, 0.8431, 0.7961; 0.6353, 0.6784, 0; ...
                1.0000, 0.8627, 0; 0.8902, 0.4471, 0.1333; ...
                0.6118, 0.0510, 0.0863];
            B = colormap_TUM(1:num_wfs, :);
        end
    end
end
