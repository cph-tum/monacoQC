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
