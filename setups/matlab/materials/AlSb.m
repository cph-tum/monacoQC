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
