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

classdef AlAsSb < ternary
    %AlAsSb represents AlAsSb material.
    %
    properties
    end
    %
    methods
        % Constructs AlAsSb.
        function obj = AlAsSb(conc)
            if (nargin == 0)
                disp('Generated a AlAsSb object with x = 0.5 of AlAs.');
                conc = 0.5;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_Al must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class AlSb', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class AlAs', ...
                        ' to generate an object!']);
                end
            else
                name = 'AlAsSb';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = AlSb();
            obj.b2 = AlAs();
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XX
            obj.C = parameter(-0.91, 0.15, 0, ...
                0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        end
    end
end
