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

classdef InAlAs < ternary
    %InAlAs represents InAlAs material.
    %
    properties
    end
    %
    methods
        % Constructs InAlAs.
        function obj = InAlAs(conc)
            if (nargin == 0)
                disp('Generated a InAlAs object with x = 0.52 of InAs.');
                conc = 0.52;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_In must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class AlAs', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class InAs', ...
                        ' to generate an object!']);
                end
            else
                name = 'InAlAs';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = AlAs();
            obj.b2 = InAs();
            obj.c_parab = 1.7;
            % Source: Hendorfer et al. 1993
            % https://doi.org/10.1103/PhysRevB.48.2328
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XIV
            obj.C = parameter(-0.0383, 0.15, 0, 0.7, 0, 0, 0, 0, 0, ...
                -1.4, -4.44, -4.81, 0.049, 0);
        end
    end
end
