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

classdef InGaSb < ternary
    %InGaSb represents InGaSb material.
    %
    properties
    end
    %
    methods
        % Constructs InGaSb.
        function obj = InGaSb(conc)
            if (nargin == 0)
                disp('Generated a InGaSb object with x = 0.5 of InSb.');
                conc = 0.5;
            end
            % Reasonability check
            if ((conc <= 0) || (conc >= 1))
                % Out of scope
                if ((conc < 0) || (conc > 1))
                    error(['Parameter conc_In must be', ...
                        'in the interval ]0, 1[!']);
                    % Check Ternary material
                elseif (conc == 0)
                    error(['Please use the binary alloy class GaSb', ...
                        ' to generate an object!']);
                elseif (conc == 1)
                    error(['Please use the binary alloy class InSb', ...
                        ' to generate an object!']);
                end
            else
                name = 'InGaSb';
            end
            % intialize base class and properties at temperature t
            obj = obj@ternary(name);
            obj.conc = conc;
            obj.b1 = GaSb();
            obj.b2 = InSb();
            % Bowing parameters
            % Source: Vurgaftman et al. 2001, Table XVIII
            % Source: Ioffe.ru
            obj.C = parameter(0.415, 0.1, 0, 0.415, 0, 0, 0, 0, 0, ...
                2.61, -6.84, 0, 0.0091, 0);
        end
    end
end
