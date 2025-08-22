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

classdef (Abstract) bisection_method < handle
    % Abstract method class for bisection method.

    methods (Static)
        function Eh = solve(FunFcnIn, x, options)
            % Finds the zeros for a given function.
            %
            % Syntax:
            %   Eh = solve(FunFcnIn, x, options)
            %
            % Input Arguments:
            %   FunFcnIn (function-handle):
            %   x (vector): 2-by-1 vector containing the bounds of search
            %     interval.
            %   options (struct): Should contain the field ``TolX``, which
            %     defines the absolute error to terminate the loop.
            %
            % Output Arguments:
            %   Eh (scalar): Returns zero.

            while (x(2) - x(1) > options.TolX)
                Eh = 0.5 * (x(1) + x(2));
                bh = FunFcnIn(Eh);
                b = FunFcnIn(x(2));
                if (bh * b < 0)
                    x(1) = Eh;
                else
                    b = bh;
                    x(2) = Eh;
                end
            end
        end
    end
end
