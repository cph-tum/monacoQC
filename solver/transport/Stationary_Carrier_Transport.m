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

classdef (Abstract) Stationary_Carrier_Transport
    % Abstract class for self-consistent simulations of stationary carrier
    % transport in QCLs.

    properties
        max_iteration % scalar | 15 (default): Max iteration of self-consistent loop.
        relTol % scalar | 0.001 (default): Max relative error of self-consistent loop.
        device % device-object: Device description.
        scenario % scenario-object: Scenario description.
    end

    methods
        function obj = Stationary_Carrier_Transport(d, s)
            % Constructs an object of type Stationary_Carrier_Transport.
            %
            % Syntax:
            %   obj = Stationary_Carrier_Transport(d, s)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   s (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.

            obj.device = d;
            obj.scenario = s;
            obj.max_iteration = 15;
            obj.relTol = 1e-3;
        end

        function obj = set_solver_options(obj, options)
            % Sets options related to the stationary carrier transport
            % solver.
            %
            % Syntax:
            %   obj = set_solver_options(obj)
            %   obj = set_solver_options(__, name, value)
            %
            % Name Value Input:
            %   max_iteration (scalar): Max iteration of self-consistent
            %     loop.
            %   relTol (scalar): Desired relative error between two
            %     subsequent iterations to terminate the loop.

            arguments
                obj
                options.max_iteration
                options.relTol
            end

            % set properties
            keys = fieldnames(options);
            for i = 1:length(keys)
                key = keys{i};
                obj.(key) = options.(key);
            end
        end
    end

    methods (Abstract)
        solve(obj) % Solves carrier transport.
    end
end
