classdef (Abstract) solver < handle
    %solver Abstract solver class for solving the Schrödinger equation.
    properties (SetAccess = protected)
        name = '' % Solver name
    end
    
    methods
        % Constructs solver.
        function obj = solver(name)
            obj.name = name;
        end
    end
    methods (Abstract)
        [eig_system, cond_profile] = solve(obj, func_carr_distr, ...
            eig_system_old);
    end
end
