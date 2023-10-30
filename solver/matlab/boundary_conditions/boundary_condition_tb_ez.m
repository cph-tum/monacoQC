classdef boundary_condition_tb_ez < boundary_condition_tb
    %boundary_condition_ez Boundary condition class for ez eigenstates.
    
    properties
        e_multiplet % Maximum energy spacing of multiplets.
        d_multiplet % Minimum dipole moment of multiplets.
    end
    
    methods
        function obj = boundary_condition_tb_ez(d, s, sim_c)
            % Creates an object with properties describing the simulation
            % boundary conditions for ez state solutions.
            obj = obj@boundary_condition_tb(d, s, sim_c);
            obj.e_multiplet = s.e_multiplet;
            obj.d_multiplet = s.d_multiplet;
        end
        
        
    end
    methods (Static)
        function bc_start = get_bc_start()
            bc_start = [0; 1];
        end
    end
end