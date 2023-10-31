classdef current_density < handle
    % current_density contains the change in current density over time and
    % return the steady state current density of the QCL.
    properties
        J; % Vector of the current density in [A/m^2]
        t; % Vector of the time in [s]
        window; % Window wich is showen in the density plot
        avring_time; % Start time to calculate the avg current [s]
    end
    
    methods
        function obj = current_density(t, J)
            obj.t = t; % Time
            obj.J = J; % Current density
            obj.set_avring_time(9e-12); % Default start time
            obj.set_window(100) % Default window
        end
        
        function set_avring_time(obj, avring_time)
            obj.avring_time = avring_time;
        end
        
        function set_window(obj, window)
            obj.window = window;
        end
        
        function average_curr_dens = moving_average(obj)
            % Sum up all densities starting at t > avring_time and average
            % it over the number of added densities.
            average_curr_dens = sum((obj.t > obj.avring_time) ...
                .*obj.J) / ...
                sum((obj.t > obj.avring_time)) / 1e7;
            % Divide by 1e7 to get the average current density in [kA/cm^2]
        end
        
        function plot_curr_dens(obj)
            I = smooth(obj.J, obj.window);
            p = (1:(length(obj.t) - obj.window)) + obj.window / 2;
            plot(obj.t(p)*1e12, I(p));
            xlabel('t/ps');
            ylabel('J/(A/m^2)');
        end
        
        function J = curr2(obj, sc, carr_dist, d)
            % Sheet density.
            n2D = d.dens_sheet;
            % Mean value of the occupations of the two middle periods.
            occ = carr_dist.occupation;
            % Number of wavefunctions in one period.
            num_wfs = length(occ);
            % Index of wavefunctions in the observed period.
            ind_wfs = (num_wfs + 1):(2 * num_wfs);
            % Scattering rates form the observed period to the right
            % (e.g. from 1 to 2).
            scat_right = sc.get_scattering_matrix('right', ind_wfs);
            % Scattering rates form the observed period to the left
            % (e.g. from 2 to 1).
            scat_left = sc.get_scattering_matrix('left', ind_wfs);
            
            % Carriers moving to the next period.
            right = scat_right .* occ;
            % Carriers moving to the previous period.
            left = scat_left .* occ;
            
            % J %in kA/m^2 if n2D is in 1/m^2.
            % right: +lvl, left: -lvl.
            J = -sum(sum(right-left)) * ...
                n2D * phys_const.e0 / 1000;
            
            
        end
    end
end
