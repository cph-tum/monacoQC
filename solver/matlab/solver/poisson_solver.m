classdef poisson_solver < handle
    %poisson_solver Solver class to solve Poisson equation.
    properties
        
    end
    
    methods (Static)
        % Calculate additional electrostatic potential arising from space
        % charge effects.
        function DV = calc_potential(sim_c, rho)
            % Get indices for one device period.
            ind_period = (sim_c.l_period - 1) / sim_c.dz_poisson;
            eps = phys_const.eps0 * sim_c.rel_permittivity;
            M = zeros(ind_period, ind_period);
            dz = sim_c.z_unit * sim_c.dz_poisson;
            % Finite Difference method
            for mz = 1:ind_period
                M(mz, mz) = -2 / dz^2;
                if (mz > 1)
                    M(mz, mz-1) = 1 / dz^2;
                end
                if (mz < ind_period)
                    M(mz, mz+1) = 1 / dz^2;
                end
            end
            % Solve Poisson equation, see Eq. 20-23 in Jirauschek and Kubis
            % 2014 (https://aip.scitation.org/doi/abs/10.1063/1.4863665).
            DVh = transpose(M \ (phys_const.e0 / eps ...
                * rho(2:ind_period+1)));
            DVh = [0, DVh];
            % Extend DVh to the length of number of periods + 1.
            DVh = repmat(DVh, 1, sim_c.num_periods_wf+2);
            % DVh vector includes all periods plus one terminating barrier.
            DVh = DVh(1:length(sim_c.vec_z_poisson));
            DVh = DVh - mean(DVh);
            % DV vector should equal obj.sim_const.vec_z_tm vector in size.
            DV = interp1(sim_c.vec_z_poisson, DVh*phys_const.e0, ...
                sim_c.vec_z_tm);
        end
        % Get space charge density rho for given eigenstates
        % and carrier distribution.
        function rho = get_rho(sim_c, fun_carr_distr, eig)
            % Electron density depends on eigenstates object, if
            % empty assume uniform distribution.
            if (isempty(eig))
                rhoe = sim_c.vec_rhod * 0 + 1;
            else
                % Get carrier distribution from function handle with
                % eigenstates as input argument.
                carr_distr = fun_carr_distr(eig);
                % Calculate electron density using eigenstates and
                % corresponding carrier distribution.
                % For the interpolation use second wavefunction period, to
                % be sure there is no wavefunction truncation on the
                % simulation boundary.
                rhoe1 = sum(((eig.psi(:, ...
                    (eig.num_wfs + 1):(2 * eig.num_wfs)).^2) ...
                    * diag(carr_distr.occupation((1): ...
                    (eig.num_wfs))))');
                rhoe = sim_c.vec_rhod * 0;
                % Consider the wavefunctions extending over
                % 2 times number of wavefunction periods
                % to assure a complete space charge density over the
                % simulation domain.
                for k = -sim_c.num_periods_wf:sim_c.num_periods_wf
                    rhoe2 = interp1(eig.z_wf', rhoe1, ...
                        sim_c.vec_z_poisson-k*sim_c.l_period, 'v5cubic');
                    rhoe2(isnan(rhoe2)) = 0;
                    rhoe = rhoe + rhoe2';
                end
            end
            % Normalize electron density with respect to donor density.
            rhoe = ...
                -rhoe * trapz(sim_c.vec_z_poisson, sim_c.vec_rhod) ...
                / trapz(sim_c.vec_z_poisson, rhoe');
            % Vector with total charge density in cm^(-3).
            rho = sim_c.vec_rhod + rhoe;
        end
    end
end