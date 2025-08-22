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

classdef poisson_solver < handle
    % Solver class to solve the Poisson equation.

    methods (Static)
        function DV = calc_potential(sim_c, rho)
            % Calculate additional electrostatic potential arising from
            % space charge effects according to Eq. 20-23 in
            % https://doi.org/10.1063/1.4863665.
            %
            % Syntax:
            %   DV = calc_potential(sim_c, rho)
            %
            % Input Arguments:
            %   sim_c (sim_constants-object): Contains simulation constants.
            %   rho (vector): Charge density [C].
            %
            % Output Arguments:
            %   DV (vector): Electrostatic potential.

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

        function rho = get_rho(sim_c, fun_carr_distr, eig)
            % Get space charge density rho for given eigenstates
            % and carrier distribution.
            %
            % Syntax:
            %   rho = get_rho(sim_c, fun_carr_distr, [])
            %   rho = get_rho(sim_c, fun_carr_distr, eig)
            %
            % Input Arguments:
            %   sim_c (sim_constants-object): Contains simulation constants.
            %   fun_carr_distr (carrier_distribution-object): Contains the
            %     subband occupations required to calculate the mobile
            %     charge carrier density. A uniform distribution is assumed
            %     by default.
            %   eig (eigenstates-objec): Contains information about the
            %     wavefunctions, eigenenergies and effective masses.
            %
            % Output Arguments:
            %   rho (vector): Space charge density [1/m^3].

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
