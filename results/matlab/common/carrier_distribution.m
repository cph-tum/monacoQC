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

classdef carrier_distribution < handle
    %carrier_distribution Contains the carrier distributions
    % of the specified levels in a QCL system.
    properties (SetAccess = private)
        distribution; % Vector with carrier distributions for one period.
        occupation; % Vector with level occupations for one period.
        E_kin; % Vector of kinetic energy discretization in the subbands.
        T_e; % Vector of electron temperature in the subbands.
        fE; % Function handle of a carrier distribution function for fit.
        nS % Sheet density.
    end
    
    methods
        function obj = carrier_distribution(dist, occ, E, ns)
            % Constructs carrier_distribution.
            obj.distribution = dist; % Carrier distribution
            obj.occupation = occ; % Level occupation
            if nargin > 2
                obj.E_kin = E; % Kinetic energy
                obj.nS = ns; % Sheet density
                % Calculate electron temperatures.
                obj.T_e = ...
                    sum((obj.E_kin * ones(1, length(obj.occupation))) ...
                    .*obj.distribution) * ...
                    phys_const.e0 ./ sum(obj.distribution) / phys_const.kB;
            end
            
            % Set default fitting carrier distribution function.
            obj.set_fit_carr_dist('Maxwell Boltzmann');
        end
        
        function dist_i = get_distribution(obj, ind)
            % Returns carrier distribution for the given state index.
            n_wf = length(obj.T_e);
            % Find corresponding state in the first period.
            ni = mod(ind-1, n_wf) + 1;
            dist_i = obj.distribution(:, ni);
        end
        
        function occ_i = get_occupation(obj, ind)
            % Returns level occupation for the given state index.
            n_wf = length(obj.T_e);
            % Find corresponding state in the first period.
            ni = mod(ind-1, n_wf) + 1;
            occ_i = obj.occupation(ni);
        end
        
        function te_i = get_electron_temperature(obj, ind)
            % Returns electron temperature for the given state index.
            n_wf = length(obj.T_e);
            % Find corresponding state in the first period.
            ni = mod(ind-1, n_wf) + 1;
            te_i = obj.T_e(ni);
        end
        function set_fit_carr_dist(obj, name_fit_dist)
            % Sets carrier distribution function used for fitting.
            name_fit_dist = validatestring(name_fit_dist, ...
                {'Maxwell Boltzmann', 'Fermi Dirac'}, ...
                'set_fit_carr_dist', 'Fit carrier distribution');
            if strcmp(name_fit_dist, 'Maxwell Boltzmann')
                obj.fE = @(x, E) Maxwell_Boltzmann.calc(x(1), x(2), E);
            elseif strcmp(name_fit_dist, 'Fermi Dirac')
                obj.fE = @(x, E) Fermi_Dirac.calc(x(1), x(2), x(3), E);
            end
        end
        
        function [A, T_mb, mu] = fit_distribution(obj, ind_i)
            % Fits carrier distribution of level i.
            % Fits to predefined distribution and returns electron
            % temperature of the fitted curve.
            v = 2:length(obj.E_kin);
            [x, ~] = lsqcurvefit(obj.fE, ...
                [obj.distribution(1, ind_i), obj.T_e(ind_i), 0], ...
                obj.E_kin(v), obj.distribution(v, ind_i));
            A = x(1);
            T_mb = x(2);
            mu = x(3);
        end
        
        function plot_carr_distr(obj)
            % Plots carrier distributions of levels over kinetic energy.
            hold on;
            num_wfs = size(obj.distribution, 2);
            A = custom_colormap.colormap_old(num_wfs);
            set(gca, 'ColorOrder', A);
            plot(obj.E_kin, obj.distribution(:, :), 'linewidth', 2);
            xlabel('E/eV');
            ylabel('f(E)');
            for i = 1:num_wfs
                [A, T_mb, mu] = obj.fit_distribution(i);
                plot(obj.E_kin, ...
                    obj.fE([A, T_mb, mu], obj.E_kin), ':');
            end
            hold off;
        end
        
        function profile_carr_dist(obj, cond_band, eigen, d)
            % Set Nk to get the right dimensions when interpolating.
            Nk = length(obj.distribution(:, 1)) - 1;
            % Kinetic energy.
            Ek = obj.E_kin(1:Nk);
            % Probability density
            psis2 = eigen.get_psi_squared();
            % Fermi distribution.
            dis = obj.distribution(1:Nk, :);
            % Energy.
            E = min(cond_band.Vh):0.0001:max(eigen.E) + max(obj.E_kin);
            % Probability density distribution of all 4 periods.
            FE = 0;
            % Distribution of all states.
            dens = zeros(length(E), length(eigen.E));
            
            % Interpolate the density for every eigenstate (j)
            for j = 1:length(eigen.E)
                % Set window (imin - imax) in which the interpolated
                % distribution (dis) will be inserted.
                [~, imin] = min(abs(E-eigen.E(j)));
                [~, imax] = min(abs(E-eigen.E(j)-max(Ek)));
                Einterp = E(imin:imax) - eigen.E(j);
                % Inserte the interpolation at the positon of the j-th
                % eigenstate.
                dens(imin:imax, j) = interp1(Ek, dis(:, ...
                    mod(j-1, length(obj.distribution(1, :)))+1), ...
                    Einterp, 'spline');
                % Neglect negative interpolations to avoide complex values
                % in the logarithmic mapping.
                dens = dens .* (dens >= 0);
                % Combine the distribution of all states (dens) with the
                % probability density of all states as a vector
                % multiplication and sum them up to get the overall
                % probability density distribution (FE).
                FE = FE + (dens(:, j) * (psis2(:, j) ...
                    -eigen.E(j))');
            end
            
            % Logarithmic mapping.
            rho0 = 1e18;
            n2D = d.dens_sheet;
            c = 4 * n2D / (E(2) - E(1)) / (eigen.z_wf(2) - ...
                eigen.z_wf(1)) * 1e8 / sum(sum(FE));
            x1 = 0.01;
            x2 = 5.0;
            nmax = 64;
            a1 = (nmax - 1) / log(x2/x1);
            a2 = exp(1/a1) / x1;
            image(eigen.z_wf, E, max(1, a1*log(a2*c*FE/rho0)));
            num_color = ceil(max(max(1, a1*log(a2*c*FE/rho0)), [], 'all'));
            colormap(jet(num_color));
            set(gca, 'ydir', 'normal');
            % Conduction band profile.
            hold on;
            plot(cond_band.zv, cond_band.Vh/phys_const.e0, 'Color', ...
                [1, 1, 1], 'linewidth', 2);
            hold off;
            colorbar('ytick', [1, (1:6) / 6 * nmax], 'yticklabel', ...
                [x1, round(1000*exp((1:6)/6*nmax/a1)/a2) / 1000]);
            max(max(c*FE/rho0))
            % Optional:
            % hold on;
            % plot(eigen.z_wf,cond_band.psis2,'Color',[1 1 1],...
            % 'linewidth',2);
        end
    end
end
