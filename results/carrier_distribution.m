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

classdef carrier_distribution < handle
    % Contains the carrier distributions of the specified levels in one
    % QCL-period.

    properties (SetAccess = private)
        distribution % matrix: K-space distribution of occupation probabilities within each subband [-].
        occupation % vector: Occupation probabilities of each subbands [-].
        E_kin % matrix: Discretized kinetic energy vector for each subband [eV].
        T_e % vector: Electron temperature within each subband [K].
        fE % function_handle: Fit function for carrier distribution.
        nS % matrix: K-space distribution of electron density within each subband [1/m^2].
    end

    methods
        function obj = carrier_distribution(dist, occ, E, ns)
            % Constructs an object of type carrier_distribution.
            %
            % Syntax:
            %   obj = carrier_distribution(dist, occ)
            %   obj = carrier_distribution(dist, occ, E, ns)
            %
            % Input Arguments:
            %   dist (matrix): K-space distribution of occupation
            %     probabilities for each subband. The individual
            %     distributions should be specified as columns of dist.
            %   occ (vector): Occupation probabilities of each subband [-].
            %   E (column vector | matrix): Discretized kinetic energy
            %     vector for each subband [eV]. The individual Ekin vectors
            %     should be specified as columns of E.
            %   ns (matrix): K-space distribution of electron densities
            %     within each subband [1/m^2]. The individual distributions
            %     should be specified as the columns of ns.

            obj.distribution = dist; % Carrier distribution
            obj.occupation = occ; % Level occupation
            if nargin > 2
                % Kinetic energy (ensure that each subband has
                % corresponding Ekin vector)
                obj.E_kin = E .* ones(1, length(obj.occupation));
                obj.nS = ns; % Sheet density
                % Calculate electron temperatures.
                obj.T_e = ...
                    sum(obj.E_kin.*obj.distribution) * ...
                    phys_const.e0 ./ sum(obj.distribution) ./ phys_const.kB;
            end

            % Set default fitting carrier distribution function.
            obj.set_fit_carr_dist('Maxwell Boltzmann');
        end

        function dist_i = get_distribution(obj, ind)
            % Returns carrier distribution for a single subband.
            %
            % Syntax:
            %   dist_i = get_distribution(obj, ind)
            %
            % Input Arguments:
            %   ind (scalar): Subband index.
            %
            % Output Arguments:
            %   dist_i (vector): Carrier distribution of specified subband.

            n_wf = length(obj.T_e);
            % Find corresponding state in the first period.
            ni = mod(ind-1, n_wf) + 1;
            dist_i = obj.distribution(:, ni);
        end

        function occ_i = get_occupation(obj, ind)
            % Returns subband occupation for a single subband.
            %
            % Syntax:
            %   occ_i = get_occupation(obj, ind)
            %
            % Input Arguments:
            %   ind (scalar): Subband index.
            %
            % Output Arguments:
            %   occ_i (scalar): Occupation of specified subband.

            n_wf = length(obj.T_e);
            % Find corresponding state in the first period.
            ni = mod(ind-1, n_wf) + 1;
            occ_i = obj.occupation(ni);
        end

        function te_i = get_electron_temperature(obj, ind)
            % Returns electron temperature for a single subband.
            %
            % Syntax:
            %   te_i = get_electron_temperature(obj, ind)
            %
            % Input Arguments:
            %   ind (scalar): Subband index.
            %
            % Output Arguments:
            %   te_i (scalar): Electron temperature of specified subband [K].

            n_wf = length(obj.T_e);
            % Find corresponding state in the first period.
            ni = mod(ind-1, n_wf) + 1;
            te_i = obj.T_e(ni);
        end

        function set_fit_carr_dist(obj, name_fit_dist)
            % Sets carrier distribution fit function.
            %
            % Syntax:
            %   set_fit_carr_dist(obj, name_fit_dist)
            %
            % Input Arguments:
            %   name_fit_dist (char): Name of the fit function. Valid
            %     choices are ``Maxwell Boltzmann`` and ``Fermi Dirac``.

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
            % Fits carrier distribution of a single subband according to
            % the pre-defined fit-function.
            %
            % Syntax:
            %   [A, T_mb, mu] = fit_distribution(obj, ind_i)
            %
            % Input Arguments:
            %   ind_i (scalar): Subband index.
            %
            % Output Arguments:
            %   A (scalar): Fitted amplitude.
            %   T_mb (scalar): Fitted electron temperature [K].
            %   mu (scalar): Fitted chemical potential [eV].

            % Fits to predefined distribution and returns electron
            % temperature of the fitted curve.
            v = 2:length(obj.E_kin);
            [x, ~] = lsqcurvefit(obj.fE, ...
                [obj.distribution(1, ind_i), obj.T_e(ind_i), 0], ...
                obj.E_kin(v, ind_i), obj.distribution(v, ind_i));
            A = x(1);
            T_mb = x(2);
            mu = x(3);
        end

        function plot_carr_distr(obj, plot_fit)
            % Plots carrier distributions of subbands over kinetic energy.
            %
            % Syntax:
            %   plot_carr_distr(obj)
            %   plot_carr_distr(obj, plot_fit)
            %
            % Input Arguments:
            %   plot_fit (logical): Adds fit function to the plot.

            if nargin < 2
                plot_fit = 1;
            end

            figure;
            num_wfs = size(obj.distribution, 2);
            A = custom_colormap.colormap_old(num_wfs);
            set(gca, 'ColorOrder', A);

            % Plot distribution function
            plot(obj.E_kin, obj.distribution, 'linewidth', 2);
            hold on;
            xlabel('E (eV)');
            ylabel('f(E)');
            if plot_fit
                for i = 1:num_wfs
                    [A, T_mb, mu] = obj.fit_distribution(i);
                    plot(obj.E_kin, ...
                        obj.fE([A, T_mb, mu], obj.E_kin), ':');
                end
            end
            hold off;
        end

        function profile_carr_dist(obj, cond_band, eigen, d)
            % Plots carrier distribution over position and energy with
            % conduction band profile.
            %
            % Syntax:
            %   profile_carr_dist(obj, cond_band, eigen, d)
            %
            % Input Arguments:
            %   cond_band (conduction_band-object): Contains information
            %     about the conduction band profile.
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.
            %   d (device-object): Contains information about the
            %     structure/ geometry and materials of the QCL.

            % Set Nk to get the right dimensions when interpolating.
            Nk = length(obj.distribution(:, 1)) - 1;
            % Kinetic energy.
            Ek = obj.E_kin(1:Nk, :);
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
                [~, imax] = min(abs(E-eigen.E(j)-max(Ek, [], "all")));
                Einterp = E(imin:imax) - eigen.E(j);
                % Inserte the interpolation at the positon of the j-th
                % eigenstate.
                index = mod(j-1, length(obj.distribution(1, :))) + 1;
                dens(imin:imax, j) = interp1(Ek(:, index), ...
                    dis(:, index), Einterp, 'spline');
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

        function to_hdf5(obj, filename)
            % Saves carrier_distribution object in hdf5 format.
            %
            % Syntax:
            %   to_hdf5(obj, filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            write_to_hdf5(filename, "/kinetic_energy", obj.E_kin);
            write_to_hdf5(filename, "/carrier_distribution", obj.distribution);
            write_to_hdf5(filename, "/occupation", obj.occupation);
            write_to_hdf5(filename, "/carrier_sheet_density", obj.nS);
        end
    end

    methods (Static)
        function obj = from_hdf5(filename)
            % Constructs carrier_distribution object from hdf5 file.
            %
            % Syntax:
            %   obj = from_hdf5(filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            dist = h5read(filename, "/carrier_distribution");
            occ = h5read(filename, "/occupation");
            Ekin = h5read(filename, "/kinetic_energy");
            nS = h5read(filename, "/carrier_sheet_density");
            obj = carrier_distribution(dist, occ', Ekin, nS);
        end

        function obj = generate_uniform(num_wfs)
            % Constructs uniform carrier_distribution object.
            %
            % Syntax:
            %   obj = generate_uniform(num_wfs)
            %
            % Input Arguments:
            %   num_wfs (scalar): Number of subbands.

            num_E_kin = 150; % Set kinetic energy vector E_kin.
            dist = ones(num_E_kin, num_wfs) / num_E_kin / num_wfs;
            occ = ones(1, num_wfs) / num_wfs;

            obj = carrier_distribution(dist, occ);
        end

        function obj = generate_thermal(T, n2D, eig_system)
            % Constructs thermailized carrier_distribution object.
            %
            % Syntax:
            %   obj = generate_thermal(T, n2D, eig_system)
            %
            % Input Arguments:
            %   T (vector): Electron temperatures for each subband [K].
            %   n2D (vector): Sheet densities of each subband [1/m^2].
            %   eig_system (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.

            E = eig_system.E(1:eig_system.num_wfs)';
            m_star = phys_const.me * eig_system.m_eff(1:eig_system.num_wfs)';

            % Calculate the thermal distribution of levels.
            Ea = min(E);

            % Function handle sheet density.
            nsi = @(E, mu, T) ...
                (m_star * phys_const.kB * T) / (pi * phys_const.hbar^2) .* ...
                (log(1+exp(phys_const.e0*(mu - E)/(phys_const.kB * T))));
            ns = @(E, mu, T) sum(nsi(E, mu, T));
            f = @(mu) (ns(E, mu, T) / n2D - 1);

            % Setting up options for fzero method.
            opt = optimset('TolX', thermal_carrier_distribution.dE);

            % Use fzero to find chemical potential.
            mu = fzero(f, Ea, opt);

            % Calculate occupations from chemical potential.
            occ = nsi(E, mu, T)' / n2D;

            % Set kinetic energy vector E_kin.
            num_E_kin = 150;
            E_kin_max = 0.0012 * T;

            E_kin = (1:num_E_kin)' * E_kin_max / num_E_kin;

            % Calculate carrier distribution.
            dist = Fermi_Dirac.calc(1, T, mu, E+E_kin);
            dist = dist / sum(sum(dist)) * n2D;
            % Calculate sheet density.
            nS = nsi(E+E_kin, mu, T);

            obj = carrier_distribution(dist, occ, E_kin, nS);
        end
    end
end
