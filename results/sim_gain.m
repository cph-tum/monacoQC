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

classdef sim_gain
    % Calculates and plots the optical gain inside the cavity.

    properties
        eigen; % eigenstates-object: Object of class eigenstates.
        device; % device-object: Object of class device.
        dep; % dephasing_rates-object: Object of class dephasing_rates.
        c_dist; % carrier_distribution-object: Object of class carrier_distribution.
    end

    methods
        function obj = sim_gain(d, eigen, deph, carr_dist)
            % Constructs an object of type sim_gain.
            %
            % Syntax:
            %   obj = sim_gain(d, eigen, deph, carr_dist)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.
            %   deph (dephasing_rates-object): Contains information about
            %     pure dephasing rates and lifetime broadening values.
            %   carr_dist (carrier_distribution-object): Contains
            %     information about the carrier distributions of all
            %     subbands.

            obj.eigen = eigen;
            obj.device = d;
            obj.dep = deph;
            obj.c_dist = carr_dist;
        end

        function gain = calc_gain(obj, fmin, fmax, ind_i, ind_j, id)
            % Calculates the gain spectrum including spectral broadening
            % (lifetime broadening and pure dephasing) with lorentzian line
            % shape function. In the gain spectrum, all possible
            % transitions in k-space for one QCL-period (i.e. transitions
            % within the period and to the next nearest period) can be
            % considered. If no specific transition is specified, the total
            % gain, including all possible transitions, is calculated.
            %
            % Syntax:
            %   gain = calc_gain(obj, fmin, fmax, ind_i, ind_j, id)
            %   gain = calc_gain(obj, fmin, fmax, [], [], id)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of gain spectrum.
            %   fmax (scalar): Maximum frequency of gain spectrum.
            %   ind_i (scalar): Initial state index.
            %   ind_j (scalar): Final state index.
            %   id (char): Broadening mechanism, including either both
            %     lifetime broadening and pure dephasing rates
            %     (``id=total``) or only lifetime broadening
            %     (``id=ltbroad``).
            %
            % Output Arguments:
            %   gain (vector): Gain spectrum [1/m].

            % Initial.
            L_p = obj.device.l_period;
            n_eff = obj.device.n_eff;
            omega = 2 * pi * (fmin + (fmax - fmin) * (1:100) / 100);
            g = zeros(size(omega));
            nS = obj.c_dist.nS();
            if ~isempty(ind_i) && ~isempty(ind_j)
                % Iterate over every possible transition.
                x = 1;
                y = obj.eigen.num_wfs * 4;
            else
                % Iterate over every transition of the mid period.
                x = obj.eigen.num_wfs * 2 + 1;
                y = obj.eigen.num_wfs * 3;
            end
            for i = x:y
                for j = 1:obj.eigen.num_wfs * 4
                    % Dipole elements between i -> j, j -> i (matrix ixj).
                    d_ij(i, j) = obj.eigen.get_dipole_element(i, j);
                    % Initial gain of transition i -> j over frequency
                    % vector omega.
                    g_ij = zeros(size(omega));
                    % Resolve in k-space.
                    for k = 1:length(obj.c_dist.E_kin)
                        % Set option for dephasing rate:
                        % - only lifetime broadening {'ltbroad'}
                        % - lifetime broadening + pure dephasing {'total'}
                        switch id
                            case 'ltbroad'
                                gamma_ijk = ...
                                    obj.dep.get_lt_broadening(i, j, k);
                            case 'total'
                                gamma_ijk = ...
                                    obj.dep.get_dephasing_rate ...
                                    (i, j, k);
                            otherwise
                                error(['Error: Choose one of ...', ...
                                    'the following ', 'options ...', ...
                                    '{''total'', ''ltbroad''}.'])
                        end
                        % Gain of transition i -> j over frequency
                        % vector omega in k.
                        g_ij = g_ij ...
                            +sign(obj.get_transition_freq(i, j, k)) ...
                            * omega / phys_const.c0 / n_eff * d_ij(i, j)^2 ...
                            / phys_const.eps0 / phys_const.hbar ...
                            * (nS(k ...
                            , mod(i-1, obj.eigen.num_wfs)+1) ...
                            / L_p / 1e-10) ./ (1 + (omega ...
                            -abs(obj.get_transition_freq(i, j, k))).^2 ...
                            / gamma_ijk^2) / gamma_ijk;
                    end
                    if ~isempty(ind_i) && ~isempty(ind_j) && ...
                            i == ind_i && ind_j == j
                        % Return gain of transition i -> j.
                        gain = g_ij;
                        return
                    else
                        % Return overall gain.
                        g = g + g_ij;
                    end
                end
            end
            gain = g;
        end

        function gain = calc_gain_avg(obj, fmin, fmax, ind_i, ind_j, id)
            % Calculates the gain using k-averaged dephasing rates.
            %
            % Syntax:
            %   gain = calc_gain_avg(obj, fmin, fmax, ind_i, ind_j, id)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of gain spectrum.
            %   fmax (scalar): Maximum frequency of gain spectrum.
            %   ind_i (scalar): Initial state index.
            %   ind_j (scalar): Final state index.
            %   id (char): Broadening mechanism including either both
            %     lifetime broadening and pure dephasing rates
            %     (``id=total``) or only lifetime broadening
            %     (``id=ltbroad``).
            %
            % Output Arguments:
            %   gain (vector): Gain spectrum [1/m].

            L_p = obj.device.l_period * 1e-10;
            n_eff = obj.device.n_eff;
            omega = 2 * pi * (fmin + (fmax - fmin) * (1:100) / 100);
            g = zeros(size(omega));
            nS = obj.c_dist.get_occupation(1:4*obj.eigen.num_wfs) * obj.device.dens_sheet;

            x = obj.eigen.num_wfs * 2 + 1;
            y = obj.eigen.num_wfs * 3;

            for i = x:y
                for j = 1:4 * obj.eigen.num_wfs
                    % dipole moments
                    d_ij = obj.eigen.get_dipole_element(i, j);
                    % transition rates
                    switch id
                        case 'ltbroad'
                            gamma_ij = ...
                                obj.dep.get_lt_broadening(i, j);
                        case 'total'
                            gamma_ij = ...
                                obj.dep.get_dephasing_rate ...
                                (i, j);
                        otherwise
                            error(['Error: Choose one of ...', ...
                                'the following ', 'options ...', ...
                                '{''total'', ''ltbroad''}.'])
                    end
                    % lorentzian lineshape
                    omega_ij = (obj.eigen.E(i) - obj.eigen.E(j)) ./ phys_const.hbar * phys_const.e0;
                    lorentz = 1 ./ (1 + (omega - abs(omega_ij)).^2 ./ gamma_ij.^2) ./ gamma_ij;
                    % gain
                    g_ij = sign(omega_ij) * omega ./ (phys_const.hbar * phys_const.eps0 ...
                        * phys_const.c0 * n_eff * L_p) .* nS(i) .* d_ij.^2 .* lorentz;
                    g = g + g_ij;
                end
            end
            gain = g;
        end

        function w_ijk = get_transition_freq(obj, ind_i, ind_j, ind_k)
            % Returns the transition frequencies from state ik to state jk.
            %
            % Syntax:
            %   w_ijk = get_transition_freq(obj, ind_i, ind_j, ind_k)
            %
            % Input Arguments:
            %   ind_i (scalar): Initial state index.
            %   ind_j (scalar): Final state index.
            %   ind_k (scalar): Index of wavenumber k.
            %
            % Output Arguments:
            %   w_ijk (scalar): Transition frequency [Hz].

            w_ij = (obj.eigen.E(ind_i) - obj.eigen.E(ind_j)) * ...
                phys_const.e0() / phys_const.hbar;
            % Transition frequencies i -> j | j -> i in k.
            w_ijk = w_ij + phys_const.e0 ...
                * obj.c_dist.E_kin(ind_k, 1) ...
                * (1 - obj.eigen.m_eff ...
                (mod(ind_i-1, obj.eigen.num_wfs) + 1) ...
                / obj.eigen.m_eff ...
                (mod(ind_j-1, obj.eigen.num_wfs) + 1)) ...
                / phys_const.hbar;
        end

        function gain_ij = get_gain_ij(obj, fmin, fmax, ind_i, ind_j, id, avg)
            % Returns the gain of the transition from subband i to j.
            %
            % Syntax:
            %   gain_ij = get_gain_ij(obj, fmin, fmax, ind_i, ind_j, id, avg)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of gain spectrum.
            %   fmax (scalar): Maximum frequency of gain spectrum.
            %   ind_i (scalar): Initial state index.
            %   ind_j (scalar): Final state index.
            %   id (char): Broadening mechanism including either both
            %     lifetime broadening and pure dephasing rates
            %     (``id=total``) or only lifetime broadening
            %     (``id=ltbroad``).
            %   avg (logical): Use either k-dependent (``avg=1``) or
            %     k-averaged (``avg=0``) dephasing rates. Default is
            %     ``avg=0``.
            %
            % Output Arguments:
            %   gain_ij (vector): Gain spectrum [1/m].

            if avg
                gain_ij = obj.calc_gain_avg(fmin, fmax, ind_i, ind_j, id);
            else
                gain_ij = obj.calc_gain(fmin, fmax, ind_i, ind_j, id);
            end
        end

        function gain = get_gain(obj, fmin, fmax, id, avg)
            % Returns the total gain, considering all transitions in one
            % period.
            %
            % Syntax:
            %   gain = get_gain(obj, fmin, fmax, id)
            %   gain = get_gain(obj, fmin, fmax, id, avg)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of gain spectrum.
            %   fmax (scalar): Maximum frequency of gain spectrum.
            %   id (char): Broadening mechanism including either both
            %     lifetime broadening and pure dephasing rates
            %     (``id=total``) or only lifetime broadening
            %     (``id=ltbroad``).
            %   avg (logical): Use either k-dependent (``avg=1``) or
            %     k-averaged (``avg=0``) dephasing rates. Default is
            %     ``avg=0``.
            %
            % Output Arguments:
            %   gain (vector): Gain spectrum [1/m].

            if nargin < 5
                avg = 0;
            end
            if avg
                gain = obj.calc_gain_avg(fmin, fmax, [], [], id);
            else
                gain = obj.calc_gain(fmin, fmax, [], [], id);
            end
        end

        function plot_gain(obj, fmin, fmax, id)
            % Plots the total gain, considering all transitions in one
            % period.
            %
            % Syntax:
            %   plot_gain(obj, fmin, fmax, id)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of gain spectrum.
            %   fmax (scalar): Maximum frequency of gain spectrum.
            %   id (char): Broadening mechanism including either both
            %     lifetime broadening and pure dephasing rates
            %     (``id=total``) or only lifetime broadening
            %     (``id=ltbroad``).

            f = (fmin + (fmax - fmin) * (1:100) / 100);
            if nargin < 4
                id = "total";
            end
            g = obj.get_gain(fmin, fmax, id, 0);

            hold on;
            plot(f/1e12, g/100);
            title(['Overall gain (''', id, ''')'], 'FontSize', 17)
            x = xlabel('Frequency [THz]');
            x.FontSize = 17;
            x.Position(2) = x.Position(2) - 0.001;
            y = ylabel('Gain [1/cm]');
            y.FontSize = 17;
            y.Position(1) = y.Position(1) - 0.001;
            hold off;
            [maxg, imax] = max(g);
            disp(['f=', num2str(f(imax)/1e12), ' THz, FWHM=', ...
                num2str(full_width_half_max.calc(g)* ...
                (f(2) - f(1))/1e12), ' THz, max=', ...
                num2str(maxg/100), ' 1/cm']);
        end

        function plot_gain_avg(obj, fmin, fmax, id)
            % Plots the overall gain with k-averaged rates.
            %
            % Syntax:
            %   plot_gain_avg(obj, fmin, fmax)
            %   plot_gain_avg(obj, fmin, fmax, id)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of gain spectrum.
            %   fmax (scalar): Maximum frequency of gain spectrum.
            %   id (char): Broadening mechanism including either both
            %     lifetime broadening and pure dephasing rates
            %     (``id=total``) or only lifetime broadening
            %     (``id=ltbroad``).

            f = (fmin + (fmax - fmin) * (1:100) / 100);
            if nargin < 4
                id = "total";
            end
            g = obj.get_gain(fmin, fmax, id, 1);

            hold on;
            plot(f/1e12, g/100);
            title(['Overall gain (''', id, ''')'], 'FontSize', 17)
            x = xlabel('Frequency [THz]');
            x.FontSize = 17;
            x.Position(2) = x.Position(2) - 0.001;
            y = ylabel('Gain [1/cm]');
            y.FontSize = 17;
            y.Position(1) = y.Position(1) - 0.001;
            hold off;
            [maxg, imax] = max(g);
            disp(['f=', num2str(f(imax)/1e12), ' THz, FWHM=', ...
                num2str(full_width_half_max.calc(g)* ...
                (f(2) - f(1))/1e12), ' THz, max=', ...
                num2str(maxg/100), ' 1/cm']);
        end

        function plot_gain_ij(obj, fmin, fmax, ind_i, ind_j, id)
            % Plots gain of transition between subband i to j.
            %
            % Syntax:
            %   plot_gain_ij(obj, fmin, fmax, ind_i, ind_j, id)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency of gain spectrum.
            %   fmax (scalar): Maximum frequency of gain spectrum.
            %   ind_i (scalar): Initial state index.
            %   ind_j (scalar): Final state index.
            %   id (char): Broadening mechanism including either both
            %     lifetime broadening and pure dephasing rates
            %     (``id=total``) or only lifetime broadening
            %     (``id=ltbroad``).

            f = (fmin + (fmax - fmin) * (1:100) / 100);
            if nargin > 5
                g = obj.get_gain_ij(fmin, fmax, ind_i, ind_j, id);
            else
                id = 'total';
                g = obj.get_gain_ij(fmin, fmax, ind_i, ind_j, id);
            end
            hold on;
            plot(f/1e12, g/100);
            title(['Gain of transition ', num2str(ind_i), ...
                ' \rightarrow ', num2str(ind_j), ' (''', id, ''')'] ...
                , 'FontSize', 17)
            x = xlabel('Frequency [THz]');
            x.FontSize = 17;
            x.Position(2) = x.Position(2) - 0.001;
            y = ylabel('Gain [1/cm]');
            y.FontSize = 17;
            y.Position(1) = y.Position(1) - 0.001;
            hold off;
            [maxg, imax] = max(g);
            disp(['f=', num2str(f(imax)/1e12), ' THz, FWHM=', ...
                num2str(full_width_half_max.calc(g)* ...
                (f(2) - f(1))/1e12), ' THz, max=', ...
                num2str(maxg/100), ' 1/cm']);
        end
    end

    methods (Static)
        function [gain, f] = calculate_gain(d, dipoles, E, ...
                occ, gamma, fmin, fmax, num_points)
            % Static method for calculating the optical gain spectrum.
            %
            % Syntax:
            %   [gain, f] = calculate_gain(d, dipoles, E, occ, gamma, fmin, fmax)
            %   [gain, f] = calculate_gain(d, dipoles, E, occ, gamma, fmin, fmax, num_points)
            %
            % Input Arguments:
            %   d (device-object): Contains information about the
            %     structure, geometry and materials of the QCL.
            %   dipole (matrix): Dipole matrix of 4 QCL periods [Cm].
            %   E (vector): Energies of eigenstates [eV].
            %   occ (vector | matrix): Occupation probabilities (k-averaged
            %     or k-resolved).
            %   gamma (matrix | 3-d array): Dephasing rates (k-averaged or
            %     k-resolved).
            %   fmin (scalar): Minimum frequency of the gain spectrum.
            %   fmax (scalar): Maximum frequency of the gain spectrum.
            %   num_points (scalar): Number of points for the gain
            %     spectrum.
            %
            % Output Arguments:
            %   gain (vector): Gain spectrum [1/m].
            %   f (vector): Frequency vector of gain spectrum.

            if nargin < 8
                num_points = 100;
            end

            % If k-resolved occupations and dephasing rates are given,
            % calculate k-resolved optical gain.
            if length(size(gamma)) == 3
                num_k = size(gamma, 3);
                f = (fmin + (fmax - fmin) * (1:num_points) / num_points);
                gain = zeros(size(f));
                for i = 1:num_k
                    [g_k, ~] = sim_gain.calculate_gain(d, dipoles, ...
                        E, occ(i, :), gamma(:, :, i), ...
                        fmin, fmax, num_points);
                    gain = gain + g_k;
                end
                return
            end

            L_p = d.l_period * 1e-10;
            n_eff = d.n_eff;
            omega = 2 * pi * (fmin + (fmax - fmin) * (1:num_points) / num_points);
            g = zeros(size(omega));
            nS = occ * d.dens_sheet;
            num_wfs = size(dipoles, 1) / 4;

            x = num_wfs * 2 + 1;
            y = num_wfs * 3;

            for i = x:y
                for j = 1:4 * num_wfs
                    % dipole moments
                    d_ij = dipoles(i, j);
                    gamma_ij = gamma(i, j);
                    % lorentzian lineshape
                    omega_ij = (E(i) - E(j)) ./ ...
                        phys_const.hbar * phys_const.e0;
                    lorentz = 1 ./ (1 + (omega - abs(omega_ij)).^2 ./ ...
                        gamma_ij.^2) ./ gamma_ij;
                    % gain
                    g_ij = sign(omega_ij) * omega ./ ...
                        (phys_const.hbar * phys_const.eps0 ...
                        * phys_const.c0 * n_eff * L_p) .* ...
                        nS(mod(i-1, num_wfs)+1) .* d_ij.^2 .* lorentz;
                    g = g + g_ij;
                end
            end
            gain = g;
            f = omega / 2 / pi;
        end
    end
end
