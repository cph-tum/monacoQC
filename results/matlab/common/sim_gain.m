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

classdef sim_gain
    % sim_gain generates an object of class emc_sim gain to plot the
    % overall gain in the cavity or of just one transition depending on the
    % chosen dephasing mechanism.
    properties
        eigen; % Object of class eigenstates.
        device; % Object of class device.
        dep; % Object of class dephasing_rates.
        c_dist; % Object of class carrier_distribution.
    end
    methods
        % Constructs sim_gain.
        function obj = sim_gain(d, eigen, deph, carr_dist)
            obj.eigen = eigen;
            obj.device = d;
            obj.dep = deph;
            obj.c_dist = carr_dist;
        end
        
        % Return transition frequencies i -> j | j -> i (matrix ixj).
        function w_ijk = get_transition_freq(obj, ind_i, ind_j, ind_k)
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
        % Return the gain of transition i -> j.
        function gain_ij = get_gain_ij(obj, fmin, fmax, ind_i, ind_j, id)
            gain_ij = obj.calc_gain(fmin, fmax, ind_i, ind_j, id);
        end
        % Return the overall gain.
        function gain = get_gain(obj, fmin, fmax, id)
            gain = obj.calc_gain(fmin, fmax, [], [], id);
        end
        % Plot the overall gain.
        function plot_gain(obj, fmin, fmax, id)
            f = (fmin + (fmax - fmin) * (1:1000) / 1000);
            if nargin > 3
                g = obj.get_gain(fmin, fmax, id);
            else
                id = 'total';
                g = obj.get_gain(fmin, fmax, id);
            end
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
        % Plot gain of transition i -> j.
        function plot_gain_ij(obj, fmin, fmax, ind_i, ind_j, id)
            f = (fmin + (fmax - fmin) * (1:1000) / 1000);
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
    
    methods (Access = private)
        % Calculate the gain.
        function gain = calc_gain(obj, fmin, fmax, ind_i, ind_j, id)
            % Initial.
            L_p = obj.device.l_period;
            n_eff = obj.device.n_eff;
            omega = 2 * pi * (fmin + (fmax - fmin) * (1:1000) / 1000);
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
                    % Dipole elements between i -> j < j -> i (matrix ixj).
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
    end
end