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

classdef tb_base_transform < handle
    % Base transformation class for tight binding eigenstates.

    methods (Static)
        function V_t = get_V_t(ind_z_L, ind_z_R, V)
            % Calculates the tight-binding potential.
            %
            % Syntax:
            %   V_t = get_V_t(ind_z_L, ind_z_R, V)
            %
            % Input Arguments:
            %   ind_z_L (scalar): Start index of thight-binding period.
            %   ind_z_R (scalar): End index of thight-binding period.
            %   V (vector): Untransformed conduction band profile of
            %     unbiased structure.
            %
            % Output Arguments:
            %   V_t (vector): Tight-binding potential.

            V_t = zeros(size(V));
            V_t(ind_z_L:ind_z_R) = V(ind_z_L:ind_z_R);
            V_t(1:ind_z_L) = V(ind_z_L);
            V_t(ind_z_R:end) = V(ind_z_R);
        end

        function ind_period_tb = find_index_tb_period(ind_z_L, ind_z_R, ...
                z_wf, psi, num_periods)
            % Returns index of the wavefunction period, which corresponds
            % to the tight-binding period.
            %
            % Syntax:
            %   ind_period_tb = find_index_tb_period(ind_z_L, ind_z_R, z_wf, psi, num_periods)
            %
            % Input Arguments:
            %   ind_z_L (scalar): Start index of tight-binding period.
            %   ind_z_R (scalar): End index of tight-binding period.
            %   z_wf (vector): Z-grid vector.
            %   psi (matrix): Matrix containing the wavefunctions of all
            %     states as columns.
            %   num_periods (scalar): Number of periods of the device.
            %
            % Output Arguments:
            %   ind_period_tb (scalar): Index of the tight-binding period.

            % Find position of
            % tight-binding period boundaries (left, right).
            z_ind_L = z_wf(ind_z_L);
            z_ind_R = z_wf(ind_z_R);
            % Find period index of wavefunctions in the tight-binding
            % period. Therefore, we calculate the expectation value of the
            % position vector for wavefunction i.
            num_wfs = size(psi, 2) / num_periods;
            ind_period = zeros(num_wfs, 1);
            for i = 1:size(psi, 2)
                z_wf_i = trapz(z_wf, z_wf'.* ...
                    psi(:, i).^2);
                if (z_wf_i > z_ind_L && z_wf_i < z_ind_R)
                    ind_wf = mod(i - 1, num_wfs) + 1;
                    ind_period(ind_wf) = floor((i - 1)/num_wfs) + 1;
                end
            end

            ind_period_tb = unique(ind_period);
            % Warning, if one wavefunction period does not comprise
            % a tight-binding period.
            if (length(ind_period_tb) ~= 1)
                warning(['The wavefunction period does not', ...
                    ' reside in the tight-binding period!'])
            end

        end

        function ham = calc_hamiltonian(z_wf, psi, E, ind_period_tb, ...
                Vt, Vh, num_periods)
            % Calculates the tigh-binding Hamiltonian.
            %
            % Syntax:
            %   ham = calc_hamiltonian(z_wf, psi, E, ind_period_tb, Vt, Vh, num_periods)
            %
            % Input Arguments:
            %   z_wf (vector): Z-grid vector.
            %   psi (matrix): Matrix containing the wavefunctions of all
            %     states as columns.
            %   E (vector): Eigenenergies of all states.
            %   ind_period_tb (scalar): Index of tight-binding period.
            %   Vt (vector): Tight-binding potential profile [J].
            %   Vh (vector): Normal potential profile [J].
            %   num_periods (scalar): Number of periods of the device.
            %
            % Output Arguments:
            %   ham (matrix): Tight-binding Hamiltonian.

            % Number of considered wavefunction periods.
            m = size(E, 1) / num_periods;

            % Setting up Hamiltonian with eigenenergies
            % on the main diagonal.
            ham = diag(E);

            % Add Rabi energies to the off-diagonal elements
            % of the Hamiltonian.
            for ni = ((ind_period_tb - 1) * m + 1):(ind_period_tb * m)
                for nf = (ind_period_tb * m + 1):((ind_period_tb + 1) * m)

                    % Tunneling to right-neighbouring well.
                    tr = trapz(z_wf, psi(:, ni) ...
                        .*(Vh - Vt)'.*psi(:, nf));
                    % Tunneling from the left-neigbouring well.
                    tl = trapz(z_wf, psi(:, ni-m).* ...
                        (Vh - Vt)'.*psi(:, nf-m));
                    for k = (-3):3
                        if ((ni - k * m > 0) && ...
                                (ni - k * m <= 4 * m) && ...
                                (nf - k * m > 0) && (nf - k * m <= 4 * m))
                            % Calculate mean value of the Rabi energy.
                            ham(ni-k*m, nf-k*m) = sqrt(abs(tl*tr)) ...
                                / phys_const.e0;
                            ham(nf-k*m, ni-k*m) = ham(ni-k*m, nf-k*m);

                        end
                    end
                end
            end
        end
    end
end
