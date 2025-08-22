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

classdef eigenstates < handle
    % Contains all information about the eigenstates of the quantum system.

    properties %(SetAccess = private)
        psi % matrix: Wavefunctions of the defined QCL [1/m^0.5)].
        E % vector: Eigenenergies of the given eigenstates [eV].
        z_wf % vector: Position vector of wavefunctions [Angstrom].
        num_wfs % scalar: Number of wavefunctions in one QCL-period.
        Escale; % scalar: Scaling factor for plotting.
        hamiltonian; % matrix: System Hamiltonian for 4 QCL-periods [eV].
        dipoles % matrix: Dipole matrix for 4 QCL-periods [Cm].
        m_eff % effective_mass-object: Contains effective mass for each subband [-].
        E_bound_CBO % scalar: Energy with respect to conduction band potential.
        A_test % matrix: Includes all eigenenergies and E_kin values in the tm_solver.
    end

    methods
        function obj = eigenstates(E, psi, z_wf, meff, Ekin)
            % Constructs an object of type eigenstates.
            %
            % Syntax:
            %   obj = eigenstates(E, psi, z_wf, meff)
            %   obj = eigenstates(E, psi, z_wf, meff, Ekin)
            %
            % Input Arguments:
            %   E (vector | matrix): Eigenenergies of the given eigenstates
            %     [eV]. It is also possible to specify E as a matrix with
            %     the eigenenergies on its diagonal.
            %   psi (matrix): Wavefunctions of the given eigenstates.
            %   z_wf (vector): Position vector of wavefunctions.
            %   meff (vector): Effective masses of the eigenstates.
            %   Ekin (vector): Kinetic energy vector.

            % Define hamiltonian and eigenstates by checking the size of
            % the input value E.
            if (isvector(E))
                obj.E = E; % Eigenenergies in eV.
                obj.hamiltonian = diag(E);
            else
                obj.E = diag(E);
                obj.hamiltonian = E;
            end
            obj.psi = psi; % Wavefunctions.
            obj.z_wf = reshape(z_wf, length(z_wf), 1); % Position vector.
            obj.num_wfs = length(obj.E) / 4; % Number of wavefunctions.
            obj.set_Escale(0.05); % Default scaling factor.
            obj.m_eff = meff; % Subband energy resolved effective masses.
            if (nargin < 5)
                obj.E_bound_CBO = 0;
            else
                obj.E_bound_CBO = Ekin;
            end
        end

        function psi = get_wavefct(obj, ind)
            % Returns wavefunction of a specific eigenstate.
            %
            % Syntax:
            %   psi = get_wavefct(obj, ind)
            %
            % Input Arguments:
            %   ind (scalar | vector): Index of the eigenstate.
            %
            % Output Arguments:
            %   psi (vector): Wavefunction.

            psi = obj.psi(:, ind);
        end

        function E_i = get_E_i(obj, ind)
            % Returns eigenenergy of a specific eigenstate.
            %
            % Syntax:
            %   E_i = get_E_i(obj, ind)
            %
            % Input Arguments:
            %   ind (scalar | vector): Index of the eigenstate.
            %
            % Output Arguments:
            %   E_i (vector): Eigenenergy [eV].

            E_i = obj.E(ind);
        end

        function set_Escale(obj, escale)
            obj.Escale = escale;
        end

        function d_ij = get_dipole_element(obj, ind_i, ind_j)
            % Returns dipole matrix element for the given pair of states.
            %
            % Syntax:
            %   d_ij = get_dipole_element(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_i (scalar): Index of first state.
            %   ind_j (scalar): Index of second state.
            %
            % Output Arguments:
            %   d_ij (scalar): Dipole moment [Cm].

            % Wavefunction with index ind_i.
            psi_i = obj.psi(:, ind_i) ...
                / sqrt(trapz(obj.z_wf, obj.psi(:, ind_i).^2));
            % Wavefunction with index ind_j;
            psi_j = obj.psi(:, ind_j) ...
                / sqrt(trapz(obj.z_wf, obj.psi(:, ind_j).^2));
            % Expectation value of the position vector.
            zs = trapz(obj.z_wf, obj.z_wf.*abs(psi_i.*psi_j)) ...
                / trapz(obj.z_wf, abs(psi_i.*psi_j));
            % Correction of integration center for static dipole
            % moments
            % Would be best to have zs in center of each multiplet
            if (ind_i == ind_j)
                zs = mean(obj.z_wf);
            end
            % Dipole matrix element.
            d_ij = phys_const.e0 * trapz(obj.z_wf-zs, psi_i ...
                .*(obj.z_wf - zs).*psi_j) * 1e-10;
        end

        function plot_wavefunctions(obj, cond_band)
            % Plots conduction band profile and magnitude square of all
            % wavefunctions shifted by their respective eigenenergy.
            %
            % Syntax:
            %   plot_wavefunctions(obj, cond_band)
            %
            % Input Arguments:
            %   cond_band (conduction_band-object): Contains information
            %     about the conduction band profile.

            % ColorOrder
            color = custom_colormap.colormap_old(obj.num_wfs);
            % Conductionband profile
            hold on;
            set(gca, 'ColorOrder', color);
            % Probability densities.
            psis2 = obj.get_psi_squared();
            list = cell(obj.num_wfs*4, 1);
            for c = 1:obj.num_wfs * 4
                list{c} = ['Wavefunction', int2str(obj.num_wfs*4-c+1)];
                % conductionband profile + Wavefunctions
                plot(-(obj.z_wf / 10) ...
                    , psis2(:, end-c+1), '-', 'linewidth', 2);
            end
            legend(list);
            cond_band.plot_profile();
        end

        function dipoles = get.dipoles(obj)
            % Returns complete dipole matrix.
            dipoles = zeros(4*obj.num_wfs);
            for ind_i = 1:4 * obj.num_wfs
                for ind_j = 1:4 * obj.num_wfs
                    dipoles(ind_i, ind_j) = ...
                        obj.get_dipole_element(ind_i, ind_j);
                end
            end
        end

        function mark_off_diag = tunnel_marker(obj)
            % Creates matrix specifying which transition corresponds to a
            % tunneling transition. Tunneling transitions are identified
            % by non-zero off-diagonal elements of the hamiltonian matrix.
            %
            % Syntax:
            %   mark_off_diag = tunnel_marker(obj)
            %
            % Output Arguments:
            %   mark_off_diag (matrix): Tunneling marker matrix, where only
            %       the elements corresponding to a tunneling transition
            %       are non-zero and have value 1.

            mark_off_diag = obj.hamiltonian - diag(diag(obj.hamiltonian));
            mark_off_diag = mark_off_diag ~= 0;
        end

        function E_ij = get_anticrossing_energy(obj, ind_i, ind_j)
            % Returns anticrossing energy for selected pair of states.
            %
            % Syntax:
            %   E_ij = get_anticrossing_energy(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_j (scalar): Index of first state.
            %   ind_i (scalar): Index of second state.
            %
            % Output Arguments:
            %   E_ij (scalar): Anticrossing energy [eV].

            E_ij = obj.hamiltonian(ind_i, ind_j);
        end

        function E_ij = get_resonance_energy(obj, ind_i, ind_j)
            % Returns resonance energy (transition energy) for selected
            % pair of states.
            %
            % Syntax:
            %   E_ij = get_resonance_energy(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_j (scalar): Index of first state.
            %   ind_i (scalar): Index of second state.
            %
            % Output Arguments:
            %   E_ij (scalar): Resonance energy [eV].

            E_ij = abs(diag(obj.hamiltonian(ind_i, ind_i))- ...
                diag(obj.hamiltonian(ind_j, ind_j)));
        end

        function omega_ij = get_resonance_freq(obj, ind_i, ind_j)
            % Returns resonance frequency (transition frequency) for
            % selected pair of states.
            %
            % Syntax:
            %   omega_ij = get_resonance_freq(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_j (scalar): Index of first state.
            %   ind_i (scalar): Index of second state.
            %
            % Output Arguments:
            %   omega_ij (scalar): Resonance frequency [rad/s].

            E_ij = obj.get_resonance_energy(ind_i, ind_j); % in eV.
            omega_ij = E_ij * phys_const.e0 / phys_const.hbar;
        end

        function omega_rabi_ij = get_rabi_freq(obj, ind_i, ind_j)
            % Returns Rabi frequency for selected pair of states.
            %
            % Syntax:
            %   omega_rabi_ij = get_rabi_freq(obj, ind_i, ind_j)
            %
            % Input Arguments:
            %   ind_j (scalar): Index of first state.
            %   ind_i (scalar): Index of second state.
            %
            % Output Arguments:
            %   omega_rabi_ij (scalar): Rabi frequency [rad/s].

            omega_rabi_ij = obj.hamiltonian(ind_i, ind_j) ...
                * phys_const.e0 / phys_const.hbar;
        end

        function wt = calc_tunnel_rates(obj)
            % Calculates tunneling rates (Rabi frequency) for all
            % transitions, for which the off-diagonal elements of the
            % hamiltonian is not zero.
            %
            % Syntax:
            %   wt = calc_tunnel_rates(obj)
            %
            % Output Arguments:
            %   wt (matrix): Tunneling rates.

            mark_off_diag = obj.tunnel_marker();
            wt = (mark_off_diag .* obj.hamiltonian) ...
                * phys_const.e0 / phys_const.hbar;
        end

        function wd = calc_resonance_energies(obj, mark_off_diag)
            % Calculates resonance energies for selected transitions.
            %
            % Syntax:
            %   wd = calc_resonance_energies(obj)
            %   wd = calc_resonance_energies(obj, mark_off_diag)
            %
            % Input Arguments:
            %   mark_off_diag (matrix): Boolean matrix specifying for which
            %     transitions in the 4 QCL-periods the resonance energy
            %     should be calculated. If mark_off_diag is not provided as
            %     input, the resonance energies for all transitions,
            %     for which the off-diagonal elements of the hamiltonian is
            %     not zero, are calculated.
            %
            % Output Arguments:
            %   wd (matrix): Resonance energies for the selected states [eV].

            if nargin < 2
                mark_off_diag = obj.tunnel_marker();
            end
            wd = zeros(4*obj.num_wfs);
            for i = 1:4 * obj.num_wfs
                for j = 1:4 * obj.num_wfs
                    if mark_off_diag(i, j) == 1 %i~=j%
                        wd(i, j) = abs(diag(obj.hamiltonian(i, i))- ...
                            diag(obj.hamiltonian(j, j))) * ...
                            phys_const.e0 / phys_const.hbar;
                    end
                end
            end
        end

        function write_wavefcts(obj, ind_wfs, filename)
            % Write wavefunctions of selected states into filename.csv file.
            %
            % Syntax:
            %   write_wavefcts(obj, ind_wfs, filename)
            %
            % Input Arguments:
            %   ind_wfs (scalar | vector): Indices of the states.
            %   filename (char): Specifies absolute or relativ path of
            %     csv-file (without extension).

            psi2_write = obj.get_psi_squared(ind_wfs);
            T_wavefcts = array2table([(-1) * obj.z_wf / 10, psi2_write]);
            header = cell(1, length(ind_wfs));
            header{1} = 'z';
            for i = 2:length(ind_wfs) + 1
                header{i} = append('wf_', num2str(ind_wfs(i-1)));
            end
            T_wavefcts.Properties.VariableNames(1:size(header, 2)) ...
                = header;
            writetable(T_wavefcts, filename);
        end

        function psi_squared = get_psi_squared(obj, ind_wfs)
            % Get probability densities for selected states.
            %
            % Syntax:
            %   psi_squared = get_psi_squared(obj)
            %   psi_squared = get_psi_squared(obj, ind_wfs)
            %
            % Input Arguments:
            %   ind_wfs (scalar | vector): Indices of the states. If
            %      ind_wfs is not provided as input, psi^2 of all states
            %      in the QCL system is returned.
            %
            % Output Arguments:
            %   psi_squared (vector | matrix): Wavefunction squared of the
            %     selected states.

            if (nargin <2)
                % Calculate probability densities.
                for j = 1:4 * obj.num_wfs
                    psi_squared(:, j) = obj.psi(:, j).^2 / ...
                        max(max(obj.psi.^2)) * obj.Escale + obj.E(j);
                end
            else
                for i = 1:length(ind_wfs)
                    psi_squared(:, i) = obj.psi(:, ind_wfs(i)).^2 / ...
                        max(max(obj.psi.^2)) * obj.Escale + ...
                        obj.E(ind_wfs(i));
                end
            end
        end

        function to_hdf5(obj, filename)
            % Saves eigenstates object in hdf5 format.
            %
            % Syntax:
            %   to_hdf5(obj, filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            write_to_hdf5(filename, "/wave_function", obj.psi(end:-1:1, :));
            write_to_hdf5(filename, "/z_position", -obj.z_wf(end:-1:1));
            write_to_hdf5(filename, "/hamiltonian", obj.hamiltonian');
            write_to_hdf5(filename, "/effective_mass", obj.m_eff);
        end
    end

    methods (Static)
        function obj = from_hdf5(filename)
            % Constructs eigenstates object from hdf5 file.
            %
            % Syntax:
            %   obj = from_hdf5(filename)
            %
            % Input Arguments:
            %   filename (string): Name of hdf5 file.

            psi = h5read(filename, "/wave_function");
            z_wf = h5read(filename, "/z_position");
            H = h5read(filename, "/hamiltonian");
            m_eff = h5read(filename, "/effective_mass");
            obj = eigenstates(H', psi(end:-1:1, :), -z_wf(end:-1:1), m_eff);
        end
    end
end
