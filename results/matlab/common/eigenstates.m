classdef eigenstates < handle
    %eigenstates Contains all information about the eigenstates of the
    %quantum system.
    properties %(SetAccess = private)
        psi; % Wavefunctions of the defined QCL.
        E; % Eigenenergies of the given eigenstates in eV.
        z_wf; % Position vector of wavefunctions.
        num_wfs; % Number of wavefunctions in one period.
        Escale; % Scaling factor.
        hamiltonian; % System Hamiltonian for 4 energy periods in eV.
        dipoles; % Dipole matrix for 4 energy periods.
        m_eff; % Object of of subband energy resolved effective masses.
        E_bound_CBO; % Energy with respect to conduction band potential V.
        A_test; % Includes all Eigenenergies and E_kin values
        % in the tm solver.
    end
    
    methods
        function obj = eigenstates(E, psi, z_wf, meff, Ekin)
            % Constructs eigenstates.
            
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
            % Returns wavefunction of given index.
            psi = obj.psi(:, ind);
        end
        
        function E_i = get_E_i(obj, ind)
            % Returns eigenenergy of given index.
            E_i = obj.E(ind);
        end
        
        function set_Escale(obj, escale)
            obj.Escale = escale;
        end
        
        
        function psi_squared = get_psi_squared(obj)
            % Calculate probability densities.
            for j = 1:4 * obj.num_wfs
                psi_squared(:, j) = obj.psi(:, j).^2 / ...
                    max(max(obj.psi.^2)) * obj.Escale + obj.E(j);
            end
        end
        
        function d_ij = get_dipole_element(obj, ind_i, ind_j)
            % Returns dipole matrix element for the given pair of states.
            % Wavefunction with index ind_i.
            psi_i = obj.psi(:, ind_i) ...
                / sqrt(trapz(obj.z_wf, obj.psi(:, ind_i).^2));
            % Wavefunction with index ind_j;
            psi_j = obj.psi(:, ind_j) ...
                / sqrt(trapz(obj.z_wf, obj.psi(:, ind_j).^2));
            % Expectation value of the position vector.
            zs = trapz(obj.z_wf, obj.z_wf.*abs(psi_i.*psi_j)) ...
                / trapz(obj.z_wf, abs(psi_i.*psi_j));
            %Correction of integration center for static dipole
            %moments
            %Would be best to have zs in center of each multiplet
            if (ind_i == ind_j)
                zs = mean(obj.z_wf);
            end
            % Dipole matrix element.
            d_ij = phys_const.e0 * trapz(obj.z_wf-zs, psi_i ...
                .*(obj.z_wf - zs).*psi_j) * 1e-10;
        end
        
        function plot_wavefunctions(obj, cond_band)
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
            mark_off_diag = obj.hamiltonian - diag(diag(obj.hamiltonian));
            mark_off_diag = mark_off_diag ~= 0;
        end
        
        function E_ij = get_anticrossing_energy(obj, ind_i, ind_j)
            % Returns anticrossing energy for the given pair of states.
            E_ij = obj.hamiltonian(ind_i, ind_j);
        end
        % Returns resonance energy for the given pair of states.
        function E_ij = get_resonance_energy(obj, ind_i, ind_j)
            % E_ij in eV.
            E_ij = abs(diag(obj.hamiltonian(ind_i, ind_i))- ...
                diag(obj.hamiltonian(ind_j, ind_j)));
        end
        % Returns resonance frequency for the given pair of states.
        function omega_ij = get_resonance_freq(obj, ind_i, ind_j)
            % E_ij in eV.
            E_ij = obj.get_resonance_energy(ind_i, ind_j);
            % omega_ij in rad/s.
            omega_ij = E_ij * phys_const.e0 / phys_const.hbar;
        end
        % Returns Rabi frequency for the given pair of states.
        function omega_rabi_ij = get_rabi_freq(obj, ind_i, ind_j)
            % omega_rabi_ij in rad/s.
            omega_rabi_ij = obj.hamiltonian(ind_i, ind_j) ...
                * phys_const.e0 / phys_const.hbar;
        end
        % Calculate tunneling rates.
        function wt = calc_tunnel_rates(obj)
            mark_off_diag = obj.tunnel_marker();
            wt = (mark_off_diag .* obj.hamiltonian) ...
                * phys_const.e0 / phys_const.hbar;
        end
        
        function wd = calc_resonance_energies(obj, mark_off_diag)
            if (~exist("mark_off_diag", 'var'))
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
        
    end
    
    
end