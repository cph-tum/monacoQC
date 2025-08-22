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

classdef (Abstract) ScatteringRate < handle
    % Base class for all scattering mechanisms.

    properties
        Name % string: Name of scattering mechanism.
        tau_inv % matrix: Transistion rates between subbands.
        E % vector: Energy levels of subbands [J].
        mEff % vector: Effective masses of subbands [kg].
        T % scalar: Lattice temperature [K].
        T_e % vector: Electron temperatures of subbands [K].
        psi % vector: Normalized wavefunctions of subbands [1/m^(1/2)].
        z % vector: Grid in z-direction [Angstrom].
        k % vector: Grid in k-space [1/m].
        indices_i_j % vector: Indices of considered subbands.
        num_states % scalar: Number of considered subbands.
    end

    methods
        function obj = ScatteringRate(eigen, device, scenario, options)
            % Constructs an object of type ScatteringRate.
            %
            % Syntax:
            %   obj = AlloyRate(eigen, device, scenario)
            %   obj = AlloyRate(eigen, device, scenario, options)
            %
            % Input Arguments:
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.
            %   device (device-object): Contains information about the
            %     structure/ geometry and materials of the QCL.
            %   scenario (scenario-object): Contains information about the
            %     specific scenario considered for the simulation.
            %   options (cell-array): Array containing name-value pairs for
            %     changing values of default properties. Valid names are
            %     ``num_k`` and ``Te``.

            obj.E = eigen.E * phys_const.e0;
            obj.mEff = eigen.m_eff * phys_const.me;
            obj.T = scenario.T;
            obj.T_e = scenario.T * ones(1, scenario.num_wavefct*4);

            % eigenstates
            obj.z = eigen.z_wf;
            obj.psi = eigen.psi;
            norm = trapz(obj.z*1e-10, obj.psi.^2);
            obj.psi = eigen.psi ./ sqrt(norm); % normalize psi

            % k-space (default value)
            Emax = phys_const.kB * obj.T + device.get_CBO * phys_const.e0;
            obj.set_k(15, Emax);

            % overwrite default values by user defined options
            if nargin > 3
                obj.set_options(options);
            end

            % the scattering rates are computed over the range of
            % two periods: we consider periods 2-3
            p_min = 2;
            p_max = 3;
            start = (p_min - 1) * eigen.num_wfs + 1;
            stop = p_max * eigen.num_wfs;
            obj.indices_i_j = start:stop;
            obj.num_states = stop - start + 1;
            obj.tau_inv = zeros(obj.num_states);
        end

        function set_eigenstates(obj, eigen)
            % Sets new eigenenergies, effective masses and wavefunctions.
            %
            % Syntax:
            %   set_eigenstates(obj, eigen)
            %
            % Input Arguments:
            %   eigen (eigenstates-object): Contains information about
            %     eigenenergies, wavefunctions and effective masses.

            % set new eigenenergies
            obj.E = eigen.E * phys_const.e0;
            % set new effective masses
            obj.mEff = eigen.m_eff * phys_const.me;
            % set new wavefunctions
            obj.psi = eigen.psi;
            norm = trapz(obj.z*1e-10, obj.psi.^2);
            obj.psi = eigen.psi ./ sqrt(norm); % normalize psi
        end

        function set_k(obj, nk, Emax)
            % Sets new k wavevector.
            %
            % Syntax:
            %   set_k(obj, nk, Emax)
            %
            % Input Arguments:
            %   nk (scalar): Number of discretization points.
            %   Emax (scalar): Maximum kinetic energy [J].

            kMax = sqrt(2*max(obj.mEff)*Emax) / phys_const.hbar;
            kInd = 0.5:1:nk;
            obj.k = kInd / max(kInd) * kMax;
        end

        function set_options(obj, options)
            % Changes default values of some properties.
            %
            % Syntax:
            %   set_options(obj, options)
            %
            % Input Arguments:
            %   options (cell-array): Array containing name-value pairs for
            %     changing values of default properties. Valid names are
            %     ``numk`` and ``Te``.

            for i = 1:2:length(options)
                key = options{i};
                value = options{i+1};
                if key == "num_k"
                    kMax = max(obj.k);
                    obj.k = (0.5:1:(value)) / (value - 0.5) .* kMax;
                elseif key == "Te"
                    obj.T_e = value;
                end
            end
        end
    end

    methods (Abstract)
        calculate(obj) % Calculate transition rates.
        calculate_parallel(obj) % Calculate transition rates using
        % parallelization technique.
        update(obj) % Update properties which depend on the occupations.
    end
end
