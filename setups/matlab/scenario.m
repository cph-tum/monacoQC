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

classdef scenario < handle
    % Describes specific simulation parameters.

    properties
        name % char: Name of the scenario.
        fmin % scalar: Minimum frequency [Hz].
        fmax % scalar: Maximum frequency [Hz].
        T % scalar: Temperature [K].
        V % scalar: Applied bias field [kV/cm].
        t_sim % scalar: End time of EMC simulation [s].
        t_stat % scalar: EMC time for stationary state [s].
        dt % scalar: Time step of EMC simulation [s].
        basis_sp % string: Type of basis states used in Schrödinger solver.
        num_wavefct = 10 % scalar: Number of wavefunctions per period.
        dz_sp = 1 % 1 (default) | scalar: Grid spacing of the Schrödinger-Poisson solver [Angstrom].
        e_multiplet = 5e-3 % 5e-3 (default) | scalar: Max energy spacing of multiplets for ez-states.
        d_multiplet = 1e-10 % 1e-10 (default) | scalar: Min dipole moment of multiplets for ez-states.
    end

    methods
        function obj = scenario(temperature, bias, t_sim, num_wf, basis_sp)
            % Constructs an object of type scenario.
            %
            % Syntax:
            %   obj = scenario(temperature, bias, t_sim, num_wf, basis_sp)
            %
            % Input Arguments:
            %   temperature (scalar): Temperature [K].
            %   bias (scalar): Applied bias field [kV/cm].
            %   t_sim (scalar): End time of EMC simulation [s].
            %   num_wf (scalar): Number of wavefunctions per period.
            %   basis_sp (scalar): Type of basis states used for solving
            %     the Schrödinger equation. Valid choices are ``tb``
            %     (tight-binding states), ``ext`` (extended states) and
            %     ``ez``.

            obj.T = temperature;
            obj.V = bias;
            obj.t_sim = t_sim;
            obj.t_stat = t_sim;
            dtsim = t_sim * 1e-4;
            if (dtsim < 1e-14)
                dtsim = 1e-14;
            end
            obj.dt = dtsim;
            obj.num_wavefct = num_wf;
            obj.basis_sp = basis_sp;
        end

        function obj = set_grid_sp(obj, dz)
            % Sets grid spacing for the Schrödinger-Poisson solver.
            %
            % Syntax:
            %   set_grid_sp(obj, dz)
            %
            % Input Arguments:
            %   dz (scalar): Grid spacing [Angstrom].

            obj.dz_sp = dz;
        end

        function obj = set_freq_range(obj, fmin, fmax)
            % Sets frequency range for the cavity field.
            %
            % Syntax:
            %   set_freq_range(obj, fmin, fmax)
            %
            % Input Arguments:
            %   fmin (scalar): Minimum frequency [Hz].
            %   fmax (scalar): Maximum frequency [Hz].

            if fmax < fmin
                error('''fmax'' has be greater than ''fmin''!')
            end
            obj.fmin = fmin;
            obj.fmax = fmax;
        end
    end
end
