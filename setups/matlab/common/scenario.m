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

classdef scenario < handle
    %scenario Describes specific simulation parameters.
    %
    properties
        name % Name of the scenario.
        fmin % Frequency range minimum.
        fmax % Frequency range maximum.
        T % Temperature in K.
        V % Bias field in kV/cm.
        t_sim % Simulation end time in s.
        t_stat % time for stationary state in s.
        dt % simulation timestep in s.
        basis_sp % Simulation type (tightbinding or extended).
        num_wavefct = 10 % number of wave functions per period.
        dz_sp = 1; % Grid spacing of the Schrödinger-Poisson solver
        % in Angstrom.
        % only ez-states:
        e_multiplet = 5e-3; % Maximum energy spacing of multiplets.
        d_multiplet = 1e-10; % Minimum dipole moment of multiplets.
    end
    %
    methods
        % Constructs scenario.
        function obj = scenario(temperature, bias, t_sim, num_wf, basis_sp)
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
        % Sets grid spacing for the Schrödinger-Poisson solver.
        function obj = set_grid_sp(obj, dz)
            obj.dz_sp = dz;
        end
    end
end
