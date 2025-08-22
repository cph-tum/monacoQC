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

%% Script for running simulations with SP_rate_solver() class
%  Coupled Schrödinger-Poisson-Rate equations are solved self-consistently
%  for simulating static charge carrier transport

% QCL structure
[d, s] = read_input_file("benz2007");
s.V = 9; % applied bias [kV/cm]
s.T = 77; % lattice temperature [K]
s.basis_sp = "tb"; % basis state: "ext", "tb", "ez"

% Initialize Carrier Transport solver
scattering_mech = {"tunneling", "alloy disorder", "impurity", ...
    "acoustic phonon", "lo phonon emission", "lo phonon absorption", ...
    "interface roughness"}; %"electron electron"
SPR_Solver = SP_Rate_solver(d, s, scattering_mech);
SPR_Solver = SPR_Solver.set_solver_options("max_iteration", 10, "relTol", 1e-2);
SPR_Solver = SPR_Solver.set_rate_eq_options("num_k", 40, "relTol", 1e-2, ...
    "num_theta", 20, "parallel", true, "max_iteration", 10);

% Simulate static carrier transport
[cond, eigen, curr, dist, deph, sc, gain, rateSolver, m] = SPR_Solver.solve();
