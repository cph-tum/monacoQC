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

%% Script for running simulations with RateSolver class
%  Rate equations are solved self-consistently

%% Set up device and scenario
[d, s] = read_input_file("benz2007");
s.V = 9; % applied bias [kV/cm]
s.T = 77; % lattice temperature [K]
s.basis_sp = "tb"; % basis state: "ext", "tb", "ez"
s.set_freq_range(1e12, 5e12);

%% Solve Schrödinger-Poisson equation
SP_solver = tm_solver(d, s);
[eigen, cond_band] = SP_solver.solve();
cond = conduction_band(cond_band.zv, SP_solver.sim_const.vec_V_0);

%% Set up rate equation solver
% Specify parameters for rate solver as cell-array with unique key-value pairs
rateOptions = {"num_k", 40, "Te", s.T * ones(1, s.num_wavefct*4)};

% Create objects of scattering rate mechanisms
roughnessRate = InterfaceRoughnessRate(eigen, d, s, rateOptions);
impurityRate = ImpurityRate(eigen, d, s, rateOptions);
tunneling = TunnelingRate(eigen, d, s, cond, rateOptions);
loPhononRateAbsorption = LOPhononRate(eigen, d, s, "absorption", rateOptions);
loPhononRateEmission = LOPhononRate(eigen, d, s, "emission", rateOptions);
acousticRate = AcousticPhononRate(eigen, d, s, rateOptions);
alloyRate = AlloyRate(eigen, d, s, rateOptions);
opticsRate = OpticalRate(eigen, d, s, rateOptions);
eeRate = ElectronElectronRate(eigen, d, s, [rateOptions, {"num_theta", 20, "num_interp", 200}]);

rates = { ...
    roughnessRate, ...
    impurityRate, ...
    acousticRate, ...
    loPhononRateEmission, ...
    loPhononRateAbsorption, ...
    alloyRate, ...
    eeRate, ...
    tunneling, ...
    };

%% Solve rate equations
rateSolver = RateSolver(rates);
dist = rateSolver.solve("parallel", true, "relTol", 1e-2, "nmax", 10);
