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

function [d, s] = bismuto2010(basis_sp)
%bismuto2010 creates simulation scenario and device for a mid-IR setup.
%
%    bismuto2010(basis_sp) creates the mid-infrared QCL setup
%    (wavelength 8.5µm, alias: EV2103) in Bismuto et al., 2010
%    (https://doi.org/10.1063/1.3377008). It returns an object of
%    simulation device and scenario. The input argument basis_sp
%    specifies theeigenstates basis for the Schrödinger-Poisson solver.
%

%% define scenario and device
% temperature in K
temperature = 300;

% barrier material
b = InAlAs(0.52);
% well material
w = InGaAs(0.53);

% Set strain to zero for comparison to fortran files in CI.
b.set_strain('false');
w.set_strain('false');
% set up period
period = {; ...
    layer(b, 34, 1.2e17); ...
    layer(w, 31, 1.2e17); ...
    layer(b, 24); ...
    layer(w, 33); ...
    layer(b, 17); ...
    layer(w, 36); ...
    layer(b, 14); ...
    layer(w, 43); ...
    layer(b, 11); ...
    layer(w, 48); ...
    layer(b, 10); ...
    layer(w, 53); ...
    layer(b, 8); ...
    layer(w, 18); ...
    layer(b, 40); ...
    layer(w, 29); ...
    };

% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);

% bias field in kV/cm
bias = 57.0;

% Simulation time in s
t_sim = 10e-11;

% SP simulation basis
if (nargin == 0)
    basis_sp = 'tb';
end
% Number of wavefunctions
num_wavefct = 9;

% set up scenario
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);