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

function [d, s] = giorgetta2008(basis_sp)
%giorgetta2008 creates input files for a mid-infrared test setup.
%
%    giorgetta2008(folder) creates the strain-balanced
%    mid-infrared QCD setup (wavelength 4 µm, alias: N1037) described in
%    Giorgetta et al., 2008 (https://doi.org/10.1063/1.2902301).
%    It returns an object of
%    simulation device and scenario. The input argument basis_sp
%    specifies theeigenstates basis for the Schrödinger-Poisson solver.
%

%% define scenario and device
% temperature in K
temperature = 100;

% barrier material
b = InAlAs(0.45);
% well material
w = InGaAs(0.61);

% set up period
period = {; ...
    layer(b, 40); ...
    layer(w, 45, 1e18); ...
    layer(b, 60); ...
    layer(w, 10); ...
    layer(b, 50); ...
    layer(w, 13); ...
    layer(b, 40); ...
    layer(w, 16); ...
    layer(b, 35); ...
    layer(w, 19); ...
    layer(b, 30); ...
    layer(w, 22); ...
    layer(b, 30); ...
    layer(w, 27); ...
    layer(b, 26); ...
    layer(w, 33); ...
    };

% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);

% bias field in kV/cm
bias = 0;

% Simulation time in s
t_sim = 1e-10;

% SP simulation basis
if (nargin == 0)
    basis_sp = 'QCD';
end
%number of wavefunctions
num_wavefct = 9;

% set up scenario
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);