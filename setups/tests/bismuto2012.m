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

function [d, s] = bismuto2012(basis_sp)
%bismuto2012 creates input files for a mid-infrared test setup.
%
%    bismuto2012(folder) creates the strain-balanced mid-infrared QCL setup
%    (wavelength 4.3 µm, alias: EV1429) described in
%    Bismuto et al., 2012 (https://doi.org/10.1063/1.4734389), also
%    published in the PhD Thesis of J. Wolf (ETH Zürich, 2017).
%    It returns an object of simulation device and scenario.
%    The input argument basis_sp specifies theeigenstates basis
%    for the Schrödinger-Poisson solver.
%

%% define scenario and device
% temperature in K
temperature = 300;

% barrier material
b = InAlAs(0.335);
% well material
w = InGaAs(0.635);

% set up period
period = {; ...
    layer(b, 35); ...
    layer(w, 17); ...
    layer(b, 24); ...
    layer(w, 17); ...
    layer(b, 19); ...
    layer(w, 19, 0.099e+18); ...
    layer(b, 20, 0.099e+18); ...
    layer(w, 19, 0.099e+18); ...
    layer(b, 22, 0.099e+18); ...
    layer(w, 21, 0.099e+18); ...
    layer(b, 14); ...
    layer(w, 23); ...
    layer(b, 15); ...
    layer(w, 26); ...
    layer(b, 19); ...
    layer(w, 27); ...
    layer(b, 18); ...
    layer(w, 35); ...
    layer(b, 10); ...
    layer(w, 38); ...
    layer(b, 13); ...
    layer(w, 11); ...
    };

% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);

% bias field in kV/cm
bias = 79;

% Simulation time in s
t_sim = 10e-11;

% Schrödinger-Poisson solver basis.
if (nargin == 0)
    basis_sp = 'ext';
end
num_wavefct = 12;

% set up scenario
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);
