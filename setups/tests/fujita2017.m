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

function [d, s] = fujita2017(basis_sp)
% fujita2017 creates input files for a THz DFG test setup.
%
%    fujita2017(folder) creates the strain-balanced
%    THz difference frequency generation QCL setup(wavelength 6.8µm)
%    described in Fujita et al., 2017
%    (https://iopscience.iop.org/article/10.7567/APEX.10.082102).
%    It returns an object of
%    simulation device and scenario. The input argument basis_sp
%    specifies theeigenstates basis for the Schrödinger-Poisson solver.
%

%% define scenario and device
% barrier material
b = InAlAs(0.44);
% well material
w = InGaAs(0.6);

% set up period
period = {; ...
    layer(b, 37); ...
    layer(w, 26); ...
    layer(b, 28); ...
    layer(w, 27); ...
    layer(b, 24, 1.0e17); ...
    layer(w, 28, 1.0e17); ...
    layer(b, 21, 1.0e17); ...
    layer(w, 29, 1.0e17); ...
    layer(b, 18); ...
    layer(w, 30); ...
    layer(b, 16); ...
    layer(w, 32); ...
    layer(b, 15); ...
    layer(w, 36); ...
    layer(b, 12); ...
    layer(w, 45); ...
    layer(b, 11); ...
    layer(w, 49); ...
    layer(b, 9); ...
    layer(w, 60); ...
    layer(b, 26); ...
    layer(w, 24); ...
    };

% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);

% set waveguide
% effective refractive index.
neff = 3.125;
d.n_eff = neff;
r = (neff - 1) / (neff + 1);
w = waveguide(3e-3, 640, [], 0.6, r);
d.set_waveguide(w);

% temperature in K
temperature = 300;

% bias field in kV/cm
% Consolino: V = 13.5V
V = 13.7e-3;
periods = 40;
length = 1e-8 * periods * d.l_period();
bias = V / length;
% bias = 58;
% Simulation time in s
t_sim = 10e-11;

% SP simulation basis
if (nargin == 0)
    basis_sp = 'tb';
end
%number of wavefunctions
num_wavefct = 12;

% set up scenario
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);
