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

function [d, s] = ...
    forrer2021(basis_sp, bias)
%forrer2021 creates simulation scenario and device for a THz QCL test setup.
%
%    forrer2021(basis_sp, bias) generates device and scenario for the
%    THz QCL simulated in Forrer et al., 2021
%    (https://doi.org/10.1063/5.0041339). Here, basis_sp specifies the
%    wavefunction basis for the SP solver and bias gives the applied bias.
%
%    forrer2021() generates device and scenario with tight-binding basis
%    and for a bias of 49 mV per period.

%% define scenario and device
% barrier material
b = AlGaAs(0.15);
% well material
w = GaAs();
% set up period
period = {; ...
    layer(b, 58); ...
    layer(w, 183, 2.3e16); ...
    layer(b, 41); ...
    layer(w, 92); ...
    layer(b, 39); ...
    layer(w, 115); ...
    layer(b, 29); ...
    layer(w, 106); ...
    };
% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);
% set waveguide
% effective refractive index.
neff = 3.6;
d.n_eff = neff;
r = (neff - 1) / (neff + 1);
w = waveguide(5e-3, 2200, [], 1, r);
d.set_waveguide(w);

% temperature in K
temperature = 80;

% bias field in kV/cm
if (nargin < 2)
    V = 49e-6;
else
    V = bias * 1e-6;
end
l_period = d.l_period * 1e-8;
bias = V / l_period;

% Simulation time in s
t_sim = 5e-11;

% SP simulation basis
if (nargin <= 1)
    basis_sp = 'tb';
end

% Number of wavefunctions
num_wavefct = 5;

% set up scenario
% s = scenario(temperature, bias, t_sim);
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);
s.e_multiplet = 5e-3;
end
