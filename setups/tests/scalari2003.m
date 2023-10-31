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

function [d, s] = scalari2003(basis_sp)
%scalari2003 creates input files for a terahertz test setup.
%
%    scalari2003(folder) creates the input files for the terahertz QCL
%    described in Scalari et al., 2003 (https://doi.org/10.1063/1.1571653).
%    The resulting input files are saved in folder.
%
%    scalari2003(folder, basis_sp) creates the input files and specifies
%    the eigenstates basis basis_sp for the Schrödinger-Poisson solver.
%
%    scalari2003() saves the input files in a folder to be selected via
%    GUI.

%% define scenario and device
% barrier material
b = AlGaAs(0.15);
% well material
w = GaAs();
% Set strain to zero for comparison to fortran files in CI.
b.set_strain('false');
w.set_strain('false');
% set up period
period = {; ...
    layer(b, 24); ...
    layer(w, 110, 2.5e16); ...
    layer(b, 15); ...
    layer(w, 120); ...
    layer(b, 12); ...
    layer(w, 138); ...
    layer(b, 10); ...
    layer(w, 160); ...
    layer(b, 9); ...
    layer(w, 163); ...
    layer(b, 6); ...
    layer(w, 90); ...
    layer(b, 35); ...
    layer(w, 121); ...
    layer(b, 32); ...
    layer(w, 110); ...
    };

% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);

% temperature in K
temperature = 10;

% bias field in kV/cm
bias = 2.55;

% Simulation time in s
t_sim = 7e-11;

% SP simulation basis
if (nargin == 0)
    basis_sp = 'tb';
end
%number of wavefunctions
num_wavefct = 8;

% set up scenario
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);