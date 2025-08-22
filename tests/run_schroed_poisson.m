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

% QCL structure
[d, s] = read_input_file("scalari2010");
s.V = 12; % applied bias [kV/cm]
s.T = 10; % lattice temperature [K]
s.basis_sp = "tb"; % basis state: "ext", "tb", "ez"

% Run Schrödinger-Poisson solver
SP_solver = tm_solver(d, s);
[eigen, cond_band] = SP_solver.solve();
