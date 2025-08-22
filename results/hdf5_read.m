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

function [cond_band, eigen, curr_dens, carr_dist, deph, sc] = hdf5_read(filename)
% Loads simulation results from hdf5 file and generates specific post-
% processing objects.
%
% Syntax:
%   [cond_band, eigen, curr_dens, carr_dist, deph, sc] = hdf5_read(filename)
%
% Input Arguments:
%   filename (string): Name of hdf5 file.
%
% Output Arguments:
%   cond_band (conduction_band-object): Contains conduction band profile.
%   eigen (eigenstates-objct): Contains E, psi and m_eff.
%   curr_dens (current_density-object): Contains current density.
%   carr_dist (carrier_distribution-object): Contains subband occupations.
%   deph (dephasing_rates-object): Contains dephasing rates.
%   sc (scattering_rates-object): Contains scattering rates.

% Eigenstate object.
eigen = eigenstates.from_hdf5(filename);
% Conduction band profile object.
cond_band = conduction_band.from_hdf5(filename);
% Carrier distribution object.
carr_dist = carrier_distribution.from_hdf5(filename);
% Dephasing rates object.
deph = dephasing_rates.from_hdf5(filename, carr_dist);
% Current density object.
curr_dens = current_density.from_hdf5(filename);
% Scattering rates object.
sc = scattering_rates.from_hdf5(filename);
end
