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

function mb_gaspare2020(filename, outfolder)
% mb_gaspare2020 creates the mbsolve/ mb input files for the
% THz QCL design EV1116 described in Gaspare et al., 2020
% (https://doi.org/10.1515/nanoph-2020-0378)

% load simulation results
[cond_band, eigen, curr_dens, carr_dist, deph, sc] = ...
    load_hdf5_results(filename);
[d, ~] = read_input_file('gaspare2020');

% set up waveguide
l_wg = 2.5e-3; % in m
w_wg = 7e-5; % in m
h_wg = 1e-5; % in m
loss_wg = 857.33; % in 1/m
overlap_f = 1;
refl_coeff_left = 0.8366; % field reflectivity left
refl_coeff_right = 0.8366; % field reflectivity right
wg = laser_waveguide(l_wg, loss_wg, h_wg, overlap_f, ...
    refl_coeff_left, refl_coeff_right, w_wg);
d.set_waveguide(wg);

% set up scenario
bias = 7;
temperature = 25;
s = scenario(temperature, bias, 7e-11, 5, 'tb');
s.name = 'gaspare2020';
s.fmin = 1.5e12;
s.fmax = 5e12;

% Create mb_setup.
gain = sim_gain(d, eigen, deph, carr_dist);
setup = setup_mb(fullfile(outfolder, s.name), s, gain, sc);

% Customize default mbsetup.
setup.set_r_trip(50);
setup.set_shape_cavity('Fabry-Perot');
setup.set_GVD(1e-22);
setup.set_i_wf([7, 8, 9, 10, 11]);
setup.set_input_data;

% Define optical transitions.
num_opt = 1;
trans_opt = setup.mb_input_data.find_optical_trans(setup.f_c, num_opt);
for i = 1:num_opt
    ni = trans_opt{i}(1);
    nf = trans_opt{i}(2);
    setup.mb_input_data.add_dipole_pair(ni, nf);
end

% Define tunneling transitions.
num_tun = 1;
trans_tun = setup.mb_input_data.find_tunneling_trans(num_tun);
for i = 1:num_tun
    ni = trans_tun{i}(1);
    nf = trans_tun{i}(2);
    setup.mb_input_data.add_tunnel_pair(ni, nf);
end

% Generate input files for the full system.
setup.generate_mbsolve();

% Generate input files for two level system.
setup.set_rates_2lvl();
setup.dir = fullfile(outfolder, [char(s.name), '_2lvl']);
setup.generate_mbsolve_2lvl();

end
