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

function mb_fujita2017(filename, outfolder)
% mb_fujita2017 creates the input files for the terahertz
% difference frequency generation QCL comb setup described
% in Fujita et al., 2017 (https://doi.org/10.7567/APEX.10.082102).

% Investigated QCL is described by an array of wavefunctions.
period_mb = [29, 23, 28, 21, 24, 22, 20, 19, 15, 18, 13, 14];

% Load simulation results.
[d, s] = read_input_file('fujita2017');
s.name = 'fujita2017';
s.fmin = 35e12;
s.fmax = 55e12;
[cond_band, eigen, curr_dens, carr_dist, deph, sc] = ...
    load_hdf5_results(filename);
gain = sim_gain(d, eigen, deph, carr_dist);

% Generate mb_fujita2017 class object.
setup = setup_mb(fullfile(outfolder, s.name), s, gain, sc);
setup.set_i_wf(period_mb);
setup.set_A_act(2e-6*10e-6);
setup.set_input_data;

% Find optical mid-IR transition.
f_midIR = phys_const.c0 / 6.8e-6; % in Hz
trans_opt = setup.mb_input_data.find_optical_trans(f_midIR, 5);
% trans_opt = {[21, 22], [23, 24], [21, 20], [23, 22]};
num_opt = size(trans_opt, 2);
for i = 1:num_opt
    ni = trans_opt{i}(1);
    nf = trans_opt{i}(2);
    setup.mb_input_data.add_dipole_pair(ni, nf);
end

% Add dipole pairs taking part in the DFG process.
% Calc chi2
f_THz = 2.92e12; % in Hz
n_3D = 1.6915e22; % in m^{-3}
f_IR_1 = phys_const.c0 / 6.5e-6; % in Hz
f_IR_2 = phys_const.c0 / 6.9e-6; % in Hz
[chi2, trans_THz, triplet_DFG] = setup.mb_input_data. ...
    find_triplet_DFG(5, f_IR_1, f_IR_2, f_THz, n_3D);
% trans_THz = [23,21; 24,22; 22,20];
for i = 1:size(trans_THz, 2)
    ni = trans_THz{i}(1);
    nf = trans_THz{i}(2);
    setup.mb_input_data.add_dipole_pair(ni, nf);
end

% Add tunneling transition.
num_tun = 4;
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
