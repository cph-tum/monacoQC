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

function mb_fujita2017 = ...
    mb_fujita2017(name_mat, num_midIR, num_THz, num_tun, period_mb)
%mb_fujita2017 creates mbsolve python input files for a
% THz DFG QCL test setup.
%
%    mb_fujita2017(name, num_midIR, num_THz, num_tun, period_mb) creates the  
%    input files for the terahertz difference frequency generation QCL comb  
%    setup described in Fujita et al., 2017 
%    (https://doi.org/10.7567/APEX.10.082102).
%    The MAT-file name_mat includes carrier transport simulation results,
%    which are saved in class objects defined in the module results. 
%    The number of optical (mid-IR, THz) and tunneling transitions can be 
%    specified as inputs, and are selected automatically. As an alternative, 
%    one can define the transitions within a cell array and use them as 
%    input arguments specifing the corresponding transitions. 
%    In the input variable period_mb the indices of the 
%    reduced active quantum system are stored.
%
% Load simulation results.
load(name_mat);
% set waveguide
% period length in m.
l_p = 1.04e-6;
% wavelength in m.
lambda_midIR = 6.8e-6;
% effective refractive index.
neff = lambda_midIR / (2 * l_p);
d.n_eff = neff;
r = (neff - 1) / (neff + 1);
w = waveguide(3e-3, 320,[], 0.6, r);
d.set_waveguide(w);

% Investigated QCL is described by an array of wavefunctions.
% period_mb = [29, 23, 28, 21, 24, 22, 20, 19, 15, 18, 13, 14];
mb_name = append(name_mat, '_', num2str(num_midIR), 'midIR', ...
    '_', num2str(num_THz), 'THz', ...
    '_', num2str(num_tun), 'tun');
f_midIR = phys_const.c0 / lambda_midIR; % in Hz
mb_fujita2017 = mbsolve_sim(mb_name, d, period_mb, ...
f_midIR, eigen, carr_dist, deph, sc);

% Find optical mid-IR transition.
if(iscell(num_midIR))
    trans_midIR = num_midIR;
    num_midIR = length(trans_midIR);
else
    trans_midIR = mb_fujita2017.find_optical_trans(f_midIR, num_midIR);
end
% trans_midIR = {[21, 22], [23, 24], [21, 20], [23, 22]};
% num_midIR = size(trans_midIR, 2);

for i = 1:num_midIR
    ni = trans_midIR{i}(1);
    nf = trans_midIR{i}(2);
    mb_fujita2017.add_dipole_pair(ni, nf);
end
% Add dipole pairs taking part in the DFG process.
% Calc chi2
f_THz = 2.92e12; % in Hz
n_3D = 1.6915e22; % in m^{-3}
f_IR_1 = phys_const.c0 / 6.5e-6; % in Hz
f_IR_2 = phys_const.c0 / 6.9e-6; % in Hz
if(iscell(num_THz))
    trans_THz = num_THz;
    num_THz = length(trans_THz);
else
    [chi2, trans_THz, triplet_DFG] = ...
        mb_fujita2017.find_triplet_DFG(num_DFG, f_IR_1, f_IR_2, f_THz, ...
            n_3D);
end
% trans_THz = [23,21; 24,22; 22,20];
for i = 1:num_THz
    ni = trans_THz{i}(1);
    nf = trans_THz{i}(2);
    mb_fujita2017.add_dipole_pair(ni, nf);
end

% Add tunneling transition
if(iscell(num_tun))
    trans_tun = num_tun;
    num_tun = length(trans_tun);
else
    trans_tun = mb_fujita2017.find_tunneling_trans(num_tun);
end

trans_tun = mb_fujita2017.find_tunneling_trans(num_tun);
for i = 1:num_tun
    ni = trans_tun{i}(1);
    nf = trans_tun{i}(2);
    mb_fujita2017.add_tunnel_pair(ni, nf);
end

% Generate python input file.
mb_fujita2017.generate(pwd)

end