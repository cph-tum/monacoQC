function mb_forrer2021 = ...
    mb_forrer2021_5lvl(name_mat, num_opt, num_tun, period_mb)
%mb_forrer2021_5lvl creates mbsolve python input files for a
% THz harmonic comb test setup.
%
%    mb_forrer2021(name, num_opt, num_tun, period_mb) creates the input 
%    files for a terahertz harmonic comb test setup based on the active 
%    gain medium given in Forrer et al., 2021 
%    (https://doi.org/10.1063/5.0041339).
%    The MAT-file name_mat includes carrier transport simulation results,
%    which are saved in class objects defined in the module results. 
%    The number of optical and tunneling transitions can be specified  
%    as inputs, and are selected automatically. As an alternative, 
%    one can define the transitions within a cell array and use them as 
%    input arguments specifing the corresponding transitions. 
%    In the input variable period_mb the indices of the reduced active 
%    quantum system are stored.
% 
% Load simulation results.
load(name_mat);
% set waveguide
% effective refractive index.
neff = 3.6;
d.n_eff = neff;
r = (neff - 1) / (neff + 1);
w = waveguide(5e-3, 2200,[], 1, r);
d.set_waveguide(w);

mb_name = append(name_mat, '_', num2str(num_opt), 'opt', ...
    '_', num2str(num_tun), 'tun');
% Generate mb_forrer2021 class object.
f_THz = 3.5e12; % in Hz
mb_forrer2021 = mbsolve_sim(mb_name, d, period_mb, f_THz, eigen, ...
    carr_dist, deph, sc);

% Find optical THz transition
if(iscell(num_opt))
    trans_opt = num_opt;
    num_opt = length(trans_opt);
else
    trans_opt = mb_forrer2021.find_optical_trans(f_THz, num_opt);
end

for i = 1:num_opt
    ni = trans_opt{i}(1);
    nf = trans_opt{i}(2);
    mb_forrer2021.add_dipole_pair(ni, nf);
end

% Add tunneling transition
if(iscell(num_tun))
    trans_tun = num_tun;
    num_tun = length(trans_tun);
else
    trans_tun = mb_forrer2021.find_tunneling_trans(num_tun);
end

for i = 1:num_tun
    ni = trans_tun{i}(1);
    nf = trans_tun{i}(2);
    mb_forrer2021.add_tunnel_pair(ni, nf);
end

% Generate python input file.
mb_forrer2021.generate(pwd)

end