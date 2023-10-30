function [d, s] = wolf2017(basis_sp)
%wolf2017 creates input files for a mid-infrared test setup.
%
%    wolf2017(folder) creates the strain-balanced
%    mid-infrared QCL setup (wavelength 8.2 µm, alias: EV2017) described in
%    the PhD thesis of J. Wolf (ETH Zuerich, 2017). It returns an object of
%    simulation device and scenario. The input argument basis_sp
%    specifies theeigenstates basis for the Schrödinger-Poisson solver.
%

%% define scenario and device
% barrier material
b = InAlAs(0.36);
% well material
w = InGaAs(0.58);

% set up period
period = {; ...
    layer(b, 31); ...
    layer(w, 30); ...
    layer(b, 19); ...
    layer(w, 29); ...
    layer(b, 16); ...
    layer(w, 32, 1.20276e+17); ...
    layer(b, 13, 1.20276e+17); ...
    layer(w, 38, 1.20276e+17); ...
    layer(b, 13); ...
    layer(w, 45); ...
    layer(b, 10); ...
    layer(w, 50); ...
    layer(b, 7); ...
    layer(w, 58); ...
    layer(b, 12); ...
    layer(w, 25); ...
    };

% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);

% temperature in K
temperature = 300;

% bias field in kV/cm
% Values in mV/period
bias = 47;
% Simulation time in s
t_sim = 10e-11;

% SP simulation basis
if (nargin == 0)
    basis_sp = 'tb';
end

% number of wavefunctions
num_wavefct = 9;

% set up scenario
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);