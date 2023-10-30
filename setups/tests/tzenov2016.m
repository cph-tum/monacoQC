function [d, s] = tzenov2016(basis_sp)
%tzenov2016 creates input files for a terahertz test setup.
%
%    tzenov2016(folder) creates the input files for the terahertz QCL setup
%    simulated in Tzenov et al., 2016
%    (https://doi.org/10.1364/oe.24.023232).
%    The resulting input files are saved in folder. Structure is from
%    Ph.D. Thesis Qi Qin - Development of tunable terahertz
%    quantum cascade wire lasers; (4 THz, alias: FL183S)
%    (http://hdl.handle.net/1721.1/78548). In Tzenov2016, a modification in
%    the doped layer was conducted. The well layer was only partially doped
%    (50 Angstrom). It returns an object of
%    simulation device and scenario. The input argument basis_sp
%    specifies theeigenstates basis for the Schrödinger-Poisson solver.
%

%% define scenario and device
% barrier material
b = AlGaAs(0.15);
% well material
w = GaAs();

% set up period
period = {; ...
    layer(b, 48); ...
    layer(w, 92); ...
    layer(b, 33); ...
    layer(w, 56); ...
    layer(w, 50, 1.9e16); ...
    layer(w, 56); ...
    layer(b, 40); ...
    layer(w, 74); ...
    layer(b, 23); ...
    layer(w, 76); ...
    };

% number of periods
num_periods = 5;

% set up device
d = device(period, num_periods);

% temperature in K
temperature = 50;

% bias field in kV/cm
bias = 11;

% Simulation time in s
t_sim = 5e-11;

% SP simulation basis
if (nargin == 0)
    basis_sp = 'tb';
end
% number wavefunctions
num_wavefct = 5;

% set up scenario
s = scenario(temperature, bias, t_sim, num_wavefct, basis_sp);