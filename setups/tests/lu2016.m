function [d,s] =  lu2016(basis_sp, bias)
    %lu2016 creates an object of device and scenario for a THz DFG test setup.
    %
    %    lu2016(basis_sp, bias) creates generates device and scenario for the for the strain-balanced
    %    THz difference frequency generation QCL (wavelength 7.8µm) described 
    %    in Lu et al., 2016 (http://www.nature.com/articles/srep23595). 
    %    The resulting input files are saved in folder.

    % Add library to path.
    monaco_path = fullfile(pwd, '..', '..');
    addpath(genpath(fullfile(monaco_path, 'schrod-poisson')), ...
    genpath(fullfile(monaco_path, 'setups')));
    
    % barrier material
    b = InAlAs_CBO(0.37);
    % well materials
    w1 = InGaAs_CBO(0.65);
    w2 = InGaAs_CBO(0.53);
    
    % set up period
    period = {
        layer(b, 27); ...
        layer(w1, 28); ...
        layer(b, 24); ...
        layer(w1, 28); ...
        layer(b, 19); ...
        layer(w1, 28); ...
        layer(b, 17, 1.7e17); ...
        layer(w1, 13); ...
        layer(w2, 18); ...
        layer(b, 16, 1.7e17); ...
        layer(w1, 13); ...
        layer(w2, 20); ...
        layer(b, 15); ...
        layer(w1, 15); ...
        layer(w2, 22); ...
        layer(b, 15); ...
        layer(w1, 15); ...
        layer(w2, 24); ...
        layer(b, 17); ...
        layer(w1, 21); ...
        layer(w2, 30); ...
        layer(b, 9); ...
        layer(w1, 26); ...
        layer(w2, 32); ...
        layer(b, 9); ...
        layer(w1, 21); ...
        };
    
    % number of periods
    num_periods = 5;
    
    % set up device
    d = device(period, num_periods);
    % set waveguide 
    % effective refractive index.
    % neff = 3.15;
    % d.n_eff = neff; 
    % r = (neff - 1) / (neff + 1); 
    % w = waveguide_QC(3e-3, 640, 0.6, r);
    % d.set_waveguide(w);

    % temperature in K
    temperature = 293;
    
    % bias field in kV/cm
    if (nargin < 2)
        V = 48;
    else
        V = bias;
    end

    % SP simulation basis
    if (nargin <= 1)
       basis_sp = 'tb';
    end

    % Number of wavefunctions
    num_wavefct = 9;


    % Simulation time in s
    t_sim = 10e-11;


    % set up scenario
    s = scenario(temperature, V, t_sim, num_wavefct, basis_sp);
    
    
    end
    