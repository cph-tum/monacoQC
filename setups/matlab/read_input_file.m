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

function [d, s] = read_input_file(name)

%% read input file in .yaml format
input = yaml.ReadYaml(char(name+".yaml"));

%% create materials
for m = input.structure.materials
    % create material
    if strcmp(m{1}.composition, "binary")
        materials.(m{1}.name) = binary(m{1}.material);
    elseif strcmp(m{1}.composition, "ternary")
        materials.(m{1}.name) = ternary(m{1}.material, m{1}.ratio);
    end
    % set temperature
    materials.(m{1}.name).temp = input.scenario.temperature;
    % set strain information (defaults to false)
    if isfield(m{1}, 'strain')
        materials.(m{1}.name).strain = m{1}.strain;
    else
        materials.(m{1}.name).strain = false;
    end
end

%% create period
period = {};
for p = input.structure.period
    if ~isfield(p{1}, 'doping')
        p{1}.doping = 0;
    end

    period{end+1} = layer(materials.(p{1}.material), ...
        p{1}.thickness, p{1}.doping, p{1}.material);
end

%% create device
d = device(period, input.scenario.num_periods);
% set device substrate
d.substrate = materials.substrate;
% set material system (does not need custom set function)
d.material_system = input.structure.material_system;

%% create scenario
s = scenario(input.scenario.temperature, ...
    input.scenario.bias_field, ...
    input.scenario.simulation_time, ...
    input.scenario.num_wavefunctions, ...
    input.scenario.wavefunction_basis);

%% create waveguide
if isfield(input, 'waveguide')
    cav_l = input.waveguide.cavity_length;
    if isfield(input.waveguide, 'cavity_width')
        cav_w = input.waveguide.cavity_width;
    else
        cav_w = 0;
    end
    if isfield(input.waveguide, 'cavity_height')
        cav_h = input.waveguide.cavity_height;
    else
        cav_h = 0;
    end
    r_left = input.waveguide.field_reflectivity_left;
    r_right = input.waveguide.field_reflectivity_right;
    a_w = input.waveguide.power_loss;
    overlap = input.waveguide.overlap_factor;

    % create instance of class waveguide
    wg = laser_waveguide(cav_l, a_w, cav_h, overlap, ...
        r_left, r_right, cav_w);

    % add waveguide to device
    d = d.set_waveguide(wg);
end
