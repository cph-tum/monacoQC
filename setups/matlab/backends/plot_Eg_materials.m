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

function plot_Eg_materials(dir_pdf, dir_mat, T)
% Calculate Conduction band edge in eV with respect to the value of GaAs

%% Find all given materials
materials = dir(dir_mat);
for i = 1:length(materials)
    n_mat{i} = extractBefore(materials(i).name, '.m');
end
n_mat = n_mat(~cellfun(@isempty, n_mat));
x = 0.01:0.01:0.99;
binary = {};
ternary = {};
mat = {};
for i = 1:length(n_mat)
    mat{i} = eval(strcat(n_mat{1, i}, '()'));
    if (mat{i}.n_comp == 2)
        binary{end+1} = mat{i};
    elseif (mat{i}.n_comp == 3)
        ternary{end+1} = mat{i};
    end
end

%% Plot Ec of given materials
size_ternary = 12;
size_binary = 15;
size_axis = 18;
fig = figure('units', 'centimeters');
hold on;
papersize = [24, 18];
% Ternary
k = 1;
for i = 1:size(ternary, 2)
    ternary{1, i}.set_strain(0);
    
    if (ne(0.5, ternary{1, i}.conc))
        a0_ternary_lm(k) = ternary{1, i}.get_lattice_constant(T);
        Eg_ternary_lm(k) = ternary{1, i}.get_Eg0(T);
        ternary_name_lm{k} = ternary{1, i}.get_name();
        k = k + 1;
    end
    for j = 1:length(x)
        ternary{1, i}.set_conc(x(j));
        a0_sweep(j) = ternary{1, i}.get_lattice_constant(T);
        Eg_sweep(j) = ternary{1, i}.get_Eg0(T);
    end
    plot(a0_sweep, Eg_sweep, 'k--', 'LineWidth', 1);
end
% Binary
delta_a = 5e-3;
for i = 1:size(binary, 2)
    a0_binary(i) = binary{1, i}.get_lattice_constant(T);
    Eg_binary(i) = binary{1, i}.get_Eg0(T);
    plot(a0_binary(i), Eg_binary(i), '.', 'MarkerSize', 15, ...
        'Color', [0, 101, 189]/255);
    text(a0_binary(i)+5e-3, Eg_binary(i), binary{1, i}.get_name(), ...
        'VerticalAlignment', 'bottom', 'FontSize', size_binary);
end

for i = 1:length(a0_ternary_lm)
    if (strcmp(ternary_name_lm{i}, AlGaAs().get_name()))
        text(a0_ternary_lm(i)+5e-3, Eg_ternary_lm(i), ...
            ternary_name_lm{i}, 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'right', 'FontSize', size_ternary);
    elseif (strcmp(ternary_name_lm{i}, InGaP().get_name()))
        text(a0_ternary_lm(i)+5e-3, Eg_ternary_lm(i), ...
            ternary_name_lm{i}, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left', 'FontSize', size_ternary);
    elseif (strcmp(ternary_name_lm{i}, InAlAs().get_name()))
        text(a0_ternary_lm(i)+5e-3, Eg_ternary_lm(i), ...
            ternary_name_lm{i}, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left', 'FontSize', size_ternary);
    else
        text(a0_ternary_lm(i)+5e-3, Eg_ternary_lm(i), ...
            ternary_name_lm{i}, 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'left', 'FontSize', size_ternary);
    end
    plot(a0_ternary_lm(i), Eg_ternary_lm(i), '.', 'MarkerSize', 15, ...
        'Color', [227, 114, 34]/255);
end
string_T = strcat('T=', num2str(T), ' K');
text(0.1, 0.1, string_T, 'Units', 'normalized');
xlabel('Lattice constant in Å', 'FontSize', size_axis);
ylabel('Gamma-valley Eg in eV', 'FontSize', size_axis);
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize', papersize);
ax = gca;
ax.FontSize = size_binary;
pdf_name = strcat('Eg_materials_', num2str(T), 'K');
path = fullfile(dir_pdf, pdf_name);
print(fig, path, '-dpdf', '-fillpage');
end
