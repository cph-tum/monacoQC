function plot_Ec_materials(dir_pdf, dir_mat, T)
% Plot the Conduction band edge in eV with respect to the value of GaAs

%% Find materials
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

%% Calculate Ec of ternaries
fig = figure('units', 'centimeters');
hold on;
size_ternary = 12;
size_binary = 15;
size_axis = 18;
papersize = [24, 18];
k = 1;
for i = 1:size(ternary, 2)
    ternary{1, i}.set_strain(0);
    if (ne(0.5, ternary{1, i}.conc))
        a0_ternary_lm(k) = ternary{1, i}.get_lattice_constant(T);
        Ec_ternary_lm(k) = ternary{1, i}.get_Ec0();
        ternary_name_lm{k} = ternary{1, i}.get_name();
        k = k + 1;
    end
    for j = 1:length(x)
        ternary{1, i}.set_conc(x(j));
        a0_sweep(j) = ternary{1, i}.get_lattice_constant(T);
        Ec_sweep(j) = ternary{1, i}.get_Ec0();
    end
    plot(a0_sweep, Ec_sweep, 'k--', 'LineWidth', 1);
end

%% Calculate Ec of binaries
for i = 1:size(binary, 2)
    a0_binary(i) = binary{1, i}.get_lattice_constant(T);
    Ec_binary(i) = binary{1, i}.get_Ec0();
    plot(a0_binary(i), Ec_binary(i), '.', 'MarkerSize', 15, ...
        'Color', [0, 101, 189]/255);
    text(a0_binary(i)+5e-3, Ec_binary(i), binary{1, i}.get_name(), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontSize', size_binary);
end

for i = 1:length(a0_ternary_lm)
    text(a0_ternary_lm(i)+5e-3, Ec_ternary_lm(i), ...
        ternary_name_lm{i}, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'FontSize', size_ternary);
    plot(a0_ternary_lm(i), Ec_ternary_lm(i), '.', 'MarkerSize', 15, ...
        'Color', [227, 114, 34]/255);
end
string_T = strcat('T=', num2str(T), ' K');
text(0.1, .1, string_T, 'Units', 'normalized');
xlabel('Lattice constant in Å', 'FontSize', size_axis);
ylabel('Gamma-valley Ec w.r.t. GaAs in eV', 'FontSize', size_axis);
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize', papersize);
ax = gca;
ax.FontSize = size_binary;
pdf_name = strcat('Ec_materials_', num2str(T), 'K');
path = fullfile(dir_pdf, pdf_name);
print(fig, path, '-dpdf', '-fillpage');
end
