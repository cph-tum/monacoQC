function err = plot2D_CBO_strain(dir, binary, ternary, substrate, ...
    temp, orientation)

%% Conduction Band Offset Ternary material composition
x_ternary = ternary.get_conc();
x = 0.01:0.01:0.99;
binary.set_strain(0);
name_ternary = ternary.get_name();
Ec1_unstrained = binary.get_Ec(substrate, temp, orientation);
binary.set_strain(1);
Ec1_strained = binary.get_Ec(substrate, temp, orientation);
for j = 1:length(x)
    ternary.set_conc(x(j));
    Ec2_strained(j) = ternary.get_Ec(substrate, temp, orientation);
    ternary.set_strain(0);
    Ec2_unstrained(j) = ternary.get_Ec(substrate, temp, orientation);
    ternary.set_strain(1);
end
CBO_unstrained = Ec2_unstrained - Ec1_unstrained;
CBO_strained = Ec2_strained - Ec1_strained;

%% Plot CBO
fig = figure('units', 'centimeters');
papersize = [12, 9];
hold on;
plot(x, CBO_strained, 'Color', [162, 173, 0]/255, 'LineWidth', 1);
plot(x, CBO_unstrained, 'Color', [227, 114, 34]/255, 'LineWidth', 1);
i_x = min(find(x >= x_ternary, 1));
plot(x(i_x), CBO_unstrained(i_x), '.', 'MarkerSize', 15, ...
    'Color', [0, 101, 189]/255);
str_xlabel = strcat(ternary.name(1:2), '_x', ...
    ternary.name(3:4), '_{1-x}', ternary.name(5:6));
xlabel(str_xlabel);
ylabel('CBO in eV');
legend('strained', 'unstrained', 'Location', 'best');
if (CBO_unstrained(i_x+1) > CBO_unstrained(i_x))
    text(x(i_x), CBO_unstrained(i_x), name_ternary, ...
        'VerticalAlignment', 'top');
else
    text(x(i_x), CBO_unstrained(i_x), name_ternary, ...
        'VerticalAlignment', 'bottom');
end
string_T = strcat('T=', num2str(temp), ' K');
text(0.05, 1, string_T, 'units', 'normalized');

%% Generate Pdf
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize', papersize);
string = strcat('CBO_', binary.name, '_', ternary.name, '.pdf');
path = fullfile(dir, string);
print(fig, path, '-dpdf', '-fillpage');

%% Calculate error
err = abs((CBO_unstrained(i_x) - CBO_strained(i_x))/CBO_unstrained(i_x));
if (err > 0.02)
    error(['The relative error of calculated CBO ', ...
        'for the unstrained case is to high (%2.2f %%)!'], err*100);
end
end
