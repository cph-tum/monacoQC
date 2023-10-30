function err = plot_meff_comp(dir, mat, substrate, temp, ...
    orientation)

%% CI-CD: Provide strain calculations of effective mass meff
x = 0.01:0.01:0.99;
x_mat = mat.get_conc();
name_mat = mat.get_name();
meffGamma = zeros(length(x), 1);
meffL_strain = zeros(length(x), 1);
meffII_strain = zeros(length(x), 1);

%% Effective mass calculations
for i = 1:length(x)
    mat.set_conc(x(i));
    meffL_strain(i) = mat.get_meffL(substrate, temp, orientation);
    meffII_strain(i) = mat.get_meffII(substrate, temp, orientation);
    mat.set_strain(0);
    meffGamma(i) = mat.get_meffL(substrate, temp, orientation);
    mat.set_strain(1);
end

%% Plot effective masses
fig = figure('units', 'centimeters');
hold on;
papersize = [12, 9];
plot(x, meffGamma, 'LineWidth', 1, 'Color', [162, 173, 0]/255)
plot(x, meffII_strain, 'LineWidth', 1, 'Color', [227, 114, 34]/255)
plot(x, meffL_strain, 'LineWidth', 1, 'Color', [218, 215, 203]/255)
i_x = min(find(x >= x_mat, 1));
plot(x(i_x), meffGamma(i_x), '.', 'MarkerSize', 15, ...
    'Color', [0, 101, 189]/255);
str_xlabel = strcat('Composition x (', mat.name(1:2), '_x', ...
    mat.name(3:4), '_{1-x}', mat.name(5:6), ')');
xlabel(str_xlabel);
ylabel('Effective mass in 1/m_0');

%% Generate Pdf
if (meffGamma(i_x+1) > meffGamma(i_x))
    text(x(i_x), meffGamma(i_x), name_mat, 'VerticalAlignment', 'top');
else
    text(x(i_x), meffGamma(i_x), name_mat, 'VerticalAlignment', 'bottom');
end
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize', papersize);
legend({'m_{eff,\Gamma}', 'm_{eff,II}', ...
    'm_{eff,L}'}, 'Location', 'best');
string_T = strcat('T=', num2str(temp), ' K');
text(0.05, 1, string_T, 'units', 'normalized');
string = strcat('meff_', mat.name, '.pdf');
path = fullfile(dir, string);
print(fig, path, '-dpdf', '-fillpage');

%% Calculate error
err = abs(meffGamma(i_x)-meffL_strain(i_x)) / meffGamma(i_x);
if (err > 0.02)
    error(['The relative error of calculated effective masses ', ...
        'for the unstrained case is to high (%2.2f %%)!'], err*100);
end
