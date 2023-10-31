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

function err = plot3D_CBO_strain(dir, well, barr, substrate, temp, orientation)

%% Conduction Band Offset Ternary material composition
y_well = well.get_conc();
x_barr = barr.get_conc();
x = 0.01:0.01:0.99;
y = x;

%% Calculate CBO, if material system is unstrained
well.set_strain(0);
barr.set_strain(0);
CBO_unstrained = barr.get_Ec(substrate, temp, orientation) - ...
    well.get_Ec(substrate, temp, orientation);

%% Calculate CBO for all material compositions including strain
well.set_strain(1);
barr.set_strain(1);
for i = 1:length(y)
    well.set_conc(y(i));
    Ec1(i) = well.get_Ec(substrate, temp, orientation);
    for j = 1:length(x)
        barr.set_conc(x(j));
        Ec2(j) = barr.get_Ec(substrate, temp, orientation);
        CBO(i, j) = Ec2(j) - Ec1(i);
    end
end

%% Contour plot CBO
[X, Y] = meshgrid(x, y);
CBO_max = round(max(max(CBO)), 1);
CBO_min = round(min(min(CBO)), 1);
range = linspace(CBO_min, CBO_max, 20);
range = round(range, 1);
fig = figure('units', 'centimeters');
papersize = [12, 9];
contour(X, Y, CBO, range, 'ShowText', 'on', 'LineWidth', 1);
hold on;
i_y = min(find(y >= y_well, 1));
i_x = min(find(x >= x_barr, 1));
plot3(x(i_x), y(i_y), CBO_unstrained, '.', 'MarkerSize', 15, ...
    'Color', [0, 101, 189]/255);

str_ylabel = strcat(well.name(1:2), '_x', ...
    well.name(3:4), '_{1-x}', well.name(5:6));
str_xlabel = strcat(barr.name(1:2), '_x', ...
    barr.name(3:4), '_{1-x}', barr.name(5:6));
xlabel(str_xlabel);
ylabel(str_ylabel);

string_T = strcat('T=', num2str(temp), ' K');
text(0.05, 1.05, string_T, 'units', 'normalized');

%% Generate Pdf
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize', papersize);
string = strcat('CBO_', well.name, '_', barr.name, '.pdf');
path = fullfile(dir, string);
print(fig, path, '-dpdf', '-fillpage');

%% Calculate error
CBO_strained = Ec2(i_x) - Ec1(i_y);
err = abs(CBO_unstrained-CBO_strained) / CBO_unstrained;
if (err > 0.02)
    error(['The relative error of calculated CBO ', ...
        'for the unstrained case is to high (%2.2f %%)!'], err*100);
end
end
