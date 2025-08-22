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

function f = hdf5isdataset(filename, dataset_name)
% Checks if specific dataset is contained in hdf5 file.
%
% Syntax:
%   f = hdf5isdataset(filename, dataset_name)
%
% Input Arguments:
%   filename (string): Name of hdf5 file.
%   varname (string): Name of dataset.

f = false;
if ~exist(filename, 'file')
    return
end

info = h5info(filename);
group_names = strsplit(dataset_name, "/");
dataset_name = group_names(end);
group = "";

for i = 2:length(group_names)
    % In the last iteration check for the dataset name.
    if i == length(group_names)
        for j = 1:length(info.Datasets)
            dataset = info.Datasets(j);
            if ~isempty(dataset)
                if dataset_name == string(dataset.Name)
                    f = true;
                    return
                end
            end
        end
        return
    end

    % Check for subgroups.
    group = group + "/" + group_names(i);
    subgroups = info.Groups;
    for j = 1:length(subgroups)
        if group == string(subgroups(j).Name)
            info = subgroups(j);
            break
        end
        if j == length(subgroups)
            return
        end
    end
end
end
