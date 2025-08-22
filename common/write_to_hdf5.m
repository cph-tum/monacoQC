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

function write_to_hdf5(filename, dataset_name, data)
% Adds new data to an existing hdf5 file. If the file already contains a
% data entry with the same name, the new data overwrites the old data.
%
% Syntax:
%   write_to_hdf5(filename, dataset_name, data)
%
% Input Arguments:
%   filename (string): Name of hdf5 file.
%   dataset_name (string): Name of dataset.
%   data (numeric array | string): Data to be added to the file.

if isempty(data)
    return
end

if ~hdf5isdataset(filename, dataset_name)
    if isa(data, "string")
        h5create(filename, dataset_name, size(data), "Datatype", "string");
    else
        % drop useless dimensions!
        if any(size(data) == 1)
            shape = numel(data);
        else
            shape = size(data);
        end
        h5create(filename, dataset_name, shape);
    end
end

h5write(filename, dataset_name, data)

end
