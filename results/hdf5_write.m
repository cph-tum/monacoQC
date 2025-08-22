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

function hdf5_write(filename, varargin)
% Saves objects to hdf5 file.
%
% Syntax:
%   hdf5_write(filename, obj1, obj2, ...)
%
% Input Arguments:
%   filename (string): Name of hdf5 file.
%   obj (postprocessing-object): Postprocessing object for which
%     to_hdf5() method is implemented.

for i = 1:length(varargin)
    varargin{i}.to_hdf5(filename);
end
end
