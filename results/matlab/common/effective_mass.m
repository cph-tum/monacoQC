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

classdef effective_mass < handle
    %effective_mass Contains the effective masses of the specified levels
    % in a QCL system.
    properties (SetAccess = private)
        m_eff; % Vector including the effective masses of the eigenstates.
    end
    
    methods
        function obj = effective_mass(m)
            % Constructs effective_mass.
            obj.m_eff = m; % Effective mass
        end
        
        function m_eff_i = get_m_eff(obj, ind)
            % Returns effective mass of the given state index.
            n_wf = size(obj.m_eff, 1);
            ni = mod(ind-1, n_wf) + 1;
            m_eff_i = obj.m_eff(ni);
        end
        
    end
    
end
