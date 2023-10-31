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
