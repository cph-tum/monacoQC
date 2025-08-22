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

classdef ez_base_transform < handle
    % Base transformation class for ez eigenstates.

    methods (Static)
        function marker = ...
                find_multiplets(eig, e_multiplet, d_multiplet)
            % Determines multiplets (groups of energy eigenstates that each
            % differ by less than the predefined energy e_mulitplet).
            %
            % Syntax:
            %   marker = find_multiplets(eig, e_multiplet, d_multiplet)
            %
            % Input Arguments:
            %   eig (eigenstates-object): Contains information about the
            %     wavefunctions, eigenenergies and effective masses.
            %   e_multiplet (scalar): Maximum energy difference of
            %     multiplets [J].
            %   d_multiplet (scalar): Minimum dipole moment of multiplets
            %     [Cm].
            %
            % Output Arguments:
            %   marker (matrix): Contains information about multiplet
            %     states to be transformed.

            eef = eig.E;
            dipoles = eig.dipoles;

            [eef_sorted, sort_idx] = sort(eef);
            dE_bool = abs(diff(eef_sorted)) < e_multiplet;
            dE_bool(end+1) = 0;
            % Set up matrix marking states that should be transformed in ez
            marker = zeros(length(eef));
            ones_tuple = [];
            for i = 1:(length(dE_bool) - 1)
                if dE_bool(i) == true
                    ones_tuple = [ones_tuple, sort_idx(i), sort_idx(i+1)];
                end

                ones_tuple = unique(ones_tuple);

                if dE_bool(i+1) == false
                    flip_tuple = flip(ones_tuple);
                    for j = ones_tuple
                        for n = flip_tuple
                            marker(j, n) = 1;
                        end
                    end
                    ones_tuple = [];
                end
            end
            % Return markers of multiplets.
            marker = marker .* (abs(dipoles) > d_multiplet * ...
                phys_const.e0);
        end

        function [E_non_degen, psi_ez, meff_ez, d_ez] = ...
                transform(eig, marker)
            % Transforms extended eigenstates into ez states with
            % multiplets of localized states.
            %
            % Syntax:
            %   [E_non_degen, psi_ez, meff_ez, d_ez] = transform(eig, marker)
            %
            % Input Arguments:
            %   eig (eigenstates-object): Contains information about the
            %     wavefunctions, eigenenergies and effective masses.
            %   marker (matrix): Contains information about multiplet
            %     states to be transformed.
            %
            % Output Arguments:
            %   E_non_degen (vector): Transformed, non-degenerate
            %     eigenenergies.
            %   psi_ez (matrix): Transformed wavefunctions.
            %   meff_ez (vector): Transformed effective masses.
            %   d_ez (matrix): Transformed dipole matrix.

            eef = eig.E;
            psief = eig.psi;
            meff = eig.m_eff;
            dz = eig.dipoles;
            num_wfs = 4 * eig.num_wfs;
            dz_multiplet = dz .* marker;

            % transform to new basis
            d_ez = zeros(size(dz_multiplet));
            E_ez = zeros(size(dz_multiplet));
            psi_ez = zeros(size(psief));
            meff_ez = zeros(size(meff));

            i = 1;
            used = [];
            while i <= num_wfs
                inds = mod(find(marker(i, :)), num_wfs);
                inds(inds == 0) = num_wfs;
                if ~isempty(inds) && ~ismember(i, used)

                    d_red = zeros(length(inds));
                    e_red = zeros(length(inds));
                    psi_red = zeros(length(psief(:, 1)), length(inds));
                    meff_red = zeros(length(inds), 1);

                    for j = 1:length(inds)
                        for n = 1:length(inds)
                            d_red(j, n) = dz_multiplet(inds(j), inds(n));
                        end
                        e_red(j, j) = eef(inds(j));
                        psi_red(:, j) = psief(:, inds(j));
                        meff_red(j) = meff(inds(j));
                    end

                    [v, d_trans] = eigs(d_red, size(d_red, 1), 'sr');
                    e_trans = (v' * e_red) * (v')^(-1);
                    psi_trans = psi_red * (v')^(-1);
                    meff_trans = (v.^2)' * meff_red;

                    for j = 1:length(inds)
                        for n = 1:length(inds)
                            E_ez(inds(j), inds(n)) = e_trans(j, n);
                        end
                        d_ez(inds(j), inds(j)) = d_trans(j, j);
                        psi_ez(:, inds(j)) = psi_trans(:, j);
                        meff_ez(inds(j)) = meff_trans(j);
                    end

                    diff_ind = diff(inds);
                    if (max(diff_ind) > 1)
                        for j = 1:length(inds)
                            if (diff_ind(j) > 1)
                                i = 1 + inds(j);
                                used = [used, inds(j+1:end)];
                                break;
                            end
                        end
                    else
                        i = 1 + inds(end);
                    end

                else
                    if ~ismember(i, used)
                        psi_ez(:, i) = psief(:, i);
                        E_ez(i, i) = eef(i);
                        meff_ez(i) = meff(i);
                        d_ez(i, i) = dz(i);
                    end
                    i = i + 1;
                end
            end

            % handle degenerate states via uncertainty of tunneling barrier
            % abs to get rid of possible negative tunneling energies
            E_non_degen = abs(E_ez);
        end
    end
end
