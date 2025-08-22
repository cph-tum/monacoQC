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

classdef BO_sim < carrier_transport_sim
    %BO_sim Creates simulation setup for a n-dimensional
    % Bayesian optimization problem of a QCD.

    properties
        variedLayers = []; % Layer to be varied.
        var_int_layers = []; % Variation intervals of varied layers in
        % Angstrom.
        subgrid_BO = 0.1; % Investigated subgrid for the
        % specified variation interval  of the varied layers.
        basePeriod % Period of an active quantum casade (QC) device.
        numPeriod = 0; % Number of periods.
    end

    methods
        function obj = BO_sim(variedLayers, var_int_layers, subgrid_BO, ...
                basePeriod, numPeriod, numWavefunctions)
            % Constructs an instance of this class.
            obj.variedLayers = variedLayers;
            obj.var_int_layers = var_int_layers;
            obj.subgrid_BO = subgrid_BO;
            obj.basePeriod = basePeriod;
            obj.numPeriod = numPeriod;
            obj.num_wavefct = numWavefunctions;
        end

        function result = evaluate_next_device(obj, x, scenario, ...
                sim_options, dir, fun)
            % Evaluates next device with changed layer sequence.
            period = {};
            for i = 1:size(obj.basePeriod, 1)
                period{end+1} = copy(obj.basePeriod{i});
            end

            str_path = [];
            % Change length of specific layers of the QC device period.
            for m = 1:length(obj.variedLayers)
                period{obj.variedLayers(m)}.length = ...
                    period{obj.variedLayers(m)}.length ...
                    +x.(m) * obj.subgrid_BO;
                str_path = strcat(str_path, 'x', num2str(m-1), '_', ...
                    num2str(period{obj.variedLayers(m)}.length));
            end

            % Make directory.
            path = fullfile(dir, str_path);
            mkdir(path);
            % Construct device object.
            d = device(period, obj.numPeriod);
            % Generate input files.
            obj.generate(d, scenario, sim_options, path);
            % Calculate result of the given merit function.
            result = fun(path);
        end

        function results = bayes_opt(obj, fun, scenario, sim_options, ...
                dir, varargin)
            % Run Bayesian optimization. Introduce merit function @fun,
            % scenario.
            if nargin > 4
                [varargin{:}] = convertStringsToChars(varargin{:});
            end
            % Use variation intervals and take into account the
            % subgrid_BO for the specification optimizableVariable for the
            % varied layers.
            for i = 1:length(obj.variedLayers)
                VariableDesc(i) = optimizableVariable(strcat('d', ...
                    num2str(i)), [obj.var_int_layers(i, 1), ...
                    obj.var_int_layers(i, 2)]/obj.subgrid_BO, ...
                    'Type', 'integer');
            end
            % Define objective function.
            ObjectiveFcn = @(x)obj.evaluate_next_device(x, ...
                scenario, sim_options, dir, fun);
            C = bayesoptim.suppressWarnings();
            % Define options for the Bayesian optimization.
            Options = bayesoptim.BayesoptOptions(ObjectiveFcn, ...
                VariableDesc, varargin);
            % Start Bayesian optimization.
            results = BayesianOptimization(Options);
        end
    end
end
