classdef ProjectResults
    %ProjectResults Public interface for scientific result interpretation.

    methods (Static)
        function result = Temporal(temporal, resultsCfg)
            %Temporal Select temporal candidates using configured policies.

            if nargin < 2 || isempty(resultsCfg)
                resultsCfg = struct();
            end

            ProjectResults.ValidateTemporalInput(temporal);
            parameters = ProjectResults.NormalizeTemporalConfig(resultsCfg);
            summaryTable = temporal.summaryTable;

            requiredColumns = { ...
                'qValue_FDR', 'pValue', 'signConsistency'};
            if parameters.strong.requireCRLBPass
                requiredColumns{end + 1} = 'CRLB_pair_status';
            end
            ProjectResults.RequireColumns( ...
                summaryTable, requiredColumns, 'temporal.summaryTable');

            strongRows = ...
                summaryTable.qValue_FDR < parameters.strong.maxQValue & ...
                summaryTable.signConsistency >= ...
                parameters.strong.minSignConsistency;
            if parameters.strong.requireCRLBPass
                strongRows = strongRows & ...
                    summaryTable.CRLB_pair_status == "PASS";
            end

            exploratoryRows = ...
                summaryTable.pValue < parameters.exploratory.maxPValue & ...
                summaryTable.signConsistency >= ...
                parameters.exploratory.minSignConsistency;

            result = struct();
            result.name = "Temporal correlation result interpretation";
            result.parameters = parameters;
            result.strongCandidates = summaryTable(strongRows, :);
            result.exploratoryCandidates = summaryTable(exploratoryRows, :);
            result.diagnostics = struct( ...
                'nSummaryRows', height(summaryTable), ...
                'nStrongCandidates', sum(strongRows), ...
                'nExploratoryCandidates', sum(exploratoryRows));
        end
    end

    methods (Static, Access = private)
        function ValidateTemporalInput(temporal)
            if ~isstruct(temporal) || ~isscalar(temporal) || ...
                    ~isfield(temporal, 'summaryTable') || ...
                    ~istable(temporal.summaryTable)
                error('ProjectResults:InvalidTemporalInput', ...
                    'temporal.summaryTable must be a table.');
            end
        end

        function parameters = NormalizeTemporalConfig(cfg)
            if ~isstruct(cfg) || ~isscalar(cfg)
                error('ProjectResults:InvalidTemporalConfig', ...
                    'resultsCfg must be a scalar struct.');
            end

            parameters = cfg;
            if ~isfield(parameters, 'strong') || ...
                    ~isstruct(parameters.strong) || ...
                    ~isscalar(parameters.strong)
                parameters.strong = struct();
            end
            if ~isfield(parameters, 'exploratory') || ...
                    ~isstruct(parameters.exploratory) || ...
                    ~isscalar(parameters.exploratory)
                parameters.exploratory = struct();
            end

            parameters.strong = ProjectResults.SetDefault( ...
                parameters.strong, 'maxQValue', 0.05);
            parameters.strong = ProjectResults.SetDefault( ...
                parameters.strong, 'requireCRLBPass', true);
            parameters.strong = ProjectResults.SetDefault( ...
                parameters.strong, 'minSignConsistency', 0.75);
            parameters.exploratory = ProjectResults.SetDefault( ...
                parameters.exploratory, 'maxPValue', 0.05);
            parameters.exploratory = ProjectResults.SetDefault( ...
                parameters.exploratory, 'minSignConsistency', 0.75);

            parameters.strong.requireCRLBPass = logical( ...
                parameters.strong.requireCRLBPass);
        end

        function RequireColumns(T, requiredColumns, tableName)
            missingColumns = requiredColumns(~ismember( ...
                requiredColumns, T.Properties.VariableNames));
            if ~isempty(missingColumns)
                error('ProjectResults:MissingTemporalColumns', ...
                    '%s is missing required column(s): %s.', ...
                    tableName, strjoin(missingColumns, ', '));
            end
        end

        function s = SetDefault(s, fieldName, defaultValue)
            if ~isfield(s, fieldName) || isempty(s.(fieldName))
                s.(fieldName) = defaultValue;
            end
        end
    end
end
