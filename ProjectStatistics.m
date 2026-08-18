classdef ProjectStatistics
    %ProjectStatistics Public interface for project statistical analyses.

    methods (Static)
        function result = PairwiseEmpiricalVsModel(pairwiseData, statsCfg)
            %PairwiseEmpiricalVsModel Compare empirical and model correlations.
            %
            % pairwiseData is the centrally prepared analysisData.pairwise
            % view. Filtering and eligibility decisions must be completed by
            % ApplyAnalysisFilters before this method is called.

            if nargin < 2 || isempty(statsCfg)
                statsCfg = struct();
            end

            ProjectStatistics.ValidatePairwiseData(pairwiseData);
            statsCfg = ProjectStatistics.NormalizePairwiseConfig(statsCfg);

            % The validated implementation currently accepts the complete
            % analysis-data envelope. Supply only the focused pairwise view;
            % the other fields are compatibility outputs and do not take part
            % in any pairwise calculation.
            legacyInput = struct();
            legacyInput.kind = "MRSAnalysisData";
            legacyInput.pairwise = pairwiseData;
            legacyInput.crlbMajorityTable = table();
            legacyInput.filterReport = struct();

            legacy = TestPairwiseEmpiricalVsModelCorrelation( ...
                legacyInput, statsCfg);

            diagnostics = struct();
            diagnostics.method = legacy.method;
            diagnostics.nullHypothesis = legacy.nullHypothesis;
            diagnostics.patientIDs = legacy.patientIDs;
            diagnostics.metabolites = legacy.metabolites;
            diagnostics.nPairs = height(legacy.pairSummaryTable);
            diagnostics.nSuccessfulPairs = sum( ...
                legacy.pairSummaryTable.status == "ok");
            diagnostics.nFailedPairs = sum( ...
                legacy.pairSummaryTable.status == "failed");
            diagnostics.nPatientRows = height(legacy.patientPairTable);
            diagnostics.nSuccessfulPatientRows = sum( ...
                legacy.patientPairTable.status == "ok");
            diagnostics.nFailedPatientRows = sum( ...
                legacy.patientPairTable.status == "failed");

            outputFiles = struct('summaryCsv', "", 'patientCsv', "");
            if legacy.options.exportResults
                outputFiles.summaryCsv = string(fullfile( ...
                    legacy.options.outputDir, ...
                    'Pairwise_Empirical_vs_Model_Summary.csv'));
                outputFiles.patientCsv = string(fullfile( ...
                    legacy.options.outputDir, ...
                    'Pairwise_Empirical_vs_Model_PatientLevel.csv'));
            end

            result = struct();
            result.name = "Pairwise empirical vs model correlation";
            result.parameters = legacy.options;
            result.parameters.useFDR = statsCfg.useFDR;
            result.summaryTable = legacy.pairSummaryTable;
            result.patientTable = legacy.patientPairTable;
            result.diagnostics = diagnostics;
            result.outputFiles = outputFiles;
        end
    end

    methods (Static, Access = private)
        function ValidatePairwiseData(pairwiseData)
            if ~isstruct(pairwiseData) || ~isscalar(pairwiseData)
                error('ProjectStatistics:InvalidPairwiseData', ...
                    'pairwiseData must be a scalar struct from analysisData.pairwise.');
            end

            requiredFields = { ...
                'patientIDs', 'metabolites', 'pairTable', ...
                'patientEligibilityTable', 'minValidParts', ...
                'minPatientsForGroupTest'};
            missingFields = requiredFields(~isfield(pairwiseData, requiredFields));
            if ~isempty(missingFields)
                error('ProjectStatistics:InvalidPairwiseData', ...
                    'pairwiseData is missing required field(s): %s.', ...
                    strjoin(missingFields, ', '));
            end
        end

        function statsCfg = NormalizePairwiseConfig(statsCfg)
            if ~isstruct(statsCfg) || ~isscalar(statsCfg)
                error('ProjectStatistics:InvalidPairwiseConfig', ...
                    'statsCfg must be a scalar struct.');
            end

            if ~isfield(statsCfg, 'useFDR')
                statsCfg.useFDR = true;
            end
            useFDRIsTrue = (islogical(statsCfg.useFDR) || ...
                isnumeric(statsCfg.useFDR)) && isscalar(statsCfg.useFDR) && ...
                isfinite(double(statsCfg.useFDR)) && ...
                double(statsCfg.useFDR) == 1;
            if ~useFDRIsTrue
                error('ProjectStatistics:UnsupportedPairwiseConfig', ...
                    ['The validated pairwise implementation always applies ', ...
                    'Benjamini-Hochberg FDR correction; useFDR must be true.']);
            end
            statsCfg.useFDR = true;
        end
    end
end
