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

        function result = TemporalCorrelation(temporalData, statsCfg)
            %TemporalCorrelation Run temporal metabolite-correlation tests.
            %
            % temporalData is the centrally prepared analysisData.temporal
            % view. Filtering and eligibility decisions must be completed by
            % ApplyAnalysisFilters before this method is called.

            if nargin < 2 || isempty(statsCfg)
                statsCfg = struct();
            end

            ProjectStatistics.ValidateTemporalData(temporalData);
            statsCfg = ProjectStatistics.NormalizeTemporalConfig(statsCfg);

            % Keep TestTemporalMetaboliteCorrelations authoritative. Its
            % centralized-input adapter accepts this focused view inside the
            % smallest compatible analysis-data envelope.
            legacyInput = struct();
            legacyInput.kind = "MRSAnalysisData";
            legacyInput.temporal = temporalData;

            legacy = TestTemporalMetaboliteCorrelations( ...
                legacyInput, statsCfg);

            summaryTable = table();
            if isfield(legacy, 'groupFisherTable')
                summaryTable = legacy.groupFisherTable;
            end
            patientTable = table();
            if isfield(legacy, 'patientCorrelationTable')
                patientTable = legacy.patientCorrelationTable;
            end
            permutationResults = table();
            if isfield(legacy, 'circularShiftTable')
                permutationResults = legacy.circularShiftTable;
            end

            parameters = struct();
            parameterFields = { ...
                'permutationPairs', 'doFisherGroupTest', ...
                'doCircularShiftPermutation', 'nGroupPermutations', ...
                'rngSeed', 'useFDR', 'exportTables', 'outputDir'};
            for idx = 1:numel(parameterFields)
                fieldName = parameterFields{idx};
                parameters.(fieldName) = legacy.config.(fieldName);
            end
            parameters.fdrMethod = "Benjamini-Hochberg";

            diagnostics = struct();
            diagnostics.patientIDs = legacy.patientIDs;
            diagnostics.metabolites = legacy.metabList;
            diagnostics.nPatients = numel(legacy.patientIDs);
            diagnostics.nMetabolites = numel(legacy.metabList);
            diagnostics.nFisherTests = height(summaryTable);
            diagnostics.nSuccessfulFisherTests = 0;
            diagnostics.nFailedFisherTests = 0;
            if ismember('pValue', summaryTable.Properties.VariableNames)
                diagnostics.nSuccessfulFisherTests = sum(isfinite(summaryTable.pValue));
                diagnostics.nFailedFisherTests = sum(~isfinite(summaryTable.pValue));
            end
            diagnostics.nPatientRows = height(patientTable);
            diagnostics.nPatientRowsUsed = 0;
            if ismember('usedInGroupTest', patientTable.Properties.VariableNames)
                diagnostics.nPatientRowsUsed = sum(patientTable.usedInGroupTest);
            end
            diagnostics.nPermutationTests = height(permutationResults);
            diagnostics.nSuccessfulPermutationTests = 0;
            diagnostics.nFailedPermutationTests = 0;
            if ismember('groupCircularShiftPValue', ...
                    permutationResults.Properties.VariableNames)
                successful = isfinite( ...
                    permutationResults.groupCircularShiftPValue);
                diagnostics.nSuccessfulPermutationTests = sum(successful);
                diagnostics.nFailedPermutationTests = sum(~successful);
            end

            outputFiles = struct( ...
                'groupFisherCsv', "", ...
                'patientCsv', "", ...
                'circularShiftCsv', "", ...
                'crlbQualityCsv', "");
            if legacy.config.exportTables
                if legacy.config.doFisherGroupTest
                    outputFiles.groupFisherCsv = string(fullfile( ...
                        legacy.config.outputDir, 'Group_Fisher_Z_Tests.csv'));
                    outputFiles.patientCsv = string(fullfile( ...
                        legacy.config.outputDir, 'Patient_Level_Correlations.csv'));
                end
                if legacy.config.doCircularShiftPermutation
                    outputFiles.circularShiftCsv = string(fullfile( ...
                        legacy.config.outputDir, ...
                        'Circular_Shift_Permutation_Tests.csv'));
                end
                outputFiles.crlbQualityCsv = string(fullfile( ...
                    legacy.config.outputDir, 'CRLB_Reliability_Table.csv'));
            end

            result = struct();
            result.name = "Temporal metabolite correlation";
            result.parameters = parameters;
            result.summaryTable = summaryTable;
            result.patientTable = patientTable;
            result.diagnostics = diagnostics;
            result.outputFiles = outputFiles;
            result.permutationResults = permutationResults;
            result.crlbQualityTable = legacy.groupCRLBQualityTable;
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

        function ValidateTemporalData(temporalData)
            if ~isstruct(temporalData) || ~isscalar(temporalData)
                error('ProjectStatistics:InvalidTemporalData', ...
                    ['temporalData must be a scalar struct from ', ...
                    'analysisData.temporal.']);
            end

            requiredFields = { ...
                'patientIDs', 'metabolites', 'pairTable', ...
                'patientEligibilityTable', 'minValidParts', 'minPatients', ...
                'crlbThreshold', 'requiredGoodFraction', ...
                'crlbQualityTable', 'circularShiftPrepared', ...
                'empiricalMeanAbsCorrTable', ...
                'modelMeanAmplitudeCorrTable', ...
                'modelMeanAbsAmplitudeCorrTable'};
            missingFields = requiredFields(~isfield(temporalData, requiredFields));
            if ~isempty(missingFields)
                error('ProjectStatistics:InvalidTemporalData', ...
                    'temporalData is missing required field(s): %s.', ...
                    strjoin(missingFields, ', '));
            end
        end

        function statsCfg = NormalizeTemporalConfig(statsCfg)
            if ~isstruct(statsCfg) || ~isscalar(statsCfg)
                error('ProjectStatistics:InvalidTemporalConfig', ...
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
                error('ProjectStatistics:UnsupportedTemporalConfig', ...
                    ['The authoritative temporal implementation always ', ...
                    'applies Benjamini-Hochberg FDR correction; ', ...
                    'useFDR must be true.']);
            end
            statsCfg.useFDR = true;
        end
    end
end
