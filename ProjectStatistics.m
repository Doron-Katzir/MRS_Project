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

        function result = WishartCovarianceLRT(wishartData, statsCfg)
            %WishartCovarianceLRT Run the three covariance LRT modes.
            %
            % wishartData is the centrally prepared analysisData.wishart
            % view. Patient, metabolite-panel, and minimum-observation
            % decisions must be completed by ApplyAnalysisFilters.

            if nargin < 2 || isempty(statsCfg)
                statsCfg = struct();
            end

            ProjectStatistics.ValidateWishartData(wishartData);
            statsCfg = ProjectStatistics.NormalizeWishartConfig(statsCfg);

            legacyInput = struct();
            legacyInput.kind = "MRSAnalysisData";
            legacyInput.wishart = wishartData;
            legacyInput.crlbMajorityTable = wishartData.crlbMajorityTable;

            modeNames = ["modeA", "modeB", "modeC"];
            selectionModes = [ ...
                "perPatientLargestValid", "fixed", "largestCommon"];
            modeDisplayNames = [ ...
                "Per-patient largest valid subset", ...
                "Predefined fixed metabolite panel", ...
                "Largest common valid subset"];

            modes = struct();
            combinedSummary = table();
            combinedPatients = table();
            modeDiagnostics = struct();
            outputFiles = struct();

            for modeIdx = 1:numel(modeNames)
                modeName = modeNames(modeIdx);
                modeField = char(modeName);
                legacyCfg = struct();
                legacyCfg.metaboliteSelectionMode = selectionModes(modeIdx);
                legacyCfg.filterView = modeName;
                legacyCfg.runOnlyCommonPanelPatients = ...
                    statsCfg.runOnlyCommonPanelPatients;
                legacyCfg.minMetabolites = statsCfg.minMetabolites;
                legacyCfg.alpha = statsCfg.alpha;
                legacyCfg.applyNumericalRidge = ...
                    statsCfg.applyNumericalRidge;
                legacyCfg.ridgeScale = statsCfg.ridgeScale;
                legacyCfg.maxRidgeSteps = statsCfg.maxRidgeSteps;
                legacyCfg.exportResults = statsCfg.exportResults;
                legacyCfg.outputDir = statsCfg.outputDirs.(modeField);

                legacy = TestWishartCovarianceLRT(legacyInput, legacyCfg);

                modeParameters = struct();
                modeParameters.metaboliteSelectionMode = ...
                    legacy.options.metaboliteSelectionMode;
                modeParameters.runOnlyCommonPanelPatients = ...
                    legacy.options.runOnlyCommonPanelPatients;
                modeParameters.minMetabolites = legacy.options.minMetabolites;
                modeParameters.alpha = legacy.options.alpha;
                modeParameters.applyNumericalRidge = ...
                    legacy.options.applyNumericalRidge;
                modeParameters.ridgeScale = legacy.options.ridgeScale;
                modeParameters.maxRidgeSteps = legacy.options.maxRidgeSteps;
                modeParameters.exportResults = legacy.options.exportResults;
                modeParameters.outputDir = string(legacy.options.outputDir);

                diagnostics = struct();
                diagnostics.patientIDs = legacy.patientIDs;
                diagnostics.allCandidatePatientIDs = ...
                    legacy.allCandidatePatientIDs;
                diagnostics.preparedMetabolites = ...
                    wishartData.views.(modeField).metabolites;
                diagnostics.minValidParts = ...
                    wishartData.views.(modeField).minValidParts;
                diagnostics.nCandidatePatients = ...
                    numel(legacy.allCandidatePatientIDs);
                diagnostics.nPatientsTested = height(legacy.patientSummaryTable);
                diagnostics.nSuccessfulPatients = sum( ...
                    legacy.patientSummaryTable.status == "ok");
                diagnostics.nFailedPatients = sum( ...
                    legacy.patientSummaryTable.status == "failed");
                diagnostics.covarianceDimensions = unique( ...
                    legacy.patientSummaryTable.nMetabolites( ...
                    legacy.patientSummaryTable.status == "ok"));
                diagnostics.nEmpiricalRidges = sum( ...
                    legacy.patientSummaryTable.ridgeEmpirical > 0);
                diagnostics.nModelRidges = sum( ...
                    legacy.patientSummaryTable.ridgeModel > 0);
                diagnostics.nPositiveDefiniteFailures = sum(contains( ...
                    lower(legacy.patientSummaryTable.errorMessage), ...
                    "positive definite"));
                diagnostics.commonPanel = legacy.commonPanel;

                modeFiles = struct( ...
                    'patientCsv', "", ...
                    'groupCsv', "", ...
                    'commonPanelDiscoveryCsv', "", ...
                    'commonPanelCsv', "");
                if legacy.options.exportResults
                    modeFiles.patientCsv = string(fullfile( ...
                        legacy.options.outputDir, ...
                        'Wishart_LRT_PerPatient.csv'));
                    modeFiles.groupCsv = string(fullfile( ...
                        legacy.options.outputDir, 'Wishart_LRT_Group.csv'));
                    if ~isempty(legacy.panelDiscoveryTable)
                        modeFiles.commonPanelDiscoveryCsv = string(fullfile( ...
                            legacy.options.outputDir, ...
                            'Wishart_LRT_CommonPanelDiscovery.csv'));
                    end
                    if ~isempty(legacy.commonPanel)
                        modeFiles.commonPanelCsv = string(fullfile( ...
                            legacy.options.outputDir, ...
                            'Wishart_LRT_CommonPanel.csv'));
                    end
                end

                modeResult = struct();
                modeResult.name = modeDisplayNames(modeIdx);
                modeResult.parameters = modeParameters;
                modeResult.summaryTable = legacy.groupTable;
                modeResult.patientTable = legacy.patientSummaryTable;
                modeResult.diagnostics = diagnostics;
                modeResult.outputFiles = modeFiles;
                modeResult.patientResultsByID = legacy.patientResultsByID;
                modeResult.commonPanel = legacy.commonPanel;
                modeResult.panelDiscoveryTable = legacy.panelDiscoveryTable;
                modes.(modeField) = modeResult;
                modeDiagnostics.(modeField) = diagnostics;
                outputFiles.(modeField) = modeFiles;

                modeSummary = addvars(legacy.groupTable, ...
                    repmat(modeName, height(legacy.groupTable), 1), ...
                    'Before', 1, 'NewVariableNames', 'mode');
                combinedSummary = [combinedSummary; modeSummary]; %#ok<AGROW>
                modePatients = addvars(legacy.patientSummaryTable, ...
                    repmat(modeName, height(legacy.patientSummaryTable), 1), ...
                    'Before', 1, 'NewVariableNames', 'mode');
                combinedPatients = [combinedPatients; modePatients]; %#ok<AGROW>
            end

            parameters = statsCfg;
            parameters.modeOrder = modeNames;

            diagnostics = struct();
            diagnostics.method = legacy.method;
            diagnostics.nullHypothesis = legacy.nullHypothesis;
            diagnostics.testStatistic = legacy.testStatistic;
            diagnostics.dfDefinition = legacy.dfDefinition;
            diagnostics.modeOrder = modeNames;
            diagnostics.modes = modeDiagnostics;

            result = struct();
            result.name = "Wishart covariance likelihood-ratio test";
            result.parameters = parameters;
            result.summaryTable = combinedSummary;
            result.patientTable = combinedPatients;
            result.diagnostics = diagnostics;
            result.outputFiles = outputFiles;
            result.modes = modes;
            result.crlbMajorityTable = wishartData.crlbMajorityTable;
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

        function ValidateWishartData(wishartData)
            if ~isstruct(wishartData) || ~isscalar(wishartData)
                error('ProjectStatistics:InvalidWishartData', ...
                    ['wishartData must be a scalar struct from ', ...
                    'analysisData.wishart.']);
            end

            requiredFields = { ...
                'patientIDs', 'views', 'crlbMajorityTable', ...
                'patientDataByID'};
            missingFields = requiredFields(~isfield(wishartData, requiredFields));
            if ~isempty(missingFields)
                error('ProjectStatistics:InvalidWishartData', ...
                    'wishartData is missing required field(s): %s.', ...
                    strjoin(missingFields, ', '));
            end

            requiredViews = {'modeA', 'modeB', 'modeC'};
            missingViews = requiredViews(~isfield(wishartData.views, requiredViews));
            if ~isempty(missingViews)
                error('ProjectStatistics:InvalidWishartData', ...
                    'wishartData.views is missing required mode(s): %s.', ...
                    strjoin(missingViews, ', '));
            end
        end

        function statsCfg = NormalizeWishartConfig(statsCfg)
            if ~isstruct(statsCfg) || ~isscalar(statsCfg)
                error('ProjectStatistics:InvalidWishartConfig', ...
                    'statsCfg must be a scalar struct.');
            end

            statsCfg = ProjectStatistics.SetDefaultField( ...
                statsCfg, 'minMetabolites', 2);
            statsCfg = ProjectStatistics.SetDefaultField( ...
                statsCfg, 'alpha', 0.05);
            statsCfg = ProjectStatistics.SetDefaultField( ...
                statsCfg, 'applyNumericalRidge', true);
            statsCfg = ProjectStatistics.SetDefaultField( ...
                statsCfg, 'ridgeScale', 1e-8);
            statsCfg = ProjectStatistics.SetDefaultField( ...
                statsCfg, 'maxRidgeSteps', 10);
            statsCfg = ProjectStatistics.SetDefaultField( ...
                statsCfg, 'runOnlyCommonPanelPatients', true);
            statsCfg = ProjectStatistics.SetDefaultField( ...
                statsCfg, 'exportResults', false);

            if ~isfield(statsCfg, 'outputDirs') || ...
                    ~isstruct(statsCfg.outputDirs) || ...
                    ~isscalar(statsCfg.outputDirs)
                statsCfg.outputDirs = struct();
            end
            statsCfg.outputDirs = ProjectStatistics.SetDefaultField( ...
                statsCfg.outputDirs, 'modeA', fullfile( ...
                pwd, 'WishartLRTResults_ModeA_PerPatientLargestValid'));
            statsCfg.outputDirs = ProjectStatistics.SetDefaultField( ...
                statsCfg.outputDirs, 'modeB', fullfile( ...
                pwd, 'WishartLRTResults_ModeB_PredefinedPanel'));
            statsCfg.outputDirs = ProjectStatistics.SetDefaultField( ...
                statsCfg.outputDirs, 'modeC', fullfile( ...
                pwd, 'WishartLRTResults_ModeC_LargestCommon'));
        end

        function s = SetDefaultField(s, fieldName, defaultValue)
            if ~isfield(s, fieldName) || isempty(s.(fieldName))
                s.(fieldName) = defaultValue;
            end
        end
    end
end
