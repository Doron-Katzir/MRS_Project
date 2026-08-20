classdef ProjectExports
    %ProjectExports Public interface for project result exports.

    methods (Static)
        function result = Temporal(temporal, temporalResults, exportCfg)
            %Temporal Build and optionally write simplified temporal reports.

            if nargin < 3 || isempty(exportCfg)
                exportCfg = struct();
            end

            ProjectExports.ValidateTemporalInputs(temporal, temporalResults);
            parameters = ProjectExports.NormalizeTemporalConfig(exportCfg);

            mainColumns = { ...
                'metaboliteA', 'metaboliteB', 'groupMeanR', ...
                'pValue', 'qValue_FDR', 'nPatientsUsed', ...
                'nPositivePatients', 'nNegativePatients', ...
                'signConsistency', 'CRLB_pair_status', ...
                'absLCModelCorr'};
            mainColumns = ProjectExports.AvailableColumns( ...
                temporal.summaryTable, mainColumns);

            mainTable = temporal.summaryTable(:, mainColumns);
            if ismember('qValue_FDR', mainTable.Properties.VariableNames)
                mainTable = sortrows(mainTable, 'qValue_FDR', 'ascend');
            end

            strongTable = temporalResults.strongCandidates(:, mainColumns);
            if height(strongTable) > 0 && ...
                    ismember('qValue_FDR', strongTable.Properties.VariableNames)
                strongTable = sortrows(strongTable, 'qValue_FDR', 'ascend');
            end

            exploratoryTable = ...
                temporalResults.exploratoryCandidates(:, mainColumns);
            if height(exploratoryTable) > 0 && ...
                    ismember('pValue', exploratoryTable.Properties.VariableNames)
                exploratoryTable = sortrows( ...
                    exploratoryTable, 'pValue', 'ascend');
            end

            circularShiftTable = table();
            if ~isempty(temporal.permutationResults)
                shiftColumns = { ...
                    'metaboliteA', 'metaboliteB', 'observedGroupR', ...
                    'groupCircularShiftPValue', 'qValue_FDR', ...
                    'nPatientsUsed'};
                shiftColumns = ProjectExports.AvailableColumns( ...
                    temporal.permutationResults, shiftColumns);
                circularShiftTable = ...
                    temporal.permutationResults(:, shiftColumns);
                if ismember('qValue_FDR', ...
                        circularShiftTable.Properties.VariableNames)
                    circularShiftTable = sortrows( ...
                        circularShiftTable, 'qValue_FDR', 'ascend');
                    % Presentation-only compatibility rename. FDR remains
                    % authoritative in the temporal statistical result.
                    circularShiftTable.Properties.VariableNames{ ...
                        strcmp(circularShiftTable.Properties.VariableNames, ...
                        'qValue_FDR')} = 'qValueCircularShift_FDR';
                elseif ismember('groupCircularShiftPValue', ...
                        circularShiftTable.Properties.VariableNames)
                    circularShiftTable = sortrows(circularShiftTable, ...
                        'groupCircularShiftPValue', 'ascend');
                end
            end

            patientTable = table();
            if ~isempty(temporal.patientTable)
                patientColumns = { ...
                    'metaboliteA', 'metaboliteB', 'patientID', ...
                    'rValue', 'nValidParts'};
                patientColumns = ProjectExports.AvailableColumns( ...
                    temporal.patientTable, patientColumns);
                patientTable = temporal.patientTable(:, patientColumns);
            end

            crlbTable = table();
            if ~isempty(temporal.crlbQualityTable)
                crlbColumns = { ...
                    'metabolite', 'nCRLBUnder100', 'nInstances', ...
                    'fractionCRLBUnder100', 'fails90PercentRule'};
                crlbColumns = ProjectExports.AvailableColumns( ...
                    temporal.crlbQualityTable, crlbColumns);
                crlbTable = temporal.crlbQualityTable(:, crlbColumns);
                if ismember('fractionCRLBUnder100', ...
                        crlbTable.Properties.VariableNames)
                    crlbTable = sortrows( ...
                        crlbTable, 'fractionCRLBUnder100', 'ascend');
                end
            end

            tables = struct();
            tables.mainFisher = mainTable;
            tables.strongCandidates = strongTable;
            tables.exploratoryCandidates = exploratoryTable;
            tables.circularShift = circularShiftTable;
            tables.patientCorrelations = patientTable;
            tables.crlbReliability = crlbTable;

            outputFiles = ProjectExports.EmptyTemporalOutputFiles();
            if parameters.writeSimplifiedTables
                if ~isfolder(parameters.outputDir)
                    mkdir(parameters.outputDir);
                end

                outputFiles = ProjectExports.TemporalOutputFiles( ...
                    parameters.outputDir);
                writetable(mainTable, outputFiles.mainFisherCsv);
                writetable(strongTable, outputFiles.strongCandidatesCsv);
                writetable(exploratoryTable, ...
                    outputFiles.exploratoryCandidatesCsv);
                if ~isempty(circularShiftTable)
                    writetable(circularShiftTable, ...
                        outputFiles.circularShiftCsv);
                else
                    outputFiles.circularShiftCsv = "";
                end
                if ~isempty(patientTable)
                    writetable(patientTable, outputFiles.patientCsv);
                else
                    outputFiles.patientCsv = "";
                end
                if ~isempty(crlbTable)
                    writetable(crlbTable, outputFiles.crlbQualityCsv);
                else
                    outputFiles.crlbQualityCsv = "";
                end
            end

            if parameters.printSummary
                fprintf("\nFinished temporal correlation tests.\n");
                fprintf("Number of metabolite pairs tested with Fisher-z: %d\n", ...
                    height(temporal.summaryTable));
                if ~isempty(temporal.permutationResults)
                    fprintf(["Number of metabolite pairs tested with ", ...
                        "circular shift: %d\n"], ...
                        height(temporal.permutationResults));
                end
                if parameters.writeSimplifiedTables
                    fprintf("Saved simplified CSV files to:\n%s\n", ...
                        parameters.outputDir);
                end
            end

            result = struct();
            result.name = "Temporal simplified exports";
            result.parameters = parameters;
            result.tables = tables;
            result.outputFiles = outputFiles;
            result.diagnostics = struct( ...
                'nFilesWritten', sum(strlength(string( ...
                struct2cell(outputFiles))) > 0), ...
                'knownOptionalSchemaMismatches', [ ...
                "absLCModelCorr vs absLCModelAmplitudeCorr", ...
                "rValue vs r"]);
        end
    end

    methods (Static, Access = private)
        function ValidateTemporalInputs(temporal, temporalResults)
            requiredTemporalFields = { ...
                'summaryTable', 'patientTable', 'permutationResults', ...
                'crlbQualityTable'};
            if ~isstruct(temporal) || ~isscalar(temporal)
                error('ProjectExports:InvalidTemporalInput', ...
                    'temporal must be a scalar result struct.');
            end
            missingFields = requiredTemporalFields(~isfield( ...
                temporal, requiredTemporalFields));
            if ~isempty(missingFields)
                error('ProjectExports:InvalidTemporalInput', ...
                    'temporal is missing required field(s): %s.', ...
                    strjoin(missingFields, ', '));
            end

            requiredResultFields = { ...
                'strongCandidates', 'exploratoryCandidates'};
            if ~isstruct(temporalResults) || ~isscalar(temporalResults)
                error('ProjectExports:InvalidTemporalResults', ...
                    'temporalResults must be a scalar result struct.');
            end
            missingFields = requiredResultFields(~isfield( ...
                temporalResults, requiredResultFields));
            if ~isempty(missingFields)
                error('ProjectExports:InvalidTemporalResults', ...
                    'temporalResults is missing required field(s): %s.', ...
                    strjoin(missingFields, ', '));
            end
        end

        function parameters = NormalizeTemporalConfig(cfg)
            if ~isstruct(cfg) || ~isscalar(cfg)
                error('ProjectExports:InvalidTemporalConfig', ...
                    'exportCfg must be a scalar struct.');
            end
            parameters = ProjectExports.SetDefault( ...
                cfg, 'writeSimplifiedTables', true);
            parameters = ProjectExports.SetDefault(parameters, ...
                'outputDir', fullfile( ...
                pwd, 'TemporalCorrelationStats_Simplified_AllPairs'));
            parameters = ProjectExports.SetDefault( ...
                parameters, 'printSummary', true);
            parameters.writeSimplifiedTables = logical( ...
                parameters.writeSimplifiedTables);
            parameters.printSummary = logical(parameters.printSummary);
            parameters.outputDir = string(parameters.outputDir);
        end

        function columns = AvailableColumns(T, requestedColumns)
            columns = requestedColumns(ismember( ...
                requestedColumns, T.Properties.VariableNames));
        end

        function files = TemporalOutputFiles(outputDir)
            files = struct();
            files.mainFisherCsv = string(fullfile(outputDir, ...
                'Main_FisherZ_Results_Simplified.csv'));
            files.strongCandidatesCsv = string(fullfile(outputDir, ...
                'Strong_Candidates_q_FDR_CRLB_PASS.csv'));
            files.exploratoryCandidatesCsv = string(fullfile(outputDir, ...
                'Exploratory_Candidates_p_uncorrected.csv'));
            files.circularShiftCsv = string(fullfile(outputDir, ...
                'Circular_Shift_All_Pairs_Simplified.csv'));
            files.patientCsv = string(fullfile(outputDir, ...
                'Patient_Level_Correlations_Simplified.csv'));
            files.crlbQualityCsv = string(fullfile(outputDir, ...
                'CRLB_Reliability_Simplified.csv'));
        end

        function files = EmptyTemporalOutputFiles()
            files = struct( ...
                'mainFisherCsv', "", ...
                'strongCandidatesCsv', "", ...
                'exploratoryCandidatesCsv', "", ...
                'circularShiftCsv', "", ...
                'patientCsv', "", ...
                'crlbQualityCsv', "");
        end

        function s = SetDefault(s, fieldName, defaultValue)
            if ~isfield(s, fieldName) || isempty(s.(fieldName))
                s.(fieldName) = defaultValue;
            end
        end
    end
end
