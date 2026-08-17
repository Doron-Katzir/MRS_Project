%% Example caller for PlotTemporalCorrelationDiagnostics.m
% Run this after you already created:
%   covOutputs = MetabCovarianceByPatient(cfg);
%   deGraafOutputs = DeGraafAmplitudeCorrelationByPatient(cfg);

plotCfg = struct();

%% Patient selection
% Options:
%   "all"
%   ["P01", "P03", "P11"]
%   1:5
plotCfg.patientIDs = "all";

%% Matrix plots
plotCfg.doMatrixPlots = true;
plotCfg.doDifferenceMatrix = true;
plotCfg.doCRLBOverlayMatrix = true;
plotCfg.matrixMetabs = "all";

% Example subset:
% plotCfg.matrixMetabs = ["Glu", "Gln", "Glu+Gln", "GABA", "NAA+NAAG", "Cr+PCr"];

%% CRLB rule
plotCfg.crlbThreshold = 100;
plotCfg.requiredGoodFraction = 0.90;

%% Time-series plots
plotCfg.doTimeSeriesPlots = true;
plotCfg.timeSeriesGroups = { ...
    ["GPC", "PCh", "GPC+PCh"], ...
    ["Glu", "Gln", "Glu+Gln"], ...
    ["NAA", "NAAG", "NAA+NAAG"], ...
    ["Cr", "PCr", "Cr+PCr"]};

%% Per-patient scatter plots
plotCfg.doPatientScatterPlots = true;
plotCfg.scatterPairs = [ ...
    "Glu",     "GABA"; ...
    "Glu",     "Glu+Gln"; ...
    "GABA",    "Glu+Gln"];

%% Sum-metabolite scatter plots
plotCfg.doSumScatterPlots = true;
plotCfg.sumMetabs = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

%% Multi-patient pooled within-patient z-scored scatter
plotCfg.doPooledZScatter = true;
plotCfg.pooledZPairs = [ ...
    "Glu",     "GABA"; ...
    "Glu+Gln", "GABA"];

%% Ranked table and exports
plotCfg.doRankedTable = true;
plotCfg.exportRankedTable = true;
plotCfg.exportCRLBQualityTable = true;

%% Saving figures
plotCfg.saveFigures = true;
plotCfg.outputDir = fullfile(pwd, "TemporalCorrelationPlots");

%% Run plotting pipeline
plotOutputs = PlotTemporalCorrelationDiagnostics( ...
    covOutputs, ...
    deGraafOutputs, ...
    plotCfg);
