function BerkeleyAOtumblingEThreshold(options)

%% TODO
% Will need to add code to compare thresholds for the three E positions, as
% well as figure generation/saving code of various stages of the
% computation so we can visualized what happened and make figures for a
% paper.
%

% Examples:
%{
    BerkeleyAOtumblingEThreshold( ...
        'fastParams', true, ...
        'validationThresholds',[0.0283]);

    BerkeleyAOtumblingEThreshold( ...
        'fastParams', true, ...
        'temporalModulationParams_numFrame', 6, ...
        'temporalModulationParams_xShiftPerFrame', [0 0 0 0 0 0], ...
        'temporalModulationParams_yShiftPerFrame', [0 0 0 0 0 0], ...
        'temporalModulationParams_backgroundRGBPerFrame', [1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0], ...
        'temporalModulationParams_stimOnFrames', [0 1 1 1 0 0], ...
        'validationThresholds',[]);
%}

%% Pick up optional arguments
%
% A number of these get passed into t_BerkeleyAOtumblingSceneEngine.
arguments
    % Run with fast parameters overrides
    options.fastParams (1,1) logical = true;

    % Keep rng doing the same thing each time for validation
    options.rngSeed (1,1) double = 12;

    % Print out/plot  more diagnostics, or not
    options.verbose (1,1) logical = false;
    options.visualEsOnMosaic (1,1) logical = false;
    options.visualizeScene (1,1) logical = true;

    % Wavelength support
    options.wave (:,1) double = (500:5:870)';

    % Psychometric parameters
    options.letterSizesNumExamined = 9;
    options.nTest = 512;
    options.thresholdP =  0.781;

    % Optics parameters
    options.pupilDiameterMm (1,1) double = 6;
    options.accommodatedWl  (1,1) double = 840;
    options.defocusDiopters (1,1) double =  0.05;

    % Use cone contrast
    options.useConeContrast (1,1) logical = false

     % Choose noise model
    %   Choices: 'Poisson'
    %            'Gaussian'
    options.whichNoisyInstanceNre (1,:) char = 'Poisson'
    options.gaussianSigma double = [];

     % Apply temporal filter?
    %
    % The timebase of the filter is assumed to match the frame rate, so we
    % only need to specify a list of filter values.  Since these can
    % describe gain changes as well as temporal processing per se, we do
    % not require that these sum to 1.  If you want to preserve signal
    % magnitude, you can normalize the filter values yourself, so that they
    % sum to 1. This can also be set to some string, e.g.,
    % 'photocurrentImpulseResponseBased', in which case the filter values
    % are computed on the fly
    %
    % Finally, this can be 'watsonFilter'
    options.temporalFilterValues (1,:) = []

    % Choose classifier engine
    %    rcePoisson - signal known exactly Poission max likelihood
    %    rceTemplateDistance - signal known exactly nearest L2 template
    %                 distance.
    options.whichClassifierEngine (1,:) char = 'rcePoisson'

    % AO scene parameters
    options.displayNPixels (1,1) double = 128;
    options.displayFOVDeg (1,1) double = 1.413;
    options.AOPrimaryWls (1,3) double = [840 683 543]; % [700 683 54];
    options.AOPrimaryFWHM (1,3) double = [22 27 23];
    options.AOCornealPowersUW (1,3) double = [141.4 10 10];
    options.plotDisplayCharacteristics (1,1) logical = false;
    options.chromaSpecification_type (1,:) char = 'RGBsettings';
    options.chromaSpecification_backgroundRGB (1,3) double = [1 0 0];
    options.chromaSpecification_foregroundRGB (1,3) double = [0 0 0];
    options.eHeightMin (1,1) double = 30;
    options.temporalModulationParams_frameRateHz (1,1) double = 60;
    options.temporalModulationParams_numFrame (1,1) double = 3;
    options.temporalModulationParams_xShiftPerFrame (1,:) double = [0 10/60 0];
    options.temporalModulationParams_yShiftPerFrame (1,:) double = [0 0 10/60];
    options.temporalModulationParams_backgroundRGBPerFrame (:,:) double = [0 0 0; 1 0 0; 0 0 0];
    options.temporalModulationParams_stimOnFrames (:,:) double = [0 1 0];

    % Run the validation check?  This gets overridden to empty if other
    % options change the conditions so that the validation data don't
    % apply.
    options.validationThresholds (1,:) double = []

    % These options do not get passed to t_BerkeleyAOtumblingETutorial
    options.writeFigures (1,1) logical = true;
end

%% Initialize
close all;

%% Fast parameter overrides
if (options.fastParams)
    options.displayNPixels = 128;
    options.letterSizesNumExamined = 3;
    options.nTest = 64;
end

% Make sure figures and results directories exist so that output writes
% don't fail
if (options.writeFigures)
    rootPath = ISETBerkeleyAOTumblingERootPath;
    if (~exist(fullfile(rootPath,'local','figures'),'dir'))
        mkdir(fullfile(rootPath,'local','figures'));
    end
    
    % Make sure figures and results directories exist so that output writes
    % don't fail
    rootPath = ISETBerkeleyAOTumblingERootPath;
    if (~exist(fullfile(rootPath,'local','figures'),'dir'))
        mkdir(fullfile(rootPath,'local','figures'));
    end
    if (~exist(fullfile(rootPath,'local','results'),'dir'))
        mkdir(fullfile(rootPath,'local','results'));
    end

    outputFiguresDir = fullfile(rootPath, 'local', 'figures');
    outputResultsDir = fullfile(rootPath, 'local', 'results');
end

% Set up summary filename and output dir
summaryFileName = 'Foo';
%summaryFileName = sprintf('Summary_%dms.mat', round(1000*integrationTime));
outputResultsDir = fullfile(ISETBerkeleyAOTumblingERootPath,'local','results',strrep(summaryFileName, '.mat',''));
outputFiguresDir =  fullfile(ISETBerkeleyAOTumblingERootPath,'local','figures',strrep(summaryFileName, '.mat',''));
if (~exist(outputResultsDir,'dir'))
    mkdir(outputResultsDir);
end
if (~exist(outputFiguresDir,'dir'))
    mkdir(outputFiguresDir);
end

%% Do all the hard work in the CSF generator tutorial function
%
% Add all parameters to an options struct and call
fn = fieldnames(options);
for i = 1:numel(fn)
    tutorialOptions.(fn{i}) = options.(fn{i});
end
optionsTemp = options;
optionsTemp = rmfield(optionsTemp,'writeFigures');
tutorialOptionsCell = [fieldnames(optionsTemp) , struct2cell(optionsTemp)]';
[logThreshold, logMAR, questObj, psychometricFunction, fittedPsychometricParams, ...
    trialByTrialStimulusAlternatives,trialByTrialPerformance] = ...
    t_BerkeleyAOtumblingEThreshold(tutorialOptionsCell{:});

% Print the threshold estimate
fprintf('Current threshold estimate: %g\n', 10 ^ logThreshold);

% TrailByTrial data template
keys = trialByTrialStimulusAlternatives.keys;
fprintf('trialByTrialStimulusAlternatives contents:\n');
for i = 1:length(keys)
    trialByTrialStimulusAlternatives(keys{i});
end

% Plot the derived psychometric function and other things.  The lower
% level routines put this in ISETBioJandJRootPath/figures.
% pdfFileName = sprintf('Performance_Reps_%d.pdf', options.nTest);
% pdfFileName = sprintf('%s_%s_%d_%d.pdf', responseFlag, options.exportCondition, ...
%         thresholdParameters.logThreshLimitHigh, thresholdParameters.logThreshLimitLow);
% plotDerivedPsychometricFunction(questObj, threshold, fittedPsychometricParams, ... 
%     thresholdParameters, fullfile(outputFiguresDir,pdfFileName), ...
%     'xRange', [10.^-thresholdParameters.logThreshLimitLow  10.^-thresholdParameters.logThreshLimitHigh]);
% if (options.visualEsOnMosaic)
%     pdfFileName = sprintf('Simulation_Reps_%d.pdf', options.nTest);
%     visualizeSimulationResults(questObj, threshold, fittedPsychometricParams, ...
%         thresholdParameters, tumblingEsceneEngines, theNeuralEngine, ...
%         fullfile(outputFiguresDir,pdfFileName));
% end
% 
% % Export the results
% exportFileName = sprintf('Results_Reps_%d.mat', options.nTest);
% fprintf('Saving data to %s\n', fullfile(outputResultsDir,exportFileName));
% exportSimulation(questObj, threshold, fittedPsychometricParams, ...
%     thresholdParameters, classifierPara, questEnginePara, ...
%     tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
%     fullfile(outputResultsDir,exportFileName));

% Save summary,  This allows examination of the numbers and/or
% replotting.
% save(fullfile(outputResultsDir,summaryFileName),"examinedPSFDataFiles","threshold","logMAR","LCA","TCA","theConeMosaic");
end



