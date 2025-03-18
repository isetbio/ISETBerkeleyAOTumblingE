function BerkeleyAOtumblingEThreshold(options)

%% TODO
% Will need to add code to compare thresholds for the three E positions, as
% well as figure generation/saving code of various stages of the
% computation so we can visualized what happened and make figures for a
% paper.
%
% Need to use AO optics in neural engine, and turn off chromatic
% aberration.  This will model what the AO system does in producing the
% retinal image.
%
% Need to add photocurrent filtering.

%% 10/14/24, QF, DHB.
%
% When the photocurrent is turned on (via responseFlag field
% in options below, it crashes. This is because some assertion fails.  So
% we need to track that down.
%
% QF reports that this script takes about 30 minutes to run for excitations
% in its default configuration.
%
% Figure out the photocurrent issue.

% Examples:
%{
    BerkeleyAOtumblingEThreshold( ...
    ...
    );
%}

%% Pick up optional arguments
%
% A number of these get passed into t_BerkeleyAOtumblingSceneEngine.
arguments
    % Run with fast parameters overrides
    options.fastParams (1,1) logical = true;

    % Print out/plot  more diagnostics, or not
    options.verbose (1,1) logical = false;
    options.visualEsOnMosaic (1,1) logical = false;
    options.visualizeScene (1,1) logical = false;                 

    % Wavelength support
    options.wave (:,1) double = (500:5:870)';

     % Psychometric parameters
    options.letterSizesNumExamined = 9;
    options.nTest = 512;
    options.thresholdP =  0.781;

    % Scene parameters
    options.displayNPixels (1,1) double = 128   ;                 
    options.displayFOVDeg (1,1) double = 1.413*0.25;       
    options.AOPrimaryWls (1,3) double = [840 683 543]; 
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

    % Optics parameters
    options.pupilDiameterMm (1,1) double = 6;
    options.accommodatedWl  (1,1) double = 840;
    options.defocusDiopters (1,1) double =  0.05;

    % Use cone contrast
    options.useConeContrast (1,1) logical = false

     % Choose noise model
    %   Choices: 'Poisson'
    %                  'Gaussian'
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
    options.temporalFilterValues (1,:) = []
    options.exportConditionLabel (1,:) char = 'no change'; 

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
summaryFileName = sprintf('Summary_%dms.mat', round(1000*integrationTime));
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
optionsTemp = rmfield(optionsTemp,'exportConditionLabel');
optionsTemp = rmfield(optionsTemp,'writeFigures');
tutorialOptionsCell = [fieldnames(optionsTemp) , struct2cell(optionsTemp)]';
[logThreshold, logMAR, questObj, psychometricFunction, fittedPsychometricParams, ...
    trialByTrialStimulusAlternatives,trialByTrialPerformance] = ...
    t_BerkeleyAOtumblingEThreshold(tutorialOptionsCell{:});

% Print the threshold estimate
fprintf('Current threshold estimate: %g\n', 10 ^ logThreshold);

% temporary solution for trailByTrial printing
% Create a containers.Map object (dictionary equivalent in MATLAB)
keys = trialByTrialStimulusAlternatives.keys;
fprintf('trialByTrialStimulusAlternatives contents:\n');
for i = 1:length(keys)
    trialByTrialStimulusAlternatives(keys{i})
end

% Plot the derived psychometric function and other things.  The lower
% level routines put this in ISETBioJandJRootPath/figures.
% pdfFileName = sprintf('Performance_Reps_%d.pdf', options.nTest);
pdfFileName = sprintf('%s_%s_%d_%d.pdf', responseFlag, options.exportCondition, ...
        thresholdParameters.logThreshLimitHigh, thresholdParameters.logThreshLimitLow);
plotDerivedPsychometricFunction(questObj, threshold, fittedPsychometricParams, ... 
    thresholdParameters, fullfile(outputFiguresDir,pdfFileName), ...
    'xRange', [10.^-thresholdParameters.logThreshLimitLow  10.^-thresholdParameters.logThreshLimitHigh]);
if (options.visualEsOnMosaic)
    pdfFileName = sprintf('Simulation_Reps_%d.pdf', options.nTest);
    visualizeSimulationResults(questObj, threshold, fittedPsychometricParams, ...
        thresholdParameters, tumblingEsceneEngines, theNeuralEngine, ...
        fullfile(outputFiguresDir,pdfFileName));
end

% Export the results
exportFileName = sprintf('Results_Reps_%d.mat', options.nTest);
fprintf('Saving data to %s\n', fullfile(outputResultsDir,exportFileName));
exportSimulation(questObj, threshold, fittedPsychometricParams, ...
    thresholdParameters, classifierPara, questEnginePara, ...
    tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
    fullfile(outputResultsDir,exportFileName));

% Save summary,  This allows examination of the numbers and/or
% replotting.
% save(fullfile(outputResultsDir,summaryFileName),"examinedPSFDataFiles","threshold","logMAR","LCA","TCA","theConeMosaic");
end

function theNeuralEngine = createNeuralResponseEngine(responseType, paramStruct)
    % Validate responseType 
    if ~ischar(responseType)
        error('responseType must be a string.');
    end
    
    if strcmp(responseType, 'excitation')
        % This calculates isomerizations in a patch of cone mosaic with Poisson
        % noise, and includes optical blur.
        nreNoiseFreeParams = nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic;
        
        % Set optics params
        wls = paramStruct.wave;
        fieldSizeDegs = paramStruct.displayFOVDeg;
        integrationTime = 1/paramStruct.temporalModulationParams_frameRateHz;
        accommodatedWl = paramStruct.accommodatedWl;
        pupilDiameterMm = paramStruct.pupilDiameterMm;
        defocusDiopters = paramStruct.defocusDiopters;
        
        nreNoiseFreeParams.opticsParams.wls = wls;
        nreNoiseFreeParams.opticsParams.pupilDiameterMM = pupilDiameterMm;
        nreNoiseFreeParams.opticsParams.defocusAmount = defocusDiopters;
        nreNoiseFreeParams.opticsParams.accommodatedWl = accommodatedWl;
        nreNoiseFreeParams.opticsParams.zCoeffs = zeros(66,1);
        nreNoiseFreeParams.opticsParams.defeatLCA = true;
        nreNoiseFreeParams.verbose = paramStruct.verbose;
        
        % Cone params
        nreNoiseFreeParams.coneMosaicParams.wave = wls;
        nreNoiseFreeParams.coneMosaicParams.fovDegs = fieldSizeDegs;
        nreNoiseFreeParams.coneMosaicParams.timeIntegrationSeconds = integrationTime;
    
        % Create the neural response engine
        theNeuralEngine = neuralResponseEngine(@nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic, nreNoiseFreeParams);
    
    elseif strcmp(responseType, 'photocurrent')
        % This calculates photocurrent in a patch of cone mosaic with Poisson
        % noise, and includes optical blur.
        nreNoiseFreeParams = nreAOPhotocurrentWithNoEyeMovementsCMosaic;
        
        % Set optics params
        wls = paramStruct.wave;
        fieldSizeDegs = paramStruct.displayFOVDeg;
        integrationTime = 1/paramStruct.temporalModulationParams_frameRateHz;
        accommodatedWl = paramStruct.accommodatedWl;
        pupilDiameterMm = paramStruct.pupilDiameterMm;
        defocusDiopters = paramStruct.defocusDiopters;
        
        nreNoiseFreeParams.opticsParams.wls = wls;
        nreNoiseFreeParams.opticsParams.pupilDiameterMM = pupilDiameterMm;
        nreNoiseFreeParams.opticsParams.defocusAmount = defocusDiopters;
        nreNoiseFreeParams.opticsParams.accommodatedWl = accommodatedWl;
        nreNoiseFreeParams.opticsParams.zCoeffs = zeros(66,1);
        nreNoiseFreeParams.opticsParams.defeatLCA = true;
        nreNoiseFreeParams.verbose = paramStruct.verbose;
        
        % Cone params
        nreNoiseFreeParams.coneMosaicParams.wave = wls;
        nreNoiseFreeParams.coneMosaicParams.fovDegs = fieldSizeDegs;
        nreNoiseFreeParams.coneMosaicParams.timeIntegrationSeconds = integrationTime; % 1/60

        % Create the neural response engine
        theNeuralEngine = neuralResponseEngine(@nreAOPhotocurrentWithNoEyeMovementsCMosaic, nreNoiseFreeParams); 
    else
        % If responseFlag is neither 'excitation' nor 'photocurrent'
        error('Invalid response type: %s', responseType);
    end
end

