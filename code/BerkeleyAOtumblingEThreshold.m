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

%% Pick up optional arguments
%
% A number of these get passed into t_BerkeleyAOtumblingSceneEngine.
arguments
    options.visualizeScene (1,1) logical = false;          % passing this into t_berkeleyAOtumblingESceneEngine
    options.displayNPixels (1,1) double = 512*0.2;         % Scale down because thresholds are small relative to instrument field size
    options.displayFOVDeg (1,1) double = 1.413*0.2;        % Scale down because thresholds are small relative to instrument field size
    options.wave (:,1) double = (500:5:870)';
    options.AOPrimaryWls (1,3) double = [840 683 543]; % [700 683 54];
    options.AOPrimaryFWHM (1,3) double = [22 27 23];
    options.AOCornealPowersUW (1,3) double = [141.4 10 10];
    options.ambientSpd (:,1) double = zeros(size((500:5:870)'));  % Adjust this if wave changes
    options.pupilSizeMM (1,1) double = 6;
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
    options.responseFlag (1,:) char = 'excitation'; % 'excitation' or 'photocurrent'
end

% Define the AO scene parameters for the experiment we are modeling
% Imaging (840 nm) power in uW.

aoSceneParams = struct(...
    'defocusDiopters', 0.05, ...
    'pupilDiameterMm', 6, ...
    'accommodatedWl', 840, ...
    'angleList', [0 90 180 270], ...
    'testing', false, ...
    'write', true, ...
    'verbose', false ...
);

% Initialize
close all;

% Make sure figures and results directories exist so that output writes
% don't fail
if(aoSceneParams.write)
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

% Get the tumbling E scene engines.
%
% At the moment cannot vary the set of orientations, but would be easy
% enough to allow this with a key/value pair passed to the called tutorial
% function.
%
% The scene engine tutorial returns its parameters, which are used below to
% try to match things up as best as possible.

% get scene engine arguments passed down into the tumbline E
% scene generation tutorial function.
orientations = aoSceneParams.angleList;
optionsCell = [fieldnames(options), struct2cell(options)]';
optionsCell = optionsCell(:)';

[sce0,sce90,sce180,sce270,backgroundSceneEngine,sceneParams] = t_BerkeleyAOtumblingESceneEngine(optionsCell{:});
tumblingEsceneEngines = {sce0, sce90, sce180, sce270};
clear sce0 sce90 sce180 sce270

% Parameters. These control many aspects of what gets done, particularly the subject.
%
% To run a different subject or pupil size, change 'psdDataSubdir'
% field below to have the desired subject number in the directory string, and the desired
% pupil string in the name if using other than the default 4 mm.
%
% Parameter fields below allow change of integration time, lens age,
% MPD, L:M:S proportion (although you won't get many S because of the
% tritanopic zone).
%
% Note that the custom pupil diameter field does not affect the
% optical quality because we are using the PSF read from the PSF data
% file.  What it does is allow you to set a pupil diameter different
% from the one for which the PSF was computed.  What this does is allow
% independent control of retinal illuminance and PSF, if you want to
% separate the two effects.  The full data set does contain PSFs
% computed for different pupil sizes for each TCA/LCA combination which
% may be used to explore the effect of pupil size on optical quality.
% Note again that using a PSF computed for various pupil sizes will not
% affect the retinal illuminance as that is controlled by the pupil
% size set explicitly here.  The PSF file naming convention for the
% various pupil sizes is described in the README.
%
% The code in this script is reasonably clever about creating figure
% and results subdirectories to hold its ouput, that keep the separate
% conditions you might run separate.  But it may not be perfect at
% this, particularly if you start digging deeper into the code and
% customizing more things.
params = struct(...
    'letterSizesNumExamined',  9, ...                           % How many sizes to use for sampling the psychometric curve (9 used in the paper)
    'maxLetterSizeDegs', 1, ...                                 % The maximum letter size in degrees of visual angle
    'sceneUpSampleFactor', 4, ...                               % Upsample scene, so that the pixel for the smallest scene is < cone aperture
    'mosaicIntegrationTimeSeconds', 1/options.temporalModulationParams_frameRateHz, ... % Integration time, matched to frame rate
    'nTest', 512, ...                                           % Number of trial to use for computing Pcorrect
    'thresholdP', 0.781, ...                                    % Probability correct level for estimating threshold performance
    'customLensAgeYears', [], ...                               % Lens age in years (valid range: 20-80), or empty to use the default age of 32.
    'customMacularPigmentDensity', [], ...                      % Cu∂stom MPD, or empty to use the default density of 0.35; example, 0.7
    'customConeDensities', [], ...                              % Custom L-M-S ratio or empty to use default; example [0.6 0.3 0.1]
    'customPupilDiameterMM', [], ...                            % Custom pupil diameter in MM or empty to use the value from the psfDataFile
    'visualizedPSFwavelengths', [], ...                         % Vector with wavelengths for visualizing the PSF. If set to empty[] there is no visualization; example 400:20:700
    'visualizeDisplayCharacteristics', false, ...     % Flag, indicating whether to visualize the display characteristics
    'visualizeScene', false, ...            % Flag, indicating whether to visualize one of the scenes
    'visualEsOnMosaic', false, ...   % Flag, indicating whether to visualize E's against mosaic as function of their size
    'outputResultsDir', outputResultsDir, ...
    'outputFiguresDir', outputFiguresDir ...                   % directory for saving output figures
    );

% Check on size
if (params.maxLetterSizeDegs ~= 1)
    error('At present, the underlying code only allows a max parameter value of 1.');
end

% Set up summary filename and output dir
summaryFileName = sprintf('Summary_%dms.mat', round(1000*params.mosaicIntegrationTimeSeconds));
% if (~isempty(params.customMacularPigmentDensity))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_MPD_%2.2f.mat', params.customMacularPigmentDensity));
% end
% if (~isempty(params.customPupilDiameterMM))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_pupilDiamMM_%2.2f.mat', params.customPupilDiameterMM));
% end
% if (~isempty(params.customConeDensities))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_cones_%2.2f_%2.2f_%2.2f.mat', params.customConeDensities(1), params.customConeDensities(2), params.customConeDensities(3)));
% end
% if (~isempty(params.customLensAgeYears))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_lensAge_%d.mat', params.customLensAgeYears));
% end
params.outputResultsDir = fullfile(ISETBerkeleyAOTumblingERootPath,'local','results',strrep(summaryFileName, '.mat',''));
params.outputFiguresDir =  fullfile(ISETBerkeleyAOTumblingERootPath,'local','figures',strrep(summaryFileName, '.mat',''));
if (~exist(params.outputResultsDir,'dir'))
    mkdir(params.outputResultsDir);
end
if (~exist(params.outputFiguresDir,'dir'))
    mkdir(params.outputFiguresDir);
end

% Unpack simulation params
letterSizesNumExamined = params.letterSizesNumExamined;
maxLetterSizeDegs = params.maxLetterSizeDegs;
mosaicIntegrationTimeSeconds = params.mosaicIntegrationTimeSeconds;
nTest = params.nTest;
thresholdP = params.thresholdP;

%% Create neural response engine for photopigment Excitations

% add all parameters to options struct
fn = fieldnames(aoSceneParams);
for i = 1:numel(fn)
    options.(fn{i}) = aoSceneParams.(fn{i});
end

responseFlag = options.responseFlag;    % 'excitation' or 'photocurrent'. Indicating whether to create neural response engine for cone excitation or photocurrent 
theNeuralEngine = createNeuralResponseEngine(responseFlag, options);

% Poisson n-way AFC
classifierEngine = responseClassifierEngineNWay(@rcePoisson);

% Parameters associated with use of the Poisson classifier.
classifierPara = struct('trainFlag', 'none', ...
    'testFlag', 'random', ...
    'nTrain', 1, 'nTest', nTest);

%% Parameters for threshold estimation/quest engine
thresholdParameters = struct(...
    'maxParamValue', maxLetterSizeDegs, ...    % The maximum value of the examined param (letter size in degs).  At present need this to be 1.
    'logThreshLimitLow', 3, ...                % minimum log10(normalized param value).  The convention is that this means 10^-3.
    'logThreshLimitHigh', 1, ...               % maximum log10(normalized param value).  The convention is that this means 10^-1.
    'logThreshLimitDelta', 0.01, ...
    'slopeRangeLow', 1/20, ...
    'slopeRangeHigh', 500/20, ...
    'slopeDelta', 2/20, ...
    'thresholdCriterion', thresholdP, ...
    'guessRate', 1/numel(orientations), ...
    'lapseRate', [0 0.02]);

% Parameters for Quest
questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', nTest*letterSizesNumExamined, ...
    'maxTrial', nTest*letterSizesNumExamined, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

% Compute psychometric function for the 4AFC paradigm with the 4 E scenes
[threshold, questObj, psychometricFunction, fittedPsychometricParams] = computeThreshold(...
    tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdParameters, questEnginePara, ...
    'TAFC', false, ...
    'visualizeAllComponents', ~true, ...
    'beVerbose', true, ...
    'parameterIsContrast',false);

% Get the threshold estimate
fprintf('Current threshold estimate: %g\n', 10 ^ threshold);

% Plot the derived psychometric function and other things.  The lower
% level routines put this in ISETBioJandJRootPath/figures.
pdfFileName = sprintf('Performance_Reps_%d.pdf', nTest);
plotDerivedPsychometricFunction(questObj, threshold, fittedPsychometricParams, ... 
    thresholdParameters, fullfile(params.outputFiguresDir,pdfFileName), ...
    'xRange', [10.^-thresholdParameters.logThreshLimitLow  10.^-thresholdParameters.logThreshLimitHigh]);
if (params.visualEsOnMosaic)
    pdfFileName = sprintf('Simulation_Reps_%d.pdf', nTest);
    visualizeSimulationResults(questObj, threshold, fittedPsychometricParams, ...
        thresholdParameters, tumblingEsceneEngines, theNeuralEngine, ...
        fullfile(params.outputFiguresDir,pdfFileName));
end

% Export the results
exportFileName = sprintf('Results_Reps_%d.mat', nTest);
% if (~isempty(params.customMacularPigmentDensity))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_MPD_%2.2f.mat', params.customMacularPigmentDensity));
% end
% if (~isempty(params.customPupilDiameterMM))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_PupilDiamMM_%2.2f.mat', params.customPupilDiameterMM));
% end
% if (~isempty(params.customConeDensities))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_cones_%2.2f_%2.2f_%2.2f.mat', params.customConeDensities(1), params.customConeDensities(2), params.customConeDensities(3)));
% end
% if (~isempty(params.customLensAgeYears))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_lensAge_%d.mat', params.customLensAgeYears));
% end
% 
fprintf('Saving data to %s\n', fullfile(params.outputResultsDir,exportFileName));
exportSimulation(questObj, threshold, fittedPsychometricParams, ...
    thresholdParameters, classifierPara, questEnginePara, ...
    tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
    fullfile(params.outputResultsDir,exportFileName));

% logMAR(iPSF) = log10(threshold(iPSF)*60/5);

% Save summary,  This allows examination of the numbers and/or
% replotting.
% save(fullfile(params.outputResultsDir,summaryFileName),"examinedPSFDataFiles","threshold","logMAR","LCA","TCA","theConeMosaic");


end

function theNeuralEngine = createNeuralResponseEngine(responseType, paramStruct)
    % Validate responseType 
    if ~ischar(responseType)
        error('responseType must be a string.');
    end
    
    if strcmp(responseType, 'excitation')
        % This calculates isomerizations in a patch of cone mosaic with Poisson
        % noise, and includes optical blur.
        neuralParams = nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic;
        
        % Set optics params
        wls = paramStruct.wave;
        fieldSizeDegs = paramStruct.displayFOVDeg;
        integrationTime = 1/paramStruct.temporalModulationParams_frameRateHz;
        accommodatedWl = paramStruct.accommodatedWl;
        pupilDiameterMm = paramStruct.pupilDiameterMm;
        defocusDiopters = paramStruct.defocusDiopters;
        
        neuralParams.opticsParams.wls = wls;
        neuralParams.opticsParams.pupilDiameterMM = pupilDiameterMm;
        neuralParams.opticsParams.defocusAmount = defocusDiopters;
        neuralParams.opticsParams.accommodatedWl = accommodatedWl;
        neuralParams.opticsParams.zCoeffs = zeros(66,1);
        neuralParams.opticsParams.defeatLCA = true;
        neuralParams.verbose = paramStruct.verbose;
        
        % Cone params
        neuralParams.coneMosaicParams.wave = wls;
        neuralParams.coneMosaicParams.fovDegs = fieldSizeDegs;
        neuralParams.coneMosaicParams.timeIntegrationSeconds = integrationTime;
    
        % Create the neural response engine
        theNeuralEngine = neuralResponseEngine(@nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic, neuralParams);
    
    elseif strcmp(responseType, 'photocurrent')
        % This calculates photocurrent in a patch of cone mosaic with Poisson
        % noise, and includes optical blur.
        neuralParams = nreAOPhotocurrentWithNoEyeMovementsCMosaic;
        
        % Set optics params
        wls = paramStruct.wave;
        fieldSizeDegs = paramStruct.displayFOVDeg;
        integrationTime = 1/paramStruct.temporalModulationParams_frameRateHz;
        accommodatedWl = paramStruct.accommodatedWl;
        pupilDiameterMm = paramStruct.pupilDiameterMm;
        defocusDiopters = paramStruct.defocusDiopters;
        
        neuralParams.opticsParams.wls = wls;
        neuralParams.opticsParams.pupilDiameterMM = pupilDiameterMm;
        neuralParams.opticsParams.defocusAmount = defocusDiopters;
        neuralParams.opticsParams.accommodatedWl = accommodatedWl;
        neuralParams.opticsParams.zCoeffs = zeros(66,1);
        neuralParams.opticsParams.defeatLCA = true;
        neuralParams.verbose = paramStruct.verbose;
        
        % Cone params
        neuralParams.coneMosaicParams.wave = wls;
        neuralParams.coneMosaicParams.fovDegs = fieldSizeDegs;
        neuralParams.coneMosaicParams.timeIntegrationSeconds = integrationTime; % 1/60

        % Create the neural response engine
        theNeuralEngine = neuralResponseEngine(@nreAOPhotocurrentWithNoEyeMovementsCMosaic, neuralParams); 
    else
        % If responseFlag is neither 'excitation' nor 'photocurrent'
        error('Invalid response type: %s', responseType);
    end
end

