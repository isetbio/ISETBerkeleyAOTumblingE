function [logThreshold, logMAR, questObj, psychometricFunction, fittedPsychometricParams, ...
    trialByTrialStimulusAlternatives,trialByTrialPerformance,theNeuralEngine] = BerkeleyAOTumblingEThreshold(options)

% Examples:
%{
    fileSuffix = 'FastExample';
    BerkeleyAOTumblingEThreshold( ...
        'fastParams', true, ...
        'visualizeScene', true, ...
        'visualizeEsOnMosaic', true, ...
        'visualizeEsWhichFrames', 2, ... 
        'temporalModulationParams_backgroundRGBPerFrame', [0 0 0; 1 0 0; 0 0 0], ...
        'fileSuffix', fileSuffix, ...
        'validationThresholds',[0.028]);

    BerkeleyAOTumblingEThreshold( ...
        'fastParams', false, ...
        'temporalModulationParams_numFrame', 6, ...
        'temporalModulationParams_xShiftPerFrameMin', { [0 0 0 0 0 0] }, ...
        'temporalModulationParams_yShiftPerFrameMin', { [0 0 0 0 0 0] }, ...
        'temporalModulationParams_backgroundRGBPerFrame', [1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0], ...
        'temporalModulationParams_stimOnFrames', [0 1 1 1 0 0], ...
        'temporalModulationParams_frameRateHz', 30, ...
        'temporalFilterValues', 'photocurrentImpulseResponseBased', ...
        'minLetterSizeMinutes', 10, ...
        'maxLetterSizeMinutes', 10, ...
        'letterSizesNumExamined', 1, ...
        'nTest', 512, ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 200, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'validationThresholds',[]);

     BerkeleyAOTumblingEThreshold( ...
        'fastParams', false, ...
        'temporalModulationParams_numFrame', 6, ...
        'temporalModulationParams_xShiftPerFrameMin',  { [0 0 0  0  0 0] }, ...
        'temporalModulationParams_yShiftPerFrameMin',  { [0 0 2  4  0 0] }, ...
        'temporalModulationParams_backgroundRGBPerFrame', [1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0], ...
        'temporalModulationParams_stimOnFrames', [0 1 1 1 0 0], ...
        'temporalModulationParams_frameRateHz', 30, ...
        'temporalFilterValues', 'photocurrentImpulseResponseBased', ...
        'minLetterSizeMinutes', 10, ...
        'maxLetterSizeMinutes', 10, ...
        'letterSizesNumExamined', 1, ...
        'nTest', 512, ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 200, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'validationThresholds',[]);

    BerkeleyAOTumblingEThreshold( ...
        'fastParams', false, ...
        'temporalModulationParams_numFrame', 6, ...
        'temporalModulationParams_xShiftPerFrameMin',  { [0 0 0  0  0 0] }, ...
        'temporalModulationParams_yShiftPerFrameMin',  { [0 0 4  6  0 0] }, ...
        'temporalModulationParams_backgroundRGBPerFrame', [1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0], ...
        'temporalModulationParams_stimOnFrames', [0 1 1 1 0 0], ...
        'temporalModulationParams_frameRateHz', 30, ...
        'temporalFilterValues', 'photocurrentImpulseResponseBased', ...
        'minLetterSizeMinutes', 10, ...
        'maxLetterSizeMinutes', 10, ...
        'letterSizesNumExamined', 1, ...
        'nTest', 512, ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 200, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'validationThresholds',[]);
%}

%% Pick up optional arguments
%
% A number of these get passed into t_BerkeleyAOtumblingSceneEngine.
arguments
    % Run with fast parameters overrides
    options.fastParams (1,1) logical = true;

    % Keep rng doing the same thing each time for validation.
    % 0 means don't touch the seed so each run is random.
    options.rngSeed (1,1) double = 12;

    % Print out/plot  more diagnostics, or not
    options.verbose (1,1) logical = false;
    options.visualizeEsOnMosaic (1,1) logical = false;
    options.visualizeEsWhichFrames (1,:) double = 1;
    options.visualizeScene (1,1) logical = true;
    options.visualizeEachResponse (1,1) logical = false;
    options.responseVisualizationFunction = [];
    options.responseVideoFileName(1,:) char = '';
    options.maxVisualizedNoisyResponseInstances = 1;
    options.maxVisualizedNoisyResponseInstanceStimuli = 1;

    % Wavelength support
    options.wave (:,1) double = (500:5:870)';

    % Psychometric parameters
    options.letterSizesNumExamined = 9;
    options.nTest = 512;
    options.thresholdP =  0.781;
    options.minLetterSizeMinutes (1,1) double = 0.02;
    options.maxLetterSizeMinutes (1,1) double = 2;

    % Optics parameters
    options.pupilDiameterMm (1,1) double = 6;
    options.accommodatedWl  (1,1) double = 840;
    options.defocusDiopters (1,1) double = 0.05;

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
    options.watsonParams_tau = 6.25;

    % Choose classifier engine
    %    rcePoisson - signal known exactly Poission max likelihood
    %    rceTemplateDistance - signal known exactly nearest L2 template
    %                 distance.
    options.whichClassifierEngine (1,:) char = 'rcePoisson'

    % AO scene parameters
    options.displayNPixels (1,1) double = 128;
    options.displayFOVDeg (1,1) double = 1.413;
    options.eccDegs (1,2) double = [0 0];
    options.AOPrimaryWls (1,3) double = [840 683 543]; % [700 683 54];
    options.AOPrimaryFWHM (1,3) double = [22 27 23];
    options.AOCornealPowersUW (1,3) double = [141.4 10 10];
    options.plotDisplayCharacteristics (1,1) logical = false;
    options.chromaSpecification_type (1,:) char = 'RGBsettings';
    options.chromaSpecification_backgroundRGB (1,3) double = [1 0 0];
    options.chromaSpecification_foregroundRGB (1,3) double = [0 0 0];
    options.temporalModulationParams_frameRateHz (1,1) double = 60;
    options.temporalModulationParams_numFrame (1,1) double = 3;
    options.jitterMinutesX (1,1) double = 0;
    options.jitterMinutesY (1,1) double = 0;
    options.temporalModulationParams_xShiftPerFrameMin (:,:) cell = {[0 10/60 0]};
    options.temporalModulationParams_yShiftPerFrameMin (:,:) cell = {[0 0 10/60]};
    options.temporalModulationParams_backgroundRGBPerFrame (:,:) double = [0 0 0; 1 0 0; 0 0 0];
    options.temporalModulationParams_stimOnFrames (:,:) double = [0 1 0];

    % Run the validation check?  This gets overridden to empty if other
    % options change the conditions so that the validation data don't
    % apply.
    options.validationThresholds (1,:) double = []

    % These options do not get passed to t_BerkeleyAOTumblingETutorial
    options.writeFigures (1,1) logical = true;
    options.writeSummary (1,1) logical = true;
    options.rootPath (1,:) char = '';
    options.fileSuffix (1,:) char = 'Example';
end

%% Initialize
close all;

%% Fast parameter overrides
if (options.fastParams)
    options.displayNPixels = 128;
    options.letterSizesNumExamined = 3;
    options.nTest = 64;
end

%% Set root path unless passed explicitly
if (isempty(options.rootPath))
    options.rootPath = getpref('ISETBerkeleyAOTumblingE','dataPath');
end

%% Set up summary filename and output dir
%
% Saving these with an options string that the user can set to denote the
% condtions.
summaryFileName = ['BerkeleyAOTumblingEThreshold_' options.fileSuffix];
options.outputResultsDir = fullfile(options.rootPath,'results',summaryFileName);
options.outputFiguresDir =  fullfile(options.rootPath,'figures',summaryFileName);
if (~exist(options.outputResultsDir,'dir'))
    mkdir(options.outputResultsDir);
end
if (~exist(options.outputFiguresDir,'dir'))
    mkdir(options.outputFiguresDir);
end
options.scenePdfFileBase = fullfile('Scene');
options.visualizeEsFileBase = fullfile('EsOnMosaic');

%% Do all the hard work in the CSF generator tutorial function
%
% Add all parameters to an options struct and call
fn = fieldnames(options);
for i = 1:numel(fn)
    tutorialOptions.(fn{i}) = options.(fn{i});
end
optionsTemp = options;
optionsTemp = rmfield(optionsTemp,'rootPath')';
optionsTemp = rmfield(optionsTemp,'writeFigures');
optionsTemp = rmfield(optionsTemp,'writeSummary');
optionsTemp = rmfield(optionsTemp,'fileSuffix');

% We'll make these plots here if we want them, so pass false to the
% lower level call.
optionsTemp.plotPsychometric = false;

% Put structure in right form for arguments block
tutorialOptionsCell = [fieldnames(optionsTemp) , struct2cell(optionsTemp)]';

% Do the hard work
[logThreshold, logMAR, questObj, psychometricFunction, fittedPsychometricParams, thresholdPara, ...
    trialByTrialStimulusAlternatives, trialByTrialPerformance, trialByTrialWhichResponses,theNeuralEngine] = ...
    t_BerkeleyAOTumblingEThreshold(tutorialOptionsCell{:});
threshold = 10.^logThreshold;

%% Print the threshold estimate
fprintf('Current threshold estimate: %g\n', threshold);

% Plot the psychometric function
pdfFileName = fullfile(options.outputFiguresDir,sprintf('Performance_Reps_%d.pdf', options.nTest));
[stimulusLevels, pCorrect] = plotPsychometricFunction(questObj, threshold, fittedPsychometricParams, ...
        thresholdPara, pdfFileName, 'xRange', [options.minLetterSizeMinutes/60  options.maxLetterSizeMinutes/60]);

% Trial by trial data template
stimKeys = trialByTrialStimulusAlternatives.keys;
performanceKeys = trialByTrialPerformance.keys;
if (length(stimKeys) ~= length(performanceKeys))
    error('Should have same number of stimuli as perfomance');
end
whichResponseKeys = trialByTrialWhichResponses.keys;
if (length(stimKeys) ~= length(whichResponseKeys))
    error('Should have same number of stimuli as which responses');
end
for i = 1:length(stimKeys)
    tByTStimAlternatives = trialByTrialStimulusAlternatives(stimKeys{i});
    tByTPerformance = trialByTrialPerformance(performanceKeys{i});
    tByTResponseAlternatives = trialByTrialWhichResponses(whichResponseKeys{i});
    checkPerformance = double(tByTStimAlternatives == tByTResponseAlternatives);
    if (any(tByTPerformance ~= checkPerformance))
        error('Inconsistency in performance report');
    end
    if (pCorrect ~= sum(checkPerformance)/length(checkPerformance))
        error('Inconsistency in pCorrect');
    end

    % Get confusion matrix
    nAlternatives = length(unique(tByTStimAlternatives));
    pRespondAAWithStimBB = zeros(nAlternatives,nAlternatives);
    for aa = 1:nAlternatives
        for bb = 1:nAlternatives
            index = find(tByTStimAlternatives == bb);
            pRespondAAWithStimBB(aa,bb) = sum(tByTResponseAlternatives(index) == aa)/length(index);           
        end
    end

end

% Print out table of stimulus levels and pCorrect
fprintf('\nMeasured performance\n')
for ii = 1:length(stimulusLevels)
    fprintf('%0.2f min (%0.3f deg), %0.2f pCorrect\n',60*stimulusLevels(ii),stimulusLevels(ii),pCorrect(ii));
end
fprintf('\n');

% Save summary,  This allows examination of the numbers and/or
% replotting.
save(fullfile(options.outputResultsDir,[summaryFileName '.mat']),'-v7.3');

end



