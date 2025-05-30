% BerkeleyAOTumblingECalcs
%
% Code that models Tuten AO stabilized and not tumbling E experiment.
% This sets up options to match the experiment and runs the calculations.

%% Clear
clear; close all;

%% Name
calcName = 'Visualize6';

%% Visualization
visualizeScene = true;
visualizeEsOnMosaic = false;
visualizeEachResponse =  true;
responseVisualizationFunction = @nreVisualizeCMosaic;
maxVisualizedNoisyResponseInstances = 8;
maxVisualizedNoisyResponseInstanceStimuli = 4;

%% Parameters
%
% Letter size
letterSizeMinutes = 10;

% Eccentricity
eccDegs = [-1 0];

% Frames off, on , off in native experimental frame duration.
% These get expanded if we simulate faster.
baseOffFramesStart = 1;
baseOnFrames = 3;
baseOffFramesEnd = 6;

% Experimental and simulation frame rates.  The latter
% must be an integer multiple of the former.
expTemporalFrequencyHz = 30;
temporalFrequencyHz = 90;
if (rem(temporalFrequencyHz,expTemporalFrequencyHz) ~= 0)
    error('Temporal frequency must be an integer multiple of experimental temporal frequency');
end

% Compute simulation frame rate timing
frameMultiplier = temporalFrequencyHz/expTemporalFrequencyHz ;
offFramesStart = frameMultiplier*baseOffFramesStart;
onFrames = frameMultiplier*baseOnFrames;
offFramesEnd = frameMultiplier*baseOffFramesEnd;
totalFrames = offFramesStart + onFrames + offFramesEnd;

% When is stimulus on?
stimOnFrames = zeros(1,totalFrames);
stimOnFrames(offFramesStart+1:offFramesStart+onFrames) = ones(1,onFrames);

% Number of tests to simulate for each condition
nTest = 12;

% Background info
backgroundRGB = [1 0 0];
backgroundRGBPerFrame = backgroundRGB(ones(totalFrames,1),:);
foregroundRGB = [0 0 0];

% Set up y shift vectors for each of the three steps
% rawShiftsMinutes = [0 1 2 4];
rawShiftsMinutes = [0];
nShifts = length(rawShiftsMinutes);
for ss = 1:nShifts
    baseShiftMinutes = rawShiftsMinutes(ss);
    shiftIndex = 1;
    for jj = 0:baseOnFrames-1
        for ii = 1:frameMultiplier
            theShiftOn{ss}(shiftIndex) = jj*baseShiftMinutes;
            shiftIndex = shiftIndex + 1;
        end
    end
    theShift{ss} = zeros(1,totalFrames);
    theShift{ss}(offFramesStart+1:offFramesStart+onFrames) = theShiftOn{ss};
end

%% Calculations for each filter model and shift, positive y direction shifts
nReplications = 1;
filterModels = {[], 'photocurrentImpulseResponseBased', 'watsonFilter'};
nFilterModels = length(filterModels);
watsonParams_tau = 12;
noiseSds = [0.4 0.4 0.4];
jitterRangeMinutes = 2;
for ss = 1:nShifts
    for rr = 1:nReplications
        jitterMinutesX = jitterRangeMinutes*(rand-0.5);
        jitterMinutesY = jitterRangeMinutes*(rand-0.5);
        for ff = 1:nFilterModels
            fileSuffix = sprintf('%s_posYShift_%d_Rep_%d_filter_%d',calcName,ss,rr,ff);
            fprintf('%s\n',fileSuffix);
            BerkeleyAOTumblingEThreshold( ...
                'fastParams', false, ...
                'rngSeed', 0, ...
                'eccDegs', eccDegs, ...
                'chromaSpecification_backgroundRGB', [1 0 0], ...
                'chromaSpecification_foregroundRGB', [0 0 0], ...
                'temporalModulationParams_numFrame', totalFrames, ...
                'temporalModulationParams_xShiftPerFrameMin', jitterMinutesX*ones(1,totalFrames), ...
                'temporalModulationParams_yShiftPerFrameMin', theShift{ss}+jitterMinutesY, ...
                'temporalModulationParams_backgroundRGBPerFrame', backgroundRGBPerFrame, ...
                'temporalModulationParams_stimOnFrames', stimOnFrames, ...
                'temporalModulationParams_frameRateHz', temporalFrequencyHz , ...
                'temporalFilterValues', filterModels{ff}, ...
                'watsonParams_tau', watsonParams_tau, ...
                'minLetterSizeMinutes', letterSizeMinutes , ...
                'maxLetterSizeMinutes', letterSizeMinutes , ...
                'letterSizesNumExamined', 1, ...
                'nTest', nTest, ...
                'useConeContrast', true, ...
                'whichNoisyInstanceNre', 'Gaussian', ...
                'gaussianSigma', noiseSds(ff), ...
                'whichClassifierEngine', 'rceTemplateDistance', ...
                'visualizeScene', visualizeScene, ...
                'visualizeEsOnMosaic', visualizeEsOnMosaic, ...
                'visualizeEsWhichFrames', offFramesStart+1, ...
                'visualizeEachResponse', visualizeEachResponse, ...
                'responseVisualizationFunction', responseVisualizationFunction, ...
                'responseVideoFileName', fileSuffix, ...
                'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances, ...
                'maxVisualizedNoisyResponseInstanceStimuli',maxVisualizedNoisyResponseInstanceStimuli, ...
                'fileSuffix', fileSuffix, ...
                'validationThresholds',[]);
        end
    end
end

