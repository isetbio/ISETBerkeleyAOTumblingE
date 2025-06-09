% BerkeleyAOTumblingECalcsIntermixed
%
% Code that models Tuten AO stabilized and not tumbling E experiment.
% This sets up options to match the experiment and runs the calculations.

%% Clear
clear; close all;

%% Name
calcName = 'Calcs6';

%% Visualization
visualizeScene = false;
visualizeEsOnMosaic = false;

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
%
% Must be evenly dividable by number of conditions, here 64
nTest = 12800;

% Background info
backgroundRGB = [1 0 0];
backgroundRGBPerFrame = backgroundRGB(ones(totalFrames,1),:);
foregroundRGB = [0 0 0];

% Set up shift vectors for each of the three steps
rawShiftsMinutes = [0 1 2 4];
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
    theShifts{ss} = zeros(1,totalFrames);
    theShifts{ss}(offFramesStart+1:offFramesStart+onFrames) = theShiftOn{ss};
end

%% Shift directions
shiftDirectionNames = {'posYShift' 'negYShift' 'posXShift' 'negXShift'};
nDirections = length(shiftDirectionNames);
for dd = 1:nDirections
    for ss = 1:nShifts
        switch shiftDirectionNames{dd}
            case 'posYShift'
                theShiftsX{dd,ss} = zeros(1,totalFrames);
                theShiftsY{dd,ss} = theShifts{ss};
            case 'negYShift'
                theShiftsX{dd,ss} = zeros(1,totalFrames);
                theShiftsY{dd,ss} = -theShifts{ss};
            case 'posXShift'
                theShiftsX{dd,ss} = theShifts{ss};
                theShiftsY{dd,ss} = zeros(1,totalFrames);
            case 'negXShift'
                theShiftsX{dd,ss} = -theShifts{ss};
                theShiftsY{dd,ss} = zeros(1,totalFrames);
            otherwise
                error('Unknown shift direction name');
        end
    end
end

%% Calculations for each filter model and shift, positive y direction shifts
nReplications = 16;
filterModels = {[], 'photocurrentImpulseResponseBased', 'watsonFilter'};
% nReplications = 1;
% filterModels = {[]};
nFilterModels = length(filterModels);
watsonParams_tau = 12;
noiseSds = [20 20 20];
jitterRangeMinutes = 3;
%parfor rr = 1:nReplications
for rr = 8:nReplications
    for ff = 1:nFilterModels
        jitterMinutesX = jitterRangeMinutes*(rand-0.5);
        jitterMinutesY = jitterRangeMinutes*(rand-0.5);
        fileSuffix = sprintf('%s_%s_Rep_%d_filter_%d',calcName,'Intermixed',rr,ff);
        fprintf('%s\n',fileSuffix);
        BerkeleyAOTumblingEThreshold( ...
            'fastParams', false, ...
            'rngSeed', 0, ...
            'eccDegs', eccDegs, ...
            'chromaSpecification_backgroundRGB', [1 0 0], ...
            'chromaSpecification_foregroundRGB', [0 0 0], ...
            'jitterMinutesX', jitterMinutesX, ...
            'jitterMinutesY', jitterMinutesY, ...
            'temporalModulationParams_numFrame', totalFrames, ...
            'temporalModulationParams_xShiftPerFrameMin', theShiftsX, ...
            'temporalModulationParams_yShiftPerFrameMin', theShiftsY, ...
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
            'fileSuffix', fileSuffix, ...
            'validationThresholds',[]);
    end
end



