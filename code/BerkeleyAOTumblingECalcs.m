% TODO
%
% 1) Check that background scene sequence really is
% 2) Consider computing on contrast
% 3) Accumulate data
% 4) Control plotting

%% Clear
clear; close all;

%% Visualization
visualizeScene = false;

%% Parameters
%
% Letter size
letterSizeMinutes = 10;

% Frames off, on , off in native experimental frame duration.
% These get expanded if we simulate faster.
baseOffFramesStart = 1;
baseOnFrames = 3;
baseOffFramesEnd = 3;

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
nTest = 512;

% Background info
backgroundRGB = [1 0 0];
backgroundRGBPerFrame = backgroundRGB(ones(totalFrames,1),:);
foregroundRGB = [0 0 0];

% Set up y shift vectors for each of the three steps
nYShifts = 3;
rawYShiftMinutes = 2;
for ss = 1:nYShifts
    baseShiftMinutes = (ss-1)*rawYShiftMinutes;
    shiftIndex = 1;
    for jj = 0:baseOnFrames-1
        for ii = 1:frameMultiplier
            yShiftOn{ss}(shiftIndex) = jj*baseShiftMinutes;
            shiftIndex = shiftIndex + 1;
        end
    end
    yShift{ss} = zeros(1,totalFrames);
    yShift{ss}(offFramesStart+1:offFramesStart+onFrames) = yShiftOn{ss};
end

%% Calculations for each filter model and shift
nReplications = 3;
%filterModels = {[], 'photocurrentImpulseResponseBased', 'watsonFilter'};
filterModels = {'watsonFilter'};
watsonParams_tau = 12;
noiseSds = [20 20 20];
yJitterRangeMinutes = 2;
for ss = 1:nYShifts
    yJitterMinutes(ss) = yJitterRangeMinutes*(rand-0.5);
    for rr = 1:nReplications
        for ff = 1:length(filterModels)
            [logThreshold(ss,rr,ff), logMAR(ss,rr,ff), questObj{ss,rr,ff}, ...
                psychometricFunction{ss,rr,ff}, fittedPsychometricParamsParams{ss,rr,ff}, ...
                trialByTrialStimulusAlternatives{ss,rr,ff},trialByTrialPerformancePhotocurrent{ss,rr,ff}] = ...
                BerkeleyAOtumblingEThreshold( ...
                'fastParams', false, ...
                'rngSeed', 0, ...
                'visualizeScene', false, ... 
                'chromaSpecification_backgroundRGB', [1 0 0], ...
                'chromaSpecification_foregroundRGB', [0 0 0], ...
                'temporalModulationParams_numFrame', totalFrames, ...
                'temporalModulationParams_xShiftPerFrameMin', zeros(1,totalFrames), ...
                'temporalModulationParams_yShiftPerFrameMin', yShift{ss}+yJitterMinutes(ss), ...
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
                'validationThresholds',[]);

                % Thus counts on us only studying one letter size, which is
                % enforced by the options in the call above.
                keys = psychometricFunction{ss,rr,ff}.keys;
                pCorrect(ss,rr,ff) = psychometricFunction{ss,rr,ff}(keys{1});
        end
    end
end

