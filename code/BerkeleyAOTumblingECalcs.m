% BerkeleyAOTumblingECalcs
%
% Code that models Tuten AO stabilized and not tumbling E experiment.
% This sets up options to match the experiment and runs the calculations.

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
nTest = 512;

% Background info
backgroundRGB = [1 0 0];
backgroundRGBPerFrame = backgroundRGB(ones(totalFrames,1),:);
foregroundRGB = [0 0 0];

% Set up y shift vectors for each of the three steps
rawShiftsMinutes = [0 2 4 1];
nShifts = length(rawShiftsMinutes];
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
nReplications = 3;
filterModels = {[], 'photocurrentImpulseResponseBased', 'watsonFilter'};
watsonParams_tau = 12;
noiseSds = [20 20 20];
jitterRangeMinutes = 2;
for ss = 1:nShifts
    % TEMP
    if (ss == 4)

        jitterMinutes(ss) = jitterRangeMinutes*(rand-0.5);
        for rr = 1:nReplications
            for ff = 1:length(filterModels)
                fileSuffix = sprintf('Calcs_posYShift_%d_Rep_%d_filter_%d',ss,rr,ff);
                fprintf('%s\n',fileSuffix);
                [logThreshold(ss,rr,ff), logMAR(ss,rr,ff), questObj{ss,rr,ff}, ...
                    psychometricFunction{ss,rr,ff}, fittedPsychometricParamsParams{ss,rr,ff}, ...
                    trialByTrialStimulusAlternatives{ss,rr,ff},trialByTrialPerformancePhotocurrent{ss,rr,ff}] = ...
                    BerkeleyAOTumblingEThreshold( ...
                    'fastParams', false, ...
                    'rngSeed', 0, ...
                    'chromaSpecification_backgroundRGB', [1 0 0], ...
                    'chromaSpecification_foregroundRGB', [0 0 0], ...
                    'temporalModulationParams_numFrame', totalFrames, ...
                    'temporalModulationParams_xShiftPerFrameMin', zeros(1,totalFrames), ...
                    'temporalModulationParams_yShiftPerFrameMin', theShift{ss}+jitterMinutes(ss), ...
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
                    'visualizeScene', true, ...
                    'visualizeEsOnMosaic', true, ...
                    'visualizeEsWhichFrames', offFramesStart+1, ...
                    'fileSuffix', fileSuffix, ...
                    'validationThresholds',[]);

                % Thus counts on us only studying one letter size, which is
                % enforced by the options in the call above.
                keys = psychometricFunction{ss,rr,ff}.keys;
                pCorrect(ss,rr,ff) = psychometricFunction{ss,rr,ff}(keys{1});
            end
        end
    end
end

% Neg y shift
for ss = 1:nShifts
    % TEMP
    if (ss == 4)

        jitterMinutes(ss) = jitterRangeMinutes*(rand-0.5);
        for rr = 1:nReplications
            for ff = 1:length(filterModels)
                fileSuffix = sprintf('Calcs_negYShift_%d_Rep_%d_filter_%d',ss,rr,ff);
                fprintf('%s\n',fileSuffix);
                [logThreshold(ss,rr,ff), logMAR(ss,rr,ff), questObj{ss,rr,ff}, ...
                    psychometricFunction{ss,rr,ff}, fittedPsychometricParamsParams{ss,rr,ff}, ...
                    trialByTrialStimulusAlternatives{ss,rr,ff},trialByTrialPerformancePhotocurrent{ss,rr,ff}] = ...
                    BerkeleyAOTumblingEThreshold( ...
                    'fastParams', false, ...
                    'rngSeed', 0, ...
                    'chromaSpecification_backgroundRGB', [1 0 0], ...
                    'chromaSpecification_foregroundRGB', [0 0 0], ...
                    'temporalModulationParams_numFrame', totalFrames, ...
                    'temporalModulationParams_xShiftPerFrameMin', zeros(1,totalFrames), ...
                    'temporalModulationParams_yShiftPerFrameMin', -theShift{ss}+jitterMinutes(ss), ...
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
                    'visualizeScene', true, ...
                    'visualizeEsOnMosaic', true, ...
                    'visualizeEsWhichFrames', offFramesStart+1, ...
                    'fileSuffix', fileSuffix, ...
                    'validationThresholds',[]);

                % Thus counts on us only studying one letter size, which is
                % enforced by the options in the call above.
                keys = psychometricFunction{ss,rr,ff}.keys;
                pCorrect(ss,rr,ff) = psychometricFunction{ss,rr,ff}(keys{1});
            end
        end
    end
end

% Pos x shift
for ss = 1:nShifts
    % TEMP
    if (ss == 4)

        jitterMinutes(ss) = jitterRangeMinutes*(rand-0.5);
        for rr = 1:nReplications
            for ff = 1:length(filterModels)
                fileSuffix = sprintf('Calcs_posXShift_%d_Rep_%d_filter_%d',ss,rr,ff);
                fprintf('%s\n',fileSuffix);
                [logThreshold(ss,rr,ff), logMAR(ss,rr,ff), questObj{ss,rr,ff}, ...
                    psychometricFunction{ss,rr,ff}, fittedPsychometricParamsParams{ss,rr,ff}, ...
                    trialByTrialStimulusAlternatives{ss,rr,ff},trialByTrialPerformancePhotocurrent{ss,rr,ff}] = ...
                    BerkeleyAOTumblingEThreshold( ...
                    'fastParams', false, ...
                    'rngSeed', 0, ...
                    'chromaSpecification_backgroundRGB', [1 0 0], ...
                    'chromaSpecification_foregroundRGB', [0 0 0], ...
                    'temporalModulationParams_numFrame', totalFrames, ...
                    'temporalModulationParams_xShiftPerFrameMin', theShift{ss}+jitterMinutes(ss), ...
                    'temporalModulationParams_yShiftPerFrameMin', zeros(1,totalFrames), ...
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
                    'visualizeScene', true, ...
                    'visualizeEsOnMosaic', true, ...
                    'visualizeEsWhichFrames', offFramesStart+1, ...
                    'fileSuffix', fileSuffix, ...
                    'validationThresholds',[]);

                % Thus counts on us only studying one letter size, which is
                % enforced by the options in the call above.
                keys = psychometricFunction{ss,rr,ff}.keys;
                pCorrect(ss,rr,ff) = psychometricFunction{ss,rr,ff}(keys{1});
            end
        end
    end
end

% Neg x shift
for ss = 1:nShifts
    % TEMP
    if (ss == 4)

        jitterMinutes(ss) = jitterRangeMinutes*(rand-0.5);
        for rr = 1:nReplications
            for ff = 1:length(filterModels)
                fileSuffix = sprintf('Calcs_negXShift_%d_Rep_%d_filter_%d',ss,rr,ff);
                fprintf('%s\n',fileSuffix);
                [logThreshold(ss,rr,ff), logMAR(ss,rr,ff), questObj{ss,rr,ff}, ...
                    psychometricFunction{ss,rr,ff}, fittedPsychometricParamsParams{ss,rr,ff}, ...
                    trialByTrialStimulusAlternatives{ss,rr,ff},trialByTrialPerformancePhotocurrent{ss,rr,ff}] = ...
                    BerkeleyAOTumblingEThreshold( ...
                    'fastParams', false, ...
                    'rngSeed', 0, ...
                    'chromaSpecification_backgroundRGB', [1 0 0], ...
                    'chromaSpecification_foregroundRGB', [0 0 0], ...
                    'temporalModulationParams_numFrame', totalFrames, ...
                    'temporalModulationParams_xShiftPerFrameMin', -theShift{ss}+jitterMinutes(ss), ...
                    'temporalModulationParams_yShiftPerFrameMin', zeros(1,totalFrames), ...
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
                    'visualizeScene', true, ...
                    'visualizeEsOnMosaic', true, ...
                    'visualizeEsWhichFrames', offFramesStart+1, ...
                    'fileSuffix', fileSuffix, ...
                    'validationThresholds',[]);

                % Thus counts on us only studying one letter size, which is
                % enforced by the options in the call above.
                keys = psychometricFunction{ss,rr,ff}.keys;
                pCorrect(ss,rr,ff) = psychometricFunction{ss,rr,ff}(keys{1});
            end
        end
    end
end