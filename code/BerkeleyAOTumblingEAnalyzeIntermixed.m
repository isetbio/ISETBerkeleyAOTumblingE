% BerkeleyAOTumblingEAnalyzeIntermixed
%
% Collect up and analyze output of the calcuations

%% Clear
clear; close all;

%% Verbose
verbose = false;

%% Parameters
%
% These need to match up with those used at calculation time.
calcName = 'Calcs6';
theShifts = [0 1 2 4];
nShifts = length(theShifts);
nReplications = 16;
filterModels = {'None', 'Photocurrent', 'Watson w/ Adaptation'};
nFilterModels = length(filterModels);
shiftDirections = [90 270 0 180];
nDirections = length(shiftDirections);
eOrientations = [0 90 180 270];
nEOrientations = length(eOrientations);

%% Set up summary filename and output dir, where we read in data from
options.rootPath = getpref('ISETBerkeleyAOTumblingE','dataPath');

%% Loop over data and read it in
for rr = 1:nReplications
    for ff = 1:length(filterModels)
        % Full intermixed output, by filter and replication
        options.fileSuffix = sprintf('%s_%s_Rep_%d_filter_%d',calcName,'Intermixed',rr,ff);
        summaryFileName = ['BerkeleyAOTumblingEThreshold_' options.fileSuffix];
        options.outputResultsDir = fullfile(options.rootPath,'results',summaryFileName);
        theData{rr,ff} = load(fullfile(options.outputResultsDir,[summaryFileName '.mat']), ...
            'options','pCorrect','pRespondAAWithStimBB','tByTPerformance','tByTResponseAlternatives','tByTStimAlternatives','options');
        nTest(rr,ff) = length(theData{rr,ff}.tByTStimAlternatives);

        % For debugging, show confusion matrix
        if (verbose)
            figure; imagesc(theData{rr,ff}.pRespondAAWithStimBBMappedToE);
            axis('square');
            title('original');
        end

        % Get total number of unique response alternatives
        nAlternatives = length(unique(theData{rr,ff}.tByTStimAlternatives));

        % Handle 0 motion versus direction ambiguity.  0 motion is the same
        % in all directions, so the confusions need to be redistributed.
        % Assume (but check) that first shift is 0.
        if (theShifts(1) ~= 0)
            error('Analysis code assumes first shift is 0');
        end

        % Handle fact that motion direction for 0 motion is ambiguous, and
        % we want to allocate responses across a four nominal motion
        % directions equally.  Ugh. I think this is right now.
        zeroMotionIndicesBase = [1 2 3 4];
        nZeroMotionResponses(rr,ff) = 0;
        tByTResponseAlternativesNew = theData{rr,ff}.tByTResponseAlternatives;
        for aa = 1:nAlternatives
            % Find the trials that correspond to this alternative
            index = find(theData{rr,ff}.tByTStimAlternatives == aa);

            % If there are any, find those that have a response
            % corresponding to 0 motion and reallocate across motion
            % directions, preserving E identity
            if (~isempty(index))
                for kk = 1:length(index)
                    for yy = 1:nDirections
                        for ll = 1:length(zeroMotionIndicesBase)
                            if (theData{rr,ff}.tByTResponseAlternatives(index(kk)) == zeroMotionIndicesBase(ll) + (yy-1)*nShifts*nEOrientations)
                                % Yes, it's zeroMotionIndicesBase(ll).
                                % Choose reallocation direction
                                nZeroMotionResponses(rr,ff) = nZeroMotionResponses(rr,ff) + 1;
                                whichMotionIndexToUse = randi(nDirections);
                                newResponse = zeroMotionIndicesBase(ll) + (whichMotionIndexToUse-1)*nShifts*nEOrientations;
                                tByTResponseAlternativesNew(index(kk)) = newResponse;
                            end
                        end
                    end
                end
            end
        end
        theData{rr,ff}.tByTResponseAlternatives = tByTResponseAlternativesNew;
        theData{rr,ff}.tByTPerformance  = double(theData{rr,ff}.tByTStimAlternatives == theData{rr,ff}.tByTResponseAlternatives);
        theData{rr,ff}.pRespondAAWithStimBBMappedToE = zeros(nAlternatives,nAlternatives);
        for aa = 1:nAlternatives
            for bb = 1:nAlternatives
                index = find(theData{rr,ff}.tByTStimAlternatives == bb);
                if (~isempty(index))
                    theData{rr,ff}.pRespondAAWithStimBBMappedToE(aa,bb) = sum(theData{rr,ff}.tByTResponseAlternatives(index) == aa)/length(index);
                end
            end
        end
        checkSum = sum(theData{rr,ff}.pRespondAAWithStimBBMappedToE,1);
        if (any(abs(checkSum - 1) > 1e-6))
            error('Some probability went into a black hole');
        end
        if (verbose)
            figure; imagesc(theData{rr,ff}.pRespondAAWithStimBBMappedToE);
            axis('square');
            title('0 shift fix');
            title(sprintf('%d total zero motion responses, %d zero motion trials',nZeroMotionResponses(rr,ff),nTest(rr,ff)/(nDirections*nEOrientations)));
        end
     
        % Map alternatives down from N to 4 (for E direction), and get
        % confusion matrix for E orientation, averaged over direction and
        % shift
        theData{rr,ff}.tByTStimAlternativesMappedToE = rem(theData{rr,ff}.tByTStimAlternatives-1,4)+1;
        theData{rr,ff}.tByTResponseAlternativesMappedToE = rem(theData{rr,ff}.tByTResponseAlternatives-1,4)+1;
        nAlternativesE = length(unique(theData{rr,ff}.tByTStimAlternativesMappedToE));
        theData{rr,ff}.pRespondAAWithStimBBMappedToE = zeros(nAlternativesE,nAlternativesE);
        for aa = 1:nAlternativesE
            for bb = 1:nAlternativesE
                index = find(theData{rr,ff}.tByTStimAlternativesMappedToE == bb);
                if (isempty(index))
                    error('A stimulus alternative vanished');
                end
                theData{rr,ff}.pRespondAAWithStimBBMappedToE(aa,bb) = sum(theData{rr,ff}.tByTResponseAlternativesMappedToE (index) == aa)/length(index);
            end
        end
        checkSum = sum(theData{rr,ff}.pRespondAAWithStimBBMappedToE,1);
        if (max(abs(checkSum - 1) > 1e-6))
            error('Some probability went into a black hole');
        end
        if (verbose)
            figure; imagesc(theData{rr,ff}.pRespondAAWithStimBBMappedToE);
            axis('square');
        end
    end
end

%% Report mean number of zero motion responses
fprintf('Mean number of zero motion responses: %d, expect %d\n',round(mean(nZeroMotionResponses(:))),mean(nTest(:))/(nDirections*nEOrientations));

%% Do various aggregation by filter and hold onto the pieces
%
% This evolved a bit over time, and is a bit redundant with some of what we
% did above.  The function below is at the end of this file.
for ff = 1:nFilterModels
    [meanPCorrectByShift{ff}, semPCorrectByShift{ff}, ...
        meanPRespondAAWithStimBBMappedToEByShift{ff}, semPRespondAAWithStimBBMappedToEByShift{ff}, ...
        meanPIncorrect180DegByOrientationByShift{ff}, meanPIncorrect90DegByOrientationByShift{ff}, ...
        pRespondAAWithStimBBByStimDirByStimShift{ff},meanPRespondAAWithStimBBByStimDirByStimShift{ff}] = ...
        GetPCorrectByShift(theData,ff);
end

[~,index] = sort(theShifts);
index0 = find(theShifts == 0);
for ff = 1:nFilterModels
    for dd = 1:nDirections
        % Sort out relation between this shift direction and the four E
        % orientations, parallel and perpindicular motion to E bar
        % orientation.
        parEOrientations = [shiftDirections(dd) rem(shiftDirections(dd)+180,360)];
        perpEOrientations = setdiff(eOrientations,parEOrientations);
        parEOrientationIndex = [];
        perpEOrientationIndex = [];
        for oo = 1:length(parEOrientations)
            parEOrientationIndex = [parEOrientationIndex find(parEOrientations(oo) == eOrientations)];
        end
        for oo = 1:length(perpEOrientations)
            perpEOrientationIndex = [perpEOrientationIndex find(perpEOrientations(oo) == eOrientations)];
        end
        if (length(parEOrientationIndex) ~= 2 | length(perpEOrientationIndex) ~= 2)
            error('Par/perp indexing not computed correctly');
        end

        % Adjust matrices into canonical order.  First two columns/row are
        % E's parallel with direction of motion, last two are
        % perpindicular.
        for ss = 1:nShifts
            for rr  = 1:nReplications
                pRespondAAWithStimBBCanonical{ff}{dd,ss,rr} = ...
                    pRespondAAWithStimBBByStimDirByStimShift{ff}{dd,ss,rr}([parEOrientationIndex perpEOrientationIndex],[parEOrientationIndex perpEOrientationIndex]);

                % These should be understood as AND.  For example,
                % pCorrectPar means probability of correct and parallel,
                % not probability of correct given parallel
                pCorrectPar{ff}(dd,ss,rr) = (pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(1,1) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(2,2))/nEOrientations;
                pCorrectPerp{ff}(dd,ss,rr) = (pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(3,3) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(4,4))/nEOrientations;
                pIncorrect180Par{ff}(dd,ss,rr) = (pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(1,2) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(2,1))/nEOrientations;
                pIncorrect180Perp{ff}(dd,ss,rr) = (pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(3,4) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(4,3))/nEOrientations;
                pIncorrect90Par{ff}(dd,ss,rr) = (pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(3,1) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(3,2) + ...
                    pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(4,1) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(4,2))/nEOrientations;
                 pIncorrect90Perp{ff}(dd,ss,rr) = (pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(1,3) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(1,4) + ...
                    pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(2,3) + pRespondAAWithStimBBCanonical{ff}{dd,ss,rr}(2,4))/nEOrientations;
                 
                 plotErrorPar{ff}(dd,ss,rr) = pIncorrect90Par{ff}(dd,ss,rr)/(pIncorrect90Par{ff}(dd,ss,rr) + pIncorrect180Par{ff}(dd,ss,rr));
                 plotErrorPerp{ff}(dd,ss,rr) = pIncorrect90Perp{ff}(dd,ss,rr)/(pIncorrect90Perp{ff}(dd,ss,rr) + pIncorrect180Perp{ff}(dd,ss,rr));
end

            % Get mean over replications
            meanPCorrectPar{ff}(dd,ss) = mean(pCorrectPar{ff}(dd,ss,:));
            meanPCorrectPerp{ff}(dd,ss) = mean(pCorrectPerp{ff}(dd,ss,:));
            meanPCorrect{ff}(dd,ss) = meanPCorrectPar{ff}(dd,ss) + meanPCorrectPerp{ff}(dd,ss);
            meanPIncorrect180Par{ff}(dd,ss) = mean(pIncorrect180Par{ff}(dd,ss,:));
            meanPIncorrect180Perp{ff}(dd,ss) = mean(pIncorrect180Perp{ff}(dd,ss,:));
            meanPIncorrect90Par{ff}(dd,ss) = mean(pIncorrect90Par{ff}(dd,ss,:));
            meanPIncorrect90Perp{ff}(dd,ss) = mean(pIncorrect90Perp{ff}(dd,ss,:));
            if (abs(meanPCorrectPar{ff}(dd,ss) + meanPCorrectPerp{ff}(dd,ss) +  ...
                    meanPIncorrect180Par{ff}(dd,ss) +  meanPIncorrect180Perp{ff}(dd,ss) + ...
                    meanPIncorrect90Par{ff}(dd,ss) + meanPIncorrect90Perp{ff}(dd,ss)- 1) > 1e-6)
                error('Some probability went away');
            end

            meanPlotErrorPar{ff}(dd,ss) = mean(plotErrorPar{ff}(dd,ss,:));
            meanPlotErrorPerp{ff}(dd,ss) = mean(plotErrorPerp{ff}(dd,ss,:));
        end
    end

    % Average over directions and get standard errors
    for ss = 1:nShifts
        grandMeanPCorrectPar{ff}(ss) = mean(meanPCorrectPar{ff}(:,ss));
        semMeanPCorrectPar{ff}(ss) = std(meanPCorrectPar{ff}(:,ss))/sqrt(nDirections);
        grandMeanPCorrectPerp{ff}(ss) = mean(meanPCorrectPerp{ff}(:,ss));
        semMeanPCorrectPerp{ff}(ss) = std(meanPCorrectPerp{ff}(:,ss))/sqrt(nDirections);
        grandMeanPCorrect{ff}(ss) = mean(meanPCorrect{ff}(:,ss));
        semMeanPCorrect{ff}(ss) = std(meanPCorrect{ff}(:,ss))/sqrt(nDirections);
        grandMeanPIncorrect180Par{ff}(ss) = mean(meanPIncorrect180Par{ff}(:,ss));
        semMeanPIncorrect180Par{ff}(ss) = std(meanPIncorrect180Par{ff}(:,ss))/sqrt(nDirections);
        grandMeanPIncorrect180Perp{ff}(ss) = mean(meanPIncorrect180Perp{ff}(:,ss));
        semMeanPIncorrect180Perp{ff}(ss) = std(meanPIncorrect180Perp{ff}(:,ss))/sqrt(nDirections);
        grandMeanPIncorrect90Par{ff}(ss) = mean(meanPIncorrect90Par{ff}(:,ss));
        semMeanPIncorrect90Par{ff}(ss) = std(meanPIncorrect90Par{ff}(:,ss))/sqrt(nDirections);
        grandMeanPIncorrect90Perp{ff}(ss) = mean(meanPIncorrect90Perp{ff}(:,ss));
        semMeanPIncorrect90Perp{ff}(ss) = std(meanPIncorrect90Perp{ff}(:,ss))/sqrt(nDirections);

        grandMeanPlotErrorPar{ff}(ss) = mean(meanPlotErrorPar{ff}(:,ss));
        semMeanPlotErrorPar{ff}(ss) = std(meanPlotErrorPar{ff}(:,ss))/sqrt(nDirections);
        grandMeanPlotErrorPerp{ff}(ss) = mean(meanPlotErrorPerp{ff}(:,ss));
        semMeanPlotErrorPerp{ff}(ss) = std(meanPlotErrorPerp{ff}(:,ss))/sqrt(nDirections);
    end
end

%% Plot setup
theColors = [ [0.25 0.25 0.5] ; [0.75 0.75 0]; [0 0 1] ];

% Normalized summary figure
summaryFigureNormalized = figure; clf;
set(gcf,'Position',[100,100,1200,400]);
set(gca,'FontName','Helvetica','FontSize',14);
nextSummaryPlotNormalized = 1;

% Raw version (not normalized)
summaryFigureRaw = figure; clf;
set(gcf,'Position',[100,100,1200,400]);
set(gca,'FontName','Helvetica','FontSize',14);
nextSummaryPlotRaw = 1;

% Error figure
summaryFigureError = figure; clf;
set(gcf,'Position',[100,100,1200,400]);
set(gca,'FontName','Helvetica','FontSize',14);
nextSummaryPlotError = 1;

%% Plots
for ff = 1:nFilterModels
    % Raw probability plot
    figure(summaryFigureRaw);
    subplot(1,3,nextSummaryPlotRaw); hold on;
    legendStr = {};
    nextLegend = 1;

    parh = errorbar(theShifts/2,grandMeanPCorrectPar{ff}, ...
        semMeanPCorrectPar{ff},'o');
    parh.Color = theColors(1,:);
    parh.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts/2,grandMeanPCorrectPar{ff}, ...
        '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
    legendStr{nextLegend} = 'Parallel'; nextLegend = nextLegend + 1;

    perph = errorbar(theShifts/2,grandMeanPCorrectPerp{ff}, ...
        semMeanPCorrectPerp{ff},'o');
    perph.Color = theColors(2,:);
    perph.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts/2,grandMeanPCorrectPerp{ff}, ...
        '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));
    legendStr{nextLegend} = 'Orthogonal'; nextLegend = nextLegend + 1;

    % Add the sum
    plot(theShifts/2,grandMeanPCorrect{ff}, ...
        'k:','LineWidth',2);
    legendStr{nextLegend} = 'Both'; nextLegend = nextLegend + 1;

    % Tidy up
    xlim([0 4.5]/2);
    ylim([0 0.5]);
    legend(legendStr,'FontSize',14,'Location','SouthEast');
    xlabel('Imposed Motion (Bar Width/Frame');
    ylabel('Raw Fraction Correct')
    title(sprintf('Filter: %s',filterModels{ff}));
    nextSummaryPlotRaw = nextSummaryPlotRaw + 1;

    % Normalzed fraction correct by shift plot
    figure(summaryFigureNormalized);
    set(gcf,"Position",[100 100 2100 700]);
    subplot(1,3,nextSummaryPlotNormalized); hold on;
    set(gca,'FontName','Helvetica','FontSize',18);
    legendStr = {};
    nextLegend = 1;

    parh = errorbar(theShifts/2,grandMeanPCorrectPar{ff}/grandMeanPCorrectPar{ff}(index0), ...
        semMeanPCorrectPar{ff}/grandMeanPCorrectPar{ff}(index0),'o');
    parh.Color = theColors(1,:);
    parh.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts/2,grandMeanPCorrectPar{ff}/grandMeanPCorrectPar{ff}(index0), ...
        '-o','MarkerSize',12,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
    legendStr{nextLegend} = 'Parallel'; nextLegend = nextLegend + 1;

    perph = errorbar(theShifts/2,grandMeanPCorrectPerp{ff}/grandMeanPCorrectPerp{ff}(index0), ...
        semMeanPCorrectPerp{ff}/grandMeanPCorrectPerp{ff}(index0),'o');
    perph.Color = theColors(2,:);
    perph.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts/2,grandMeanPCorrectPerp{ff}/grandMeanPCorrectPerp{ff}(index0), ...
        '-o','MarkerSize',12,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));
    legendStr{nextLegend} = 'Orthogonal'; nextLegend = nextLegend + 1;

    % Add the normalized average
    % plot(theShifts/2,grandMeanPCorrect{ff}/grandMeanPCorrect{ff}(index0), ...
    %     'k:','LineWidth',2);
    % legendStr{nextLegend} = 'Average'; nextLegend = nextLegend + 1;

    % Tidy up
    xlim([0 5]/2);
    ylim([0.7 1.05]);
    xlabel('Motion (Bar Width/Frame','FontSize',20);
    ylabel({'Probability Correct' ; '(normalized)'},'FontSize',20);
    title(sprintf('Filter: %s',filterModels{ff}),'FontSize',20);
    legend(legendStr,'FontSize',14,'Location','SouthEast');
    nextSummaryPlotNormalized = nextSummaryPlotNormalized + 1;

    % Summary figure error
    figure(summaryFigureError);
    set(gcf,"Position",[100 100 2100 700]);
    subplot(1,3,nextSummaryPlotError); hold on;
    set(gca,'FontName','Helvetica','FontSize',18);
    legendStr = {};
    nextLegend = 1;

    parh = errorbar(theShifts/2,grandMeanPlotErrorPar{ff}, ...
        semMeanPlotErrorPar{ff},'o');
    parh.Color = theColors(1,:);
    parh.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts/2,grandMeanPlotErrorPar{ff}, ...
        '-o','MarkerSize',12,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
    legendStr{nextLegend} = 'Parallel'; nextLegend = nextLegend + 1;

    perph = errorbar(theShifts/2,grandMeanPlotErrorPerp{ff}, ...
        semMeanPlotErrorPerp{ff},'o');
    perph.Color = theColors(2,:);
    perph.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts/2,grandMeanPlotErrorPerp{ff}, ...
        '-o','MarkerSize',12,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));
    legendStr{nextLegend} = 'Orthogonal'; nextLegend = nextLegend + 1;

    % Tidy up
    xlim([0 4.5]/2);
    ylim([0.3 0.9]);
    legend(legendStr,'FontSize',14,'Location','SouthEast');
    xlabel('Motion (Bar Width/Frame','FontSize',20);
    ylabel('Fraction 90 deg Error','FontSize',20);
    title(sprintf('Filter: %s',filterModels{ff}),'FontSize',20);
    nextSummaryPlotError = nextSummaryPlotError + 1;
end

%% Get pCorrect and other wonderful things as a function base shift
%
% Mean taken over replications
function [meanPCorrectByShift, semPCorrectByShift, ...
    meanPRespondAAWithStimBBMappedToEByShift, semPRespondAAWithStimBBMappedToEByShift, ...
    meanPIncorrect180DegByOrientationByShift, meanPIncorrect90DegByOrientationByShift, ...
    pRespondAAWithStimBBByStimDirByStimShift,meanPRespondAAWithStimBBByStimDirByStimShift] = ...
    GetPCorrectByShift(theData,whichFilter)
%
% Outputs
%    meanPCorrectByShift - Vector by shift.  For each entry, is the mean of
%        the diagonal entries of meanPRespondAAWithStimBBByStimShift.
%    semPCorrectByShift - Will eventually be the SEM of
%        meanPCorrectByShift, but is currently a vector of NaNs.
%    pRespondAAWithStimBBByStimDirByStimShift - Cell array by direction,
%        shift, and replication.  Each is the confusion matrix by E
%        orientation.
%    meanPRespondAAWithStimBBByStimDirByStimShift - Cell array by direction
%       and shift, averaged over replication. Each is the confusion matrix
%       by E orientation.
%   meanPRespondAAWithStimBBByStimShift - Cell array by shift, averaged
%      over directions. Each is the confusion matrix by E orientation.

% We really count on this being 4 in this routine.
% Don't change without thinking hard.
nEOrientations = 4;

% Get number of filters and replications
nFilters = size(theData,2);
nReplications = size(theData,1);

% Get directions and shifts
%
% The order of the data is (inner to outer): EOrientation, shift, direction
[nDirections,nShifts] = size(theData{1,1}.options.temporalModulationParams_xShiftPerFrameMin);

% Allocate output variables
meanPCorrectByShift = zeros(1,nShifts);
semPCorrectByShift = zeros(1,nShifts);
meanPRespondAAWithStimBBMappedToEByShift = zeros(nEOrientations,nEOrientations,nShifts);
semPRespondAAWithStimBBMappedToEByShift = zeros(nEOrientations,nEOrientations,nShifts);
meanPIncorrect180DegByOrientationByShift = zeros(1,nShifts);
meanPIncorrect90DegByOrientationByShift = zeros(1,nShifts);

% Get indices in output variables for the desired direction
nPerDirection = nEOrientations*nShifts;
nPerShift = nEOrientations;

% Get variables sorted out by stimulus direction, and then collapse over
% response directions
for rr = 1:nReplications
    % For each direction of the stimulus for this replication, pull out
    % each stimulus direction into its own confusion matrix with E
    % orientation and shift size as key variables.
    for dds = 1:nDirections
        pRespondAAWithStimBBByStimDir{dds,rr} = zeros(nEOrientations*nShifts,nEOrientations*nShifts);
        directionIndexS = ((dds-1)*nPerDirection + 1):dds*nPerDirection;
        tempPRespondAAWithStimBB = theData{rr,whichFilter}.pRespondAAWithStimBB(:,directionIndexS);

        % Collapose over response directions, which we don't record in the experiment and
        % thus that we don't know about.
        for ddr = 1:nDirections
            directionIndexR = ((ddr-1)*nPerDirection + 1):ddr*nPerDirection;
            pRespondAAWithStimBBByStimDir{dds,rr} =  pRespondAAWithStimBBByStimDir{dds,rr} + tempPRespondAAWithStimBB(directionIndexR,:);
        end
        if (any(abs(sum(pRespondAAWithStimBBByStimDir{dds,rr},1)-1) > 1e-6))
            error('Failure of probabilities to sum to unity');
        end

        % Subdivide further by separating out the confusion matrices by shift size
        for sss = 1:nShifts
            pRespondAAWithStimBBByStimDirByStimShift{dds,sss,rr} = zeros(nEOrientations,nEOrientations);
            shiftIndexS = ((sss-1)*nPerShift  + 1):sss*nPerShift;
            tempPRespondAAWithStimBB = pRespondAAWithStimBBByStimDir{dds,rr}(:,shiftIndexS);
            for ssr = 1:nShifts
                shiftIndexR = ((ssr-1)*nPerShift  + 1):ssr*nPerShift;
                pRespondAAWithStimBBByStimDirByStimShift{dds,sss,rr} = ...
                    pRespondAAWithStimBBByStimDirByStimShift{dds,sss,rr}  + tempPRespondAAWithStimBB(shiftIndexR,:);
            end
            if (any(abs(sum(pRespondAAWithStimBBByStimDirByStimShift{dds,sss,rr},1)-1) > 1e-6))
                error('Failure of probabilities to sum to unity');
            end
        end
    end
end

%% Average up over replications
for dds = 1:nDirections
    for sss = 1:nShifts
        meanPRespondAAWithStimBBByStimDirByStimShift{dds,sss} = zeros(nEOrientations,nEOrientations);
        for rr = 1:nReplications
            meanPRespondAAWithStimBBByStimDirByStimShift{dds,sss} = ...
                meanPRespondAAWithStimBBByStimDirByStimShift{dds,sss} + pRespondAAWithStimBBByStimDirByStimShift{dds,sss,rr};
        end
        meanPRespondAAWithStimBBByStimDirByStimShift{dds,sss} = meanPRespondAAWithStimBBByStimDirByStimShift{dds,sss}/nReplications;
        if (any(abs(sum(meanPRespondAAWithStimBBByStimDirByStimShift{dds,sss},1)-1) > 1e-6))
            error('Failure of probabilities to sum to unity');
        end
    end
end

%% Average over stim directions
for sss = 1:nShifts
    meanPRespondAAWithStimBBByStimShift{sss} = zeros(nEOrientations,nEOrientations);
    for dds = 1:nDirections
        meanPRespondAAWithStimBBByStimShift{sss} = ...
            meanPRespondAAWithStimBBByStimShift{sss} + meanPRespondAAWithStimBBByStimDirByStimShift{dds,sss};
    end
    meanPRespondAAWithStimBBByStimShift{sss} = meanPRespondAAWithStimBBByStimShift{sss}/nDirections;
    if (any(abs(sum(meanPRespondAAWithStimBBByStimShift{sss} ,1)-1) > 1e-6))
        error('Failure of probabilities to sum to unity');
    end

    % Summarize
    meanPCorrectByShift(sss) = mean(diag(meanPRespondAAWithStimBBByStimShift{sss}));
    semPCorrectByShift(sss) = NaN;

    meanPRespondAAWithStimBBMappedToEByShift(:,:,sss) = meanPRespondAAWithStimBBByStimShift{sss};
    semPRespondAAWIthStimBBMappedToEByShift(:,:,sss) = NaN*meanPRespondAAWithStimBBByStimShift{sss};

end

meanPIncorrect180DegByOrientationByShift = NaN;
meanPIncorrect90DegByOrientationByShift = NaN;

% % Average over replications for each shift
% for ss = 1:nShifts
%     % Get omnibus pCorrect and the confusion matrix for each shift,
%     % averaged over replications.  Also SEMs of these quantities.
%     pCorrectByShiftTemp = zeros(1,nReplications);
%     pRespondAAWithStimBBMappedToEByShiftTemp = zeros(nEOrientations,nEOrientations,nReplications);
%     for rr = 1:nReplications
%         pRespondAAWithStimBBMappedToEByShiftTemp(:,:,rr) = theData{whichDirection,ss,rr,whichFilter}.pRespondAAWithStimBBMappedToE;
%         if (abs( pCorrectByShiftTemp(rr)) ~= mean(diag(pRespondAAWithStimBBMappedToEByShiftTemp(:,:,rr))) > 1e-10)
%             error('Inconsistent accounting of responses');
%         end
%     end
% 
%     pCorrectByShiftTemp(rr) = theData{rr,whichFilter}.pCorrect;
% 
%     meanPCorrectByShift(ss) = mean(pCorrectByShiftTemp);
%     semPCorrectByShift(ss) = std(pCorrectByShiftTemp)/sqrt(nReplications);
%     meanPRespondAAWithStimBBMappedToEByShift(:,:,ss) = mean(pRespondAAWithStimBBMappedToEByShiftTemp,3);
%     semPRespondAAWithStimBBMappedToEByShift(:,:,ss) = std(pRespondAAWithStimBBMappedToEByShiftTemp,[],3)/sqrt(nReplications);
% 
%     % The columns of the confusion matrix must sum to 1, because the
%     % rows correspond to the choices for each stimulus.  Check.
%     for cc = 1:size(meanPRespondAAWithStimBBMappedToEByShift,1)
%         if (abs(sum(meanPRespondAAWithStimBBMappedToEByShift(:,cc,ss)) - 1) > 1e-10)
%             error('Response probs for a given stimulus do not sum to unity');
%         end
%     end
% 
%     % Get probability of choosing a 180 deg rotated E.
%     for oo = 1:nEOrientations
%         pIncorrectByShiftTemp = [];
%         index180 = mod(oo-1+2,nEOrientations)+1;
%         pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; meanPRespondAAWithStimBBMappedToEByShift(index180,oo,ss)];
%         meanPIncorrect180DegByOrientationByShift(oo,ss) = sum(pIncorrectByShiftTemp);
%     end
% 
%     % Get probability of choosing a 90 deg rotated E
%     for oo = 1:nEOrientations
%         pIncorrectByShiftTemp = [];
%         index90 = mod(oo-1+1,nEOrientations) + 1;
%         pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; meanPRespondAAWithStimBBMappedToEByShift(index90,oo,ss)];
%         index90 = mod(oo-1+3,nEOrientations) + 1;
%         pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; meanPRespondAAWithStimBBMappedToEByShift(index90,oo,ss)];
%         meanPIncorrect90DegByOrientationByShift(oo,ss) = sum(pIncorrectByShiftTemp);
%     end
% 
%     % Sanity checks
%     if (abs(meanPCorrectByShift(ss) - mean(diag(meanPRespondAAWithStimBBMappedToEByShift(:,:,ss)))) > 1e-10)
%         error('Inconsistent accounting of responses');
%     end
% 
%     tempPCorrect = diag(meanPRespondAAWithStimBBMappedToEByShift(:,:,ss));
%     for oo = 1:nEOrientations
%         if (abs(tempPCorrect(oo) + meanPIncorrect90DegByOrientationByShift(oo,ss) + meanPIncorrect180DegByOrientationByShift(oo,ss) - 1) > 1e-10)
%             error('Prob of correct and incorrect responses do not sum to unity');
%         end
%     end
% end
end
