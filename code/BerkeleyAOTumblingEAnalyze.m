% BerkeleyAOTumblingEAnalyze
%
% Collect up and analyze output of the calcuations

%% Clear
clear; close all;

%% Parameters
calcName = 'Calcs5';
theShifts = [0 1 2 4];
nShifts = length(theShifts);
nReplications = 8;
filterModels = {'No temporal filter', 'Photocurrent', 'Watson w/ Adaptation'};
nFilterModels = length(filterModels);
shiftDirections = [90 270 0 180];
nShiftDirections = length(shiftDirections);
eOrientations = [0 90 180 270];
nEOrientations = length(eOrientations);

%% Set up summary filename and output dir, where we read in data from
options.rootPath = getpref('ISETBerkeleyAOTumblingE','dataPath');
%options.outputFiguresDir =  fullfile(options.rootPath,'figures',summaryFileName);

%% Loop over data and read it in
for ss = 1:nShifts
    for rr = 1:nReplications
        for ff = 1:length(filterModels)
            % Pos Y shift
            options.fileSuffix = sprintf('%s_posYShift_%d_Rep_%d_filter_%d',calcName,ss,rr,ff);
            summaryFileName = ['BerkeleyAOTumblingEThreshold_' options.fileSuffix];
            options.outputResultsDir = fullfile(options.rootPath,'results',summaryFileName);
            theData{1,ss,rr,ff} = load(fullfile(options.outputResultsDir,[summaryFileName '.mat']), ...
                'options','pCorrect','pRespondAAWithStimBB','tByTPerformance','tByTResponseAlternatives','tByTStimAlternatives');

            % Neg Y shift
            options.fileSuffix = sprintf('%s_negYShift_%d_Rep_%d_filter_%d',calcName,ss,rr,ff);
            summaryFileName = ['BerkeleyAOTumblingEThreshold_' options.fileSuffix];
            options.outputResultsDir = fullfile(options.rootPath,'results',summaryFileName);
            theData{2,ss,rr,ff} = load(fullfile(options.outputResultsDir,[summaryFileName '.mat']), ...
                'options','pCorrect','pRespondAAWithStimBB','tByTPerformance','tByTResponseAlternatives','tByTStimAlternatives');

            % Pos X shift
            options.fileSuffix = sprintf('%s_posXShift_%d_Rep_%d_filter_%d',calcName,ss,rr,ff);
            summaryFileName = ['BerkeleyAOTumblingEThreshold_' options.fileSuffix];
            options.outputResultsDir = fullfile(options.rootPath,'results',summaryFileName);
            theData{3,ss,rr,ff} = load(fullfile(options.outputResultsDir,[summaryFileName '.mat']), ...
                'options','pCorrect','pRespondAAWithStimBB','tByTPerformance','tByTResponseAlternatives','tByTStimAlternatives');

            % Neg X shift
            options.fileSuffix = sprintf('%s_negXShift_%d_Rep_%d_filter_%d',calcName,ss,rr,ff);
            summaryFileName = ['BerkeleyAOTumblingEThreshold_' options.fileSuffix];
            options.outputResultsDir = fullfile(options.rootPath,'results',summaryFileName);
            theData{4,ss,rr,ff} = load(fullfile(options.outputResultsDir,[summaryFileName '.mat']), ...
                'options','pCorrect','pRespondAAWithStimBB','tByTPerformance','tByTResponseAlternatives','tByTStimAlternatives');
        end
    end
end

%% Aggregate data
theColors = [ [0.25 0.25 0.5] ; [0.75 0.75 0]; [0 0 1] ];

% Keep directions separate montage
montageByDirectionFigureNormalized = figure; clf;
set(gcf,'Position',[100,100,1500,750]);
set(gca,'FontName','Helvetica','FontSize',14);
nextMontageByDirectionPlot = 1;

% Error by direction montage
errorMontageByDirectionFigure = figure; clf;
set(gcf,'Position',[100,100,1500,750]);
set(gca,'FontName','Helvetica','FontSize',14);
nextErrorMontageByDirectionNormalizedPlot = 1;

summaryFigureNormalized = figure; clf;
set(gcf,'Position',[100,100,1200,400]);
set(gca,'FontName','Helvetica','FontSize',14);
nextSummaryPlotNormalized = 1;

% Raw version (not normalized)
rawVersionPlot = false;
if (rawVersionPlot)
    summaryFigureRaw = figure; clf;
    set(gcf,'Position',[100,100,1200,400]);
    set(gca,'FontName','Helvetica','FontSize',14);
    nextSummaryPlotRaw = 1;
end

[~,index] = sort(theShifts);
index0 = find(theShifts == 0);
for ff = 1:nFilterModels
    for dd = 1:nShiftDirections
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

        % Get the data organized by shift for this motion direction and filter, averaged over replications
        [meanPCorrectByShift{dd,ff}, semPCorrectByShift{dd,ff}, ...
            meanPRespondAAWithStimBBByShift{dd,ff}, semPRespondAAWithStimBBByShift{dd,ff}, ...
            meanPIncorrect180DegByOrientationByShift{dd,ff}, ...
            meanPIncorrect90DegByOrientationByShift{dd,ff}] = ...
            GetPCorrectByShift(theData,dd,ff);
        meanPIncorrect180DegByShift{dd,ff} = sum(meanPIncorrect180DegByOrientationByShift{dd,ff},1);
        meanPIncorrect90DegByShift{dd,ff} = sum(meanPIncorrect90DegByOrientationByShift{dd,ff},1);

        % Get parallel and orthogonal pCorrect
        %
        % This gets done by analyzing the confusion matrix.
        for ss = 1:length(meanPCorrectByShift{dd,ff})
            % Average performance over the two parallel directions
            parMeanPCorrectByShift{dd,ff}(ss) = (meanPRespondAAWithStimBBByShift{dd,ff}(parEOrientationIndex(1),parEOrientationIndex(1),ss) + ...
                meanPRespondAAWithStimBBByShift{dd,ff}(parEOrientationIndex(2),parEOrientationIndex(2),ss) )/2;
            parMeanPIncorrect90DegByShift{dd,ff}(ss) = (meanPIncorrect90DegByOrientationByShift{dd,ff}(parEOrientationIndex(1),ss) + ...
                meanPIncorrect90DegByOrientationByShift{dd,ff}(parEOrientationIndex(2),ss) )/2;
            parMeanPIncorrect180DegByShift{dd,ff}(ss) = (meanPIncorrect180DegByOrientationByShift{dd,ff}(parEOrientationIndex(1),ss) + ...
                meanPIncorrect180DegByOrientationByShift{dd,ff}(parEOrientationIndex(2),ss) )/2;

            % And the two perpindicular directions
            perpMeanPCorrectByShift{dd,ff}(ss) = (meanPRespondAAWithStimBBByShift{dd,ff}(perpEOrientationIndex(1),perpEOrientationIndex(1),ss) + ...
                meanPRespondAAWithStimBBByShift{dd,ff}(perpEOrientationIndex(2),perpEOrientationIndex(2),ss))/2;
            perpMeanPIncorrect90DegByShift{dd,ff}(ss) = (meanPIncorrect90DegByOrientationByShift{dd,ff}(perpEOrientationIndex(1),ss) + ...
                meanPIncorrect90DegByOrientationByShift{dd,ff}(perpEOrientationIndex(2),ss) )/2;
            perpMeanPIncorrect180DegByShift{dd,ff}(ss) = (meanPIncorrect180DegByOrientationByShift{dd,ff}(perpEOrientationIndex(1),ss) + ...
                meanPIncorrect180DegByOrientationByShift{dd,ff}(perpEOrientationIndex(2),ss) )/2;

            % As a check, make sure the parallel and perpindicular average
            % to the omnibus number
            meanPCorrectByShiftCheck{dd,ff}(ss) = (parMeanPCorrectByShift{dd,ff}(ss) + perpMeanPCorrectByShift{dd,ff}(ss))/2;
            if (abs(meanPCorrectByShift{dd,ff}(ss) - meanPCorrectByShiftCheck{dd,ff}(ss)) > 1e-10)
                error('Inconsistent accounting of responses');
            end

            % Rough estimate of SEMs. Average parallel and perpindicular
            % and then divide by sqrt 2.
            parSemPCorrectByShift{dd,ff}(ss) = (semPRespondAAWithStimBBByShift{dd,ff}(parEOrientationIndex(1),parEOrientationIndex(1),ss) + ...
                semPRespondAAWithStimBBByShift{dd,ff}(parEOrientationIndex(2),parEOrientationIndex(2),ss) )/(2*sqrt(2));
            perpSemPCorrectByShift{dd,ff}(ss) = (semPRespondAAWithStimBBByShift{dd,ff}(perpEOrientationIndex(1),perpEOrientationIndex(1),ss) + ...
                semPRespondAAWithStimBBByShift{dd,ff}(perpEOrientationIndex(2),perpEOrientationIndex(2),ss))/(2*sqrt(2));
        end

        % Make plot for this direction and filter and add to the  montage
        figure(montageByDirectionFigureNormalized);
        subplot(3,4,nextMontageByDirectionPlot); hold on;
        legendStr = {};
        nextLegend = 1;

        parh = errorbar(theShifts(index),parMeanPCorrectByShift{dd,ff}(index)/parMeanPCorrectByShift{dd,ff}(index0), ...
            parSemPCorrectByShift{dd,ff}(index)/parMeanPCorrectByShift{dd,ff}(index0),'o');
        parh.Color = theColors(1,:);
        parh.LineWidth = 2;
        legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
        plot(theShifts(index),parMeanPCorrectByShift{dd,ff}(index)/parMeanPCorrectByShift{dd,ff}(index0), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
        legendStr{nextLegend} = 'Parallel'; nextLegend = nextLegend + 1;

        perph = errorbar(theShifts(index),perpMeanPCorrectByShift{dd,ff}(index)/perpMeanPCorrectByShift{dd,ff}(index0), ...
            perpSemPCorrectByShift{dd,ff}(index)/perpMeanPCorrectByShift{dd,ff}(index0),'o');
        perph.Color = theColors(2,:);
        perph.LineWidth = 2;
        legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
        plot(theShifts(index),perpMeanPCorrectByShift{dd,ff}(index)/perpMeanPCorrectByShift{dd,ff}(index0), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));
        legendStr{nextLegend} = 'Orthogonal'; nextLegend = nextLegend + 1;

        xlim([0 5]);
        ylim([0.7 1.1]);
        xlabel('Shift Per Frame (minutes)','FontSize',14);
        ylabel({'Probability Correct' ; '(normalized)'},'FontSize',14);
        title({sprintf('Filter: %s',filterModels{ff}) ;  sprintf('Shift Direction; %d)',shiftDirections(dd))},'FontSize',16);
        legend(legendStr,'FontSize',14,'Location','SouthEast');

        % Move to next subplot
        nextMontageByDirectionPlot = nextMontageByDirectionPlot + 1;

        % Error analysis figure
        figure(errorMontageByDirectionFigure);
        subplot(3,4,nextErrorMontageByDirectionNormalizedPlot); hold on;
        legendStr = {};
        nextLegend = 1;

        plot(theShifts(index),parMeanPIncorrect90DegByShift{dd,ff}(index)./(parMeanPIncorrect90DegByShift{dd,ff}(index)+parMeanPIncorrect180DegByShift{dd,ff}(index)), ...
            '-o','MarkerSize',10,'LineWidth',3,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
        legendStr{nextLegend} = 'Parallel'; nextLegend = nextLegend + 1;

        plot(theShifts(index),perpMeanPIncorrect90DegByShift{dd,ff}(index)./(perpMeanPIncorrect90DegByShift{dd,ff}(index)+perpMeanPIncorrect180DegByShift{dd,ff}(index)), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));
        legendStr{nextLegend} = 'Orthogonal'; nextLegend = nextLegend + 1;

        plot(theShifts,0.67*ones(size(theShifts)),'k:');

        xlim([0 5]);
        ylim([0.25 0.75]);
        xlabel('Shift Per Frame (minutes)','FontSize',14);
        ylabel({'Probability 90 deg error' ; ''},'FontSize',14);
        title({sprintf('Filter: %s',filterModels{ff}) ;  sprintf('Shift Direction; %d)',shiftDirections(dd))},'FontSize',16);
        legend(legendStr,'FontSize',14,'Location','SouthEast');

        nextErrorMontageByDirectionNormalizedPlot = nextErrorMontageByDirectionNormalizedPlot + 1;
    end

    % Get average performance over directions
    parMeanPCorrectByShiftByDirection{ff} = zeros(nShifts,nShiftDirections);
    perpMeanPCorrectByShiftByDirection{ff} = zeros(nShifts,nShiftDirections);
    for dd = 1:nShiftDirections
        parMeanPCorrectByShiftByDirection{ff}(:,dd) = parMeanPCorrectByShift{dd,ff};
        perpMeanPCorrectByShiftByDirection{ff}(:,dd) = perpMeanPCorrectByShift{dd,ff};
    end
    parMeanPCorrectByShiftMeanOverDirection{ff} = mean(parMeanPCorrectByShiftByDirection{ff},2);
    perpMeanPCorrectByShiftMeanOverDirection{ff} = mean(perpMeanPCorrectByShiftByDirection{ff},2);
    parMeanPCorrectByShiftSemOverDirection{ff} = std(parMeanPCorrectByShiftByDirection{ff},[],2)/sqrt(nShiftDirections);
    perpMeanPCorrectByShiftSemOverDirection{ff} = std(perpMeanPCorrectByShiftByDirection{ff},[],2)/sqrt(nShiftDirections);

    % Make plot for this  filter and add to the  montage
    % This version normalized
    figure(summaryFigureNormalized);
    subplot(1,3,nextSummaryPlotNormalized); hold on;
    legendStr = {};
    nextLegend = 1;

    parh = errorbar(theShifts(index),parMeanPCorrectByShiftMeanOverDirection{ff}(index)/parMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        parMeanPCorrectByShiftSemOverDirection{ff} (index)/parMeanPCorrectByShiftMeanOverDirection{ff}(index0),'o');
    parh.Color = theColors(1,:);
    parh.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts(index),parMeanPCorrectByShiftMeanOverDirection{ff}(index)/parMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
    legendStr{nextLegend} = 'Parallel'; nextLegend = nextLegend + 1;

    perph = errorbar(theShifts(index),perpMeanPCorrectByShiftMeanOverDirection{ff}(index)/perpMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        perpMeanPCorrectByShiftSemOverDirection{ff} (index)/perpMeanPCorrectByShiftMeanOverDirection{ff}(index0),'o');
    perph.Color = theColors(2,:);
    perph.LineWidth = 2;
    legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
    plot(theShifts(index),perpMeanPCorrectByShiftMeanOverDirection{ff}(index)/perpMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));
    legendStr{nextLegend} = 'Orthogonal'; nextLegend = nextLegend + 1;

    % Add the normalized average
    plot(theShifts(index),(perpMeanPCorrectByShiftMeanOverDirection{ff}(index)/perpMeanPCorrectByShiftMeanOverDirection{ff}(index0) + ...
        parMeanPCorrectByShiftMeanOverDirection{ff}(index)/parMeanPCorrectByShiftMeanOverDirection{ff}(index0))/2, ...
        'k:','LineWidth',2);
    legendStr{nextLegend} = 'Average'; nextLegend = nextLegend + 1;

    xlim([0 5]);
    ylim([0.7 1.1]);
    xlabel('Shift Per Frame (minutes)','FontSize',14);
    ylabel({'Probability Correct' ; '(normalized)'},'FontSize',14);
    title(sprintf('Filter: %s',filterModels{ff}),'FontSize',16);
    legend(legendStr,'FontSize',14,'Location','SouthEast');
    nextSummaryPlotNormalized = nextSummaryPlotNormalized + 1;

    % Make plot for this  filter and add to the  montage
    % This version raw.  It doesn't add that much information in the end,
    % so can turn it off.
    if (rawVersionPlot)
        figure(summaryFigureRaw);
        subplot(1,3,nextSummaryPlotRaw); hold on;
        legendStr = {};
        nextLegend = 1;

        parh = errorbar(theShifts(index),parMeanPCorrectByShiftMeanOverDirection{ff}(index), ...
            parMeanPCorrectByShiftSemOverDirection{ff} (index),'o');
        parh.Color = theColors(1,:);
        parh.LineWidth = 2;
        legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
        plot(theShifts(index),parMeanPCorrectByShiftMeanOverDirection{ff}(index), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
        legendStr{nextLegend} = 'Parallel'; nextLegend = nextLegend + 1;

        perph = errorbar(theShifts(index),perpMeanPCorrectByShiftMeanOverDirection{ff}(index), ...
            perpMeanPCorrectByShiftSemOverDirection{ff} (index),'o');
        perph.Color = theColors(2,:);
        perph.LineWidth = 2;
        legendStr{nextLegend} = ''; nextLegend = nextLegend + 1;
        plot(theShifts(index),perpMeanPCorrectByShiftMeanOverDirection{ff}(index), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));
        legendStr{nextLegend} = 'Orthogonal'; nextLegend = nextLegend + 1;

        % Add the average
        plot(theShifts(index),(perpMeanPCorrectByShiftMeanOverDirection{ff}(index) + ...
            parMeanPCorrectByShiftMeanOverDirection{ff}(index))/2, ...
            'k:','LineWidth',2);
        legendStr{nextLegend} = 'Average'; nextLegend = nextLegend + 1;

        xlim([0 5]);
        ylim([0.25 1]);
        xlabel('Shift Per Frame (minutes)','FontSize',14);
        ylabel({'Probability Correct' ; '(raw)'},'FontSize',14);
        title(sprintf('Filter: %s',filterModels{ff}),'FontSize',16);
        legend(legendStr,'FontSize',14,'Location','SouthEast');
        nextSummaryPlotRaw = nextSummaryPlotRaw + 1;
    end
end

%% Get pCorrect and other wonderful things as a function base shift
%
% Mean taken over replications
function [meanPCorrectByShift, semPCorrectByShift, ...
    meanPRespondAAWithStimBBByShift, semPRespondAAWithStimBBByShift, ...
    meanPIncorrect180DegByOrientationByShift, meanPIncorrect90DegByOrientationByShift] = ...
    GetPCorrectByShift(theData,whichShiftDirection,whichFilter)

% Get number of filters and replications
nShifts = size(theData,2);
nReplications = size(theData,3);

% We really count on this being 4 in this routine.
% Don't change without thinking hard.
nEOrientations = 4;

% Allocate output variables
meanPCorrectByShift = zeros(1,nShifts);
semPCorrectByShift = zeros(1,nShifts);
meanPRespondAAWithStimBBByShift = zeros(nEOrientations,nEOrientations,nShifts);
semPRespondAAWithStimBBByShift = zeros(nEOrientations,nEOrientations,nShifts);
meanPIncorrect180DegByOrientationByShift = zeros(1,nShifts);
meanPIncorrect90DegByOrientationByShift = zeros(1,nShifts);

% Sanity check
if (size(theData{whichShiftDirection,1,1,whichFilter}.pRespondAAWithStimBB,1) ~= nEOrientations | ...
        size(theData{whichShiftDirection,1,1,whichFilter}.pRespondAAWithStimBB,2) ~= nEOrientations)
    errorr('Saved confusion matrix is not the expected size');
end

% Average over replications for each shift
for ss = 1:nShifts
    % Get omnibus pCorrect and the confusion matrix for each shift,
    % averaged over replications.  Also SEMs of these quantities.
    pCorrectByShiftTemp = zeros(1,nReplications);
    pRespondAAWithStimBBByShiftTemp = zeros(nEOrientations,nEOrientations,nReplications);
    for rr = 1:nReplications
        pCorrectByShiftTemp(rr) = theData{whichShiftDirection,ss,rr,whichFilter}.pCorrect;
        pRespondAAWithStimBBByShiftTemp(:,:,rr) = theData{whichShiftDirection,ss,rr,whichFilter}.pRespondAAWithStimBB;
        if (abs( pCorrectByShiftTemp(rr)) ~= mean(diag(pRespondAAWithStimBBByShiftTemp(:,:,rr))) > 1e-10)
            error('Inconsistent accounting of responses');
        end
    end
    meanPCorrectByShift(ss) = mean(pCorrectByShiftTemp);
    semPCorrectByShift(ss) = std(pCorrectByShiftTemp)/sqrt(nReplications);
    meanPRespondAAWithStimBBByShift(:,:,ss) = mean(pRespondAAWithStimBBByShiftTemp,3);
    semPRespondAAWithStimBBByShift(:,:,ss) = std(pRespondAAWithStimBBByShiftTemp,[],3)/sqrt(nReplications);

    % The columns of the confusion matrix must sum to 1, because the
    % rows correspond to the choices for each stimulus.  Check.
    for cc = 1:size(meanPRespondAAWithStimBBByShift,1)
        if (abs(sum(meanPRespondAAWithStimBBByShift(:,cc,ss)) - 1) > 1e-10)
            error('Response probs for a given stimulus do not sum to unity');
        end
    end

    % Get probability of choosing a 180 deg rotated E.
    for oo = 1:nEOrientations
        pIncorrectByShiftTemp = [];
        index180 = mod(oo-1+2,nEOrientations)+1;
        pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; meanPRespondAAWithStimBBByShift(index180,oo,ss)];
        meanPIncorrect180DegByOrientationByShift(oo,ss) = sum(pIncorrectByShiftTemp);
    end

    % Get probability of choosing a 90 deg rotated E
    for oo = 1:nEOrientations
        pIncorrectByShiftTemp = [];
        index90 = mod(oo-1+1,nEOrientations) + 1;
        pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; meanPRespondAAWithStimBBByShift(index90,oo,ss)];
        index90 = mod(oo-1+3,nEOrientations) + 1;
        pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; meanPRespondAAWithStimBBByShift(index90,oo,ss)];
        meanPIncorrect90DegByOrientationByShift(oo,ss) = sum(pIncorrectByShiftTemp);
    end

    % Sanity checks
    if (abs(meanPCorrectByShift(ss) - mean(diag(meanPRespondAAWithStimBBByShift(:,:,ss)))) > 1e-10)
        error('Inconsistent accounting of responses');
    end

    tempPCorrect = diag(meanPRespondAAWithStimBBByShift(:,:,ss));
    for oo = 1:nEOrientations
        if (abs(tempPCorrect(oo) + meanPIncorrect90DegByOrientationByShift(oo,ss) + meanPIncorrect180DegByOrientationByShift(oo,ss) - 1) > 1e-10)
            error('Prob of correct and incorrect responses do not sum to unity');
        end
    end
end
end
