% BerkeleyAOTumblingEAnalyze
%
% Collect up and analyze output of the calcuations

%% Clear
clear; close all;

%% Parameters
calcName = 'Calcs';
theShifts = [0 2 4 0.5];
nShifts = length(theShifts);
nReplications = 3;
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
montageFigure = figure; clf; 
set(gcf,'Position',[100,100,1500,750]);
set(gca,'FontName','Helvetica','FontSize',14);
nextMontagePlot = 1;

errorMontageFigure = figure; clf; 
set(gcf,'Position',[100,100,1500,750]);
set(gca,'FontName','Helvetica','FontSize',14);
nextErrorMontagePlot = 1;

summaryFigure = figure; clf; 
set(gcf,'Position',[100,100,1200,400]);
set(gca,'FontName','Helvetica','FontSize',14);
nextSummaryPlot = 1;

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
            meanPIncorrectParByShift{dd,ff}, semPIncorrectParByShift{dd,ff}, ...
            meanPIncorrectPerpByShift{dd,ff},semPIncorrectPerpByShift{dd,ff}] = ...
            GetPCorrectByShift(theData,dd,ff);

        %% Get parallel and orthogonal pCorrect
        %
        % This gets done by analyzing the confusion matrix.
        for ss = 1:length(meanPCorrectByShift{dd,ff})
            % Average performance over the two parallel directions
            parMeanPCorrectByShift{dd,ff}(ss) = (meanPRespondAAWithStimBBByShift{dd,ff}(parEOrientationIndex(1),parEOrientationIndex(1),ss) + ...
                meanPRespondAAWithStimBBByShift{dd,ff}(parEOrientationIndex(2),parEOrientationIndex(2),ss) )/2;

            % And the two perpindicular directions
            perpMeanPCorrectByShift{dd,ff}(ss) = (meanPRespondAAWithStimBBByShift{dd,ff}(perpEOrientationIndex(1),perpEOrientationIndex(1),ss) + ...
                meanPRespondAAWithStimBBByShift{dd,ff}(perpEOrientationIndex(2),perpEOrientationIndex(2),ss))/2;

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
        figure(montageFigure);
        subplot(3,4,nextMontagePlot); hold on;

        parh = errorbar(theShifts(index),parMeanPCorrectByShift{dd,ff}(index)/parMeanPCorrectByShift{dd,ff}(index0), ...
            parSemPCorrectByShift{dd,ff}(index)/parMeanPCorrectByShift{dd,ff}(index0),'o');
            parh.Color = theColors(1,:);
            parh.LineWidth = 2;
        plot(theShifts(index),parMeanPCorrectByShift{dd,ff}(index)/parMeanPCorrectByShift{dd,ff}(index0), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));

        perph = errorbar(theShifts(index),perpMeanPCorrectByShift{dd,ff}(index)/perpMeanPCorrectByShift{dd,ff}(index0), ...
            perpSemPCorrectByShift{dd,ff}(index)/perpMeanPCorrectByShift{dd,ff}(index0),'o');
        perph.Color = theColors(2,:);
        perph.LineWidth = 2;
        plot(theShifts(index),perpMeanPCorrectByShift{dd,ff}(index)/perpMeanPCorrectByShift{dd,ff}(index0), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));

        xlim([0 5]);
        ylim([0.7 1.1]);
        xlabel('Shift Per Frame (minutes)','FontSize',14);
        ylabel({'Probability Correct' ; '(normalized)'},'FontSize',14);
        title({sprintf('Filter: %s',filterModels{ff}) ;  sprintf('Shift Direction; %d)',shiftDirections(dd))},'FontSize',16);

        % Move to next subplot
        nextMontagePlot = nextMontagePlot + 1;

        % Error analysis figure
        figure(errorMontageFigure);
        subplot(3,4,nextErrorMontagePlot); hold on;
        plot(theShifts(index),meanPIncorrectParByShift{dd,ff}(index), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));
        plot(theShifts(index),2*meanPIncorrectPerpByShift{dd,ff}(index), ...
            '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));

        xlim([0 5]);
        ylim([0 0.4]);
        xlabel('Shift Per Frame (minutes)','FontSize',14);
        ylabel({'Probability Correct' ; '(normalized)'},'FontSize',14);
        title({sprintf('Filter: %s',filterModels{ff}) ;  sprintf('Shift Direction; %d)',shiftDirections(dd))},'FontSize',16);
     
        nextErrorMontagePlot = nextErrorMontagePlot + 1;
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
    figure(summaryFigure);
    subplot(1,3,nextSummaryPlot); hold on;

    parh = errorbar(theShifts(index),parMeanPCorrectByShiftMeanOverDirection{ff}(index)/parMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        parMeanPCorrectByShiftSemOverDirection{ff} (index)/parMeanPCorrectByShiftMeanOverDirection{ff}(index0),'o');
    parh.Color = theColors(1,:);
    parh.LineWidth = 2;
    plot(theShifts(index),parMeanPCorrectByShiftMeanOverDirection{ff}(index)/parMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(1,:),'MarkerFaceColor',theColors(1,:));

     perph = errorbar(theShifts(index),perpMeanPCorrectByShiftMeanOverDirection{ff}(index)/perpMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        perpMeanPCorrectByShiftSemOverDirection{ff} (index)/perpMeanPCorrectByShiftMeanOverDirection{ff}(index0),'o');
    perph.Color = theColors(2,:);
    perph.LineWidth = 2;
    plot(theShifts(index),perpMeanPCorrectByShiftMeanOverDirection{ff}(index)/perpMeanPCorrectByShiftMeanOverDirection{ff}(index0), ...
        '-o','MarkerSize',8,'LineWidth',2,'Color',theColors(2,:),'MarkerFaceColor',theColors(2,:));

    xlim([0 5]);
    ylim([0.7 1.1]);
    xlabel('Shift Per Frame (minutes)','FontSize',14);
    ylabel({'Probability Correct' ; '(normalized)'},'FontSize',14);
    title(sprintf('Filter: %s',filterModels{ff}),'FontSize',16);
    nextSummaryPlot = nextSummaryPlot + 1;
end

%% Get pCorrect and other wonderful things as a function base shift
function [meanPCorrectByShift, semPCorrectByShift, ...
    meanPRespondAAWithStimBBByShift, semPRespondAAWithStimBBByShift, ...
    meanPIncorrectParByShift, semPIncorrectParByShift, ...
    meanPIncorrectPerpByShift,semPIncorrectPerpByShift] = ...
    GetPCorrectByShift(theData,whichShiftDirection,whichFilter)

    % Get number of filters and replications
    nShifts = size(theData,2);
    nReplications = size(theData,3);
    nEOrientations = 4;

    % Allocate output variables
    meanPCorrectByShift = zeros(1,nShifts);
    semPCorrectByShift = zeros(1,nShifts);
    meanPRespondAAWithStimBBByShift = zeros(nEOrientations,nEOrientations,nShifts);
    semPRespondAAWithStimBBByShift = zeros(nEOrientations,nEOrientations,nShifts);
    meanPIncorrectParByShift = zeros(1,nShifts);
    semPIncorrectParByShift = zeros(1,nShifts);
    meanPIncorrectPerpByShift = zeros(1,nShifts);
    semPIncorrectPerpByShift = zeros(1,nShifts);

    % Sanity check
    if (size(theData{whichShiftDirection,1,1,whichFilter}.pRespondAAWithStimBB,1) ~= nEOrientations | ...
        size(theData{whichShiftDirection,1,1,whichFilter}.pRespondAAWithStimBB,2) ~= nEOrientations)
        errorr('Saved confusion matrix is not the expected size');
    end

    % Average over replications for each shift
    for ss = 1:nShifts
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
        semPRespondAAWithStimBBByShift(:,:,ss) = std(pRespondAAWithStimBBByShiftTemp,[],3);

        pIncorrectByShiftTemp = [];
        for zz = [2 4]
            pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; diag(meanPRespondAAWithStimBBByShift(:,:,ss),zz)];
            pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; diag(meanPRespondAAWithStimBBByShift(:,:,ss),-zz)];
        end
        meanPIncorrectParByShift(ss) = mean(pIncorrectByShiftTemp);
        semPIncorrectParByShift = std(pIncorrectByShiftTemp)/length(pIncorrectByShiftTemp);

        pIncorrectByShiftTemp = [];
        for zz = [1 3]
            pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; diag(meanPRespondAAWithStimBBByShift(:,:,ss),zz)];
            pIncorrectByShiftTemp = [pIncorrectByShiftTemp ; diag(meanPRespondAAWithStimBBByShift(:,:,ss),-zz)];
        end
        meanPIncorrectPerpByShift(ss) = mean(pIncorrectByShiftTemp);
        semPIncorrectPerpByShift = std(pIncorrectByShiftTemp)/length(pIncorrectByShiftTemp);

        % Sanity check
        if (abs(meanPCorrectByShift(ss) - mean(diag(meanPRespondAAWithStimBBByShift(:,:,ss)))) > 1e-10)
            error('Inconsistent accounting of responses');
        end
    end
end
