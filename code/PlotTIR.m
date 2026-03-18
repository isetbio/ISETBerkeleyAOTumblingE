%% TIR plots
%
% Data filenames are somewhat hard coded, but easy enough to adjust

%% Clear and close
clear; close all;

%% Set title string
titleStrs = {'Photocurrent' 'Watson Filter'};
calcName = 'Calcs6';

%% Point at data
rootPath = getpref('ISETBerkeleyAOTumblingE','dataPath');
resultsDir = 'results';
figsDir = 'figures';

%% Load in some data
%
% Any offset will do, just need to match the filter number
for ss = 1:length(titleStrs)
    titleStr = titleStrs{ss};
    switch (titleStr)
        case 'Photocurrent'
            load(fullfile(rootPath,resultsDir, ...
                ['/BerkeleyAOTumblingEThreshold_' calcName '_Intermixed_Rep_1_filter_2/BerkeleyAOTumblingEThreshold_' ...
                calcName '_Intermixed_Rep_1_filter_2.mat']));
        case 'Watson Filter'
            load(fullfile(rootPath,resultsDir, ...
                ['/BerkeleyAOTumblingEThreshold_' calcName '_Intermixed_Rep_1_filter_3/BerkeleyAOTumblingEThreshold_' ...
                calcName '_Intermixed_Rep_1_filter_3.mat']));
    end

    % Normalizer
    normalizer = max(theNeuralEngine.noiseFreeComputeParams.temporalFilter.filterValues(:));

    %% Plot and save
    theFig = figure; clf; hold on;
    set(gca,'FontName','Helvetica','FontSize',14);
    plot(theNeuralEngine.noiseFreeComputeParams.temporalFilter.temporalSupport, ...
        theNeuralEngine.noiseFreeComputeParams.temporalFilter.filterValues/normalizer, ...
        'r','LineWidth',3);
    plot([0 theNeuralEngine.noiseFreeComputeParams.temporalFilter.temporalSupport(end)],[0 0],'k:','LineWidth',2);
    xlabel('Time (sec)','FontSize',18);
    ylabel('Filter Response (arb units)','FontSize',18);
    title(titleStr,'FontSize',20);
    switch (titleStr)
        case 'Photocurrent'
            ylim([-1 1]);
            print(theFig,fullfile(rootPath,figsDir,[calcName '_PhotocurrentTIR.tiff']),'-dtiff','-r300');

        case 'Watson Filter'
            ylim([-1 1]);
            print(theFig,fullfile(rootPath,figsDir,[calcName '_WatsonFilterTIR.tiff']),'-dtiff','-r300');
    end

end
