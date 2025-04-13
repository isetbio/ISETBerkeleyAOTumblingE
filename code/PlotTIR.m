
%% Load in the results for which you want to plot the TIR

%% Set title string
titleStr = 'Photocurrent';
titleStr = 'Watson Filter';

% Normalizer
normalizer = max(  theNeuralEngine.noiseFreeComputeParams.temporalFilter.filterValues(:));

%%
figure; clf; hold on;
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
    case 'Watson Filter'
        ylim([-1 1]);
end
