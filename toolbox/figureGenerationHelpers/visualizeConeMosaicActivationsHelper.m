function visualizeConeMosaicActivationsHelper(theConeMosaic, ...
        theConeMosaicModulation, theNoisyConeMosaicModulations, ...
         domainVisualizationLimits, ...
         domainVisualizationTicks)

     domainVisualizationTicks.y = [-0.1 0 0.1];
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 4, ...
       'heightMargin',  0.0, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.00);

    hFig = figure(4);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);


    for iOri = 1:4
        ax = subplot('Position', subplotPosVectors(1,iOri).v);
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'activation', theConeMosaicModulation{iOri}, ...
            'activationRange', 20*[-1 1], ...
            'verticalActivationColorBarInside', true, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'noYLabel', (iOri>1), ...
            'plotTitle', sprintf('cone modulation (max: %2.1f%%)',max(abs(theConeMosaicModulation{iOri}(:)))), ...
            'fontSize', 16, ...
            'backgroundColor', [0 0 0]);
        xtickangle(ax, 0);
        ytickangle(ax, 0);
    end

    projectBaseDir = ISETBioJandJRootPath();
    pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'MeanConeMosaicActivations.pdf'];
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);


    hFig = figure(5);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);

    for iInstance = 1:4
        ax = subplot('Position', subplotPosVectors(1,iInstance).v);
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'activation', squeeze(theNoisyConeMosaicModulations{1}(iInstance,:,:)), ...
            'activationRange', 20*[-1 1], ...
            'verticalActivationColorBarInside', true, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'noYLabel', (iOri>1), ...
            'plotTitle', sprintf('noisy cone mosaic activation\n(instance %d)', iInstance),...
            'fontSize', 16, ...
            'backgroundColor', [0 0 0]);
        xtickangle(ax, 0);
        ytickangle(ax, 0);
    end

    projectBaseDir = ISETBioJandJRootPath();
    pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'NoisyConeMosaicActivations.pdf'];
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);


    hFig = figure(6);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);

    ax = subplot('Position', subplotPosVectors(1,4).v);
    theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'plotTitle', 'cone mosaic',...
            'fontSize', 16, ...
            'backgroundColor', [0 0 0]);

    wave = theConeMosaic.wave;
    photopigment = theConeMosaic.pigment;
    macular = theConeMosaic.macular();
    lens = Lens();

    ax = subplot('Position', subplotPosVectors(1,1).v);
    plot(ax, wave, photopigment.quantalEfficiency(:,3), 'bo-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
    xlabel(ax, 'wavelength (nm)');
    ylabel(ax, 'quantal efficiency');
    axis(ax, 'square'); grid(ax, 'on');
    set(ax, 'YLim', [0 0.5], 'YTick', 0:0.1:0.5, 'FontSize', 16)
    title(ax, 'S-cone');

    ax = subplot('Position', subplotPosVectors(1,2).v);
    plot(ax, wave, photopigment.quantalEfficiency(:,2), 'go-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [0.5 1 0.5], 'MarkerEdgeColor', [0.5 0 0], 'LineWidth', 1.5);
    xlabel(ax, 'wavelength (nm)');
    axis(ax, 'square'); grid(ax, 'on');
    set(ax, 'YLim', [0 0.5], 'YTick', 0:0.1:0.5, 'FontSize', 16)
    title(ax, 'M-cone');

    ax = subplot('Position', subplotPosVectors(1,3).v);
    plot(ax, wave, photopigment.quantalEfficiency(:,1), 'ro-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5],  'LineWidth', 1.5);
    xlabel(ax, 'wavelength (nm)');
    axis(ax, 'square'); grid(ax, 'on');
    set(ax, 'YLim', [0 0.5], 'YTick', 0:0.1:0.5, 'FontSize', 16)
    title(ax, 'L-cone');

    projectBaseDir = ISETBioJandJRootPath();
    pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'QuantalEfficiencies.pdf'];
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);


    hFig = figure(100);
    clf;
    set(hFig, 'Position', [10 10 500 1000], 'Color', [1 1 1]);

   
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    plot(ax, wave, log10(photopigment.quantalEfficiency(:,3) .* macular.transmittance .* lens.transmittance), 'bo-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
    hold(ax, 'on');
    plot(ax, wave, log10(photopigment.quantalEfficiency(:,2) .* macular.transmittance .* lens.transmittance), 'go-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [0.5 1 0.5], 'MarkerEdgeColor', [0.5 0 0], 'LineWidth', 1.5);
    
    plot(ax, wave, log10(photopigment.quantalEfficiency(:,1) .* macular.transmittance .* lens.transmittance), 'ro-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5],  'LineWidth', 1.5);
    
    xlabel(ax, 'wavelength (nm)');
    grid(ax, 'on');
    set(ax, 'YLim', [-7 0], 'YTick', -10:1:0, 'FontSize', 16);

    projectBaseDir = ISETBioJandJRootPath();
    pdfFile = fullfile(projectBaseDir, 'figures', 'ConeFundamentalsLog.pdf');
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);
end


