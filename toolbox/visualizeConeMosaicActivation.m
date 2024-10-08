function visualizeConeMosaicActivation(theConeMosaic, theOptics, excitations, modulations, activationLabel)
    
    colsNum = 3;
    rowsNum = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                 'rowsNum', rowsNum, ...
                 'colsNum', colsNum, ...
                 'heightMargin',  0.08, ...
                 'widthMargin',    0.02, ...
                 'leftMargin',     0.03, ...
                 'rightMargin',    0.00, ...
                 'bottomMargin',   0.07, ...
                 'topMargin',      0.05);


    hFig = figure();
    set(hFig, 'Position', [10 10 1400 800]);
    ax = subplot('Position', subplotPosVectors(1,1).v);
    visualizedRangeDegs = max(theConeMosaic.sizeDegs)*0.5*[-1 1 -1 1];
    theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'domain', 'degrees', ...
             'domainVisualizationLimits', visualizedRangeDegs, ...
             'domainVisualizationTicks', struct('x', -0.5:0.1:0.5, 'y', -0.5:0.1:0.5), ...
             'visualizedConeAperture', 'geometricArea', ...
             'crossHairsOnMosaicCenter', true, ...
             'noXLabel', true, ...
             'crossHairsColor', [0 0 0], ...
             'backgroundColor', [0.75 0.75 0.75], ...
             'fontSize', 16, ...
             'plotTitle', 'cone mosaic' ...
             );
    


   

    
    targetWavelength = 500;
    thePSFdata = psfDataFromOI(theOptics, targetWavelength);

    ax = subplot('Position', subplotPosVectors(2,1).v);
    coVisualizeMosaicAndPSF(hFig, ax, theConeMosaic, theOptics, thePSFdata, ...
        visualizedRangeDegs, targetWavelength);

    targetWavelength = 550;
    thePSFdata = psfDataFromOI(theOptics, targetWavelength);

    ax = subplot('Position', subplotPosVectors(2,2).v);
    coVisualizeMosaicAndPSF(hFig, ax, theConeMosaic, theOptics, thePSFdata, ...
        visualizedRangeDegs, targetWavelength);

    targetWavelength = 600;
    thePSFdata = psfDataFromOI(theOptics, targetWavelength);

    ax = subplot('Position', subplotPosVectors(2,3).v);
    coVisualizeMosaicAndPSF(hFig, ax, theConeMosaic, theOptics, thePSFdata, ...
        visualizedRangeDegs, targetWavelength);


    
    
    nTrials = size(excitations,1);

    for iTrial = 1:min([10 nTrials])

            cMap = brewermap(1024, '*greys');
            ax = subplot('Position', subplotPosVectors(1,2).v);
            theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'domain', 'degrees', ...
             'domainVisualizationLimits', visualizedRangeDegs, ...
             'domainVisualizationTicks', struct('x', -0.5:0.1:0.5, 'y', -0.5:0.1:0.5), ...
             'visualizedConeAperture', 'geometricArea', ...
             'crossHairsOnMosaicCenter', true, ...
             'crossHairsColor', [1 1 1], ...
             'activation', excitations(iTrial,:,:), ...
             'verticalActivationColorBarInside', true, ...
             'activationColorMap', cMap, ...
             'activationRange', prctile(excitations(:), [10 90]), ...
             'backgroundColor', cMap(1,:), ...
             'fontSize', 16, ...
             'colorbarFontSize', 12, ...
             'noXLabel', true, ...
             'plotTitle', sprintf('%s (%d) \nexcitations', activationLabel, iTrial) ...
             );


        cMap = brewermap(1024, '*RdBu');
        ax = subplot('Position', subplotPosVectors(1,3).v);

        theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
             'domain', 'degrees', ...
             'domainVisualizationLimits', visualizedRangeDegs, ...
             'domainVisualizationTicks', struct('x', -0.5:0.1:0.5, 'y', -0.5:0.1:0.5), ...
             'visualizedConeAperture', 'geometricArea', ...
             'crossHairsOnMosaicCenter', true, ...
             'crossHairsColor', [1 1 1], ...
             'activation', modulations(iTrial,:,:), ...
             'verticalActivationColorBarInside', true, ...
             'activationColorMap', cMap, ...
             'activationRange', prctile(abs(modulations(:)), 97)*[-1 1], ...
             'backgroundColor', [1 1 1], ...
             'fontSize', 16, ...
             'colorbarFontSize', 12, ...
             'noXLabel', true, ...
             'plotTitle', sprintf('%s (%d) \nmodulations', activationLabel, iTrial) ...
             );
    
    end

end

function coVisualizeMosaicAndPSF(hFig, ax, theConeMosaic, theOI, thePSFdata, visualizedRangeDegs, targetWavelength)

    theConeMosaic.visualize('figureHandle', hFig, ...
        'axesHandle', ax, ...
        'domain', 'degrees', ...
        'visualizedConeAperture', 'geometricArea', ...
        'domainVisualizationLimits', visualizedRangeDegs, ...
        'domainVisualizationTicks', struct('x', -0.50:0.1:0.5, 'y', -0.50:0.1:0.5), ...
        'crossHairsOnMosaicCenter', true, ...
        'labelCones', true, ...
        'plotTitle', sprintf('%d nm', targetWavelength), ...
        'fontSize', 16);

    
    % Add contour plot of the PSF
    hold(ax, 'on');
    cmap = brewermap(10,'greys');
    alpha = 'matchZLevel';
    contourLineColor = [0.0 0.0 0.0];
    visualizedPSF = squeeze(thePSFdata.psf);
    visualizedPSF = visualizedPSF / max(visualizedPSF(:));

    cMosaic.semiTransparentContourPlot(ax, ...
        thePSFdata.supportArcMin/60, thePSFdata.supportArcMin/60, ...
        visualizedPSF, [0.1:0.2:0.99], cmap, alpha, contourLineColor);

    
    % Add horizontal slice through the PSF
    m = (size(visualizedPSF,1)-1)/2+1;
    visualizedPSFslice = squeeze(visualizedPSF(m,:));
    
    xx = thePSFdata.supportArcMin/60;
    yy = min(visualizedRangeDegs)*0.95 + visualizedPSFslice*max(visualizedRangeDegs)*0.5*0.9;

    hL = plot(ax,xx, yy, '-', 'LineWidth', 4.0);
    hL.Color = [1,1,0.8,0.7];
    plot(ax,xx, yy, 'k-', 'LineWidth', 2);
    
end

function thePSFdata = psfDataFromOI(theOI, targetWavelength)

    thePSFdata.supportWavelengthNM = oiGet(theOI, 'wave');
    theOptics = oiGet(theOI, 'optics');
    % Get PSF slice at target wavelength
    thePSFdata.psf = opticsGet(theOptics,'psf data',targetWavelength);

    % Extract support in arcmin
    psfSupportMicrons = opticsGet(theOptics,'psf support','um');
    if (isfield(theOptics, 'micronsPerDegree'))
        micronsPerDegree = theOptics.micronsPerDegree;
    else
        focalLengthMeters = opticsGet(theOptics, 'focalLength');
        focalLengthMicrons = focalLengthMeters * 1e6;
        micronsPerDegree = focalLengthMicrons * tand(1);
    end
    
    xGridMinutes = 60*psfSupportMicrons{1}/micronsPerDegree;
    yGridMinutes = 60*psfSupportMicrons{2}/micronsPerDegree;
    xSupportMinutes = xGridMinutes(1,:);
    ySupportMinutes = yGridMinutes(:,1);

    thePSFdata.supportArcMin = xSupportMinutes ;
end
