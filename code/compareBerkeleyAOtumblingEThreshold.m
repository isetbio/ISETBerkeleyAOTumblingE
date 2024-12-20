%% Visualize AO Tumbling E Scene for checking
t_BerkeleyAOtumblingESceneEngine('visualizeScene', true, ...
    'temporalModulationParams_numFrame', 8, ...
    'temporalModulationParams_xShiftPerFrame', [0 0 0 0 10/60 0 0 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0 0 0 0 0 0], ...
    'temporalModulationParams_backgroundRGBPerFrame', [0 0 0; 0 0 0; 0 0 0; 1 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0], ...
    'responseFlag', 'excitation')

t_BerkeleyAOtumblingESceneEngine('visualizeScene', true, ...
    'temporalModulationParams_numFrame', 2, ...
    'temporalModulationParams_xShiftPerFrame', [0 10/60], ...
    'temporalModulationParams_yShiftPerFrame', [0 0], ...
    'temporalModulationParams_backgroundRGBPerFrame', [1 0 0; 1 0 0], ...
    'responseFlag', 'excitation')

%% photopigment excitation
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 8, ...
    'temporalModulationParams_xShiftPerFrame', [0 0 0 0 10/60 0 0 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0 0 0 0 0 0], ...
    'temporalModulationParams_backgroundRGBPerFrame', [0 0 0; 0 0 0; 0 0 0; 1 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0], ...
    'responseFlag', 'excitation');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 2, ...
    'temporalModulationParams_xShiftPerFrame', [0 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0], ...
    'responseFlag', 'excitation', 'exportCondition', 'no change');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 2, ...
    'temporalModulationParams_xShiftPerFrame', [0 10/60], ...
    'temporalModulationParams_yShiftPerFrame', [0 0], ...
    'responseFlag', 'excitation','exportCondition', 'shift 1');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 2, ...
    'temporalModulationParams_xShiftPerFrame', [0 20/60], ...
    'temporalModulationParams_yShiftPerFrame', [0 0], ...
    'responseFlag', 'excitation','exportCondition', 'shift 2');

%% photocurrent
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 0 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'photocurrent','exportCondition', 'no change');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 10/60 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'photocurrent','exportCondition', 'shift 1');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 20/60 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'photocurrent','exportCondition', 'shift 2');