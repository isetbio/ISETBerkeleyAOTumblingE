% photopigment excitation
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 0 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'excitation');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 10/60 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'excitation');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 20/60 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'excitation');

% photocurrent
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 0 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'photocurrent');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 10/60 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'photocurrent');
BerkeleyAOtumblingEThreshold('temporalModulationParams_numFrame', 3, ...
    'temporalModulationParams_xShiftPerFrame', [0 20/60 0], ...
    'temporalModulationParams_yShiftPerFrame', [0 0 0], ...
    'responseFlag', 'photocurrent');