% The core code is in, well, the code subdirectory.

% Run BerkeleyAOTumblingECalcsIntermixed for our model of the experment,
% where all E orientations and directions were run intermixed.  Use
% 'Calcs6' as set.  This takes a long time, and can cause Matlab to crash
% even on quite a large machine (M2 Studio with > 192 GB ram).  It writes
% output into the directory specified by the preferences.

% Run BerkeleyAOTumblingECalcsIntermixedAnalyze to produce the key figures.
% This will also have the data in the form that the figures use to plot it,
% so you can snag that manually if desired.  Use 'Calcs6' and match
% nReplications to the number you computed for.
%
% TODO: Save plots for paper.  Save out tabular data.

% Run BerkeleyAOTumblingECalcsVisualize to produce some movies of the
% stimulus.
%
% TODO: Check that this still runs and document what it produces.

% Run PlotTIR for a plot of the temporal impulse response.
%
% TODO: Check that this still runs and document what it produces.

% Run validation/ieValidateISETBerkeleyAOTumblingE to run some basic
% validations on the operation of the code.  Checks that the examples
% provided in core routine BerkeleyAOTumblingEThreshold run without
% crashing, and in one case that it produces the same answer that it did
% when we developed the code. These passed in March 2026.  You will need
% the full ISETBio suite setup for this, including isetvalidate.

% See the configuration subdirectory for dependencies in ToolboxToolbox json
% format, and a local hook template that shows the preferences you need to
% set up.  See github.com/toolboxhub/toolboxtoolbox.git for the
% ToolboxToolbox.