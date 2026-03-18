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
% This also saves out the data that was plotted in three tables, one for
% each of the plots.  I have not checked in detail that these match up with
% the plots, but I spotted checked a few numbers that do. If you use these,
% plot from the tables and verify matching up with the Matlab generated
% plots.

% 2026-03-18, DHB: Reran the calculations above (4 replications) and
% verfied that the plots have the key features as the data run and saved in
% June 2025 with 16 replications, with slight differences I attribute to
% run-to-run noise. Regenerated videos. So everything seems stable and replicable.

% Run BerkeleyAOTumblingECalcsVisualize to produce some movies of the
% stimulus. These come out in directories with Visualize6 in their names,
% under the same results and figures directories as the main calculations.
% There are movies both with and without noise, for each filter and spatial
% displacement.

% Run PlotTIR for plots of the temporal impulse responses (photocurrent and Watson).

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