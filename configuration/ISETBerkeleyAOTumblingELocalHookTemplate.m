
%% Define project
projectName = 'ISETBerkeleyAOTumblingE';

%% Delete any old prefs
if (ispref(projectName))
    rmpref(projectName);
end

%% Specify base paths for materials and data.
%
% The switch statement below allows definition of
% user/machine specific directory locations. You may
% not have all the routines that are called to make
% this work, but in the end all you need to do is
% make sure the locations are set correctly for your
% machine.
[~, userID] = system('whoami');
userID = strtrim(userID);

% Define the base dir for everything in the project. In the
% brainard lab, this starts at the lab dropbox tree, and the
% code below finds that on any Mac OS machine.  You can
% hardcode baseDir to be something else on your machine, if you
% want.  And, you don't really even need it as long as you get
% the locations set below to be apprpriate for your machine.
switch userID  
    otherwise
        if ismac
            dbJsonConfigFile = '~/.dropbox/info.json';
            fid = fopen(dbJsonConfigFile);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            val = jsondecode(str);
            baseDir = val.business.path;
        end
end

% Set 'dropboxPath' preference to the top level directory where files will be
% stored.  In our lab, it is inside of dropbox so we always begin
% by pointing to that.  But there is nothing magic about having
% it available under dropbox.
setpref(projectName,'dropboxPath',baseDir);

% Set project path 
setpref(projectName,'dataPath',fullfile(baseDir,projectName));