%save_model_workspace_20120425.m

%written by the Senger Research Group, 4/25/12

%this program saves only the variables specified to the workspace
%specified.  This is handy for doing quick saves.

%this program also clears the workspace of un-saved variables.

path1='/Users/rsenger/Documents/Research VT/TCU/E_coli_isobutanol';
file1='model workspace.mat';

save(fullfile(path1,file1),'model','t','t1','r','FBAsolution');

clearvars -except file1 path1

load(fullfile(path1,file1));

clear file1 path1

fprintf(['Workspace saved.' '\n']);

