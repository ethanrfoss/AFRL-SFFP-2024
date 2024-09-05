
% Add Path to GMAT
addpath('C:\Users\ethan\OneDrive - Stanford\Documents\GMAT\bin');

% Load GMAT
load_gmat();

% Load Script:
GMATAPI.LoadScript("C:\Users\ethan\OneDrive - Stanford\Documents\GMAT\Instruction Manuals\GMAT\GMAT_examples\CAPSTONE.script");

% Set Output File:
ReportFile = gmat.gmat.GetObject("CAPSTONE_Ephem");
ReportFile.SetField("Filename","C:\Users\ethan\OneDrive - Stanford\Documents\MATLAB\AFOSR SFFP\GMAT\Data\CAPSTONE_EarthMoonSwitch_ephem2.txt")

% Run Script:
ScriptRan = gmat.gmat.RunScript();
if ~ScriptRan
    disp('Script Failed to Run!');
end