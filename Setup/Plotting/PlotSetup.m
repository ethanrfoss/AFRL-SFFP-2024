%% Plotter Setup
% Function Sets the Plot Settings
function S = PlotSetup(S)

% Latex Plot Settings:
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% Set Font:
set(groot,'DefaultTextFontname','CMU Serif');
set(groot,'DefaultAxesFontname','CMU Serif');

end