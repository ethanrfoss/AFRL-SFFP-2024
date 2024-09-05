%% SundmanTransformation - Regularization Value for the CR3BP
% 
% Inputs:
%    System - Structure Containing the System variables
%    x - Current State
% 
% Outputs:
%    g - Sundman Value
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function g = SundmanTransformation(System,x)

% Primary Distances:
r1 = ((x(1)+System.mu)^2+x(2)^2+x(3)^2)^(1/2); % Distance to First Primary
r2 = ((x(1)+System.mu-1)^2+x(2)^2+x(3)^2)^(1/2); % Distance to Second Primary

% Sundman Transform:
g = (1-System.mu)/r1^(3/2) + System.mu/r2^(3/2);

end