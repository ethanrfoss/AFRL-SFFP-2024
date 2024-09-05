%% State Transition Matrix Equations of Motion for CR3BP:
function Phidot = STMCR3BP(System,t,x,Phi,Options)
% STMCR3BP - Equations of Motion for State Transition Matrix Evolution for
% the CR3BP
% 
% Inputs:
%    System - Structure Containing the System variables
%    t - Current Time (Does not matter)
%    x - Current State
%    Phi - Current State Transition Matrix (flattened)
% 
% Outputs:
%    Phidot - Flattened STM Derivative
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

% Options:
DefaultOptions.Regularize = false;
DefaultOptions.Planar = false;
if nargin < 5
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

% Jacobian
A = JacobianCR3BP(System,t,x,Options);

% Reshape State Transition Matrix
if Options.Planar
    Phi = reshape(Phi,4,4);
else
    Phi = reshape(Phi,6,6);
end

% State Transition Matrix Change
Phidot = A*Phi;

% Flatten
if Options.Planar
    Phidot = reshape(Phidot,16,1);
else
    Phidot = reshape(Phidot,36,1);
end

end