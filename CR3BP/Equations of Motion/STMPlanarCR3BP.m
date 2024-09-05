%% State Transition Matrix Equations of Motion for Planar CR3BP:
function Phidot = STMPlanarCR3BP(System,t,x,Phi)
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

% Relative Distances
rho1 = (x(1)+System.mu).^2+x(2).^2;
rho2 = (x(1)-1+System.mu).^2+x(2).^2;

% Reshape State Transition Matrix
Phi = reshape(Phi,4,4);

% Construct A
Uxx = [(System.mu - 1)/rho1^(3/2) - System.mu/rho2^(3/2) + (3*System.mu*(2*System.mu + 2*x(1) - 2)*(System.mu + x(1) - 1))/(2*rho2^(5/2)) - (3*(2*System.mu + 2*x(1))*(System.mu + x(1))*(System.mu - 1))/(2*rho1^(5/2)) + 1,(3*x(2)*System.mu*(System.mu + x(1) - 1))/rho2^(5/2) - (3*x(2)*(System.mu + x(1))*(System.mu - 1))/rho1^(5/2);
        -x(2)*((3*(2*System.mu + 2*x(1))*(System.mu - 1))/(2*rho1^(5/2)) - (3*System.mu*(2*System.mu + 2*x(1) - 2))/(2*rho2^(5/2))),(System.mu - 1)/rho1^(3/2) - System.mu/rho2^(3/2) + x(2)*((3*x(2)*System.mu)/rho2^(5/2) - (3*x(2)*(System.mu - 1))/rho1^(5/2)) + 1];

A = [zeros(2,2),eye(2,2)
     Uxx,[0 2;-2 0]];

% State Transition Matrix Change
Phidot = A*Phi;

% Flatten
Phidot = reshape(Phidot,16,1);

end