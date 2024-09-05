%% Computation of the Lagrange Points of the CR3BP:
function LagrangePoints = LagrangePoints(System)
% LagrangePoints - Computes the Lagrange Points of the CR3BP
% 
% Inputs:
%    System - Structure Containing the System variables
% 
% Outputs:
%    LagrangePoints - Structure Containing the locations of the L1-5
%    Lagrange Points
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

% Xi Axis Root Function:
f = EOMCR3BP(System);
f = @(xi) [0 0 0 1 0 0]*f([xi;0;0;0;0;0],0);

LagrangePoints.L1 = [fsolve(f,.5,optimoptions('fsolve', 'Display', 'off','OptimalityTolerance', 1e-8));0];
LagrangePoints.L2 = [fsolve(f,1.5,optimoptions('fsolve', 'Display', 'off','OptimalityTolerance', 1e-8));0];
LagrangePoints.L3 = [fsolve(f,-1.0,optimoptions('fsolve', 'Display', 'off','OptimalityTolerance', 1e-8));0];
LagrangePoints.L4 = [.5-System.mu;sqrt(3)/2];
LagrangePoints.L5 = [.5-System.mu;-sqrt(3)/2];

end