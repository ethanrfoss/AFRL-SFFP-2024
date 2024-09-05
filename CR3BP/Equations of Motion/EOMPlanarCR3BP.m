%% Equations of Motion for Planar CR3BP:
function xdot = EOMPlanarCR3BP(System,t,x)
% EOMCR3BP - Equations of Motion for the Planar CR3BP
% 
% Inputs:
%    System - Structure Containing the System variables
%    t - Current Time (Does not matter)
%    x - Current State
% 
% Outputs:
%    xdot - State Derivative
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

rho13 = ((x(1)+System.mu)^2+x(2)^2)^(3/2);
rho23 = ((x(1)+System.mu-1)^2+x(2)^2)^(3/2);

xdot(1:4,1) = [x(3)
             x(4)
             -(1-System.mu)*(x(1)+System.mu)/rho13-System.mu*(x(1)+System.mu-1)/rho23+2*x(4)+x(1)     
             -((1-System.mu)/rho13+System.mu/rho23)*x(2)-2*x(3)+x(2)];

end