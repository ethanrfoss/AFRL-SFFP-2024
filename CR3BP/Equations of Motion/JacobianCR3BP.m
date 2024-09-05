%% Jacobian of the CR3BP:
function A = JacobianCR3BP(System,Options)
% JacobianCR3BP - Jacobian for the CR3BP
% 
% Inputs:
%    System - Structure Containing the System variables
%    t - Current Time (Does not matter)
%    x - Current State
% 
% Outputs:
%    A - Evaluated jacobian
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

% Options:
DefaultOptions.Regularize = false;
DefaultOptions.Planar = false;
if nargin < 2
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

if Options.Planar
    A = @(x,t) JacobianCR3BPPlanar(System,x);
else
    A = @(x,t) JacobianCR3BPNonPlanar(System,x);
end

% Regularization:
if Options.Regularize
    g = @(x) SundmanTransformation(System,x);
    dgdx = @(x) [(3*(2*System.mu + 2*x(1))*(System.mu - 1))/(4*((System.mu + x(1))^2 + x(2)^2 + x(3)^2)^(7/4)) - (3*System.mu*(2*System.mu + 2*x(1) - 2))/(4*((System.mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(7/4)), (3*x(2)*(System.mu - 1))/(2*((System.mu + x(1))^2 + x(2)^2 + x(3)^2)^(7/4)) - (3*x(2)*System.mu)/(2*((System.mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(7/4)), (3*x(3)*(System.mu - 1))/(2*((System.mu + x(1))^2 + x(2)^2 + x(3)^2)^(7/4)) - (3*System.mu*x(3))/(2*((System.mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(7/4)), 0, 0, 0];
    f = EOMCR3BP(System,Options);
    A = @(x,t) 1/g(x)*A(x,t)-1/g(x)^2*f(x,t)*dgdx(x);
end

end

function A = JacobianCR3BPPlanar(System,x)

% Relative Distances (squared)
r1 = (x(1)+System.mu).^2+x(2).^2;
r2 = (x(1)-1+System.mu).^2+x(2).^2;

% Construct A
Uxx = [(System.mu - 1)/r1^(3/2) - System.mu/r2^(3/2) + (3*System.mu*(2*System.mu + 2*x(1) - 2)*(System.mu + x(1) - 1))/(2*r2^(5/2)) - (3*(2*System.mu + 2*x(1))*(System.mu + x(1))*(System.mu - 1))/(2*r1^(5/2)) + 1,(3*x(2)*System.mu*(System.mu + x(1) - 1))/r2^(5/2) - (3*x(2)*(System.mu + x(1))*(System.mu - 1))/r1^(5/2)
-x(2)*((3*(2*System.mu + 2*x(1))*(System.mu - 1))/(2*r1^(5/2)) - (3*System.mu*(2*System.mu + 2*x(1) - 2))/(2*r2^(5/2))),(System.mu - 1)/r1^(3/2) - System.mu/r2^(3/2) + x(2)*((3*x(2)*System.mu)/r2^(5/2) - (3*x(2)*(System.mu - 1))/r1^(5/2)) + 1];

A = [zeros(2,2),eye(2,2)
     Uxx,[0 2;-2 0]];

end

function A = JacobianCR3BPNonPlanar(System,x)

% Relative Distances (squared)
r1 = (x(1)+System.mu).^2+x(2).^2+x(3).^2;
r2 = (x(1)-1+System.mu).^2+x(2).^2+x(3).^2;

% Construct A
Uxx = [(System.mu - 1)/r1^(3/2) - System.mu/r2^(3/2) + (3*System.mu*(2*System.mu + 2*x(1) - 2)*(System.mu + x(1) - 1))/(2*r2^(5/2)) - (3*(2*System.mu + 2*x(1))*(System.mu + x(1))*(System.mu - 1))/(2*r1^(5/2)) + 1,(3*x(2)*System.mu*(System.mu + x(1) - 1))/r2^(5/2) - (3*x(2)*(System.mu + x(1))*(System.mu - 1))/r1^(5/2),(3*System.mu*x(3)*(System.mu + x(1) - 1))/r2^(5/2) - (3*x(3)*(System.mu + x(1))*(System.mu - 1))/r1^(5/2)
-x(2)*((3*(2*System.mu + 2*x(1))*(System.mu - 1))/(2*r1^(5/2)) - (3*System.mu*(2*System.mu + 2*x(1) - 2))/(2*r2^(5/2))),(System.mu - 1)/r1^(3/2) - System.mu/r2^(3/2) + x(2)*((3*x(2)*System.mu)/r2^(5/2) - (3*x(2)*(System.mu - 1))/r1^(5/2)) + 1,x(2)*((3*System.mu*x(3))/r2^(5/2) - (3*x(3)*(System.mu - 1))/r1^(5/2))
-x(3)*((3*(2*System.mu + 2*x(1))*(System.mu - 1))/(2*r1^(5/2)) - (3*System.mu*(2*System.mu + 2*x(1) - 2))/(2*r2^(5/2))),x(3)*((3*x(2)*System.mu)/r2^(5/2) - (3*x(2)*(System.mu - 1))/r1^(5/2)), x(3)*((3*System.mu*x(3))/r2^(5/2) - (3*x(3)*(System.mu - 1))/r1^(5/2)) + (System.mu - 1)/r1^(3/2) - System.mu/r2^(3/2)];

A = [zeros(3,3),eye(3,3)
     Uxx,[0 2 0;-2 0 0;0 0 0]];

end

