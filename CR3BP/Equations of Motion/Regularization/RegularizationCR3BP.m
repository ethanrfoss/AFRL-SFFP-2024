%% Regualrization of CR3BP Data
function [s] = RegularizationCR3BP(System,t,x,s0)
% RegularizationCR3BP - Determines Regularized Time from State and time
% data
% 
% Inputs:
%    System - Structure Containing the System variables
%    t - Time Vector
%    x - State Vectors
%    s0 - Initial Regularized Time
% 
% Outputs:
%    s - Regularized time vector
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024


% Initial Time:
if nargin < 4
    s0 = 0;
end

% Integrate Through:
s = zeros(1,length(t));
s(1) = s0;
g = zeros(1,length(t));
g(1) = SundmanTransformation(System,x(:,1));
for i = 2:length(t)
    s(i) = trapz(t(1:i),g(1:i))+s0;
end

end
