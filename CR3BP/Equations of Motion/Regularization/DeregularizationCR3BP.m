%% Deregularization of CR3BP Data
function [t] = DeregularizationCR3BP(System,s,x,t0)
% DeregularizationCR3BP - Determines Standard Time from State and 
% Regularized time data
% 
% Inputs:
%    System - Structure Containing the System variables
%    t - Regularized Time Vector
%    x - State Vectors
%    t0 - Initial Time
% 
% Outputs:
%    t - Time vector
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024


% Initial Time:
if nargin < 4
    t0 = 0;
end

% Integrate Through:
t = zeros(1,length(s));
t(1) = t0;
g = zeros(1,length(t));
g(1) = SundmanTransformation(System,x(:,1));
for i = 2:length(s)
    t(i) = trapz(s(1:i),1./g(1:i))+t0;
end

end
