%% JacobiConstant - Jacobi Constant for the CR3BP
% 
% Inputs:
%    System - Structure Containing the System variables
%    x - Current State
% 
% Outputs:
%    C - Jacobian Constant
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function C = JacobiConstant(System,x)

% Adjust for Possible Planarity:
if size(x,1) == 4
    x = [x(1:2,:);zeros(1,size(x,2));x(3:4,:);zeros(1,size(x,2))];
end

% Loop Through to compute JacobiConstant:
C = zeros(1,size(x,2));
for i = 1:size(x,2)
    
    % Distances:
    r1 = ((x(1,i)+System.mu)^2+x(2,i)^2+x(3,i)^2)^(1/2);
    r2 = ((x(1,i)+System.mu-1)^2+x(2,i)^2+x(3,i)^2)^(1/2);

    % Jacobi:
    U = (1-System.mu)/r1 + System.mu/r2 + 1/2*(x(1,i)^2+x(2,i)^2);
    C(i) = 2*U - norm(x(4:6,i))^2;

end

end