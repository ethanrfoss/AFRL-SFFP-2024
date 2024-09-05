%% Gauss Newton Method - Gauss Newton Root Finding
% 
% Inputs:
%    f - Function to find the root of
%    x0 - Initial Guess of the root
%    Options - Structure Array for Gauss Newton Options. The Options and 
%    their default values are:
%        Maximum Iteration (50) - Maximum Iterations of Descent
%        Tolerance (10^-6) - Root Tolerance
%        ConditionThreshold (10^10) - Breaking Threshold for Matrix
%        Condition
%        Printing (false)
%        MaxTol (1000) - Breaking Tolerance Threshold
% 
% Outputs:
%    x - Converged root
%    Converged - Convergence boolean
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function [x,Converged] = GaussNewton(f,x0,Options)

% Set up Options
DefaultOptions.MaximumIterations = 50;
DefaultOptions.Tolerance = 10^-6;
DefaultOptions.CheckCondition = false;
DefaultOptions.ConditionThreshold = 10^10;
DefaultOptions.Printing = false;
DefaultOptions.MaxTol = 1000;
if nargin < 3
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

% Initialize:
x = x0;
i = 0;
Converged = false;

% Loop:
while i < Options.MaximumIterations
    
    % Evaluate
    [fval,Jval] = f(x);

    % Check Convergence
    if norm(fval) < Options.Tolerance
        Converged = true;
        if Options.Printing
            disp(['Gauss Newton Converged. Iteration: ' num2str(i) ' Current Tolerance: ' num2str(norm(fval))])
        end
        break;
    end

    if norm(fval) > Options.MaxTol
        Converged = false;
        if Options.Printing
            disp('Gauss Newton Unconverged, Maximum Tolerance Reached');
        end
        break;
    end

    % Print:
    if Options.Printing
        disp(['Gauss Newton Iteration: ' num2str(i) ' Current Tolerance: ' num2str(norm(fval))]);
    end

    % Condition:
    if Options.CheckCondition && cond(Jval) >= Options.ConditionThreshold
        error('Condition Number Threshold Exceeded');
    end

    % Compute Step:
    dx = -Jval\fval;
    
    % Update:
    x = x + dx;

    % Iterate:
    i = i + 1;

end

end