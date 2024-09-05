%% Levenberg Marquadt Method
function [x,Converged] = LevenbergMarquadt(f,x0,Options)
% LevenbergMarquadt - Levenberg Marquadt Root Finding
% 
% Inputs:
%    f - Function to find the root of
%    x0 - Initial Guess of the root
%    Options - Structure Array for Gauss Newton Options. The Options and 
%    their default values are:
%        Maximum Iteration (50) - Maximum Iterations of Descent
%        Tolerance (10^-6) - Root Tolerance
%        ConditionThreshold (10^10) - Breaking Threshold for Matrix
%        CheckCondition (false) - Condition Check
%        betash (10) - Shrinking Parameter for Trust Region
%        betagr (10) - Growth Parameter for Trust Region
%        lambda (1) - Initial Trust Region Size
%        Printing (false)
%        MaxLambda (10^8) - Breaking Thresold for Trust Region
% 
% Outputs:
%    x - Converged root
%    Converged - Convergence boolean
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

% Set up Options
DefaultOptions.MaximumIterations = 50;
DefaultOptions.Tolerance = 10^-6;
DefaultOptions.ConditionThreshold = 10^10;  
DefaultOptions.lambda = 1;
DefaultOptions.betash = 10;
DefaultOptions.betagr = 10;
DefaultOptions.Printing = false;
DefaultOptions.CheckCondition = false;
DefaultOptions.MaxLambda = 10^8;
if nargin < 3
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

% Initialize:
x = x0;
i = 0;
Converged = false;
lambda = Options.lambda;

% Evaluate
[fval,Jval] = f(x);

% Loop:
while i < Options.MaximumIterations
    

    % Normal Equations
    if size(Jval,1) >= size(Jval,2)
        A = Jval'*Jval;
        b = -Jval'*fval;
        C = eye(size(Jval,2));
    else
        C = Jval';
        A = Jval*Jval';
        b = -fval;
    end

    % Check Condition Number
    if Options.CheckCondition && cond(A + lambda*diag(diag(A))) > Options.ConditionThreshold
        break;
    end

    % Compute Step:
    dx = C*((A + lambda*diag(diag(A)))\b);
    
    % Update:
    xnew = x + dx;
    [fnew,Jnew] = f(xnew); 

    % Check Convergence
    if norm(fnew) < Options.Tolerance
        Converged = true;
        x = xnew;
        if Options.Printing
            disp(['Levenberg Marquadt Converged. Iteration: ' num2str(i) ' Current Tolerance: ' num2str(norm(fnew)) ' Trust Region Size: ' num2str(lambda)])
        end
        break;
    end

    % Trust Region Update:
    if norm(fnew) < norm(fval)
        lambda = lambda/Options.betash;
        x = xnew;
        fval = fnew;
        Jval = Jnew; 
    else
        lambda = lambda*Options.betagr;
    end

    if lambda >= Options.MaxLambda
        break;
    end

    % Iterate:
    i = i + 1;

    % Print:
    if Options.Printing
        disp(['Levenberg Marquadt Iteration: ' num2str(i) ' Current Tolerance: ' num2str(norm(fval)) ' Trust Region Size: ' num2str(lambda)]);
    end

end

end