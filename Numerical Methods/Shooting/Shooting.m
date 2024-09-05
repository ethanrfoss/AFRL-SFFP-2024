%% Shooting Method:
function [Solution] = Shooting(System,TPBVP,InitialGuess,Options)
% LevenbergMarquadt - Levenberg Marquadt Root Finding
% 
% Inputs:
%    System - Structure of System Parameters
%    TPBVP - Structure of TPBVP Problem Functions. These function handles
%    are required:
%        f(x) - Dynamics Function 
%        g(x0,xf,p,T) - Boundary Conditions
%    These functions can be optionally specified. If they are not
%    specified, then they are computed:
%        A(x) - Dynamics Jacobian
%        H0(x0,xf,p,T) - Initial Value Boundary Jacobian
%        Hf(x0,xf,p,T) - Final Value Boundary Jacobian
%        I(x0,xf,p,T) - Parameter Boundary Jacobian
%        Jf(x0,xf,p,T) - Time Boundary Jacobian
%    InitialGuess - Initial Guess structure for the shooting method. The
%    structure requires:
%        x - Vectors of States, corresponding to initial conditions for
%        each shooting arc
%        p - Parameters
%        T - Vector of Time segments, corresponding to each shooting arc
%    Options - Structure Array for Gauss Newton Options. The Options and 
%    their default values are:
%        Maximum Iteration (50) - Maximum Iterations
%        Tolerance (10^-12) - Root Tolerance
%        AbsoluteTolerance (10^-15) - Integration Absolute Tolerance
%        RelativeTolerance (10^-12) - Integration Relative Tolerance
%        FreeTime (false) - Free or Fixed Final Time problem
%        RootFinder (LevenbergMarquadt) - Root Finder for the shooting
%        method. The options are GaussNewton, LevenbergMarquadt, or fsolve
% 
% Outputs:
%    Solution - Solution Structure of the problem which includes:
%        t - Time Vector of continuous solution
%        x - State vector of continuous solution
%        Converged - Convergence boolean
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

% Set up Options
DefaultOptions.MaximumIterations = 20;
DefaultOptions.Tolerance = 10^-12;
DefaultOptions.AbsoluteTolerance = 10^-15;
DefaultOptions.RelativeTolerance = 10^-12;
DefaultOptions.FreeTime = false;
DefaultOptions.RootFinder = 'LevenbergMarquadt';
if nargin < 4
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end 

% Compute Jacobians if Necessary
if ~isfield(TPBVP,'A')
    [TPBVP] = ComputeTPBVP(System,TPBVP);
end

% Shooting Arc Count:
TPBVP.N = length(InitialGuess.T);
if TPBVP.N ~= size(InitialGuess.x,2)
    error('Incommensurate Arc Count');
end
InitialGuess.x = reshape(InitialGuess.x,[TPBVP.N*TPBVP.nx,1]);
InitialGuess.T = reshape(InitialGuess.T,[TPBVP.N,1]);

% Find Root:
if Options.FreeTime
    if isequal(Options.RootFinder,'LevenbergMarquadt')
        [z0,Converged] = LevenbergMarquadt(@(z)ShootingRootFunction(z,TPBVP,Options),[InitialGuess.x;InitialGuess.p;InitialGuess.T],Options);
    elseif isequal(Options.RootFinder,'GaussNewton')
        [z0,Converged] = GaussNewton(@(z)ShootingRootFunction(z,TPBVP,Options),[InitialGuess.x;InitialGuess.p;InitialGuess.T],Options);
    elseif isequal(Options.RootFinder,'fsolve')
        [z0,~,flag] = fsolve(@(z)ShootingRootFunction(z,TPBVP,Options),[InitialGuess.x;InitialGuess.p;InitialGuess.T],optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',Options.MaximumIterations));
        Converged = flag > 0;
    else
        error('Invalid Root Finder Specified');
    end
else
    if isequal(Options.RootFinder,'LevenbergMarquadt')
        [z0,Converged] = LevenbergMarquadt(@(z)ShootingRootFunction(z,TPBVP,Options,InitialGuess.T),[InitialGuess.x;InitialGuess.p],Options);
    elseif isequal(Options.RootFinder,'GaussNewton')
        [z0,Converged] = GaussNewton(@(z)ShootingRootFunction(z,TPBVP,Options,InitialGuess.T),[InitialGuess.x;InitialGuess.p],Options);
    elseif isequal(Options.RootFinder,'fsolve')
        [z0,~,flag] = fsolve(@(z)ShootingRootFunction(z,TPBVP,Options,InitialGuess.T),[InitialGuess.x;InitialGuess.p],optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','off','MaxIterations',Options.MaximumIterations));
        Converged = flag > 0;
    else
        error('Invalid Root Finder Specified');
    end
end

% Integrate to Get Trajectory
if Options.FreeTime
    [Solution.t,Solution.x] = ode89(@(t,x)TPBVP.f(x),[0 sum(z0(end-TPBVP.N+1:end))],z0(1:TPBVP.nx),odeset('AbsTol',Options.AbsoluteTolerance,'RelTol',Options.RelativeTolerance));
else
    [Solution.t,Solution.x] = ode89(@(t,x)TPBVP.f(x),[0 sum(InitialGuess.T)],z0(1:TPBVP.nx),odeset('AbsTol',Options.AbsoluteTolerance,'RelTol',Options.RelativeTolerance));
end
Solution.x = Solution.x';
Solution.Converged = Converged;

end