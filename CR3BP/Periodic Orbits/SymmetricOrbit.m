%% SymmetricOrbit - Determines a symettric orbit about the x-z plane
% crossing
%
% Inputs: 
%    System - Structure Containing the System variables, which is
%    constructed using System = InitializeCR3BP(System).
%    InitialGuess - Initial Guess structure for the periodic orbit. It must
%    contain:
%       x - Initial Guess for the state at the x-z plane crossing
%       P - Initial Guess for the period of the orbit
%    Options - Structure Containing Orbit Options. The Options and their
%    default values are:
%       Plot (true) - Plot the periodic orbit
%       PlotEigen (false) - Plot the eigenvalues of the orbit on the
%       complex plane
%       Save (false) - Save the resulting periodic orbit
%       Regularize (false) - Regularize the dynamics
%       N (1) - Shooting Arcs
%       Shooting - Structure Containing shooting method options:
%           AbsoluteTolerance (10^-15) - Absolute Integration Tolerance
%           RelativeTolerance (10^-12) - Relative Integration Tolerance
%           Tolerance (10^-10) - Solution Tolerance
%           Printing (false) - Print Diagnostics
%           RootFinder (GaussNewton) - Root finding method
%       SaveName (Optional) - Save Name for the Orbit. If none is specified
%       the save name is prompted
%       Axis (Optional) - Axis of the Plot. If none is specified a default
%       axis is set 
% Outputs:
%    f - Figure handle
%
% Example Usage:
%    System = CislunarSystem;
%    System = InitializeCR3BP(System);
%    InitialGuess.P = 2.726;
%    InitialGuess.x = [.8255;0;0;.1048];
%    SymmetricOrbit(System,InitialGuess)
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function Orbit = SymmetricOrbit(System,InitialGuess,Options)

% Options:
DefaultOptions.Plot = true; % Plotting
DefaultOptions.PlotEigen = false;
DefaultOptions.PlotInertial = false;
DefaultOptions.Save = false; 
DefaultOptions.SaveDirectory = fileparts(mfilename('fullpath'));
DefaultOptions.Regularize = false;
DefaultOptions.N = 1; % Shooting arc
DefaultShooting.AbsoluteTolerance = 10^-15;
DefaultShooting.RelativeTolerance = 10^-12;
DefaultShooting.Tolerance = 10^-10;
DefaultShooting.MaximumIterations = 200;
DefaultShooting.Printing = false;
DefaultShooting.RootFinder = 'GaussNewton';
DefaultOptions.Shooting = DefaultShooting;
DefaultOptions.FreeTime = false;
if nargin < 3
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.Planar = length(InitialGuess.x) == 4;
% TPBVP.nx = length(InitialGuess.x);

% Setup TPBVP:
TPBVP = SymmetricOrbitTPBVP(System,Options);

% Construct Initial Guess
% if Options.Regularize
%     % Convert Time:
%     [t,x] = ode89(@(t,x)EOMCR3BP(System,t,x,Options),[0 InitialGuess.P/2],InitialGuess.x,odeset('AbsTol',Options.Shooting.AbsoluteTolerance,'RelTol',Options.Shooting.RelativeTolerance));
%     [s] = RegularizationCR3BP(System,t,x',0);
%     InitialGuess.P = s(end)*2; % Regularized Period
% end
if Options.N == 1
    InitialGuess.T = InitialGuess.P/2;
else
    InitialGuess.T = ones(1,Options.N)*InitialGuess.P/(2*Options.N);
    [~,InitialGuess.x] = ode89(@(t,x) TPBVP.f(x,t),[0 cumsum(InitialGuess.T(1:end-1))],InitialGuess.x,odeset('RelTol',Options.Shooting.RelativeTolerance,'AbsTol',Options.Shooting.AbsoluteTolerance));
    InitialGuess.x = InitialGuess.x';
    if Options.N == 2
        InitialGuess.x = [InitialGuess.x(:,1) InitialGuess.x(:,end)];
    end
end
InitialGuess.p = [];

% Solve the shooting problem
[Solution] = Shooting(System,TPBVP,InitialGuess,Options.Shooting);
Orbit.ShootingSolution = Solution;
Orbit.Converged = Solution.Converged;
if ~Orbit.Converged
    Orbit.t = [];
    Orbit.x = [];
    Orbit.Period = [];
    Orbit.STM = [];
    Orbit.Monodromy = [];
    Orbit.EigenVectors = [];
    Orbit.EigenValues = [];
    Orbit.StabilityIndex = [];
    Orbit.JacobiConstant = [];
    Orbit.InertialState = [];
    return;MMR
end

% Integrate for full period:
Orbit.Period = 2*Solution.t(end);
[Orbit.t,xphi] = ode89(@(t,xphi) [TPBVP.f(xphi(1:TPBVP.nx),t);reshape(TPBVP.A(xphi(1:TPBVP.nx),t)*reshape(xphi(TPBVP.nx+1:TPBVP.nx+TPBVP.nx^2),[TPBVP.nx,TPBVP.nx]),[TPBVP.nx^2,1])],[0 Orbit.Period],[Solution.x(:,1);reshape(eye(TPBVP.nx),[TPBVP.nx^2,1])],odeset('RelTol',Options.Shooting.RelativeTolerance,'AbsTol',Options.Shooting.AbsoluteTolerance));
Orbit.x = xphi(:,1:TPBVP.nx)';
Orbit.STM = reshape(xphi(:,TPBVP.nx+1:TPBVP.nx+TPBVP.nx^2)',[TPBVP.nx,TPBVP.nx,length(Orbit.t)]);
Orbit.Monodromy = Orbit.STM(:,:,end);

% Determine Eigenstructure and Stability Index
[V,E] = eigs(Orbit.Monodromy);
Orbit.EigenValues = diag(E);
Orbit.EigenVectors = V;
Orbit.StabilityIndex = 1/2*(abs(Orbit.EigenValues)+1./abs(Orbit.EigenValues));
% Orbit.StabilityIndex = uniquetol(Orbit.StabilityIndex,10^-8);
Orbit.StabilityIndex = max(Orbit.StabilityIndex);

% Inertial State:
Orbit.InertialState = CR3BPFrameTransform(System,Orbit.t,Orbit.x,["Barycenter","Earth"],["Rotating","Inertial"]);

% Jacobi Constant:
Orbit.JacobiConstant = JacobiConstant(System,Orbit.x(:,1));

% Plotting
if Options.Plot
    f = TabFigures(CR3BPOrbitFigure(System,Orbit,Options), "Rotating Frame Orbit");
end
if Options.PlotInertial
    f = TabFigures(InertialOrbitFigure(System,Orbit,Options), "Inertial Frame Orbit",f);
end

% Saving Data and Figure:
if Options.Save
    SaveData(struct('Orbit',Orbit),Options);
    SaveFigure(f,Options);
end

end