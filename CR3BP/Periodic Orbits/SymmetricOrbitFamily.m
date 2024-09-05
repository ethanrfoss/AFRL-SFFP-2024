%% SymmetricOrbitFamily - Determines a family of symettric orbit about the 
% x-z plane crossing
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
%       Nm (10) - Number of Backwards orbit continuations
%       Np (10) - Number of Forwards orbit continuations
%       Plot (true) - Plot the periodic orbit
%       PlotEigen (false) - Plot the eigenvalues of the orbit on the
%       complex plane
%       dS (.01) - Pseudo Arclength step
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
%    Options.Nm = 50;
%    Options.Np = 50;
%    Option.dS = .001;
%    SymmetricOrbitFamily(System,InitialGuess,Options)
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function Orbit = SymmetricOrbitFamily(System,InitialGuess,Options)

% Options:
DefaultOptions.Nm = 10; % Backward Steps
DefaultOptions.Np = 10; % Forward Steps
DefaultOptions.Plot = true; % Plotting
DefaultOptions.Colorbar = true; % Jacobi Constant colorbar
DefaultOptions.PlotEigen = false;
DefaultOptions.dS = .01; % Pseudostep
DefaultOptions.N = 1; % Shooting Arcs
DefaultOptions.Save = false; 
DefaultOptions.SaveDirectory = fileparts(mfilename('fullpath'));
DefaultOptions.Regularize = false;
DefaultOptions.Colormap = cool;
DefaultShooting.AbsoluteTolerance = 10^-15;
DefaultShooting.RelativeTolerance = 10^-12;
DefaultShooting.Tolerance = 10^-10;
DefaultShooting.MaximumIterations = 200;
DefaultShooting.Printing = false;
DefaultShooting.RootFinder = 'GaussNewton';
DefaultOptions.Shooting = DefaultShooting;
if nargin < 3
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.Planar = size(InitialGuess.x,1) == 4;

% Construct Orbit Options:
OrbitOptions.N = Options.N;
OrbitOptions.Regularize = Options.Regularize;
OrbitOptions.Shooting = Options.Shooting;
OrbitOptions.Planar = Options.Planar;
OrbitOptions.Plot = false;

% Solve for Starting Periodic Orbit if necessary:
if isfield(InitialGuess,'JacobiConstant') % Initial Guess was given as a periodic orbit
    Orbit(Options.Nm+1) = InitialGuess;
else
    Orbit(Options.Nm+1) = SymmetricOrbit(System,InitialGuess,OrbitOptions);
end

% Set to Free time:
OrbitOptions.Shooting.FreeTime = true;

% Loop Though Negative Orbits:
for i = Options.Nm:-1:1
    % Get Tangent Direction of Previous Orbit:
    [~,J] = ShootingRootFunction([Orbit(i+1).x(:,1);Orbit(i+1).Period/2],SymmetricOrbitTPBVP(System,struct('FreeTime',true,'N',1,'Planar',OrbitOptions.Planar)),OrbitOptions.Shooting);
    DX = null(J);
    if exist('DXprev')
        if norm(DX-DXprev) >= norm(-DX-DXprev)
            DX = -DX/norm(DX);
        else
            DX = DX/norm(DX);
        end
    end
    DXprev = DX;

    % Get Pseudoarclength BCs:
    BC = PseudoArclengthBC(DX,-Options.dS,Orbit(i+1).x(:,1),Orbit(i+1).Period/2);
    OrbitOptions = MergeStructs(BC,OrbitOptions);

    % Construct Initial Guess for Next Orbit:
    InitialGuess.x = Orbit(i+1).x(:,1);
    InitialGuess.P = Orbit(i+1).Period;

    % Get Orbit:
    Orbit(i) = SymmetricOrbit(System,InitialGuess,OrbitOptions);
    if ~Orbit(i).Converged
        Orbit(1:i) = [];
        Options.Nm = Options.Nm-i;
        disp(['Pseudo Arclength Failed to Converge at step ' num2str(i) ', terminating early.']);
        break;
    end

    % Print:
    disp(['Pseudoarclength Converged at step ' num2str(i)]);

end
clear("DXprev");

% Loop Though Positive Orbits:
for i = Options.Nm+2:Options.Nm+Options.Np
    % Get Tangent Direction of Previous Orbit:
    [~,J] = ShootingRootFunction([Orbit(i-1).x(:,1);Orbit(i-1).Period/2],SymmetricOrbitTPBVP(System,struct('FreeTime',true,'N',1,'Planar',OrbitOptions.Planar)),OrbitOptions.Shooting);
    DX = null(J);
    if exist('DXprev')
        if norm(DX-DXprev) >= norm(-DX-DXprev)
            DX = -DX/norm(DX);
        else
            DX = DX/norm(DX);
        end
    end
    DXprev = DX;

    % Get Pseudoarclength BCs:
    BC = PseudoArclengthBC(DX,Options.dS,Orbit(i-1).x(:,1),Orbit(i-1).Period/2);
    OrbitOptions = MergeStructs(BC,OrbitOptions);

    % Construct Initial Guess for Next Orbit:
    InitialGuess.x = Orbit(i-1).x(:,1);
    InitialGuess.P = Orbit(i-1).Period;

    % Get Orbit:
    Orbit(i) = SymmetricOrbit(System,InitialGuess,OrbitOptions);
    if ~Orbit(i).Converged
        Options.Np = i-Options.Nm-2;
        disp(['Pseudo Arclength Failed to Converge at step ' num2str(i) ', terminating early.']);
        break;
    end

    % Print:
    disp(['Pseudo Arclength Converged at step ' num2str(i)]);

end

% Plotting
Options.HighlightIndices = Options.Nm+1;
if Options.Plot
    f = TabFigures(CR3BPFamilyFigure(System,Orbit,Options), "Rotating Frame Orbit");
end
if Options.PlotInertial
    f = TabFigures(InertialFamilyFigure(System,Orbit,Options), "Inertial Frame Orbit",f);
end

% Saving Data and Figure:
if Options.Save
    SaveData(struct('Orbit',Orbit),Options);
    SaveFigure(f,Options);
end

end

