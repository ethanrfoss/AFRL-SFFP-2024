
function [QPO] = QPOFamily(System,Orbit,Options)

% Options:
DefaultOptions.Plot = true; % Plotting
DefaultOptions.PlotInertial = true;
DefaultOptions.Save = false;
DefaultOptions.SaveDirectory = fileparts(mfilename('fullpath'));
DefaultOptions.Regularize = false;
DefaultOptions.N1 = 25; % Invariant Circle Points
DefaultOptions.N0 = 1; % Shooting Arcs
DefaultOptions.N = 50;
DefaultOptions.dS = .01;
DefaultOptions.R = .01; % Initial Guess Radius
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

% QPO Options:
QPOOptions = Options;
QPOOptions.Plot = false;
QPOOptions.PlotInertial = false;
QPOOptions.Save = false;

% Get Initial QPO:
if isfield(Orbit,'u') % Converged QPO Given
    QPO(1) = Orbit;
    Options.Planar = size(Orbit.u,1) == 4;
    Options.nx = size(Orbit.u,1);
else
    Options.Planar = size(Orbit.x,1) == 4;
    Options.nx = size(Orbit.x,1);
    QPO(1) = GMOS(System,Orbit,QPOOptions);
end
if QPO(1).Converged
    disp(['QPO Converged at iteration ' num2str(1)]);
else
    error(['QPO failed to converged at iteration ' num2str(1)]);
end

% Get Second QPO:
QPOOptions.QPO1 = QPO(1);
QPO(2) = GMOS(System,Orbit,QPOOptions);
if QPO(2).Converged
    disp(['QPO Converged at iteration ' num2str(2)]);
else
    error(['QPO failed to converged at iteration ' num2str(2)]);
end

% Loop Through:
for i = 3:Options.N
    
    % Set Reference:
    QPOOptions.QPO1 = QPO(i-1);
    QPOOptions.QPO2 = QPO(i-2);

    % Compute QPO:
    QPO(i) = GMOS(System,QPO(i-1),QPOOptions);

    % Print
    if QPO(i).Converged
        disp(['QPO Converged at iteration ' num2str(i)]);
    else
        error(['QPO failed to converged at iteration ' num2str(i)]);
    end
end

% Plotting
if Options.Plot
    f = TabFigures(CR3BPQPOFamilyFigure(System,QPO,Options), "Rotating Frame QPO");
end
% if Options.PlotInertial
%     f = TabFigures(InertialQPOFamilyFigure(System,QPO,Options), "Inertial Frame QPO",f);
% end

% Saving Data and Figure:
if Options.Save
    SaveData(struct('QPO',QPO),Options);
    SaveFigure(f,Options);
end

end