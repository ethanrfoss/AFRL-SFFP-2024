%% CR3BPFigure - Creates a Figure of the CR3BP
%
% Inputs: 
%    System - Structure Containing the System variables, which is
%    constructed using System = InitializeCR3BP(System).
%    Options - Structure Containing Plotting Options. The Options and their
%    default values are:
%       Labels (true) - Adds axis labels
%       LagrangePoints (true) - Plots Lagrange Points
%       LagrangePointsLabels (true) - Labels the Lagrange Points
%       Ticks (true) - Adds ticks to the axes
%       TickSpacing (.5) - Tick Spacing
%       Planar (false) - Planar or Non-Planar CR3BP
%       ZeroVelocityCurves (false) - Plots Zero Velocity Curves
%       JacobiConstant (3) - Jacobi Constant Used for Zero Velocity Curves
%       Axis ([-1.4 1.4 -1.4 1.4 -.2 .2]) - Axis Values of the plot
%       PrimaryScale (1) - Size Scale of the Primary Body
%       SecondaryScale (1) - Size Scale of the Secondary Body
%       PrimaryColor ([0,.4,.8]) - Color of the Primary Body
%       SecondaryColor ([.5,.5,.5]) - Color of the Secondary Body
%       ZeroVelocityCurvesColor ([0.9294    0.6941    0.1255]) - Color of
%       the Zero Velocity Curves
%       Box (false) - Box On/Off
%       Grid (false) - Grid On/Off
%
% Outputs:
%    f - Figure handle
%
% Example Usage:
%    System = CislunarSystem;
%    System = InitializeCR3BP(System);
%    FigOptions.Planar = true;
%    FigOptions.ZeroVelocityCurves = true;
%    FigOptions.JacobiConstant = 3.1722;
%    CR3BPFigure(System,FigOptions);
%    FigOptions.Planar = false;
%    FigOptions.Axis = [.8 1.2 -.2 .2 -.15 .15];
%    FigOptions.Box = true;
%    FigOptions.LagrangePointLabels = true;
%    CR3BPFigure(System,FigOptions);
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function f = InertialFigure(System,Options)

% Options:
DefaultOptions.Labels = true;
DefaultOptions.Ticks = true;
DefaultOptions.TickSpacing = .5;
DefaultOptions.Planar = false;
DefaultOptions.Axis = [-1.4 1.4 -1.4 1.4 -.2 .2];
DefaultOptions.PrimaryScale = 1;
DefaultOptions.SecondaryScale = 1;
DefaultOptions.PrimaryColor = [0,.4,.8];
DefaultOptions.SecondaryColor = [.5,.5,.5];
DefaultOptions.Box = false;
DefaultOptions.Grid = false;
% DefaultOptions.NewFigure = true;
if nargin < 2
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

% Create Figure:
% if Options.NewFigure
%     figure; 
% end
f = figure; hold on; axis equal;

% Set to Latex:
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% Grid and Box
grid(Options.Grid);
box(Options.Box);

% Set Ticks:
ax = gca;
ax.Units = 'normalized';
if Options.Ticks
    ax.XTick = -5:Options.TickSpacing:5;
    ax.YTick = -5:Options.TickSpacing:5;
    if ~Options.Planar
        ax.ZTick = -5:Options.TickSpacing:5;
    end
else
    ax.XTick = [];
    ax.YTick = [];
    if ~Options.Planar
        ax.ZTick = [];
    end
end

% Axis Labels:
if Options.Labels
    xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex'); 
    if ~Options.Planar
        zlabel('$z$','interpreter','latex');
    end
end

% Plot Bodies:
if Options.Planar
    X1 = cos(linspace(0,2*pi-2*pi/50,50))*Options.PrimaryScale*System.R1;
    Y1 = sin(linspace(0,2*pi-2*pi/50,50))*Options.PrimaryScale*System.R1;
    P1 = polyshape(X1,Y1);
    plot(P1,"EdgeColor",Options.PrimaryColor,"FaceColor",Options.PrimaryColor,"FaceAlpha",1);
else
    [Xsphere,Ysphere,Zsphere] = sphere(50);
    surf(-System.mu+Xsphere*Options.PrimaryScale*System.R1,Ysphere*Options.PrimaryScale*System.R1,Zsphere*Options.PrimaryScale*System.R1,'FaceColor',Options.PrimaryColor,'EdgeColor',Options.PrimaryColor);
end
X2 = cos(linspace(0,2*pi,50));
Y2 = sin(linspace(0,2*pi,50));
plot(X2,Y2,"Color",Options.SecondaryColor,"LineWidth",2,'LineStyle','--');

% Axis:
if Options.Planar
    axis(Options.Axis(1:4));
else
    axis(Options.Axis);
end

end