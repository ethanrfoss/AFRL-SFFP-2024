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
function f = CR3BPFigure(System,Options)

% Options:
DefaultOptions.Labels = true;
DefaultOptions.LagrangePoints = true;
DefaultOptions.LagrangePointLabels = false;
DefaultOptions.Ticks = true;
DefaultOptions.TickSpacing = .5;
DefaultOptions.Planar = false;
DefaultOptions.ZeroVelocityCurves = false;
DefaultOptions.JacobiConstant = 3;
DefaultOptions.Axis = [-1.4 1.4 -1.4 1.4 -.2 .2];
DefaultOptions.PrimaryScale = 1;
DefaultOptions.SecondaryScale = 1;
DefaultOptions.PrimaryColor = [0,.4,.8];
DefaultOptions.SecondaryColor = [.5,.5,.5];
DefaultOptions.ZeroVelocityCurvesColor = [0.9294    0.6941    0.1255];
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
    xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); 
    if ~Options.Planar
        zlabel('$\zeta$','interpreter','latex');
    end
end

% Zero Velocity Curves
if Options.ZeroVelocityCurves
    Phi = @(X,Y,Z) 1/2*(X.^2+Y.^2) + (1-System.mu)./sqrt((X+System.mu).^2+Y.^2 + Z.^2) + System.mu./sqrt((X-1+System.mu).^2+Y.^2 + Z.^2);
    if Options.Planar
        fimplicit(@(x,y) 2*Phi(x,y,0)-Options.JacobiConstant,Options.Axis(1:4),'Color','k');
        X = linspace(Options.Axis(1), Options.Axis(2), 500); Y = linspace(Options.Axis(3), Options.Axis(4), 500); [X, Y] = meshgrid(X, Y);
        contourf(X, Y, -2*Phi(X,Y,zeros(size(X)))+Options.JacobiConstant, [0 0], 'LineColor', 'none'); % Filling the regions
        colormap([1 1 1; Options.ZeroVelocityCurvesColor]); % White for regions not filled
        caxis([-0.001 0.001]); % Set the color axis to highlight the zero level
    else
        fimplicit3(@(x,y,z) 2*Phi(x,y,z) - Options.JacobiConstant,Options.Axis,'FaceAlpha',.5,'FaceColor',Options.ZeroVelocityCurvesColor,'EdgeColor','none','MeshDensity',100);
        X = linspace(Options.Axis(1),Options.Axis(2),100);
        Y = linspace(Options.Axis(3),Options.Axis(4),100);
        Z = linspace(Options.Axis(5),Options.Axis(6),100);
        [X,Y,Z] = meshgrid(X,Y,Z);
        V = 2*Phi(X,Y,Z) - Options.JacobiConstant;
        t1 = contourslice(X,Y,Z,V,Options.Axis(1),[],[],[0 0]); for i = 1:length(t1) t1(i).EdgeColor = 'k'; end
        t2 = contourslice(X,Y,Z,V,Options.Axis(2),[],[],[0 0]); for i = 1:length(t2) t2(i).EdgeColor = 'k'; end
        t3 = contourslice(X,Y,Z,V,[],Options.Axis(3),[],[0 0]); for i = 1:length(t3) t3(i).EdgeColor = 'k'; end
        t4 = contourslice(X,Y,Z,V,[],Options.Axis(4),[],[0 0]); for i = 1:length(t4) t4(i).EdgeColor = 'k'; end
        t5 = contourslice(X,Y,Z,V,[],[],Options.Axis(5),[0 0]); for i = 1:length(t5) t5(i).EdgeColor = 'k'; end
        t6 = contourslice(X,Y,Z,V,[],[],Options.Axis(6),[0 0]); for i = 1:length(t6) t6(i).EdgeColor = 'k'; end
    end
end

% Plot Bodies:
if Options.Planar
    X1 = cos(linspace(0,2*pi-2*pi/50,50))*Options.PrimaryScale*System.R1-System.mu;
    Y1 = sin(linspace(0,2*pi-2*pi/50,50))*Options.PrimaryScale*System.R1;
    X2 = cos(linspace(0,2*pi-2*pi/50,50))*Options.SecondaryScale*System.R2+1-System.mu;
    Y2 = sin(linspace(0,2*pi-2*pi/50,50))*Options.SecondaryScale*System.R2;
    P1 = polyshape(X1,Y1);
    P2 = polyshape(X2,Y2);
    plot(P1,"EdgeColor",Options.PrimaryColor,"FaceColor",Options.PrimaryColor,"FaceAlpha",1);
    plot(P2,"EdgeColor",Options.SecondaryColor,"FaceColor",Options.SecondaryColor,"FaceAlpha",1);
else
    [Xsphere,Ysphere,Zsphere] = sphere(50);
    surf(-System.mu+Xsphere*Options.PrimaryScale*System.R1,Ysphere*Options.PrimaryScale*System.R1,Zsphere*Options.PrimaryScale*System.R1,'FaceColor',Options.PrimaryColor,'EdgeColor',Options.PrimaryColor);
    surf(1-System.mu+Xsphere*Options.SecondaryScale*System.R2,Ysphere*Options.SecondaryScale*System.R2,Zsphere*Options.SecondaryScale*System.R2,'FaceColor',Options.SecondaryColor,'EdgeColor',Options.SecondaryColor);
end

% Lagrange Points:
if Options.LagrangePoints
    plot(System.LagrangePoints.L1(1),System.LagrangePoints.L1(2),'r*');
    plot(System.LagrangePoints.L2(1),System.LagrangePoints.L2(2),'r*');
    plot(System.LagrangePoints.L3(1),System.LagrangePoints.L3(2),'r*');
    plot(System.LagrangePoints.L4(1),System.LagrangePoints.L4(2),'g*');
    plot(System.LagrangePoints.L5(1),System.LagrangePoints.L5(2),'g*');
    if Options.LagrangePointLabels
        if (System.LagrangePoints.L1(1) >= Options.Axis(1) && System.LagrangePoints.L1(1) <= Options.Axis(2)) text(System.LagrangePoints.L1(1),System.LagrangePoints.L1(2),'$L_1$','VerticalAlignment', 'top', 'HorizontalAlignment', 'left','Color','r'); end
        if (System.LagrangePoints.L2(1) >= Options.Axis(1) && System.LagrangePoints.L2(1) <= Options.Axis(2)) text(System.LagrangePoints.L2(1),System.LagrangePoints.L2(2),'$L_2$','VerticalAlignment', 'top', 'HorizontalAlignment', 'left','Color','r'); end
        if (System.LagrangePoints.L3(1) >= Options.Axis(1) && System.LagrangePoints.L3(1) <= Options.Axis(2)) text(System.LagrangePoints.L3(1),System.LagrangePoints.L3(2),'$L_3$','VerticalAlignment', 'top', 'HorizontalAlignment', 'left','Color','r'); end
        if (System.LagrangePoints.L4(1) >= Options.Axis(1) && System.LagrangePoints.L4(1) <= Options.Axis(2) && System.LagrangePoints.L4(2) >= Options.Axis(3) && System.LagrangePoints.L4(2) <= Options.Axis(4)) text(System.LagrangePoints.L4(1),System.LagrangePoints.L4(2),'$L_4$','VerticalAlignment', 'top', 'HorizontalAlignment', 'left','Color','g'); end
        if (System.LagrangePoints.L5(1) >= Options.Axis(1) && System.LagrangePoints.L5(1) <= Options.Axis(2) && System.LagrangePoints.L5(2) >= Options.Axis(3) && System.LagrangePoints.L5(2) <= Options.Axis(4)) text(System.LagrangePoints.L5(1),System.LagrangePoints.L5(2),'$L_5$','VerticalAlignment', 'top', 'HorizontalAlignment', 'left','Color','g'); end
    end
end

% Axis:
if Options.Planar
    axis(Options.Axis(1:4));
else
    axis(Options.Axis);
end

end