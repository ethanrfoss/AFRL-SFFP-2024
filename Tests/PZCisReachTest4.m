
%% Cis Reach

% PZ Cis Reach Test 3:
System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

Options.Planar = true;
Options.Spacecraft = true;
Options.TensorOrder = 1;
Options.Ephemeris = false;
f1 = EOMCR3BP(System,Options);
Options.TensorOrder = 3;
% f = @(x,u,t) EOMCR3BP(System,0,x,Options) + System.Spacecraft.Thrust/System.Spacecraft.WetMass*[0;0;u(1);u(2)];
f = EOMCR3BP(System,Options);
load('MMR1to1Family.mat'); % Initial Orbit
x0 = Orbit(300).x(:,1);
u0 = [0;0];
T = 0:.1:3;

Options.STTOrder = 3;
Options.MaxPZOrder = 100;
Options.ReductionOrder = 100;
% R = Reachability(f,T,x0,u0,x0,PolynomialZonotope(Ellipsoid([1 0;0 1],[0;0]),struct('Order',2,'Method','Squircle','ID','Unique')),Options);
R = Reachability(f,T,x0,u0,x0,PolynomialZonotope([1 0;0 1],[1 0;0 1]),Options);

CR3BPFigure(System,Options);
hold on;
color = autumn(length(T));
for i = 2:length(T)
    plot(Project(R(i),[1,2]));
    i
end