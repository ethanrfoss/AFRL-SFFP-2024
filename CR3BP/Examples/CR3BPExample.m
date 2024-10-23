
%% CR3BP Example:

%% Comparison of CR3BP and Ephemeris Model:

% Setup Cislunar System:
System = CislunarSystem; % Get parameters of cislunar system
System = InitializeCR3BP(System); % Initialize System (compute Lagrange points, mass parameter, etc)

% Get Equations of Motion for CR3BP:
Options.Planar = false;
Options.Ephemeris = false;

f1 = EOMCR3BP(System,Options); % Function Handle

% Get Equations of Motion for Ephemeris:
Options.Planar = false;
Options.Ephemeris = true;
Options.Epoch = '2024-09-06 12:00:00 TDB';
Options.Primary = {'Earth','Moon'};
Options.Perturber = {'Sun','Venus'};

f2 = EOMCR3BP(System,Options); % Function Handle

% Get an Initial Condition an Integrate:
x0 = [.8;0;0;0;.2;0];
[t1,x1] = ode113(@(t,x) f1(x,t),[0 10],x0); x1 = x1';
[t2,x2] = ode113(@(t,x) f2(x,t),[0 10],x0); x2 = x2';

% Plot Results:
Options.Labels = true;
Options.LagrangePoints = true;
Options.LagrangePointLabels = true;
Options.ZeroVelocityCurves = true;
Options.JacobiConstant = JacobiConstant(System,x0);
CR3BPFigure(System,Options);
p1 = plot3(x1(1,:),x1(2,:),x1(3,:),'LineWidth',2);
p2 = plot3(x2(1,:),x2(2,:),x2(3,:),'LineWidth',2);
legend([p1,p2],'CR3BP','Ephemeris');
view(-15,45);

%% Periodic Orbit Generation in the CR3BP:

clear Options;

% Set Options:
Options.Plot = true;
Options.Save = false; 
Options.FreeTime = false;

% % Create Initial Guess:
% InitialGuess.x = [1.08;0;.202;0;-.201;0];
% InitialGuess.P = 2.5;
% 
% % Get Periodic Orbit:
% Orbit = SymmetricOrbit(System,InitialGuess,Options);

% Example 2:
InitialGuess.x = [.9985;0;-.01488;0;1.144387;-.00121];
InitialGuess.P = 6.537;

% Get Periodic Orbit:
Orbit = SymmetricOrbit(System,InitialGuess,Options);

%% Generate a Family About Orbit:

% Options:
Options.Nm = 50; % Backward Steps
Options.Np = 150; % Forward Steps
Options.dS = .01; % Pseudostep

% Generate:
SymmetricOrbitFamily(System,Orbit,Options);
xlim([.8 1.2]); ylim([-.25 .25]);
view(30,17);

%% Generate a Quasi-Periodic Orbit:

% Create Initial Guess:
InitialGuess.x = [1.08;0;.202;0;-.201;0];
InitialGuess.P = 2.5;

% Get Periodic Orbit:
Orbit = SymmetricOrbit(System,InitialGuess,Options);
view(40,17);

% Set Options:
Options.Plot = true; % Plotting
Options.PlotInertial = false;
Options.N1 = 25; % Invariant Circle Points
Options.N0 = 1; % Shooting Arcs
Options.R = .1; % Initial Guess Radius

% Generate with GMOS:
QPO = GMOS(System,Orbit,Options);
view(40,17);

 