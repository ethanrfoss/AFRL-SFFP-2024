
%% STT Test With Inputs:
System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);
Options.Planar = true;
f = @(t,x) EOMCR3BP(System,t,x(1:4),Options);
fu = @(t,x) [f(t,x) + System.Spacecraft.Thrust/System.Spacecraft.WetMass*[0 0;0 0;1 0;0 1]*x(5:6);0;0];

% Compute:
Options.N = 3;
T = 0:1;
x0 = [.2765;0;0;2.1379;0;0];
Phi = StateTransitionTensors(f,T,x0,Options);

% DSTTs:
[Phi,R] = DirectionalStateTransitionTensors(f,T,x0,Options);

