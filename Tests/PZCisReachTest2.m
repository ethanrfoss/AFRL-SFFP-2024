
System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

Options.Planar = true;
f = @(x,u) EOMCR3BP(System,0,x,Options) + System.Spacecraft.Thrust/System.Spacecraft.WetMass*[zeros(2,2);eye(2,2)]*u;
x0 = [.2765;0;0;2.1379];

sys = nonlinearSys(f);

params.tFinal = 5; %final time
params.R0 = polyZonotope(x0,[],[],[]);
params.U = polyZonotope([0;0],[0 1 0 -4.9348 0 4.0587 0;0 0 3.1416 0 -5.1677 0 2.5502],[],[0 1 1 1 1 1 1;0 0 1 2 3 4 5]);


options.timeStep = 0.1;
options.zonotopeOrder = 50; %zonotope order
options.taylorTerms = 4;
options.intermediateOrder = 50;
options.errorOrder = 20;
options.lagrangeRem.simplify = 'simplify';

% reachability algorithm
options.alg = 'poly';
options.tensorOrder = 3;

% special settings for polynomial zonotopes
polyZono.maxDepGenOrder = 50;
polyZono.maxPolyZonoRatio = 0.01;
polyZono.restructureTechnique = 'reducePca';

R = reach(sys, params, options);