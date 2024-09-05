
%% Reach Test for cislunar:

Options.Planar = true;
Options.N = 4;
Options.AutomaticDifferentiation = false;
System = CislunarSystem;
System = InitializeCR3BP(System);
f = @(t,x) EOMCR3BP(System,0,x,Options);
T = 0:1:5;
x0 = [.2765;0;0;2.1379];
Phi = StateTransitionTensors(f,T,x0,Options);

N = 1000;
delta = .01;
E = delta*[zeros(length(x0)-2,N);-ones(1,N/4) linspace(-1,1,N/4) ones(1,N/4) linspace(1,-1,N/4);linspace(1,-1,N/4) -ones(1,N/4) linspace(-1,1,N/4) ones(1,N/4)]+x0;
% E = delta*[zeros(length(x0)-2,N);cos(linspace(0,2*pi,N));sin(linspace(0,2*pi,N))]+x0;

%% Reach:
dR = PolynomialZonotope([[0;0;delta;0] [0;0;0;delta]],[1 0;0 1]);
% dR = CartesianProduct(PolynomialZonotope([0;0],zeros(0,1)),PolynomialZonotope(Ellipsoid(delta^2*eye(2),[0;0])));
R(1) = ExactSum(PolynomialZonotope(x0,zeros(0,1)),dR);
for i = 2:length(T)
    R(i) = PolynomialZonotope(Phi{1,i},zeros(0,1));
    R(i) = ExactSum(R(i),LinearMap(dR,Phi{2,i}));
    R(i) = ExactSum(R(i),TensorMap(dR,1/2*Phi{3,i}));
    R(i) = ExactSum(R(i),TensorMap(dR,1/6*Phi{4,i}));
end

% Propagate Points:
for i = 1:size(Phi,1)
    XSTT{i} = zeros(length(x0),length(T),length(E));
    for j = 1:length(E)
        XSTT{i}(:,:,j) = TensorPropagate(E(:,j),Phi(1:i,:));
    end
end
X = zeros(length(x0),length(T),length(E));
for i = 1:length(E)
    [t,x] = ode89(f,T,E(:,i),odeset('AbsTol',10^-15,'RelTol',10^-12));
    X(:,:,i) = x';
end

% Plot:
figure; hold on;
plot(Project(R(end),[1 2]));
plot(squeeze(XSTT{4}(1,end,:)),squeeze(XSTT{4}(2,end,:)),'r');
plot(squeeze(X(1,end,:)),squeeze(X(2,end,:)),'k');
legend('Polynomial Zonotope','STT Trajectories','Exact Trajectories');

%% With Inputs:

System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

R(1) = PolynomialZonotope(x0,zeros(0,1));

f = @(t,z) [EOMCR3BP(System,0,z(1:4),Options) + System.Spacecraft.Thrust/System.Spacecraft.WetMass*[0;0;z(5);z(6)];0;0];
T = 0:.1:.2;
x0 = [.2765;0;0;2.1379;0;0];
Options.Broken = true;
Phi = StateTransitionTensors(f,T,x0,Options);

dR = PolynomialZonotope([0;0;0;0],zeros(0,1));
for i = 2:length(T)
    U = PolynomialZonotope(Ellipsoid([1 0;0 1],[0;0]));
    X = CartesianProduct(dR,U);
    dR = Project(ExactSum(ExactSum(LinearMap(X,Phi{2,i}),TensorMap(X,1/2*Phi{3,i})),TensorMap(X,1/6*Phi{4,i})),[1:4]);
    R(i) = ExactSum(PolynomialZonotope(Phi{1,i}(1:4,:),zeros(0,1)),dR);
end
