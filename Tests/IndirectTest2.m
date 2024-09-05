
System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

B = [zeros(2,2);System.Spacecraft.Thrust/System.Spacecraft.WetMass*eye(2,2)];

Options.Planar = true;
Options.TensorOrder = 1;
Options.Spacecraft = false;
Options.Ephemeris = false;
f = EOMCR3BP(System,Options);
A = JacobianCR3BP(System,Options);

lambdadot = @(z) -A(z(1:4),0).'*z(5:8);
g = @(z,t) [f(z(1:4),t)-B*B'*z(5:8)/max(norm(B'*z(5:8)),10^-8);lambdadot(z)/norm(z(5:8))-z(5:8)*z(5:8)'*lambdadot(z)/norm(z(5:8))^3];
% g = @(z,t) [f(z(1:4),t)-B*B'*z(5:8)/max(norm(B'*z(5:8)),10^-8);-A(z(1:4),t)'*z(5:8)];

N = 1000;
T = 2.5;
load('MMR1to1Family.mat');
x0 = Orbit(50).x(:,1);

% lambda0 = 2*[cos(linspace(0,2*pi,N));sin(linspace(0,2*pi,N))];

%% ODE:

xf = zeros(4,N);
lambdaf = zeros(4,N);
for i = 1:N
    lambda0 = [rand(4,1)];
    lambda0 = lambda0/norm(lambda0);
    z0 = [x0;lambda0];
    [t,z] = ode89(@(t,z)g(z,t),[0,T],z0);
    xf(:,i) = z(end,1:4)';
    x{i} = z(:,1:4)';
    lambdaf(:,i) = z(end,5:8)';
end
% for i = N/2+1:N
%     lambda0 = [0;0;rand(2,1)];
%     lambda0 = lambda0/norm(lambda0);
%     z0 = [x0;lambda0];
%     [t,z] = ode89(@(t,z)g(z),[0,T],z0);
%     xf(:,i) = z(end,1:4)';
%     x{i} = z(:,1:4)';
%     lambdaf(:,i) = z(end,5:8)';
% end

figure; hold on;
% color = jet(N);
% subplot(1,2,1); hold on;
% for i = 1:length(xf)-1
%     plot(lambda0(1,i:i+1),lambda0(2,i:i+1),'Color',color(i,:),'LineWidth',3);
% end
% subplot(1,2,2); hold on;
% for i = 1:length(xf)-1
%     plot(x{i}(1,:),x{i}(2,:),'Color',[0 0 0 .2]);
%     % plot(xf(1,i:i+1),xf(2,i:i+1),'Color',color(i,:),'LineWidth',3);
% end
plot(xf(1,:),xf(2,:),'ro');


%% TPBVP:
System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

B = [zeros(2,2);System.Spacecraft.Thrust/System.Spacecraft.WetMass*eye(2,2)];

Options.Planar = true;
Options.TensorOrder = 1;
Options.Spacecraft = false;
f = EOMCR3BP(System,Options);
A = JacobianCR3BP(System,Options);

lambda0 = ones(4,1)/2;
TPBVP.Variables = @(S) deal(sym('z',[8,1]),sym([]));
TPBVP.Dynamics = @(S,z) [f(z(1:4),0)-B*B.'*z(5:8)/(sqrt(z(5:8).'*B*B.'*z(5:8))+10^-6);-A(z(1:4),0).'*z(5:8)];
TPBVP.BoundaryConditions = @(S,z0,zf,p,t) [z0(1:4)-x0;zf(5:8)-lambda0];
TPBVP = ComputeTPBVP(System,TPBVP);

N = 200;
lambda0 = [cos(linspace(0,2*pi,N));sin(linspace(0,2*pi,N));zeros(1,N);zeros(1,N)];
Options.RootFinder = 'GaussNewton';
Options.MaximumIterations = 1000;

i = 1;
while true
    % lambda0 = randn(4,1); lambda0 = lambda0/norm(lambda0);
    TPBVP.g = @(z0,zf,p,t) [z0(1:4)-x0;zf(5:8)-lambda0(:,i)];
    InitialGuess.x = [x0;randn(4,1)];
    InitialGuess.T = T;
    InitialGuess.p = [];

    Solution = Shooting(System,TPBVP,InitialGuess,Options);

    if ~Solution.Converged
        continue;
    end

    Trajectory(i).t = Solution.t;
    Trajectory(i).x = Solution.x(1:4,:);
    Trajectory(i).lambda = Solution.x(5:8,:);
    Trajectory(i).u = zeros(2,length(Trajectory(i).t));
    for j = 1:length(Trajectory(i).t)
        Trajectory(i).u(:,j) = -B'*Trajectory(i).lambda(:,j)/(sqrt(Trajectory(i).lambda(:,j)'*B*B'*Trajectory(i).lambda(:,j))+10^-6);
    end
    break;
end

i = 2;
Options.RootFinder = 'fsolve';
while i < N
    % lambda0 = randn(4,1); lambda0 = lambda0/norm(lambda0);
    TPBVP.g = @(z0,zf,p,t) [z0(1:4)-x0;zf(5:8)-lambda0(:,i)];
    InitialGuess.x = [x0;Trajectory(i-1).lambda(:,1)];
    InitialGuess.T = T;
    InitialGuess.p = [];

    Solution = Shooting(System,TPBVP,InitialGuess,Options);

    if ~Solution.Converged
        continue;
    end

    Trajectory(i).t = Solution.t;
    Trajectory(i).x = Solution.x(1:4,:);
    Trajectory(i).lambda = Solution.x(5:8,:);
    Trajectory(i).u = zeros(2,length(Trajectory(i).t));
    for j = 1:length(Trajectory(i).t)
        Trajectory(i).u(:,j) = -B'*Trajectory(i).lambda(:,j)/(sqrt(Trajectory(i).lambda(:,j)'*B*B'*Trajectory(i).lambda(:,j))+10^-6);
    end
    i = i + 1
end

Options.Axis = [-.25 2.25 -1 1.5];
Options.LagrangePointLabels = true;
CR3BPFigure(System,Options);
hold on;
for i = 1:length(Trajectory)
    plot(Trajectory(i).x(1,:),Trajectory(i).x(2,:),'Color',[0 0 0 .2]);
    plot(Trajectory(i).x(1,end),Trajectory(i).x(2,end),'g.');
end

%% FFOE:

System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

Options.Planar = true;
Options.TensorOrder = 1;
Options.Spacecraft = true;
f = EOMCR3BP(System,Options);

T = 0:.01:2.5;
Options.Iterations = 1;
Options.N = 100;
Trajectory = DiscreteExtremals(f,T,x0,Options);

Options.Axis = [-.25 2.25 -1 1.5];
Options.LagrangePointLabels = true;
CR3BPFigure(System,Options);
hold on;
for i = 1:length(Trajectory)
    plot(Trajectory(i).x(1,:),Trajectory(i).x(2,:),'Color',[0 0 0 .2]);
    plot(Trajectory(i).x(1,end),Trajectory(i).x(2,end),'b.');
end