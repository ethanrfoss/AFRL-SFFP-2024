
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
x0 = [.2765;0;0;2.1379];
u0 = [0;0];
T = 0:.1:3;

%% ODE Extremals:
Trajectory = Extremals(f1,T,x0,struct('InputSet','UnitSphere','N',100,'BVP',false));

Options.Axis = [-.25 2.25 -1 1.5];
Options.LagrangePointLabels = true;
CR3BPFigure(System,Options);
hold on;
for i = 1:length(Trajectory)
    plot(Trajectory(i).x(1,:),Trajectory(i).x(2,:),'Color',[0 0 0 .2]);
    plot(Trajectory(i).x(1,end),Trajectory(i).x(2,end),'r.');
end

%% FFOE:
% x = sym('x',[4,1]);
% u = sym('u',[2,1]);
% Variables = @(S) deal(x,u);
% FFOE.Seg = 30;
% FFOE.dt = .1;
% FFOE.N = 100;
% Dynamics = @(S,x,u) f1(x,u,0);
% System.x0 = x0;
% System.t0 = 0;
% FFOE.IntegrationAbsTol = 10^-9;
% FFOE.IntegrationRelTol = 10^-6;
% FFOE.ReachStates = [1;1;0;0];
% [Xref,X,~] = FastFirstOrderEstimateReachability(FFOE,System,Variables,Dynamics);
Trajectory3 = DiscreteExtremals(f1,T,x0,struct('InputSet','UnitSphere','N',100,'CostateSampling',[eye(2,2);zeros(2,2)]));

Options.Axis = [-.25 2.25 -1 1.5];
Options.LagrangePointLabels = true;
CR3BPFigure(System,Options);
hold on;
for i = 1:length(Trajectory3)
    plot(Trajectory3(i).x(1,:),Trajectory3(i).x(2,:),'Color',[0 0 0 .2]);
    plot(Trajectory3(i).x(1,end),Trajectory3(i).x(2,end),'b.');
end

%% Indirect Extremals:
Trajectory2 = Extremals(f1,T,x0,struct('InputSet','UnitSphere','N',100,'BVP',true,'CostateSampling',[eye(2,2);zeros(2,2)]));

Options.Axis = [-.25 2.25 -1 1.5];
Options.LagrangePointLabels = true;
CR3BPFigure(System,Options);
hold on;
for i = 1:length(Trajectory2)
    plot(Trajectory2(i).x(1,:),Trajectory2(i).x(2,:),'Color',[0 0 0 .2]);
    plot(Trajectory2(i).x(1,end),Trajectory2(i).x(2,end),'g.');
end

%% Compute Reachable Sets:
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

%% Compute Reachable Sets (CORA):
% Options.STTOrder = 3;
% Options.MaxPZOrder = 100;
% Options.ReductionOrder = 100;
% % R = Reachability(f,T,x0,u0,x0,PolynomialZonotope(Ellipsoid([1 0;0 1],[0;0]),struct('Order',2,'Method','Squircle','ID','Unique')),Options);
% U = PolynomialZonotope(Ellipsoid([1 0;0 1],[0;0]));
% R = Reachability(f,T,x0,u0,x0,polyZonotope([0;0],U.G,[],U.E),Options);
% 
% CR3BPFigure(System,Options);
% hold on;
% color = autumn(length(T));
% for i = 2:length(T)
%     plot(project(R(i),[1,2]));
%     i
% end

%% Combine:
figure; hold on;
hold on;
for i = 1:length(Trajectory)
    p1 = plot(Trajectory(i).x(1,end),Trajectory(i).x(2,end),'r.');
end
for i = 1:length(Trajectory3)
    p2 = plot(Trajectory3(i).x(1,end),Trajectory3(i).x(2,end),'b.');
end
for i = 1:length(Trajectory2)
    p3 = plot(Trajectory2(i).x(1,end),Trajectory2(i).x(2,end),'g.');
end
p4 = plot(Project(R(end),[1,2]));
xlabel('$\xi$');
ylabel('$\eta$');
legend([p1,p2,p3,p4],'ODE','$\texttt{FFOE}$','Indirect','Algorithm 1');
grid on;
axis equal;

figure; hold on;
hold on;
for i = 1:length(Trajectory)
    p1 = plot(Trajectory(i).x(3,end),Trajectory(i).x(4,end),'r.');
end
for i = 1:size(X,3)
    p2 = plot(X(3,end,i),X(4,end,i),'b.');
end
for i = 1:length(Trajectory2)
    p3 = plot(Trajectory2(i).x(3,end),Trajectory2(i).x(4,end),'g.');
end
p4 = plot(Project(R(end),[3,4]));
xlabel('$\dot{\xi}$');
ylabel('$\dot{\eta}$');
legend([p1,p2,p3,p4],'ODE','$\texttt{FFOE}$','Indirect','Algorithm 1');
grid on;
axis equal;

%% Ephemeris:

System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

Options.Ephemeris = true;
Options.TensorOrder = 2;
Options.Spacecraft = false;
f = EOMCR3BP(System,Options);

EOM.f = @(x,t) [eye(6,6) zeros(6,36)]*f(x,t);
EOM.g = @(x,t) System.Spacecraft.Thrust/System.Spacecraft.WetMass*[zeros(3,3);eye(3,3)];
EOM.A = @(x,t) reshape([zeros(36,6) eye(36,36)]*f(x,t),[6,6]);

Options.TensorOrder = 2;
Options.Spacecraft = true;
f = EOMCR3BP(System,Options);

x0 = [.9133975094440237363357937283581;
      -0.01;
                          -0.16;
                                                            .03;
                         -0.08;
                                                            .01];

% x0 = [0.9136   -0.0477   -0.0845   -0.0065   -0.0421    0.3460]';
u0 = [0;0;0];
T = 0:.05:1;

Trajectory = Extremals(EOM,[0 1],x0,struct('InputSet','UnitCube','N',100,'BVP',false,'Integrator',@ode113,'ODESettings',odeset));

CR3BPFigure(System,Options);
hold on;
for i = 1:length(Trajectory)
    plot3(Trajectory(i).x(1,:),Trajectory(i).x(2,:),Trajectory(i).x(3,:),'Color',[0 0 0 .2]);
    plot3(Trajectory(i).x(1,end),Trajectory(i).x(2,end),Trajectory(i).x(3,end),'g.');
end

Options.STTOrder = 2;
Options.MaxPZOrder = 50;
Options.ReductionOrder = 50;
ReachSet = Reachability(f,T,x0,u0,x0,PolynomialZonotope(eye(3,3),eye(3,3)),Options);
% ReachSet = Reachability(f,0:.05:1,x0,u0,x0,PolynomialZonotope(Ellipsoid(eye(3,3),[0;0;0])),Options);

Options.Planar = true;
Options.Axis = [.8 1.1 -.15 .15];
Options.TickSpacing = .05;
CR3BPFigure(System,Options); hold on;
for i = 2:length(T)
    plot(Project(ReachSet(i),[1,2]));
    i
end
for i = 1:length(Trajectory)
    plot(Trajectory(i).x(1,:),Trajectory(i).x(2,:),'Color',[0 0 0 .1]);
end

% Options.Planar = false;
% CR3BPFigure(System,Options); hold on;
figure; hold on;
X2 = cos(linspace(0,2*pi-2*pi/50,50))*System.R2;
Y2 = sin(linspace(0,2*pi-2*pi/50,50))*System.R2;
P2 = polyshape(X2,Y2);
plot(P2,"EdgeColor","none","FaceColor",[.5 .5 .5],"FaceAlpha",1);
for i = 2:length(T)
    plot(Project(ReachSet(i),[2,3]));
    i
end
for i = 1:length(Trajectory)
    plot(Trajectory(i).x(2,:),Trajectory(i).x(3,:),'Color',[0 0 0 .1]);
end
axis equal;
xlabel('$\eta$');
ylabel('$\zeta$');
