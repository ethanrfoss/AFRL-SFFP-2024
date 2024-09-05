
function STTNRHOComparison

% Setup System
System = CislunarSystem;
System.Spacecraft = LunarIceCube;
System = InitializeCR3BP(System);

% Set up EOM:
Options.Planar = false;
Options.TensorOrder = 4;
% Options.Ephemeris = true;
f = EOMCR3BP(System,Options);

% Initial Condition and Displacements:
x0 = [1.0133975094440237363357937283581;
      -0.0000000000000000000000000004572425006398181505229021818393;
                          -0.17537156531341002718704658036586;
                                                            0;
                         -0.083683586737119008969187916591181;
                                                            0];
delta = .005;
T = 1.396*[0:4];

% Initial set:
N = 500;
t1 = (0:N-1)*2*pi/((1+sqrt(5))/2);
t2 = acos(1-2*(0:N-1 + .5)/N);
X0 = x0 + delta*[cos(t1).*sin(t2);sin(t1).*sin(t2);cos(t2);zeros(3,N)];
TT = convhull(X0(1,:),X0(2,:),X0(3,:));

% Propagate State Transition Tensors:
Options.N = Options.TensorOrder;
Options.ODESettings = odeset('AbsTol',10^-9,'RelTol',10^-6);
Options.Integrator = @ode113;
Phi = StateTransitionTensors2(f,T,x0,Options);

% Propagate Each Point at Each Order:
XSTT{1} = zeros(6,length(T),N);
XSTT{2} = zeros(6,length(T),N);
XSTT{3} = zeros(6,length(T),N);
for i = 2:Options.N
    for j = 1:N
        XSTT{i-1}(:,:,j) = TensorPropagate(X0(:,j),Phi(1:i,:));
    end
end

% Propagate ODEs:
X = zeros(6,length(T),N);
Options.TensorOrder = 1;
f = EOMCR3BP(System,Options);
for j = 1:N
    [t,x] = ode89(@(t,x) f(x,t),T,X0(:,j),Options.ODESettings);
    X(:,:,j) = x';
end

% Plot:
colors = hsv(4);
for i = 1:length(T)
    figure; hold on;
    trisurf(TT,squeeze(X(1,i,:)),squeeze(X(2,i,:)),squeeze(X(3,i,:)),'EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(1,:));
    trisurf(TT,squeeze(XSTT{1}(1,i,:)),squeeze(XSTT{1}(2,i,:)),squeeze(XSTT{1}(3,i,:)),'EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(2,:));
    trisurf(TT,squeeze(XSTT{1}(1,i,:)),squeeze(XSTT{2}(2,i,:)),squeeze(XSTT{2}(3,i,:)),'EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(3,:));
    trisurf(TT,squeeze(XSTT{1}(1,i,:)),squeeze(XSTT{3}(2,i,:)),squeeze(XSTT{3}(3,i,:)),'EdgeColor','none','FaceAlpha',.5,'FaceColor',colors(4,:));
    legend('Exact','First Order','Second Order','Third Order');
end

end