
%% State Transition Tensors Test:

% System 1: Linear
% f = @(t,x) [0 x(2);x(1) 0];
% Options.N = 2;
% T = [0 1 2 3];
% Phi = StateTransitionTensors(f,T,[1;1],Options);

% Van der Pol:
mu = 1; 
f = @(t,x) [x(2);mu*(1-x(1)^2)*x(2)-x(1)];
Options.N = 5;
T = 0:5;
x0 = [0;.5];
Phi = StateTransitionTensors(f,T,x0,Options);
E = .2*[cos(linspace(0,2*pi,1000));sin(linspace(0,2*pi,1000))]+x0;
% E = .1*[-linspace(-1,1,500) ones(1,500) linspace(1,-1,500) -ones(1,500);-ones(1,500) linspace(-1,1,500) ones(1,500) linspace(1,-1,500)]+x0;
EllipsoidComparison(f,x0,T,E,Phi);

% Propagate the Set as a Zonotope:
% figure; hold on;
% Gi = diag([.1,.1]);
R0 = conPolyZono(ellipsoid([.2 0;0 .2],x0));
for i = 1:size(Phi,2)
    Q{1} = squeeze(1/2*Phi{3,i}(1,:,:));
    Q{2} = squeeze(1/2*Phi{3,i}(2,:,:));
    R(i) = mtimes(Phi{2,i},R0) + quadMap(R0,R0,Q);
    % plot(R(i));
end


% Duffing: 
% f = @(t,x) [x(2);-x(1)^3-x(1)];
% Options.N = 4;
% T = 0:.7:4;
% x0 = [2;2];
% Phi = StateTransitionTensors(f,T,x0,Options);
% EllipsoidComparison(f,x0,T,100,.2,Phi);

% 3BP:
% Options.Planar = true;
% Options.N = 4;
% Options.AutomaticDifferentiation = false;
% System = CislunarSystem;
% System = InitializeCR3BP(System);
% f = @(t,x) EOMCR3BP(System,0,x,Options);
% T = 0:.25:5;
% x0 = [.2765;0;0;2.1379];
% Phi = StateTransitionTensors(f,T,x0,Options);
% E = .02*[zeros(length(x0)-2,1000);cos(linspace(0,2*pi,1000));sin(linspace(0,2*pi,1000))]+x0;
% 
% CR3BPFigure(System);
% EllipsoidComparison(f,x0,T,E,Phi);

function EllipsoidComparison(f,x0,T,E,Phi)

X = zeros(length(x0),length(T),length(E));
XSTT = cell(size(Phi,1));
for i = 1:length(E)
    [t,x] = ode89(f,T,E(:,i),odeset('AbsTol',10^-15,'RelTol',10^-12));
    X(:,:,i) = x';
end
for i = 1:size(Phi,1)
    XSTT{i} = zeros(size(X));
    for j = 1:length(E)
        XSTT{i}(:,:,j) = TensorPropagate(E(:,j),Phi(1:i,:));
    end
end

figure; hold on;
Colors = ['b','g','r','y','m'];
for i = 1:length(T)
    p1 = plot(squeeze(X(1,i,:)),squeeze(X(2,i,:)),'k');
    for j = 2:3%size(Phi,1)
        p2(j-1) = plot(squeeze(XSTT{j}(1,i,:)),squeeze(XSTT{j}(2,i,:)),Colors(j));
    end
end
Legends = {'Truth','STM','2nd Order STT','3rd Order STT','4th Order STT'};
legend([p1,[p2(:)]'],Legends{1:3});

end



