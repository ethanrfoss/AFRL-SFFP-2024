
%% State Transition Tensors Test:

% 3BP:
% Options.N = 3;
% Options.AutomaticDifferentiation = false;
% System = CislunarSystem;
% System = InitializeCR3BP(System);
% System.mu = .0121505856;
% f = @(t,x) EOMCR3BP(System,0,x);
% T = [0:.3:.9]*1.396265;
% x0 = [1.013417655693384;0;-.175374764978708;0;-.083721347178432;0];
% [Phi,R] = DirectionalStateTransitionTensors(f,T,x0,Options);
% E = .001*[zeros(3,100);cos(linspace(0,2*pi,100));zeros(1,100);sin(linspace(0,2*pi,100))]+x0;
% EllipsoidComparison(f,x0,T,E,Phi,R);

% 3BP:
Options.Planar = true;
Options.N = 4;
Options.AutomaticDifferentiation = false;
System = CislunarSystem;
System = InitializeCR3BP(System);
f = @(t,x) EOMCR3BP(System,0,x,Options);
T = 0:.25:5;
x0 = [.2765;0;0;2.1379];
% Options.EigenCutoff = 10^-100;
[Phi,R] = DirectionalStateTransitionTensors(f,T,x0,Options);
E = .02*[zeros(length(x0)-2,100);cos(linspace(0,2*pi,100));sin(linspace(0,2*pi,100))]+x0;
CR3BPFigure(System,struct('LagrangePointLabels',true));
EllipsoidComparison(f,x0,T,E,Phi,R);

function EllipsoidComparison(f,x0,T,E,Phi,R)

X = zeros(length(x0),length(T),length(E));
XSTT = cell(size(Phi,1));
for i = 1:length(E)
    [t,x] = ode89(f,T,E(:,i),odeset('AbsTol',10^-15,'RelTol',10^-12));
    X(:,:,i) = x';
end
for i = 1:size(Phi,1)
    XSTT{i} = zeros(size(X));
    for j = 1:length(E)
        XSTT{i}(:,:,j) = DirectionalTensorPropagate(E(:,j),Phi(1:i,:),R);
    end
end

% figure; hold on;
Colors = ['b','g','r','y','m'];
for i = 1:length(T)
    p1 = plot(squeeze(X(1,i,:)),squeeze(X(2,i,:)),'k');
    for j = 2:size(Phi,1)
        p2(j-1) = plot(squeeze(XSTT{j}(1,i,:)),squeeze(XSTT{j}(2,i,:)),Colors(j));
    end
end
Legends = {'Truth','STM','2nd Order STT','3rd Order STT','4th Order STT'};
legend([p1,[p2(:)]'],Legends{1:size(Phi,1)});

end



