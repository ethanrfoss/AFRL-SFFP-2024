
%% State Transition Tensors Test:

% System 1: Linear
% f = @(t,x) [0 x(2);x(1) 0];
% Options.N = 2;
% T = [0 1 2 3];
% Phi = StateTransitionTensors(f,T,[1;1],Options);

% Van der Pol:
% mu = 1; 
% f = @(t,x) [x(2);mu*(1-x(1)^2)*x(2)-x(1)];
% Options.N = 4;
% T = 0:5;
% x0 = [2;2];
% Phi = StateTransitionTensors(f,T,x0,Options);
% EllipsoidComparison(f,x0,T,100,.2,Phi);

% Duffing: 
% f = @(t,x) [x(2);-x(1)^3-x(1)];
% Options.N = 4;
% T = 0:.7:4;
% x0 = [2;2];
% Phi = StateTransitionTensors(f,T,x0,Options);
% EllipsoidComparison(f,x0,T,100,.2,Phi);

% 3BP:
Options.Planar = true;
Options.N = 4;
Options.AutomaticDifferentiation = false;
System = CislunarSystem;
System = InitializeCR3BP(System);
f = @(t,x) EOMCR3BP(System,0,x,Options);
T = 0:.25:5;
x0 = [.2765;0;0;2.1379];
Phi = StateTransitionTensors(f,T,x0,Options);
E = .02*[zeros(length(x0)-2,1000);cos(linspace(0,2*pi,1000));sin(linspace(0,2*pi,1000))]+x0;

CR3BPFigure(System);
EllipsoidComparison(f,x0,T,E,Phi);

% HR3BP:
f = @(t,x) [x(3);x(4);2*x(4)-x(1)/(x(1)^2+x(2)^2)^(3/2);-2*x(3)-x(2)/(x(1)^2+x(2)^2)^(3/2)];
x0 = [.69010031015662;-.06716709529872;-.11045639526249;.03184084790390];
Phi = StateTransitionTensors(f,T,x0,Options);
E = [5.1*10^-4*cos(linspace(0,2*pi,1000));5.1*10^-4*sin(linspace(0,2*pi,1000));]+x0;
CR3BPFigure(System);
EllipsoidComparison(f,x0,T,E,Phi);


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

Colors = ['b','g','r','y'];
for i = 1:length(T)
    p1 = plot(squeeze(X(1,i,:)),squeeze(X(2,i,:)),'k');
    for j = 2:size(Phi,1)
        p2(j-1) = plot(squeeze(XSTT{j}(1,i,:)),squeeze(XSTT{j}(2,i,:)),Colors(j));
    end
end
Legends = {'Truth','STM','2nd Order STT','3rd Order STT','4th Order STT'};
legend([p1,[p2(:)]'],Legends{1:size(Phi,1)});

end



