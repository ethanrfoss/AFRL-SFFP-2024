
% STT Indirect Test:

clear all; close all;

% Set Up Van Der Pol:
mu = 0; 
ep = 10^-6;
syms x1 x2;
x = [x1;x2];
f = [x2;mu*(1-x1^2)*x2-x1];
A = jacobian(f,x);
B = [0;.03];
syms lambda1 lambda2;
lambda = [lambda1;lambda2];
z = [x;lambda];
g = [f-B*B'*lambda/(sqrt(lambda'*B*B'*lambda+10^-6));-A'*lambda];
Ag = jacobian(g,z);
g = matlabFunction(g,"Vars",{z});
Ag = matlabFunction(Ag,"Vars",{z});
f = matlabFunction(f,"Vars",{x});
clear x lambda z;

% Propagate Exact Reachable Set:
N = 1000;
T = 5;
dt = .1;
t = 0:dt:T;
x0 = [1;1];

lambda0 = [1*cos(linspace(0,2*pi,N));1*sin(linspace(0,2*pi,N))];

x = zeros(2,length(t),N);
lambda = zeros(2,length(t),N);
u = zeros(1,length(t),N);
for i = 1:N
    z0 = [x0;lambda0(:,i)];
    [~,z] = ode89(@(t,z) g(z(1:4)),t,z0);
    x(:,:,i) = z(:,1:2)';
    lambda(:,:,i) = z(:,3:4)';
    u(:,:,i) = B'*z(:,3:4)'/norm(B'*z(:,3:4)');
end

% Costates Plot:
% figure; hold on
% for i = 1:N
%     plot(lambda(1,:,i),lambda(2,:,i),'Color',[0 0 0 .1])
% end

% Determine State Transition Tensors about a reference:
Options.N = 4;
Phi = StateTransitionTensors(@(t,z) g(z),t,[x0;1;0],Options);
for i = 1:N
    z0 = [x0;lambda0(:,i)];
    z = TensorPropagate(z0,Phi);
    xstt(:,:,i) = z(1:2,:);
    lambdastt(:,:,i) = z(3:4,:);
    ustt(:,:,i) = B'*z(3:4,:)/norm(B'*z(3:4,:));
end

% Plot Exact and STT Reach sets:
figure; hold on;
for i = 1:length(t)
    plot(squeeze(x(1,i,:)),squeeze(x(2,i,:)),'k');
    plot(squeeze(xstt(1,i,:)),squeeze(xstt(2,i,:)),'b');
end
