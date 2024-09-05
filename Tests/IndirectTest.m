
clear all; close all;

A = [0 1;-1 0];
B = [0;.3];

g = @(z) [A*z(1:2)-B*B'*z(3:4)/max(norm(B'*z(3:4)+10^-6));-A'*z(3:4)];

N = 100;
T = 5;
x0 = [1;1];

lambda0 = [1*cos(linspace(0,2*pi,N));1*sin(linspace(0,2*pi,N))];

xf = zeros(2,N);
lambdaf = zeros(2,N);
for i = 1:N
    z0 = [x0;lambda0(:,i)];
    [t{i},z] = ode89(@(t,z)g(z),[0,T],z0);
    xf(:,i) = z(end,1:2)';
    x{i} = z(:,1:2)';
    u{i} = B'*z(:,3:4)';
    u{i} = u{i}./vecnorm(u{i},2,1);
    lambdaf(:,i) = z(end,3:4)';
end

figure;
color = jet(N);
subplot(2,2,1); hold on;
for i = 1:N-1
    plot(lambda0(1,i:i+1),lambda0(2,i:i+1),'Color',color(i,:),'LineWidth',3);
end
subplot(2,2,2); hold on;
for i = 1:N-1
    plot(x{i}(1,:),x{i}(2,:),'Color',[0 0 0 .2]);
    plot(xf(1,i:i+1),xf(2,i:i+1),'Color',color(i,:),'LineWidth',3);
end
subplot(2,2,[3,4]); hold on;
for i = 1:N-1
    plot(t{i},u{i},'Color',color(i,:),'LineWidth',1);
end

%% Confirm Reachability with yalmip:
N = 100;
t = linspace(0,T,N);
dt = t(2)-t(1);

R = expm([A B;0 0 0]*dt);
Ad = R(1:2,1:2);
Bd = R(1:2,3);

AA = Ad^N;
BB = zeros(2,N);
for i = 1:N
    BB(:,i) = Ad^(N-i)*Bd;
end

u = sdpvar(N,1);
x = sdpvar(2,1);
figure; hold on;
plot([x == Ad^(N-1)*x0+BB*u, u>=-ones(N,1), u<=ones(N,1)],x);
for i = 1:N-1
    plot(xf(1,i:i+1),xf(2,i:i+1),'Color',color(i,:),'LineWidth',3);
end

%% Propagate a trajectory
% optimize([[3.25438;0.160948] == BB*u, u>=-ones(N,1), u<=ones(N,1)],[],sdpsettings('solver','mosek'));
% ut = value(u);
% xt = zeros(2,length(ut));
% xt(:,1) = x0;
% for i = 1:length(ut)-1
%     [t,xr] = ode89(@(t,x) A*x+B*ut(i),[0,dt],xt(:,i),odeset('AbsTol',10^-15,'RelTol',10^-12));
%     xt(:,i+1) = xr(end,:)';
% end
% 
% plot(xt(1,:),xt(2,:),'Color','k','LineWidth',3);