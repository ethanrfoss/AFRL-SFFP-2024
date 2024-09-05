
%% Zonotope Reachability:

Ac = [0 -.5;1 0];
Bc = [0;.1];
x0 = [1;1];
dt = .1;
T = 5;
t = 0:dt:T;

% Discretize Dynamics:
AB = expm(dt*[Ac Bc;zeros(1,3)]);
A = AB(1:2,1:2); B = AB(1:2,3);

%% Polynomial Zonotopes:
R(1) = polyZonotope(x0);
U = polyZonotope(0,0,1,[]);

for i = 1:length(t)-1
    R(i+1) = A*R(i) + B*U;
end

% Plot:
figure; hold on;
for i = 1:length(t)
    plot(R(i));
end

%% Constrained Polynomial Zonotope:
R(1) = polyZonotope(x0);
% U = conPolyZono(polyZonotope(0,0,1,[]));

for i = 1:length(t)-1
    R(i+1) = A*R(i) + B*U;
end

%% Compute Exact Reachable set through method of characteristics:
g = @(x,lambda) [Ac*x-Bc*Bc'*lambda/(sqrt(lambda'*Bc*Bc'*lambda+10^-6));-Ac'*lambda];
N = 1000;
lambda0 = [1*cos(linspace(0,2*pi,N));1*sin(linspace(0,2*pi,N))];

xf = zeros(2,N);
lambdaf = zeros(2,N);
for i = 1:N
    z0 = [x0;lambda0(:,i)];
    [~,z] = ode89(@(t,z) g(z(1:2),z(3:4)),t,z0);
    x{i} = z(:,1:2)';
    u{i} = B'*z(:,3:4)';
    plot(x{i}(1,:),x{i}(2,:),'Color',[0 0 0 .01]);
end

