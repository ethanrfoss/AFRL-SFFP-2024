function [xr,x,u] = FastFirstOrderEstimateReachability(FFOE,System,Variables,Dynamics)

% % Get Prolem Dimensions:
[x,u] = Variables(System);
dims.nx = length(x);
dims.nu = length(u);
% 
% % Linearize:
[LinearizedDynamics] = Linearize(System,Variables,Dynamics);
% LinearizedDynamics.f

% Preallocate
xr = zeros(dims.nx,FFOE.Seg);
x = zeros(dims.nx,FFOE.Seg,FFOE.N);
u = zeros(dims.nu,FFOE.Seg,FFOE.N);
Fx = zeros(dims.nx,dims.nx,FFOE.Seg-1);
Fu = zeros(dims.nx,dims.nu,FFOE.Seg-1);

% Initial Conditions
xr(:,1) = System.x0;
t = System.t0;

% Simulated Open Loop Trajectory with no thrust to get Partials
% textprogressbar('Propogating Open Loop Trajectory: ');
for i = 1:FFOE.Seg
    [xr(:,i+1),Fx(:,:,i),Fu(:,:,i)] = StepFirstOrder(FFOE,xr(:,i),zeros(dims.nu,1),t,t+FFOE.dt,LinearizedDynamics,dims);
    t = t + FFOE.dt;
    % textprogressbar(min(100*i/FFOE.N,100));
end
% textprogressbar('Done');

% Use First Order Estimate to Backpropogate Thrust Control law
% textprogressbar('Back Propogating Costates: ');
for j = 1:FFOE.N
    p = randn(dims.nx,1).*FFOE.ReachStates;
    p = p/norm(p); % Initial Costate with Random Unit Vector
    for i = FFOE.Seg:-1:1
        u(1:dims.nu-1,i,j) = -1*Fu(1:dims.nx,1:dims.nu-1,i)'*p/norm(Fu(1:dims.nx,1:dims.nu-1,i)'*p,2);
        u(dims.nu,i,j) = 1;
        p = Fx(:,:,i)'*p;
    end
    % textprogressbar(min(100*j/FFOE.N,100));
end
% textprogressbar('Done');

% Recompute Trajectories Subject to Control Law
% textprogressbar('Simulating Trajectories: ');
for j = 1:FFOE.N
    t = System.t0;
    x(:,1,j) = System.x0;
    for i = 1:FFOE.Seg
        x(:,i+1,j) = StepFlow(FFOE,x(:,i,j),u(:,i,j),t,t+FFOE.dt,LinearizedDynamics);
        t = t + FFOE.dt;
    end
    % textprogressbar(min(100*j/FFOE.N,100));
end
% textprogressbar('Done\n');

end

function [x1] = StepFlow(FFOE,x0,u0,t0,t1,LinearizedDynamics)

tspan = [t0 t1];
[~,x] = ode89(@(t,x)xdot(x,u0,LinearizedDynamics),tspan,x0,odeset('AbsTol',FFOE.IntegrationAbsTol,'RelTol',FFOE.IntegrationRelTol)); % Perform Integration

x1 = x(end,:)';

end

function [x1,Fx1,Fu1] = StepFirstOrder(FFOE,x0,u0,t0,t1,LinearizedDynamics,dims)

tspan = [t0 t1];
xfxfu0 = [x0;reshape(eye(dims.nx,dims.nx),[dims.nx^2,1]);reshape(zeros(dims.nx,dims.nu),[dims.nx*dims.nu,1])];
[~,xfxfu] = ode89(@(t,xfxfu)[xdot(xfxfu(1:dims.nx),u0,LinearizedDynamics); Fxdot(xfxfu(1:dims.nx),u0,xfxfu(dims.nx+1:dims.nx+dims.nx^2),LinearizedDynamics,dims); Fudot(xfxfu(1:dims.nx),u0,xfxfu(dims.nx+dims.nx^2+1:dims.nx+dims.nx^2+dims.nx*dims.nu),LinearizedDynamics,dims)],tspan,xfxfu0,odeset('AbsTol',FFOE.IntegrationAbsTol,'RelTol',FFOE.IntegrationRelTol)); % Perform Integration

xfxfuf = xfxfu(end,:)';

x1 = xfxfuf(1:dims.nx);
Fx1 = reshape(xfxfuf(dims.nx+1:dims.nx+dims.nx^2),[dims.nx,dims.nx]);
Fu1 = reshape(xfxfuf(dims.nx+dims.nx^2+1:dims.nx+dims.nx^2+dims.nx*dims.nu),[dims.nx,dims.nu]);

end

function xdot = xdot(x,u,LinearizedDynamics)

xdot = LinearizedDynamics.f(x,u);

end

function Fxdot = Fxdot(x,u,Fx,LinearizedDynamics,dims)

A = LinearizedDynamics.A(x,u);

Fx = reshape(Fx,[dims.nx,dims.nx]);

Fxdot = A*Fx;

Fxdot = reshape(Fxdot,[dims.nx^2,1]);

end

function Fudot = Fudot(x,u,Fu,LinearizedDynamics,dims)

A = LinearizedDynamics.A(x,u);
B = LinearizedDynamics.B(x,u);

Fu = reshape(Fu,[dims.nx,dims.nu]);

Fudot = A*Fu + B;

Fudot = reshape(Fudot,[dims.nx*dims.nu,1]);

end