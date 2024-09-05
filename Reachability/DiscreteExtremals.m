
function Trajectory = DiscreteExtremals(EOM,T,x0,Options)

% Options:
DefaultOptions.N = 1; % Number of Trajectories
DefaultOptions.Iterations = 1;
DefaultOptions.Integrator = @ode45;
DefaultOptions.ODESettings = odeset('AbsTol',10^-15,'RelTol',10^-12);
DefaultOptions.InputSet = 'UnitSphere';
DefaultOptions.BVP = false;
DefaultOptions.RootFinder = 'GaussNewton';
DefaultOptions.CostateSampling = 'UnitSphere';
DefaultOptions.Tolerance = 10^-6;
if nargin < 4
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.nx = length(x0);

% Linearize:
dims = InputDimensions(EOM);
dims = struct('nx',dims(1),'nu',dims(2));
[A,B] = Linearize(EOM,dims);

for i = 1:Options.N

    % Preallocate:
    Trajectory(i).lambda = zeros(dims.nx,length(T));
    Trajectory(i).lambda(:,end) = RandomSample(Options.CostateSampling,Options);
    uprev = zeros(dims.nu,length(T));
    [xprev,PhiA,PhiB] = Discretize(EOM,A,B,T,uprev,x0,dims,Options);
    Trajectory(i).u = zeros(dims.nu,length(T)); 
    
    for j = 1:Options.Iterations
        
        % Discretize:
        if ~all(all(uprev==Trajectory(i).u))
            uprev = Trajectory(i).u;
            [xprev,PhiA,PhiB] = Discretize(EOM,A,B,T,uprev,x0,dims,Options);
        end

        % Back Propagate Costates:
        for k = length(T)-1:-1:1
            Trajectory(i).lambda(:,k) = PhiA(:,:,k)'*Trajectory(i).lambda(:,k+1);
            Trajectory(i).u(:,k) = -PhiB(:,:,k)'*Trajectory(i).lambda(:,k)/sqrt(Trajectory(i).lambda(:,k)'*PhiB(:,:,k)*PhiB(:,:,k)'*Trajectory(i).lambda(:,k));
        end

        % Propagate:
        Trajectory(i).x = Propagate(EOM,T,Trajectory(i).u,x0,dims,Options);

        if norm(xprev(:,end)-Trajectory(i).x(:,end)) <= Options.Tolerance
            break;
        end

    end
end

end

function [A,B] = Linearize(EOM,dims)

x = sym('x',[dims.nx,1]);
u = sym('u',[dims.nu,1]);
t = sym('t',1);
f = EOM(x,u,t);
A = matlabFunction(jacobian(f,x),"Vars",{x,u,t});
B = matlabFunction(jacobian(f,u),"Vars",{x,u,t});

end

function [x,PhiA,PhiB] = Discretize(f,A,B,T,u,x0,dims,Options)

x = zeros(dims.nx,length(T));
x(:,1) = RandomSample(x0);
PhiA = zeros(dims.nx,dims.nx,length(T)-1);
PhiB = zeros(dims.nx,dims.nu,length(T)-1);
for i = 1:length(T)-1
    z0 = [x(:,i);reshape(eye(dims.nx),[dims.nx^2,1]);zeros(dims.nx*dims.nu,1)];
    [~,z] = Options.Integrator(@(t,z) [f(z(1:dims.nx),u(:,i),t);reshape(A(z(1:dims.nx),u(:,i),t)*reshape(z(dims.nx+1:dims.nx+dims.nx^2),[dims.nx,dims.nx]),[dims.nx^2,1]);reshape(A(z(1:dims.nx),u(:,i),t)*reshape(z(dims.nx+dims.nx^2+1:dims.nx+dims.nx^2+dims.nx*dims.nu),[dims.nx,dims.nu])+B(z(1:dims.nx),u(:,i),t),[dims.nx*dims.nu,1])],[T(i) T(i+1)],z0,Options.ODESettings);
    x(:,i+1) = z(end,1:dims.nx)';
    PhiA(:,:,i) = reshape(z(end,dims.nx+1:dims.nx+dims.nx^2)',[dims.nx,dims.nx]);
    PhiB(:,:,i) = reshape(z(end,dims.nx+dims.nx^2+1:dims.nx+dims.nx^2+dims.nx*dims.nu)',[dims.nx,dims.nu]); 
end

end

function [x] = Propagate(f,T,u,x0,dims,Options)

x = zeros(dims.nx,length(T));
x(:,1) = RandomSample(x0,Options);
for i = 1:length(T)-1
    [~,z] = Options.Integrator(@(t,x) f(x,u(:,i),t),[T(i) T(i+1)],x(:,i),Options.ODESettings);
    x(:,i+1) = z(end,1:dims.nx)';
end

end

function x = RandomSample(X,Options)

if isnumeric(X)
    if size(X,2) == 1
        x = X;
    else
        x = X*normalize(randn(size(X,2),1),'norm',2);
    end
elseif isa(X,'function_handle')
    x = X();
elseif isa(X,'char')
    if isequal(X,'UnitSphere')
        x = randn(Options.nx,1);
        x = x/norm(x);
    else
        error('Invalid Random Sampling');
    end
else
    error('Invalid Random Sampling');
end

end