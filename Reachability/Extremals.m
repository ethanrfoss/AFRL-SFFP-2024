
function Trajectory = Extremals(EOM,T,x0,Options)

% Options:
DefaultOptions.N = 1; % Number of Trajectories
DefaultOptions.Integrator = @ode45;
DefaultOptions.ODESettings = odeset('AbsTol',10^-15,'RelTol',10^-12);
DefaultOptions.InputSet = 'UnitSphere';
DefaultOptions.BVP = false;
DefaultOptions.RootFinder = 'GaussNewton';
DefaultOptions.CostateSampling = 'UnitSphere';
if nargin < 4
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.nx = length(x0);

% Decompose Dynamics into drift and input drift:
if isa(EOM,'struct')
    f = EOM.f;
    g = EOM.g;
    A = EOM.A;
    dims = InputDimensions(f);
    dims = struct('nx',dims(1),'nu',size(g(x0,0),2));

    u = ExtremalInput(g,dims,Options);
    F = @(z,t) [f(z(1:dims.nx),t) + g(z(1:dims.nx))*u(z(1:dims.nx),z(dims.nx+1:2*dims.nx),t);-A(z(1:dims.nx),t)'*z(dims.nx+1:2*dims.nx)];
else
    dims = InputDimensions(EOM);
    if dims(3) ~= 1
        error('Wrong Input Dimensions or Order for EOM');
    end
    dims = struct('nx',dims(1),'nu',dims(2));
    [f,g] = ControlAffineSystem(EOM,dims);

    % Compute Dynamics: 
    u = ExtremalInput(g,dims,Options);
    F = AugmentedDynamics(f,g,u,dims);
end

% Compute Trajectories:
if Options.BVP
    Trajectory = ExtremalBVP(F,u,T,x0,dims,Options);
else
    Trajectory = ExtremalIVP(F,u,T,x0,dims,Options);
end

end

function [f,g] = ControlAffineSystem(f,dims)

x = sym('x',[dims.nx,1]);
u = sym('u',[dims.nu,1]);
t = sym('t');
f = f(x,u,t);

try
    g = matlabFunction(jacobian(f,u),"Vars",{x,t});
catch
    error('System Must Be Control Affine');
end
f = matlabFunction(subs(f,u,zeros(dims.nu,1)),"Vars",{x,t});

end

function Trajectory = ExtremalIVP(F,u,T,x0,dims,Options)

Trajectory = repmat(struct('x',[],'lambda',[],'u',[],'t',[]),Options.N,1);
for i = 1:Options.N
    
    z0 = [RandomSample(x0);RandomSample(Options.CostateSampling,Options)];
    [t,z] = Options.Integrator(@(t,z) F(z,t),T,z0,Options.ODESettings);

    Trajectory(i).t = t;
    Trajectory(i).x = z(:,1:dims.nx)';
    Trajectory(i).lambda = z(:,dims.nx+1:2*dims.nx)';
    Trajectory(i).u = zeros(dims.nu,length(t));
    for j = 1:length(t)
        Trajectory(i).u(:,j) = u(Trajectory(i).x(:,j),Trajectory(i).lambda(:,j),t(j));
    end
    % i

end

end

function Trajectory = ExtremalBVP(F,u,T,X0,dims,Options)

% Compute TPBVP:
TPBVP.Variables = @(S) deal(sym('z',[2*Options.nx,1]),sym([]));
TPBVP.Dynamics = @(S,z) F(z,0);
TPBVP.BoundaryConditions = @(S,z0,zf,p,t) [z0(1:Options.nx)-RandomSample(X0);zf(Options.nx+1:2*Options.nx)];
TPBVP = ComputeTPBVP(struct(),TPBVP);

% Solve Series of TPBVPs:
Options.FreeTime = false;
Trajectory = repmat(struct('x',[],'lambda',[],'u',[],'t',[]),Options.N,1);
Options.RootFinder = 'fsolve';
Options.MaximumIterations = 100;
i = 1;
maxRuns = 10;
runs = 1;
while i < Options.N
    
    x0 = RandomSample(X0);
    if size(Options.CostateSampling,2) == Options.N
        TPBVP.g = @(z0,zf,p,t) [z0(1:Options.nx)-x0;zf(Options.nx+1:2*Options.nx)-Options.CostateSampling(:,i)];
        if i == 1 || ~Solution.Converged
            InitialGuess.x = [x0;RandomSample('UnitSphere',Options)];
        else
            InitialGuess.x = [x0;Trajectory(i-1).lambda(:,1)];
        end
    else
        TPBVP.g = @(z0,zf,p,t) [z0(1:Options.nx)-x0;zf(Options.nx+1:2*Options.nx)-RandomSample(Options.CostateSampling)];
        InitialGuess.x = [x0;RandomSample('UnitSphere',Options)];
    end
    InitialGuess.p = [];
    InitialGuess.T = T(end)-T(1);
    [Solution] = Shooting(struct(),TPBVP,InitialGuess,Options);
    % if ~Solution.Converged && i~=1
    %     % disp(['Failed to Converge at ' num2str(i) 'th iteration, breaking.']);
    %     % Trajectory(i:end) = [];
    %     % break;
    % 
    % end
    % if ~Solution.Converged && i==1
    %     continue;
    % end
    if ~Solution.Converged
        runs = runs + 1;
        if runs > maxRuns
            runs = 1;
            i = i + 1;
        end
        continue;
    end

    Trajectory(i).t = Solution.t;
    Trajectory(i).x = Solution.x(1:dims.nx,:);
    Trajectory(i).lambda = Solution.x(dims.nx+1:2*dims.nx,:);
    Trajectory(i).u = zeros(dims.nu,length(Trajectory(i).t));
    for j = 1:length(Trajectory(i).t)
        Trajectory(i).u(:,j) = u(Trajectory(i).x(:,j),Trajectory(i).lambda(:,j),Trajectory(i).t);
    end
    % i
    i = i + 1;
    runs = 1;

end
end

function u = ExtremalInput(g,dims,Options)

if isequal(Options.InputSet,'UnitSphere')    
    u = @(x,lambda,t) -g(x,t)'*lambda/(norm(g(x,t)'*lambda)+10^-6);
elseif isequal(Options.InputSet,'UnitCube')
    us = [zeros(dims.nu,1) (dec2bin(0:2^(dims.nu)-1)-'0')'*2-1];
    % us = [zeros(dims.nu,1) (dec2base(0:3^(dims.nu)-1,3)-'0')'-1];
    u = @(x,lambda,t) us(:,find(min(lambda'*g(x,t)*us) == lambda'*g(x,t)*us,1));
else
    error('Invalid Input Set Specified');
end

end

function F = AugmentedDynamics(f,g,u,dims)

% Compute Partials:
lambda = sym('lambda',[dims.nx,1]);
x = sym('x',[dims.nx,1]);
t = sym('t');
fc = matlabFunction(jacobian(lambda.'*f(x,t),x).',"Vars",{x,lambda,t});
gc = matlabFunction(jacobian(lambda.'*g(x,t),x).',"Vars",{x,lambda,t});

% Dynamics:
F = @(z,t) [f(z(1:dims.nx),t) + g(z(1:dims.nx),t)*u(z(1:dims.nx),z(dims.nx+1:2*dims.nx),t);-fc(z(1:dims.nx),z(dims.nx+1:2*dims.nx),t)-gc(z(1:dims.nx),z(dims.nx+1:2*dims.nx),t)*u(z(1:dims.nx),z(dims.nx+1:2*dims.nx),t)];

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
        x = normalize(randn(Options.nx,1),'norm',2);
    else
        error('Invalid Random Sampling');
    end
else
    error('Invalid Random Sampling');
end

end