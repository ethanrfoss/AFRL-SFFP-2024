
function QPO = GMOS(System,Orbit,Options)

% Options:
DefaultOptions.Plot = true; % Plotting
DefaultOptions.PlotInertial = true;
DefaultOptions.Save = false;
DefaultOptions.SaveDirectory = fileparts(mfilename('fullpath'));
DefaultOptions.Regularize = false;
DefaultOptions.N1 = 25; % Invariant Circle Points
DefaultOptions.N0 = 1; % Shooting Arcs
DefaultOptions.R = .01; % Initial Guess Radius
DefaultShooting.AbsoluteTolerance = 10^-15;
DefaultShooting.RelativeTolerance = 10^-12;
DefaultShooting.Tolerance = 10^-10;
DefaultShooting.MaximumIterations = 200;
DefaultShooting.Printing = false;
DefaultShooting.RootFinder = 'GaussNewton';
DefaultOptions.Shooting = DefaultShooting;
if nargin < 3
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

% Get Initial Guess of QPO:
if isfield(Orbit,'u') % Converged QPO Given
    QPO = Orbit;
    Options.Planar = size(Orbit.u,1) == 4;
    Options.nx = size(Orbit.u,1);
else
    Options.Planar = size(Orbit.x,1) == 4;
    Options.nx = size(Orbit.x,1);
    QPO = QPOInitialGuess(System,Orbit,Options);
end

% Compute QPO:
z0 = [reshape(permute(QPO.u,[1,3,2]),[Options.nx*Options.N0*Options.N1,1]);QPO.T;QPO.rho];
if isequal(Options.Shooting.RootFinder,'LevenbergMarquadt')
    [z0,Converged] = LevenbergMarquadt(@(z)GMOSRootFunction(z,System,Options),z0,Options);
elseif isequal(Options.Shooting.RootFinder,'GaussNewton')
    [z0,Converged] = GaussNewton(@(z)GMOSRootFunction(z,System,Options),z0,Options);
elseif isequal(Options.Shooting.RootFinder,'fsolve')
    [z0,~,flag] = fsolve(@(z)GMOSRootFunction(z,System,Options),z0,optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',Options.MaximumIterations));
    Converged = flag > 0;
else
    error('Invalid Root Finder Specified');
end

% Convergence:
QPO.Root = z0;
QPO.Converged = Converged;
if ~QPO.Converged
    QPO.u = [];
    QPO.T = [];
    QPO.rho = [];
    QPO.uv = [];
    QPO.w0 = [];
    QPO.w1 = [];
    QPO.Trajectory = [];
    QPO.EigenValues = [];
    QPO.StabilityIndex = [];
    QPO.EigenVectors = [];
    QPO.Du0 = [];
    QPO.Du1 = [];
    QPO.JacobiConstant = [];
    return;
end

% Decompose root:
QPO.u = permute(reshape(z0(1:Options.nx*Options.N0*Options.N1),[Options.nx,Options.N1,Options.N0]),[1,3,2]);
QPO.T = z0(Options.nx*Options.N0*Options.N1+1:Options.nx*Options.N0*Options.N1+Options.N0);
QPO.rho = z0(end);
QPO.uv = reshape(QPO.u,[Options.nx*Options.N0*Options.N1,1]);
QPO.w0 = 2*pi/sum(QPO.T);
QPO.w1 = QPO.rho/sum(QPO.T);

% Integrate Trajectories and STM for a full period:
D = kron(exp(1i*2*pi/Options.N1 * (-(Options.N1-1)/2:(Options.N1-1)/2)'*(0:(Options.N1-1))),eye(Options.nx,Options.nx));
Dinv = 1/Options.N1*kron(exp(-1i*2*pi/Options.N1*(0:(Options.N1-1))'*(-(Options.N1-1)/2:(Options.N1-1)/2)),eye(Options.nx,Options.nx));
Q = kron(diag(exp(-1i*-QPO.rho*(-(Options.N1-1)/2:(Options.N1-1)/2))),eye(Options.nx,Options.nx));
for i = 1:Options.N0
    Ju = zeros(Options.nx*Options.N1,Options.nx*Options.N1);
    for j = 1:Options.N1
        % Integrate Point:
        xphi0 = [QPO.u(:,i,j); reshape(eye(Options.nx),[Options.nx^2,1])];
        tspan = [0 sum(QPO.T)];
        [QPO.Trajectory{i}{j}.t,xphi] = ode89(@(t,xphi)[EOMCR3BP(System,t,xphi(1:Options.nx),Options); STMCR3BP(System,t,xphi(1:Options.nx),xphi(Options.nx+1:Options.nx+Options.nx^2),Options)],tspan,xphi0,odeset('AbsTol',Options.Shooting.AbsoluteTolerance,'RelTol',Options.Shooting.RelativeTolerance));
        QPO.Trajectory{i}{j}.x = xphi(:,1:Options.nx)';
        QPO.STM{i}{j} = reshape(xphi(end,Options.nx+1:Options.nx+Options.nx^2),[Options.nx,Options.nx]);
        Ju(Options.nx*(j-1)+1:Options.nx*j,Options.nx*(j-1)+1:Options.nx*j) = QPO.STM{i}{j};

        % Inertial State:
        QPO.Trajectory{i}{j}.InertialState = CR3BPFrameTransform(System,QPO.Trajectory{i}{j}.t,QPO.Trajectory{i}{j}.x,["Barycenter","Earth"],["Rotating","Inertial"]);
    end

    % Stability Properties:
    [V,E] = eig(real(Dinv*Q*D)*Ju);
    QPO.EigenValues{i} = diag(E);
    QPO.StabilityIndex = max(1/2*(abs(QPO.EigenValues{i})+1./abs(QPO.EigenValues{i})));
    % QPO.UnstableEigenVector{i} = reshape(V(:,find(max(abs(QPO.EigenValues{i}))==abs(QPO.EigenValues{i}),1)),[Options.nx,Options.N1]);
    % QPO.StableEigenVector{i} = V(:,find(min(abs(QPO.EigenValues{i}))==abs(QPO.EigenValues{i}),1));
    QPO.EigenVectors{i} = V;

end

% Post Processing:
% Jacobi Constant:
QPO.JacobiConstant = JacobiConstant(System,QPO.u(:,1,1));

% Partials:
[DFTMatrices] = DFT(Options.N1);
DFTMatrices.DH = kron(real(DFTMatrices.D*DFTMatrices.H),eye(Options.nx));
QPO.Du1 = zeros(size(QPO.u));
QPO.Du0 = zeros(size(QPO.u));
for i = 1:Options.N0
    QPO.Du1(:,i,:) = reshape(DFTMatrices.DH*reshape(QPO.u(:,i,:),[Options.nx*Options.N1,1]),[Options.nx,1,Options.N1]);
    for j = 1:Options.N1
        QPO.Du0(:,i,j) = 1/QPO.w0*(EOMCR3BP(System,0,QPO.u(:,i,j),Options)-QPO.w1*QPO.Du1(:,i,j));
    end
end

% Plotting
if Options.Plot
    f = TabFigures(CR3BPQPOFigure(System,QPO,Options), "Rotating Frame QPO");
end
if Options.PlotInertial
    f = TabFigures(InertialQPOFigure(System,QPO,Options), "Inertial Frame QPO",f);
end

% Saving Data and Figure:
if Options.Save
    SaveData(struct('QPO',QPO),Options);
    SaveFigure(f,Options);
end

end

function [F,J] = GMOSRootFunction(z,System,Options)

% Deconstruct Root:
u0 = permute(reshape(z(1:Options.nx*Options.N0*Options.N1),[Options.nx,Options.N1,Options.N0]),[1,3,2]);
T = z(Options.nx*Options.N0*Options.N1+1:Options.nx*Options.N0*Options.N1+Options.N0);
rho = z(end)*T/sum(T);

% Numerical Integration:
uf = zeros(Options.nx,Options.N0,Options.N1);
Ju = zeros(Options.nx*Options.N1,Options.nx*Options.N1,Options.N0);
fuf = zeros(Options.nx,Options.N0,Options.N1);
for i = 1:Options.N0
    for j = 1:Options.N1

        % Integrate Point:
        xphi0 = [u0(:,i,j); reshape(eye(Options.nx),[Options.nx^2,1])];
        tspan = [sum(T(1:i-1)) sum(T(1:i))];
        [~,xphi] = ode89(@(t,xphi)[EOMCR3BP(System,t,xphi(1:Options.nx),Options); STMCR3BP(System,t,xphi(1:Options.nx),xphi(Options.nx+1:Options.nx+Options.nx^2),Options)],tspan,xphi0,odeset('AbsTol',Options.Shooting.AbsoluteTolerance,'RelTol',Options.Shooting.RelativeTolerance));
        
        % Populate:
        uf(:,i,j) = xphi(end,1:Options.nx)';
        Ju(Options.nx*(j-1)+1:Options.nx*j,Options.nx*(j-1)+1:Options.nx*j,i) = reshape(xphi(end,Options.nx+1:Options.nx+Options.nx^2),[Options.nx,Options.nx]);
        fuf(:,i,j) = EOMCR3BP(System,0,uf(:,i,j),Options);

    end
end

% Debug Plot:
% figure; hold on;
% for i = 1:Options.N0
%     for j = 1:Options.N1
%         plot(u0(1,i,j),u0(2,i,j),'go');
%         plot(uf(1,i,j),uf(2,i,j),'ro');
%     end
% end

% Rotation Matrices:
D = kron(exp(1i*2*pi/Options.N1 * (-(Options.N1-1)/2:(Options.N1-1)/2)'*(0:(Options.N1-1))),eye(Options.nx,Options.nx));
Dinv = 1/Options.N1*kron(exp(-1i*2*pi/Options.N1*(0:(Options.N1-1))'*(-(Options.N1-1)/2:(Options.N1-1)/2)),eye(Options.nx,Options.nx));
Q = zeros(Options.N1*Options.nx,Options.N1*Options.nx,Options.N0);
for i = 1:Options.N0
    Q(:,:,i) = kron(diag(exp(-1i*-rho(i)*(-(Options.N1-1)/2:(Options.N1-1)/2))),eye(Options.nx,Options.nx));
end

% Construct Root Function and Jacobian:
F1 = zeros(Options.nx*Options.N0*Options.N1,1);
J11 = zeros(Options.nx*Options.N0*Options.N1,Options.nx*Options.N0*Options.N1);
if Options.N0 == 1
    F1 = real(Dinv*Q*D)*reshape(uf(:,1,:),[Options.nx*Options.N1,1]) - reshape(u0(:,1,:),[Options.nx*Options.N1,1]);
    J11 = real(Dinv*Q*D)*Ju-eye(Options.nx*Options.N1);
else
    for i = 1:Options.N0-1
        F1(Options.nx*Options.N1*(i-1)+1:Options.nx*Options.N1*i) = real(Dinv*Q(:,:,i)*D)*reshape(uf(:,i,:),[Options.nx*Options.N1,1]) - reshape(u0(:,i+1,:),[Options.nx*Options.N1,1]);
        J11(Options.nx*Options.N1*(i-1)+1:Options.nx*Options.N1*i,Options.nx*Options.N1*(i-1)+1:Options.nx*Options.N1*i) = real(Dinv*Q(:,:,i)*D)*Ju(:,:,i);
        J11(Options.nx*Options.N1*(i-1)+1:Options.nx*Options.N1*i,Options.nx*Options.N1*(i)+1:Options.nx*Options.N1*(i+1)) = -eye(Options.nx*Options.N1);
    end
    F1(Options.nx*Options.N1*(Options.N0-1)+1:Options.nx*Options.N1*Options.N0) = real(Dinv*Q(:,:,end)*D)*reshape(uf(:,end,:),[Options.nx*Options.N1,1]) - reshape(u0(:,1,:),[Options.nx*Options.N1,1]);
    J11(Options.nx*Options.N1*(Options.N0-1)+1:Options.nx*Options.N1*Options.N0,Options.nx*Options.N1*(Options.N0-1)+1:Options.nx*Options.N1*Options.N0) = real(Dinv*Q(:,:,Options.N0)*D)*Ju(:,:,Options.N0);
    J11(Options.nx*Options.N1*(Options.N0-1)+1:Options.nx*Options.N1*Options.N0,1:Options.nx*Options.N1) = -eye(Options.nx*Options.N1);
end
J12 = zeros(Options.nx*Options.N0*Options.N1,Options.N0);
for i = 1:Options.N0
    J12(Options.nx*Options.N1*(i-1)+1:Options.nx*Options.N1*i,i) = real(Dinv*Q(:,:,i)*D)*reshape(fuf(:,i,:),[Options.nx*Options.N1,1]);
end
J13 = zeros(Options.nx*Options.N0*Options.N1,1);
for i = 1:Options.N0
    J13(Options.nx*Options.N1*(i-1)+1:Options.nx*Options.N1*i,1) = T(i)/sum(T)*real(Dinv*Q(:,:,i)*kron(diag(1i*(-(Options.N1-1)/2:(Options.N1-1)/2)),eye(Options.nx,Options.nx))*D*reshape(uf(:,i,:),[Options.nx*Options.N1,1]));
end
F = F1;
J = [J11 J12 J13];

% Additional Root Functions:
if isfield(Options,'QPO1')
    F = [F;dot(z(1:Options.nx*Options.N0*Options.N1),reshape(permute(Options.QPO1.Du1,[1,3,2]),[Options.nx*Options.N0*Options.N1,1]))];
    J = [J; reshape(permute(Options.QPO1.Du1,[1,3,2]),[Options.nx*Options.N0*Options.N1,1])' zeros(1,Options.N0+1)];
    F = [F;dot(z(1:Options.nx*Options.N0*Options.N1)-reshape(permute(Options.QPO1.u,[1,3,2]),[Options.nx*Options.N0*Options.N1,1]),reshape(permute(Options.QPO1.Du0,[1,3,2]),[Options.nx*Options.N0*Options.N1,1]))];
    J = [J; reshape(permute(Options.QPO1.Du0,[1,3,2]),[Options.nx*Options.N0*Options.N1,1])' zeros(1,Options.N0+1)];
end
if isfield(Options,'JacobiConstant')
    Jacobi = 0;
    J21 = zeros(1,Options.nx*Options.N0*Options.N1);
    for i = 1:Options.N0
        for j = 1:Options.N1
            Jacobi = Jacobi + 1/(Options.N0*Options.N1)*JacobiConstant(System,u0(:,i,j));
            J21(1,Options.nx*(Options.N1*(i-1)+(j-1))+1:Options.nx*(Options.N1*(i-1)+j)) = 1/(Options.N0*Options.N1)*JacobiConstantPartial(System,u0(:,i,j));
        end
    end
    F = [F;Jacobi-Options.JacobiConstant];
    J = [J;J21 zeros(1,Options.N0+1)];
end
if isfield(Options,'Period')
    F = [F;sum(T)-Options.Period];
    J = [J;zeros(1,Options.nx*Options.N0*Options.N1) ones(1,Options.N0) 0];
end
if isfield(Options,'QPO2')
    DU = reshape(permute(Options.QPO1.u-Options.QPO2.u,[1,3,2]),[Options.nx*Options.N0*Options.N1,1]);
    F = [F;dot(reshape(permute(u0-Options.QPO1.u,[1,3,2]),[Options.nx*Options.N0*Options.N1,1]),DU/norm(DU))-Options.dS];
    J = [J;DU'/norm(DU) zeros(1,Options.N0+1)];
end

end