function QPO = QPOInitialGuess(System,Orbit,Options)

% Sanity Checks
if mod(Options.N0,2) ~= 1 || mod(Options.N1,2) ~= 1
    error('An odd Number of QPO points must be given');
end

% Simulate for One Period
if Options.N0 == 1
    T = Orbit.Period;
    x = Orbit.x(:,1);
    Phi = eye(Options.nx);
else
    T = ones(1,Options.N0)*Orbit.Period/Options.N0;
    xphi0 = [Orbit.x(:,1); reshape(eye(Options.nx),[Options.nx^2,1])];
    [t,xphi] = ode89(@(t,xphi)[EOMCR3BP(System,t,xphi(1:Options.nx),Options); STMCR3BP(System,t,xphi(1:Options.nx),xphi(Options.nx+1:Options.nx+Options.nx^2),Options)],[0 cumsum(T(1:end-1))],xphi0,odeset('AbsTol',Options.Shooting.AbsoluteTolerance,'RelTol',Options.Shooting.RelativeTolerance)); %, 'Events',@(t,x)PoincareMapIntersection(t,x,[0;-1;0],[0;0;0])));
    x = xphi(:,1:Options.nx)';
    Phi = reshape(xphi(:,Options.nx+1:Options.nx+Options.nx^2)',[Options.nx,Options.nx,length(t)]);
end

% theta0, theta1
theta0 = 2*pi*(0:Options.N0-1)/Options.N0;
theta1 = 2*pi*(0:Options.N1-1)/Options.N1;

% State grid u is N0 by N1 by 6 matrix
QPO.u = zeros(Options.nx,Options.N0,Options.N1);

% Determine Monodromy Matrix
[V,E] = eig(Orbit.Monodromy);

% Decompose Eigenstructure
Es = cplxpair(diag(E));
[~,ind] = ismember(diag(E),Es(1));
w10 = V(:,find(ind));
rho = angle(Es(1));

% Determine Invariant Circles along Each Time grid point
for i = 1:length(T)
    QPO.u(:,i,:) = reshape((Phi(:,:,i)*Options.R*(real(w10)*cos(theta1-rho*sum(T(1:i-1))/Orbit.Period)-imag(w10)*sin(theta1-rho*sum(T(1:i-1))/Orbit.Period)) + x(:,i)),[Options.nx,1,Options.N1]);
end

% Reshape u to vectoral form:
QPO.uv = reshape(QPO.u,Options.nx*Options.N0*Options.N1,1);

% Initial Rotation Rates:
QPO.w0 = 2*pi/Orbit.Period; 
QPO.w1 = rho/Orbit.Period;
QPO.rho = rho;
QPO.T = T';

end