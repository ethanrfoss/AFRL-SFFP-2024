%% Shooting Root Function
% Gives the Evaluated root and jacobian for a shooting method
function [F,J] = ShootingRootFunction(z0,TPBVP,Options,T)

% Initial Conditions
x0 = reshape(z0(1:TPBVP.nx*TPBVP.N),[TPBVP.nx,TPBVP.N]); % Extract Initial Conditions
p = z0(TPBVP.N*TPBVP.nx+1:TPBVP.N*TPBVP.nx+TPBVP.np); % Extract Parameters
if Options.FreeTime
    T = z0(TPBVP.N*TPBVP.nx+TPBVP.np+1:end);
end

% Integrate Each Arc:
xf = zeros(TPBVP.nx,TPBVP.N); % Final Points
Phif = zeros(TPBVP.nx,TPBVP.nx,TPBVP.N); % Final STMs
for i = 1:TPBVP.N
    phi0 = eye(TPBVP.nx,TPBVP.nx);
    xphi0 = [x0(:,i); reshape(phi0,[TPBVP.nx^2,1])];
    % Integrate
    [~,xphi] = ode89(@(t,xphi)[TPBVP.f(xphi(1:TPBVP.nx)); reshape(TPBVP.A(xphi(1:TPBVP.nx))*reshape(xphi(TPBVP.nx+1:TPBVP.nx^2+TPBVP.nx),[TPBVP.nx,TPBVP.nx]),[TPBVP.nx^2,1])],[sum(T(1:i-1)) sum(T(1:i))],xphi0,odeset('AbsTol',Options.AbsoluteTolerance,'RelTol',Options.RelativeTolerance)); % Perform Integration
    
    % Get Final States:
    xf(:,i) = xphi(end,1:TPBVP.nx)';
    Phif(:,:,i) = reshape(xphi(end,TPBVP.nx+1:end),[TPBVP.nx,TPBVP.nx]); % Terminal STM
end

% Get Linearized Matrices:
g = TPBVP.g(x0(:,1),xf(:,end),p,sum(T));
H0 = TPBVP.H0(x0(:,1),xf(:,end),p,sum(T));
Hf = TPBVP.Hf(x0(:,1),xf(:,end),p,sum(T));
I = TPBVP.I(x0(:,1),xf(:,end),p,sum(T));
Jf = TPBVP.Jf(x0(:,1),xf(:,end),p,sum(T));

% Root:
F = [g; reshape(x0(:,2:end)-xf(:,1:end-1),[TPBVP.nx*(TPBVP.N-1),1])];

% Jacobian:
if TPBVP.N == 1
    J11 = H0+Hf*Phif(:,:,TPBVP.N);
else
    J11 = [H0 zeros(size(g,1),(TPBVP.N-2)*TPBVP.nx) Hf*Phif(:,:,TPBVP.N)];
end
J21 = zeros(TPBVP.nx*(TPBVP.N-1),TPBVP.nx*TPBVP.N);
for i = 1:TPBVP.N-1
    J21((i-1)*TPBVP.nx+1:i*TPBVP.nx,(i-1)*TPBVP.nx+1:i*TPBVP.nx) = -Phif(:,:,i);
    J21((i-1)*TPBVP.nx+1:i*TPBVP.nx,i*TPBVP.nx+1:(i+1)*TPBVP.nx) = eye(TPBVP.nx);
end
J12 = I;
J22 = zeros(TPBVP.nx*(TPBVP.N-1),TPBVP.np);
J = [J11 J12;J21 J22];
if Options.FreeTime
    J13 = [kron(ones(1,TPBVP.N-1),Jf) Jf + Hf*TPBVP.f(xf(:,end))];
    J23 = zeros(TPBVP.nx*(TPBVP.N-1),TPBVP.N);
    for i = 1:TPBVP.N-1
        J23((i-1)*TPBVP.nx+1:i*TPBVP.nx,i) = -TPBVP.f(xf(:,i));
    end
    J = [J [J13;J23]];
end

end