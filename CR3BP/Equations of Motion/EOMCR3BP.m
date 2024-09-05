%% EOMCR3BP - Equations of Motion for the CR3BP
% 
% Inputs:
%    System - Structure Containing the System variables
%    t - Current Time (Does not matter)
%    x - Current State
%    Options - Structure Array for EOM Options. The Options and their
%    default values are:
%        Regularize (false) - Dynamic Regularization
% 
% Outputs:
%    f - Matlab Function Handle for Equations of Motion
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function f = EOMCR3BP(System,Options)

% Options:
DefaultOptions.Regularize = false;
DefaultOptions.Planar = false;
DefaultOptions.Ephemeris = false;
DefaultOptions.Spacecraft = false;
DefaultOptions.TensorOrder = 1;
if nargin < 2
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

if Options.TensorOrder ~= 1

    if Options.Ephemeris
        
        % Options:
        DefaultOptions.Epoch = '2024-05-01 12:00:00 TDB';
        DefaultOptions.Primary = {'Earth','Moon'};
        DefaultOptions.Perturber = {'Sun','Venus'};
        Options = MergeStructs(Options,DefaultOptions);
        
        % Load Kernels:
        LoadKernels;
        
        % Get Gravitational Data for Each Body:
        System.GM = [arrayfun(@(B) GravitationalParameter(B),Options.Primary) arrayfun(@(B) GravitationalParameter(B),Options.Perturber)]*(System.TU^2)/(System.DU/1000)^3;

        if Options.Spacecraft
            f = @(x,t) STTEphemerisInputs(System,x,t(1),Options);
        else
            f = @(x,t) STTEphemeris(System,x,t(1),Options);
        end
    elseif Options.Planar
        if Options.Spacecraft
            f = @(x,t) STTCR3BPPlanarInputs(System,x,Options);
        else
            f = @(x,t) STTCR3BPPlanar(System,x,Options);
        end
    else
        if Options.Spacecraft
            f = @(x,t) STTCR3BPNonplanarInputs(System,x,Options);
        else
            f = @(x,t) STTCR3BPNonplanar(System,x,Options);
        end
    end
else

    if Options.Ephemeris
        
        % Options:
        DefaultOptions.Epoch = '2024-05-01 12:00:00 TDB';
        DefaultOptions.Primary = {'Earth','Moon'};
        DefaultOptions.Perturber = {'Sun','Venus'};
        Options = MergeStructs(Options,DefaultOptions);
        
        % Load Kernels:
        LoadKernels;
    
        % Get Gravitational Data for Each Body:
        System.GM = [arrayfun(@(B) GravitationalParameter(B),Options.Primary) arrayfun(@(B) GravitationalParameter(B),Options.Perturber)]*(System.TU^2)/(System.DU/1000)^3;
        
        % Function Handle:
        f = @(x,t) EOMEphemeris(System,x,t(1),Options);

    elseif Options.Planar
        f = EOMCR3BPPlanar(System);
    else
        f = EOMCR3BPNonplanar(System);
    end
    
    if Options.Spacecraft
        if Options.Planar
            f = @(x,u,t) f(x,t) + System.Spacecraft.Thrust/System.Spacecraft.WetMass*[0;0;u(1);u(2)];
        else
            f = @(x,u,t) f(x,t) + System.Spacecraft.Thrust/System.Spacecraft.WetMass*[0;0;0;u(1);u(2);u(3)];
        end
    end
    
    if Options.Regularize
        f = @(x,varargin) 1/SundmanTransformation(System,x)*f(x,varargin{:});
    end

end

end

function xdot = EOMEphemeris(System,x,t,Options)

% Time Varying Quantities:
[M,Mdot,Mddot,R,Rddot,r] = TimeVaryingQuantities(System,t,Options);

% State Varying Quantities:
xdot = StateVaryingQuantities(System,x,M,Mdot,Mddot,R,Rddot,r,Options);

end

function D = STTEphemeris(System,x,t,Options)

% Time Varying Quantities:
[M,Mdot,Mddot,R,Rddot,r] = TimeVaryingQuantities(System,t,Options);

% 1st Order:
A{1} = StateVaryingQuantities(System,x,M,Mdot,Mddot,R,Rddot,r,Options);

if Options.TensorOrder >= 2
    
    d = r;
    for i = 1:length(r)
        d{i} = x(1:3)-M'*(r{i}-R);
    end
    A{2} = zeros(6,6);
    A{2}(1:3,4:6) = eye(3,3);
    for i = 1:length(r)
        A{2}(4:6,1:3) = A{2}(4:6,1:3)+System.GM(i)*(3/norm(d{i})^5*d{i}*d{i}' - 1/norm(d{i})^3*eye(3,3));
    end
    A{2}(4:6,1:3) = A{2}(4:6,1:3)-M'*Mddot;
    A{2}(4:6,4:6) = -2*M'*Mdot;

end

if Options.TensorOrder >= 3

    A{3} = zeros(6,6,6);
    for i = 1:length(r)
        A{3}(4:6,1:3,1:3) = A{3}(4:6,1:3,1:3) + System.GM(i)*(-15/norm(d{i})^7*tensorprod(d{i}*d{i}',d{i}',[3],[1])+3/norm(d{i})^5*(tensorprod(eye(3,3),d{i}',[3],[1]) + permute(tensorprod(eye(3,3),d{i}',[3],[1]),[3 1 2]) + permute(tensorprod(eye(3,3),d{i}',[3],[1]),[2 3 1])));
    end

end

if Options.TensorOrder >= 4

    A{4} = zeros(6,6,6,6);
    for i = 1:length(r)
        T = tensorprod(tensorprod(eye(3,3),d{i}',[3],[1]),d{i}',[4],[1]);
        P(:,:,1,1) = [3 0 0;0 1 0;0 0 1]; P(:,:,2,1) = [0 1 0;1 0 0;0 0 0]; P(:,:,3,1) = [0 0 1;0 0 0;1 0 0];
        P(:,:,1,2) = [0 1 0;1 0 0;0 0 0]; P(:,:,2,2) = [1 0 0;0 3 0;0 0 1]; P(:,:,3,2) = [0 0 0;0 0 1;0 1 0];
        P(:,:,1,3) = [0 0 1;0 0 0;1 0 0]; P(:,:,2,3) = [0 0 0;0 0 1;0 1 0]; P(:,:,3,3) = [1 0 0;0 1 0;0 0 3];
        A{4}(4:6,1:3,1:3,1:3) = A{4}(4:6,1:3,1:3,1:3) + System.GM(i)*(105/norm(d{i})^9*tensorprod(tensorprod(d{i}*d{i}',d{i}',[3],[1]),d{i}',[4],[1]) - 15/norm(d{i})^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(d{i})^5*P);
    end

end

D = sparse(Detensify(A));

end

function D = STTEphemerisInputs(System,z,t,Options)

x = z(1:6);
u = z(7:9);
B = [zeros(3,3);eye(3,3)]*System.Spacecraft.Thrust/System.Spacecraft.WetMass;

% Time Varying Quantities:
[M,Mdot,Mddot,R,Rddot,r] = TimeVaryingQuantities(System,t,Options);

% 1st Order:
A{1} = StateVaryingQuantities(System,x,M,Mdot,Mddot,R,Rddot,r,Options);
A{1} = [A{1} + B*u;0;0;0];

if Options.TensorOrder >= 2
    
    d = r;
    for i = 1:length(r)
        d{i} = x(1:3)-M'*(r{i}-R);
    end
    A{2} = zeros(9,9);
    A{2}(1:3,4:6) = eye(3,3);
    for i = 1:length(r)
        A{2}(4:6,1:3) = A{2}(4:6,1:3)+System.GM(i)*(3/norm(d{i})^5*d{i}*d{i}' - 1/norm(d{i})^3*eye(3,3));
    end
    A{2}(4:6,1:3) = A{2}(4:6,1:3)-M'*Mddot;
    A{2}(4:6,4:6) = -2*M'*Mdot;
    A{2}(1:6,7:9) = B;

end

if Options.TensorOrder >= 3

    A{3} = zeros(9,9,9);
    for i = 1:length(r)
        A{3}(4:6,1:3,1:3) = A{3}(4:6,1:3,1:3) + System.GM(i)*(-15/norm(d{i})^7*tensorprod(d{i}*d{i}',d{i}',[3],[1])+3/norm(d{i})^5*(tensorprod(eye(3,3),d{i}',[3],[1]) + permute(tensorprod(eye(3,3),d{i}',[3],[1]),[3 1 2]) + permute(tensorprod(eye(3,3),d{i}',[3],[1]),[2 3 1])));
    end

end

if Options.TensorOrder >= 4

    A{4} = zeros(9,9,9,9);
    for i = 1:length(r)
        T = tensorprod(tensorprod(eye(3,3),d{i}',[3],[1]),d{i}',[4],[1]);
        P(:,:,1,1) = [3 0 0;0 1 0;0 0 1]; P(:,:,2,1) = [0 1 0;1 0 0;0 0 0]; P(:,:,3,1) = [0 0 1;0 0 0;1 0 0];
        P(:,:,1,2) = [0 1 0;1 0 0;0 0 0]; P(:,:,2,2) = [1 0 0;0 3 0;0 0 1]; P(:,:,3,2) = [0 0 0;0 0 1;0 1 0];
        P(:,:,1,3) = [0 0 1;0 0 0;1 0 0]; P(:,:,2,3) = [0 0 0;0 0 1;0 1 0]; P(:,:,3,3) = [1 0 0;0 1 0;0 0 3];
        A{4}(4:6,1:3,1:3,1:3) = A{4}(4:6,1:3,1:3,1:3) + System.GM(i)*(105/norm(d{i})^9*tensorprod(tensorprod(d{i}*d{i}',d{i}',[3],[1]),d{i}',[4],[1]) - 15/norm(d{i})^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(d{i})^5*P);
    end

end

D = sparse(Detensify(A));

end

% function V = Tensors(System,t,x,V,Options)
% 
% [M,Mdot,Mddot,R,Rddot,r] = TimeVaryingQuantities(System,t,Options);
% V = V(x,M,Mdot,Mddot,R,Rddot,r{:});
% 
% end

function [M,Mdot,Mddot,R,Rddot,r] = TimeVaryingQuantities(System,t,Options)

% Get Initial Time:
EpochTime = cspice_str2et(Options.Epoch);
TimeSeconds = t*System.TU + EpochTime;

% Pos, Vel, Acc, and Jerk:
PrimaryPVAJ{1} = SolarSystemInertialPVAJ(Options.Primary{1},TimeSeconds);
PrimaryPVAJ{2} = SolarSystemInertialPVAJ(Options.Primary{2},TimeSeconds);
for i = 1:length(Options.Perturber)
    PerturberState{i} = SolarSystemInertialState(Options.Perturber{i},TimeSeconds);
end

% Normalize:
PrimaryPVAJ{1}(1:3) = PrimaryPVAJ{1}(1:3)/(System.DU/1000); PrimaryPVAJ{2}(1:3) = PrimaryPVAJ{2}(1:3)/(System.DU/1000);
PrimaryPVAJ{1}(4:6) = PrimaryPVAJ{1}(4:6)*(System.TU)/(System.DU/1000); PrimaryPVAJ{2}(4:6) = PrimaryPVAJ{2}(4:6)*(System.TU)/(System.DU/1000);
PrimaryPVAJ{1}(7:9) = PrimaryPVAJ{1}(7:9)*(System.TU)^2/(System.DU/1000); PrimaryPVAJ{2}(7:9) = PrimaryPVAJ{2}(7:9)*(System.TU)^2/(System.DU/1000);
PrimaryPVAJ{1}(10:12) = PrimaryPVAJ{1}(10:12)*(System.TU)^3/(System.DU/1000); PrimaryPVAJ{2}(10:12) = PrimaryPVAJ{2}(10:12)*(System.TU)^3/(System.DU/1000);
r{1} = PrimaryPVAJ{1}(1:3);
r{2} = PrimaryPVAJ{2}(1:3);
for i = 1:length(PerturberState)
    PerturberState{i}(1:3) = PerturberState{i}(1:3)/(System.DU/1000);
    PerturberState{i}(4:6) = PerturberState{i}(4:6)*(System.TU)/(System.DU/1000);
    r{i+2} = PerturberState{i}(1:3);
end

% Coordinate Frame Origin Position:
R = (System.GM(1)*PrimaryPVAJ{1}(1:3) + System.GM(2)*PrimaryPVAJ{2}(1:3))/(System.GM(1)+System.GM(2));
% Rdot = (PrimaryGM(1)*PrimaryPVAJ{1}(4:6) + PrimaryGM(2)*PrimaryPVAJ{2}(4:6))/(PrimaryGM(1)+PrimaryGM(2));
Rddot = (System.GM(1)*PrimaryPVAJ{1}(7:9) + System.GM(2)*PrimaryPVAJ{2}(7:9))/(System.GM(1)+System.GM(2));

% Coordinate Frame Orientation:
% Functions for Unit Vector Rates of Change:
vhat = @(v) v/norm(v);
vhatdot = @(v,vdot) vdot/norm(v) - v'*vdot/norm(v)^3*v;
vhatddot = @(v,vdot,vddot) vddot/norm(v) -2*v'*vdot/norm(v)^3*vdot -(vdot'*vdot+v'*vddot)/norm(v)^3*v+3*(v'*vdot)^2/norm(v)^5*v;

% Unit Vectors:
xhat = vhat(PrimaryPVAJ{2}(1:3) - PrimaryPVAJ{1}(1:3));
xhatdot = vhatdot(PrimaryPVAJ{2}(1:3) - PrimaryPVAJ{1}(1:3),PrimaryPVAJ{2}(4:6) - PrimaryPVAJ{1}(4:6));
xhatddot = vhatddot(PrimaryPVAJ{2}(1:3) - PrimaryPVAJ{1}(1:3),PrimaryPVAJ{2}(4:6) - PrimaryPVAJ{1}(4:6),PrimaryPVAJ{2}(7:9) - PrimaryPVAJ{1}(7:9));
zhat = vhat(cross(xhat,PrimaryPVAJ{2}(4:6)-PrimaryPVAJ{1}(4:6)));
zhatdot = vhatdot(cross(xhat,PrimaryPVAJ{2}(4:6)-PrimaryPVAJ{1}(4:6)),cross(xhatdot,PrimaryPVAJ{2}(4:6)-PrimaryPVAJ{1}(4:6))+cross(xhat,PrimaryPVAJ{2}(7:9)-PrimaryPVAJ{1}(7:9)));
zhatddot = vhatddot(cross(xhat,PrimaryPVAJ{2}(4:6)-PrimaryPVAJ{1}(4:6)),cross(xhatdot,PrimaryPVAJ{2}(4:6)-PrimaryPVAJ{1}(4:6))+cross(xhat,PrimaryPVAJ{2}(7:9)-PrimaryPVAJ{1}(7:9)),cross(xhatddot,PrimaryPVAJ{2}(4:6)-PrimaryPVAJ{1}(4:6))+2*cross(xhatdot,PrimaryPVAJ{2}(7:9)-PrimaryPVAJ{1}(7:9))+cross(xhat,PrimaryPVAJ{2}(10:12)-PrimaryPVAJ{1}(10:12)));
yhat = cross(zhat,xhat);
yhatdot = cross(zhatdot,xhat) + cross(zhat,xhatdot);
yhatddot = cross(zhatddot,xhat) + 2*cross(zhatdot,xhatdot) + cross(zhat,xhatddot);

% Transformation from Rotating to Inertial Frame:
M = [xhat yhat zhat];
Mdot = [xhatdot yhatdot zhatdot];
Mddot = [xhatddot yhatddot zhatddot];

end

function xdot = StateVaryingQuantities(System,x,M,Mdot,Mddot,R,Rddot,r,Options)

% Get State Vector in Inertial Frame:
RotatingPosition = x(1:3);
RotatingVelocity = x(4:6);
InertialPosition = M*RotatingPosition + R;

% Get Accelerations due to bodies:
A = @(GM,r) GM/sqrt(r.'*r)^3*r;
% InertialAcceleration = A(PrimaryGM(1),PrimaryPVAJ{1}(1:3)-InertialPosition(1:3));
% InertialAcceleration = InertialAcceleration + A(PrimaryGM(2),PrimaryPVAJ{2}(1:3)-InertialPosition(1:3));
% for i = 1:length(Options.Perturber)
%     InertialAcceleration = InertialAcceleration + A(PerturberGM(i),PerturberState{i}(1:3)-InertialPosition(1:3));
% end
InertialAcceleration = zeros(3,1);
for i = 1:length(Options.Perturber)+2
    InertialAcceleration = InertialAcceleration + A(System.GM(i),r{i}-InertialPosition);
end

% Get Acceleration in Rotating Frame:
RotatingAcceleration = M.'*InertialAcceleration-M.'*Rddot-M.'*Mddot*RotatingPosition-2*M.'*Mdot*RotatingVelocity;

% Rate of Change:
xdot = [RotatingVelocity;RotatingAcceleration];

end

function f = EOMCR3BPPlanar(System)

r1 = @(x) ((x(1)+System.mu)^2+x(2)^2)^(1/2);
r2 = @(x) ((x(1)+System.mu-1)^2+x(2)^2)^(1/2);
f = @(x,t) [x(3)
            x(4)
            -(1-System.mu)*(x(1)+System.mu)/r1(x)^3-System.mu*(x(1)+System.mu-1)/r2(x)^3+2*x(4)+x(1)     
            -((1-System.mu)/r1(x)^3+System.mu/r2(x)^3)*x(2)-2*x(3)+x(2)];

end

function f = EOMCR3BPNonplanar(System)

r1 = @(x) ((x(1)+System.mu)^2+x(2)^2+x(3)^2)^(1/2);
r2 = @(x) ((x(1)+System.mu-1)^2+x(2)^2+x(3)^2)^(1/2);
f = @(x,t) [x(4)
            x(5)
            x(6)
            -(1-System.mu)*(x(1)+System.mu)/r1(x)^3-System.mu*(x(1)+System.mu-1)/r2(x)^3+2*x(5)+x(1)     
            -((1-System.mu)/r1(x)^3+System.mu/r2(x)^3)*x(2)-2*x(4)+x(2)
            -((1-System.mu)/r1(x)^3+System.mu/r2(x)^3)*x(3)];

end

function D = STTCR3BPPlanar(System,x,Options)

% 1st Order:
% r1 = ((x(1)+System.mu)^2+x(2)^2)^(1/2);
% r2 = ((x(1)+System.mu-1)^2+x(2)^2)^(1/2);
r1 = [x(1)+System.mu;x(2)];
r2 = [x(1)+System.mu-1;x(2)];
A{1} =  [x(3)
            x(4)
            -(1-System.mu)*(x(1)+System.mu)/norm(r1)^3-System.mu*(x(1)+System.mu-1)/norm(r2)^3+2*x(4)+x(1)     
            -((1-System.mu)/norm(r1)^3+System.mu/norm(r2)^3)*x(2)-2*x(3)+x(2)];

if Options.TensorOrder >= 2
    
    A{2} = zeros(4,4);
    A{2}(1:2,3:4) = eye(2,2);
    A{2}(3:4,1:2) = A{2}(3:4,1:2)+(1-System.mu)*(3/norm(r1)^5*r1*r1' - 1/norm(r1)^3*eye(2,2));
    A{2}(3:4,1:2) = A{2}(3:4,1:2)+(System.mu)*(3/norm(r2)^5*r2*r2' - 1/norm(r2)^3*eye(2,2));
    A{2}(3:4,1:2) = A{2}(3:4,1:2)+[1 0;0 1];
    A{2}(3:4,3:4) = [0 2;-2 0];

end

if Options.TensorOrder >= 3

    A{3} = zeros(4,4,4);
    A{3}(3:4,1:2,1:2) = A{3}(3:4,1:2,1:2) + (1-System.mu)*(-15/norm(r1)^7*tensorprod(r1*r1',r1',[3],[1])+3/norm(r1)^5*(tensorprod(eye(2,2),r1',[3],[1]) + permute(tensorprod(eye(2,2),r1',[3],[1]),[3 1 2]) + permute(tensorprod(eye(2,2),r1',[3],[1]),[2 3 1])));
    A{3}(3:4,1:2,1:2) = A{3}(3:4,1:2,1:2) + (System.mu)*(-15/norm(r2)^7*tensorprod(r2*r2',r2',[3],[1])+3/norm(r2)^5*(tensorprod(eye(2,2),r2',[3],[1]) + permute(tensorprod(eye(2,2),r2',[3],[1]),[3 1 2]) + permute(tensorprod(eye(2,2),r2',[3],[1]),[2 3 1])));

end

if Options.TensorOrder >= 4

    A{4} = zeros(4,4,4,4);
    P(:,:,1,1) = [3 0;0 1]; P(:,:,2,1) = [0 1;1 0];
    P(:,:,1,2) = [0 1;1 0]; P(:,:,2,2) = [1 0;0 3];
    T = tensorprod(tensorprod(eye(2,2),r1',[3],[1]),r1',[4],[1]);
    A{4}(3:4,1:2,1:2,1:2) = A{4}(3:4,1:2,1:2,1:2) + (1-System.mu)*(105/norm(r1)^9*tensorprod(tensorprod(r1*r1',r1',[3],[1]),r1',[4],[1]) - 15/norm(r1)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r1)^5*P);
    T = tensorprod(tensorprod(eye(2,2),r2',[3],[1]),r2',[4],[1]);
    A{4}(3:4,1:2,1:2,1:2) = A{4}(3:4,1:2,1:2,1:2)+ (System.mu)*(105/norm(r2)^9*tensorprod(tensorprod(r2*r2',r2',[3],[1]),r2',[4],[1]) - 15/norm(r2)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r2)^5*P);

end

D = (Detensify(A));

end

function D = STTCR3BPNonplanar(System,x,Options)

% 1st Order:
% r1 = ((x(1)+System.mu)^2+x(2)^2)^(1/2);
% r2 = ((x(1)+System.mu-1)^2+x(2)^2)^(1/2);
r1 = [x(1)+System.mu;x(2);x(3)];
r2 = [x(1)+System.mu-1;x(2);x(3)];
A{1} =  [x(4)
            x(5)
            x(6)
            -(1-System.mu)*(x(1)+System.mu)/norm(r1)^3-System.mu*(x(1)+System.mu-1)/norm(r2)^3+2*x(5)+x(1)     
            -((1-System.mu)/norm(r1)^3+System.mu/norm(r2)^3)*x(2)-2*x(4)+x(2)
            -((1-System.mu)/norm(r1)^3+System.mu/norm(r2)^3)*x(3)];

if Options.TensorOrder >= 2
    
    A{2} = zeros(6,6);
    A{2}(1:3,4:6) = eye(3,3);
    A{2}(4:6,1:3) = A{2}(4:6,1:3)+(1-System.mu)*(3/norm(r1)^5*r1*r1' - 1/norm(r1)^3*eye(3,3));
    A{2}(4:6,1:3) = A{2}(4:6,1:3)+(System.mu)*(3/norm(r2)^5*r2*r2' - 1/norm(r2)^3*eye(3,3));
    A{2}(4:6,1:3) = A{2}(4:6,1:3)+[1 0 0;0 1 0;0 0 0];
    A{2}(4:6,4:6) = [0 2 0;-2 0 0;0 0 0];

end

if Options.TensorOrder >= 3

    A{3} = zeros(6,6,6);
    A{3}(4:6,1:3,1:3) = A{3}(4:6,1:3,1:3) + (1-System.mu)*(-15/norm(r1)^7*tensorprod(r1*r1',r1',[3],[1])+3/norm(r1)^5*(tensorprod(eye(3,3),r1',[3],[1]) + permute(tensorprod(eye(3,3),r1',[3],[1]),[3 1 2]) + permute(tensorprod(eye(3,3),r1',[3],[1]),[2 3 1])));
    A{3}(4:6,1:3,1:3) = A{3}(4:6,1:3,1:3) + (System.mu)*(-15/norm(r2)^7*tensorprod(r2*r2',r2',[3],[1])+3/norm(r2)^5*(tensorprod(eye(3,3),r2',[3],[1]) + permute(tensorprod(eye(3,3),r2',[3],[1]),[3 1 2]) + permute(tensorprod(eye(3,3),r2',[3],[1]),[2 3 1])));

end

if Options.TensorOrder >= 4

    A{4} = zeros(6,6,6,6);
    P(:,:,1,1) = [3 0 0;0 1 0;0 0 1]; P(:,:,2,1) = [0 1 0;1 0 0;0 0 0]; P(:,:,3,1) = [0 0 1;0 0 0;1 0 0];
    P(:,:,1,2) = [0 1 0;1 0 0;0 0 0]; P(:,:,2,2) = [1 0 0;0 3 0;0 0 1]; P(:,:,3,2) = [0 0 0;0 0 1;0 1 0];
    P(:,:,1,3) = [0 0 1;0 0 0;1 0 0]; P(:,:,2,3) = [0 0 0;0 0 1;0 1 0]; P(:,:,3,3) = [1 0 0;0 1 0;0 0 3];
    T = tensorprod(tensorprod(eye(3,3),r1',[3],[1]),r1',[4],[1]);
    A{4}(4:6,1:3,1:3,1:3) = A{4}(4:6,1:3,1:3,1:3) + (1-System.mu)*(105/norm(r1)^9*tensorprod(tensorprod(r1*r1',r1',[3],[1]),r1',[4],[1]) - 15/norm(r1)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r1)^5*P);
    T = tensorprod(tensorprod(eye(3,3),r2',[3],[1]),r2',[4],[1]);
    A{4}(4:6,1:3,1:3,1:3) = A{4}(4:6,1:3,1:3,1:3) + (System.mu)*(105/norm(r2)^9*tensorprod(tensorprod(r2*r2',r2',[3],[1]),r2',[4],[1]) - 15/norm(r2)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r2)^5*P);

end

D = sparse(Detensify(A));

end

function D = STTCR3BPPlanarInputs(System,z,Options)

% 1st Order:
% r1 = ((x(1)+System.mu)^2+x(2)^2)^(1/2);
% r2 = ((x(1)+System.mu-1)^2+x(2)^2)^(1/2);
x = z(1:4);
u = z(5:6);
B = [zeros(2,2);eye(2,2)]*System.Spacecraft.Thrust/System.Spacecraft.WetMass;

r1 = [x(1)+System.mu;x(2)];
r2 = [x(1)+System.mu-1;x(2)];
A{1} =  [x(3)
            x(4)
            -(1-System.mu)*(x(1)+System.mu)/norm(r1)^3-System.mu*(x(1)+System.mu-1)/norm(r2)^3+2*x(4)+x(1)     
            -((1-System.mu)/norm(r1)^3+System.mu/norm(r2)^3)*x(2)-2*x(3)+x(2)];
A{1} = [A{1} + B*u;0;0];

if Options.TensorOrder >= 2
    
    A{2} = zeros(6,6);
    A{2}(1:2,3:4) = eye(2,2);
    A{2}(3:4,1:2) = A{2}(3:4,1:2)+(1-System.mu)*(3/norm(r1)^5*r1*r1' - 1/norm(r1)^3*eye(2,2));
    A{2}(3:4,1:2) = A{2}(3:4,1:2)+(System.mu)*(3/norm(r2)^5*r2*r2' - 1/norm(r2)^3*eye(2,2));
    A{2}(3:4,1:2) = A{2}(3:4,1:2)+[1 0;0 1];
    A{2}(3:4,3:4) = [0 2;-2 0];
    A{2}(1:4,5:6) = B;

end

if Options.TensorOrder >= 3

    A{3} = zeros(6,6,6);
    A{3}(3:4,1:2,1:2) = A{3}(3:4,1:2,1:2) + (1-System.mu)*(-15/norm(r1)^7*tensorprod(r1*r1',r1',[3],[1])+3/norm(r1)^5*(tensorprod(eye(2,2),r1',[3],[1]) + permute(tensorprod(eye(2,2),r1',[3],[1]),[3 1 2]) + permute(tensorprod(eye(2,2),r1',[3],[1]),[2 3 1])));
    A{3}(3:4,1:2,1:2) = A{3}(3:4,1:2,1:2) + (System.mu)*(-15/norm(r2)^7*tensorprod(r2*r2',r2',[3],[1])+3/norm(r2)^5*(tensorprod(eye(2,2),r2',[3],[1]) + permute(tensorprod(eye(2,2),r2',[3],[1]),[3 1 2]) + permute(tensorprod(eye(2,2),r2',[3],[1]),[2 3 1])));

end

if Options.TensorOrder >= 4

    A{4} = zeros(6,6,6,6);
    P(:,:,1,1) = [3 0;0 1]; P(:,:,2,1) = [0 1;1 0];
    P(:,:,1,2) = [0 1;1 0]; P(:,:,2,2) = [1 0;0 3];
    T = tensorprod(tensorprod(eye(2,2),r1',[3],[1]),r1',[4],[1]);
    A{4}(3:4,1:2,1:2,1:2) = A{4}(3:4,1:2,1:2,1:2) + (1-System.mu)*(105/norm(r1)^9*tensorprod(tensorprod(r1*r1',r1',[3],[1]),r1',[4],[1]) - 15/norm(r1)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r1)^5*P);
    T = tensorprod(tensorprod(eye(2,2),r2',[3],[1]),r2',[4],[1]);
    A{4}(3:4,1:2,1:2,1:2) = A{4}(3:4,1:2,1:2,1:2)+ (System.mu)*(105/norm(r2)^9*tensorprod(tensorprod(r2*r2',r2',[3],[1]),r2',[4],[1]) - 15/norm(r2)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r2)^5*P);

end

D = (Detensify(A));

end

function D = STTCR3BPNonplanarInputs(System,z,Options)

% 1st Order:
% r1 = ((x(1)+System.mu)^2+x(2)^2)^(1/2);
% r2 = ((x(1)+System.mu-1)^2+x(2)^2)^(1/2);
x = z(1:6);
u = z(7:9);
B = [zeros(3,3);eye(3,3)]*System.Spacecraft.Thrust/System.Spacecraft.WetMass;

r1 = [x(1)+System.mu;x(2);x(3)];
r2 = [x(1)+System.mu-1;x(2);x(3)];
A{1} =  [x(4)
            x(5)
            x(6)
            -(1-System.mu)*(x(1)+System.mu)/norm(r1)^3-System.mu*(x(1)+System.mu-1)/norm(r2)^3+2*x(5)+x(1)     
            -((1-System.mu)/norm(r1)^3+System.mu/norm(r2)^3)*x(2)-2*x(4)+x(2)
            -((1-System.mu)/norm(r1)^3+System.mu/norm(r2)^3)*x(3)];
A{1} = [A{1} + B*u;0;0;0];

if Options.TensorOrder >= 2
    
    A{2} = zeros(9,9);
    A{2}(1:3,4:6) = eye(3,3);
    A{2}(4:6,1:3) = A{2}(4:6,1:3)+(1-System.mu)*(3/norm(r1)^5*r1*r1' - 1/norm(r1)^3*eye(3,3));
    A{2}(4:6,1:3) = A{2}(4:6,1:3)+(System.mu)*(3/norm(r2)^5*r2*r2' - 1/norm(r2)^3*eye(3,3));
    A{2}(4:6,1:3) = A{2}(4:6,1:3)+[1 0 0;0 1 0;0 0 0];
    A{2}(4:6,4:6) = [0 2 0;-2 0 0;0 0 0];
    A{2}(1:6,7:9) = B;

end

if Options.TensorOrder >= 3

    A{3} = zeros(9,9,9);
    A{3}(4:6,1:3,1:3) = A{3}(4:6,1:3,1:3) + (1-System.mu)*(-15/norm(r1)^7*tensorprod(r1*r1',r1',[3],[1])+3/norm(r1)^5*(tensorprod(eye(3,3),r1',[3],[1]) + permute(tensorprod(eye(3,3),r1',[3],[1]),[3 1 2]) + permute(tensorprod(eye(3,3),r1',[3],[1]),[2 3 1])));
    A{3}(4:6,1:3,1:3) = A{3}(4:6,1:3,1:3) + (System.mu)*(-15/norm(r2)^7*tensorprod(r2*r2',r2',[3],[1])+3/norm(r2)^5*(tensorprod(eye(3,3),r2',[3],[1]) + permute(tensorprod(eye(3,3),r2',[3],[1]),[3 1 2]) + permute(tensorprod(eye(3,3),r2',[3],[1]),[2 3 1])));

end

if Options.TensorOrder >= 4

    A{4} = zeros(9,9,9,9);
    P(:,:,1,1) = [3 0 0;0 1 0;0 0 1]; P(:,:,2,1) = [0 1 0;1 0 0;0 0 0]; P(:,:,3,1) = [0 0 1;0 0 0;1 0 0];
    P(:,:,1,2) = [0 1 0;1 0 0;0 0 0]; P(:,:,2,2) = [1 0 0;0 3 0;0 0 1]; P(:,:,3,2) = [0 0 0;0 0 1;0 1 0];
    P(:,:,1,3) = [0 0 1;0 0 0;1 0 0]; P(:,:,2,3) = [0 0 0;0 0 1;0 1 0]; P(:,:,3,3) = [1 0 0;0 1 0;0 0 3];
    T = tensorprod(tensorprod(eye(3,3),r1',[3],[1]),r1',[4],[1]);
    A{4}(4:6,1:3,1:3,1:3) = A{4}(4:6,1:3,1:3,1:3) + (1-System.mu)*(105/norm(r1)^9*tensorprod(tensorprod(r1*r1',r1',[3],[1]),r1',[4],[1]) - 15/norm(r1)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r1)^5*P);
    T = tensorprod(tensorprod(eye(3,3),r2',[3],[1]),r2',[4],[1]);
    A{4}(4:6,1:3,1:3,1:3) = A{4}(4:6,1:3,1:3,1:3) + (System.mu)*(105/norm(r2)^9*tensorprod(tensorprod(r2*r2',r2',[3],[1]),r2',[4],[1]) - 15/norm(r2)^7*(permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2])) + 3/norm(r2)^5*P);

end

D = sparse(Detensify(A));

end