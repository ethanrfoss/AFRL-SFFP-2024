
function f = EOMCR3BPEphemeris(System,Options)

% Options:
DefaultOptions.Epoch = '2024-05-01 12:00:00 TDB';
DefaultOptions.Primary = {'Earth','Moon'};
DefaultOptions.Perturber = {'Sun','Venus'};
if nargin < 2
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end

% Load Kernels:
LoadKernels;

% Function Handle:
f = @(x,t) EOM(System,x,t,Options);

end

function xdot = EOM(System,x,t,Options)

% Get Initial Time:
EpochTime = cspice_str2et(Options.Epoch);
TimeSeconds = t*System.TU + EpochTime;

% Get Gravitational Data for Each Body:
PrimaryGM = arrayfun(@(B) GravitationalParameter(B),Options.Primary);
PerturberGM = arrayfun(@(B) GravitationalParameter(B),Options.Perturber);

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
PerturberState{1}(1:3) = PerturberState{1}(1:3)/(System.DU/1000); PerturberState{2}(1:3) = PerturberState{2}(1:3)/(System.DU/1000);
PerturberState{1}(4:6) = PerturberState{1}(4:6)*(System.TU)/(System.DU/1000); PerturberState{2}(4:6) = PerturberState{2}(4:6)*(System.TU)/(System.DU/1000);
PrimaryGM = PrimaryGM*(System.TU^2)/(System.DU/1000)^3;
PerturberGM = PerturberGM*(System.TU^2)/(System.DU/1000)^3;

% Coordinate Frame Origin Position:
R = (PrimaryGM(1)*PrimaryPVAJ{1}(1:3) + PrimaryGM(2)*PrimaryPVAJ{2}(1:3))/(PrimaryGM(1)+PrimaryGM(2));
Rdot = (PrimaryGM(1)*PrimaryPVAJ{1}(4:6) + PrimaryGM(2)*PrimaryPVAJ{2}(4:6))/(PrimaryGM(1)+PrimaryGM(2));
Rddot = (PrimaryGM(1)*PrimaryPVAJ{1}(7:9) + PrimaryGM(2)*PrimaryPVAJ{2}(7:9))/(PrimaryGM(1)+PrimaryGM(2));

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

% Get State Vector in Inertial Frame:
RotatingPosition = x(1:3);
RotatingVelocity = x(4:6);
InertialPosition = M*RotatingPosition + R;

% Get Accelerations due to bodies:
A = @(GM,r) GM/norm(r)^3*r;
InertialAcceleration = A(PrimaryGM(1),PrimaryPVAJ{1}(1:3)-InertialPosition(1:3));
InertialAcceleration = InertialAcceleration + A(PrimaryGM(2),PrimaryPVAJ{2}(1:3)-InertialPosition(1:3));
for i = 1:length(Options.Perturber)
    InertialAcceleration = InertialAcceleration + A(PerturberGM(i),PerturberState{i}(1:3)-InertialPosition(1:3));
end

% Get Acceleration in Rotating Frame:
RotatingAcceleration = M'*InertialAcceleration-M'*Rddot-M'*Mddot*RotatingPosition-2*M'*Mdot*RotatingVelocity;

% Rate of Change:
xdot = [RotatingVelocity;RotatingAcceleration];

end