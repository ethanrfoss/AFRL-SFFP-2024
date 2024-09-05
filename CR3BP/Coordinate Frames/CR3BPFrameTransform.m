%% CR3BPFrameTransformation: Function for transforming state data to
% different coordinate systems.
%
% Inputs: 
%    System - Structure Containing the System variables, which is
%    constructed using System = InitializeCR3BP(System).
%    t - Time Data
%    x - State Data:
%    Origins - String Vector containing Origin of current system and Origin 
%    of transforming frame, in that order. Options are Earth, Barycenter,
%    and Moon
%    Frames - String vector contataing current frame and transforming
%    frame, in that order. Options are Rotating and Inertial. 
%
% Outputs:
%    X - Transformed data
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024
function X = CR3BPFrameTransform(System,t,x,Origins,Frames)

% Standard frame, defined as rotating barycentric:
Q.Rotating = @(t) [1;0;0;0];
R.Barycenter = @(t) [0;0;0];

% Additional frames and origins:
Q.Inertial = @(t) [cos(t/2);0;0;sin(t/2)];
R.Earth = @(t) [-System.mu;0;0];
R.Moon = @(t) [1-System.mu;0;0];

% Checks:
assert(all(isfield(Q,Frames)));
assert(all(isfield(R,Origins)));

% Construct Quaternion and Displacement vectors:
QT = @(t) QuaternionMultiply(QuaternionConjugate(Q.(Frames(1))(t)),Q.(Frames(2))(t));
RT = @(t) QuaternionRotate(Q.(Frames(1))(t),R.(Origins(2))(t)-R.(Origins(1))(t));

% Perform Conversion:
nx = size(x,1);
if nx == 4 % Planar Case
    x = [x(1:2,:);zeros(1,length(t));x(3:4,:);zeros(1,length(t))];
end
X = CoordinateTransformation(t,x,RT,QT);
if nx == 4
    X = [X(1:2,:);X(4:5,:)];
end

end