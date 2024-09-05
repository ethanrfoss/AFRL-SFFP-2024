%% CoordinateTransformation - Transforms a sequence of time varying states 
% to another frame represented by time functional displacements and
% rotations
%
% Inputs:
%    t - Time Vector
%    x - State vector
%    R - Time function that represents the displacement vector to the
%    transformed frame origin expressed in the original frame
%    Q - Time function that represents the quaternion rotation to the
%    transformed frame
%
% Outputs:
%    X - Transformed states
%
% Example Usage:
%    t = 0;
%    x = [.8;0;0;0;.6;0];
%    R = @(t) [-.012;0;0];
%    Q = @(t) [cos(t/2);0;0;-sin(t/2)];
%    X = CoordinateTransformation(t,x,R,Q)
function X = CoordinateTransformation(t,x,R,Q)

% Derive position and quaternion rate of change:
syms tt;
V = matlabFunction(diff(R(tt),tt),"Vars",tt);
O = matlabFunction(diff(Q(tt),tt),"Vars",tt);

% Loop Through State vector:
X = zeros(size(x));
for i = 1:length(t)
    X(1:3,i) = QuaternionRotate(Q(t(i)),x(1:3,i)-R(t(i)));
    X(4:6,i) = [zeros(3,1) eye(3)]*(QuaternionMultiply(QuaternionMultiply(Q(t(i)),[0;x(1:3,i)-R(t(i))]),QuaternionConjugate(O(t(i)))) + QuaternionMultiply(QuaternionMultiply(O(t(i)),[0;x(1:3,i)-R(t(i))]),QuaternionConjugate(Q(t(i)))) + QuaternionMultiply(QuaternionMultiply(Q(t(i)),[0;x(4:6,i)-V(t(i))]),QuaternionConjugate(Q(t(i))))); 
end  

end