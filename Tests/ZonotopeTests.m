
%% Zonotope Tests:

% Ellipsoidal Zonotope:
% G = [1 0 0;0 1 0];
% E = [1 0 0;0 1 0;0 0 0];
% A = [1 1 -.5];
% b = .5;
% R = [2 0 0;0 2 0;0 0 1];
% E = conPolyZono(c,G,E,A,b,R);
% 
% % Equivalent Ellipsoid:
% E = conPolyZono(ellipsoid([1 0;0 1],[0;0]));
% 
% % Quadratic Maps:
% Q{1} = [2 1;1 -1];
% Q{2} = [-1 1;1 -1];
% 
% EQ = quadMap(E,Q);
% 
% % Plot:
% N = 10000;
% S = [cos(linspace(0,2*pi,N));sin(linspace(0,2*pi,N))];
% ES = zeros(size(S));
% for i = 1:length(S)
%     ES(:,i) = EQ.G*prod(rand*S(:,i).^EQ.E(1:2,:))';
% end
% plot(ES(1,:),ES(2,:),'.');
% 
% % Circular Mesh:
% N1 = 100;
% N2 = 100;
% R = linspace(0,1,N1);
% theta = linspace(0,2*pi,N2);
% [RR,TT] = meshgrid(R,theta);

% Test Minkowski of two circularly constrained poly zonos:
% Sets:
c1 = [0;0];
A1 = [1 2 0;0 2 1];
E1 = [2 1 0;0 1 2];
c2 = [0;0];
A2 = [2 0 2;0 1 2];
E2 = [3 1 0;0 0 1];


function plotZono(c,A,E)

% Circular Mesh:
N1 = 100;
N2 = 100;
R = linspace(0,1,N1);
theta = linspace(0,2*pi,N2);
[RR,TT] = meshgrid(R,theta);

% P = zeros(size(A,1),N1*N2);
figure; hold on;
for i = 1:N1
    for j = 1:N2
        a = [R(i)*cos(theta(j));R(i)*sin(theta(j))];
        P = c + A*prod(a.^E)';
        plot(P(1),P(2),'k.');
    end
    i
end
end

