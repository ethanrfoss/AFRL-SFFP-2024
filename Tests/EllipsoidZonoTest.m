
% Ellipsoid Test:

N = 2;

G = zeros(2,2*N+3);

G(1,2:2:2*N+2) = (-1).^(0:N).*pi.^(2*(0:N))./factorial(2*(0:N));
G(2,3:2:2*N+3) = (-1).^(0:N).*pi.^(2*(0:N)+1)./factorial(2*(0:N)+1);

E = [0 0 1:(2*N+1);0 ones(1,2*N+2)];

PZ = PolynomialZonotope(G,E);

P = RandomPoints(PZ,100000);

figure; hold on;
scatter(P(1,:),P(2,:),'Marker','.','MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',.1,'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',.1);
plot(polyshape(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100))));

% Ellipsoid function:
EP = Ellipsoid([1 0;0 1],[0;0]);
PZ = PolynomialZonotope(EP);

P = RandomPoints(PZ,100000);

figure; hold on;
scatter(P(1,:),P(2,:),'Marker','.','MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',.1,'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',.1);
% plot(polyshape(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100))));

% 3D:
EP = Ellipsoid(eye(3),[0;0;0]);
PZ = PolynomialZonotope(EP);

P = RandomPoints(PZ,100000);

figure; hold on;
scatter3(P(1,:),P(2,:),P(3,:),'Marker','.','MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',.1,'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',.1);


% Gc(1,2:2:2*N+2) = (-1).^(0:N).*pi.^(2*(0:N))./factorial(2*(0:N));
% Gs(2,3:2:2*N+3) = (-1).^(0:N).*pi.^(2*(0:N)+1)./factorial(2*(0:N)+1);
% Ec = 0:(2*N+1);
% Es = 0:(2*N+1);
% 
% E = kron(Ec)