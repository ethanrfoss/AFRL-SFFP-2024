
function OE = OrbitalElements(System,x,Elements)

% Element Options:
Orbital.a = @SemimajorAxis;
Orbital.e = @Eccentricity;
Orbital.emag = @(S,x) norm(Eccentricity(S,x));
Orbital.i = @Inclination;
Orbital.Omega = @RAAN;
Orbital.omega = @ArgumentOfPeriapsis;
Orbital.theta = @TrueAnomaly;
Orbital.h = @AngularMomentum;
Orbital.hmag = @(S,x) norm(AngularMomentum(S,x));
Orbital.P = @OrbitalPeriod;
Orbital.Keplerian = @(S,x) [Orbital.a(S,x);Orbital.emag(S,x);Orbital.i(S,x);Orbital.Omega(S,x);Orbital.omega(S,x);Orbital.theta(S,x)];

% Planarity:
if size(x,1) == 4
    x = [x(1:2,:);zeros(1,size(x,2));x(3:4,:);zeros(1,size(x,2))];
end

% Get Index Lengths
ind = zeros(length(Elements),1);
for i = 1:length(Elements)
    ind(i) = length(Orbital.(Elements(i))(System,x(:,1)));
end

% Solve Orbital Elements
OE = zeros(sum(ind),size(x,2));
for i = 1:size(x,2)
    for j = 1:length(Elements)
        OE(sum(ind(1:j-1))+1:sum(ind(1:j)),i) = Orbital.(Elements(j))(System,x(:,i));
    end
end

end

%% Semi-Major Axis
function a = SemimajorAxis(System,x)
a = -System.mu/(2*(norm(x(4:6))^2/2-System.mu/norm(x(1:3))));
end

%% Angular Momentum
function h = AngularMomentum(System,x)
h = cross(x(1:3),x(4:6));
end

%% Eccentricity
function e = Eccentricity(System,x)
e = cross(x(4:6),AngularMomentum(System,x))/System.mu-x(1:3)/norm(x(1:3));
end

%% Inclination
function i = Inclination(System,x)
h = AngularMomentum(System,x);
i = acos(h(3)/norm(h));
end

%% Right Ascention of the Ascending Node
function Omega = RAAN(System,x)
w = cross([0;0;1],AngularMomentum(System,x));
if w(2) >= 0
    Omega = acos(w(1)/norm(w));
else
    Omega = 2*pi - acos(w(1)/norm(w));
end
end

%% Argument of Periapsis
function omega = ArgumentOfPeriapsis(System,x)
w = cross([0;0;1],AngularMomentum(System,x));
e = Eccentricity(System,x);
if e(3)>=0
    omega = acos(dot(e,w)/(norm(e)*norm(w)));
else
    omega = 2*pi-acos(dot(e,w)/(norm(e)*norm(w)));
end
end

%% True Anomaly
function theta = TrueAnomaly(System,x)
e = Eccentricity(System,x);
if dot(e,x(1:3)) >= 0
    theta = acos(dot(e,x(1:3))/(norm(e)*norm(x(1:3))));
else
    theta = 2*pi-acos(acos(dot(e,x(1:3))/(norm(e)*norm(x(1:3)))));
end
end

%% Orbital Period
function P = OrbitalPeriod(System,x)
P = 2*pi*sqrt(SemimajorAxis(System,x)^3/System.mu);
end