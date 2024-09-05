
clear all; close all;

mu = 1; 
ep = 10^-6;
syms x1 x2;
x = [x1;x2];
f = [x2;mu*(1-x1^2)*x2-x1];
A = jacobian(f,x);
B = [0;.1];
syms lambda1 lambda2;
lambda = [lambda1;lambda2];
z = [x;lambda];
g = [f-B*B'*lambda/(sqrt(lambda'*B*B'*lambda+10^-6));-A'*lambda];
Ag = jacobian(g,z);
g = matlabFunction(g,"Vars",{z});
Ag = matlabFunction(Ag,"Vars",{z});
f = matlabFunction(f,"Vars",{x});
clear x lambda z;

N = 1000;
T = 10;
dt = .1;
TT = 0:dt:T;
x0 = [1;1];

lambda0 = [1*cos(linspace(0,2*pi,N));1*sin(linspace(0,2*pi,N))];

xf = zeros(2,N);
lambdaf = zeros(2,N);
for i = 1:N
    z0 = [x0;lambda0(:,i);0];
    [t{i},z] = ode89(@(t,z)[g(z(1:4));z(3:4)'*f(z(1:2))-norm(B'*z(3:4))],TT,z0);
    xf(:,i) = z(end,1:2)';
    x{i} = z(:,1:2)';
    u{i} = B'*z(:,3:4)';
    lambdaf(:,i) = z(end,3:4)';
    V{i} = z(:,5)';
end

% figure;
% color = jet(N);
% subplot(2,2,1); hold on;
% for i = 1:N-1
%     plot(lambda0(1,i:i+1),lambda0(2,i:i+1),'Color',color(i,:),'LineWidth',3);
% end
% subplot(2,2,2); hold on;
% for i = 1:N-1
%     plot(x{i}(1,:),x{i}(2,:),'Color',[0 0 0 .2]);
%     plot(xf(1,i:i+1),xf(2,i:i+1),'Color',color(i,:),'LineWidth',3);
% end
% subplot(2,2,[3,4]); hold on;
% for i = 1:N-1
%     plot(t{i},u{i},'Color',color(i,:),'LineWidth',1);
% end

%% Propagate an ellipsoidal set:
% P0 = [zeros(2,2) zeros(2,2);zeros(2,2) 1*eye(2,2)];
% z0 = [x0;0;1];
% Ad = [0 1;-1 0];
% % J = [Ad zeros(2,1) -B;zeros(2,2) -Ad'];
% [~,zp] = ode89(@(t,zp) [g(zp(1:4));reshape(Ag(zp(1:4))*reshape(zp(5:20),[4,4]),[16,1])],TT,[z0;reshape(eye(4,4),[16,1])]);
% Phi = reshape(zp(:,5:20)',[4,4,length(TT)]);
% z = zp(:,1:4)';
% P = zeros(4,4,length(TT));
% for i = 1:length(TT)
%     P(:,:,i) = Phi(:,:,i)*P0*Phi(:,:,i)';
% end

%% Ellipsoid:
function E = Ellipsoid(x,P)

[V,D] = eig(P);

E = V*sqrt(D)*[cos(linspace(0,2*pi,100));sin(linspace(0,2*pi,100))]+x;

end

%%

% Create a new figure
f = uifigure('Name', 'Van der Pol Reachability', 'NumberTitle', 'off');

% Create a UI axes for the plot
a = uiaxes('Parent', f, 'Position', [50, 100, 400, 300]); hold(a,'on');

% Axis:
X = [x{:}];
set(a,'XLim',[min(X(1,:)) max(X(1,:))]);
set(a,'YLim',[min(X(2,:)) max(X(2,:))]);

% Initial data for the plot
color = jet(N);
for i = 1:N-1
    p1(i) = plot(a,x{i}(1,:),x{i}(2,:),'Color',[0 0 0 .2]);
    p2(i) = plot(a,xf(1,i:i+1),xf(2,i:i+1),'Color',color(i,:),'LineWidth',3);
end
% E = Ellipsoid(z(1:2,end),P(1:2,1:2,end));
% p3 = plot(a,z(1,end),z(2,end),'bx');
% p4 = plot(a,E(1,:),E(2,:),'Color','b','LineWidth',2);

% Add a slider control
sld = uislider(f, 'Orientation', 'vertical', ...
               'ValueChangedFcn', @(i1,i2)slider(p1,p2,TT,x,i1,i2), ...
               'Limits', [0, T], ...
               'MajorTicks', TT, ...
               'Position', [500, 100, 3, 300]);

% Callback function for the slider
function slider(p1,p2,TT,x,hObject, ~)
    % Get the current value of the slider
    t = hObject.Value;

    ind = find(min(abs(t-TT))==abs(t-TT));

    % Update the plot based on the slider value
    for i = 1:length(x)-1
        set(p1(i),'XData',x{i}(1,1:ind),'YData',x{i}(2,1:ind));
        set(p2(i),'XData',[x{i}(1,ind) x{i+1}(1,ind)],'YData',[x{i}(2,ind) x{i+1}(2,ind)]);
    end
    % E = Ellipsoid(z(1:2,ind),P(1:2,1:2,ind));
    % set(p3,'XData',z(1,ind),'YData',z(2,ind));
    % set(p4,'XData',E(1,:),'YData',E(2,:));

end
