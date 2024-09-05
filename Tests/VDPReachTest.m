
% Van Der Pol Test:
f = @(x,u,t) [x(2);1*(1-x(1)^2)*x(2)-x(1)+.05*u(1)];
T = 0:.1:3;
x0 = [.3;0];

f1 = @(z,t) [f(z(1:2),z(3:4),t);0;0];

% Extremal Trajectories:
Options.InputSet = 'UnitCube';
Options.N = 100;
Trajectory = Extremals(f,T,x0,Options);

% Reachable Sets:
u0 = [0;0];
Options.STTOrder = 3;
Options.MaxPZOrder = 200;
Options.ReductionOrder = 200;
Options.DirectionalSTT = false;
Options.EigenCutoff = 10^-2;
% [ReachSet,InputSet] = Reachability(f1,T,x0,u0,x0,PolynomialZonotope(Ellipsoid(eye(2,2),[0,0]),struct('Order',2,'Method','Squircle','ID','Unique')),Options);
[ReachSet,InputSet] = Reachability(f1,T,x0,u0,x0,PolynomialZonotope(eye(2,2),eye(2,2)),Options);

figure; hold on;
color = autumn(length(T));
for i = 2:length(T)
    p1 = plot(ReachSet(i),struct('Splits',10));
    i
end
for i = 1:length(Trajectory)
    p2 = plot(Trajectory(i).x(1,:),Trajectory(i).x(2,:),'Color',[0 0 0 .2]);
end
% for i = 1:length(T)
%     plot(arrayfun(@(t) t.x(1,i),Trajectory),arrayfun(@(t) t.x(2,i),Trajectory),'Color',color(i,:),'LineWidth',2);
% end
legend([p1,p2],'Reachable Sets','Extremal Trajectories');
xlabel('$x_1$'); ylabel('$x_2$');
grid on;

%% Test Sampling of Unit Circle:

% Get Trajectories:
% Options.InputSet = 'UnitSphere';
% Trajectory2 = Extremals(f,T,x0,Options);
% 
% % Sample reachable sets:
% N = 1000000;
% alpha = rand(ReachSet(end).dims.p,N)*2-1;
% % for i = 1:length(InputSet)
% %     [~,ind] = ismember(InputSet(i).id,ReachSet(end).id);
% %     r = sqrt(rand(1,N));
% %     t = 2*pi*rand(1,N);
% %     alpha(ind,:) = [r.*cos(t);r.*sin(t)];
% % end
% 
% G = ReachSet(end).G;
% E = ReachSet(end).E;
% P = zeros(2,N);
% for i = 1:N
%     P(:,i) = G*prod(alpha(:,i).^E,1)';
% end
% 
% figure; hold on
% plot(arrayfun(@(t) t.x(1,end),Trajectory2),arrayfun(@(t) t.x(2,end),Trajectory2),'r.');
% plot(arrayfun(@(t) t.x(1,end),Trajectory),arrayfun(@(t) t.x(2,end),Trajectory),'r.');
% plot(P(1,:),P(2,:),'b.');
% 
