
function f = CR3BPOrbitFigure(System,Orbit,Options)

% Create Figure:
Options.JacobiConstant = Orbit.JacobiConstant;
Options.ZeroVelocityCurves = true;
Options.LagrangePointLabels = true;
if ~isfield(Options,'Axis') 
    Options.Axis = [min(Orbit.x(1,:)) max(Orbit.x(1,:)); min(Orbit.x(2,:)) max(Orbit.x(2,:)); min(Orbit.x(3,:)) max(Orbit.x(3,:))]';
    Options.Axis = reshape((Options.Axis-mean(Options.Axis))*1.5+mean(Options.Axis),[1 6]);
end
f = CR3BPFigure(System,Options);

% Plot Orbits:
if Options.Planar
    p = plot(Orbit.x(1,:),Orbit.x(2,:),'Color','k','LineWidth',2); uistack(p,'bottom');
else
    p = plot3(Orbit.x(1,:),Orbit.x(2,:),Orbit.x(3,:),'Color','k','LineWidth',2); uistack(p,'bottom');
end
% axis tight;
axis(Options.Axis);

end