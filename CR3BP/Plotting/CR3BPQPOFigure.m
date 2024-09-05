
function f = CR3BPQPOFigure(System,QPO,Options)

% Create Figure:
Options.JacobiConstant = QPO.JacobiConstant;
Options.ZeroVelocityCurves = true;
Options.LagrangePointLabels = true;
if ~isfield(Options,'Axis') 
    X = cellfun(@(x) getfield(x,'x'),QPO.Trajectory{1},'UniformOutput',false); X = [X{:}];
    Options.Axis = [min(X(1,:)) max(X(1,:)); min(X(2,:)) max(X(2,:)); min(X(3,:)) max(X(3,:))]';
    Options.Axis = reshape((Options.Axis-mean(Options.Axis))*1.5+mean(Options.Axis),[1 6]);
end
f = CR3BPFigure(System,Options);

% Plot Orbits:
if Options.Planar
    for i = 1:size(QPO.u,3)
        p = plot(QPO.Trajectory{1}{i}.x(1,:),QPO.Trajectory{1}{i}.x(2,:),'Color',[0 0 0 .5]); uistack(p,'bottom');
    end
else
    for i = 1:size(QPO.u,3)
        p = plot3(QPO.Trajectory{1}{i}.x(1,:),QPO.Trajectory{1}{i}.x(2,:),QPO.Trajectory{1}{i}.x(3,:),'Color',[0 0 0 .5]); uistack(p,'bottom');
    end
end
axis tight;

end