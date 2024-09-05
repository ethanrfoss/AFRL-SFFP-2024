
function f = CR3BPQPOFamilyFigure(System,QPO,Options)

% Create Figure:
if ~isfield(Options,'JacobiConstant')
    Options.ZeroVelocityCurves = true;
    Options.Colorbar = false;
end
Options.LagrangePointLabels = true;
if ~isfield(Options,'Axis') 
    X = cellfun(@(x) getfield(x,'x'),QPO(:).Trajectory{1},'UniformOutput',false); X = [X{:}];
    Options.Axis = [min(X(1,:)) max(X(1,:)); min(X(2,:)) max(X(2,:)); min(X(3,:)) max(X(3,:))]';
    Options.Axis = reshape((Options.Axis-mean(Options.Axis))*1.5+mean(Options.Axis),[1 6]);
end
f = CR3BPFigure(System,Options);

% Plot Orbits:
for i = 1:length(QPO)
    for j = 1:size(QPO(i).u,3)
        if Options.Colorbar
            color = interp1(linspace(min([QPO(:).JacobiConstant]),max([QPO(:).JacobiConstant]),length(Options.Colormap)),Options.Colormap,QPO(i).JacobiConstant);
        else
            color = [0 0 0 .5];
        end
        if Options.Planar
            p = plot(QPO(i).Trajectory{1}{j}.x(1,:),QPO(i).Trajectory{1}{j}.x(2,:),'Color',color); uistack(p,'bottom');
        else
            p = plot3(QPO(i).Trajectory{1}{j}.x(1,:),QPO(i).Trajectory{1}{j}.x(2,:),QPO(i).Trajectory{1}{j}.x(3,:),'Color',color); uistack(p,'bottom');
        end
    end
end
if Options.Colorbar
    colormap(Options.Colormap);
    cbar = colorbar; ylabel(cbar,'Jacobi Constant');
    clim([min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant])]);
end

axis tight;

end