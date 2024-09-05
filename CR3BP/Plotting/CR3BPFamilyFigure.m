
function f = CR3BPFamilyFigure(System,Orbit,Options)

Options.LagrangePointLabels = true;
if ~isfield(Options,'Axis') 
    X = [Orbit(:).x];
    Options.Axis = [min(X(1,:)) max(X(1,:)); min(X(2,:)) max(X(2,:)); min(X(3,:)) max(X(3,:))]';
    Options.Axis = reshape((Options.Axis-mean(Options.Axis))*2+mean(Options.Axis),[1 6]);
end
f = CR3BPFigure(System,Options);
if Options.Planar
    if isfield(Options,'HighlightIndices')
        for i = 1:length(Options.HighlightIndices)
            p = plot(Orbit(Options.HighlightIndices(i)).x(1,:),Orbit(Options.HighlightIndices(i)).x(2,:),'Color','k','LineWidth',1); uistack(p,'bottom');
        end
    end
    for i = 1:length(Orbit)
        if Orbit(i).Converged
            color = interp1(linspace(min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant]),length(Options.Colormap)),Options.Colormap,Orbit(i).JacobiConstant);
            p = plot(Orbit(i).x(1,:),Orbit(i).x(2,:),'Color',color,'LineWidth',1); uistack(p,'bottom');
        end
    end
else
    if isfield(Options,'HighlightIndices')
        for i = 1:length(Options.HighlightIndices)
            p = plot3(Orbit(Options.HighlightIndices(i)).x(1,:),Orbit(Options.HighlightIndices(i)).x(2,:),Orbit(Options.HighlightIndices(i)).x(3,:),'Color','k','LineWidth',2); uistack(p,'bottom');
        end
    end
    for i = 1:length(Orbit)
        if Orbit(i).Converged
            color = interp1(linspace(min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant]),length(Options.Colormap)),Options.Colormap,Orbit(i).JacobiConstant);
            p = plot3(Orbit(i).x(1,:),Orbit(i).x(2,:),Orbit(i).x(3,:),'Color',color,'LineWidth',1);  uistack(p,'bottom');
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