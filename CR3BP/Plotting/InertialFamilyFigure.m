
function f = InertialFamilyFigure(System,Orbit,Options)

f = InertialFigure(System);
if Options.Planar
    if isfield(Options,'HighlightIndices')
        for i = 1:length(Options.HighlightIndices)
            p = plot(Orbit(Options.HighlightIndices(i)).InertialState(1,:),Orbit(Options.HighlightIndices(i)).InertialState(2,:),'Color','k','LineWidth',2); uistack(p,'bottom');
        end
    end
    for i = 1:length(Orbit)
        color = interp1(linspace(min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant]),length(Options.Colormap)),Options.Colormap,Orbit(i).JacobiConstant);
        p = plot(Orbit(i).InertialState(1,:),Orbit(i).InertialState(2,:),'Color',color,'LineWidth',1); uistack(p,'bottom');
    end
else
    if isfield(Options,'HighlightIndices')
        for i = 1:length(Options.HighlightIndices)
            p = plot3(Orbit(Options.HighlightIndices(i)).InertialState(1,:),Orbit(Options.HighlightIndices(i)).InertialState(2,:),Orbit(Options.HighlightIndices(i)).InertialState(3,:),'Color','k','LineWidth',2); uistack(p,'bottom');
        end
    end
    for i = 1:length(Orbit)
        color = interp1(linspace(min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant]),length(Options.Colormap)),Options.Colormap,Orbit(i).JacobiConstant);
        p = plot3(Orbit(i).InertialState(1,:),Orbit(i).InertialState(2,:),Orbit(i).InertialState(3,:),'Color',color,'LineWidth',1); uistack(p,'bottom');
    end
end
if Options.Colorbar
    colormap(Options.Colormap);
    cbar = colorbar; ylabel(cbar,'Jacobi Constant');
    clim([min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant])]);
end
axis tight;

end