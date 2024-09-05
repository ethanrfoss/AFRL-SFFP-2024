function [Figures] = OrbitalElementsFigure(Orbit,Elements,Options)

Names = ["Semi-Major Axis","Eccentricity","Inclination","Right Ascension of the Ascending Node","Argument of Periapsis","True Anomaly"];

for i = 1:length(Elements)
    Figures(i) = figure; hold on; box on; grid on;

    if isfield(Options,'HighlightIndices')
        for j = 1:length(Options.HighlightIndices)
            p = plot(Orbit(Options.HighlightIndices(j)).t,Orbit(Options.HighlightIndices(j)).Elements(Elements(i),:),'Color','k','LineWidth',2); uistack(p,'bottom');
        end
    end
    for j = 1:length(Orbit)
        if Orbit(i).Converged
            color = interp1(linspace(min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant]),length(Options.Colormap)),Options.Colormap,Orbit(j).JacobiConstant);
            p = plot(Orbit(j).t,Orbit(j).Elements(Elements(i),:),'Color',color,'LineWidth',1); uistack(p,'bottom');
        end
    end
    xlabel('Time');
    ylabel(['Osculating ' Names(i)]);
    if Options.Colorbar
        colormap(Options.Colormap);
        cbar = colorbar; ylabel(cbar,'Jacobi Constant');
        clim([min([Orbit(:).JacobiConstant]),max([Orbit(:).JacobiConstant])]);
    end
    axis tight;
end

end