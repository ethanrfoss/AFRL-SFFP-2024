function Saved = SaveFigure(Figure,Options)

if ~isfield(Options,'SaveName')
    Options.SaveName = input('Enter Save Name for Figure: ',"s");
end
if ~isfield(Options,'SaveDirectory')
    Options.SaveDirectory = cd;
end
if isfile([Options.SaveDirectory '\SavedImages\' strrep(Options.SaveName,' ','') '.fig'])
    if input('Figure Already Exists, Override? (y/n): ',"s") ~= 'y'
        disp('Discarding Figure'); 
        Saved = false;
        return;
    end
end
Saved = true;
savefig(Figure,[Options.SaveDirectory '\SavedImages\' strrep(Options.SaveName,' ','') '.fig']);
if any(arrayfun(@(c) isa(c, 'matlab.ui.container.TabGroup'), Figure.Children))
    Children = Figure.Children;
    TabGroup = findobj(Children,'Type','uitabgroup');
    if ~isempty(TabGroup)
        Tabs = TabGroup.Children;
        for i = 1:length(Tabs)
            exportgraphics(Tabs(i),[Options.SaveDirectory '\SavedImages\' strrep([Options.SaveName Tabs(i).Title],' ','')  '.png'],'Resolution',300);
        end
    end
else
    exportgraphics(Figure,[Options.SaveDirectory '\SavedImages\' strrep([Options.SaveName],' ','')  '.png'],'Resolution',300);
end

end