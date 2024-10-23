function TabFig = TabFigures(Figures, Titles, TabFig)
    
    % Create a new figure for the tab group if one not specified
    if nargin < 3
        TabFig = figure;
        TabGroup = uitabgroup('Parent', TabFig);
    else
        TabGroup = TabFig.Children;
    end
   
    % Iterate through each figure handle and create tabs
    for i = 1:length(Figures)
        % Create a tab for each figure
        Tab = uitab(TabGroup, 'Title', Titles(i));
        
        % Copy the contents of the figure into the tab
        copyobj([findobj(Figures(i).Children,'type','axes'), findobj(Figures(i).Children,'type','colorbar'), findobj(Figures(i).Children,'type','legend'), findobj(Figures(i).Children,'type','UIControl')], Tab);
        if ~isempty(findobj(Figures(i).Children,'type','colorbar'))
            colormap(findobj(Tab.Children,'type','axes'),colormap(Figures(i)));
        end

        % Close figure
        close(Figures(i));
    end
    
    % Make the tab group visible
    TabGroup.Visible = 'on';
    drawnow;

end
