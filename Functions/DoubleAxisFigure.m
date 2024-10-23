function f = DoubleAxisFigure(XData,YData,XLabel,YLabel)

% f = figure;
% 
% yyaxis left;
% plot(XData,YData{1}, '-b','LineWidth',2);
% xlabel(XLabel);
% ylabel(YLabel(1));
% 
% yyaxis right;
% plot(XData,YData{2},'-r','LineWidth',2);
% ylabel(YLabel(2));

f = figure;

ax1=axes('yaxislocation','left','YColor','b'); hold on;
ax2=axes('yaxislocation','right','XColor','none','YColor','r'); hold on;
set([ax1 ax2],'color','none');

plot(ax1,XData,YData{1}, '-b','LineWidth',2);
xlabel(ax1,XLabel);
ylabel(ax1,YLabel(1));

plot(ax2,XData,YData{2}, '-r','LineWidth',2);
ylabel(ax2,YLabel(2));

linkaxes([ax1 ax2],'x');

end