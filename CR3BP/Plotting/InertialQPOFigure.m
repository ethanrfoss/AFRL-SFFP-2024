
function f = InertialQPOFigure(System,QPO,Options)

f = InertialFigure(System);
if Options.Planar
    for i = 1:size(QPO.u,3)
        p = plot(QPO.Trajectory{1}{i}.InertialState(1,:),QPO.Trajectory{1}{i}.InertialState(2,:),'Color','k','LineWidth',1); uistack(p,'bottom');
    end
else
    for i = 1:size(QPO.u,3)
        p = plot(QPO.Trajectory{1}{i}.InertialState(1,:),QPO.Trajectory{1}{i}.InertialState(2,:),QPO.Trajectory{1}{i}.InertialState(3,:),'Color','k','LineWidth',1); uistack(p,'bottom');
    end
end
axis tight;

end