
function f = InertialOrbitFigure(System,Orbit,Options)

f = InertialFigure(System);
if Options.Planar
    p = plot(Orbit.InertialState(1,:),Orbit.InertialState(2,:),'Color','k','LineWidth',2); uistack(p,'bottom');
else
    p = plot3(Orbit.InertialState(1,:),Orbit.InertialState(2,:),Orbit.InertialState(3,:),'Color','k','LineWidth',2); uistack(p,'bottom');
end
axis tight;

end