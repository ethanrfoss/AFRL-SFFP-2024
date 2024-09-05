
function State = SolarSystemInertialState(Body,Time)

[State] = cspice_spkezr(upper(Body), Time, 'J2000', 'NONE', 'SSB');

end