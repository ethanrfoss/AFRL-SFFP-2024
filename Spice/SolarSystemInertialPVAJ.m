
function [PVAJ] = SolarSystemInertialPVAJ(Body,Time)

P1 = cspice_spkezr(upper(Body), Time-1, 'J2000', 'NONE', 'SSB');
P2 = cspice_spkezr(upper(Body), Time, 'J2000', 'NONE', 'SSB');
P3 = cspice_spkezr(upper(Body), Time+1, 'J2000', 'NONE', 'SSB');
PVAJ = [P2;1/2*(P3(4:6)-P1(4:6));P3(4:6)-2*P2(4:6)+P1(4:6)]; % Position, Velocity, Acceleration, Jerk

end