%% Lunar Icecube Spacecraft System Parameters
% This function initializes Earth-Moon system parameters
function Spacecraft = LunarIceCube

% Mass, Thrust, Specific Impulse:
Spacecraft.Thrust = 1.1*10^-3; % Maximum Thrust Capability [N]
Spacecraft.WetMass = 14; % Initial Mass [kg]
Spacecraft.ExhaustVelocity = 9.81*2800; % Exhaust Velocity [m/s]

end