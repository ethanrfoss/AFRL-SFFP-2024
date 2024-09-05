%% Deep Space 1 Spacecraft System Parameters
% This function initializes Earth-Moon system parameters
function Spacecraft = DeepSpaceOne

% Mass, Thrust, Specific Impulse:
Spacecraft.Thrust = .092; % Maximum Thrust Capability [N]
Spacecraft.WetMass = 486.3; % Initial Mass [kg]
Spacecraft.ExhaustVelocity = 9.81*3200; % Exhaust Velocity [m/s]

end