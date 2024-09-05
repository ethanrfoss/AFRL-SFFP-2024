%% Initialize CR3BP System
function System = InitializeCR3BP(System)
% InitializeCR3BP - Computes CR3BP Parameters from System Data
% 
% Inputs:
%    System - Structure Containing the System variables
% 
% Outputs:
%    System - Structure Containing the System variables with CR3BP
%    Parameters
%
% Example Usage:
%    System = CislunarSystem;
%    System = InitializeCR3BP(System);
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

% Compute Gravitational Parameter:
System.mu = (System.PrimaryBody.Mass/System.SecondaryBody.Mass+1)^(-1); % Non-Dimensional Gravity Constant

% Construct Non-Dimensionalized Parameters:
System.m1 = System.PrimaryBody.Mass/(System.MU);
System.m2 = System.SecondaryBody.Mass/(System.MU);
System.D = System.OrbitalDistance/System.DU; % Should be 1
System.R1 = System.PrimaryBody.Radius/System.DU;
System.R2 = System.SecondaryBody.Radius/System.DU;

% Spacecraft Normalization:
if isfield(System,"Spacecraft")
    System.Spacecraft.Thrust = System.Spacecraft.Thrust*System.TU^2/(System.MU*System.DU);
    System.Spacecraft.WetMass = System.Spacecraft.WetMass/System.MU;
    System.Spacecraft.ExhaustVelocity = System.Spacecraft.ExhaustVelocity*System.TU/System.DU;
end

% Lagrange Points of Earth Moon System:
System.LagrangePoints = LagrangePoints(System);

end