%% Earth Moon System Parameters
function System = CislunarSystem(System)
% CislunarSystem - Values for the Cislunar System
% 
% Inputs:
%    System (optional) - Structure Containing the System variables
% 
% Outputs:
%    System - Structure Containing the System variables
%
% Author: Ethan Foss (erfoss@stanford.edu)
% Date: 06/28/2024

% System Parameters ( Earth-Moon ):
System.PrimaryBody.Mass = 5.9722*10^24; % Earth Mass (kg)
System.SecondaryBody.Mass = 7.34767309*10^22; % Moon Mass (kg)
System.PrimaryBody.Radius = 6370*1000; % Earth Radius (m)
System.SecondaryBody.Radius = 1740*1000; % Moon Radius (m)
System.OrbitalDistance = 382500*1000; % Earth Moon Distance (m)
System.OrbitalRate = 1/sqrt(System.OrbitalDistance^3/(GravitationalConstant*(System.PrimaryBody.Mass+System.SecondaryBody.Mass))); % Earth Moon Rotation Rate (rad/s)
System.OrbitalEccentricity = .0549; % Eccentricity of Earth-Moon System

% Canonicalization:
System.TU = 1/System.OrbitalRate; % Time Unit (s)
System.DU = System.OrbitalDistance; % Distance Unit (m)
System.MU = System.PrimaryBody.Mass+System.SecondaryBody.Mass; % Mass Unit (kg)

end