
function MMRGeneration

% Setup System:
System = CislunarSystem;
System = InitializeCR3BP(System);

% Set up Options for Plotting:
Options.PlotInertial = true;
Options.SaveDirectory = fileparts(mfilename('fullpath'));
Options.Save = true;

% 1:1 Family
% MMR.p = 1; MMR.q = 1;
% MMR.e = .85;
% Options.Nm = 400;
% Options.Np = 200;
% Options.N = 1;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 1:2 Family
% MMR.p = 1; MMR.q = 2;
% MMR.e = .85;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 1;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 1:3 Family
% MMR.p = 1; MMR.q = 3;
% MMR.e = .85;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 1;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 2:1 Family
% MMR.p = 2; MMR.q = 1;
% MMR.e = 0;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 3;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 2:3 Family
% MMR.p = 2; MMR.q = 3;
% MMR.e = 0.8;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 3;
% MMR.x = [.5488;0;0;1.139512];
% MMR.P = 18.497;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 2:5 Family
% MMR.p = 2; MMR.q = 5;
% MMR.e = .4;
% Options.Nm = 400;
% Options.Np = 90;
% Options.N = 1;
% Options.Apoapsis = true;
% MMRFamily(System,MMR,Options);

% 3:1 Family
% MMR.p = 3; MMR.q = 1;
% MMR.e = .9;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 5;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 3:2 Family
% MMR.p = 3; MMR.q = 2;
% MMR.e = .2;
% Options.Nm = 150;
% Options.Np = 250;
% Options.N = 1;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 3:4 Family
% MMR.p = 3; MMR.q = 4;
% MMR.e = 0.9;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 1;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 4:1 Family
% MMR.p = 4; MMR.q = 1;
% MMR.e = 0.6;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 3;
% Options.Apoapsis = false;
% MMRFamily(System,MMR,Options);

% 4:3 Family
MMR.p = 4; MMR.q = 3;
MMR.e = 0.8;
Options.Nm = 120;
Options.Np = 400;
Options.N = 1;
Options.Apoapsis = true;
MMRFamily(System,MMR,Options);

% 4:5 Family
% MMR.p = 4; MMR.q = 5;
% MMR.e = -.605;
% Options.Nm = 400;
% Options.Np = 400;
% Options.N = 1;
% Options.Apoapsis = true;
% MMRFamily(System,MMR,Options);

end