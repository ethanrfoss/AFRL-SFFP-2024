
function Orbit = MMRFamily(System,MMR,Options)

% Options:
DefaultOptions.Plot = true; % Plotting
DefaultOptions.PlotInertial = true;
DefaultOptions.Save = false; 
DefaultOptions.SaveDirectory = fileparts(mfilename('fullpath'));
DefaultOptions.Colorbar = true;
DefaultOptions.Colormap = cool;
if nargin < 3
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.Planar = true;

% Initial Guess for MMR:
a = (MMR.q^2/MMR.p^2)^(1/3); % Semimajor Axis
if Options.Apoapsis
    r = -a*(1+MMR.e); % Argument of apoapsis
    V = -sqrt(2*(1-System.mu)*(1/-r-1/(2*a)));
else
    r = a*(1-MMR.e);
    V = sqrt(2*(1-System.mu)*(1/r-1/(2*a)));
end
if ~isfield(MMR,'x')
    InitialGuess.x = CR3BPFrameTransform(System,0,[r;0;0;V],["Earth","Barycenter"],["Inertial","Rotating"]);
else
    InitialGuess.x = MMR.x;
end 
if ~isfield(MMR,'P')
    InitialGuess.P = 2*pi*MMR.q;
else
    InitialGuess.P = MMR.P;
end

% Solve For First Orbit Using fsolve:
OrbitOptions.Shooting.AbsoluteTolerance = 10^-15;
OrbitOptions.Shooting.RelativeTolerance = 10^-12;
OrbitOptions.Shooting.Tolerance = 10^-10;
OrbitOptions.Shooting.MaximumIterations = 500;
OrbitOptions.Shooting.Printing = false;
OrbitOptions.Shooting.RootFinder = 'fsolve';
OrbitOptions.Save = false;
OrbitOptions.Plot = Options.Plot;
OrbitOptions.PlotInertial = Options.PlotInertial;
OrbitOptions.N = Options.N;
Orbit = SymmetricOrbit(System,InitialGuess,OrbitOptions);
OrbitOptions = rmfield(OrbitOptions,'Shooting');

% Save Data
if Options.Save
    Options.SaveName = ['MMR' num2str(MMR.p) 'to' num2str(MMR.q)];
    SaveData(struct('Orbit',Orbit),Options);
    f = gcf;
    SaveFigure(f,Options);
end

% Compute Families:
OrbitOptions.Nm = Options.Nm;
OrbitOptions.Np = Options.Np;
Orbit = SymmetricOrbitFamily(System,Orbit,OrbitOptions);

% Compute Extra MMR Data:
TBP.mu = 1-System.mu;
for i = 1:length(Orbit)
    Orbit(i).Elements = OrbitalElements(TBP,Orbit(i).InertialState,"Keplerian");
    Orbit(i).OsculatingPeriod = OrbitalElements(TBP,Orbit(i).InertialState,"P");
    disp(['Computing Elements for Orbit ' num2str(i)]);
end

% Plot Elements and Stability:
if Options.Plot
    f = gcf;
    f = TabFigures(DoubleAxisFigure([Orbit(:).JacobiConstant],{[Orbit(:).StabilityIndex],[Orbit(:).Period]},"Jacobi Constant",["Stability Index","Period"]), "Period and Stability",f);
    Options.HighlightIndices = Options.Nm+1;
    f = TabFigures(OrbitalElementsFigure(Orbit,[1,2],Options), ["Semimajor Axis","Eccentricity"],f);
end

% Save Data
if Options.Save
    Options.SaveName = ['MMR' num2str(MMR.p) 'to' num2str(MMR.q) 'Family'];
    SaveData(struct('Orbit',Orbit),Options);
    if Options.Plot
        SaveFigure(f,Options);
    end
end

end