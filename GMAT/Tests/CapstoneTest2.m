
%% Initialize:
clear all;
close all;

%% Directories:
GMATDirectory = 'C:\Users\ethan\OneDrive - Stanford\Documents\GMAT';

%% Add Path to GMAT
addpath([GMATDirectory '\bin']);

%% Load GMAT
load_gmat();
gmat.gmat.Clear();

%% Solar System User-Modified Values:
GMATAPI.LoadScript([cd '\..\Scripts\Ephemeris\DE440.script']);

%% Coordinate System:

EarthMoonBary = GMATAPI.Construct("Barycenter","EarthMoonBary");
EarthMoonBary.SetField("OrbitColor","Gold");
EarthMoonBary.SetField("OrbitColor","DarkGray");
EarthMoonBary.SetField("BodyNames","{Earth, Luna}");

%% Spacecraft:

Capstone = GMATAPI.Construct("Spacecraft","Capstone");
Capstone.SetField("DateFormat","TDBGregorian");
Capstone.SetField("Epoch","19 Nov 2022 00:00:00.000");
Capstone.SetField("CoordinateSystem","EarthMJ2000Ec");
Capstone.SetField("DisplayStateType","Cartesian");
Capstone.SetField("X",-3.998540475751689E+05);
Capstone.SetField("Y",6.126049756976774E+04);
Capstone.SetField("Z",-2.695661234835022E+04);
Capstone.SetField("VX",-1.767694283860905E-02);
Capstone.SetField("VY",-9.729575033742976E-01);
Capstone.SetField("VZ",1.417589603935756E-01);
Capstone.SetField("DryMass",100);
Capstone.SetField("Cd",2.2);
Capstone.SetField("Cr",1.5);
Capstone.SetField("DragArea",15);
Capstone.SetField("SRPArea",1);
Capstone.SetField("SPADDragScaleFactor",1);
Capstone.SetField("SPADSRPScaleFactor",1);
Capstone.SetField("AtmosDensityScaleFactor",1);
Capstone.SetField("ExtendedMassPropertiesModel","None");
Capstone.SetField("NAIFId",-10011001);
Capstone.SetField("NAIFIdReferenceFrame",-9011001);
Capstone.SetField("OrbitColor","Red");
Capstone.SetField("TargetColor","Teal");
% P = gmat.Rmatrix(6,6);
% P.SetElement(0,0,1e70); P.SetElement(1,1,1e70); P.SetElement(2,2,1e70); P.SetElement(3,3,1e70);  P.SetElement(4,4,1e70);P.SetElement(5,5,1e70);
% Capstone.SetField("OrbitErrorCovariance",gmat.Rmatrix(6,6));
Capstone.SetField("Attitude","CoordinateSystemFixed");
Capstone.SetField("SPADSRPInterpolationMethod","Bilinear");
Capstone.SetField("SPADSRPScaleFactorSigma",1e70);
Capstone.SetField("SPADDragInterpolationMethod","Bilinear");
Capstone.SetField("SPADDragScaleFactorSigma",1e70);
Capstone.SetField("AtmosDensityScaleFactorSigma",1e70);
Capstone.SetField("ModelFile","aura.3ds");
Capstone.SetField("ModelOffsetX",0);
Capstone.SetField("ModelOffsetY",0);
Capstone.SetField("ModelOffsetZ",0);
Capstone.SetField("ModelRotationX",0);
Capstone.SetField("ModelRotationY",0);
Capstone.SetField("ModelRotationZ",0);
Capstone.SetField("ModelScale",1);
Capstone.SetField("AttitudeDisplayStateType","Quaternion");
Capstone.SetField("AttitudeRateDisplayStateType","AngularVelocity");
Capstone.SetField("AttitudeCoordinateSystem","EarthMJ2000Eq");
Capstone.SetField("EulerAngleSequence","321");

%% Force Models:

% Gravity Fields:
MoonGravityField = GMATAPI.Construct("GravityField");
MoonGravityField.SetField("BodyName","Luna");
MoonGravityField.SetField("Degree",50);
MoonGravityField.SetField("Order",50);
MoonGravityField.SetField("StmLimit",100);
MoonGravityField.SetField("PotentialFile",[GMATDirectory '\data\gravity\luna\LP165P.cof']);
MoonGravityField.SetField("TideModel","None");

% Point Masses:
EarthGravity = GMATAPI.Construct("PointMassForce");
EarthGravity.SetField("BodyName","Earth");
MoonGravity = GMATAPI.Construct("PointMassForce");
MoonGravity.SetField("BodyName","Luna");
SunGravity = GMATAPI.Construct("PointMassForce");
SunGravity.SetField("BodyName","Sun");
MercuryGravity = GMATAPI.Construct("PointMassForce");
MercuryGravity.SetField("BodyName","Mercury");
VenusGravity = GMATAPI.Construct("PointMassForce");
VenusGravity.SetField("BodyName","Venus");
MarsGravity = GMATAPI.Construct("PointMassForce");
MarsGravity.SetField("BodyName","Mars");
JupiterGravity = GMATAPI.Construct("PointMassForce");
JupiterGravity.SetField("BodyName","Jupiter");
SaturnGravity = GMATAPI.Construct("PointMassForce");
SaturnGravity.SetField("BodyName","Saturn");
UranusGravity = GMATAPI.Construct("PointMassForce");
UranusGravity.SetField("BodyName","Uranus");
NeptuneGravity = GMATAPI.Construct("PointMassForce");
NeptuneGravity.SetField("BodyName","Neptune");
PlutoGravity = GMATAPI.Construct("PointMassForce");
PlutoGravity.SetField("BodyName","Pluto");

% Solar Radiation Pressure:
SRP = GMATAPI.Construct("SolarRadiationPressure");
SRP.SetField("Flux",1367);
SRP.SetField("SRPModel","Spherical");
SRP.SetField("Nominal_Sun",149597870.691);

% Relativity:
Relativity = GMATAPI.Construct("RelativisticCorrection");

% Earth Moon J2 with SRP Force Model:
EarthMoonJ2_SRP = GMATAPI.Construct("ForceModel","EarthMoonJ2_SRP");
EarthMoonJ2_SRP.SetField("CentralBody","Earth");
%EarthMoonJ2_SRP.SetField("PrimaryBodies","{Luna}");
%EarthMoonJ2_SRP.SetField("PointMasses","{Sun, Venus, Mars, Jupiter, Saturn, Uranus, Neptune}");
%EarthMoonJ2_SRP.SetField("Drag","None");
%EarthMoonJ2_SRP.SetField("SRP","On");
% EarthMoonJ2_SRP.SetField("RelativisticCorrection",true);
% EarthMoonJ2_SRP.SetField("ErrorControl","RSSStep");
% EarthMoonJ2_SRP.AddForce(MoonGravityField);
EarthMoonJ2_SRP.AddForce(SunGravity);
EarthMoonJ2_SRP.AddForce(VenusGravity);
EarthMoonJ2_SRP.AddForce(MarsGravity);
EarthMoonJ2_SRP.AddForce(JupiterGravity);
EarthMoonJ2_SRP.AddForce(SaturnGravity);
EarthMoonJ2_SRP.AddForce(UranusGravity);
EarthMoonJ2_SRP.AddForce(NeptuneGravity);
EarthMoonJ2_SRP.AddForce(SRP);
EarthMoonJ2_SRP.AddForce(Relativity);

% Moon J2:
MoonJ2 = GMATAPI.Construct("ForceModel","MoonJ2");
MoonJ2.SetField("CentralBody","Moon");
MoonJ2.AddForce(SunGravity);
MoonJ2.AddForce(EarthGravity);
MoonJ2.AddForce(MercuryGravity);
MoonJ2.AddForce(VenusGravity);
MoonJ2.AddForce(MarsGravity);
MoonJ2.AddForce(JupiterGravity);
MoonJ2.AddForce(SaturnGravity);
MoonJ2.AddForce(UranusGravity);
MoonJ2.AddForce(NeptuneGravity);
MoonJ2.AddForce(PlutoGravity);
MoonJ2.AddForce(MoonGravityField);
MoonJ2.AddForce(Relativity);

%% Propagators:

Prop = GMATAPI.Construct("Propagator","Prop");
Prop.SetField("Type","RungeKutta89");
Prop.SetField("FM","MoonJ2");
Prop.SetField("InitialStepSize",60);
Prop.SetField("Accuracy",1e-14);
Prop.SetField("MinStep",60);
Prop.SetField("MaxStep",60);
%Prop.SetField("StopIfAccuracyIsViolated",true);

%% Propagation:

% State Manager:
PSM = gmat.PropagationStateManager();
PSM.SetObject(Capstone);
PSM.BuildState();
MoonJ2.SetPropStateManager(PSM);
MoonJ2.SetState(PSM.GetState());

% Test The Model:
GMATAPI.Initialize();

% Setup Direct Access to Force Model:
MoonJ2.BuildModelFromMap();
MoonJ2.UpdateInitialData();

% Test Derivative Data:
MoonJ2.GetDerivatives(Capstone.GetState().GetState());

% Create Propagator:
Prop.AddPropObject(Capstone);
Prop.PrepareInternals();
Propagator = Prop.GetPropagator();

% Propagate:
t = 0:60:33*24*60*60; % 30 Days
x = zeros(6,length(t));
xr = zeros(6,length(t));
x(:,1) = Capstone.GetState().GetState();
xr(:,1) = Propagator.GetState();
for i = 1:length(t)-1
    Propagator.Step(t(i+1)-t(i));
    Propagator.UpdateSpaceObject();
    xr(:,i+1) = Propagator.GetState();
    x(:,i+1) = Capstone.GetState().GetState();
end
