
% Ephemeris STT Test:
x0 = [.8;0;0;0;0;0];
Options.Ephemeris = true;
Options.TensorOrder = 3;
Options.N = 3;
f = EOMCR3BP(System,Options);

Phi = StateTransitionTensors2(f,0:5,x0,Options);