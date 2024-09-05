
function LoadKernels

Directory = fileparts(mfilename('fullpath'));

% Add Paths:
addpath([Directory '\mice\lib']);
addpath([Directory '\mice\src\mice']);

% Load Kernels:
cspice_furnsh([Directory '\Kernels\pck00010.tpc']);
cspice_furnsh([Directory '\Kernels\de-403-masses.tpc']);
cspice_furnsh([Directory '\Kernels\naif0012.tls']);
cspice_furnsh([Directory '\Kernels\de440s.bsp']);

end