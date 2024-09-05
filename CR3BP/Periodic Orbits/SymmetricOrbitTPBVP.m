function TPBVP = SymmetricOrbitTPBVP(System,Options)

TPBVP.f = EOMCR3BP(System,Options);
TPBVP.A = JacobianCR3BP(System,Options);
TPBVP.N = Options.N;
if Options.Planar
    TPBVP.g = @(x0,xf,p,T) [x0(2);x0(3);xf(2);xf(3)];
    TPBVP.H0 = @(x0,xf,p,T) [0 1 0 0;0 0 1 0;0 0 0 0;0 0 0 0];
    TPBVP.Hf = @(x0,xf,p,T) [0 0 0 0;0 0 0 0;0 1 0 0;0 0 1 0];
    TPBVP.I = @(x0,xf,p,T) zeros(4,0);
    TPBVP.Jf = @(x0,xf,p,T) zeros(4,1);
    TPBVP.nx = 4;
    TPBVP.np = 0;
else
    TPBVP.g = @(x0,xf,p,T) [x0(2);x0(4);x0(6);xf(2);xf(4);xf(6)];
    TPBVP.H0 = @(x0,xf,p,T) [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0];
    TPBVP.Hf = @(x0,xf,p,T) [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1];
    TPBVP.I = @(x0,xf,p,T) zeros(6,0);
    TPBVP.Jf = @(x0,xf,p,T) zeros(6,1);
    TPBVP.nx = 6;
    TPBVP.np = 0;
end
if isfield(Options,'g')
    TPBVP.g = @(x0,xf,p,T) [TPBVP.g(x0,xf,p,T);Options.g(x0,xf,p,T)];
    TPBVP.H0 = @(x0,xf,p,T) [TPBVP.H0(x0,xf,p,T);Options.H0(x0,xf,p,T)];
    TPBVP.Hf = @(x0,xf,p,T) [TPBVP.Hf(x0,xf,p,T);Options.Hf(x0,xf,p,T)];
    TPBVP.I = @(x0,xf,p,T) [TPBVP.I(x0,xf,p,T);Options.I(x0,xf,p,T)];
    TPBVP.Jf = @(x0,xf,p,T) [TPBVP.Jf(x0,xf,p,T);Options.Jf(x0,xf,p,T)];
end

end