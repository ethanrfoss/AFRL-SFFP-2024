function TPBVP = PseudoArclengthBC(DX,dS,x0prev,Tprev)

nx = length(x0prev);
TPBVP.g = @(x0,xf,p,T) [x0-x0prev;T-Tprev]'*DX-dS;
TPBVP.H0 = @(x0,xf,p,T) DX(1:nx)';
TPBVP.Hf = @(x0,xf,p,T) zeros(1,nx);
TPBVP.I = @(x0,xf,p,T) zeros(1,0);
TPBVP.Jf = @(x0,xf,p,T) DX(end);

end