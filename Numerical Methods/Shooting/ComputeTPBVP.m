function [TPBVP] = ComputeTPBVP(System,TPBVP)

[x,p] = TPBVP.Variables(System);
x = sym(x);
p = sym(p);
syms t;
for i = 1:length(x)
    x0(i,:) = sym([char(x(i)) '0']);
    xf(i,:) = sym([char(x(i)) 'f']);
end
TPBVP.nx = length(x);
TPBVP.np = length(p);

% Get Dynamics and BCs
f = TPBVP.Dynamics(System,x);
g = TPBVP.BoundaryConditions(System,x0,xf,p,t);

% Linearize
TPBVP.A = matlabFunction(jacobian(f,x),"Vars",{x});
TPBVP.H0 = matlabFunction(jacobian(g,x0),"Vars",{x0,xf,p,t});
TPBVP.Hf = matlabFunction(jacobian(g,xf),"Vars",{x0,xf,p,t});
TPBVP.I = matlabFunction(jacobian(g,p),"Vars",{x0,xf,p,t});
TPBVP.Jf = matlabFunction(jacobian(g,t),"Vars",{x0,xf,p,t});
TPBVP.f = matlabFunction(f,"Vars",{x});
TPBVP.g = matlabFunction(g,"Vars",{x0,xf,p,t});

end