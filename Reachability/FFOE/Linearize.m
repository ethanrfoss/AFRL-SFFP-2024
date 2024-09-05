function [LinearizedDynamics] = Linearize(System,Variables,Dynamics)

[x,u] = Variables(System);
[f] = Dynamics(System,x,u);

LinearizedDynamics.A = matlabFunction(Jacobian(f,x),"Vars",{x,u});
LinearizedDynamics.B = matlabFunction(Jacobian(f,u),"Vars",{x,u});
LinearizedDynamics.f = matlabFunction(f,"Vars",{x,u});

end

function fx = Jacobian(f,x)

if isempty(f)
    fx = sym(zeros(0,length(x)));
    return;
end

fx = [];

for i = 1:length(x)
    fx = [fx diff(f,x(i))];
end

end