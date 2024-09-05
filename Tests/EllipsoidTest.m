
mu = 0; 
syms x1 x2;
x = [x1;x2];
f = [2*x2;mu*(1-x1^2)*x2-x1];
A = jacobian(f,x);
f = matlabFunction(f,"Vars",{x});
A = matlabFunction(A,"Vars",{x});

T = 5;
dt = 1;
TT = 0:dt:T;

P0 = [.2 .1;.1 .1];
x0 = [1;1];

[V,D] = eig(P0);
E0 = V*sqrt(D)*[cos(linspace(0,2*pi,100));sin(linspace(0,2*pi,100))]+x0;

figure; hold on;
for i = 1:length(TT)
    Phi = expm(A(0)*TT(i));
    P = Phi*P0*Phi';
    [V,D] = eig(P);
    E = V*sqrt(D)*[cos(linspace(0,2*pi,100));sin(linspace(0,2*pi,100))]+Phi*x0;
    plot(E(1,:),E(2,:),'b','LineWidth',4);
    ET = Phi*E0;
    plot(ET(1,:),ET(2,:),'r','LineWidth',2);
end


