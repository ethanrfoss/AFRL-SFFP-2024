
%% Derivative Tests:
% Test different computation strategies for higher order derivative
% tensors
System = CislunarSystem;
System = InitializeCR3BP(System);
Options.nx = 6;
f = @(t,x) EOMCR3BP(System,t,x);

%% Matlab Symbolic Analytic Derivatives:
for i = 1:5
    Options.N = i;
    tic;
    D = DerivativeTensors(f,Options);
    T= toc;
    disp(['Computation Time for Matlab Symbolic Derivatives of Order ' num2str(i) ': ' num2str(T)]);
end

N = 1000;
x = rand(6,N);
for i = 1:5
    tic;
    for j = 1:N
        D{i}(0,x(:,i));
    end
    T = toc;
    disp(['Average Evaluation Time for Matlab Symbolic Derivatives of Order ' num2str(i) ': ' num2str(T/N)]);
end

%% Casadi Symbolic Analytic Derivatives:
for i = 1:5
    Options.N = i;
    tic;
    D = AutomaticDerivativeTensors(f,Options);
    T= toc;
    disp(['Computation Time for Casadi Automatic Derivatives of Order ' num2str(i) ': ' num2str(T)]);
end

for i = 1:5
    tic;
    for j = 1:N
        full(D{i}(0,x(:,i)));
    end
    T = toc;
    disp(['Average Evaluation Time for Casadi Automatic Derivatives of Order ' num2str(i) ': ' num2str(T/N)]);
end

%% Casadi Numerical Derivatives:

