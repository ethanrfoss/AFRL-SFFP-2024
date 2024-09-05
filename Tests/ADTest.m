clear all;
import casadi.*

System = CislunarSystem;
System = InitializeCR3BP(System);

Options.Planar = true;
Options.nx = 4;
Options.N = 3;
EOM = @(t,x) EOMCR3BP(System,t,x,Options);


% Create Symbolic Variables in CasADi:
t = MX.sym('t',1);
x = MX.sym('x',[Options.nx,1]);
A{1} = EOM(t,x);

% Loop Through:
for i = 2:Options.N
    
    % Allocate:
    n = repmat(Options.nx,1,i);
    % A{i} = MX.sym(['A' num2str(i)],repmat(Options.nx,1,i));
    A{i}  = MX.zeros(repmat(Options.nx,1,i));

    % Compute gradients:
    NN = num2cell(dec2base(0:Options.nx^(i-1)-1,Options.nx)-'0'+1);
    for j = 1:size(NN,1)
        
        % Compute Gradient from previous gradients:
        A{i}(NN{j,:},:) = gradient(A{1}(NN{j,:}),x);

    end

end

for i = 1:Options.N
    A{i} = Function(['A' num2str(i)],{t,x},{A{i}});
end