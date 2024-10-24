clear all;
import casadi.*

System = CislunarSystem;
System = InitializeCR3BP(System);

Options.Planar = true;
Options.nx = 4;
Options.N = 3;
EOM = @(t,x) EOMCR3BP(System,t,x,Options);
A = cell(Options.N);

% Create Symbolic Variables in CasADi:
t = MX.sym('t',1);
x = MX.sym('x',[Options.nx,1]);
A{1} = EOM(t,x);

% Loop Through:
for i = 2:Options.N
    
    % Allocate:
    A{i}  = cell(repmat(Options.nx,1,i));

    % Compute gradients:
    NN = num2cell(dec2base(0:Options.nx^(i-1)-1,Options.nx)-'0'+1);
    for j = 1:size(NN,1)
        
        % Compute Gradient from previous gradients:
        G = gradient(A{i-1}{NN{j,:}},x);
        for k = 1:Options.nx
            A{i}{NN{j,:},k} = G(k);
        end

    end

end

for i = 1:Options.N
    NN = num2cell(dec2base(0:Options.nx^i-1,Options.nx)-'0'+1);
    D{i} = cell(size(A{i}));
    for j = 1:size(NN,1)
        D{i}{NN{j,:}} = Function(['A' num2str(i) strrep(num2str([NN{j,:}]),' ','')],{t,x},{A{i}{NN{j,:}}});
    end
end