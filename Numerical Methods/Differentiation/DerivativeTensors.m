
function [A] = DerivativeTensors(EOM,Options)

% Allocate Arrays:
A = cell(Options.N,1);

if ~isa(EOM,'sym')
    t = sym('t',1);
    x = sym('x',[Options.nx,1]);
    A{1} = EOM(x,t);
    Symbolic = false;
else
    x = Options.Variable;
    A{1} = EOM;
    Symbolic = true;
end

% Loop Through:
for i = 2:Options.N
    
    % Allocate:
    A{i}  = sym(zeros(repmat(Options.nx,1,i)));

    % Compute gradients:
    if Options.nx == 1
        NN = {1};
    else
        NN = num2cell(dec2base(0:Options.nx^(i-1)-1,Options.nx)-'0'+1);
    end
    for j = 1:size(NN,1)
        
        % Compute Gradient from previous gradients:
        A{i}(NN{j,:},:) = gradient(A{i-1}(NN{j,:}),x);

    end

end

if ~Symbolic
    for i = 1:Options.N
        A{i} = matlabFunction(A{i},"Vars",{t,x});
    end
end

end