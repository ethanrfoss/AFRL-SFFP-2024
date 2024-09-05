
function x = TensorPropagate(x0,Phi)

% Get Initial Displacement:
dx0 = x0 - Phi{1,1};

% Use STTs to Propagate:
x = [Phi{1,:}];
for i = 1:size(Phi,2)
    for j = 2:size(Phi,1)
        x(:,i) = x(:,i) + 1/factorial(j-1)*TensorVectorMultiplication(Phi{j,i},dx0);
    end
end

end