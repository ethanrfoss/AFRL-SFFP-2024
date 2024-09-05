
function x = DirectionalTensorPropagate(x0,Phi,R)

% Get Initial Displacement:
dx0 = x0 - Phi{1,1};

% Use STTs to Propagate:
x = [Phi{1,:}];
for i = 1:size(Phi,2)
    if size(Phi,1) > 1
        x(:,i) = x(:,i) + Phi{2,i}*dx0;
    end
    for j = 3:size(Phi,1)
        x(:,i) = x(:,i) + 1/factorial(j-1)*TensorVectorMultiplication(Phi{j,i},R(1:size(Phi{j,i},ndims(Phi{j,i})),:)*dx0);
    end
end

end