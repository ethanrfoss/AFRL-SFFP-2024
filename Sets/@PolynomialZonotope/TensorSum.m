
function PZ = TensorSum(PZ,Phi,c)

% Determine the Exponent Matrix:
if any(all(PZ.E==0,1))
    G = [sum(PZ.G(:,all(PZ.E==0,1)),2) PZ.G(:,~all(PZ.E==0,1))];
    E = [sum(PZ.E(:,all(PZ.E==0,1)),2) PZ.E(:,~all(PZ.E==0,1))];
    h = PZ.dims.h;
else
    G = [zeros(PZ.dims.n,1) PZ.G];
    E = [zeros(PZ.dims.p,1) PZ.E];
    h = PZ.dims.h + 1;
end
E0 = E;
G0 = G;
for i = 1:ndims(Phi{end})-2
    E = kron(E,ones(1,h))+kron(ones(1,h^i),E0);
end
[E,~,ind] = unique(E','rows','stable'); E = E';
G = zeros(size(Phi{1},1),size(E,2));

for k = 1:length(Phi)
    Gk = Phi{k};
    for i = ndims(Phi{k}):-1:2
        Gk = tensorprod(Gk,G0,i,1);
    end
    % Gk = reshape(Gk,[size(Gk,1),PZ.dims.h^(ndims(Phi{k})-1)]);
    Gk = Gk*1/factorial(ndims(Phi{k})-1);

    for i = 1:(h)^(ndims(Phi{k})-1)
        % G(:,ind(i)) = G(:,ind(i)) + Gk(:,ind2sub(repmat(PZ.dims.h,[ndims(Phi{k})-1,1]),i));
        G(:,ind(i)) = G(:,ind(i)) + Gk(:,i);
    end
    % G(:,1:(PZ.dims.h)^(ndims(Phi{k})-1)) = G(:,1:(PZ.dims.h)^(ndims(Phi{k})-1)) + Gk(:,ind1(1:(PZ.dims.h)^(ndims(Phi{k})-1)));
end
if nargin > 2
    G(:,1) = G(:,1) + c;
end
PZ = PolynomialZonotope(G,E,PZ.id);

end