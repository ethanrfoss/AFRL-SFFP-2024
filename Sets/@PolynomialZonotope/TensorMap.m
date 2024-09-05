
function PZ = TensorMap(PZ,T)

N = ndims(T);
G = T;
E = PZ.E;
for i = N:-1:2
    G = tensorprod(G,PZ.G,i,1);
end
% E = zeros(PZ.dims.p,PZ.dims.h^(N-1));
for i = 1:N-2
    E = kron(E,ones(1,PZ.dims.h))+kron(ones(1,PZ.dims.h^i),PZ.E);
end
% Eind = (dec2base(0:5^(3)-1,5)-'0')+1;
% for i = 1:PZ.dims.h^(N-1)
%     E(:,i) = sum(PZ.E(:,Eind(i,:)),2);
% end
G = reshape(G,[size(G,1),PZ.dims.h^(N-1)]);
PZ = PolynomialZonotope(G,E,PZ.id);

end