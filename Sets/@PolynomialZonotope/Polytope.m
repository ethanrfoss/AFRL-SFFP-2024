
function PT = Polytope(PZ)

ZeroInd = all(PZ.E==0,1);
H = any(PZ.E>1,1) & ~ZeroInd;
K = ~H & ~ZeroInd;
Z = Zonotope(PolynomialZonotope([PZ.G(:,ZeroInd) PZ.G(:,H)],[PZ.E(:,ZeroInd) PZ.E(:,H)],PZ.id));
G = [Z.c PZ.G(:,K) Z.G];
E = [zeros(PZ.dims.p+Z.dims.l,1) [PZ.E(:,K);zeros(Z.dims.l,sum(K))] [zeros(PZ.dims.p,Z.dims.l);eye(Z.dims.l)]];
alpha = (dec2bin(0:2^(size(E,1))-1)-'0')*2-1;
Vertices = zeros(size(alpha,1),PZ.dims.n);
for i = 1:size(E,1)
    Vertices(i,:) = (G*prod(alpha(i,:)'.^(E),1)')';
end
PT = Polytope(Vertices);

end