function PZ = Compact(PZ)

[E,~,ind] = unique(PZ.E','rows');
E = E';
G = zeros(PZ.dims.n,size(E,2));
for i = 1:size(E,2)
    G(:,i) = sum(PZ.G(:,ind==i),2);
end
E(:,all(G==0,1)) = [];
G(:,all(G==0,1)) = [];
PZ = PolynomialZonotope(G,E,PZ.id);

end