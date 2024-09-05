
function PZ = Reduce(PZ,rho,Options)

% Fix so that Zero Index not included!

DefaultOptions.Method = 'PCA';
if nargin < 3
    Options = DefaultOptions;
end

a = max(0,min(PZ.dims.h,ceil(PZ.dims.h-PZ.dims.n*(rho-1))));
[~,SortInd] = sort(vecnorm(PZ.G,2,1));
G = PZ.G(:,SortInd);
E = PZ.E(:,SortInd);
K = 1:a;
Kbar = a+1:size(PZ.E,2);
Z = Reduce(Zonotope(PolynomialZonotope(G(:,K),E(:,K),PZ.id)),1,Options);

Ebar = E(:,Kbar);
N = any(Ebar~=0,2);
Ebar = Ebar(N,:);
PZ = Compact(PolynomialZonotope([Z.c G(:,Kbar) Z.G],blkdiag([zeros(size(Ebar,1),1) Ebar],eye(size(Z.G,2))),[PZ.id(N) Z.id]));

end