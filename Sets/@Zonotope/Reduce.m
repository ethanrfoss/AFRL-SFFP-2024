
function Z = Reduce(Z,rho,Options)

DefaultOptions.Method = 'PCA';
if nargin < 3
    Options = DefaultOptions;
end

if isequal(Options.Method,'Girard')
    [~,SortInd] = sort(vecnorm(Z.G,1,1)-vecnorm(Z.G,inf,1));
    G = Z.G(:,SortInd);
    H = 1:ceil(max(0,Z.dims.l-Z.dims.n*(rho-1)));
    K = ceil(max(0,Z.dims.l-Z.dims.n*(rho-1)))+1:Z.dims.l;
    Z = Zonotope(Z.c,[G(:,K) diag(sum(abs(G(:,H)),2))]);
elseif isequal(Options.Method,'PCA')
    [~,SortInd] = sort(vecnorm(Z.G,1,1)-vecnorm(Z.G,inf,1));
    G = Z.G(:,SortInd);
    H = 1:ceil(max(0,Z.dims.l-Z.dims.n*(rho-1)));
    K = ceil(max(0,Z.dims.l-Z.dims.n*(rho-1)))+1:Z.dims.l;
    U = PrincipleComponentAnalysis([-G(:,H) G(:,H)]);
    Z = Zonotope(Z.c,[G(:,K) U*diag(sum(abs(U'*G(:,H)),2))]);
else
    error('Invalid Reduction Method');
end


end