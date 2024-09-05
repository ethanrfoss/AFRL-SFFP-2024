
function U = PrincipleComponentAnalysis(X)

[U,~,~] = svd(cov(X'));

end