function PZ = LinearMap(PZ,M)

if size(M,2) ~= PZ.dims.n
    error('Incommensurate Dimensions for Linear Mapping');
end
PZ = PolynomialZonotope(M*PZ.G,PZ.E,PZ.id);

end