
function PZ = CartesianProduct(PZ1,PZ2)

PZ = Compact(PolynomialZonotope(blkdiag(PZ1.G,PZ2.G),blkdiag(PZ1.E,PZ2.E),UniqueID(PZ1.dims.p+PZ2.dims.p)));

end