
function PZ = Project(PZ,ind)

PZ = Compact(PolynomialZonotope(PZ.G(ind,:),PZ.E,PZ.id));

end