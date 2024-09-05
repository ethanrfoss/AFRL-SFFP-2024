
function PZ = ExactSum(PZ1,PZ2)

[PZ1,PZ2] = MergeID(PZ1,PZ2);
PZ = PolynomialZonotope([PZ1.G PZ2.G],[PZ1.E PZ2.E],PZ1.id);
% PZ = Compact(PZ);

end