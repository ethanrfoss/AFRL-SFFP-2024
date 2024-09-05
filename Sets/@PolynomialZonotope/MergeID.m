function [PZ1,PZ2] = MergeID(PZ1,PZ2)

% Merge and Get unique identifiers
id = unique([PZ1.id PZ2.id]);

% Edit PZs:
E1 = zeros(length(id),PZ1.dims.h);
E1(ismember(id,PZ1.id),:) = PZ1.E;
PZ1 = PolynomialZonotope(PZ1.G,E1,id);
E2 = zeros(length(id),PZ2.dims.h);
E2(ismember(id,PZ2.id),:) = PZ2.E;
PZ2 = PolynomialZonotope(PZ2.G,E2,id);

end