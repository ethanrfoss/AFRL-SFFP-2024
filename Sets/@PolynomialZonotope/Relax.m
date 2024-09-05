
function PZ = Relax(PZ,i,j)

ZeroInd = find(all(PZ.E==0,1));
OneInd = find(all(sum(PZ.E,1)==1) & all(sum(~PZ.E,1)==dims.p-1));

k0 = find(PZ.E(:,i)-PZ.E(:,j));
if length(k0) ~= 1
    return;
end

etaE = exp(PZ.E(k0,i)*log(PZ.E(k0,j)/PZ.E(k0,i))/(PZ.E(k0,i)-PZ.E(k0,j)));
if all(mod(PZ.E(:,j),2)==0)
    etaG = etaE/2;
else
    etaG = etaE;
end
etac = etaE-etaG;

G = [etac*PZ.G(:,i) etaG*PZ.G(:,i) PZ.G(:,setdiff(1:PZ.dims.h,[i,j])) PZ.G(:,j)+PZ.G(:,i)];
E = [zeros(PZ.dims.p+1,1) [zeros(PZ.dims.p,1);1] [PZ.E(:,setdiff(1:PZ.dims.h,[i,j])) PZ.E(:,j);zeros(1,PZ.dims.h-1)]];

% PZ = Compact(PolynomialZonotope(G,E,[PZ.id UniqueID(1)]));
PZ = PolynomialZonotope(G,E,[PZ.id UniqueID(1)]);

end