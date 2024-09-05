
function PZ = Subset(PZ,r,l,u)

if l < -1
    error('Invalid Lower Bound');
end
if u > 1
    error('Invalid Upper Bound');
end
if r > PZ.dims.p
    error('Invalid Generator Index');
end

% Get monomial coefficients:
PosInd = find(PZ.E(r,:)>0); % K
NegInd = find(PZ.E(r,:)<=0); % H

% Get Matrices:
E = PZ.E(:,NegInd);
G = PZ.G(:,NegInd);
for i = 1:length(PosInd)
    E = [E [PZ.E(1:r-1,PosInd(i))*ones(1,PZ.E(r,PosInd(i))+1);0:PZ.E(r,PosInd(i));PZ.E(r+1:PZ.dims.p,PosInd(i))*ones(1,PZ.E(r,PosInd(i))+1)]];
    % Monomial Coefficients:
    b = zeros(1,PZ.E(r,PosInd(i))+1);
    for j = 0:PZ.E(r,PosInd(i))
        b(PZ.E(r,PosInd(i))-j+1) = nchoosek(PZ.E(r,PosInd(i)),j)*(.5*(u-l))^(PZ.E(r,PosInd(i))-j)*(.5*(u+l))^j;
    end
    G = [G PZ.G(:,PosInd(i))*b];
end

PZ = Compact(PolynomialZonotope(G,E));

end