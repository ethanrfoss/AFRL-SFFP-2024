
function Z = Zonotope(PZ)

ZeroInd = all(PZ.E==0,1);
EvenInd = all(mod(PZ.E,2)==0,1) & ~ZeroInd;
OddInd = ~EvenInd & ~ZeroInd;
Z = Zonotope(sum(PZ.G(:,ZeroInd),2)+1/2*sum(PZ.G(:,EvenInd),2),[PZ.G(:,OddInd) 1/2*PZ.G(:,EvenInd)]);
% c = zeros(PZ.dims.n,1);
% G = [];
% for i = 1:PZ.dims.h
%     if all(mod(PZ.E(:,i),2)==0)
%         c = c + 1/2*PZ.G(:,i);
%         G = [G 1/2*PZ.G(:,i)];
%     else
%         G = [G PZ.G(:,i)];
%     end
% end
% 
% Z = Zonotope(c,G);

end