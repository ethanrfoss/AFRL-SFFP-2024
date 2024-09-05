
function r = MaxSplitDirection(PZ)

ZeroInd = all(PZ.E==0,1);
% OneInd = sum(PZ.E,1)==1;
[~,ind] = max(sum((PZ.G.*~ZeroInd).^2,1));
[~,r] = max(PZ.E(:,ind));

end