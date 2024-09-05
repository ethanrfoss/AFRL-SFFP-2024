
function Poly = Polygon(Z)

% alpha = (dec2bin(0:2^(Z.dims.l)-1)-'0')*2-1;
% P = (Z.c + Z.G*alpha')';
% if Z.dims.n == 2
%     P = P(convhull(P(:,1),P(:,2)),:);
%     Poly = Polygon2D(P);
% elseif Z.dims.n == 3
%     P = P(convhull(P(:,1),P(:,2),P(:,3)),:);
%     Poly = Polygon3D(P);
% else
%     error('Invalid Polygon Conversion Dimension');
% end

if Z.dims.n == 2
    Poly = Poly2D(Z);
elseif Z.dims.n == 3
    error('This functionality does not exist yet');
else
    error('Invalid Polygon Conversion Dimension');
end

end

function Poly = Poly2D(Z)

G = Z.G;
G(:,G(2,:)<0) = -G(:,G(2,:)<0);
[~,ind] = sort(atan2(G(2,:),G(1,:)));
MaxG = sum(G,2);

P = zeros(2,size(G,2));
P(:,1) = 2*G(:,ind(1));
for i = 2:size(P,2)
    P(:,i) = P(:,i-1) + 2*G(:,ind(i));
end

P = P - MaxG;
P = [P -P];

Poly = Polygon2D((P+Z.c)');

end