
function Z = RegularPolygon(N)

if mod(N,2) ~= 0
    error('N must be even');
end
angles = (0:N/2-1)*360/N;
L = tand(180/N);

Z = Zonotope([0;0],[cosd(angles);sind(angles)]*L);

end