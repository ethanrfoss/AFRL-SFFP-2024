
function Poly = Polygon(PT)

if size(PT.Vertices,2) == 2
    Poly = Polygon2D(PT.Vertices);
else
    Poly = Polygon3D(PT.Vertices,PT.Faces);
end

end