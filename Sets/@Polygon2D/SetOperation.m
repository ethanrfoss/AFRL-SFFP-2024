
function Poly = SetOperation(PolyList,Operation)

warning('off');

polyshapes = repmat(polyshape(),length(PolyList),1);
for i = 1:length(PolyList)
    polyshapes(i) = polyshape(PolyList(i).Vertices(:,1),PolyList(i).Vertices(:,2));
end

if isequal(Operation,'union')
    Poly = union(polyshapes);
elseif isequal(Operation,'intersect')
    Poly = intersect(polyshapes);
else
    error('Unrecognized 2D Polygon Set Operation');
end

Poly = Polygon2D(Poly.Vertices);

warning('on');

end