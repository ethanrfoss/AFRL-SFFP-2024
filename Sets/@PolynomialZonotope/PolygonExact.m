
function P = PolygonExact(PZ)



end

function M = MatrixBlocks(E)

ind = 1;
i = 1;
while i < size(E,2)
    j = 1;
    Col = logical(sum(E(:,i),2));
    while i+j <= size(E,2) && all(logical(E(:,i+j)) == Col)
        j = j + 1;
    end
    M{ind} = E(Col,i:i+j-1);
    i = i + j;
    ind = ind + 1;
end

end