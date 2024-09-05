
function P = RandomPoints(PZ,N)

P = zeros(PZ.dims.n,N);
for i = 1:N
    alpha = rand(1,PZ.dims.p)*2-1;
    P(:,i) = PZ.G*prod(alpha.^(PZ.E'),2);
end

end