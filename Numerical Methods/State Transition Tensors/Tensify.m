function T = Tensify(V,nx,N,nk)

T = cell(N,size(V,2));
if nargin > 3
    ind = [1,nx*(nk.^(0:N-1))];
else
    ind = nx.^(0:N);
end
for i = 1:N
    for j = 1:size(V,2)
        if i >= 2
            if nargin > 3
                shape = [nx,repmat(nk(i),1,i-1)];
            else
                shape = repmat(nx,1,i);
            end
            T{i,j} = reshape(V(sum(ind(1:i)):sum(ind(1:i+1))-1,j),shape);
        else
            T{i,j} = V(sum(ind(1:i)):sum(ind(1:i+1))-1,j);
        end
    end
end

end