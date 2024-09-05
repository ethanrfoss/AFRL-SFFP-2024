function V = Detensify(T)

ind = 1;
V = zeros(sum(cellfun(@numel,T)),1);
if isa(T{ind},'sym')
    V = sym(V);
end
for i = 1:length(T)
    V(ind:ind+numel(T{i})-1) = reshape(T{i},[numel(T{i}),1]);
    ind = ind+numel(T{i});
end

end