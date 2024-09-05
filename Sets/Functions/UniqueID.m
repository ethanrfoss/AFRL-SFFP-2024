
function id = UniqueID(n)

persistent Counter;
if nargin == 0
    Counter = 1;
    return;
end
if isempty(Counter)
    Counter = 1;
end

id = Counter:Counter+n-1;
Counter = Counter + n;

end