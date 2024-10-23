
function b = dec2bij(m,k)

b = '';
q0 = m;
q1 = ceil(q0/k)-1;
if q1 == 0
    b = num2str(m);
    return;
end
while true
    b = [num2str(q0-q1*k) b];
    if ceil(q1/k)-1 == 0
        b = [num2str(q1) b];
        return;
    end
    q0 = q1;
    q1 = ceil(q1/k)-1;
end
    

end