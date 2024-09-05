
A = rand(4,4,4);
B = rand(4,4);
C = rand(4,4);

ABC = zeros(4,4,4);

for i = 1:4
    for a = 1:4
        for b = 1:4
            for alpha = 1:4
                for beta = 1:4
                    ABC(i,a,b) = ABC(i,a,b) + A(i,alpha,beta)*B(alpha,a)*C(beta,b);
                end
            end
        end
    end
end

ABC - tensorprod(tensorprod(A,B,[2],[1]),C,[3],[1])

A = rand(4,4);
B = rand(4,4,4);

AB = zeros(4,4,4);

for i = 1:4
    for a = 1:4
        for b = 1:4
            for alpha = 1:4
                AB(i,a,b) = AB(i,a,b) + A(i,alpha)*B(alpha,a);
            end
        end
    end
end