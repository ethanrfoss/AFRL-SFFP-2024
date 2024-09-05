
%% N = 2

h = 5;
n = 2;

A = rand(n,n,n);
G = rand(n,h);

Gb = zeros(n,h^2);

% for i = 1:h
%     for j = 1:n
%         Gb(j,h*(i-1)+1:h*i) = G(:,i)'*squeeze(A(j,:,:))*G;
%     end
% end
for l = 1:n
    for i = 1:h
        for j = 1:h
            for k1 = 1:n
                for k2 = 1:n
                    Gb(l,h*(i-1)+j) = Gb(l,h*(i-1)+j) + A(l,k1,k2)*G(k1,i)*G(k2,j);
                end
            end
        end
    end
end

Gb2 = reshape(tensorprod(tensorprod(A,G,[3],[1]),G,[2],[1]),[n,h^2]);

%% N = 3

h = 5;
n = 2;

A = rand(n,n,n,n);
G = rand(n,hW);

% Gb = zeros(n,h^3);

for l = 1:n
    for i = 1:h
        for j = 1:h
            for k = 1:h
                for k1 = 1:n
                    for k2 = 1:n
                        for k3 = 1:n
                            Gb(l,h^2*(i-1)+h*(j-1)+k) = Gb(l,h^2*(i-1)+h*(j-1)+k) + A(l,k1,k2,k3)*G(k1,i)*G(k2,j)*G(k3,k);
                        end
                    end
                end
            end
        end
    end
end


Gb2 = reshape(tensorprod(tensorprod(tensorprod(A,G,[4],[1]),G,[3],[1]),G,[2],[1]),[n,h^3]);