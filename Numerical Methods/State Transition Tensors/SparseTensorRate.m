
function Vdot = SparseTensorRate(f,Options)

if Options.N == 1
    Vdot = @(V,t) SparseTensorRate1D(f,V);
elseif Options.N == 2
    Vdot = @(V,t) SparseTensorRate2D(f,V);
elseif Options.N == 3
    Vdot = @(V,t) SparseTensorRate3D(f,V);
elseif Options.N == 4
    Vdot = @(V,t) SparseTensorRate4D(f,V);
else
    error('Invalid Order for Sparse Tensor Rate');
end

end

function Vdot = SparseTensorRate1D(f,V,t)

Vdot = f(V,t);

end

function Vdot = SparseTensorRate2D(f,V,t,nx)

[rows,~,fx] = find(f(V(1:nx),t));
Vdot = zeros(nx+nx^2,1);
ind = [0,nx];


for j = 1:length(rows)
    Inds = dec2base(row(j),nx)-'0';
    if length(Inds) <= 1
        Vdot(Inds) = fx(j);
    else
        for a = 1:nx
            Vdot(ind(2)+Inds(1)*a) = Vdot(ind(2)+Inds(1)*a) + fx(j)*V(ind(2)+Inds(2)*a);
        end
    end

end

end

function Vdot = SparseTensorRate3D(f,V,t,nx)

[rows,~,fx] = find(f(V(1:nx),t));
Vdot = zeros(nx+nx^2,1);
ind = [0,nx,nx^2];


for j = 1:length(rows)
    Inds = dec2base(row(j),nx)-'0';
    if length(Inds) <= 1
        Vdot(Inds) = fx(j);
    elseif length(Inds) <= 2
        for a = 1:nx
            Vdot(ind(2)+Inds(1)*a) = Vdot(ind(2)+Inds(1)*a) + fx(j)*V(ind(2)+Inds(2)*a);
            for b = 1:nx
                Vdot(ind(3)+Inds(1)*a*b) = Vdot(ind(3)+Inds(1)*a*b) + fx(j)*V(ind(3)+Inds(2)*a*b);
            end
        end
    else
        for a = 1:nx
            for b = 1:nx
                Vdot(ind(3)+Inds(1)*a*b) = Vdot(ind(3)+Inds(1)*a*b) + fx(j)*V(ind(2)+Inds(2)*a)*V(ind(2)+Ind(3)*b);
            end
        end
    end

end

end

function Vdot = SparseTensorRate4D(f,V,t,nx)

[rows,~,fx] = find(f(V(1:nx),t));
Vdot = zeros(nx+nx^2,1);
ind = [0,nx,nx^2,nx^3];


for j = 1:length(rows)
    Inds = dec2base(row(j),nx)-'0';
    if length(Inds) <= 1
        Vdot(Inds) = fx(j);
    elseif length(Inds) <= 2
        for a = 1:nx
            Vdot(ind(2)+Inds(1)*a) = Vdot(ind(2)+Inds(1)*a) + fx(j)*V(ind(2)+Inds(2)*a);
            for b = 1:nx
                Vdot(ind(3)+Inds(1)*a*b) = Vdot(ind(3)+Inds(1)*a*b) + fx(j)*V(ind(3)+Inds(2)*a*b);
                for c = 1:nx
                    Vdot(ind(4)+Inds(1)*a*b*c) = Vdot(ind(4)+Inds(1)*a*b*c) + fx(j)*V(ind(4)+Inds(2)*a*b*c);
                end
            end
        end
    elseif length(Inds) <= 3
        for a = 1:nx
            for b = 1:nx
                Vdot(ind(3)+Inds(1)*a*b) = Vdot(ind(3)+Inds(1)*a*b) + fx(j)*V(ind(2)+Inds(2)*a)*V(ind(2)+Ind(3)*b);
                for c = 1:nx
                    Vdot(ind(4)+Inds(1)*a*b*c) = Vdot(ind(4)+Inds(1)*a*b*c) + fx(j)*(V(ind(2)+Inds(2)*a)*V(ind(3)+Inds(3)*b*c)+V(ind(3)+Inds(1)*a*b)*V(ind(2)+Inds(3)*c)+V(ind(3)+Inds(2)*a*c)*V(ind(2)+Inds(3)*b));
                end
            end
        end
    else
        for a = 1:nx
            for b = 1:nx
                for c = 1:nx
                    Vdot(ind(4)+Inds(1)*a*b*c) = Vdot(ind(4)+Inds(1)*a*b*c) + fx(j)*V(ind(2)+Inds(2)*a)*V(ind(2)+Inds(3)*b)*V(ind(2)+Inds(4)*c);
                end
            end
        end
    end

end

end
