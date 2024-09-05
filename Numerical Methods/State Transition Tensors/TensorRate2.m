
function Vdot = TensorRate2(f,Options)

% Vdot = @(V,t) StandardTensorRate(t,V,f,Options.nx,Options.N);
if isa(f,'cell')
    Vdot = @(V,t) StandardTensorRate(t,V,f,Options.nx,Options.N);
elseif isfield(Options,'Sparse')
    if Options.N == 1
        Vdot = @(V,t) SparseTensorRate1D(f,V,t);
    elseif Options.N == 2
        Vdot = @(V,t) SparseTensorRate2D(f,V,t,Options.nx);
    elseif Options.N == 3
        Vdot = @(V,t) SparseTensorRate3D(f,V,t,Options.nx,Options.nk);
    elseif Options.N == 4
        Vdot = @(V,t) SparseTensorRate4D(f,V,t,Options.nx,Options.nk);
    else
        error('Invalid Order for Sparse Tensor Rate');
    end
else
    error('Invalid Derivative Tensor');
end

end


function Vdot = StandardTensorRate(t,V,D,nx,N)

disp(['Current Time: ' num2str(t)]);

% Convert Vector to Tensor Form:
Phi = Tensify(V,nx,N);
% f0 = full(D(V(1:nx),t));

Phidot = cell(N,1);
if N >= 1
    % Phidot{1} = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{1},'UniformOutput',false));
    Phidot{1} = D{1}(t,Phi{1});
    % Phidot{1} = f0(1:nx);
end
if N >= 2
    % A2 = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{2},'UniformOutput',false));
    A2 = D{2}(t,Phi{1});
    % A2 = reshape(f0(nx+1:nx+nx^2),[nx,nx]);
    % Phidot{2} = A2*Phi{2};
    Phidot{2} = zeros(repmat(nx,1,2));
    for i = 1:nx
        for a = 1:nx
            for alpha = 1:nx
                Phidot{2}(i,a) = Phidot{2}(i,a) + A2(i,alpha)*Phi{2}(alpha,a);
            end
        end
    end
end
if N >= 3
    % A3 = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{3},'UniformOutput',false));
    A3 = D{3}(t,Phi{1});
    % A3 = reshape(f0(nx+nx^2+1:nx+nx^2+nx^3),[nx,nx,nx]);
    Phidot{3} = zeros(repmat(nx,1,3));
    for i = 1:nx
        for a = 1:nx
            for b = 1:nx
                for alpha = 1:nx
                    Phidot{3}(i,a,b) = Phidot{3}(i,a,b) + A2(i,alpha)*Phi{3}(alpha,a,b);
                    for beta = 1:nx
                        Phidot{3}(i,a,b) = Phidot{3}(i,a,b) + A3(i,alpha,beta)*Phi{2}(alpha,a)*Phi{2}(beta,b);
                    end
                end
            end
        end
    end
end
if N >= 4
    % A4 = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{4},'UniformOutput',false));
    A4 = D{4}(t,Phi{1});
    % A4 = reshape(f0(nx+nx^2+nx^3+1:nx+nx^2+nx^3+nx^4),[nx,nx,nx,nx]);
    Phidot{4} = zeros(repmat(nx,1,4));
    for i = 1:nx
        for a = 1:nx
            for b = 1:nx
                for c = 1:nx
                    for alpha = 1:nx
                        Phidot{4}(i,a,b,c) = Phidot{4}(i,a,b,c) + A2(i,alpha)*Phi{4}(alpha,a,b,c);
                        for beta = 1:nx
                            Phidot{4}(i,a,b,c) = Phidot{4}(i,a,b,c) + A3(i,alpha,beta)*(Phi{2}(alpha,a)*Phi{3}(beta,b,c)+Phi{3}(alpha,a,b)*Phi{2}(beta,c)+Phi{3}(alpha,a,c)*Phi{2}(beta,b));
                            for gamma = 1:nx
                                Phidot{4}(i,a,b,c) = Phidot{4}(i,a,b,c) + A4(i,alpha,beta,gamma)*Phi{2}(alpha,a)*Phi{2}(beta,b)*Phi{2}(gamma,c);
                            end
                        end
                    end
                end
            end
        end
    end
end
if N >= 5
    % A5 = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{5},'UniformOutput',false));
    A5 = D{5}(t,Phi{1});
    Phidot{5} = zeros(repmat(nx,1,5));
    for i = 1:nx
        for a = 1:nx
            for b = 1:nx
                for c = 1:nx
                    for d = 1:nx
                        for alpha = 1:nx
                            Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A2(i,alpha)*Phi{5}(alpha,a,b,c,d);
                            for beta = 1:nx
                                Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A3(i,alpha,beta)*(Phi{4}(alpha,a,b,c)*Phi{2}(beta,d)+Phi{4}(alpha,a,b,d)*Phi{2}(beta,c)+Phi{4}(alpha,a,c,d)*Phi{2}(beta,b)+Phi{3}(alpha,a,b)*Phi{3}(beta,c,d)+Phi{3}(alpha,a,c)*Phi{3}(beta,b,d)+Phi{3}(alpha,a,d)*Phi{3}(beta,b,c)+Phi{2}(alpha,a)*Phi{4}(beta,b,c,d));
                                for gamma = 1:nx
                                    Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A4(i,alpha,beta,gamma)*(Phi{3}(alpha,a,b)*Phi{2}(beta,c)*Phi{2}(gamma,d)+Phi{3}(alpha,a,c)*Phi{2}(beta,b)*Phi{2}(gamma,d)+Phi{3}(alpha,a,d)*Phi{2}(beta,b)*Phi{2}(gamma,c)+Phi{2}(alpha,a)*Phi{3}(beta,b,c)*Phi{2}(gamma,d)+Phi{2}(alpha,a)*Phi{3}(beta,b,d)*Phi{2}(gamma,c)+Phi{2}(alpha,a)*Phi{2}(beta,b)*Phi{3}(gamma,c,d));
                                    for delta = 1:nx
                                        Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A5(i,alpha,beta,gamma,delta)*Phi{2}(alpha,a)*Phi{2}(beta,b)*Phi{2}(gamma,c)*Phi{3}(delta,d);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if N >= 6
    error('This order state transition tensor has not been programmed');
end

% Convert Tensor to Vector Form:
Vdot = Detensify(Phidot);

end

function Vdot = SparseTensorRate1D(f,V,t)

disp(['Current Time: ' num2str(t)]);
Vdot = f(V,t);

end

function Vdot = SparseTensorRate2D(f,V,t,nx)

disp(['Current Time: ' num2str(t)]);

[row,~,fx] = find(f(V(1:nx),t));
Vdot = zeros(nx+nx^2,1);

% f0 = f(V(1:nx));
% A2 = reshape(f0(nx+1:nx+nx^2),[nx,nx]);

for j = 1:length(row)
    % Inds = dec2bij(row(j),nx)-'0';
    Inds = [];
    q0 = row(j);
    q1 = ceil(q0/nx)-1;
    if q1 == 0
        Inds = row(j);
    else
        while true
            Inds = [Inds (q0-q1*nx)];
            if ceil(q1/nx)-1 == 0
                Inds = [Inds q1];
                break;
            end
            q0 = q1;
            q1 = ceil(q1/nx)-1;
        end
    end
    if length(Inds) <= 1
        Vdot(Inds) = fx(j);
    else
        for a = 1:nx
            % Vdot(ind(2)+(Inds(1)-1)*nx+a) = Vdot(ind(2)+(Inds(1)-1)*nx+a) + fx(j)*V(ind(2)+(Inds(2)-1)*nx+a);
            % Vdot(bij2base(Inds(1),a,nx)) = Vdot(bij2base(Inds(1),a,nx)) + fx(j)*V(bij2base(Inds(2),a,nx));
            Vdot((Inds(1))+a*nx) = Vdot((Inds(1))+a*nx) + fx(j)*V(Inds(2)+a*nx); % Fastest
        end
    end

end

% fx = f(V(1:nx),t);
% Vdot2 = [fx(1:nx);reshape(reshape(fx(nx+1:nx+nx^2),[nx,nx])*reshape(V(nx+1:nx+nx^2),[nx,nx]),[nx^2,1])];

end

function Vdot = SparseTensorRate3D(f,V,t,nx,nk)

disp(['Current Time: ' num2str(t)]);

f0 = f(V(1:nx),t);
[row,~,fx] = find([zeros(nx,1);f0(nx+1:end)]);
Vdot = zeros(nx+nx^2+nx*nk^2,1);
Vdot(1:nx) = f0(1:nx);
Vdot(nx+1:nx+nx^2) = reshape(reshape(f0(nx+1:nx+nx^2),[nx,nx])*reshape(V(nx+1:nx+nx^2),[nx,nx]),[nx^2,1]);

% A3 = reshape(V(nx+nx^2+1:end),[nx,nx,nx]);

for j = 1:length(row)
    Inds = [];
    q0 = row(j);
    q1 = ceil(q0/nx)-1;
    if q1 == 0
        Inds = row(j);
    else
        while true
            Inds = [Inds (q0-q1*nx)];
            if ceil(q1/nx)-1 == 0
                Inds = [Inds q1];
                break;
            end
            q0 = q1;
            q1 = ceil(q1/nx)-1;
        end
    end
    if length(Inds) == 2
        for a = 1:nk
            for b = 1:nk
                Vdot(Inds(1)+a*nx+b*nx^2) = Vdot(Inds(1)+a*nx+b*nx^2) + fx(j)*V(Inds(2)+a*nx+b*nx^2);
                % A3(Inds(2),a,b) - V(b*nx^2+a*nx+Inds(2))
            end
        end
    else
        for a = 1:nk
            for b = 1:nk
                Vdot(Inds(1)+a*nx+b*nx^2) = Vdot(Inds(1)+a*nx+b*nx^2) + fx(j)*V(Inds(2)+a*nx)*V(Inds(3)+b*nx);
            end
        end
    end

end

end

function Vdot = SparseTensorRate4D(f,V,t,nx,nk)

disp(['Current Time: ' num2str(t)]);

f0 = f(V(1:nx),t);
[row,~,fx] = find([zeros(nx,1);f0(nx+1:end)]);
Vdot = zeros(nx+nx^2+nx*nk^2+nx*nk^3,1);
Vdot(1:nx) = f0(1:nx);
Vdot(nx+1:nx+nx^2) = reshape(reshape(f0(nx+1:nx+nx^2),[nx,nx])*reshape(V(nx+1:nx+nx^2),[nx,nx]),[nx^2,1]);

for j = 1:length(row)
    Inds = [];
    q0 = row(j);
    q1 = ceil(q0/nx)-1;
    if q1 == 0
        Inds = row(j);
    else
        while true
            Inds = [Inds (q0-q1*nx)];
            if ceil(q1/nx)-1 == 0
                Inds = [Inds q1];
                break;
            end
            q0 = q1;
            q1 = ceil(q1/nx)-1;
        end
    end
    if length(Inds) == 2
        for a = 1:nk
            for b = 1:nk
                Vdot(Inds(1)+a*nx+b*nx^2) = Vdot(Inds(1)+a*nx+b*nx^2) + fx(j)*V(Inds(2)+a*nx+b*nx^2);
                for c = 1:nk
                    Vdot(Inds(1)+a*nx+b*nx^2+c*nx^3) = Vdot(Inds(1)+a*nx+b*nx^2+c*nx^3) + fx(j)*V(Inds(2)+a*nx+b*nx^2+c*nx^3);
                end
            end
        end
    elseif length(Inds) == 3
        for a = 1:nk
            for b = 1:nk
                Vdot(Inds(1)+a*nx+b*nx^2) = Vdot(Inds(1)+a*nx+b*nx^2) + fx(j)*V(Inds(2)+a*nx)*V(Inds(3)+b*nx);
                for c = 1:nk
                    Vdot(Inds(1)+a*nx+b*nx^2+c*nx^3) = Vdot(Inds(1)+a*nx+b*nx^2+c*nx^3) + fx(j)*(V(Inds(2)+a*nx)*V(Inds(3)+b*nx+c*nx^2)+V(Inds(2)+a*nx+b*nx^2)*V(Inds(3)+c*nx)+V(Inds(2)+a*nx+c*nx^2)*V(Inds(3)+b*nx));
                end
            end
        end
    else
        for a = 1:nk
            for b = 1:nk
                for c = 1:nk
                    Vdot(Inds(1)+a*nx+b*nx^2+c*nx^3) = Vdot(Inds(1)+a*nx+b*nx^2+c*nx^3) + fx(j)*V(Inds(2)+a*nx)*V(Inds(3)+b*nx)*V(Inds(4)+c*nx);
                end
            end
        end
    end

end

end