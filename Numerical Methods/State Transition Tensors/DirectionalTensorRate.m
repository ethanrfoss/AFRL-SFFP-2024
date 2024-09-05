
function Vdot = DirectionalTensorRate(t,V,D,R,nx,N,nk)

disp(['Current Time: ' num2str(t)]);

% Convert Vector to Tensor Form:
Phi = Tensify(V,nx,N,nk);

Phidot = cell(N,1);
if N >= 1
    % Phidot{1} = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{1},'UniformOutput',false));
    Phidot{1} = D{1}(t,Phi{1});
end
if N >= 2
    % A2 = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{2},'UniformOutput',false));
    A2 = D{2}(t,Phi{1});
    % Phidot{2} = A2*Phi{2};
    Phidot{2} = zeros(nx,nk(2));
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
    Phi2 = Phi{2}*R(1:nk(3),:)';
    Phidot{3} = zeros([nx,repmat(nk(3),1,2)]);
    for i = 1:nx
        for a = 1:nk(3)
            for b = 1:nk(3)
                for alpha = 1:nx
                    Phidot{3}(i,a,b) = Phidot{3}(i,a,b) + A2(i,alpha)*Phi{3}(alpha,a,b);
                    for beta = 1:nx
                        Phidot{3}(i,a,b) = Phidot{3}(i,a,b) + A3(i,alpha,beta)*Phi2(alpha,a)*Phi2(beta,b);
                    end
                end
            end
        end
    end
end
if N >= 4
    % A4 = cell2mat(cellfun(@(f) full(f(t,Phi{1})),D{4},'UniformOutput',false));
    A4 = D{4}(t,Phi{1});
    Phi2 = Phi{2}*R(1:nk(4),:)';
    Phidot{4} = zeros([nx,repmat(nk(4),1,3)]);
    for i = 1:nx
        for a = 1:nk(4)
            for b = 1:nk(4)
                for c = 1:nk(4)
                    for alpha = 1:nx
                        Phidot{4}(i,a,b,c) = Phidot{4}(i,a,b,c) + A2(i,alpha)*Phi{4}(alpha,a,b,c);
                        for beta = 1:nx
                            Phidot{4}(i,a,b,c) = Phidot{4}(i,a,b,c) + A3(i,alpha,beta)*(Phi2(alpha,a)*Phi{3}(beta,b,c)+Phi{3}(alpha,a,b)*Phi2(beta,c)+Phi{3}(alpha,a,c)*Phi2(beta,b));
                            for gamma = 1:nx
                                Phidot{4}(i,a,b,c) = Phidot{4}(i,a,b,c) + A4(i,alpha,beta,gamma)*Phi2(alpha,a)*Phi2(beta,b)*Phi2(gamma,c);
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
    Phi2 = Phi{2}*R(1:nk(5),:)';
    Phidot{5} = zeros([nx,repmat(nk(5),1,4)]);
    for i = 1:nx
        for a = 1:nk(5)
            for b = 1:nk(5)
                for c = 1:nk(5)
                    for d = 1:nk(5)
                        for alpha = 1:nx
                            Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A2(i,alpha)*Phi{5}(alpha,a,b,c,d);
                            for beta = 1:nx
                                Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A3(i,alpha,beta)*(Phi{4}(alpha,a,b,c)*Phi2(beta,d)+Phi{4}(alpha,a,b,d)*Phi2(beta,c)+Phi{4}(alpha,a,c,d)*Phi2(beta,b)+Phi{3}(alpha,a,b)*Phi{3}(beta,c,d)+Phi{3}(alpha,a,c)*Phi{3}(beta,b,d)+Phi{3}(alpha,a,d)*Phi{3}(beta,b,c)+Phi2(alpha,a)*Phi{4}(beta,b,c,d));
                                for gamma = 1:nx
                                    Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A4(i,alpha,beta,gamma)*(Phi{3}(alpha,a,b)*Phi2(beta,c)*Phi2(gamma,d)+Phi{3}(alpha,a,c)*Phi2(beta,b)*Phi2(gamma,d)+Phi{3}(alpha,a,d)*Phi2(beta,b)*Phi2(gamma,c)+Phi2(alpha,a)*Phi{3}(beta,b,c)*Phi2(gamma,d)+Phi2(alpha,a)*Phi{3}(beta,b,d)*Phi2(gamma,c)+Phi{2}(alpha,a)*Phi2(beta,b)*Phi{3}(gamma,c,d));
                                    for delta = 1:nx
                                        Phidot{5}(i,a,b,c,d) = Phidot{5}(i,a,b,c,d) + A5(i,alpha,beta,gamma,delta)*Phi2(alpha,a)*Phi2(beta,b)*Phi2(gamma,c)*Phi{3}(delta,d);
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