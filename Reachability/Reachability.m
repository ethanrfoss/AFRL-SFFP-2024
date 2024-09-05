
function [ReachSet,InputSet] = Reachability(EOM,T,x0,u0,X0,U,Options)

DefaultOptions.STTOrder = 4;
DefaultOptions.DirectionalSTT = false;
DefaultOptions.EigenCutoff = 10^-6;
DefaultOptions.MaxPZOrder = 100;
DefaultOptions.ReductionOrder = 100;
DefaultOptions.ReductionOptions.Method = 'PCA';
if nargin < 7
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.nx = length(x0);
Options.nu = length(u0);

% Integrate STTs:
STTOptions.Broken = true;
STTOptions.N = Options.STTOrder;
STTOptions.EigenCutoff = Options.EigenCutoff;
if Options.DirectionalSTT
    [Phi,G] = DirectionalStateTransitionTensors(@(t,z)[EOM(z(1:Options.nx),z(Options.nx+1:Options.nx+Options.nu),t);zeros(Options.nu,1)],T,[x0;u0],STTOptions);
else
    % Phi = StateTransitionTensors(@(t,z)[EOM(z(1:Options.nx),z(Options.nx+1:Options.nx+Options.nu),t);zeros(Options.nu,1)],T,[x0;u0],STTOptions);
    Phi = StateTransitionTensors2(EOM,T,[x0;u0],STTOptions);
    G = eye(Options.nx+Options.nu);
end
for i = 1:STTOptions.N
    ind = repmat({':'},1,i);
    ind{1} = Options.nx+1:Options.nx+Options.nu;
    for j = 1:length(T)
        Phi{i,j}(ind{:}) = [];
    end
end
% for i = 1:Options.nx+Options.nu
%     for j = 1:length(T)
%         Q{i,j} = squeeze(Phi{3,j}(i,:,:));
%     end
% end

% Initialize and Propagate Reachable Set:
ReachSet(1) = PolynomialZonotope(X0);
for i = 2:length(T)

    % Propagate:
    InputSet(i-1) = PolynomialZonotope(U);
    dZ = [ReachSet(i-1)-Phi{1,i-1};InputSet(i-1)-u0];
    % R(i) = Phi{1,i} + LinearMap(dZ,Phi{2,i});
    % dZ = LinearMap(dZ,G);
    % ReachSet(i) = Phi{1,i} + Phi{2,i}*dZ;
    % dZ = G*dZ;
    % for j = 3:Options.STTOrder
    %     % R(i) = R(i) + TensorMap(dZ,1/factorial(j-1)*Phi{j,i});
    %     % R(1:size(Phi{j,i},ndims(Phi{j,i})),:)
    %     ReachSet(i) = ReachSet(i) + 1/factorial(j-1)*Phi{j,i}*(G(1:size(Phi{j,i},ndims(Phi{j,i})),:)*dZ);
    % end
    ReachSet(i) = TensorSum(dZ,Phi(2:end,i),Phi{1,i});

    % Compact:
    % ReachSet(i) = Compact(ReachSet(i));

    % Reduce:
    if ReachSet(i).dims.h >= ReachSet(i).dims.n*Options.MaxPZOrder
        % ReachSet(i) = GreedyRelax(ReachSet(i),2);
        ReachSet(i) = Reduce(ReachSet(i),Options.ReductionOrder,Options.ReductionOptions);
    end

    disp(['Current Time: ' num2str(T(i))]);

end

% % Cora Version:
% ReachSet(1) = polyZonotope(X0);
% for i = 2:length(T)
% 
%     % Propagate:
%     InputSet(i-1) = polyZonotope(U.c,U.G,U.GI,U.E,UniqueID(2)');
%     dZ = cartProd_(ReachSet(i-1)-Phi{1,i-1}(1:Options.nx),InputSet(i-1)-u0);
%     % R(i) = Phi{1,i} + LinearMap(dZ,Phi{2,i});
%     % dZ = LinearMap(dZ,G);
%     ReachSet(i) = Phi{1,i} + Phi{2,i}*dZ;
%     % dZ = G*dZ;
%     ReachSet(i) = ReachSet(i) + 1/factorial(2)*quadMap(dZ,Q(:,i));
%     ReachSet(i) = project(ReachSet(i),[1:Options.nx]);
% 
%     % Compact:
%     ReachSet(i) = compact(ReachSet(i));
% 
%     % Reduce:
%     if size(ReachSet(i).E,2) >= size(ReachSet(i).G,1)*Options.MaxPZOrder
%         % ReachSet(i) = GreedyRelax(ReachSet(i),2);
%         % ReachSet(i) = Reduce(ReachSet(i),Options.ReductionOrder,Options.ReductionOptions);
%         ReachSet(i) = reduce(ReachSet(i),'pca',Options.ReductionOrder);
%     end
% 
%     disp(['Current Time: ' num2str(T(i))]);
% 
% end

end