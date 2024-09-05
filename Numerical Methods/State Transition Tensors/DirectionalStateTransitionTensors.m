
function [Phi,R] = DirectionalStateTransitionTensors(EOM,T,x0,Options)

% Options:
DefaultOptions.N = 3; % Tensor Order
DefaultOptions.AutomaticDifferentiation = false;
DefaultOptions.Integrator = @ode45;
DefaultOptions.ODESettings = odeset('AbsTol',10^-15,'RelTol',10^-12);
DefaultOptions.EigenCutoff = 10^-6;
if nargin < 4
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.nx = length(x0);

% Get Derivative Tensors:
if isfield(Options,'D')
    D = Options.D;
elseif Options.AutomaticDifferentiation
    D = AutomaticDerivativeTensors(EOM,Options);
else
    D = DerivativeTensors(EOM,Options);
end

% Get Cauchy Green Tensor:
V0 = StateTransitionTensorInitialCondition(x0,Options.nx,2);
[~,V] = Options.Integrator(@(t,V) TensorRate(t,V,D,Options.nx,2),[T(1) (T(1)+T(end))/2 T(end)],V0,Options.ODESettings);
V = V';
Phi = Tensify(V,Options.nx,2);
[R,S] = eigs(Phi{2,end}'*Phi{2,end});
R = R';
if ~isfield(Options,'nk')
    Options.nk = [1,Options.nx,repmat(sum(diag(S)/max(S(:))>Options.EigenCutoff),1,Options.N-2)];
end

% Integrate:
V0 = StateTransitionTensorInitialCondition(x0,Options.nx,Options.N);
V0 = V0(1:sum(Options.nx*(Options.nk.^(0:Options.N-1))));
if length(T) == 2
    [~,V] = Options.Integrator(@(t,V) DirectionalTensorRate(t,V,D,R,Options.nx,Options.N,Options.nk),[T(1) (T(1)+T(2))/2 T(2)],V0,Options.ODESettings);
    V = [V(1,:)' V(3,:)'];
else
    [~,V] = Options.Integrator(@(t,V) DirectionalTensorRate(t,V,D,R,Options.nx,Options.N,Options.nk),T,V0,Options.ODESettings);
    V = V';
end

% Decompose:
Phi = Tensify(V,Options.nx,Options.N,Options.nk);

end


