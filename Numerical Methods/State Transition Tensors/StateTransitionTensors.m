
function [Phi] = StateTransitionTensors(EOM,T,x0,Options)

% Options:
DefaultOptions.N = 2; % Tensor Order
DefaultOptions.AutomaticDifferentiation = false;
DefaultOptions.Integrator = @ode45;
DefaultOptions.ODESettings = odeset('AbsTol',10^-15,'RelTol',10^-12);
DefaultOptions.Broken = false;
if nargin < 4
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.nx = length(x0);

% Get Derivative Tensors:
if isfield(Options,'D')
    D = Options.D;
    if issparse(D(x0,T(1)))
        Options.Sparse = true;
    end
elseif Options.AutomaticDifferentiation
    D = AutomaticDerivativeTensors(EOM,Options);
else
    D = DerivativeTensors(EOM,Options);
end

% Initial Conditions:
V0 = StateTransitionTensorInitialCondition(x0,Options.nx,Options.nx,Options.N);

% Integrate:
if ~Options.Broken
    if length(T) == 2
        [~,V] = Options.Integrator(@(t,V) TensorRate(t,V,D,Options.nx,Options.N),[T(1) (T(1)+T(2))/2 T(2)],V0,Options.ODESettings);
        V = [V(1,:)' V(3,:)'];
    else
        [~,V] = Options.Integrator(@(t,V) TensorRate(t,V,D,Options.nx,Options.N),T,V0,Options.ODESettings);
        V = V';
    end
else
    V = zeros(length(V0),length(T));
    V(:,1) = V0;
    for i = 1:length(T)-1
        V0 = StateTransitionTensorInitialCondition(V(1:Options.nx,i),Options.nx,Options.N);
        [~,R] = Options.Integrator(@(t,V) TensorRate(t,V,D,Options.nx,Options.N),[T(i) (T(i)+T(i+1))/2 T(i+1)],V0,Options.ODESettings);
        V(:,i+1) = R(end,:)';
    end
end


% Decompose:
Phi = Tensify(V,Options.nx,Options.N);

end


