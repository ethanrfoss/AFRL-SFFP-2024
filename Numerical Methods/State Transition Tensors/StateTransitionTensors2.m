
function [Phi] = StateTransitionTensors2(EOM,T,x0,Options)

% Options:
DefaultOptions.N = 2; % Tensor Order
DefaultOptions.AutomaticDifferentiation = false;
DefaultOptions.DirectionalSTT = false;
DefaultOptions.EigenCutoff = 10^-6;
DefaultOptions.Integrator = @ode45;
DefaultOptions.ODESettings = odeset('AbsTol',10^-6,'RelTol',10^-3);
DefaultOptions.Broken = false;
if nargin < 4
    Options = DefaultOptions;
else
    Options = MergeStructs(Options,DefaultOptions);
end
Options.nx = length(x0);

% Check Dynamics:
[EOM,Options] = DynamicsChecks(EOM,x0,Options);

% Get Derivative Tensors:
if ~Options.Sparse
    if Options.AutomaticDifferentiation
        EOM = AutomaticDerivativeTensors(EOM,Options);
    else
        EOM = DerivativeTensors(EOM,Options);
    end
end

% Direction:
if Options.DirectionalSTT
    [R,Options] = DirectionalRotation(f,Options);
else
    R = eye(Options.nx);
    Options.nk = Options.nx;
end

% Initial Conditions:
V0 = StateTransitionTensorInitialCondition(x0,Options.nx,Options.nk,Options.N);

% Integrate:
Vdot = TensorRate2(EOM,Options);
if ~Options.Broken
    if length(T) == 2
        [~,V] = Options.Integrator(@(t,V) Vdot(V,t),[T(1) T(2)],V0,Options.ODESettings);
        V = [V(1,:)' V(end,:)'];
    else
        [~,V] = Options.Integrator(@(t,V) Vdot(V,t),T,V0,Options.ODESettings);
        V = V';
    end
else
    V = zeros(length(V0),length(T));
    V(:,1) = V0;
    for i = 1:length(T)-1
        % V0 = StateTransitionTensorInitialCondition(V(1:Options.nx,i),Options.nx,Options.N);
        V0 = StateTransitionTensorInitialCondition(V(1:Options.nx,i),Options.nx,Options.nk,Options.N);
        [~,R] = Options.Integrator(@(t,V) Vdot(V,t),[T(i) T(i+1)],V0,Options.ODESettings);
        V(:,i+1) = R(end,:)';
    end
end

% Decompose:
Phi = Tensify(V,Options.nx,Options.N);

end

function [f,Options] = DynamicsChecks(f,x0,Options)

InputDims = InputDimensions(f,sum(Options.nx.^(1:5)));
if InputDims(1) ~= Options.nx
    error('Incommensurate State Dimensions');
end
if isscalar(InputDims)
    f = @(x,t) f(x);
else
    if InputDims(end) ~= 1
        error('Incommensurate Time Dimension');
    end
end
if length(InputDims) == 2 % No Input
    f0 = f(x0,0);
    if length(f0) == Options.nx && Options.N ~= 1
        Options.Sparse = false;
    elseif length(f0) == Options.nx && Options.N == 1
        Options.Sparse = true;
    elseif length(f0) == sum(Options.nx.^(1:Options.N))
        Options.Sparse = true;
    else
        error('Incommensurate EOM Dimensions');
    end
elseif length(InputDims) == 3
    f = @(x,t) [f(x(1:Options.nx),x(Options.nx+1:Options.nx+InputDims(2)),t);zeros(InputDims(2),1)];
    Options.nx = Options.nx+InputDims(2);
    x0 = [x0;zeros(InputDims(2),1)];
    f0 = f(x0,0);
    if length(f0) ~= Options.nx
        error('Incommensurate EOM Dimensions');
    end
    Options.Sparse = false;
else
    error('Invalid Input Dimensions to EOM');
end

end

function [R,Options] = DirectionalRotation(f,Options)

% Get Cauchy Green Tensor:
V0 = StateTransitionTensorInitialCondition(x0,Options.nx,2);
[~,V] = Options.Integrator(@(t,V) TensorRate(t,V,D,Options.nx,2),[T(1) (T(1)+T(end))/2 T(end)],V0,Options.ODESettings);
V = V';
Phi = Tensify(V,Options.nx,2);
[R,S] = eigs(Phi{2,end}'*Phi{2,end});
R = R';
if ~isfield(Options,'nk')
    % Options.nk = [1,Options.nx,repmat(sum(diag(S)/max(S(:))>Options.EigenCutoff),1,Options.N-2)];
    Options.nk = sum(diag(S)/max(S(:))>Options.EigenCutoff);
end

end

