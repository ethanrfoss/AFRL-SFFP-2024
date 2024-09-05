
function PZ = PolynomialZonotope(EP,Options)

DefaultOptions.Order = 3;
DefaultOptions.ID = 'Unique';
DefaultOptions.Method = 'Squircle';
if nargin < 2
    Options = DefaultOptions;
end
Order = Options.Order;

% Validity Checks:
N = length(EP.c);
if size(EP.Q,1) ~= N || size(EP.Q,2) ~= N
    error('Invalid Ellipsoid Dimensions');
end
try
    M = chol(EP.Q);
catch
    error('Q must be Positive Definite');
end

if isequal(Options.Method,'Squircle')
    Coefs = [1 -1/4 -1/32 -1/128 -5/2048 -7/8192];
    Exps = [0 2 4 6 8 10];
    G = kron(eye(N),Coefs(1:Order));
    E = kron(eye(N),ones(1,Order)) + kron(diag(ones(N-1,1),1)+diag(ones(N-1,1),-1),Exps(1:Order));
elseif isequal(Options.Method,'Taylor')
    % Ellipsoid Functions:
    Gcos = @(Order) kron((-1).^(0:Order).*pi.^(2*(0:Order))./factorial(2*(0:Order)),[1 0]);
    Gsin = @(Order) kron((-1).^(0:Order).*pi.^(2*(0:Order)+1)./factorial(2*(0:Order)+1),[0 1]);
    
    % Construct:
    E = zeros(N,1+(2*Order+2)^(N-1));
    E(1,2:end) = 1;
    E(2:end,2:end) = (dec2base(0:(2*Order+2)^(N-1)-1,2*Order+2)-'0')';
    G = zeros(N,(2*Order+2)^(N-1));
    for i = 1:N
        if i == N
            Gr = Gsin(Order);
            for j = 1:N-2
                Gr = kron(Gr,Gsin(Order));
            end
        else
            Gr = Gcos(Order);
            for j = 1:i-1
                Gr = kron(Gr,Gsin(Order));
            end
        end
        G(i,2:length(Gr)) = Gr;

    end
else
    error('Invalid Conversion Method Specified');
end

% Ellipsoid Mapping:
G = [EP.c M*G];
E = [zeros(N,1) E];

% Construct:
if isequal(Options.ID,'Unique')
    PZ = PolynomialZonotope(G,E);
else
    PZ = PolynomialZonotope(G,E,Options.ID);
end

end




