
function PZ = vertcat(varargin)

PZ = PolynomialZonotope([],[]); % Initialize With Empty Zonotope;    
for i = 1:length(varargin)
    if isa(varargin{i},'PolynomialZonotope')
        PZ = ExactCartesianProduct(PZ,varargin{i});
    elseif isnumeric(varargin{i})
        PZ = ExactCartesianProduct(PZ,PolynomialZonotope(varargin{i},zeros(0,1)));
    else
        try
            PZ = ExactCartesianProduct(PZ,PolynomialZonotope(varargin{i}));
        catch
            error('Invalid Input for Concatenation');
        end
    end
end