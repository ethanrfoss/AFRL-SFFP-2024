
function PZ = plus(P1,P2)

if isa(P1,'PolynomialZonotope') && isa(P2,'PolynomialZonotope')
    PZ = ExactSum(P1,P2);
elseif isa(P1,'PolynomialZonotope') && isnumeric(P2)
    if P1.dims.n ~= size(P2,1)
        error('Incorrect Dimensions for subtraction');
    end
    PZ = Compact(PolynomialZonotope([P2 P1.G],[zeros(P1.dims.p,1) P1.E],P1.id));
elseif isa(P2,'PolynomialZonotope') && isnumeric(P1)
    if P2.dims.n ~= size(P1,1)
        error('Incorrect Dimensions for subtraction');
    end
    PZ = Compact(PolynomialZonotope([P1 P2.G],[zeros(P2.dims.p,1) P2.E],P2.id));
else
    error('Invalid Inputs for polynomial zonotope subtraction');
end

end