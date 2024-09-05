
function PZ = mtimes(P1,P2)

if isa(P1,'PolynomialZonotope')
    if isscalar(P2)
        PZ = PolynomialZonotope(P2*P1.G,P1.E,P1.id);
    elseif isnumeric(P2)
        error('Incorrect Order for Matrix-Set multiplication');
    else 
        error('Invalid Input for Polynomial Zonotope Multiplication');
    end
else
    if isscalar(P1)
        PZ = PolynomialZonotope(P1*P2.G,P2.E,P2.id);
    elseif isnumeric(P1)
        if ndims(P1) > 2
            PZ = TensorMap(P2,P1);
        else
            PZ = LinearMap(P2,P1);
        end
    else
        error('Invalid Input for Polynomial Zonotope Multiplication');
    end
end

end