
%% Relaxation Tests:

PZ = PolynomialZonotope(Ellipsoid([1 0;0 1],[0;0]));
T = randn(2,2,2,2);
PZ2 = Compact(TensorMap(PZ,T));

n = [1:8];
for i = 1:length(n)
    PZ2r(i) = GreedyRelax(PZ2,n(i));
    figure(i); hold on;
    plot(PZ2,struct('Splits',15));
    plot(PZ2r(i),struct('Splits',15));
    % plot(Zonotope(PZ2));
    % plot(Zonotope(PZ2r(i)));
    i
end



