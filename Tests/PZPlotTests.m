
%% Polynomial Zonotope Plotting Tests:

% Test the Effect of Splitting:
PZ = PolynomialZonotope(randn(2,10),randi([0 2],10,10));
Splits = 25;
colors = autumn(length(Splits));
for i = 1:length(Splits)
    figure;
    plot(PZ,struct('Splits',Splits(i)));
end

% With large p, plotting struggles. h not important. Order of exponentials
% does have small effect