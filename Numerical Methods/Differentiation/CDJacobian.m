function J = CDJacobian(z0,f)

epsilon = 10^-6;
F = f(z0);
J = zeros(length(F),length(z0));

for i = 1:length(z0)
    z1 = z0; z1(i) = z1(i) + epsilon;
    z2 = z0; z2(i) = z2(i) - epsilon;
    J(:,i) = (f(z1)-f(z2))/(2*epsilon);
    disp([num2str(i) ' of ' num2str(length(z0))]);
end

end