
function DC = JacobiConstantPartial(System,x)

if length(x) == 4
    DC = [2*x(1) + ((2*System.mu + 2*x(1))*(System.mu - 1))/((System.mu + x(1))^2 + x(2)^2)^(3/2) - (System.mu*(2*System.mu + 2*x(1) - 2))/((System.mu + x(1) - 1)^2 + x(2)^2)^(3/2), 2*x(2) - (2*x(2)*System.mu)/((System.mu + x(1) - 1)^2 + x(2)^2)^(3/2) + (2*x(2)*(System.mu - 1))/((System.mu + x(1))^2 + x(2)^2)^(3/2), -2*x(3), -2*x(4)];
else
    DC = [2*x(1) + ((2*System.mu + 2*x(1))*(System.mu - 1))/((System.mu + x(1))^2 + x(2)^2 + x(3)^2)^(3/2) - (System.mu*(2*System.mu + 2*x(1) - 2))/((System.mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(3/2), 2*x(2) - (2*x(2)*System.mu)/((System.mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(3/2) + (2*x(2)*(System.mu - 1))/((System.mu + x(1))^2 + x(2)^2 + x(3)^2)^(3/2), (2*x(3)*(System.mu - 1))/((System.mu + x(1))^2 + x(2)^2 + x(3)^2)^(3/2) - (2*System.mu*x(3))/((System.mu + x(1) - 1)^2 + x(2)^2 + x(3)^2)^(3/2), -2*x(4), -2*x(5), -2*x(6)];
end

end