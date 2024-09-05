
%% Polynomial Zonotope Tests:

% Test Zonotope:
PZ = PolynomialZonotope([-.5 1 1 1 1;-.5 1 0 -1 1],[0 1 0 1 2;0 0 1 1 0]);
RP = RandomPoints(PZ,10000);
figure; hold on;
plot(PZ);
plot(RP(1,:),RP(2,:),'k.');

% Linear Map Tests:
M = [2 0;0 3];
PZM = LinearMap(PZ,M);
RPM = M*RP;

figure; hold on;
plot(PZM);
plot(RPM(1,:),RPM(2,:),'k.');

% Tensor Map Tests:
M = randn(2,2,2);
for i = 1:length(RP)
    RPM(:,i) = TensorVectorMultiplication(M,RP(:,i));
end
PZM = TensorMap(PZ,M);
RPM2 = RandomPoints(PZM,100000);

figure; hold on;
plot(PZM);
plot(RPM2(1,:),RPM2(2,:),'b.');
plot(RPM(1,:),RPM(2,:),'k.');

% Third Order Tensor:
M = randn(2,2,2,2);
for i = 1:length(RP)
    RPM(:,i) = TensorVectorMultiplication(M,RP(:,i));
end
PZM = TensorMap(PZ,M);
RPM2 = RandomPoints(PZM,100000);

figure; hold on;
plot(RPM(1,:),RPM(2,:),'k.');
plot(RPM2(1,:),RPM2(2,:),'b.');

% Splitting:
PZ = PolynomialZonotope([2 0 1;0 2 1],[1 0 3;0 1 1]);
temp = SplitMax(PZ);
PZ1s = SplitMax(temp(1));
PZ2s = SplitMax(temp(2));
PZ1 = PZ1s(1);
PZ2 = PZ1s(2);
PZ3 = PZ2s(1);
PZ4 = PZ2s(2);

figure; hold on;
plot(PZ);
plot(PZ1);
plot(PZ2);
plot(PZ3);
plot(PZ4);

% Sums
PZ1 = PolynomialZonotope([2 0 1;1 2 1],[1 0 1;0 1 3],[1 2]);
M = [-.5 .2;-.1 .6];
PZ2 = LinearMap(PZ1,M);

figure; hold on;
PZMS = MinkowskiSum(PZ1,PZ2);
PZES = ExactSum(PZ1,PZ2);
plot(PZMS);
plot(PZES);