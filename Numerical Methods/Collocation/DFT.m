
function [DFTMatrices] = DFT(N)

DFTMatrices.H = 1/sqrt(N)*exp(-2*pi*1i/N*[-(N-1)/2:(N-1)/2]'*[0:N-1]);
DFTMatrices.D = 1/sqrt(N)*1i*exp(2*pi*1i/N*[0:N-1]'*[-(N-1)/2:(N-1)/2])*diag([-(N-1)/2:(N-1)/2]);

end