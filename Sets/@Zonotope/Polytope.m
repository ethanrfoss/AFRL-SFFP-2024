
function PT = Polytope(Z)

alpha = (dec2bin(0:2^(Z.dims.l)-1)-'0')*2-1;
PT = Polytope((Z.c + Z.G*alpha')');

end