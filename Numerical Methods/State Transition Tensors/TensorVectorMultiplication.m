
function c = TensorVectorMultiplication(A,b)

c = A;
for i = 1:ndims(A)-1
    c = tensorprod(c,b,ndims(c),1);
end

end