
% State Transition Tensor Derivations:
r = sym('r',[2,1]);

A1(:,:,1) = r*r'*r(1);
A1(:,:,2) = r*r'*r(2);

A2(:,:,1) = [3*r(1) r(2) r(3);r(2) r(1) 0;r(3) 0 r(1)];
A2(:,:,2) = [r(2) r(1) 0;r(1) 3*r(2) r(3);0 r(3) r(2)];
A2(:,:,3) = [r(3) 0 r(1);0 r(3) r(2);r(1) r(2) 3*r(3)];

dA1(:,:,:,1) = diff(A1,r(1));
dA1(:,:,:,2) = diff(A1,r(2));
dA1(:,:,:,3) = diff(A1,r(3));

dA2(:,:,:,1) = A2*r(1);
dA2(:,:,:,2) = A2*r(2);
dA2(:,:,:,3) = A2*r(3);

AA = dA1 + dA2;

t = [1;2;3];
T = tensorprod(tensorprod(eye(3,3),t',[3],[1]),t',[4],[1]);

TT = permute(T,[1 2 3 4]) + permute(T,[1 3 2 4]) + permute(T,[1 3 4 2]) + permute(T,[3 1 2 4]) + permute(T,[3 1 4 2]) + permute(T,[3 4 1 2]);

TT - subs(AA,{r(1),r(2),r(3)},{t(1),t(2),t(3)})