
function V0 = StateTransitionTensorInitialCondition(x0,nx,nk,N)

if N >= 2
    V0 = zeros(sum([nx nx^2 nx*(nk).^(2:N-1)]),1);
    V0(nx+1:nx+nx^2) = reshape(eye(nx),[nx^2,1]);
    V0(1:nx) = x0;
else
    V0 = x0;
end

end