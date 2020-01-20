% dCdA = dCdA_K1(A,Ca,t,C,x,w) computes the derivative with respect
% to A of the solution C to C'=A*C+[Ca;0] with C(0)=[0;0].
%
% A is a 2x2 matrix.
% Ca is a function handle accepting a vector as argument and returning a
% vector of the same size.
% t is a vector.
% C is the solution at the points t to C'=A*C+[Ca;0] with C(0)=[0;0]
% It is given by the function concentration_K1.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
%
% C is a matrix of size 2 x length(t).
%
% C is given by C(t)=int(expm((t-u)A)*Ca(u),u=0..t).
%
% dCdA is an array of size 4 x 4 x length(t).
% Each 4 x 4 matrix gives the values of the derivative at times t.
%
% V=dC/dA(t).H, that is dC/dA(t) acting on a matrix H, is the solution to
% V'=A*V+H*C with V(0)=[0;0].
% (dC/dA(t)).H=int(exp((t-u)A)*H*C(u),u=0..t)
% That is
% (dC/dA(t)).H=int(kron(C(u).',exp((t-u)A))*H(:),u=0..t)

function dCdA = dCdA_K1(A,Ca,t,C,x,w)

nt=length(t);
dCdA=zeros(2,4,nt);

% Integral between 0 and t(1) to have dCdA(:,:,1).
Cu=@(u)(concentration_K1(A,Ca,0,[0;0],u,x,w));
f=@(u)(kron(Cu(u).',expm((t(1)-u)*A)));
dCdA(:,:,1)=quadglv(f,0,t(1),x,w);

% dCdA.H is the solution to (dCdA.H)'=A*(dCdA.H)+HC with initial condition
% dCdAC(n).H at t(n)
for n=2:nt;
    Cu=@(u)(concentration_K1(A,Ca,t(n-1),C(:,n-1),u,x,w));
    f=@(u)(kron(Cu(u).',expm((t(n)-u)*A)));
    dCdA(:,:,n)=expm((t(n)-t(n-1))*A)*dCdA(:,:,n-1)+quadglv(f,t(n-1),t(n),x,w);
end

end