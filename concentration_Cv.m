% Cv = concentration_Cv(A,Ca,t,C,x,w) computes the solution to
% Cv'=-Ksg*Cv+Ksg*Cg with Cv(0)=0, where Cg is the first component of the
% solution to C'=A*C+[Kga*Ca;0] with C(0)=0.
%
% A is a 2x2 matrix.
% Ca is a function handle accepting a vector as argument and returning a
% vector of the same size.
% t is a vector.
% C is the solution at the points t to C'=A*C+[Kga*Ca;0] with C(0)=[0;0]
% It is given by the function concentration.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
%
% Cv is a vector of length length(t).
%
% Cv is given by Cv(t)=Ksg*int(expm(-Ksg(t-u))*Cg(u),u=0..t).

function Cv = concentration_Cv(A,Ca,t,C,x,w)

nt=length(t);
Cv=zeros(1,nt);
Ksg=-(A(1,1)+A(2,1));

% Integral between 0 and t(1) to have Cv(1).
f=@(u)(exp(-Ksg*(t(1)-u))*getrows(concentration_K1(A,Ca,0,[0;0],u,x,w),1));
Cv(1)=quadglv(f,0,t(1),x,w);

% Cv is the solution to Cv'=-Ksg*Cv+Ksg*Cg with initial condition Cv(n) at
% t(n).
for n=2:nt;
    f=@(u)(expm(-Ksg*(t(n)-u))*getrows(concentration_K1(A,Ca,t(n-1),C(:,n-1),u,x,w),1));
    Cv(n)=exp(-Ksg*(t(n)-t(n-1)))*Cv(n-1)+quadglv(f,t(n-1),t(n),x,w);
end

Cv=Ksg*Cv;

end