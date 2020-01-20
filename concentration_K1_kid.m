% C = concentration_K1_kid(A,Ca,ei,t0,C0,t,x,w) computes the solution to
% C'=A*C+Ca*ei with ei=[0;..;1;0;..0] and initial condition C(t0)=C0.
%                           i-mo
% Whereas the problem is such that t0=0 and C0=[0;..;0], the condition
% C(t0)=C0 is useful to in the function dCdA_K1.
%
% A is a 'nc x nc' matrix, with 'nc' the number of model compartments.                        
% Ca is a function handle accepting a vector as argument and returning a vector of the same size.
% t0 is a scalar.
% C0 is a vector of length nc.
% t is a vector.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
%
% C is a matrix of size nc x length(t).
%
% C is given by C(t)=expm((t-t0)A)*C0+int( expm((t-u)A)*Ca(u) ,u=t0..t) [ei]
% on the 'ei' direction.

function C = concentration_K1_kid(A,Ca,ei,t0,C0,t,x,w)

nt=length(t); nc=length(A); ind=find(ei==1);
C=zeros(nc,nt); 

% Integral from t0 to t(1) to have C(:,1)
f=@(u)(Ca(u)*getcols(expm((t(1)-u)*A),ind));
C(:,1)=expm((t(1)-t0)*A)*C0(:)+quadglv(f,t0,t(1),x,w);

% C is the solution to C'=A*C+[Ca;0] with initial condition C(:,n) at t(n)
for n=2:nt;
    f=@(u)(Ca(u)*getcols(expm((t(n)-u)*A),ind));
    C(:,n)=expm((t(n)-t(n-1))*A)*C(:,n-1)+quadglv(f,t(n-1),t(n),x,w);
end

end