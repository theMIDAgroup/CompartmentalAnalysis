% M = mat_derivative(alpha,A,K,Ca,t,C,x,w) compute the matrix of the
% derivative of dot(alpha,Cc) with respect to K and A where Cc is the
% solution Cc to Cc'=A*Cc+[K*Ca;0] with Cc(0)=[0;0] and alpha a vector.
%
% alpha is a vector of length 2.
% A is a 2x2 matrix.
% Ca is a function handle accepting a vector as argument and returning a
% vector of the same size.
% t is a vector.
% C=Cc/K is the solution at the points t to C'=A*C+[Ca;0] with C(0)=[0;0]
% It is given by the function concentration_K1.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
%
% M is a matrix of size 5 x length(t).
% Each row gives the value of the partial derivative a time t.
% The first column is for the derivative with respect to K.
% The last 4 columns are for the derivative with respect to A.
%
% If opt is present, it has to be 'M', 'A' or 'K'.
% opt='M' does not change anything.
% opt='A' only computes the derivative with respect to A and put zeros
% elsewhere.
% opt='K' only computes the derivative with respect to K and put zeros
% elsewhere.

function M = mat_derivative(alpha,A,K,Ca,t,C,x,w,opt)

if nargin==8
    opt='M';
end

nt=length(t);
alpha=alpha(:).';

% Derivative with respect to K. Since K -> C is linear, the derivative of C
% with respect to K is simply given by concentration_K1.
if opt=='M'|| opt=='K'
    dalphaCdK=(alpha*C).';
else
    dalphaCdK=zeros(size(C,2),1);
end

% Derivative with respect to A, see function dCdA_K1 for details
if opt=='M'|| opt=='A'
    dalphaCdA=K*reshape(alpha*reshape(dCdA_K1(A,Ca,t,C,x,w),2,4*nt),4,nt).';
else
    dalphaCdA=zeros(nt,4);
end

% Concatenate the column vector of length nt and the matrix of size 4 x nt.
M=[dalphaCdK,dalphaCdA];

end