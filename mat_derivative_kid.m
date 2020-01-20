% M = mat_derivative_kid(alpha,A,K,Ca,ei,t,C,x,w) compute the matrix of the
% derivative of dot(alpha,C) with respect to K and A where C is the
% solution C to C'=A*C+K*Ca*ei with ei=[0;..;1;0;..0] and initial condition C(0)=[0;..;0].
%                                          i-mo
% alpha is a vector.
% A is a 'nc x nc' matrix, with 'nc' the number of model compartments.
% Ca is a function handle accepting a vector as argument and returning a
% vector of the same size.
% t is a vector.
% C is the solution at the points t to C'=A*C+Ca*ei with C(0)=[0;..;0].
% It is given by the function concentration_K1.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
%
% M is a matrix of size length(t) x (nc^2+1).
% Each row gives the value of the partial derivative at a time t.
% The first column is for the derivative with respect to K.
% The last nc^2 columns are for the derivative with respect to A.

function M = mat_derivative_kid(alpha,A,K,Ca,ei,t,C,x,w)

nt=length(t); nc=length(A);
alpha=alpha(:).';

% Derivative with respect to K. Since K -> C is linear, the derivative
% of C with respect to K is simply given by concentration_K1.
dalphaCdK=(alpha*C).';

% Derivative with respect to A, see function dCdA_K1 for details
dalphaCdA=K*reshape(alpha*reshape(dCdA_K1_kid(A,Ca,ei,t,C,x,w),nc,(nc^2)*nt),(nc^2),nt).';

% Concatenate the column vector of length nt and the matrix of size nt x 4 (or 9).
M=[dalphaCdK,dalphaCdA]; 

end