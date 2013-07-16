function X = mypinv(A,my_tol)
%PINV   Pseudoinverse.
%   X = myPINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS(class(A)).
%
%   myPINV(A,TOL) uses the tolerance TOL instead of the default.

if ~exist('my_tol','var')
   my_tol=eps*1024;
end

if isempty(A)     % quick return
  X = zeros(size(A'),class(A));  
  return  
end

% m = number of rows of A
% n = number of columns of A
[m,n] = size(A);

% if columns is parger than rows
if n > m
    % Apply the same algorithm to A^T and traspose the result
   X = pinv(A',my_tol)';
else
    % Obtain singular values of A
    %   [U,S,V] = SVD(X) produces a diagonal matrix S, of the same 
    %   dimension as X and with nonnegative diagonal elements in
    %   decreasing order, and unitary matrices U and V so that
    %   X = U*S*V'.
    %
    %   [U,S,V] = SVD(X,0) produces the "economy size"
    %   decomposition. If X is m-by-n with m > n, then only the
    %   first n columns of U are computed and S is n-by-n.
    %   For m <= n, SVD(X,0) is equivalent to SVD(X).
        [U,S,V] = svd(A,0);
   
   %v = diag(X) returns the main diagonal of X as a vector
   if m > 1 
       s = diag(S);
   elseif m == 1
       s = S(1);
   else
       s = 0;
   end
   
   if my_tol > 0
      tol = my_tol;
   else
       %d = eps(X) is the positive distance from abs(X) to the next larger
       %in magnitude floating point number of the same precision as X. 
       %X may be either double precision or single precision. For all X, 
      tol = max(m,n) * eps(max(s));
   end
   %B = sum(A) returns sums along different dimensions of an array.
   r = sum(s > tol);
   if (r == 0)
      X = zeros(size(A'),class(A));
   else
       %s(1:r) because they come in order
       %s is overwriten and in the end it has size r x r
      s = diag(ones(r,1)./s(1:r));
      %(:,1:r) because a is r x r
      X = V(:,1:r)*s*U(:,1:r)';
   end
end