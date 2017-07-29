function [Areturn] = chordalExt(A)

%given a negative definite matrix find its chordal extension by
%first finding a minimum degree ordering of the nodes and 
%then performing a symbolic Cholesky factorisation.

n = size(A,1);
A = abs(A)>0;

B = 0.1*rand(n,n);
tmp = B'*B;
tmp = tmp.*A;
tmp = tmp + max(eig(tmp))*eye(n);

order = symamd(tmp);
Aperm = tmp(order,order);
L = chol(Aperm);
Aext = L + L';

u = zeros(n,1);
for i =1:n
    u(i) = find(order==i);
end

Areturn = abs(Aext(u,u))>0;




