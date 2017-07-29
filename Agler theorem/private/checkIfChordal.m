function [y,alpha] = checkIfChordal(R)

%input is a real n by n Adjacency matrix, function checks if it is chordal by first finding an
%elimination ordering using the MCS function and then seeing if it is a
%perfect elimination ordering.

%output is 1 if it is chordal 0 if it is not.

n = size(R,1);
R = abs(R)>0;
R0 = R -diag(diag(R));
R = R0+ diag(ones(n,1));
alpha = maxCardinalitySearch(R);
beta = zeros(n,1);

for i =1:n
    u =find(R0(alpha(i),:)==1); 
    if all(R(u,u)==1)
        beta(i) =1;
        R0(alpha(i),:) =0;
        R0(:,alpha(i)) = 0;
    else
        beta(i) = 0;
    end
end

if sum(beta)==n
   y = 1;
else
   y = 0;
end

