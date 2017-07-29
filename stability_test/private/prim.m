function [MCO,Adj,MC] = prim(A)

%Section 2.1.3 

% an implementation of Prim's algorithm

% the function  returns the maximal cliques of the graph arranged in
% topological order so thast the first column of MCO is C1 ....the mth
% column of MC if Cm

% have noticed that this function behaves incorrectly when the graph is not
% a single connected component - need to be careful with this.

[m,n] = size(A);

if m==n
    
    MC = maximalCliques(A);
    
    [p,q] = size(MC);
    
    T = zeros(q,q);
    
    for i =1:q;
        for j =i:q
            T(i,j) = sum(and(MC(:,i),MC(:,j)));
        end
    end
    
    T=T- diag(diag(T));
    T= T+T';
    
    K = zeros(q,1);
    K(1) = true;
    Adj = zeros(q,q);
    
    for i =1:q-1
        u = find(K);
        v = find(~K);
        [a,b]  =find(T(u,v)==max(max(T(u,v))));
        Adj(u(a(1)),v(b(1))) =1;
        Adj(v(b(1)),u(a(1))) =1;
        K(v(b(1))) =true;
    end
    
    numbered = zeros(q,1);
    numbered(1) =true;
    toporder = zeros(q,1);
    toporder(1,1) = q;
    
    for i =q-1:-1:1
        u = find(numbered);
        v = find(~numbered);
        [a,b] = find(Adj(u,v)==1);
        numbered(v(b(1))) =true;
        toporder(v(b(1))) = i;
    end
    
    MCO = zeros(p,q);
    for i =1:q
        MCO(:,toporder(i)) =MC(:,i);
    end  
  
    
else
    error('myApp:prim', 'Matrix is not square')
    
end

