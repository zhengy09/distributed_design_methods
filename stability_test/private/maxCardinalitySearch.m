function alpha = maxCardinalitySearch(A)
%find an elimination ordering alpha using the maximum cardinality search


initNode = 1;
G     = logical(A); 
d     = size(G, 1); 
order = zeros(1, d); 
numbered    = false(1, d); 
order(1)    = initNode; 
numbered(initNode) = true; 
for i=2:d
  % For each un-numbered node, find the one with the greatest
  % number of numbered neighbors and pick it as next in order
  score = zeros(1,d); 
  U = find(~numbered);% unnumbered verticies
  for j=1:numel(U)
    k = U(j); u = k;
    score(u) = sum((G(k, :) | G(:, k)') & numbered); 
  end
  tmpmax = max(score);
  u = find((score==tmpmax) & ~numbered);
  order(i) = u(1);
  numbered(u(1)) = true;
  
end

alpha = fliplr(order);

