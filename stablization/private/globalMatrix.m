function Xk = globalMatrix(Xk,Mc,n)
% transform the local matrix Xk into the global form

    NodeIndex = find(Mc == 1);
    E = sparse(length(NodeIndex),n);
    for i = 1:length(NodeIndex)
        E(i,NodeIndex(i)) = 1;
    end
    Xk = E'*Xk*E;
end

