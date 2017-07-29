function  Edge = edgeMultipliers(EdgeIndex,Xk,Edge,Mc,opts)
%   Detailed explanation goes here

    for i = 1:length(Edge.clique)
        Clique     = Mc(:,Edge.clique(i));                  % the nodes that this clique contains
        lNodeIndexi = sum(Clique(1:EdgeIndex(1)));          % local position of this node in Node.clique(i)
        lNodeIndexj = sum(Clique(1:EdgeIndex(2)));          % local position of this node in Node.clique(i)
        Edge.LocalMultipliers(i)  = Edge.LocalMultipliers(i) + ... 
                        opts.mu*(Xk{Edge.clique(i)}(lNodeIndexi,lNodeIndexj) - Edge.LocalVariables(i));
    end

end

