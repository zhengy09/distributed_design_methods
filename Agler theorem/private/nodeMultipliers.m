function Node = nodeMultipliers(nodeIndex,Xk,Node,Mc,opts)
% Update the multipliers in each node

    for i = 1:length(Node.clique)
        Clique     = Mc(:,Node.clique(i));              % the nodes that this clique contains
        lNodeIndex = sum(Clique(1:nodeIndex));          % local position of this node in Node.clique(i)
        Node.LocalMultipliers(i)  = Node.LocalMultipliers(i) ...
                                + opts.mu*(Xk{Node.clique(i)}(lNodeIndex,lNodeIndex) - Node.LocalVariables(i));
    end
    
end

