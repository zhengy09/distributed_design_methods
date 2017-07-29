function Node = nodeProblem(nodeIndex,Xk,Node,X,Mc,opts)
% Update the local variables in each node

NonLocalVariables = length(Node.clique);
LocalVariables    = sdpvar(NonLocalVariables,1);

%% define the cost function
Cost = 0;
for i = 1:NonLocalVariables
    Clique     = Mc(:,Node.clique(i));              % the nodes that this clique contains
    lNodeIndex = sum(Clique(1:nodeIndex));          % local position of this node in Node.clique(i)
    Cost = Cost + (LocalVariables(i) - Xk{Node.clique(i)}(lNodeIndex,lNodeIndex) - 1/opts.mu*Node.LocalMultipliers(i))^2;
end

%% define the constraints
Constraints = [sum(LocalVariables) - X(nodeIndex,nodeIndex) == 0];

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
optimize(Constraints,Cost,options);
Node.LocalVariables  = value(LocalVariables);




end

