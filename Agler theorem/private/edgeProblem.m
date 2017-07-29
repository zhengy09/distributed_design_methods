function Edge = edgeProblem(EdgeIndex,Xk,Edge,X,Mc,opts)
% Update the local variables in each edge


NonLocalVariables = length(Edge.clique);
LocalVariables    = sdpvar(NonLocalVariables,1);

%% define the cost function
Cost = 0;
for i = 1:NonLocalVariables
    Clique     = Mc(:,Edge.clique(i));              % the nodes that this clique contains
    lNodeIndexi = sum(Clique(1:EdgeIndex(1)));          % local position of this node in Node.clique(i)
    lNodeIndexj = sum(Clique(1:EdgeIndex(2)));          % local position of this node in Node.clique(i)
    Cost = Cost + (LocalVariables(i) - Xk{Edge.clique(i)}(lNodeIndexi,lNodeIndexj) - 1/opts.mu*Edge.LocalMultipliers(i))^2;
end

%% define the constraints
Constraints = [sum(LocalVariables) - X(EdgeIndex(1),EdgeIndex(2)) == 0];

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
optimize(Constraints,Cost,options);
Edge.LocalVariables  = value(LocalVariables);




end

