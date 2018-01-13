function [Xk,Node] = cliqueProblemDC(CliqueIndex,Xk,Mc,McOverlap,X,Node,Edge,EdgeRep,opts)
%  Subproblem in each clique
% 

Xk = sdpvar(size(Xk,1));    % define a yalmip variable

%% define the cost function
Cost = 0;
gOverlapInd = find(McOverlap == 1);                          % Overlapping nodes in the gloabl index
gCliqueNode = find(Mc == 1);                                 % Clique nodes in the gloabl index
[~,lOverlapInd,~] = intersect(gCliqueNode,gOverlapInd);      % Overlapping nodes in the clique index

% count the local variables in nodes
for k = 1:length(gOverlapInd)
    tmpNode   = Node{gOverlapInd(k)};                               % the node that is included in this clique
    Cost      = Cost + tmpNode.LocalMultipliers*Xk(lOverlapInd(k),lOverlapInd(k));
end

% count the local variables in edges
for i = 1:length(gOverlapInd)
    for j = i+1:length(gOverlapInd)
        tmpEdge   = Edge{gOverlapInd(i),gOverlapInd(j)};                            % the edge that is included in this clique
        if length(tmpEdge.clique) > 1 
            Cost      = Cost + tmpEdge.LocalMultipliers*Xk(lOverlapInd(i),lOverlapInd(j));
        end
    end
end

Cost = Cost + norm(Xk,'fro').^2;

%% define the constraint
Constraints = [Xk >=0];

for i = 1: length(gCliqueNode)
    if length(Node{gCliqueNode(i)}.clique) == 1 &&  Node{gCliqueNode(i)}.clique== CliqueIndex
        Constraints = [Constraints, Xk(i,i) == Node{gCliqueNode(i)}.value];%X(gCliqueNode(i),gCliqueNode(i))];
    end
    for j = i+1:length(gCliqueNode)
        tmpEdge = Edge{gCliqueNode(i),gCliqueNode(j)};
        if length(tmpEdge.clique) == 1 &&  tmpEdge.clique == CliqueIndex
            Constraints = [Constraints, Xk(i,j) == tmpEdge.value];%X(gCliqueNode(i),gCliqueNode(j));];
        end
    end
end

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
optimize(Constraints,Cost,options);
Xk      = value(Xk);

%% Set the solution to each node and edges
% the local variables in nodes
for k = 1:length(gOverlapInd)
    tmpClique = find(Node{gOverlapInd(k)}.clique == CliqueIndex); 
    Node{gOverlapInd(k)}.LocalVariables(tmpClique) = Xk(lOverlapInd(k),lOverlapInd(k));    % the node that is included in this clique
end

% count the local variables in edges
for i = 1:length(gOverlapInd)
    for j = i+1:length(gOverlapInd)
            tmpClique = find(Edge{gOverlapInd(i),gOverlapInd(j)}.clique == CliqueIndex); 
           Edge{gOverlapInd(i),gOverlapInd(j)}.LocalVariables(tmpClique) = Xk(lOverlapInd(i),lOverlapInd(j));
    end
end

end

