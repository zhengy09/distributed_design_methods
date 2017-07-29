function Xk = cliqueProblem(CliqueIndex,Xk,Mc,McOverlap,X,Node,Edge,EdgeRep,opts)
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
    tmpClique = find(tmpNode.clique == CliqueIndex); 
    Cost      = Cost + (Xk(lOverlapInd(k),lOverlapInd(k)) - tmpNode.LocalVariables(tmpClique) ...
                   + 1/opts.mu*tmpNode.LocalMultipliers(tmpClique)).^2;
end
% count the local variables in edges
for i = 1:length(gOverlapInd)
    for j = i+1:length(gOverlapInd)
        tmpEdge   = Edge{gOverlapInd(i),gOverlapInd(j)};                            % the edge that is included in this clique
        if length(tmpEdge.clique) > 1 
            tmpClique = find(tmpEdge.clique == CliqueIndex); 
            Cost      = Cost + (Xk(lOverlapInd(k),lOverlapInd(k)) - tmpEdge.LocalVariables(tmpClique) ...
                       + 1/opts.mu*tmpEdge.LocalMultipliers(tmpClique)).^2;
        end
    end
end

%% define the constraint
McNonOverlap         = Mc - McOverlap;                                % nonoverlapping elements in this clique
gNonOverlapInd       = find(McNonOverlap == 1);                       % nonoverlapping nodes in the gloabl index
[~,lNonOverlapInd,~] = intersect(gCliqueNode,gNonOverlapInd);         % Overlapping nodes in the clique index
Constraints = [Xk >=0];
% for i = 1: length(gNonOverlapInd)
%     for j = i:size(Xk,1)
%         Constraints = [Constraints, Xk(lNonOverlapInd(i),j) == X(gNonOverlapInd(i),gCliqueNode(j))];
%     end
% end

for i = 1: length(gCliqueNode)
    if length(Node{gCliqueNode(i)}.clique) == 1 &&  Node{gCliqueNode(i)}.clique== CliqueIndex
        Constraints = [Constraints, Xk(i,i) == X(gCliqueNode(i),gCliqueNode(i))];
    end
    for j = i+1:length(gCliqueNode)
        tmpEdge = Edge{gCliqueNode(i),gCliqueNode(j)};
        if length(tmpEdge.clique) == 1 &&  tmpEdge.clique == CliqueIndex
            Constraints = [Constraints, Xk(i,j) == X(gCliqueNode(i),gCliqueNode(j))];
        end
    end
end


%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
optimize(Constraints,Cost,options);
Xk      = value(Xk);

end

