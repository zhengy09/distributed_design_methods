function [Node,Edge,Clique] = cliqueProblem(CliqueIndex,Clique,Node,Edge,A,opts,iter)
%  Subproblem in each clique
%  Update the variables that belong to each clique

epsilon = 1e-2;

%% Define virables
Nk = length(Clique.node);                % how many nodes are there in this clique?
P  = cell(Nk,1);                         % Define sdp variables Pi
Xk = sdpvar(Clique.sdpsize);
for i = 1:Nk
   P{i} = sdpvar(Clique.nodsize(i));     % some are global, and others are local to each clique
end

%% define the cost function
Cost = 0;

% count the local variables in nodes
nodeC  = Clique.node;                             % Nodes in this clique
nodeCR = Clique.nodeR;                            % Overlapping nodes in this clique
accDim = cumsum([1;Clique.nodsize]);              % accumulative dimension Index
[~,lOverlapInd,~] = intersect(nodeC,nodeCR);      % Overlapping nodes in the clique index

for i = 1:length(nodeCR)  
    tmpNode   = Node{nodeCR(i)};
    tmpClique = find(tmpNode.clique == CliqueIndex); 
    lIndex    = accDim(lOverlapInd(i)):accDim(lOverlapInd(i)+1)-1;
    Cost      = Cost + norm(Xk(lIndex,lIndex) - tmpNode.LocalVariables{tmpClique} ...
                   + 1/opts.mu*tmpNode.LocalMultipliers{tmpClique},'fro').^2;
    gIndex    = find(nodeC == nodeCR(i));
    Cost      = Cost + norm(P{gIndex} - tmpNode.P ...
                   + 1/opts.mu*tmpNode.PiMultipliers{tmpClique},'fro').^2;
end

% count the local variables in edges
for i = 1:length(nodeCR)
    for j = i+1:length(nodeCR)
        tmpEdge   = Edge{nodeCR(i),nodeCR(j)};              % the edge that is included in this clique
        if length(tmpEdge.clique) > 1 
            tmpClique = find(tmpEdge.clique == CliqueIndex); 
            lIndexi    = accDim(lOverlapInd(i)):accDim(lOverlapInd(i)+1)-1;
            lIndexj    = accDim(lOverlapInd(j)):accDim(lOverlapInd(j)+1)-1;
            Cost      = Cost + norm(Xk(lIndexi,lIndexj) - tmpEdge.LocalVariables{tmpClique} ...
                       + 1/opts.mu*tmpEdge.LocalMultipliers{tmpClique},'fro').^2;
        end
    end
end

%% define the constraint
Constraints = [];
for i = 1: length(nodeC)
    if length(Node{nodeC(i)}.clique) == 1 && Node{nodeC(i)}.clique == CliqueIndex
        lIndex    = accDim(i):accDim(i+1)-1;
        Constraints = [Constraints, Xk(lIndex,lIndex) == - A{nodeC(i),nodeC(i)}'*P{i} ...
                                        - P{i}*A{nodeC(i),nodeC(i)} -epsilon*eye(Clique.nodsize(i))];
        Constraints = [Constraints, P{i}-epsilon*eye(Clique.nodsize(i)) >=0];
    end
    for j = i+1:length(nodeC)
        tmpEdge = Edge{nodeC(i),nodeC(j)};
        if length(tmpEdge.clique) == 1 &&  tmpEdge.clique == CliqueIndex
            lIndexi    = accDim(i):accDim(i+1)-1;
            lIndexj    = accDim(j):accDim(j+1)-1;
            Constraints = [Constraints, Xk(lIndexi,lIndexj) == - A{nodeC(j),nodeC(i)}'*P{j} ...
                                        - P{i}*A{nodeC(i),nodeC(j)}];
        end
    end
end

Constraints = [Constraints,Xk >=0];

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
sol     = optimize(Constraints,Cost,options);
Xk      = value(Xk);

Clique.time(iter) = sol.solvertime; 

%% set value to the correspond nodes and edges
for i = 1:length(nodeC)
    if length(Node{nodeC(i)}.clique) == 1 && Node{nodeC(i)}.clique == CliqueIndex
        Node{nodeC(i)}.P = value(P{i});
    else
        tmpClique = find(Node{nodeC(i)}.clique == CliqueIndex); 
        lIndex    = accDim(i):accDim(i+1)-1;
        Node{nodeC(i)}.CliqueVariables{tmpClique} = Xk(lIndex,lIndex);
        
        Node{nodeC(i)}.Pi{tmpClique} = value(P{i});
    end
    for j = i+1:length(nodeC)
        tmpEdge   = Edge{nodeC(i),nodeC(j)};              % the edge that is included in this clique
        if length(tmpEdge.clique) > 1 
            tmpClique = find(tmpEdge.clique == CliqueIndex); 
            lIndexi    = accDim(i):accDim(i+1)-1;
            lIndexj    = accDim(j):accDim(j+1)-1;
            Edge{nodeC(i),nodeC(j)}.CliqueVariables{tmpClique} = Xk(lIndexi,lIndexj);
        end
    end
end

end

