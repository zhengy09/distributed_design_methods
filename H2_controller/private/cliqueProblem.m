function [Node,Edge,Clique,sol] = cliqueProblem(CliqueIndex,Clique,Node,Edge,A,opts,iter)
%  Subproblem in each clique
%  Update the variables that belong to each clique

epsilon = 1e-2;

%% Define virables
Nk = length(Clique.node);                % how many nodes are there in this clique?
X  = cell(Nk,1);                         % Define sdp variables Xi
Xk = sdpvar(Clique.sdpsize);
for i = 1:Nk
   X{i} = sdpvar(Clique.nodsize(i));     % some are local to each node (shared by other clique), and others are local to each clique
end

%% define the cost function
Cost = 0;

% count the local variables in nodes
nodeC  = Clique.node;                             % Nodes in this clique
nodeCR = Clique.nodeR;                            % Overlapping nodes in this clique
accDim = cumsum([1;Clique.nodsize]);              % accumulative dimension Index
[~,lOverlapInd,~] = intersect(nodeC,nodeCR);      % Overlapping nodes in the clique index

for i = 1:length(nodeC)
    if length(Node{nodeC(i)}.clique) == 1 %&& Node{nodeC(i)}.clique == CliqueIndex
        [ti,tj] = size(Node{nodeC(i)}.Y);
        Node{nodeC(i)}.Yi = sdpvar(ti,tj,'full');
        
        [ti,tj] = size(Node{nodeC(i)}.Z);
        Node{nodeC(i)}.Zi = sdpvar(ti,tj,'full');     % for computation
        Cost = Cost + trace(Node{nodeC(i)}.Q*X{i})+trace(Node{nodeC(i)}.R*Node{nodeC(i)}.Yi);
    end
end

for i = 1:length(nodeCR)  
    tmpNode   = Node{nodeCR(i)};
    tmpClique = find(tmpNode.clique == CliqueIndex); 
    lIndex    = accDim(lOverlapInd(i)):accDim(lOverlapInd(i)+1)-1;
    Cost      = Cost + opts.mu/2*norm(Xk(lIndex,lIndex) - tmpNode.LocalVariables{tmpClique} ...
                   + tmpNode.LocalMultipliers{tmpClique},'fro').^2;
    gIndex    = find(nodeC == nodeCR(i));
    Cost      = Cost + opts.mu/2*norm(X{gIndex} - tmpNode.X ...
                   + tmpNode.XiMultipliers{tmpClique},'fro').^2;
end

% count the local variables in edges
for i = 1:length(nodeCR)
    for j = i+1:length(nodeCR)
        tmpEdge   = Edge{nodeCR(i),nodeCR(j)};              % the edge that is included in this clique
        if length(tmpEdge.clique) > 1 
            tmpClique = find(tmpEdge.clique == CliqueIndex); 
            lIndexi    = accDim(lOverlapInd(i)):accDim(lOverlapInd(i)+1)-1;
            lIndexj    = accDim(lOverlapInd(j)):accDim(lOverlapInd(j)+1)-1;
            Cost      = Cost + opts.mu/2*norm(Xk(lIndexi,lIndexj) - tmpEdge.LocalVariables{tmpClique} ...
                       + tmpEdge.LocalMultipliers{tmpClique},'fro').^2;
                   
           % added for edges
           Cost       = Cost + opts.mu/2*norm(X{find(Clique.node == nodeCR(i))} - tmpEdge.Xi ...
                       + tmpEdge.XiMultipliers{tmpClique},'fro').^2;
           Cost       = Cost + opts.mu/2*norm(X{find(Clique.node == nodeCR(j))} - tmpEdge.Xj ...
                       + tmpEdge.XjMultipliers{tmpClique},'fro').^2;
        end
    end
end

%% define the constraint
Constraints = [];
for i = 1: length(nodeC)
    if length(Node{nodeC(i)}.clique) == 1 && Node{nodeC(i)}.clique == CliqueIndex
        lIndex    = accDim(i):accDim(i+1)-1;
        Constraints = [Constraints, Xk(lIndex,lIndex) == - (A{nodeC(i),nodeC(i)}*X{i} - Node{nodeC(i)}.B*Node{nodeC(i)}.Zi) ...
                                        - (A{nodeC(i),nodeC(i)}*X{i} - Node{nodeC(i)}.B*Node{nodeC(i)}.Zi)' - Node{nodeC(i)}.M*Node{nodeC(i)}.M'];
        Constraints = [Constraints, X{i}-epsilon*eye(Clique.nodsize(i)) >=0];
        Constraints = [Constraints, [Node{nodeC(i)}.Yi Node{nodeC(i)}.Zi;Node{nodeC(i)}.Zi' X{i}] >=0];
    end
    for j = i+1:length(nodeC)
        tmpEdge = Edge{nodeC(i),nodeC(j)};
        if length(tmpEdge.clique) == 1 &&  tmpEdge.clique == CliqueIndex
            lIndexi    = accDim(i):accDim(i+1)-1;
            lIndexj    = accDim(j):accDim(j+1)-1;
            Constraints = [Constraints, Xk(lIndexi,lIndexj) == - X{i}*A{nodeC(j),nodeC(i)}' ...
                                        - A{nodeC(i),nodeC(j)}*X{j}];
        end
    end
end

Constraints = [Constraints,Xk >=0];

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi','cachesolvers',1);
sol     = optimize(Constraints,Cost,options);
Xk      = value(Xk);

Clique.time(iter) = sol.solvertime; 

%% set value to the correspond nodes and edges
for i = 1:length(nodeC)
    if length(Node{nodeC(i)}.clique) == 1 && Node{nodeC(i)}.clique == CliqueIndex
        Node{nodeC(i)}.X = value(X{i});
        Node{nodeC(i)}.Y = value(Node{nodeC(i)}.Yi);
        Node{nodeC(i)}.Z = value(Node{nodeC(i)}.Zi);
    else
        tmpClique = find(Node{nodeC(i)}.clique == CliqueIndex); 
        lIndex    = accDim(i):accDim(i+1)-1;
        Node{nodeC(i)}.CliqueVariables{tmpClique} = Xk(lIndex,lIndex);
        
        Node{nodeC(i)}.Xi{tmpClique} = value(X{i});
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

