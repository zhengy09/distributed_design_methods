function [Node,Edge,Clique] = cliqueProblem(CliqueIndex,Clique,Node,Edge,A,B,Gc,opts,iter)
%  Subproblem in each clique
%  Update the variables that belong to each clique

epsilon = 1e-2;

% count the local variables in nodes
nodeC  = Clique.node;                             % Nodes in this clique
nodeCR = Clique.nodeR;                            % Overlapping nodes in this clique
accDim = cumsum([1;Clique.nodsize]);              % accumulative dimension Index
[~,lOverlapInd,~] = intersect(nodeC,nodeCR);      % Overlapping nodes in the clique index

%% Define virables
Nk = length(Clique.node);                % how many nodes are there in this clique?
X  = cell(Nk,1);                         % Define sdp variables Xi
Z  = cell(Nk,Nk);                        % Define variables Zij -- feedback gains Kij = Zij Xi^{-1}
Xk = sdpvar(Clique.sdpsize);
for i = 1:Nk
   X{i} = sdpvar(Clique.nodsize(i));     % some are global, and others are local to each clique
   for j = 1:Nk  % some are not necessary 
       if i == j || Gc(nodeC(i),nodeC(j)) == 1
            Z{i,j} = sdpvar(Clique.nodinpt(i),Clique.nodsize(j));   % belongs to node & edge
       end
   end
end

%% define the cost function
Cost = 0;

for i = 1:length(nodeCR)  
    tmpNode   = Node{nodeCR(i)};
    tmpClique = find(tmpNode.clique == CliqueIndex); 
    lIndex    = accDim(lOverlapInd(i)):accDim(lOverlapInd(i)+1)-1;
    Cost      = Cost + norm(Xk(lIndex,lIndex) - tmpNode.LocalVariables{tmpClique} ...
                   + 1/opts.mu*tmpNode.LocalMultipliers{tmpClique},'fro').^2;
    gIndex    = find(nodeC == nodeCR(i));
    Cost      = Cost + norm(X{gIndex} - tmpNode.X ...
                   + 1/opts.mu*tmpNode.XiMultipliers{tmpClique},'fro').^2;
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
        Constraints = [Constraints, Xk(lIndex,lIndex) == - A{nodeC(i),nodeC(i)}*X{i} ...
                                        - X{i}*A{nodeC(i),nodeC(i)}' - B{nodeC(i)}*Z{i,i} ...
                                        - Z{i,i}'*B{nodeC(i)}' - epsilon*eye(Clique.nodsize(i))];
        Constraints = [Constraints, X{i}-epsilon*eye(Clique.nodsize(i)) >=0];
    end
    for j = i+1:length(nodeC)
        tmpEdge = Edge{nodeC(i),nodeC(j)};
        if length(tmpEdge.clique) == 1 &&  tmpEdge.clique == CliqueIndex
            lIndexi    = accDim(i):accDim(i+1)-1;
            lIndexj    = accDim(j):accDim(j+1)-1;
            if Gc(nodeC(i),nodeC(j)) == 0 && Gc(nodeC(j),nodeC(i)) == 0           % case of no communications (ij) and (ji)
                Constraints = [Constraints, Xk(lIndexi,lIndexj) == - X{i}*A{nodeC(j),nodeC(i)}' ...
                                        - A{nodeC(i),nodeC(j)}*X{j}];
            elseif Gc(nodeC(i),nodeC(j)) == 1 && Gc(nodeC(j),nodeC(i)) == 0        % communication (j -> i) but no (i -> j)
                Constraints = [Constraints, Xk(lIndexi,lIndexj) == - X{i}*A{nodeC(j),nodeC(i)}' ...
                                        - A{nodeC(i),nodeC(j)}*X{j} - B{nodeC(i)}*Z{i,j}];
            elseif Gc(nodeC(i),nodeC(j)) == 0 && Gc(nodeC(j),nodeC(i)) == 1        % communication (i -> j) but no (j -> i)
                Constraints = [Constraints, Xk(lIndexi,lIndexj) == - X{i}*A{nodeC(j),nodeC(i)}' ...
                                        - A{nodeC(i),nodeC(j)}*X{j} - Z{j,i}'*B{nodeC(j)}'];     
            else
                Constraints = [Constraints, Xk(lIndexi,lIndexj) == - X{i}*A{nodeC(j),nodeC(i)}' ...
                     - A{nodeC(i),nodeC(j)}*X{j} - B{nodeC(i)}*Z{i,j} - Z{j,i}'*B{nodeC(j)}'];
            end
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
        Node{nodeC(i)}.X = value(X{i});
        Node{nodeC(i)}.Z = value(Z{i,i});
    else
        tmpClique = find(Node{nodeC(i)}.clique == CliqueIndex); 
        lIndex    = accDim(i):accDim(i+1)-1;
        Node{nodeC(i)}.CliqueVariables{tmpClique} = Xk(lIndex,lIndex);
        
        Node{nodeC(i)}.Xi{tmpClique} = value(X{i});
    end
    
    for j = 1:length(nodeC)
        if Gc(nodeC(i),nodeC(j)) == 1 && i~=j
            Edge{nodeC(i),nodeC(j)}.Z = value(Z{i,j});
        end
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

