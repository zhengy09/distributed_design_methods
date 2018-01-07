function [Xk, Mc, info] = dcd(X)
% Distributed Computation of Agler's Theorem (chordal decomposition) using ADMM
% find X1,X2,...Xp
%  s.t. \sum Xk = X
%        Xk is positive semidefinte

%% parameters
opts.mu      = 1;
opts.maxIter = 100;
opts.global  = true;
opts.verbose = true;
opts.subbose = false;
opts.eps     = 1e-6;

time = zeros(opts.maxIter,1);

%% finding cliques
n     = size(X,1);            % Dimension
G     = spones(X);            % Sparsity pattern
Ge    = chordalExt(G);        % Chordal extension
Mc    = maximalCliques(Ge);   % Maximal cliques, each column is a clique  
P     = size(Mc,2);           % Number of cliques
SizeC = sum(Mc,1);            % Size of each clique 

%% Overlapping information
Node        = cell(n,1);             % membership of nodes
NodeRep     = sum(Mc,2);             % repetition of nodes in different cliques
NodeOverlap = zeros(n,1);            % Overlapping nodes
NodeOverlap(find(NodeRep>1)) = 1;    %
NodeOverInd = find(NodeRep>1);
for i = 1:n
    Node{i}.clique = find(Mc(i,:) == 1);
end

Edge    = cell(n,n);              % membership of edges
EdgeRep = zeros(n,n);             % repetition of edges in different cliques
for i = 1:n
    for j = i+1:n
         for k = 1:P
             if Mc(i,k) == 1 && Mc(j,k) == 1
                 if isfield(Edge{i,j},'clique')
                    Edge{i,j}.clique = [Edge{i,j}.clique,k];
                    EdgeRep(i,j) = EdgeRep(i,j) + 1;
                 else
                     Edge{i,j}.clique = k;
                     EdgeRep(i,j) = EdgeRep(i,j)+1;
                 end
             end
         end
    end
end

[EdgeRepi, EdgeRepj] = find(EdgeRep > 1);

McOverlap = Mc & repmat(NodeOverlap,1,P); % overlapping nodes in clique k
McNonOver = Mc - McOverlap ;              % Non-overlapping nodes in clique k

%% Initialization 
Xk = cell(P,1);               % variables in cliques
for k = 1:P
    Xk{k} = zeros(SizeC(k)); 
end

for i = 1:n                   % local variables in nodes
    if NodeRep(i) > 1  % set local variables for this node
        Node{i}.LocalVariables   = zeros(NodeRep(i),1);
        Node{i}.LocalMultipliers = zeros(NodeRep(i),1);
    end
end

for i = 1:n                  % local variables in edges
    for j = i+1:n
         if EdgeRep(i,j) > 1
             Edge{i,j}.LocalVariables   = zeros(EdgeRep(i,j),1);
             Edge{i,j}.LocalMultipliers = zeros(EdgeRep(i,j),1);
         end
    end
end

% ----------------------
% Print information
% ----------------------
if opts.verbose == true
    fprintf('-----------------------------------\n')
    fprintf(' Number of Nodes    : %8d\n',n)
    fprintf(' Number of Cliques    %8d\n',P)
    fprintf('-----------------------------------\n')
    fprintf(' Iter.    Error    Time (s)\n')
end

tic
for iter = 1:opts.maxIter
    
    %% Step 1: Each clique solve a subproblem in parallel     
    for k = 1:P
        Xk{k} = cliqueProblem(k,Xk{k},Mc(:,k),McOverlap(:,k),X,Node,Edge,EdgeRep,opts);
    end
    
    %% Step 2: Overlapping nodes and edges solve their own subproblem to update
    % local virables, which can be computed in parallel
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        Node{nodeIndex} = nodeProblem(nodeIndex,Xk,Node{nodeIndex},X,Mc,opts);
    end
    
    for i = 1:length(EdgeRepi)
        %for j = 1:length(EdgeRepj)
            EdgeIndex = [EdgeRepi(i),EdgeRepj(i)];
            Edge{EdgeIndex(1),EdgeIndex(2)} = edgeProblem(EdgeIndex,Xk,Edge{EdgeIndex(1),EdgeIndex(2)},X,Mc,opts);
        %end
    end
    
    %% Step 3: Overlapping nodes and edges update the dual virables
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        Node{nodeIndex} = nodeMultipliers(nodeIndex,Xk,Node{nodeIndex},Mc,opts);
    end
    
    for i = 1:length(EdgeRepi)
        %for j = 1:length(EdgeRepj)
            EdgeIndex = [EdgeRepi(i),EdgeRepj(i)];
            Edge{EdgeIndex(1),EdgeIndex(2)} = edgeMultipliers(EdgeIndex,Xk,Edge{EdgeIndex(1),EdgeIndex(2)},Mc,opts);
        %end
    end
    %% check convergence & output
    [Stop, Info] = ConverCheck(Xk,X,Node,Edge,Mc,EdgeRepi,EdgeRepj,opts);
    time(iter) = toc;
    if opts.verbose == true
        fprintf('%5d   %8.4f   %6.2f\n', iter,Info.presi,time(iter));
    end
    
    if Stop == true
        break;
    end
end

info.iter = iter;
info.time = time(iter);

if opts.global
    tmpX   = zeros(n);
    tmpEig = zeros(P,1); 
    for k = 1:P
        info.Xk{k} = globalMatrix(Xk{k},Mc(:,k),n);
        tmpEig(k)  = min(eig(Xk{k}));
        tmpX       = tmpX + info.Xk{k};
    end
    info.tmpX = tmpX;
end

%% Summaries
if opts.verbose == true
    fprintf('-----------------------------------\n')
    fprintf(' Error norm(sum Xk -X)   :  %8.4e \n',norm(tmpX-X))
    fprintf(' Minimal eigenvalue of Xk:  %8.4f \n',min(tmpEig))
    fprintf(' Minimal eigenvalue of X :  %8.4f \n',min(eig(X)))
    fprintf('-----------------------------------\n')
end

end

