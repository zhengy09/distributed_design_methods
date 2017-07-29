function [P, info] = dst(A,G)
% Distributed Stability Test using chordal decomposition and ADMM
% find P1,P2,...Pn
%  s.t. A^TP + PA < 0
%       P = diag(P1, P2, .., Pn) is positive definte
%
% Input: graph G  -- the dynamical interconnection graph, a matrix
%        matrix A -- cell format

%% parameters
opts.mu      = 1;
opts.maxIter = 10;
opts.verbose = true;          % display output information
opts.subbose = false;          % display output information for each subproblem
opts.eps     = 1e-4;
opts.global  = true;

time = zeros(opts.maxIter,1);

%% finding cliques
n     = size(G,1);            % Dimension number of nodes
G     = spones(G);            % Sparsity pattern
Ge    = chordalExt(G);        % Chordal extension
Mc    = maximalCliques(Ge);   % Maximal cliques, each column is a clique  
Nc    = size(Mc,2);           % Number of cliques

%% process the model data such that Ge and A are consistent
subDimension = zeros(n,1);
for i = 1:n
    subDimension(i) = size(A{i,i},1);
end

for i = 1:n
    for j = 1:n
        if Ge(i,j) == 1 && isempty(A{i,j})
            A{i,j} = zeros(subDimension(i),subDimension(j));
        end
    end
end

%% Overlapping information
tic;
Node        = cell(n,1);             % Nodes
NodeRep     = sum(Mc,2);             % repetition of nodes in different cliques
NodeOverInd = find(NodeRep > 1);     % nodes that belong to more than clique

% Node information
for i = 1:n
    Node{i}.Id = i;
    Node{i}.clique = find(Mc(i,:) == 1);   % the cliques that contain this node
end

% Clique information
Clique = cell(Nc,1);                 
for k = 1:Nc
    Clique{k}.node    = find(Mc(:,k) == 1);                     % the nodes in this clique
    Clique{k}.nodsize = subDimension(Clique{k}.node);           % the size of local dynamics in each node
    Clique{k}.sdpsize = sum(Clique{k}.nodsize);                 % the dimension of sdp constraint
    Clique{k}.nodeR   = intersect(Clique{k}.node,NodeOverInd);  % Repeated nodes in this clique
    Clique{k}.edgeR   = [];                                     % Repeated edges in this clique
    Clique{k}.time    = zeros(opts.maxIter,1);                  % Time consumption for solving this clique problem
end

Edge    = cell(n,n);              % Edges
EdgeRep = zeros(n,n);             % repetition of edges in different cliques
for i = 1:n
    for j = i+1:n
         for k = 1:Nc
             if Mc(i,k) == 1 && Mc(j,k) == 1
                 Edge{i,j}.Id = [i,j];
                 if isfield(Edge{i,j},'clique')
                    Edge{i,j}.clique = [Edge{i,j}.clique,k]; % the cliques that contain this edge
                    EdgeRep(i,j) = EdgeRep(i,j) + 1;         
                 else
                     Edge{i,j}.clique = k;
                     EdgeRep(i,j) = EdgeRep(i,j) + 1;
                 end
             end
         end
    end
end

[EdgeRepi, EdgeRepj] = find(EdgeRep > 1);  % edges that belong to more than one clique

timeChordal = toc;

%% Initialization
for i = 1:n                   % local variables in nodes
    Node{i}.P  = zeros(size(A{i,i}));    % Lyapunov variables
    if NodeRep(i) > 1                    % set local variables for this node
        for k = 1:NodeRep(i)
            % for Mi
            Node{i}.CliqueVariables{k}  = zeros(subDimension(i));   % this variable belong to clique in node.clique(k)
            Node{i}.LocalVariables{k}   = zeros(subDimension(i));   % Local variable in this node
            Node{i}.LocalMultipliers{k} = zeros(subDimension(i));   % Corresponding multipliers
            
            % for Pi
            Node{i}.Pi{k}   = zeros(size(A{i,i}));               % This Pi belongs to node.clique(k) 
            Node{i}.PiMultipliers{k} = zeros(size(A{i,i}));      % Corresponding multipliers
            
            % for time record
            Node{i}.time   = zeros(opts.maxIter,2);
        end
    end
end

for i = 1:n                  % local variables in edges
    for j = i+1:n
         if EdgeRep(i,j) > 1
             for k = 1:EdgeRep(i,j)
                 % for Mi
                 Edge{i,j}.CliqueVariables{k}  = zeros(subDimension(i),subDimension(j));
                 Edge{i,j}.LocalVariables{k}   = zeros(subDimension(i),subDimension(j));
                 Edge{i,j}.LocalMultipliers{k} = zeros(subDimension(i),subDimension(j));
                 
                 % for Pi and Pj
                 Edge{i,j}.Pi = zeros(subDimension(i));   % this is equal to node{i}.P
                 Edge{i,j}.Pj = zeros(subDimension(j));     
                 
                 % for time record
                 Edge{i}.time = zeros(opts.maxIter,2);
             end
         end
    end
end

% ----------------------
% Print information
% ----------------------
if opts.verbose == true
    fprintf('-------------------------------------\n')
    fprintf(' Number of Nodes    :     %8d\n',n)
    fprintf(' Number of Cliques        %8d\n',Nc)
    fprintf(' Max. Clique Size         %8d\n',max(sum(Mc,1)))
    fprintf(' Min. Clique Size         %8d\n',min(sum(Mc,1)))
    fprintf(' Chordal processsing (s): %8.4f\n', timeChordal)
    fprintf('-------------------------------------\n')
    fprintf(' Iter.    Error    Time (s)\n')
end

tic
for iter = 1:opts.maxIter
    
    %% Step 1: Each clique solve a subproblem in parallel     
    for k = 1:Nc
        [Node,Edge,Clique{k}] = cliqueProblem(k,Clique{k},Node,Edge,A,opts,iter);
    end
    
    %% Step 2: Overlapping nodes and edges solve their own subproblem to update
    % local virables, which can be computed in parallel
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        Node{nodeIndex} = nodeProblem(Node{nodeIndex},A{nodeIndex,nodeIndex},opts,iter);
    end
    
    for k = 1:length(EdgeRepi)
        %for j = 1:length(EdgeRepj)
            Eind = [EdgeRepi(k),EdgeRepj(k)];
            A12  = A{Eind(1),Eind(2)};
            A21  = A{Eind(2),Eind(1)};
            P1   = Node{Eind(1)}.P;
            P2   = Node{Eind(2)}.P;
            Edge{Eind(1),Eind(2)} = edgeProblem(Edge{Eind(1),Eind(2)},A12,A21,P1,P2,opts,iter);
        %end
    end
    
    %% Step 3: Overlapping nodes and edges update the dual virables
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        Node{nodeIndex} = nodeMultipliers(Node{nodeIndex},opts,iter);
    end
    
    for k = 1:length(EdgeRepi)
        %for j = 1:length(EdgeRepj)
            Eind = [EdgeRepi(k),EdgeRepj(k)];
            Edge{Eind(1),Eind(2)} = edgeMultipliers(Edge{Eind(1),Eind(2)},opts,iter);
        %end
    end
    %% check convergence & output
    [Stop, Info] = ConverCheck(Node,Edge,NodeOverInd,EdgeRepi,EdgeRepj,opts);
    time(iter) = toc;
    if opts.verbose == true
        fprintf('%5d   %8.4f   %6.2f\n', iter,Info.presi,time(iter));
    end
    
    if Stop == true
        break;
    end
end

%% time post processing
timeClique   = zeros(iter,Nc); 
timesubtotal = 0;
timeparal    = 0;
for k = 1:Nc
    timesubtotal    = timesubtotal + sum(Clique{k}.time); 
    timeClique(:,k) = Clique{k}.time(1:iter); 
end
timeparal  = timeparal + sum(max(timeClique,[],2));

timeNode   = zeros(iter,length(NodeOverInd)*2); 
for i = 1:length(NodeOverInd)
    timesubtotal               = timesubtotal + sum(sum(Node{NodeOverInd(i)}.time)); 
    timeNode(:,(i-1)*2+1:i*2)  = Node{NodeOverInd(i)}.time(1:iter,:); 
end
timeparal  = timeparal + sum(max(timeNode,[],2));

timeEdge   = zeros(iter,length(EdgeRepi)*2); 
for i = 1:length(EdgeRepi)
    timesubtotal               = timesubtotal + sum(sum(Edge{EdgeRepi(i),EdgeRepj(i)}.time)); 
    timeEdge(:,(i-1)*2+1:i*2)  = Edge{EdgeRepi(i),EdgeRepj(i)}.time(1:iter,:); 
end
if ~isempty(EdgeRepi)
    timeparal  = timeparal + sum(max(timeEdge,[],2));
end

info.Mc = Mc;
info.timeTotal    = [time(iter),timesubtotal,timeparal];  % Total time; time for SeDuMi only; estimated time for parallel computation?
info.time.chordal = timeChordal;
info.time.clique  = timeClique; 
info.time.node    = timeNode; 
info.time.edge    = timeEdge; 
%%
P       = cell(n,1);
tmpEig  = zeros(n,1);
for i = 1:n
    P{i} = Node{i}.P;
    tmpEig(i) = min(eig(P{i}));
end

if opts.global
    accDimen  = [cumsum([1;subDimension])]; 
    globalA   = zeros(sum(subDimension));  % gloabl state space model
    globalP   = zeros(sum(subDimension)); 
    for i = 1:n
        globalP(accDimen(i):accDimen(i+1)-1,accDimen(i):accDimen(i+1)-1) = P{i};
        for j = 1:n
            if G(i,j) ~= 0
                globalA(accDimen(i):accDimen(i+1)-1,accDimen(j):accDimen(j+1)-1) = A{i,j};
            end
        end
    end
    info.globalA = globalA;
    info.globalP = globalP;
    
    % test the result
    flag = min(globalA'*globalP + globalP*globalA);
    info.minEig = max(flag);
    if info.minEig  < 0
        info.flag = true;
    else
        info.flag = false;
    end
end

%% Summaries
if opts.verbose == true
    fprintf('--------------------------------------\n')
    fprintf(' Residual (Pi - P):      :   %8.2e \n',Info.presi)
    fprintf(' Minimal eigenvalue of Pk:   %8.2e \n',min(tmpEig))
    fprintf('--------------------------------------\n')
end

end

