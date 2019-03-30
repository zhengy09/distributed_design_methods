function [K,X,Z,info] = dsd(A,B,Gp,Gc)
% Distributed stabilization design using chordal decomposition and ADMM
% find K \in S(G)
%  s.t. A + BK is stable
%
% Input: graph Gp  -- the dynamical interconnection graph, a matrix
%        graph Gc  -- the communication graph
%        matrix A  -- cell format  dot(x) = Ax + Bu
%        matrix B  -- cell format

%% parameters
opts.mu      = 1;
opts.maxIter = 20;
opts.verbose = true;          % display output information
opts.subbose = false;          % display output information for each subproblem
opts.eps     = 1e-4;
opts.global  = true;

time = zeros(opts.maxIter,1);

%% finding cliques
n     = size(Gp,1);            % Dimension number of nodes
Gs    = Gp | Gp' | Gc | Gc';   % super-graph, connections in the whole system
Gs    = spones(Gs);            % Sparsity pattern of the LMI (A + BK)'P + P(A + BK)
Ge    = chordalExt(Gs);        % Chordal extension
Mc    = maximalCliques(Ge);    % Maximal cliques, each column is a clique  
Nc    = size(Mc,2);            % Number of cliques

%% process the model data such that Ge and A are consistent
subDimState = zeros(n,1);
subDimInput = zeros(n,1);
for i = 1:n
    subDimState(i) = size(A{i,i},1);
    subDimInput(i) = size(B{i},2);
end

for i = 1:n
    for j = 1:n
        if Ge(i,j) == 1 && isempty(A{i,j})
            A{i,j} = zeros(subDimState(i),subDimState(j));
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
    Clique{k}.nodsize = subDimState(Clique{k}.node);            % the size of local dynamics in each node
    Clique{k}.nodinpt = subDimInput(Clique{k}.node);            % the input dimension of each node
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
    Node{i}.X  = zeros(subDimState(i));                       % Lyapunov variables   ZX^{-1} = Kii
    Node{i}.Z  = zeros(subDimInput(i),subDimState(i));        % Feeback gain variables
    if NodeRep(i) > 1                                         % set local variables for this node
        for k = 1:NodeRep(i)
            % for Mi
            Node{i}.CliqueVariables{k}  = zeros(subDimState(i));   % this variable belong to clique in node.clique(k)
            Node{i}.LocalVariables{k}   = zeros(subDimState(i));   % Local variable in this node
            Node{i}.LocalMultipliers{k} = zeros(subDimState(i));   % Corresponding multipliers
            
            % for Pi
            Node{i}.Xi{k}   = zeros(subDimState(i));                 % This Xi belongs to node.clique(k) 
            Node{i}.XiMultipliers{k} = zeros(subDimState(i));        % Corresponding multipliers
            
            % for time record
            Node{i}.time   = zeros(opts.maxIter,2);
        end
    end
end

for i = 1:n                  
    for j = 1:n   % feedback variables
        if Gc(i,j) == 1 && i ~= j   %% communication? have a feedback gain here
             Edge{i,j}.Z = zeros(subDimInput(i),subDimState(j));
        end
    end
    
    % local variables for consensus
    for j = i+1:n
         if EdgeRep(i,j) > 1
             for k = 1:EdgeRep(i,j)
                 % for Mi
                 Edge{i,j}.CliqueVariables{k}  = zeros(subDimState(i),subDimState(j));
                 Edge{i,j}.LocalVariables{k}   = zeros(subDimState(i),subDimState(j));
                 Edge{i,j}.LocalMultipliers{k} = zeros(subDimState(i),subDimState(j));
                 
                 % for Xi and Xj
                 Edge{i,j}.Xi = zeros(subDimState(i));   % this is equal to node{i}.X
                 Edge{i,j}.Xj = zeros(subDimState(j));     
                 
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
        [Node,Edge,Clique{k}] = cliqueProblem(k,Clique{k},Node,Edge,A,B,Gc,opts,iter);
    end
    
    %% Step 2: Overlapping nodes and edges solve their own subproblem to update
    % local virables, which can be computed in parallel
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        Node{nodeIndex} = nodeProblem(Node{nodeIndex},A{nodeIndex,nodeIndex},B{nodeIndex},opts,iter);
    end
    
    for k = 1:length(EdgeRepi)
        %for j = 1:length(EdgeRepj)
            Eind = [EdgeRepi(k),EdgeRepj(k)];
            A12  = A{Eind(1),Eind(2)};  B1 = A{Eind(1)};
            A21  = A{Eind(2),Eind(1)};  B2 = A{Eind(2)};
            X1   = Node{Eind(1)}.X;
            X2   = Node{Eind(2)}.X;
           % Edge{Eind(1),Eind(2)} = edgeProblem(Edge{Eind(1),Eind(2)},A12,A21,B1,B2,X1,X2,opts,iter);
        %end
    end
    
    %% Step 3: Overlapping nodes and edges update the dual virables
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        Node{nodeIndex} = nodeMultipliers(Node{nodeIndex},opts,iter);
    end
    
    for k = 1:length(EdgeRepi)
            Eind = [EdgeRepi(k),EdgeRepj(k)];
            %Edge{Eind(1),Eind(2)} = edgeMultipliers(Edge{Eind(1),Eind(2)},opts,iter);
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
X       = cell(n,1);
Z       = cell(n,n);
K       = cell(n,n);
tmpEig  = zeros(n,1);
for i = 1:n
    X{i} = Node{i}.X;
    tmpEig(i) = min(eig(X{i}));
    Z{i,i} = Node{i}.Z;
    K{i,i} = Z{i,i}*X{i}^(-1);
    for j = 1:n
        if Gc(j,i) == 1 && i~=j
            Z{j,i} = Edge{j,i}.Z;
            K{j,i} = Z{j,i}*X{i}^(-1);
        end
    end
end

if opts.global
    accDimenState  = [cumsum([1;subDimState])]; 
    accDimenInput  = [cumsum([1;subDimInput])]; 
    globalA   = zeros(sum(subDimState));                 % gloabl state space model
    globalB   = zeros(sum(subDimState),sum(subDimInput));
    globalX   = zeros(sum(subDimState)); 
    globalZ   = zeros(sum(subDimInput),sum(subDimState)); 
    globalK   = zeros(sum(subDimInput),sum(subDimState)); 
    for i = 1:n
        globalX(accDimenState(i):accDimenState(i+1)-1,accDimenState(i):accDimenState(i+1)-1) = X{i};
        globalB(accDimenState(i):accDimenState(i+1)-1,accDimenInput(i):accDimenInput(i+1)-1) = B{i};
        for j = 1:n
            if Gp(i,j) ~= 0
                globalA(accDimenState(i):accDimenState(i+1)-1,accDimenState(j):accDimenState(j+1)-1) = A{i,j};
            end
            if Gc(i,j) ~=0 || i == j
                globalZ(accDimenInput(i):accDimenInput(i+1)-1,accDimenState(j):accDimenState(j+1)-1) = Z{i,j};
                globalK(accDimenInput(i):accDimenInput(i+1)-1,accDimenState(j):accDimenState(j+1)-1) = K{i,j};
            end
        end
    end
    info.globalA = globalA;
    info.globalB = globalB;
    info.globalX = globalX;
    info.globalZ = globalZ;
    info.globalK = globalK;
    
    % test the result
    Closedloop = globalA + globalB*globalK;
    flag = eig(Closedloop*globalX + globalX*Closedloop');
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
    fprintf(' Residual (Pi - P):        %8.2e \n',Info.presi);
    fprintf(' Max. eigenvalue A+BK      %8.2e \n',info.minEig);
    fprintf(' Max. eigenvalue A         %8.2e \n',max(eig(globalA)));
    fprintf(' Closed-loop stable        %8d   \n',info.flag);
    fprintf('--------------------------------------\n')
end

end

