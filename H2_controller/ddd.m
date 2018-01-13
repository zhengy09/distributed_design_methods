function [K, Cost, info] = ddd(G,A,B,M,Q,R,userOpts)
% Distributed Design of Decentralized Controoler using chordal decomposition and ADMM
% dot(x) = Ax + Bu + Md
% u = Kx, K is block diagonal
% Q,R cell structure, performance index

%% parameters
opts.mu      = 10;
opts.maxIter = 500;
opts.verbose = true;          % display output information
opts.subbose = false;          % display output information for each subproblem
opts.eps     = 1e-4;
opts.global  = true;

if(nargin >= 7)
    fnames = fieldnames(userOpts);
    for n=1:length(fnames)
        if isfield(opts,fnames{n})
            opts.(fnames{n}) = userOpts.(fnames{n});
        end
    end
end


time = zeros(opts.maxIter,1);

%% finding cliques
n     = size(G,1);            % Dimension number of nodes
G     = spones(G+eye(n));            % Sparsity pattern
Ge    = chordalExt(G);        % Chordal extension
Mc    = maximalCliques(Ge);   % Maximal cliques, each column is a clique  
Nc    = size(Mc,2);           % Number of cliques

%% process the model data such that Ge and A are consistent
subDim  = zeros(n,1);         % State dimension
subDimI = zeros(n,1);         % Input dimension 
subDimd = zeros(n,1);         % disturbance dimension
for i = 1:n
    subDim(i) = size(A{i,i},1);
    subDimI(i) = size(B{i},2);
    subDimd(i) = size(M{i},2);
end

for i = 1:n
    for j = 1:n
        if Ge(i,j) == 1 && isempty(A{i,j})
            A{i,j} = zeros(subDim(i),subDim(j));
        end
    end
end

%% Overlapping information
tic;
Node        = cell(n,1);             % Nodes
NodeRep     = sum(Mc,2);             % repetition of nodes in different cliques
NodeOverInd = find(NodeRep > 1);     % nodes that belong to more than one clique

% Node information
for i = 1:n
    Node{i}.Id = i;
    Node{i}.clique = find(Mc(i,:) == 1);   % the cliques that contain this node
end

% Clique information
Clique = cell(Nc,1);                 
for k = 1:Nc
    % statistic 
    Clique{k}.node    = find(Mc(:,k) == 1);                     % the nodes in this clique
    Clique{k}.nodsize = subDim(Clique{k}.node);                 % the size of local dynamics in each node
    Clique{k}.sdpsize = sum(Clique{k}.nodsize);                 % the dimension of sdp constraint
    Clique{k}.nodeR   = intersect(Clique{k}.node,NodeOverInd);  % Repeated nodes in this clique
    % Clique{k}.edgeR   = [];                                     % Repeated edges in this clique
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
    % dynamics data
    Node{i}.B = B{i};
    Node{i}.M = M{i};
    Node{i}.Q = Q{i};
    Node{i}.R = R{i};
    
    % Lyapunov variables
    Node{i}.X  = zeros(subDim(i));               % Lyapunov variables
    Node{i}.Y  = zeros(subDimI(i));              % Lyapunov variables
    Node{i}.Z  = zeros(subDimI(i),subDim(i));    % Lyapunov variables
    if NodeRep(i) > 1                    % set local variables for this node
        for k = 1:NodeRep(i)
            % for Ji
            Node{i}.CliqueVariables{k}  = zeros(subDim(i));   % this variable belong to clique in node.clique(k)
            Node{i}.LocalVariables{k}   = zeros(subDim(i));   % Local variable in this node
            Node{i}.LocalMultipliers{k} = zeros(subDim(i));   % Corresponding multipliers
            
            % for Xi
            Node{i}.Xi{k}   = zeros(subDim(i));               % This Pi belongs to node.clique(k) 
            Node{i}.XiMultipliers{k} = zeros(subDim(i));      % Corresponding multipliers
            
            % for time record
            Node{i}.time   = zeros(opts.maxIter,2);
        end
    end
end

for i = 1:n                  % local variables in edges
    for j = i+1:n
         if EdgeRep(i,j) > 1
             for k = 1:EdgeRep(i,j)
                 % for Ji
                 Edge{i,j}.CliqueVariables{k}  = zeros(subDim(i),subDim(j));
                 Edge{i,j}.LocalVariables{k}   = zeros(subDim(i),subDim(j));
                 Edge{i,j}.LocalMultipliers{k} = zeros(subDim(i),subDim(j));
                 
                 % for Xi and Xj
                  Edge{i,j}.Xi = zeros(subDim(i));   % this is for consensus
                  Edge{i,j}.Xj = zeros(subDim(j));     
                  
                  Edge{i,j}.XiMultipliers{k} = zeros(subDim(i));      % Corresponding multipliers
                  Edge{i,j}.XjMultipliers{k} = zeros(subDim(j));      % Corresponding multipliers
                 
                 % for time record
                 %Edge{i}.time = zeros(opts.maxIter,2);
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
    fprintf(' Iter.    presi    dresi    Time (s)  Cost\n')
end

tic
for iter = 1:opts.maxIter
    
    NodeOld = Node;  % for convergence checking
    %% Step 1: Each clique solve a subproblem in parallel 
    CliqueInfea = 0;
    for k = 1:Nc
        [Node,Edge,Clique{k},sol] = cliqueProblem(k,Clique{k},Node,Edge,A,opts,iter);
        if sol.problem == 1
            CliqueInfea = CliqueInfea + 1;
        end
    end
    
    %% Step 2: Overlapping nodes and edges solve their own subproblem to update
    % local virables, which can be computed in parallel
    NodeInfea = 0;
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        [Node{nodeIndex},sol] = nodeProblem(Node{nodeIndex},A{nodeIndex,nodeIndex},opts,iter);
        if sol.problem == 1
            NodeInfea = NodeInfea + 1;
        end
    end
    
    EdgeInfea = 0;
    for k = 1:length(EdgeRepi)
        Eind = [EdgeRepi(k),EdgeRepj(k)];
        A12  = A{Eind(1),Eind(2)};
        A21  = A{Eind(2),Eind(1)};
        X1   = Node{Eind(1)}.Xi;
        X2   = Node{Eind(2)}.Xi;
        [Edge{Eind(1),Eind(2)},sol] = edgeProblem(Edge{Eind(1),Eind(2)},A12,A21,X1,X2,opts,iter);
        if sol.problem == 1
            EdgeInfea = EdgeInfea + 1;
        end
    end
    
    %% Step 3: Overlapping nodes and edges update the dual virables
    for i = 1:length(NodeOverInd)
        nodeIndex = NodeOverInd(i);
        Node{nodeIndex} = nodeMultipliers(Node{nodeIndex},opts,iter);
    end
    
    for k = 1:length(EdgeRepi)
        Eind = [EdgeRepi(k),EdgeRepj(k)];
        X1   = Node{Eind(1)}.Xi;
        X2   = Node{Eind(2)}.Xi;
        Edge{Eind(1),Eind(2)} = edgeMultipliers(Edge{Eind(1),Eind(2)},X1,X2,opts,iter);
    end
    %% check convergence & output
    [Stop, Info] = ConverCheck(NodeOld,Node,Edge,NodeOverInd,EdgeRepi,EdgeRepj,opts);
    time(iter) = toc;
    if opts.verbose == true
        fprintf('%5d   %8.4f  %8.4f   %6.2f  %6.3f %4d %4d %4d\n', iter,Info.presi,Info.dresi,time(iter),Info.cost,CliqueInfea,NodeInfea,EdgeInfea);
    end
    
    if Stop == true
        %break;
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

info.Node = Node;
info.Edge = Edge;
%%
X       = [];
Y       = [];
Z       = [];
Cost    = 0;
tmpEig  = zeros(n,1);
for i = 1:n
    X = blkdiag(X,Node{i}.X);
    Y = blkdiag(Y,Node{i}.Y);
    Z = blkdiag(Z,Node{i}.Z);
    Cost = Cost + trace(Q{i}*Node{i}.X)+trace(R{i}*Node{i}.Y);
end

info.X = value(X);
info.Y = value(Y);
info.Z = value(Z);

K = Z*X^(-1);



gB   = [];
gM   = [];
gQ = [];
gR = [];

%if opts.global
    accDimen  = [cumsum([1;subDim])]; 
    gA   = zeros(sum(subDim));  % gloabl state space model
    for i = 1:n
        gB = blkdiag(gB,B{i});
        gM = blkdiag(gM,M{i});
        gQ = blkdiag(gQ,Q{i});
        gR = blkdiag(gR,R{i});
        for j = 1:n
            if G(i,j) ~= 0
                gA(accDimen(i):accDimen(i+1)-1,accDimen(j):accDimen(j+1)-1) = A{i,j};
            end
        end
    end
    info.gA = gA;
    info.gB = gB;
    info.gM = gM;
    
    ClosedSys = ss(gA - gB*K,gM,[gQ^(1/2);gR^(1/2)*K],[]);
    info.h2  = norm(ClosedSys,2);
    
    % test the result
%     flag = min(eig(globalA'*globalP + globalP*globalA));
%     info.minEig = max(flag);
%     if info.minEig  < 0
%         info.flag = true;
%     else
%         info.flag = false;
%     end
%end

%% Summaries
if opts.verbose == true
    fprintf('--------------------------------------\n')
    fprintf(' Residual (Pi - P):      :   %8.2e \n',Info.presi)
    %fprintf(' Minimal eigenvalue of Pk:   %8.2e \n',min(tmpEig))
    fprintf('--------------------------------------\n')
end

end

