
clc;clear
%% a network of 4 nodes 

N = 4;
A = cell(N,N);
B = cell(N,1);

Gp = [1 0 0 0;1 1 0 0;0 1 1 1; 1 1 0 1];
A{1,1} = [1];A{2,1} = [1];A{4,1} = [1];
A{2,2} = [2];A{3,2} = [2];A{4,2} = [2];
A{3,3} = [3];
A{3,4} = [4];A{4,4} = [4];
B{1} = [1];B{2} = [1];B{3} = [1];B{4} = [1];
M = B;
Q = B;
R = B;

    %% solutions via different methods
    
    % ADMM 
    % centralized design
    [K1, Cost1, info1] = cdd(Gp,A,B,M,Q,R);
    
    if info1.h2 ~= inf
        opts.mu = 50;
        [K,Cost, info] = ddd(Gp,A,B,M,Q,R,opts);
    else
        info = info1;
        info.iter = 1;
        K = K1;
    end

    
    gA = info1.gA;
    gB = info1.gB;
    gM = info1.gM;
    gQ = eye(size(gA,1));
    gR = eye(size(gB,2));
    
    
    % sequential design
    [K2,~,~,Tdata,Tsolver,Tgraph,Ttotal] = scSeq(A,B,Gp,Gc);
    ClosedSys = ss(gA + gB*K2,gM,[gQ^(1/2);gR^(1/2)*K2],[]);
    h20  = norm(ClosedSys,2);
    
    % trancated LQR
    K3 = lqr(gA,gB,gQ,gR);
    K3 = K3.*spones(K1);
    ClosedSys = ss(gA - gB*K3,gM,[gQ^(1/2);gR^(1/2)*K3],[]);
    h21  = norm(ClosedSys,2);

    % loalized LQR
    K41 = lqr(A{1,1},B{1},Q{1},R{1});
    K4 = [];
    for i = 1:N
        K4 = blkdiag(K4,K41);
    end
    ClosedSys = ss(gA - gB*K4,gM,[gQ^(1/2);gR^(1/2)*K4],[]);
    h22  = norm(ClosedSys,2);
    
    
    %% statistics
    Dynamic{Index}   = A;
    Iter(Index)      = info.iter;
    H2(Index,:)      = [info.h2,info1.h2,h20,h21,h22];
    Feedback{Index,1} = K;
    Feedback{Index,2} = K1;
    Feedback{Index,3} = K2;
    Feedback{Index,4} = K3;
    Feedback{Index,5} = K4;
    SolInfo{Index,1} = info;
    SolInfo{Index,2} = info1;

