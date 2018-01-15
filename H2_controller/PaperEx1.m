
clc;clear
%% a chain of 5 nodes 

% N = 5;
% 
% Gp = zeros(N,N);            % Plant graph: a line 
% Gp(1,1) = 1; Gp(1,2) = 1;
% Gp(N,N-1) = 1; Gp(N,N) = 1;
% for i = 2:N-1
%     Gp(i,i-1) = 1;
%     Gp(i,i) = 1;
%     Gp(i,i+1) = 1;
% end
% Gc = eye(N,N); % fully decentralized controller
% ew3
% 
% NumSample = 100;
% 
% H2      = zeros(NumSample,5);
% Iter    =  zeros(NumSample,1);
% SolInfo = cell(NumSample,2);
% Dynamic = cell(NumSample,1);
% 
% Feedback = cell(NumSample,5);

load PaperEx2

for Index = 97:NumSample
    
    rand('state',sum(100*clock));

    fprintf('Iteration number: %d\n', Index);
    A = cell(N,N);
    B = cell(N,1);
    M = cell(N,1);
    Q = cell(N,1);
    R = cell(N,1);
    
    n = ones(1,N)*2;
    m = ones(1,N);
    d = ones(1,N);
    a = -0.5; b = 0.5;
    for i = 1:N
       A{i,i} = [1 1;1 2];
       B{i} = [0;1];
       M{i} = [0;1];
       Q{i} = 1*eye(n(i));
       R{i} = 1*eye(m(i));
       for j = 1:N
           if Gp(i,j) == 1 && i ~= j
               A{i,j} = a + (b-a)*rand(n(i),n(j));
           end
       end
    end
    
    %% solutions via different methods
    
    % ADMM 
    % centralized design
    [K1, Cost1, info1] = cdd(Gp,A,B,M,Q,R);
    
    if info1.h2 ~= inf
        [K,Cost, info] = ddd(Gp,A,B,M,Q,R);
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
    
    save PaperEx2
end   

