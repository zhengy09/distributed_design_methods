
%function Test
   clc;
%% Generating data
 %   N      = 10;
  %  thresh = 5;
  %  G    = chordalGen(N,thresh);

    N = 5;
     G = [1 1 1 1 0 ;1 1 1 1 1;1 1 1 1 1; 1 1 1 1 1; 0 1 1 1 1 ];
    
    A = cell(N,N);
    B = cell(N,1);
    M = cell(N,1);
    Q = cell(N,1);
    R = cell(N,1);
    
    n = ones(1,N)*2;
    m = ones(1,N);
    d = ones(1,N);
    %n = randi([3,5],1,N);
    %m = randi([1,2],1,N);
    %d = randi([1,2],1,N);
    for i = 1:N
       A{i,i} = [1 2;1 2];%rand(n(i));
       B{i} = [0;1];%rand(n(i),m(i));
       M{i} = [0;1];%rand(n(i),d(i));
       Q{i} = 0*eye(n(i));
       R{i} = 0.1*eye(m(i));
       for j = 1:N
           if G(i,j) == 1 && i ~= j
               A{i,j} = 0.1*rand(n(i),n(j));
           end
       end
    end
    
    %% Test
    [K,Cost, info1] = ddd(G,A,B,M,Q,R);
    
    [K1,Cost1, info1] = cdd(G,A,B,M,Q,R);
   
    [Cost, Cost1]
    %norm(K)
%     try
%     [P, info] = dst(A,G);
%     catch
%         k = 1;
%     end
%     
    
   
   
   
   
   
