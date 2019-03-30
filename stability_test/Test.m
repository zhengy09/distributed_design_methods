
%function Test
   
%% Generating data
    N      = 20;
    thresh = 5;
    G    = chordalGen(N,thresh);

%     N = 3;
%     G = [1 1 0;1 1 1;0 1 1];
    
    A = cell(N,N);
    n = randi(3,1,N);
    %n = ones(1,N);
    for i = 1:N
       A{i,i} = rand(n(i));
       %A{i,i} = (A{i,i} + A{i,i})./2;
       A{i,i} = A{i,i} + (-max(real(eig(A{i,i})))-2).*eye(n(i));
       for j = 1:N
           if G(i,j) == 1 && i ~= j
               A{i,j} = 0.1*rand(n(i),n(j));
           end
       end
    end
    
    %% Test
    [P1, info1] = cst(G,A);
    try
    [P, info] = dst(A,G);
    catch
        k = 1;
    end
    
    
    %% Centralized solution
    
    
%     R = 10;
%     XY = zeros(N,2);
%     Degree = 0:2*pi/N:2*pi;
%     for i = 1:N
%         XY(i,:) = R.*[sin(Degree(i)),cos(Degree(i))];
%     end
%     
%     gplot(G,XY,'-*')
%     axis square

%end
   
   
   
   
   
