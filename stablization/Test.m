
%function Test
   
%% Generating data
%     N      = 30;
%     thresh = 5;
%     G    = chordalGen(N,thresh);
    
clc;clear
    N = 3;
    n = 2; m = 1;
    Gp = [1 0 0;1 1 0;0 1 1];
    Gc = Gp';
    A = cell(N,N);
    B = cell(N,1);
    %n = ones(1,N);
    for i = 1:N
        B{i} = rand(n,m);%[1;1];
        for j = 1:N
            if Gp(i,j) == 1
                A{i,j} = rand(n,n);
            end
        end
    end
    
    dsd(A,B,Gp,Gc);
   

%end
   
   
   
   
   
