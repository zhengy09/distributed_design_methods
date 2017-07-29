
function Test
   
%% Generating data
    N      = 10;
    thresh = 5;
    G    = chordalGen(N,thresh);

    A = cell(N,N);
    n = randi(5,1,N);
    %n = ones(1,N);
    for i = 1:N
       A{i,i} = rand(n(i));
       A{i,i} = (A{i,i} + A{i,i})./2;
       A{i,i} = A{i,i} + (-max(eig(A{i,i}))-5).*eye(n(i));
       for j = 1:N
           if G(i,j) == 1 && i ~= j
               A{i,j} = 0.01*rand(n(i),n(j));
           end
       end
    end
    
    %% Test
    [P, info] = dst(A,G);
    
    
    R = 10;
    XY = zeros(N,2);
    Degree = 0:2*pi/N:2*pi;
    for i = 1:N
        XY(i,:) = R.*[sin(Degree(i)),cos(Degree(i))];
    end
    
    gplot(G,XY,'-*')
    axis square

end
   
   
   
   
   
