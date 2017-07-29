function Test
   
close all;
    


%% random generate a sparse positive semidefinite matrix
    n      = 5;
    thresh = 4;
    R    = chordalGen(n,thresh);
    
    X = rand(n).*R;
    X = (X+X')/2;
    X = X + (-min(eig(X))+0.1)*eye(n);
    
    % figure 1
    %

%% Decomposition
    [Xk, Mc, info] = dcd(X);
    
%% figure
    figure; 
    subplot(2,length(Xk),1);spy(X); hold on;
    for i = 1:length(Xk)
        hold on;
        subplot(2,length(Xk),length(Xk)+i);spy(info.Xk{i});
    end
end
