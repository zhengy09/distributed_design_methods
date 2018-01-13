%function Test
   
close all;
    


%% random generate a sparse positive semidefinite matrix
    n      = 5;
    % thresh = 3;
    % R    = chordalGen(n,thresh);
    
%     R = [1 1 1 1 0;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1; 0 1 1 1 1];
%     
%     X = rand(n).*R;
%     X = (X+X')/2;
%     X = X + (-min(eig(X))+0.1)*eye(n);
%     
    %X = [4 1 1 0;1 2 1 1;1 1 3 1; 0 1 1 4];
    
    X = [4 1 0;1 3 1;0 1 3]
    
    % figure 1
    %

%% Decomposition
    [Xk, Mc, info] = dcd(X);
    [Xk1, Mc1, info1] = adc(X);
    
%% figure
%     figure; 
%     subplot(2,length(Xk),1);spy(X); hold on;
%     for i = 1:length(Xk)
%         hold on;
%         subplot(2,length(Xk),length(Xk)+i);spy(info.Xk{i});
%     end
%end
