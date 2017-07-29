function [Stop, Info] = ConverCheck(Xk,X,Node,Edge,Mc,EdgeRepi,EdgeRepj,opts)
% Check convergence

Stop = false;
n = size(X,1);
presi = 0;
dresi = 0;

% % consider the variables in nodes 
% for i = 1:n
%     Cnode   = Node{i};
%     if length(Cnode.clique) > 1
%         tmpresi = 0;
%         for li = 1:length(Cnode.clique)
%             Clique     = Mc(:,Cnode.clique(li));              % the nodes that this clique contains
%             lNodeIndex = sum(Clique(1:i));                   % local position of this node in Node.clique(i)
%             presi      = presi + (Xk{Cnode.clique(li)}(lNodeIndex,lNodeIndex) - Cnode.LocalVariables(li))^2;
%             tmpresi    = tmpresi +  Cnode.LocalVariables(li);
%         end
%         presi = presi + (tmpresi - X(i,i))^2;
%     end
% end
% 
% % consider the variables in edges
% for i = 1:length(EdgeRepi)
%         for j = 1:length(EdgeRepj)
%             EdgeIndex = [EdgeRepi(i),EdgeRepj(j)];
%             CEdge     = Edge{EdgeIndex(1),EdgeIndex(2)};
%             tmpresi = 0;
%             for locali = 1:length(CEdge.clique)
%                 Clique      = Mc(:,CEdge.clique(locali));                  % the nodes that this clique contains
%                 lNodeIndexi = sum(Clique(1:EdgeIndex(1)));          % local position of this node in Node.clique(i)
%                 lNodeIndexj = sum(Clique(1:EdgeIndex(2)));          % local position of this node in Node.clique(i)
%                 presi  =  presi + (Xk{CEdge.clique(locali)}(lNodeIndexi,lNodeIndexj) - CEdge.LocalVariables(locali))^2;
%                 tmpresi = tmpresi + CEdge.LocalVariables(locali);
%             end
%             presi = presi + (tmpresi - X(EdgeIndex(1),EdgeIndex(2)))^2;
%         end
% end

for k = 1:length(Xk)
   gXk{k} = globalMatrix(Xk{k},Mc(:,k),n);
end

tmpX = zeros(n);
for i = 1:length(Xk)
        tmpX = tmpX + gXk{i};
end
presi = norm(tmpX-X);


Info.presi = presi;
if presi < opts.eps
    Stop = true;
end



