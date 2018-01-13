function [Stop, Info] = ConverCheck(NodeOld,Node,Edge,NodeOverInd,EdgeRepi,EdgeRepj,opts)
% Check convergence

Stop = false;
presi = 0;

dresi = 0;

for i = 1:length(NodeOverInd)
    tmpNode = Node{NodeOverInd(i)};
    for k = 1:length(tmpNode.clique)
        presi = presi + norm(tmpNode.X - tmpNode.Xi{k},'fro').^2;
        presi = presi + norm(tmpNode.CliqueVariables{k} - tmpNode.LocalVariables{k},'fro').^2;
    end
end

for i = 1:length(EdgeRepi)
    
    Eind = [EdgeRepi(i),EdgeRepj(i)];  
    presi = presi + norm(Edge{Eind(1),Eind(2)}.Xi - Node{Eind(1)}.X,'fro').^2;
    presi = presi + norm(Edge{Eind(1),Eind(2)}.Xj - Node{Eind(2)}.X,'fro').^2;
    
end

cost = 0;
for i = 1:length(Node)
    cost = cost + trace(Node{i}.Q*Node{i}.X) + trace(Node{i}.R*Node{i}.Y);
    dresi = dresi + norm(Node{i}.X - NodeOld{i}.X,'fro').^2;
    dresi = dresi + norm(Node{i}.Y - NodeOld{i}.Y,'fro').^2;
    dresi = dresi + norm(Node{i}.Z - NodeOld{i}.Z,'fro').^2;
end
Info.cost = cost;
Info.presi = presi;
Info.dresi = dresi;
if presi < opts.eps && dresi < opts.eps
    Stop = true;
end



