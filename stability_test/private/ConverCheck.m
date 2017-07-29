function [Stop, Info] = ConverCheck(Node,Edge,NodeOverInd,EdgeRepi,EdgeRepj,opts)
% Check convergence

Stop = false;
presi = 0;

for i = 1:length(NodeOverInd)
    tmpNode = Node{NodeOverInd(i)};
    for k = 1:length(tmpNode.clique)
        presi = presi + norm(tmpNode.P - tmpNode.Pi{k},'fro').^2;
        presi = presi + norm(tmpNode.CliqueVariables{k} - tmpNode.LocalVariables{k},'fro').^2;
    end
end

Info.presi = presi;
if presi < opts.eps
    Stop = true;
end



