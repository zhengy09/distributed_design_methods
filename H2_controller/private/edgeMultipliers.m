function  Edge = edgeMultipliers(Edge,X1,X2,opts,iter)
%   Detailed explanation goes here

    tstart = tic;
    for i = 1:length(Edge.clique)
        Edge.LocalMultipliers{i}  = Edge.LocalMultipliers{i} + ... 
                        (Edge.CliqueVariables{i} - Edge.LocalVariables{i});
                    
        Edge.XiMultipliers{i}    = Edge.XiMultipliers{i} + (X1{i} - Edge.Xi);
        Edge.XjMultipliers{i}    = Edge.XjMultipliers{i} + (X2{i} - Edge.Xj);
    end
    Edge.time(iter,2) = toc(tstart); 
end

