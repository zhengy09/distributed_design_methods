function  Edge = edgeMultipliers(Edge,opts,iter)
%   Detailed explanation goes here

    tstart = tic;
    for i = 1:length(Edge.clique)
        Edge.LocalMultipliers{i}  = Edge.LocalMultipliers{i} + ... 
                        opts.mu*(Edge.CliqueVariables{i} - Edge.LocalVariables{i});
    end
    Edge.time(iter,2) = toc(tstart); 
end

