function Node = nodeMultipliers(Node,opts,iter)
% Update the multipliers in each node

    tstart = tic;
    for i = 1:length(Node.clique)
        Node.LocalMultipliers{i}  = Node.LocalMultipliers{i} ...
                                + opts.mu*(Node.CliqueVariables{i} - Node.LocalVariables{i});
        Node.XiMultipliers{i}  = Node.XiMultipliers{i} + opts.mu*(Node.Xi{i} - Node.X);
    end
    
    Node.time(iter,2) = toc(tstart); 
end

