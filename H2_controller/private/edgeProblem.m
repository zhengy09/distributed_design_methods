function Edge = edgeProblem(Edge,A12,A21,X1,X2,opts,iter)
% Update the local variables in each edge

NonLocalVariables   = length(Edge.clique);
LocalVariables      = cell(NonLocalVariables,1);
Xi   = sdpvar(size(Edge.Xi,1));
Xj   = sdpvar(size(Edge.Xj,1));
[ni,mi]  = size(A12);
for i = 1:NonLocalVariables
    LocalVariables{i} = sdpvar(ni,mi);
end

%% define the cost function
Cost = 0;
for i = 1:NonLocalVariables
    Cost = Cost + norm(Edge.CliqueVariables{i} - LocalVariables{i} + 1/opts.mu*Edge.LocalMultipliers{i},'fro').^2;
    Cost = Cost + norm(X1{i} - Xi + 1/opts.mu*Edge.XiMultipliers{i},'fro').^2;
    Cost = Cost + norm(X2{i} - Xj + 1/opts.mu*Edge.XjMultipliers{i},'fro').^2;
end

%% define the constraints
tmp = LocalVariables{1};
for i = 2:NonLocalVariables
    tmp = tmp + LocalVariables{i};
end
Constraints = [tmp == -(Xi*A21'+A12*Xj)];

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
sol     = optimize(Constraints,Cost,options);
Edge.time(iter) = sol.solvertime; 
%% set values
Edge.Xi = value(Xi);
Edge.Xj = value(Xj);
for i = 1:NonLocalVariables
    Edge.LocalVariables{i} = value(LocalVariables{i});
end

end
