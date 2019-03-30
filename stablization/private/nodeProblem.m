function Node = nodeProblem(Node,A,B,opts,iter)
% Update the local variables in each node

epsilon = 1e-2;

NonLocalVariables = length(Node.clique);
LocalVariables    = cell(NonLocalVariables,1);
for i = 1:NonLocalVariables
    LocalVariables{i} = sdpvar(size(A,1));
end
X = sdpvar(size(A,1));
Z = sdpvar(size(B,2),size(B,1));

%% define the cost function
Cost = 0;
for i = 1:NonLocalVariables
    Cost = Cost + norm(Node.CliqueVariables{i} - LocalVariables{i} + 1/opts.mu*Node.LocalMultipliers{i},'fro').^2;
    Cost = Cost + norm(Node.Xi{i} - X + 1/opts.mu*Node.XiMultipliers{i},'fro').^2;
end

%% define the constraints
tmp = LocalVariables{1};
for i = 2:NonLocalVariables
    tmp = tmp + LocalVariables{i};
end
Constraints = [tmp == -(A*X+X*A'+B*Z+ (B*Z)' + epsilon*eye(size(A)))];
Constraints = [Constraints, X-epsilon*eye(size(A)) >=0];

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
sol     = optimize(Constraints,Cost,options);

Node.time(iter,1) = sol.solvertime; 

%% set values
Node.X = value(X);
Node.Z = value(Z);
for i = 1:NonLocalVariables
    Node.LocalVariables{i} = value(LocalVariables{i});
end

end

