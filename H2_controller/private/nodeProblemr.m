function [Node,sol] = nodeProblemr(Node,A,opts,iter)
% Update the local variables in each node

epsilon = 1e-2;%1e-2;

NumLocalVariables = length(Node.clique);               % number of local variables
LocalVariables    = cell(NumLocalVariables,1);
for i = 1:NumLocalVariables
    LocalVariables{i} = sdpvar(size(A,1));
end
X = sdpvar(size(A,1));
Y = sdpvar(size(Node.Y,1),size(Node.Y,2),'full');
Z = sdpvar(size(Node.Z,1),size(Node.Z,2),'full');

%% define the cost function
Cost = trace(Node.Q*X+Node.R*Y);
for i = 1:NumLocalVariables
    Cost = Cost + opts.mu/2*norm(Node.CliqueVariables{i} - LocalVariables{i} + Node.LocalMultipliers{i},'fro').^2;
    Cost = Cost + opts.mu/2*norm(Node.Xi{i} - X + Node.XiMultipliers{i},'fro').^2;
end


%% define the constraints
tmp = LocalVariables{1};
for i = 2:NumLocalVariables
    tmp = tmp + LocalVariables{i};
end
%Constraints = [tmp == -(A'*X+X*A+epsilon*eye(size(A)))];
Constraints = [tmp == -(A*X - Node.B*Z) - (A*X - Node.B*Z)' - Node.M*Node.M'];
Constraints = [Constraints, X-epsilon*eye(size(A)) >=0, [Y Z;Z' X] >=0];

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi','cachesolvers',1);
sol     = optimize(Constraints,Cost,options);

Node.time(iter,1) = sol.solvertime; 

%% set values
    if sol.problem ~= 1
        Node.X = value(X);
        Node.Y = value(Y);
        Node.Z = value(Z);
        for i = 1:NumLocalVariables
            Node.LocalVariables{i} = value(LocalVariables{i});
        end
    end
end

