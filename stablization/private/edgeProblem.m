function Edge = edgeProblem(Edge,A12,A21,B1, B2, X1,X2,Gc,opts,iter)
% Update the local variables in each edge

NonLocalVariables = length(Edge.clique);
LocalVariables    = cell(NonLocalVariables,1);
[ni,mi]  = size(A12);
for i = 1:NonLocalVariables
    LocalVariables{i} = sdpvar(ni,mi);
end

Z    = cell(2,1);
Z{1} = sdpvar(size(B1,2),size());
Z{2} = sdpvar();

%% define the cost function
Cost = 0;
for i = 1:NonLocalVariables
    Cost = Cost + norm(Edge.CliqueVariables{i} - LocalVariables{i} - 1/opts.mu*Edge.LocalMultipliers{i},'fro').^2;
end

%% define the constraints
tmp = LocalVariables{1};
for i = 2:NonLocalVariables
    tmp = tmp + LocalVariables{i};
end

if Gc(1) == 0 && Gc(2) == 0
    Constraints = [tmp == -(X2*A21'+X1*A12)];
elseif Gc(1) == 1 && Gc(2) == 0
     Constraints = [tmp == -(X2*A21'+X1*A12)];
elseif Gc(1) == 0 && Gc(2) == 1
     Constraints = [tmp == -(X2*A21'+X1*A12)];
else
     Constraints = [tmp == -(X2*A21'+X1*A12)];
end

%% Get solutions
options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
sol     = optimize(Constraints,Cost,options);
Edge.time(iter) = sol.solvertime; 
%% set values
for i = 1:NonLocalVariables
    Edge.LocalVariables{i} = value(LocalVariables{i});
end

end

