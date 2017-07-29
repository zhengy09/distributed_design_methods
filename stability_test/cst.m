function [P, info] = cst(G,A)
% Centralized Stability Test
%

opts.subbose = false;
epsilon      = 1e-2;

tic
n     = size(G,1);            % Dimension number of nodes
subDimension = zeros(n,1);
for i = 1:n
    subDimension(i) = size(A{i,i},1);
end

% global form
accDimen  = [cumsum([1;subDimension])]; 
globalA   = zeros(sum(subDimension));  % gloabl state space model
P         = [];
%globalP   = zeros(sum(subDimension)); 
for i = 1:n
    P = blkdiag(P,sdpvar(subDimension(i)));
    %globalP(accDimen(i):accDimen(i+1)-1,accDimen(i):accDimen(i+1)-1) = P{i};
    for j = 1:n
        if G(i,j) ~= 0
            globalA(accDimen(i):accDimen(i+1)-1,accDimen(j):accDimen(j+1)-1) = A{i,j};
        end
    end
end


%% define Cost and constraints
Constraints = [-globalA'*P-P*globalA - epsilon*eye(sum(subDimension)) >=0];
Constraints = [Constraints, P - epsilon*eye(sum(subDimension)) >=0];
Cost = 0;

options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
sol     = optimize(Constraints,Cost,options);

timeTotal    = toc;
info.time    = [timeTotal,sol.solvertime]; 
info.globalA = globalA;

%% set values
P = value(P);

end

