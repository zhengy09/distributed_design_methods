function [K,Cost, info] = cdd(G,A,B,M,Q,R)
% Centralized Design of decentralized optimal controller
% dot(x) = Ax + Bu + Md
% u = Kx, K is block diagonal
% Q,R cell structure, performance index

opts.subbose = true;
epsilon      = 1e-2;

tic
n            = size(G,1);            % Dimension number of nodes
subDimension  = zeros(n,1);          % state dimension
subDimensionI = zeros(n,1);          % input dimension
subDimensiond = zeros(n,1);          % input dimension
for i = 1:n
    subDimension(i) = size(A{i,i},1);
    subDimensionI(i) = size(B{i},2);
    subDimensiond(i) = size(M{i},2);
end

% global form
accDimen  = [cumsum([1;subDimension])]; 
gA   = zeros(sum(subDimension));  % gloabl state space model
X         = [];
Y         = [];
Z         = [];  
gB   = [];
gM   = [];
gQ = [];
gR = [];
%globalP   = zeros(sum(subDimension)); 
for i = 1:n
    X = blkdiag(X,sdpvar(subDimension(i)));
    Y = blkdiag(Y,sdpvar(subDimensionI(i)));
    Z = blkdiag(Z,sdpvar(subDimensionI(i),subDimension(i)));
    gB = blkdiag(gB,B{i});
    gM = blkdiag(gM,M{i});
    gQ = blkdiag(gQ,Q{i});
    gR = blkdiag(gR,R{i});
    for j = 1:n
        if G(i,j) ~= 0
            gA(accDimen(i):accDimen(i+1)-1,accDimen(j):accDimen(j+1)-1) = A{i,j};
        end
    end
end

%Y = sdpvar(sum(subDimensionI));


%% define Cost and constraints
Constraints = [(gA*X - gB*Z) + (gA*X - gB*Z)'+ gM*gM' <=0];
Constraints = [Constraints, X - epsilon*eye(sum(subDimension)) >=0];

Constraints = [Constraints, [Y Z;Z' X] >= 0];

Cost = trace(gQ*X) + trace(gR*Y);

options = sdpsettings('verbose',opts.subbose,'solver','sedumi');
sol     = optimize(Constraints,Cost,options);

timeTotal    = toc;
info.time    = [timeTotal,sol.solvertime]; 
info.gA = gA;
info.gB = gB;
info.gM = gM;

%% set values
X = value(X);
Cost = value(Cost);
K = value(Z)*X^(-1);

end

