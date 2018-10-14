
%% Counter example
%  a decentralized controller exists for a full Lyapunov function
%  but a decentralized controller fails for a diagonal Lyapunov function 

clc;clear

n = 3;

A = rand(n);

A = A - (0.01+max(real(eig(A))))*eye(n);

Q = 10*eye(n)

X = lyap(A',Q)

eig(X)

A'*X + X*A


A = [1 2; -1 0.5]
B = [0 0;0 1]

K = [0 0; 0 2]

A - B*K

eig(A - B*K)

X = [];
Z = [];
for i = 1:2
    X = blkdiag(X,sdpvar(1));
    Z = blkdiag(Z,sdpvar(1));
end

Objective = 0;
epsilon = 1e-2;
Constraint = [(A*X - B*Z) + (A*X - B*Z)' + epsilon*eye(2)<=0, X- epsilon *eye(2) >= 0];

options = sdpsettings('solver','sedumi');
sol     = optimize(Constraint,Objective,options);

K1 = value(Z)*value(X)^(-1);

K2 = K;
K2(2,2) = K1(2,2)
eig(A - B*K2)



X1 = sdpvar(1);X2 = sdpvar(1);
Z2 = sdpvar(1);


