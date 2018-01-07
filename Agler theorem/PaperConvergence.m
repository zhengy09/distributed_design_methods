%function Test
   
close all;
clc;clear;

Num = 100;

n = 50;
thresh = 10;

ConIter   = zeros(Num,1);
TimeTotal = zeros(Num,1); 

for i = 1:Num


%% random generate a sparse positive semidefinite matrix
    R = chordalGen(n,thresh);
    X = rand(n).*R;
    X = (X+X')/2;
    X = X + (-min(eig(X))+0.1)*eye(n);
    [Xk, Mc, info] = dcd(X);
    
    ConIter(i) = info.iter;
    TimeTotal(i) = info.time;    
end
 