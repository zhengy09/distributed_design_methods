
clc;clear;close all

%% paramters

% distributed generator units
% Rt1 = 0.2; Ct1 = 2.2 * 10^(-3); Lt1 = 1.8 * 10^(-3);
% Rt2 = 0.3; Ct2 = 1.9 * 10^(-3); Lt2 = 2.0 * 10^(-3);
% Rt3 = 0.1; Ct3 = 1.7 * 10^(-3); Lt3 = 2.2 * 10^(-3);
% Rt4 = 0.5; Ct4 = 2.5 * 10^(-3); Lt4 = 3.0 * 10^(-3);
% Rt5 = 0.4; Ct5 = 2.0 * 10^(-3); Lt5 = 1.2 * 10^(-3);

Rt1 = 0.2; Ct1 = 2.2 * 10^(0); Lt1 = 1.8 * 10^(0);
Rt2 = 0.3; Ct2 = 1.9 * 10^(0); Lt2 = 2.0 * 10^(0);
Rt3 = 0.1; Ct3 = 1.7 * 10^(0); Lt3 = 2.2 * 10^(0);
Rt4 = 0.5; Ct4 = 2.5 * 10^(0); Lt4 = 3.0 * 10^(0);
Rt5 = 0.4; Ct5 = 2.0 * 10^(0); Lt5 = 1.2 * 10^(0);

% Transmission lines
R12 = 0.05;  L12 = 2.1 * 10^(-6); 
R13 = 0.07;  L13 = 1.8 * 10^(-6);
R34 = 0.06;  L34 = 1.0 * 10^(-6);
R24 = 0.04;  L24 = 2.3 * 10^(-6);
R45 = 0.08;  L45 = 1.8 * 10^(-6);

%% Dynmaics
A11 = [-1/(R12*Ct1)-1/(R13*Ct1) 1/Ct1; -1/Lt1 -Rt1/Lt1];
A12 = [1/(R12*Ct1) 0;0 0]; A13 = [-1/(R13*Ct1) 0; 0 0];
B1 = [0; 1/Lt1]; 
M1 = [-1/Ct1; 0]; 
H1 = [1 0];

A22 = [-1/(R12*Ct2)-1/(R24*Ct2) 1/Ct2; -1/Lt2 -Rt2/Lt2];
A21 = [1/(R12*Ct2) 0;0 0]; A24 = [-1/(R24*Ct2) 0; 0 0];
B2 = [0; 1/Lt2]; 
M2 = [-1/Ct2; 0]; 
H2 = [1 0];

A33 = [-1/(R13*Ct3)-1/(R34*Ct3) 1/Ct3; -1/Lt3 -Rt3/Lt3];
A31 = [1/(R13*Ct3) 0;0 0]; A34 = [-1/(R34*Ct3) 0; 0 0];
B3 = [0; 1/Lt3]; 
M3 = [-1/Ct3; 0]; 
H3 = [1 0];

A44 = [-1/(R24*Ct4)-1/(R34*Ct4)-1/(R45*Ct4) 1/Ct4; -1/Lt4 -Rt4/Lt4];
A42 = [1/(R24*Ct4) 0;0 0]; A43 = [-1/(R34*Ct4) 0; 0 0]; A45 = [-1/(R45*Ct4) 0; 0 0];
B4 = [0; 1/Lt4]; 
M4 = [-1/Ct2; 0]; 
H4 = [1 0];

A55 = [-1/(R45*Ct5) 1/Ct5; -1/Lt5 -Rt3/Lt5];
A54 = [1/(R45*Ct5) 0;0 0]; 
B5 = [0; 1/Lt5]; 
M5 = [-1/Ct5; 0]; 
H5 = [1 0];

%% 
hA = cell(5,5);
hB = cell(5,1);
hM = cell(5,1);

%% Dynamics with integrator
hA11 = [A11 zeros(2,1); -H1 0];  
hA12 = blkdiag(A12,0); hA13 = blkdiag(A13,0);
hB1 = [B1;0]; 
hM1 = blkdiag(M1,1); 
hH1 = [H1,0];

hA{1,1} = hA11; hA{1,2} = hA12; hA{1,3} = hA13;
hB{1} = hB1;
hM{1} = hM1;

hA22 = [A22 zeros(2,1); -H2 0];
hA21 = blkdiag(A21,0); hA24 = blkdiag(A24,0);
hB2 = [B2;0]; 
hM2 = blkdiag(M2,1); 
hH2 = [H2,0];

hA{2,2} = hA22; hA{2,1} = hA21; hA{2,4} = hA24;
hB{2} = hB2;
hM{2} = hM2;

hA33 = [A33 zeros(2,1); -H3 0];
hA31 = blkdiag(A31,0); hA34 = blkdiag(A34,0);
hB3 = [B3;0]; 
hM3 = blkdiag(M3,1); 
hH3 = [H3,0];

hA{3,3} = hA33; hA{3,1} = hA31; hA{3,4} = hA34;
hB{3} = hB3;
hM{3} = hM3;

hA44 = [A44 zeros(2,1); -H4 0];
hA42 = blkdiag(A42,0); hA43 = blkdiag(A43,0); hA45 = blkdiag(A45,0);
hB4 = [B4;0]; 
hM4 = blkdiag(M4,1); 
hH4 = [H4,0];

hA{4,4} = hA44; hA{4,2} = hA42; hA{4,3} = hA43;hA{4,5} = hA45;
hB{4} = hB4;
hM{4} = hM4;

hA55 = [A55 zeros(2,1); -H5 0];
hA54 = blkdiag(A54,0); 
hB5 = [B5;0]; 
hM5 = blkdiag(M5,1); 
hH5 = [H5,0];

hA{5,5} = hA44; hA{5,4} = hA54; 
hB{5} = hB5;
hM{5} = hM5;

%% global dynamics

Ag = [hA11 hA12 hA13 zeros(3) zeros(3);
      hA21 hA22 zeros(3) hA24 zeros(3);
      hA31 zeros(3) hA33 hA34 zeros(3);
      zeros(3) hA42 hA43 hA44 hA45;
      zeros(3) zeros(3) zeros(3) hA54 hA55];
Bg = blkdiag(hB1,hB2,hB3,hB4,hB5);
Mg = blkdiag(hM1,hM2,hM3,hM4,hM5);


Q1 = eye(3); R1 = 0.1;
Q = cell(5,1);R = cell(5,1);
 Q{1} = eye(3);Q{2} = eye(3);Q{3} = eye(3);Q{4} = eye(3);Q{5} = eye(3);
 R{1} = R1; R{2} = R1; R{3} = R1; R{4} = R1; R{5} = R1;
 
 G = [1 1 1 0 0;
      1 1 0 1 0;
      1 0 1 1 0
      0 1 1 1 1;
      0 0 0 1 1];
  
 
 [K1,Cost1, info1] = cdd(G,hA,hB,hM,Q,R);
 
 opts.eps = 1.0e-3;
 opts.MaxIter = 500;
 opts.mu  = 10;
 [K,Cost, info] = ddd(G,hA,hB,hM,Q,R,opts);
 
 [Cost, Cost1]

%% decentralized H2 controller by Yalmip

% Q  = blkdiag(Q1,Q1,Q1,Q1,Q1);R = blkdiag(R1,R1,R1,R1,R1);
% 
% % variables
% X1 = sdpvar(3);
% X  = blkdiag(X1,X1,X1,X1,X1);
% Z1 = sdpvar(1,3);
% Z  = blkdiag(Z1,Z1,Z1,Z1,Z1);
% Y1 = sdpvar(1);
% Y  = blkdiag(Y1,Y1,Y1,Y1,Y1);
% 
% % cost
% Obj = trace(Q*X) + trace(R*Y);
% % constraints
% Const = (Ag*X - Bg*Z) + (Ag*X - Bg*Z)' + Mg*Mg' <= 0;
% Const = [Const, [Y Z; Z' X]>=0, X-0.01*eye(15) >=0];
% sol = optimize(Const,Obj);
% Ko  = value(Z)*value(X)^(-1);
% 
% [Cost, value(Obj)]
% 
% Aco = Ag - Bg*Ko;
% 
% %% decentralized controller
% K1 = lqr(hA11,hB1,eye(3),0.1);
% K2 = lqr(hA22,hB2,eye(3),0.1);
% K3 = lqr(hA33,hB3,eye(3),0.1);
% K4 = lqr(hA44,hB4,eye(3),0.1);
% K5 = lqr(hA55,hB5,eye(3),0.1);
% 
% K = blkdiag(K1,K2,K3,K4,K5);
% 
% Ac = Ag - Bg*K;
% 
% C1 = [1 0 0];
% C = blkdiag(C1,C1,C1,C1,C1);
% 
% %% response
% 
% % lsim(sys,u,t,x0)
% t = 0:0.1:400;
% V1 = 50; V2 = 50; V3 = 50; V4 = 50; V5 = 50;
% d = [0*ones(1,length(t));V1*ones(1,length(t));
%     0*ones(1,length(t));V2*ones(1,length(t));
%     0*ones(1,length(t));V3*ones(1,length(t));
%     [10*ones(1,2000),5+0*sin(2*pi/10*t(2000+1:end))];V4*ones(1,length(t));
%     0*ones(1,length(t));V5*ones(1,length(t));];
% 
% %10+10*sin(2*pi/10*t(2000+1:end))
% 
% % sys1 -- optimal 
% sys1 = ss(Aco,Mg,C,[]);
% sys2 = ss(Ac,Mg,C,[]);
% 
% figure;
% lsim(sys1,d',t)
% 
% figure;
% lsim(sys2,d',t)



%max(real(eig(Ac)))







