%% Figures for the third examples

clc;clear

load Ex1Microgrid

%% response

gA = info3.gA; gB = info3.gB; gM = info3.gM;

C1 = [1 0 0];
C = blkdiag(C1,C1,C1,C1,C1);

Ac = gA - gB*K;


% lsim(sys,u,t,x0)
t = 0:0.1:600;
V1 = 50; V2 = 50; V3 = 50; V4 = 50; V5 = 50;
% d = [0*ones(1,length(t));V1*ones(1,length(t));
%     0*ones(1,length(t));V2*ones(1,length(t));
%     0*ones(1,length(t));V3*ones(1,length(t));
%     0*ones(1,length(t));V4*ones(1,length(t));
%     0*ones(1,length(t));V5*ones(1,length(t));];

d = [0*ones(1,length(t));V1*ones(1,length(t));
    0*ones(1,length(t));V2*ones(1,length(t));
    0*ones(1,length(t));[V3*ones(1,4000),48*ones(1,length(t)-4000)];
    0*ones(1,length(t));V4*ones(1,length(t));
    0*ones(1,length(t));V5*ones(1,length(t));];

%10+10*sin(2*pi/10*t(2000+1:end))

% sys1 -- optimal 
sys1 = ss(Ac,gM,C,[]);
figure;
lsim(sys1,d',t)
