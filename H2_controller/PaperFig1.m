
%% Figures for the first example

clc;clear; close all;
load PaperEx2.mat

IndPro = find(H2(:,3)>50);

IndPro1 = find(H2(IndPro,3) <  inf);
IndPro = IndPro(IndPro1);

mH2 = H2;
mH2(IndPro,3) = inf;

%% Unstable feedback for each method

NumSuc = zeros(5,1);    % the number of succesful instances
IndSuc = cell(5,1);   % the corresponding indices
for i = 1:5
    IndSuc{i} = find(mH2(:,i) < inf);
    NumSuc(i) = length(IndSuc{i});
end

fprintf('%d\n',NumSuc);


%% Performance improvement for each method

H2sel = cell(3,1);
H2mea = zeros(3,2);
PerIm = zeros(3,1);
for i = 1:3
    H2sel{i} = mH2(IndSuc{i+2},[1,i+2]);
    H2mea(i,:) = mean(mH2(IndSuc{i+2},[1,i+2]));
    PerIm(i)   = (H2mea(i,2) - H2mea(i,1))./H2mea(i,2);
end
fprintf('%3.4f\n',PerIm);

[mean(mH2(IndSuc{1})),mean(mH2(IndSuc{3})),mean(mH2(IndSuc{4})),mean(mH2(IndSuc{5}))]

%% Performance improvement for the common instances
C = intersect(IndSuc{3},IndSuc{4});
C = intersect(C,IndSuc{5});
length(C)

H2com = H2(C,:);
H2comM = mean(H2com);

for i = 1:3
    PercomIm(i) = (H2comM(1) - H2comM(i+2))./H2comM(i+2);
end

H2comM
PercomIm

%% Figures
figure;
histogram(Iter,30,'Normalization','cdf','FaceColor','[0.8, 0.8, 0.8]');
set(gca,'YTick',[0:0.2:1],'TickLabelInterpreter','latex','fontsize',9)
xlabel('Iterations to convergence','Interpreter','latex'); grid on;
set(gcf,'Position',[250 150 350 250]);

folderAddress = 'F:/Goolge Drive/1 Paper writing/2 Decentralized stabilization of linear systems/Journal paper/Figures/';
print(gcf,[folderAddress,'ChainIterations'],'-painters','-depsc2','-r600')






