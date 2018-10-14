
%% Figure -- convergence results

clc;clear;close all
load ConData.mat

figure;
histogram(ConIterR,10);
xlabel('Iteration'); ylabel('Frequency')

figure;
histogram(ConIterR,30,'Normalization','cdf');
set(gca,'YTick',[0:0.1:1],'TickLabelInterpreter','latex','fontsize',9)
xlabel('Iteractions to convergence','Interpreter','latex'); grid on;
set(gcf,'Position',[250 150 350 250]);

folderAddress = 'F:/Goolge Drive/1 Paper writing/2 Decentralized stabilization of linear systems/Journal paper/Figures/';
print(gcf,[folderAddress,'AglerIterations'],'-painters','-depsc2','-r600')
