clc;clear;close all

load Ex2Inverted1;

t = 0:0.1:19.5; 

%alpha = 1;

% x0 = [0.05; 0; 0.2; 0;0.05; 0; 0.2; 0;0.05; 0; 0.2; 0;];
x0 = [5/180*2*pi; 0; 0.2; 0;0.0; 0; 0; 0;0; 0; 0; 0;];
x1 = zeros(length(t),length(x0));
x1n = zeros(length(t),1);
x2 = zeros(size(x1));
x2n = zeros(length(t),1);
x3 = zeros(length(t),1);
for i = 1:length(t)
    x1(i,:) = expm((info2.gA-info2.gB*K1)*t(i))*x0;
    x2(i,:) = expm((info2.gA-info2.gB*K2)*t(i))*x0;
    x1n(i) = norm(x1(i,:),'fro'); 
    x2n(i) = norm(x2(i,:),'fro');
end

figure; h1 = plot(t,x2n,'m','linewidth',1); hold on
h2 = plot(t,x1n,'b','linewidth',1); 

h = legend([h1,h2],'Centralized','Sequentail','Location','Northwest');
set(h,'box','off')
set(gcf,'Position',[250 150 350 350]);box off, grid off; %xlim([10 200]); ylim([1,10^4])
%set(gca,'XTick',[10,20,50,100,150,200]) 
xlabel('Time (s)');ylabel('Norm of error |x(t)|')

figure;
%plot(t,x2(:,1:4:12)*180/2/pi,'linewidth',2); 
plot(t,x2(:,1)*180/2/pi,'k-.','linewidth',1); hold on 
plot(t,x2(:,5)*180/2/pi,'k--','linewidth',1); hold on
plot(t,x2(:,9)*180/2/pi,'k-','linewidth',1); hold on
xlim([0 10]); 
set(gcf,'Position',[250 150 300 350]);box off, grid off;
xlabel('time $(s)$','interpreter','latex','fontsize',12);
ylabel('Vertical degree $\theta_i (^o)$','interpreter','latex','fontsize',12);
h = legend('Pendulum 1','Pendulum 2','pendulum 3','Location','Northeast');
set(h,'box','off')


folderAddress = 'F:/Goolge Drive/1 Paper writing/2 Decentralized stabilization of linear systems/Journal paper/Figures/';
print(gcf,[folderAddress,'theta'],'-painters','-depsc2','-r600')

figure;
%plot(t,x2(:,3:4:12),'linewidth',2);
plot(t,x2(:,3),'k-.','linewidth',1); hold on 
plot(t,x2(:,7),'k--','linewidth',1); hold on
plot(t,x2(:,11),'k-','linewidth',1); hold on
xlim([0 15])
set(gcf,'Position',[250 150 300 350]);box off, grid off;
xlabel('time $(s)$','interpreter','latex','fontsize',12)
ylabel('Horizontal displacement $y_i (m)$ ','interpreter','latex','fontsize',12)
h = legend('Pendulum 1','Pendulum 2','Pendulum 3','Location','Northeast');
set(h,'box','off')

print(gcf,[folderAddress,'displacement'],'-painters','-depsc2','-r600')

set(h,'box','off')