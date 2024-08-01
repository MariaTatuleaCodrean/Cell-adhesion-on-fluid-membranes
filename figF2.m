clear; close all; clc;

% Color scheme
pink = [232 51 210]/255;
green = [121 251 77]/255;
cyan = [107 231 233]/255;
blue = [86 105 211]/255;
red = [216 13 60]/255;
ochre = [219 198 116]/255;
darkgreen = [51 153 102]/255;
grey = [1 1 1]*116/255;

%% Definitions of vectors and functions

% parameters
chi0 = 15;
Delta = 15;
f0 = 50; % pN
f1 = 2; % pN
f = 2.1;
myaxislim = [0 1 -12 2];

% vector of phib concentrations
phib=linspace(0.0001,0.9999);

% chemical potential [x = chi0]
mub=@(x) log(phib./(1-phib))-x*phib./(1+exp(Delta*(1-f/f1)))+f/f0;

%% Drawing figure
figure('Units','centimeters','Position',[10,10,17,7])

% plot of fixed epsilon and changing f
subplot(1,2,1)
hold on
plot(phib,mub(2),'Color',red,'linewidth',1.5);
plot(phib,mub(11.75),'Color',blue,'linewidth',1.5);
plot(phib,mub(chi0),'k','linewidth',1.5);
plot(phib,-5*phib./phib,'k--','linewidth',1.5);

xlabel('bond density $\phi_b$','interpreter','latex','FontName','times');
ylabel('chemical potential $\mu_b$','interpreter','latex','FontName','times');    

box on
axis square
axis(myaxislim)

% plot of fixed f and changing epsilon
subplot(1,2,2)
hold on
plot(phib,-10*phib./phib,'--','Color',red,'linewidth',1.5);
plot(phib,-6.95*phib./phib,'--','Color',blue,'linewidth',1.5);
plot(phib,-5*phib./phib,'k--','linewidth',1.5);
plot(phib,-3.15*phib./phib,'--','Color',blue,'linewidth',1.5);
plot(phib,mub(chi0),'k','linewidth',1.5);
plot(phib,0*phib./phib,'--','Color',red,'linewidth',1.5);

xlabel('bond density $\phi_b$','interpreter','latex','FontName','times');
ylabel('chemical potential $\mu_b$','interpreter','latex','FontName','times');    

box on
axis square
axis(myaxislim)