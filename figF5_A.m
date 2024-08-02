% clear; close all; clc;
%--------------------------------------------------------------------------
% Universal constants
%--------------------------------------------------------------------------
kT = 4.11; % pN*nm

%--------------------------------------------------------------------------
% Model parameters
%--------------------------------------------------------------------------
P.f1 = 2;  % force to activate adhesome proteins [pN]
P.f0 = 50; % force to break bond [pN]
P.chi0 = 15; % interaction strength, units of kT
P.Pi = 15; % activation energy, units of kT
P.r = 25;  % ratio of spring constants, r=kb/k0
P.fb = 20; % elastic force on stretched bond  [pN]
P.Eel = 200; % elastic energy, units of kT
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Input variables
%--------------------------------------------------------------------------
allf = linspace(0,5,1000);       % vertical force 
allmub = linspace(-10,15,1000);  % effective chemical potential 
[F,MU] = meshgrid(allf,allmub);
% Reminder: size(F) = [length(allmub) length(allf)]

%--------------------------------------------------------------------------
% Output variables
%--------------------------------------------------------------------------
mub_crit = NaN(size(allf));  % critical clustering threshold

%--------------------------------------------------------------------------
% Main routine
%--------------------------------------------------------------------------
tic
numberofchiplots = 0;
mywaitbar = waitbar(0,'Entering for loop');
% For range of forces
for ii = 1:length(allf)
    f = allf(ii);

     % For each force, find turning points of g(x) (= chemical potential)
    [phibintervals, mubintervals] = findturningpoints(P,f);

    % Find local minimum chemical potential (mub) after removing -Inf and +Inf
    % values from vector epsilintervals
    mubintervals = mubintervals(isfinite(mubintervals));

    if ~isempty(mubintervals)
        mub_crit(ii) = min(mubintervals);
    end
    waitbar(ii/length(allf),mywaitbar,[num2str(ii/length(allf)*100) '% done']);
end
toc
close(mywaitbar)

%% Plot
% Color scheme
pink = [232 51 210]/255;
green = [121 251 77]/255;
cyan = [107 231 233]/255;
blue = [86 105 211]/255;
red = [216 13 60]/255;
ochre = [219 198 116]/255;
darkgreen = [51 153 102]/255;
grey = [1 1 1]*116/255;

% Size (enlarge by factor of 3)
myFontSize = 10*3;
myLineWidth = 2*3;
myThinWidth = 1*3;

% Generate figure to correct size
figure('Units','centimeters','Position',[5 5 7*3 6*3])
hold on

% Define corners
myXlim = [min(allf) max(allf)];
myYlim = [min(allmub) max(allmub)];

% Define clustering region
% curved edge
X = allf; 
Y = mub_crit;
% add top right corner
X = [X(:); myXlim(2)]; 
Y = [Y(:); myYlim(2)];
% add top left corner
X = [X(:); myXlim(1)]; 
Y = [Y(:); myYlim(2)];

% Define background
Xb = [myXlim(1); myXlim(2); myXlim(2); myXlim(1)];
Yb = [myYlim(1); myYlim(1); myYlim(2); myYlim(2)];

% Plot
patch(Xb,Yb,grey,'EdgeColor','none');               % nay, clusters
patch(X,Y,darkgreen,'EdgeColor','none');            % yea, clusters
plot(allf,mub_crit,'w-','LineWidth',myLineWidth);   % clustering threshold

% Add horizontal lines for special values of mub
yline(-8,'Color','k','LineWidth',myThinWidth)
yline(2,'Color',pink,'LineWidth',myThinWidth)
yline(12,'Color',red,'LineWidth',myThinWidth)

% Add vertical lines for special values of f
xline(4.5,'Color','k','LineWidth',myThinWidth)
xline(2.5,'Color',pink,'LineWidth',myThinWidth)
xline(0.5,'Color',red,'LineWidth',myThinWidth)

xlabel('$f~(pN)$','Interpreter','latex')
ylabel('$\mu_b~(k_BT)$','Interpreter','latex')   
xlim(myXlim)
ylim(myYlim)
box on
set(gca,'FontSize',myFontSize)