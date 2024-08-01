clear; clc; close all; 
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
allf = [0.5 2.5 4.5];               % vertical force 
allmub = linspace(-10,15,1000);     % effective chemical potential 
[F,MU] = meshgrid(allf,allmub);
% Reminder: size(F) = [length(allmub) length(allf)]

%--------------------------------------------------------------------------
% Output variables [SS = steady state]
%--------------------------------------------------------------------------
numberSS = zeros(length(allmub),length(allf));
steadystates = zeros(length(allmub),length(allf),3); % at most 3 SS

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

    % For range of affinities
    for jj = 1:length(allmub)
        mub = allmub(jj);

        % Find steady states
        SS = zeros(3,1); % 3 SS, 1 feature (bond density, phib)

        % For each interval, check if mub is in the correct range
        for kk = 1:(length(phibintervals)-1)
            % This mub interval
            thismubinterval = [mubintervals(kk) mubintervals(kk+1)];

            if mub > min(thismubinterval) && mub < max(thismubinterval)
                % If mub just right, there is a solution in this interval
                numberSS(jj,ii) = numberSS(jj,ii)+1;

                % This phib interval
                thisphibinterval = [phibintervals(kk) phibintervals(kk+1)];

                % Find steady state
                phib = findsteadystates(P,f,mub,thisphibinterval);

                % Save into SS
                SS(numberSS(jj,ii)) = phib;
            end
        end

        % Record steady state for this mub and f
        steadystates(jj,ii,:) = SS;
    end
    waitbar(ii/length(allf),mywaitbar,[num2str(ii/length(allf)*100) '% done']);
end
toc
close(mywaitbar)

%%
xvar = allmub;
lnwdth = 1.5;

% Color scheme
pink = [232 51 210]/255;
green = [121 251 77]/255;
cyan = [107 231 233]/255;
blue = [86 105 211]/255;
red = [216 13 60]/255;
ochre = [219 198 116]/255;
darkgreen = [51 153 102]/255;
grey = [1 1 1]*116/255;

c = cat(1,red,pink,grey);

figure('Units','centimeters','Position',[15 15 10.93 6])
hold on

% Select value of f
for jj = 1:length(allf)

    % First solution (dilute - background)
    yvar = steadystates(:,jj,1,1);
    plot(xvar,yvar,'-','Color',c(jj,:),'MarkerFaceColor','k','LineWidth',lnwdth);
    
    % Second solution (unstable)
    yvar = steadystates(:,jj,2,1);
    yvar(yvar ==0) = NaN;
    plot(xvar,yvar,'--','Color',c(jj,:),'MarkerFaceColor','r','LineWidth',lnwdth)

    % Third solution (dense - clusters)
    yvar = steadystates(:,jj,3,1);
    yvar(yvar ==0) = NaN;
    plot(xvar,yvar,'-','Color',c(jj,:),'MarkerFaceColor','b','LineWidth',lnwdth);
    
    % Identify saddle-node bifurcations
    yvar = cat(1,steadystates(1,jj,2,1),steadystates(:,jj,2,1),steadystates(end,jj,2,1));
    nonzero = yvar>0;
    % Current value must be non-zero, but either the value immediately
    % before or immediately after it myst be zero
    bifurcation = nonzero(2:end-1) & (~nonzero(1:end-2) | ~nonzero(3:end));

    % Indicate saddle-node bifurcations with circle
    yvar = steadystates(:,jj,3,1);
    plot(xvar(bifurcation),yvar(bifurcation),'o','Color',c(jj,:),'MarkerFaceColor',[1 1 1],'LineWidth',lnwdth,'HandleVisibility','off')
end

% Determine legend
mylegend = cell(1,3*length(allf));
for ii = 1:length(allf)
    mylegend{3*ii-2} = ['$f= ' num2str(allf(ii)) '$~pN'];
    mylegend{3*ii-1} = '';
    mylegend{3*ii} = '';
end

legend(mylegend,'Location','southwest','Interpreter','latex')
xlabel('$\mu_b~(k_BT)$','Interpreter','latex')  
ylabel('$\phi_b$','Interpreter','latex')  
box on
axis([allmub(1) allmub(end) -0.05 1.05])
set(gca,'FontSize',11)

box on