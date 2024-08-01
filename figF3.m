clear; close all; clc;
%--------------------------------------------------------------------------
% Universal constants
%--------------------------------------------------------------------------
kT = 4.11; % pN*nm

%--------------------------------------------------------------------------
% Model parameters
%--------------------------------------------------------------------------
P.f1 = 2;  % force to activate adhesome proteins [pN]
P.fb = 20; % elastic force on stretched bond  [pN]
P.f0 = 50; % force to break bond [pN]
P.Pi = 15; % activation energy, units of kT
chi0 = 15; % interaction strength, units of kT
Eel = 200; % elastic energy, units of kT
P.r = 25;  % ratio of spring constants, r=kb/k0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot phase diagram over a range of forces
%--------------------------------------------------------------------------
fvec = linspace(0,2*(P.f1)); % input: range of forces between 0 and 5*f0
mubvec = [-10 15];  % and a range of chemical potentials

%--------------------------------------------------------------------------
% Main routine
%--------------------------------------------------------------------------
P.chi0  = 0; % interaction strength, units of kT
P.Eel   = Eel;  % elastic energy,       units of kT
mub_min = NaN(size(fvec));    % output: minimum mub for clustering
mub_max = NaN(size(fvec));    % output: maximum mub for clustering

% For range of forces
for ii = 1:length(fvec)
    f = fvec(ii);

    % For each force, find turning points of g(x) (= chemical potential)
    [phibintervals, mubintervals] = findturningpoints(P,f);

    % Find local minimum chemical potential (mub) after removing -Inf and +Inf
    % values from vector epsilintervals
    mubintervals = mubintervals(isfinite(mubintervals));

    if ~isempty(mubintervals)
        mub_min(ii) = min(mubintervals);
        mub_max(ii) = max(mubintervals);
    end
end

%--------------------------------------------------------------------------
% Add membrane-mediated interactions
%--------------------------------------------------------------------------
P.chi0  = chi0; % interaction strength, units of kT
P.Eel   = Eel;  % elastic energy,       units of kT
mub2_min = NaN(size(fvec));    % output: minimum mub for clustering
mub2_max = NaN(size(fvec));    % output: maximum mub for clustering

% For range of forces
for ii = 1:length(fvec)
    f = fvec(ii);

    % For each force, find turning points of g(x) (= chemical potential)
    [phibintervals, mubintervals] = findturningpoints(P,f);

    % Find local minimum chemical potential (mub) after removing -Inf and +Inf
    % values from vector epsilintervals
    mubintervals = mubintervals(isfinite(mubintervals));

    if ~isempty(mubintervals)
        mub2_min(ii) = min(mubintervals);
        mub2_max(ii) = max(mubintervals);
    end
end

%--------------------------------------------------------------------------
%% Draw phase diagram - patch
%--------------------------------------------------------------------------
Xlim = [min(fvec) max(fvec)];
Ylim = [mubvec(1) mubvec(2)];

% Color scheme
pink = [232 51 210]/255;
green = [121 251 77]/255;
cyan = [107 231 233]/255;
blue = [86 105 211]/255;
red = [216 13 60]/255;
ochre = [219 198 116]/255;
darkgreen = [51 153 102]/255;

pink_pastel = [250 214 246]/255;
grey = [1 1 1]*217/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero chi0, define edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = fvec; % curved edge
Y = mub_min;

X = [X(:); Xlim(2)]; % add top right corner
Y = [Y(:); Ylim(2)];

X = [X(:); Xlim(1)]; % add top left corner
Y = [Y(:); Ylim(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite chi0, define edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X2 = fvec; % curved edge
Y2 = mub2_min;

X2 = [X2(:); Xlim(2)]; % add top right corner
Y2 = [Y2(:); Ylim(2)];

X2 = [X2(:); Xlim(1)]; % add top left corner
Y2 = [Y2(:); Ylim(2)];

% Draw phase diagram using patch function
figure('Units','centimeters','Position',[10,10,17,9])
patch(X2,Y2,pink_pastel,'EdgeColor','none');
hold on
patch(X,Y,grey,'EdgeColor','none');
plot(X2(1:end-2),Y2(1:end-2),'-','Color',pink,'LineWidth',2);
plot(X(1:end-2),Y(1:end-2),'--','Color','k','LineWidth',2);
plot([Xlim(1) Xlim(2) Xlim(2) Xlim(1) Xlim(1)], [Ylim(1) Ylim(1) Ylim(2) Ylim(2) Ylim(1)],'k-')

% Add mub estimates for RGD 0.1%, RGD 2% and invasin
yline(-2,'k:','LineWidth',2)
yline(1,'k:','LineWidth',2)
yline(5,'k:','LineWidth',2)

set(gca,'fontname','times')
set(gca,'fontsize',11)

xlabel('force $f~(pN)$','Interpreter','latex')
ylabel('chemical potential $\mu_b~(k_BT)$','Interpreter','latex')
box on
axis square
axis([Xlim(1) Xlim(2) Ylim(1) Ylim(2)])