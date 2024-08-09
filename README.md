# Cell-adhesion-on-fluid-membranes
This public repository contains the MATLAB codes to make theoretical predictions and to generate the figures in the Supplemental Material for our article:

**Cell adhesion and spreading on fluid membranes through microtubules-dependent mechanotransduction**

Oleg Mikhajlov, Ram M. Adar, Maria Tătulea-Codrean, Anne-Sophie Macé, John Manzi, Fanny Tabarin, Aude Battistella, Fahima di Federico, Jean-François Joanny, Guy Tran van Nhieu, Patricia Bassereau

doi: https://doi.org/10.1101/2022.09.12.507658

The folder contains these MATLAB functions:

* **figF2.m** - plot Figure F2 from Supplemental Material (corresponds to Supp. Fig. S13 from main manuscript)
* **FigF3.m** - plot Figure F3 from Supplemental Material (corresponds to Fig. 4F from main manuscript)
* **FigF5_X.m** - plot Figure F5 panel X (X being from A to H) (corresponds to Supp. Fig. S14 from main manuscript)


* **findturningpoints.m**
This routine takes as input the parameters of the theoretical model and a given cytoskeletal force f (measured in pN) and returns: 

1) the intervals of phib between 0 and 1 where the effective chemical potential is monotonically increasing or decreasing, and
2) the values of the chemical potential at the ends of those intervals.


% INPUT

% P [structure] model parameters

%   P.f1        force to activate adhesome proteins [pN]

%   P.f0        force to break bond [pN]

%   P.chi0      interaction strength, units of kT

%   P.Pi        activation energy, units of kT

%   P.r         ratio of spring constants, r=kb/k0

%   P.fb        elastic force on stretched bond  [pN]

%   P.Eel       elastic energy, units of kT

%

%   f           vertical force applied by microtubules to each bond

%

% Units

% all forces:   pN

% all lengths:  nm

% all times:    s

% all energies in units of kT

%

% OUTPUT

% phibintervals [vector]    phib values between which the chemical

%                           potential is monotonically

%                           increasing/decreasing

% mubintervals [vector]   values of chemical potential at the end of

%                           those intervals

% by construction, phibintervals(1)  = 0,    phibintervals(end)  = 1

%                  mubintervals(1) = -Inf, mubintervals(end) = +Inf

%--------------------------------------------------------------------------


* %FINDSTEADYSTATES
% This routine takes as input the parameters of the theoretical model,
% a given cytoskeletal force f (measured in pN) and chemical potential 
% (measured in kBT), and it returns the solution for the bond density
% within the interval indicated by the input parameter phibinterval.
%
% INPUT
% P [structure] model parameters
%   P.f1        force to activate adhesome proteins [pN]
%   P.f0        force to break bond [pN]
%   P.chi0      interaction strength, units of kT
%   P.Pi        activation energy, units of kT
%   P.r         ratio of spring constants, r=kb/k0
%   P.fb        elastic force on stretched bond  [pN]
%   P.Eel       elastic energy, units of kT
%
%   f           vertical force applied by microtubules to each bond
%
% Units
% all forces:   pN
% all lengths:  nm
% all times:    s
% all energies in units of kT
%
% OUTPUT
%   phib        bond density at steady state
%--------------------------------------------------------------------------
