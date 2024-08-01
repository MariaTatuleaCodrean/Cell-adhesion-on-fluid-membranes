function [phib] = findsteadystates(P,f,mub,phibinterval)
%FINDSTEADYSTATES
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

% Effective interaction parameter
chi = (P.chi0)/(1+exp((P.Pi)*(1-f/(P.f1))));

% Define function representing (roughly) the chemical potential
fchemicalpotential = @(x) log(x./(1-x))+f/(P.f0)-chi*x+(P.Eel)*(1+(P.r)*x*f/(P.fb)).^2./(1+(P.r)*x).^2 - mub;

% Use tangent trick to remove infinite values
fun = @(x) atan(fchemicalpotential(x));

% Find solution in this interval
phib = fzero(fun,phibinterval);

end