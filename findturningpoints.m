function [phibintervals, mubintervals] = findturningpoints(P,f)
%FINDTURNINGPOINTS
% This routine takes as input the parameters of the theoretical model 
% and a given cytoskeletal force f (measured in pN) and returns:
% 1) the intervals of phib between 0 and 1 where the effective
% chemical potential is monotonically increasing or decreasing, and
% 2) the values of the chemical potential at the ends of those intervals.
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
% phibintervals [vector]    phib values between which the chemical
%                           potential is monotonically
%                           increasing/decreasing
% mubintervals [vector]   values of chemical potential at the end of
%                           those intervals
% by construction, phibintervals(1)  = 0,    phibintervals(end)  = 1
%                  mubintervals(1) = -Inf, mubintervals(end) = +Inf
%--------------------------------------------------------------------------

% Effective interaction parameter
chi = (P.chi0)/(1+exp((P.Pi)*(1-f/(P.f1))));

% Quintic coefficients
r = P.r;                                % frequently used parameter
prefactor = 2*r*(P.Eel)*(1-f/(P.fb));   % frequently used parameter
a0 = 1;
a1 = 3*r   - chi*1                 - prefactor;
a2 = 3*r^2 - chi*3*r   + chi       + prefactor*(1-r*f/(P.fb));
a3 = r^3   - chi*3*r^2 + chi*3*r   + prefactor*r*f/(P.fb);
a4 =       - chi*r^3   + chi*3*r^2;
a5 =                     chi*r^3;

% Chemical potential as a function of phib is given by
mub = @(phib) log(phib./(1-phib))+f/(P.f0)-chi*phib+(P.Eel)*(1+r*phib*f/(P.fb)).^2./(1+r*phib).^2;

% For reference, the derivative of mub w.r.t. phib is given by
% dmubdphib = @(phib) (a5*phib^5 + a4*phib^4 + a3*phib^3 + a2*phib^2 + ...
%     a1*phib + a0)./(phib*(1-phib)*(1+r*phib)^3);

% Find roots of a5*x^5 + a4*x^4 + ... + a0 = 0
% i.e. find where $\frac{\partial\mu_b}{\partial\phi_b} = 0$
rts = roots([a5 a4 a3 a2 a1 a0]);

% Select only real roots
rts = rts(~imag(rts));

% Sort in ascending order
rts = sort(rts);

% Only roots between 0 and 1 are of interest
rts = rts(rts<1);
rts = rts(rts>0);

switch length(rts)
    case 4 
        % Four roots between 0 and 1, we have two max TPs and two min TPs
        phibintervals = [0 rts(1) rts(2) rts(3) rts(4) 1];
    case 2
        % Two roots between 0 and 1, we have both max and min TPs
        phibintervals = [0 rts(1) rts(2) 1];

    case 0
        % No roots, we have a single interval (entire domain)
        phibintervals = [0 1];

    otherwise
        error('findturningpoints: number of turning points is not standard')
end
% Debug
if length(rts) > 2
    disp(['number of turning points = ' num2str(length(rts))])
end

% Determine affinity interval
mubintervals = mub(rts);
mubintervals = [-Inf; mubintervals; +Inf];

end