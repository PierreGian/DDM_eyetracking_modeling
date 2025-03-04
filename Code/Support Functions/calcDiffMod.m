function [t1, t2, p] = calcDiffMod(coh, theta, subj_eq)
% calculates a hitting time function for a particle in Brownian motion to
% barrier at height +A or -B.  The model transforms the intensity parmeters in
% vector x, to an expectation, k*x, and assumes that the variance is 1 at each
% time step.  
%
%  send theta = [k u0 A B t0_t1 t0_t2] as a vector

% 5/17/03 mns wrote it
% 9/12/03 mns made all coh terms into abs(k*x). Now handles negative coh
% 9/13/03 mns allows vectors for parameter args
% 5/13/04 tdh modified for stim fits
% 7-05 modified again by TDH 

if nargin < 3
    subj_eq = 0;
end

% Extract parameters
k = theta(1);
u0 = theta(2);
A = theta(3);
B = theta(4);
t0_t1 = theta(5);
t0_t2 = theta(6);

% make nan vectors for return values
t1 = nan(size(coh));
t2 = nan(size(coh));
p = nan(size(coh));

% this keeps u0 in units of coherence
x = coh + u0;
u = k.*x;

% Deal with limiting case where drift is near zero
L0 = abs(u) <= eps;
if any(L0)
    t1(L0) = (A * (A + 2*B) / 3) + t0_t1;
    t2(L0) = (B * (B + 2*A) / 3) + t0_t2;
    p(L0) = B ./ (A + B);
end


% Calculate choices
p(~L0) = (exp(2*u(~L0)*B)-1) ./ (exp(2*u(~L0)*B)  - exp(-2*u(~L0)*A));
% Deal with nans at extreme values
p((isnan(p) & u > 0) | (p>1)) = 1;
p((isnan(p) & u < 0) | (p<0)) = 0;

% find point to split RTs either at subjective equality or smallest coh
if subj_eq
    % split at coh with pT1 closest to 0.5
    [foo, eq_pt] = min(abs(p-0.5));
else
    [foo, eq_pt] = min(abs(coh));
end

% indices to RTs we will actually calculate
if length(eq_pt) == 1
    % if there's one "equality" point, overlap T1 and T2 RTs there
    T1 = logical(coh >= coh(eq_pt));
    T2 = logical(coh <= coh(eq_pt));
elseif length(eq_pt) == 2
    % if there's two "equality" points, don't overlap T1 and T2 RTs 
    T1 = logical(coh > min(coh(eq_pt)));
    T2 = logical(coh < max(coh(eq_pt)));
else
    error('More than two points of equality')
end


% exclude limiting case from RT calculation
T1 = T1 & ~L0;
T2 = T2 & ~L0;
    
% Calculate RTs
t1(T1) = (( -B*coth(B*u(T1)) + (A+B)*coth((A+B)*u(T1)) ) ./ u(T1)) + t0_t1;
t2(T2) = (( -A*coth(A*u(T2)) + (A+B)*coth((A+B)*u(T2)) ) ./ u(T2)) + t0_t2;


 