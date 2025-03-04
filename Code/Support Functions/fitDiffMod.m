function err = fitDiffMod(theta, data, flags)
% err = fitDiffMod(theta) returns minus log likelihood of obtaining y
% values in last three columns of data given a fit parameterizeed by
% theta. 
% 
% 
% Arguments:
%  - theta      ... row vector of parameters of the model [k u0 A S.P./B T_res]
%                   where the number of each must match necessary number
%                   given flags and S.P./B refers to the starting point
%                   offset or B barrier depending on flags
%  - data       ... matrix with columns 
%                   [coh T1_rt_mean T1_rt_se T2_rt_mean T2_rt_se pT1 n] with last six
%                   entries repeated for each condition (e.g., ustim or
%                   bias)
%  - flags      ... structure that specifies the type of fit and number of conditions
%     .conds    ... number of conditions (e.g., 2 ustim conds or 3 bias ones)
%     .type     ... 0 to fit RT and choice, 1 for just RT, and 2 for just choice
%     .u0       ... -1 don't include drift offsets, 0 don't include them
%                   for 1st cond, 1 include them for all conds
%     .barr     ... -1 don't allow S.P./B changes, 0 allow S.P. changes for
%                   N+1 conds, 1 allow S.P. changes for all conds, 2 allow B 
%                   changes for N+1 conds, 3 allow B changes for all conds
%     .resid    ... 1 for same T1/T2 residual time, 2 for different
%     .nT1      ... 1 to pass nT1 rather than pT1, 0 for usual method
%     .subj_eq  ... 0 to split at 0% coh, 1 to split where pT1 = 0.5  
%     .log_trans... NOT YET IMPLEMENTED
%c
% Returns:
%  - err        ... minus log likelihood of obtaining values in last N-1
%                   columns of data given a fit parameterized by theta
%
% 
% 

% 7-05 Created by TDH


%% Check and set flags
if nargin < 3
    flags = [];
end
% set nonexistent fields to defaults
if ~isfield(flags, 'conds')
    flags.conds = 1;
end
if ~isfield(flags, 'type')
    flags.type = 0;
end
if ~isfield(flags, 'u0')
    flags.u0 = 0;
end
if ~isfield(flags, 'barr')
    flags.barr = 0;
end
if ~isfield(flags, 'resid')
    flags.resid = 1;
end
if ~isfield(flags, 'nT1')
    flags.nT1 = 0;
end
if ~isfield(flags, 'subj_eq')
    flags.subj_eq = 0;
end


%% Extract parameter values
i = 1;  % index of where we are in theta
% coh-to-drift conversion
k = theta(i);
% drift offset
if flags.u0 == -1
    u0 = 0;
else
    if (flags.conds == 1 & flags.u0 == 0)
        error('Not valid flag setting: one condition, drift offset changes on n+1 conditions')
    end
    u0 = zeros(1,flags.conds);
    n_u0 = flags.conds - 1 + flags.u0;
    i = i+1;
    j = i+n_u0-1;
    u0(end-n_u0+1:end) = theta(i:j);
    i = j;
end
% barrier heights
if flags.barr == -1
    i = i+1;
    A = ones(1,flags.conds)*theta(i);
    B = A;
elseif flags.barr < 2
    if (flags.conds == 1 & flags.barr == 0)
        error('Not valid flag setting: one condition, starting point changes on n+1 conditions')
    end
    i = i+1;
    A = ones(1,flags.conds)*theta(i);
    B = A;
    n_SP = flags.conds - 1 + flags.barr;
    i = i+1;
    j = i+n_SP-1;
    % starting point offset implemented as symmetric changes to both bounds
    A(end-n_SP+1:end) = A(end-n_SP+1:end) - theta(i:j);
    B(end-n_SP+1:end) = B(end-n_SP+1:end) + theta(i:j);
    i = j;
else
    if (flags.conds == 1 & flags.barr == 2)
        error('Not valid flag setting: one condition, B barrier changes on n+1 conditions')
    end
    n_A = flag.conds;
    i = i+1;
    j = i+n_A;
    A = theta(i:j);
    i = j;
    n_B = flags.conds - 3 + flags.barr;
    % this will be overwritten if flags.barr == 3
    B = A;
    j = i+n_B-1;
    B(end-n_B+1:end) = theta(i:j);
    i = j;
end
% residual times
i = i+1;
t1_res = theta(i);
if flags.resid == 2
    i = i+1;
    t2_res = theta(i);
else
    t2_res = t1_res;
end


% Extract coherence
coh = data(:,1);

%% Calculate predictions and minus summed log likelihood

err = 0;
% Loop through each condition
for i=1:flags.conds
    
    % Assign appropriate values to theta_cond
    theta_cond = [k u0(i) A(i) B(i) t1_res t2_res];
    
    
    %% Calculate predicted RT and pT1
    [t1pred, t2pred, ppred] = calcDiffMod(coh, theta_cond, flags.subj_eq);
    
    %% Calculate minus summed log likelihood of data based on predictions
    
    % Extract data for current condition
    offset = 1 + (i-1)*6;
    t1obs   = data(:,offset+1);
    t1se    = data(:,offset+2);
    t2obs   = data(:,offset+3);
    t2se    = data(:,offset+4);
    np_obs  = data(:,offset+5);
    n       = data(:,offset+6);
        
    % use RTs for combined or RT only fits
    if flags.type ~= 2
        % find indices where predictions were made
        t1_ind = ~isnan(t1pred) & ~isnan(t1se);
        t2_ind = ~isnan(t2pred) & ~isnan(t2se);
        % add minus log likelihood of RT data based on Gaussian error
        err = err + sum(0.5*(((t1obs(t1_ind)-t1pred(t1_ind))./t1se(t1_ind)).^2))  + sum(log(t1se(t1_ind))) + 0.5*sum(t1_ind)*log(2*pi);
        err = err + sum(0.5*(((t2obs(t2_ind)-t2pred(t2_ind))./t2se(t2_ind)).^2))  + sum(log(t2se(t2_ind))) + 0.5*sum(t2_ind)*log(2*pi);
        % same as: err - nansum(log(normpdf(tobs,tpred,tse)))
    end
    % use choices for combined or choice only fits
    if flags.type ~= 1
        % estimate n_obs if p_obs is given
        if flags.nT1
            n_obs = np_obs;
        else
            n_obs = round(np_obs.*n);
        end
        % add minus log likelihood based on binomial distribution
        % add 'eps' before taking log(0), which results in err = Nan
        n_diff = n-n_obs+1;
        n_diff(n_diff<0) = 0;
        err = err - sum(gammaln(n+1) - gammaln(n_obs+1) - gammaln(n_diff) + n_obs.*log(ppred+eps) + (n-n_obs).*log(1-ppred+eps));
        % same as: err - sum(log(binopdf(n_obs,n,ppred)))
    end
end
%disp(num2str(err))
