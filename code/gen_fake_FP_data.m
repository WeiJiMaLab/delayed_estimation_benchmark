% function data = gen_fake_FP_data(parameters,ktype,n_trials,N_vec)
%
% Generates synthetic data using the FP model.
%
% INPUT
%  ktype:
%   1: Infinite (K=N for each set size N)
%   2: Fixed    (K is a constant)
%   3: Uniform  (K is drawn from a uniform distribution on each trial)
%   4: Poisson  (K is drawn from a Poisson distribution on each trial)
%  modelpars: parameter vector [kappa_r, Kpar, NT_slope]
%  n_trials : total number of trials to generate
%  n_vec    : set sizes for which to produce data (uniformly sampled)
%
% OUTPUT
%  data: set of synthetic data
%
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2013.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com

function data = gen_fake_FP_data(parameters,ktype,n_trials,N_vec)

% parse model pars
kappa_r  = parameters(1);
Kpar     = parameters(2);
NT_slope = parameters(3);

% loop over trials for clarity (this can be speeded up considerably by vectorizing it)
for ii=1:n_trials
    % draw random set size for this trial
    N = N_vec(randi(numel(N_vec)));
    
    % draw stimuli for this trial
    stims = rand(1,N)*2*pi-pi;
    
    % compute probability to report target
    p_target = 1-(N-1)*NT_slope;
    
    % decide which stimulus to report (if in set of encoded items)
    if rand<p_target
        resp_idx = 1; % report target
    else
        resp_idx = randi(N-1)+1; % report one of the distractors
    end

    % decide how many stimuli to encode
    if ktype==1
        K = Inf; % no capacity limit
    elseif ktype==2 
        K = Kpar; % fixed capacity limit
    elseif ktype==3
        K = poissrnd(Kpar); % draw K from Poisson distribution
    elseif ktype==4
        K = randi(Kpar+1)-1;  % draw K from uniform on [0,Kpar]
    else
        error('invalid ktype');
    end
    
    % decide which stimuli to encode
    if K>=N
        x = stims(resp_idx);
    else
        enc_idx = randperm(N);
        enc_idx = enc_idx(1:K);
        if sum(enc_idx==resp_idx)==0
            x = rand*2*pi-pi;
        elseif sum(enc_idx==resp_idx)==1
            x = stims(resp_idx);
        end
    end
    
    % add response noise
    resp = circ_vmrnd(x,kappa_r);
            
    % compute error w.r.t. to target and distractors and store data
    data.error_vec(ii) = circ_dist(resp,stims(1));
    data.dist_error_vec{ii} = circ_dist(resp,stims(2:end));
    data.N(ii) = N;
end
