% function data = gen_fake_SA_data(parameters,ktype,n_trials,N_vec)
%
% Generates synthetic data using the SA model.
%
% INPUT
%  ktype:
%   1: Infinite (K=N for each set size N)
%   2: Fixed    (K is a constant)
%   3: Uniform  (K is drawn from a uniform distribution on each trial)
%   4: Poisson  (K is drawn from a Poisson distribution on each trial)
%  modelpars: parameter vector [J1, kappa_r, Kpar, NT_slope
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

function data = gen_fake_SA_data(parameters,ktype,n_trials,N_vec)
% parse model pars
J1       = parameters(1);
kappa_r  = parameters(2);
Kpar     = parameters(3);
NT_slope = parameters(4);

% precompute mapping kappa <--> J
kmap = linspace(0,700,1e5);
Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);

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
        K = Inf;
    elseif ktype==2
        K = Kpar;
    elseif ktype==3
        K = poissrnd(Kpar);
    elseif ktype==4
        K = randi(Kpar+1)-1;  % draw K from uniform on [0,Kpar]
    else
        error('invalid ktype');
    end
    
    % decide which stimuli to encode        
    if K>=N
        p_high = mod(K,N)/N;
        if rand<p_high
            J_high = min(max(Jmap),J1*(floor(K/N)+1));
            kappa_high = interp1(Jmap,kmap,J_high,'linear','extrap');
            if isinf(K)
                x = stims(resp_idx); % encode with infinite precision
            else
                x = circ_vmrnd(stims(resp_idx),kappa_high); % encode with kappa_high
            end
        else
            J_low  = min(max(Jmap),J1*(floor(K/N)));
            kappa_low = interp1(Jmap,kmap,J_low,'linear','extrap');
            if isinf(K)
                x = stims(resp_idx); % encode with infinite precision
            else
                x = circ_vmrnd(stims(resp_idx),kappa_low); % encode with kappa_low
            end
        end
    else
        enc_idx = randperm(N);
        enc_idx = enc_idx(1:K);
        if sum(enc_idx==resp_idx)==0
            x = rand*2*pi-pi; % random guess
        elseif sum(enc_idx==resp_idx)==1
            kappa1 = interp1(Jmap,kmap,J1,'linear','extrap');
            x = circ_vmrnd(stims(resp_idx),kappa1); % encode with kappa corresponding to 1 chunk of resource
        end
    end
    
    % add response noise
    resp = circ_vmrnd(x,kappa_r);
    
    % compute error w.r.t. to target and distractors and store data
    data.error_vec(ii) = circ_dist(resp,stims(1));
    data.dist_error_vec{ii} = circ_dist(resp,stims(2:end));
    data.N(ii) = N;
end