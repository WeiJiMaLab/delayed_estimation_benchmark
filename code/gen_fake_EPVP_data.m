% function data  = gen_fake_EPVP_data(parameters,ktype,n_trials,N_vec)
%
% Generates synthetic data using the EP or VP model. To generate EP data,
% set parameter tau (i.e., modelpars(2)) to 0. 
%
% INPUT
%  ktype:
%   1: Infinite (K=N for each set size N)
%   2: Fixed    (K is a constant)
%   3: Uniform  (K is drawn from a uniform distribution on each trial)
%   4: Poisson  (K is drawn from a Poisson distribution on each trial)
%  modelpars: parameter vector [J1bar power tau kappa_r Kpar NT_slope]
%  n_trials : total number of trials to generate
%  n_vec    : set sizes for which to produce data (uniformly sampled)
%
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2013.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com

function data  = gen_fake_EPVP_data(parameters,ktype,n_trials,N_vec)

% parse model pars & flags
J1bar    = parameters(1);
power    = parameters(2);
tau      = parameters(3);
kappa_r  = parameters(4);
Kpar     = parameters(5);
NT_slope = parameters(6);

% precompute mapping between J and kappa
kmap = linspace(0,700,1e5);
Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);

% loop over trials for clarity (this can be speeded up considerably by vectorizing it)
for ii=1:n_trials
    
    % draw a random set size
    N = N_vec(randi(numel(N_vec)));
    
    % compute probability to report target (1 minus probability to report a nontarget)
    p_target = 1-(N-1)*NT_slope;
    
    % draw stimuli 
    stims = rand(1,N)*2*pi-pi;
    
    % decide which stimulus to report 
    if rand<p_target
        resp_idx = 1; % report target
    else
        resp_idx = randi(N-1)+1; % report one of the distractors
    end

    % decide how many stimuli to encode
    if ktype==1
        K = N; % everything is encoded
    elseif ktype==2
        K = Kpar; % fixed number of items are encoded
    elseif ktype==3
        K = poissrnd(Kpar); % draw K from Poisson pdf
    elseif ktype==4
        K = randi(Kpar+1)-1; % draw K from uniform on [0,Kpar]
    else
        error('invalid ktype');
    end

    % draw encoding precision of reported item and compute corresponding kappa
    Jbar = J1bar*min(K,N)^power; 
    if tau==0
        J = Jbar;  % EP
    else
        J = gamrnd(Jbar/tau,tau); % VP
    end
    J = min(J,max(Jmap));
    kappa = interp1(Jmap,kmap,J,'linear','extrap');
    
    % draw a measurement (target+noise or nontarget+noise)
    if K>=N
        % report an encoded item
        x = circ_vmrnd(stims(resp_idx),kappa);
    else
        % produce a guess if item not encoded
        enc_idx = randperm(N);
        enc_idx = enc_idx(1:K);
        if sum(enc_idx==resp_idx)==0
            x = rand*2*pi-pi; % report random guess
        elseif sum(enc_idx==resp_idx)==1
            x = circ_vmrnd(stims(resp_idx),kappa); % report item+noise
        end
    end
    
    % add response noise
    resp = circ_vmrnd(x,kappa_r); 
            
    % compute errors w.r.t. to target and nontargets and store everything
    data.error_vec(ii) = circ_dist(resp,stims(1));
    data.dist_error_vec{ii} = circ_dist(resp,stims(2:end));       
    data.N(ii) = N;
end