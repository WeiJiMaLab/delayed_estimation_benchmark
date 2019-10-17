% function [fitpars, AIC, BIC, parnames] = fit_FP_model(data,modelflags)
%
% This function fits theF FP model to a given dataset.
%
% INPUT
%  data: see one of the gen_fake_*_data.m files for structuring of this
%        variable 
%  modelflags: [dim1 dim2 dim2], where dim1: 1=FP, 2=SA, 3=EP, 4=VP
%                                      dim2: 1=A,  2=F,  3=P,  4=U
%                                      dim3: 0=no NT, 1=NT
%
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2013.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com
function [fitpars, AIC, BIC, parnames] = fit_FP_model(data,modelflags)

% init
gvar = get_gvar();
gvar.n_par = 3;                 % number of parameters: kappa_r, Kpar, NT_slope
parnames = {'kappa_r','Kpar','NT_slope'};
if modelflags(2)~=3
    gvar.samplers = [2 4 2];   % 1=normal, 2=log normal, 3=poisson, 4=uniform on [X-2, X+2]
else
    gvar.samplers = [2 2 2];   
end
gvar.ls = zeros(1,gvar.n_par);

% precompute indices of target and nontarget errors (error space is discretized)
unique_N = unique(data.N);
for ii=1:length(unique_N)
    trial_idx = find(data.N==unique_N(ii));
    data.error_idx{ii} = interp1(gvar.error_range,1:length(gvar.error_range),abs(data.error_vec(trial_idx)),'nearest','extrap');
    data.dist_error_idx{ii} = [];
    if unique_N(ii)==1
        data.dist_error_idx{ii} = [];
    else
        for jj=1:length(trial_idx)
            data.dist_error_idx{ii}(:,jj) = interp1(gvar.error_range,1:length(gvar.error_range),abs(data.dist_error_vec{trial_idx(jj)}),'nearest','extrap');
        end
    end
end

% initialize matrix with parameter samples 
gvar.popSize = gvar.popSizeStart;
X_mat = zeros(gvar.popSize,gvar.n_par);
X_mat(:,1) = 1 + rand(gvar.popSize,1)*10;    % kappa_r
X_mat(:,2) = randi(10,gvar.popSize,1);       % Kpar
X_mat(:,3) = rand(gvar.popSize,1)*.05;       % NT slope

% if K=Inf, set all Kpar samples to Inf
if modelflags(2)==1
    X_mat(:,2) = Inf;
end
% if model without non-target responses, set all NT_slope samples to 0
if modelflags(3)==0
    X_mat(:,3) = 0;
end

% compute fitness of each member in initial generation 
fprintf('\nFitting model...');
LLH = zeros(1,gvar.popSize);
for jj=1:gvar.popSize
    LLH(jj) = compute_LLH(X_mat(jj,:), data, gvar, modelflags);
end

% Run genetic algorithm
tic
nEval=gvar.popSize;
for ii=1:gvar.nGenerations

    % reproduce!
    X_mat_new = zeros(gvar.popSize*gvar.nKids,gvar.n_par);
    start_idx = 1;
    for jj=1:gvar.popSize
        
        for kk=1:gvar.n_par
            if modelflags(2)==1 && kk==2
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = Inf;
            elseif modelflags(3)==0 && kk==3
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = 0;
            else
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = reproduce(X_mat(jj,kk),gvar.ls(kk),gvar.samplers(kk),gvar.nKids);
            end
        end
        start_idx = start_idx + gvar.nKids;
    end
    
    % compute fitness of each member in new generation
    LLH_new = zeros(1,gvar.popSize*gvar.nKids);
    for jj=1:gvar.popSize*gvar.nKids
        LLH_new(jj) = compute_LLH(X_mat_new(jj,:), data, gvar, modelflags);
        nEval=nEval+1;
    end
    
    % merge old and new 
    X_mat_all = [X_mat; X_mat_new];
    LLH_all = [LLH LLH_new];

    % update population size
    gvar.popSize = max(gvar.popSizeMin,round(gvar.popSizeRate*gvar.popSize));
    
    % update population by selecting fittest members
    [LLH_all I] = sort(LLH_all,'descend');    
    LLH = LLH_all(1:gvar.popSize);
    X_mat = X_mat_all(I(1:gvar.popSize),:);
%     fprintf('%d: mean/max LLH = %2.2f/%2.2f (ETL=%2.1f min)\n',ii,mean(LLH),max(LLH),toc/ii*(gvar.nGenerations-ii)/60);
    fprintf('.');

    gvar.ls = gvar.ls + gvar.dls;
end

% compute final results
n_free_pars = 1;   % kappa_r
if modelflags(2)>1
    n_free_pars = n_free_pars+1;  % Kpar
end
if modelflags(3)==1
    n_free_pars = n_free_pars+1;  % NT_slope
end

% find ML parameter estimates
[maxLLH max_idx] = max(LLH);
fitpars = X_mat(max_idx,:);

% compute BIC and estimate marginal model likelihood
AIC = -2*maxLLH + n_free_pars*2;
BIC = -2*maxLLH + n_free_pars*log(numel(data.N));

%-%-%-%-%-%-%-%-%-%-%-%
% Likelihood function %
%-%-%-%-%-%-%-%-%-%-%-%
function LLH = compute_LLH(pars, data, gvar, modelflags)
kappa_r = pars(1);
kappa_r = min(kappa_r,500);
Kpar = pars(2);
NT_slope = pars(3);

if modelflags(2)==3
    poissW = poisspdf(0:20,Kpar);
    poissW = poissW/sum(poissW);
end

LLH = 0;
unique_N = unique(data.N);
for ii=1:length(unique_N)  
    N = unique_N(ii);
    
    % compute error distribution under K=Inf (motor noise is the only source of error)
    p_error = 1/2/pi/besseli0_fast(kappa_r,0) .* exp(kappa_r*cos(gvar.error_range));
    
    % compute probability of errors under item limit
    if modelflags(2)==1
        % no item limit - don't need to do anything
    elseif modelflags(2)==2
        % fixed item limit 
        pEncode = min(Kpar,N)/N;
        p_error = pEncode*p_error + (1-pEncode)*1/2/pi;
    elseif modelflags(2)==3
        % variable item limit, K ~ Poisson(Kmean)
        w = poissW;
        w(N+1) = sum(w((N+1):end));  % all K>N give same prediction as K=N; sum the weights
        w = w(1:(N+1));
        p_error_original = p_error;
        p_error = p_error*0;
        for K=0:N
            pEncode = min(K,N)/N;
            p_error = p_error + w(K+1)*(pEncode*p_error_original + (1-pEncode)*1/2/pi);
        end        
    elseif modelflags(2)==4
        % variable item limit, K ~ U(0,Kmax)
        p_error_original = p_error;
        p_error = p_error*0;
        for K=0:Kpar
            pEncode = min(K,N)/N;
            w = 1/(Kpar+1);
            p_error = p_error + w*(pEncode*p_error_original + (1-pEncode)*1/2/pi);
        end
    end
    
    % make sure p_error integrates to 0.5 (we're considering only the positive half of the pdf)
    p_error = p_error/sum(p_error) * 1/diff(gvar.error_range(1:2))/2;

    % incorporate non-target responses    
    if modelflags(3)==0 || N==1
        p_resp = p_error(data.error_idx{ii});
    else
        nDist=N-1;
        p_NT = min(NT_slope*nDist,1);
        p_resp = (1-p_NT)*p_error(data.error_idx{ii});
        for jj=1:nDist
            p_resp = p_resp + p_NT/nDist*p_error(data.dist_error_idx{ii}(jj,:));
        end
    end
    
    % sum log likelihood of response errors over trials set sizes
    LLH = LLH + sum(log(p_resp));       
end

% this should never happen
if isnan(LLH)
    LLH = -Inf;
end
