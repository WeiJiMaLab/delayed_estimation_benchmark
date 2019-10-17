% function [fitpars, AIC, BIC, parnames] = fit_SA_model(data,modelflags)
%
% This function fits the SA model to a given dataset.
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

function [fitpars, AIC, BIC, parnames] = fit_SA_model(data,modelflags)

% init
gvar = get_gvar();
gvar.nPar = 4;                   % number of parameters (J1, kappa_r, Kpar, NT_slope)
parnames = {'J1','kappa_r','Kpar','NT_slope'};
if modelflags(2)~=3
    gvar.samplers = [2 2 4 2];   % 1=normal, 2=log normal, 3=poisson, 4=uniform on [X-2, X+2]
elseif modelflags(2)==4
    gvar.samplers = [2 2 2 2];  
end
gvar.ls = zeros(1,gvar.nPar)+0;

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
X_mat = zeros(gvar.popSize,gvar.nPar);
X_mat(:,1) = 5 + rand(gvar.popSize,1)*5;    % J1
X_mat(:,2) = 40 + rand(gvar.popSize,1)*20;  % kappa_r
X_mat(:,3) = randi(8,gvar.popSize,1);       % Kpar
X_mat(:,4) = rand(gvar.popSize,1)*.05;      % NT slope

% if K=Inf, set all J1 and Kpar samples to Inf 
if modelflags(2)==1
    X_mat(:,1) = Inf;
    X_mat(:,3) = Inf;
end
% if model without non-target responses, set all NT_slope samples to 0
if modelflags(3)==0
    X_mat(:,4) = 0;
end

% compute fitness of each member in initial generation
fprintf('\nFitting model...');
LLH = zeros(1,gvar.popSize);
for jj=1:gvar.popSize
    LLH(jj) = compute_LLH(X_mat(jj,:), data, gvar,modelflags);
end
fprintf('\n');

% Run genetic algorithm
tic
for ii=1:gvar.nGenerations

    % reproduce!
    X_mat_new = zeros(gvar.popSize*gvar.nKids,gvar.nPar);
    start_idx = 1;
    for jj=1:gvar.popSize
        
        for kk=1:gvar.nPar
            if modelflags(2)==1 && kk==1
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = Inf; % J1 in -A model is Inf
            elseif modelflags(2)==1 && kk==3
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = Inf; % K in -A models if Inf
            elseif modelflags(3)==0 && kk==4
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = 0; % NT slope in noNT models is 0                
            else
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = reproduce(X_mat(jj,kk),gvar.ls(kk),gvar.samplers(kk),gvar.nKids);
            end
        end
        start_idx = start_idx + gvar.nKids;
    end
    
    % compute fitness of each member in new generation
    LLH_new = zeros(1,gvar.popSize*gvar.nKids);
    for jj=1:gvar.popSize*gvar.nKids
        LLH_new(jj) = compute_LLH(X_mat_new(jj,:), data, gvar,modelflags);
    end
    
    % merge old and new
    X_mat_all = [X_mat; X_mat_new];
    LLH_all = [LLH LLH_new];

    % update population size
    gvar.popSize = max(gvar.popSizeMin,round(gvar.popSizeRate*gvar.popSize));
  
    % select fittest members 
    [LLH_all, I] = sort(LLH_all,'descend');   
    LLH = LLH_all(1:gvar.popSize);
    X_mat = X_mat_all(I(1:gvar.popSize),:);
%     fprintf('%d: mean/max LLH = %2.2f/%2.2f (ETL=%2.1f min)\n',ii,mean(LLH),max(LLH),toc/ii*(gvar.nGenerations-ii)/60);
    fprintf('.');
    
    gvar.ls = gvar.ls + gvar.dls;
end

% compute final results
if modelflags(2)==1
    nFreePars = 1;   % kappa_r (J=Inf, Kpar=Inf)
else
    nFreePars = 3;   % J, kappa_r, Kpar
end
if modelflags(3)==1
    nFreePars = nFreePars+1;  % NT_slope
end

% find ML parameter estimates
[maxLLH max_idx] = max(LLH);
fitpars = X_mat(max_idx,:);

% compute BIC and estimate marginal model likelihood
AIC = -2*maxLLH + nFreePars*2;
BIC = -2*maxLLH + nFreePars*log(numel(data.N));

%-%-%-%-%-%-%-%-%-%-%-%
% Likelihood function %
%-%-%-%-%-%-%-%-%-%-%-%
function LLH = compute_LLH(pars, data, gvar, modelflags)
J1 = pars(1);
kappa_r = pars(2);
Kpar = pars(3);
NT_slope = pars(4);
Kmax = 20;  % maximum K to compute predictions for (in Poisson model, K can become arbitrarily large)
if modelflags(2)==4
    poissW = poisspdf(0:Kmax,Kpar);
    poissW = poissW/sum(poissW);
end

LLH = 0;
if modelflags(2)==1
    K_range = Inf;
elseif modelflags(2)==2
    K_range = Kpar;
elseif modelflags(2)==3
    K_range = 0:Kmax;
else
    K_range = 0:Kpar;
end

if isinf(J1)
    k1 = Inf;
else
    k1 = interp1(gvar.J_map,gvar.kappa_map,J1,'linear','extrap');
end

unique_N = unique(data.N);
for ii=1:length(unique_N)
    N = unique_N(ii);
    
    if modelflags(2)>1
        % precompute error distribution under each possible value of K
        p_error_mat = zeros(length(K_range),length(gvar.error_range));
        for jj=1:length(K_range)
            if N<=K_range(jj)
                % Case N<=K
                J_high = J1*(floor(K_range(jj)/N)+1);
                J_low = J1*floor(K_range(jj)/N);
                kappa_high = min(interp1(gvar.J_map,gvar.kappa_map,J_high),gvar.kappa_max);
                kappa_low = min(interp1(gvar.J_map,gvar.kappa_map,J_low),gvar.kappa_max);
                p_high = mod(K_range(jj),N)/N;
                kc_high = sqrt(kappa_high^2 + kappa_r^2 + 2*kappa_high*kappa_r*cos(gvar.error_range));
                kc_low = sqrt(kappa_low^2 + kappa_r^2 + 2*kappa_low*kappa_r*cos(gvar.error_range));
                p_error_high = besseli0_fast(kc_high,1)./(2*pi*besseli0_fast(kappa_high,1).*besseli0_fast(kappa_r,1)) .* exp(kc_high-kappa_high-kappa_r);
                p_error_low = besseli0_fast(kc_low,1)./(2*pi*besseli0_fast(kappa_low,1).*besseli0_fast(kappa_r,1)) .* exp(kc_low-kappa_low-kappa_r);
                p_error_mat(jj,:) = p_high * p_error_high + (1-p_high) * p_error_low;
            else
                % Case N>K
                p_guess = 1-K_range(jj)/N;
                kc = sqrt(k1^2 + kappa_r^2 + 2*k1*kappa_r*cos(gvar.error_range));
                p_error_mat(jj,:) = p_guess*(1/2/pi) + (1-p_guess)*besseli0_fast(kc,1)./(2*pi*besseli0_fast(k1,1)*besseli0_fast(kappa_r,1)) .* exp(kc-k1-kappa_r);
            end
        end
    end
        
    % compute probability of errors under distribution of K imposed by current model
    if modelflags(2)==1
        % infinite K - response noise is the only source of error
        p_error = 1/2/pi/besseli0_fast(kappa_r,0) .* exp(kappa_r*cos(gvar.error_range));
    elseif modelflags(2)==2
        % fixed K
        p_error = p_error_mat(1,:);
    elseif modelflags(2)==3
        % variable item limit, K ~ Poisson(Kmean)
        w = poissW;
        p_error = zeros(1,length(gvar.error_range));
        for K=0:Kmax
            p_error = p_error + w(K+1)*p_error_mat(K+1,:);
        end        
    elseif modelflags(2)==4
        % variable item limit, K ~ U(0,Kmax)
        p_error = zeros(1,length(gvar.error_range));
        for K=0:Kpar
            w = 1/(Kpar+1);
            p_error = p_error + w*p_error_mat(K+1,:);
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
