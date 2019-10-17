% function [fitpars, AIC, BIC, parnames] = fit_EPVP_model(data,modelflags)
%
% This function fits the EP or VP model to a given dataset.
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

function [fitpars, AIC, BIC, parnames] = fit_EPVP_model(data,modelflags)

% init
gvar = get_gvar();
gvar.n_par = 6;       % number of parameters (J1bar, power, tau, kappa_r, Kpar, NTslope)
if modelflags(2)~=3
    gvar.samplers = [2 1 2 2 4 2];   % 1=normal, 2=log normal, 3=poisson, 4=uniform on [X-2, X+2]
elseif modelflags(2)==3
    gvar.samplers = [2 1 2 2 2 2]; 
end
gvar.ls = zeros(1,gvar.n_par);

% get indices of errors
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

% initialize matrix with parameter samples with random values
gvar.popSize = gvar.popSizeStart;
X_mat = zeros(gvar.popSize,gvar.n_par);
X_mat(:,1) = 40 + rand(gvar.popSize,1)*20; % J1bar
X_mat(:,2) = -2 + rand(gvar.popSize,1);    % power
X_mat(:,3) = 20 + rand(gvar.popSize,1)*10; % tau
X_mat(:,4) = 40 + rand(gvar.popSize,1)*20; % kappa_r
X_mat(:,5) = randi(15,gvar.popSize,1);     % Kpar
X_mat(:,6) = rand(gvar.popSize,1)*.05;     % NT slope

% if EP, set all tau samples to 0
if modelflags(1)==3
    X_mat(:,3) = 0;
end
% if K=Inf, set all Kpar samples to Inf
if modelflags(2)==1
    X_mat(:,5) = Inf;
end
% if model without non-target responses, set all NT_slope samples to 0
if modelflags(3)==0
    X_mat(:,6) = 0;
end

% compute fitness of each member in initial generation 
fprintf('\nFitting model...');
LLH = zeros(1,gvar.popSize);
for jj=1:gvar.popSize
    LLH(jj) = compute_LLH(X_mat(jj,:), data, gvar, modelflags);
end
fprintf('\n');

% Run genetic algorithm
tic
for ii=1:gvar.nGenerations

    % reproduce!
    X_mat_new = zeros(gvar.popSize*gvar.nKids,gvar.n_par);
    start_idx = 1;
    for jj=1:gvar.popSize
        
        for kk=1:gvar.n_par
            if modelflags(1)==3 && kk==3
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = 0;   % tau in EP is 0
            elseif modelflags(2)==1 && kk==5
                X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = Inf; % K in -A models if Inf
            elseif modelflags(3)==0 && kk==6
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
        LLH_new(jj) = compute_LLH(X_mat_new(jj,:), data, gvar, modelflags);
    end
    
    % merge old and new 
    X_mat_all = [X_mat; X_mat_new];
    LLH_all = [LLH LLH_new];

    % update population size
    gvar.popSize = max(gvar.popSizeMin,round(gvar.popSizeRate*gvar.popSize));

    % select fittest members (keep population size constant)
    [LLH_all I] = sort(LLH_all,'descend');    
    LLH = LLH_all(1:gvar.popSize);
    X_mat = X_mat_all(I(1:gvar.popSize),:);
    
    % print some progress information
%     fprintf('%d: mean/max LLH = %2.2f/%2.2f (ETL=%2.1f min)\n',ii,mean(LLH),max(LLH),toc/ii*(gvar.nGenerations-ii)/60);
%     fprintf('    %2.1f, %2.1f, %2.1f, %2.1f, %2.1f, %2.2f\n',X_mat(1,1),X_mat(1,2),X_mat(1,3),X_mat(1,4),X_mat(1,5),X_mat(1,6));
    fprintf('.');
    
    gvar.ls = gvar.ls + gvar.dls;
end

% compute final results
n_free_pars = 3;   % Jbar, power, kappa_r
if modelflags(1)==4
    n_free_pars = n_free_pars+1;  % tau
end
if modelflags(2)>1
    n_free_pars = n_free_pars+1;  % Kpar
end
if modelflags(3)==1
    n_free_pars = n_free_pars+1;  % NT_slope
end

% find ML parameter estimates
maxLLH = max(LLH);
fitpars = X_mat(LLH==maxLLH,:);

% compute AIC and BIC
BIC = -2*maxLLH + n_free_pars*log(numel(data.N));
AIC = -2*maxLLH + 2*n_free_pars;

if modelflags(1)==3
    parnames = {'J1','power','tau','kappa_r','Kpar','NTslope'};
else
    parnames = {'J1bar','power','tau','kappa_r','Kpar','NTslope'};
end

%-%-%-%-%-%-%-%-%-%-%-%
% Likelihood function %
%-%-%-%-%-%-%-%-%-%-%-%
function LLH = compute_LLH(pars, data, gvar, modelflags)
J1bar = pars(1);
power = pars(2);
tau = pars(3);
kappa_r = pars(4);
Kpar = pars(5);
NT_slope = pars(6);

if modelflags(2)==1           % no limit; we need predictions for K=unique(data.N) encoded items
    K_range = unique(data.N);  
elseif modelflags(2)==2       % fixed limit; we need predictions for K=unique(data.N) and for K=Kpar encoded items, but not for K>N or K>Kpar
    K_range = unique([data.N Kpar]);            
    K_range = K_range(K_range<=max(data.N));
    K_range = K_range(K_range<=Kpar);
elseif modelflags(2)==3       % uniformly distributed K; we need predictions for K=0...min(Kpar,maxN)
    K_range = 0:max(data.N);   
    poissW = poisspdf(0:(max(K_range)-1),Kpar);
    poissW(end+1) = 1-sum(poissW);    
elseif modelflags(2)==4       % Poisson distributed K; we need predictions for K=0...max(N)
    K_range = 0:min(Kpar,max(data.N));
end
LLH = 0;

% (pre)compute predictions for each possible number of encoded items
p_error_mat = zeros(length(K_range),length(gvar.error_range));
for ii=1:length(K_range)
    if K_range(ii)==0
        % 0 items encoded --> response distribution is flat
        p_error_mat(ii,:) = 1/2/pi;
    else
        % draw values for J and compute corresponding kappas
        if tau==0
            J = J1bar*K_range(ii)^power;
        else
            Jbar = ones(1,gvar.nMCSamples)*J1bar*K_range(ii)^power;
            J = gamrnd(Jbar/tau,tau);
        end
        J = min(J,max(gvar.J_map));
        kappa = interp1(gvar.J_map,gvar.kappa_map,J);
        
        % compute response distribution (add motor noise, marginalize over J)
        kappa_sq = kappa.^2;
        k_c = zeros(length(gvar.error_range),numel(J));
        for jj=1:length(gvar.error_range)
            k_c(jj,:) = sqrt(kappa_r^2 + kappa_sq + 2*kappa_r*kappa*cos(gvar.error_range(jj)));
        end
        p_error_mat(ii,:) = mean(bsxfun(@rdivide,besseli0_fast(k_c,1),2*pi*besseli0_fast(kappa,1)*besseli0_fast(kappa_r,1)).*exp(bsxfun(@minus,k_c,kappa+kappa_r)),2);        
    end
end
unique_N = unique(data.N);
for ii=1:length(unique_N)       
    N = unique_N(ii);
    p_error = zeros(1,length(gvar.error_range));
        
    % compute probability of errors under item limit
    if modelflags(2)==1
        % no item limit 
        p_error = p_error_mat(K_range==N,:);
    elseif modelflags(2)==2
        % fixed item limit 
        pEncode = min(Kpar,N)/N;
        p_error = pEncode*p_error_mat(K_range==min(Kpar,N),:) + (1-pEncode)*1/2/pi;
    elseif modelflags(2)==3
        % variable item limit, K ~ Poisson(Kmean)   
        for jj=1:length(K_range)
            if K_range(jj)<N                
                pEncode = K_range(jj)/N;
                p_error = p_error + poissW(jj)*(pEncode*p_error_mat(jj,:) + (1-pEncode)*1/2/pi);
            else
                p_error = p_error + poissW(jj)*p_error_mat(K_range==N,:);
            end
        end        
    elseif modelflags(2)==4
        % variable item limit, K ~ U(0,Kmax)
        w_vec = ones(1,Kpar+1)/(Kpar+1);    % set weight to 1/(Kpar+1) for all possible K
        w_vec = w_vec(1:min(Kpar+1,N+1));   % predictions for K>N are all equal to prediction for K=N 
        w_vec(end) = 1-sum(w_vec(1:end-1));
        for jj=1:length(w_vec)
            K=jj-1;
            if K<N
                pEncode = K/N;
                p_error = p_error + w_vec(jj)*(pEncode*p_error_mat(jj,:) + (1-pEncode)*1/2/pi);
            else
                p_error = p_error + w_vec(jj)*p_error_mat(K_range==N,:);
            end
        end        
    end
       
    % make sure p_error integrates to 0.5 (we're considering only the positive half of the pdf)
    p_error = p_error/sum(p_error) * 1/diff(gvar.error_range(1:2))/2;

    % incorporate non-target responses and compute probability of response
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
        
    LLH = LLH + sum(log(p_resp));       
end

% this should never happen
if isnan(LLH)
    LLH = -Inf;
end