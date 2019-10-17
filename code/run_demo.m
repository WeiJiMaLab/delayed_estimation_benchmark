% function run_demo()
%
% The main purpose of this demo is to serve as a reference that allows you
% to figure out how to generate synthetic data and how to fit models to data. 
%
% The demo will first ask you what model you want to use to generate a set
% of synthetic data. Possible model codes are: 
%
% FP-A, FP-F, FP-P, FP-U, FP-A-NT, FP-F-NT, FP-P-NT, FP-U-NT,
% SA-A, SA-F, SA-P, SA-U, SA-A-NT, SA-F-NT, SA-P-NT, SA-U-NT,
% EP-A, EP-F, EP-P, EP-U, EP-A-NT, EP-F-NT, EP-P-NT, EP-U-NT,
% VP-A, VP-F, VP-P, VP-U, VP-A-NT, VP-F-NT, VP-P-NT, VP-U-NT.
%
% Thereafter it will ask you which model you want to fit to the generated
% dataset.
% 
% After fitting the model, the maximum-likelihood parameter estimates, AIC, 
% and BIC estimates will be shown. In addition, a plot will be generated 
% that shows the estimation error distribution of the synthetic data and 
% model fit (in this plot, errors and fit are pooled across set sizes).
%
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2013.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com

function run_demo

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Ask user input (what model to use for generating data)  %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
fprintf('------------------------------------------------------\n');
fprintf('From which model would you like to generate fake data?\n');
fprintf('------------------------------------------------------\n');
fprintf('Possible model codes are: \n');
fprintf('FP-A, FP-F, FP-P, FP-U, FP-A-NT, FP-F-NT, FP-P-NT, FP-U-NT,\n');
fprintf('SA-A, SA-F, SA-P, SA-U, SA-A-NT, SA-F-NT, SA-P-NT, SA-U-NT,\n');
fprintf('EP-A, EP-F, EP-P, EP-U, EP-A-NT, EP-F-NT, EP-P-NT, EP-U-NT,\n');
fprintf('VP-A, VP-F, VP-P, VP-U, VP-A-NT, VP-F-NT, VP-P-NT, VP-U-NT.\n\n');

% get input
gen_modelflags=NaN;
while isnan(gen_modelflags)
    gen_modelstr = upper(input('Model code (to generate data): ','s'));
    gen_modelflags = get_modelflags(gen_modelstr);
end

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Create vector with parameter values for the specified model %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
if gen_modelflags(2)==1
    Kpar = Inf;  % set capacity parameter to infinite
    kname = 'K';
elseif gen_modelflags(2)==2
    Kpar = randi(3)+2; % set fixed capacity limit K to a random integer between 2 and 5
    kname = 'K';
elseif gen_modelflags(2)==3
    Kpar = rand*3+2; % set mean of Poisson distribution on K to a random value between 2 and 5
    kname = 'Kmean';
elseif gen_modelflags(2)==4
    Kpar = randi(6)+4; % set maximum of uniform distribution on K to a random integer between 5 and 10
    kname = 'Kmax';
end

if gen_modelflags(3)==0
    NT_slope = 0; % no non-target responses
else
    NT_slope = rand*.05; % set non-target response parameter to a random value between 0 and 0.05
end

kappa_r = rand*180+20; % draw random value for kappa_r from range [20,200]

if gen_modelflags(1)==1
    % FP model: kappa_r, Kpar, NT_slope
    pars_in = [kappa_r Kpar NT_slope];
    parnames_gen = {'kappa_r',kname,'NT_slope'};
elseif gen_modelflags(1)==2
    % SA model: J1, kappa_r, Kpar, NT_slope
    J1 = rand*5+2; 
    pars_in = [J1, kappa_r, Kpar, NT_slope];
    parnames_gen = {'J1','kappa_r',kname,'NT_slope'};
elseif gen_modelflags(1)==3
    % EP model: J1, power, tau=0, kappa_r, Kpar, NTslope
    J1 = rand*10+10; 
    power = -0.5-rand;
    tau = 0;
    pars_in = [J1, power, tau, kappa_r, Kpar, NT_slope];
    parnames_gen = {'J1','power','tau','kappa_r',kname,'NT_slope'};
elseif gen_modelflags(1)==4
    % VP model: J1, power, tau, kappa_r, Kpar, NTslope
    J1bar = rand*30+20; 
    power = -0.5-rand;
    tau = rand*10+10;
    pars_in = [J1bar, power, tau, kappa_r, Kpar, NT_slope];
    parnames_gen = {'J1','power','tau','kappa_r',kname,'NT_slope'};
end

% show input parameters
fprintf('\nGenerating data from %s model with the following randomly drawn parameter values:\n\n',gen_modelstr);
for ii=1:length(parnames_gen)
    fprintf('%s = %2.2f\n',parnames_gen{ii},pars_in(ii));
end

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Generate synthetic dataset  %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
n_trials = 2000; % number of trials to generate
N_vec = [1 2 3 4 6 8]; % set sizes included in synthetic data
if gen_modelflags(1)==1
    % FP model
    data = gen_fake_FP_data(pars_in,gen_modelflags(2),n_trials,N_vec);
elseif gen_modelflags(1)==2
    % SA model
    data = gen_fake_SA_data(pars_in,gen_modelflags(2),n_trials,N_vec);
elseif gen_modelflags(1)==3
    % EP model
    data = gen_fake_EPVP_data(pars_in,gen_modelflags(2),n_trials,N_vec);
elseif gen_modelflags(1)==4
    % VP model
    data = gen_fake_EPVP_data(pars_in,gen_modelflags(2),n_trials,N_vec);
end
fprintf('\nSynthetic data generated.\n');

%-%-%-%-%-%-%-%-%
% Fit model(s)  %
%-%-%-%-%-%-%-%-%
done=0;
while ~done
    fprintf('----------------------------------------------\n');
    fprintf('Which model would you like to fit to the data?\n');
    fprintf('----------------------------------------------\n');
    fprintf('Possible model codes are: \n');
    fprintf('FP-A, FP-F, FP-P, FP-U, FP-A-NT, FP-F-NT, FP-P-NT, FP-U-NT,\n');
    fprintf('SA-A, SA-F, SA-P, SA-U, SA-A-NT, SA-F-NT, SA-P-NT, SA-U-NT,\n');
    fprintf('EP-A, EP-F, EP-P, EP-U, EP-A-NT, EP-F-NT, EP-P-NT, EP-U-NT,\n');
    fprintf('VP-A, VP-F, VP-P, VP-U, VP-A-NT, VP-F-NT, VP-P-NT, VP-U-NT.\n\n');
    
    % get input
    fit_modelflags=NaN;
    while isnan(fit_modelflags)
        fit_modelstr = upper(input('Model code (to FIT data): ','s'));
        fit_modelflags = get_modelflags(fit_modelstr);
    end

    % fit model and generate a synthetic data using ML estimate (to get a fitted error distribution for plotting)
    if fit_modelflags(1)==1
        [fitpars, AIC, BIC, parnames_fit] = fit_FP_model(data,fit_modelflags);
        fprintf('\nDone. Computing predicted response distribution...\n');    
        data_fit = gen_fake_FP_data(fitpars,fit_modelflags(2),10000,N_vec);
    elseif fit_modelflags(1)==2
        [fitpars, AIC, BIC, parnames_fit] = fit_SA_model(data,fit_modelflags);
        fprintf('\nDone. Computing predicted response distribution...\n');    
        data_fit = gen_fake_SA_data(fitpars,fit_modelflags(2),10000,N_vec);
    elseif fit_modelflags(1)==3
        [fitpars, AIC, BIC, parnames_fit] = fit_EPVP_model(data,fit_modelflags);
        fprintf('\nDone. Computing predicted response distribution...\n');    
        data_fit = gen_fake_EPVP_data(fitpars,fit_modelflags(2),10000,N_vec);
    elseif fit_modelflags(1)==4
        [fitpars, AIC, BIC, parnames_fit] = fit_EPVP_model(data,fit_modelflags);
        fprintf('\nDone. Computing predicted response distribution...\n');    
        data_fit = gen_fake_EPVP_data(fitpars,fit_modelflags(2),10000,N_vec);
    end
    
    fprintf('\n-------------RESULTS---------------\n');
    % show input parameters
    fprintf('Data were generated from %s model with the following parameter values:\n',gen_modelstr);
    for ii=1:length(parnames_gen)
        fprintf(' %s = %2.2f\n',parnames_gen{ii},pars_in(ii));
    end
    fprintf('Results from fitting these data with %s model:\n',fit_modelstr);
    for ii=1:length(parnames_fit)
        fprintf(' %s = %2.2f\n',parnames_fit{ii},fitpars(ii));
    end
    fprintf('\nAIC=%2.2f, BIC=%2.2f\n',AIC,BIC);
    fprintf('-----------------------------------\n');
            
    % plot fit  
    figure
    X = linspace(-pi,pi,52);
    X = X(1:end-1)+diff(X(1:2))/2;
    Y_emp = hist(data.error_vec,X);
    Y_emp = Y_emp/sum(Y_emp)/diff(X(1:2));
    Y_fit = hist(data_fit.error_vec,X);
    Y_fit = Y_fit/sum(Y_fit)/diff(X(1:2));
    bar(X,Y_emp,'k');
    hold on
    plot(X,Y_fit,'r-','Linewidth',3)
    legend('Data','Fit')    
    xlabel('Response error');
    ylabel('Probability');
    xlim([-pi pi]);
    title([gen_modelstr ' data fitted with ' fit_modelstr ' model; AIC=' num2str(AIC,4)]);
    
    % ask whether to fit another model    
    another = input('\nFit another model to this dataset (Y/N)? ','s');
    fprintf('\n\n');
    done = ~strcmpi(another,'Y');    
end
%----------------- HELPER FUNCTIONS --------------------%

function modelflags = get_modelflags(modelstr)
% check length of string
if length(modelstr)<4
    fprintf('Invalid model code\n');
    modelflags=NaN;
    return;
end    

% set first model dimension
if strcmp(modelstr(1:3),'FP-')
    modelflags(1)=1; % fixed-precision model
elseif strcmp(modelstr(1:3),'SA-')
    modelflags(1)=2; % slots-plus-averaging model
elseif strcmp(modelstr(1:3),'EP-')
    modelflags(1)=3; % equal-precision model
elseif strcmp(modelstr(1:3),'VP-')
    modelflags(1)=4; % variable-precision model
else
    fprintf('Invalid model code\n');
    modelflags=NaN;
    return;
end

% set second model dimension
if strcmp(modelstr(4),'A')
    modelflags(2)=1; % model in which all items are remembered
elseif strcmp(modelstr(4),'F')
    modelflags(2)=2; % model in which a fixed number of items is remembered
elseif strcmp(modelstr(4),'P')
    modelflags(2)=3; % model with Poisson-distributed number of remembered items
elseif strcmp(modelstr(4),'U')
    modelflags(2)=4; % model with uniformly distributed number of remembered items
else
    fprintf('Invalid model code\n');
    modelflags=NaN;
    return;
end

% set third model dimension
if length(modelstr)==4
    modelflags(3)=0;  % model without non-target responses
elseif strcmp(modelstr(end-2:end),'-NT')
    modelflags(3)=1;  % model with non-target responses
else
    fprintf('Invalid model code\n');
    modelflags=NaN;
    return;
end

