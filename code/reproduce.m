% function X_new = reproduce(X_old,ls,sampler_id,n_samples)
%
% This function is used by the evolutionary optimization algorithm to draw parameter values
% from certain distributions.
%
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2013.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com

function X_new = reproduce(X_old,ls,sampler_id,n_samples)
if sampler_id==1
    % normal distribution
    X_new = normrnd(X_old,exp(ls),1,n_samples);
elseif sampler_id==2
    % log normal distribution
    X_new = exp(normrnd(log(X_old),exp(ls),1,n_samples));
elseif sampler_id==3
    % Poisson distribution
    X_new = poissrnd(X_old,1,n_samples);    
elseif sampler_id==4
    % Uniform distribution on [X_old-1, X_old+1]
    X_new = max(X_old + randi(5,1,n_samples)-3,0);
end

