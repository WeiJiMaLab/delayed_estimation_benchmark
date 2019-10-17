% This folder contains a collection of Matlab files that contain code to accompany
% the paper "Factorial comparison of working memory models" by Van den Berg, Awh, and Ma (Psychological Review, 2014).
%
% The code is quite extensively annotated. However, if you have any questions/bug reports/etc,
% please contact me at nronaldvdberg@gmail.com.
%
% To get started, I'd recommend to run and look at run_demo.m.


%--- FILE DESCRIPTIONS ---%

% run_demo.m     : this walks you through generating synthetic data and fitting models to those data
%
% fit_*_model.m  : used to fit the 32 main models described in the the paper
% gen_fake_*_data: used to generate synthetic data from the 32 main models described in the paper (note that these functions were (re)written for clarity, not speed; considerable speeding up can be achieved by vectorizing the code)
% 
% get_gvar.m     : return general settings for model fitting (mostly settings of the evolutionary optimization method)
% reproduce.m    : draw parameter values from a specific distribution (used by the optimization algorithm)
%
% besseli0_fast.m: evaluation I0(.) in a way faster than Matlab's besseli()
% 
% circ_*.m       : these files are part of the Circular Statistics Toolbox by P. Berens and J. Velasco
% randi.m        : should do the same as Matlab's randi.m (some older version of Matlab don't appear to have this function)


%--- DATA STRUCTURE ---%

% The fake data that are returned by gen_fake_*_data and the input data for fit_*_model are structures with 
% three fields: error_vec, dist_error_vec, and N. 
%
% N gives the set size on each trial.
%
% The error_vec field contains the errors on all trials, i.e., the (circular) difference between the response 
% and the target value. The domain is [-pi, pi], because that is the domain of the Von Mises distribution. 
% If your errors are for example in the range [-180, 180], then you should rescale them to [-pi, pi], 
% by multiplying them by pi/180.
%
% The dist_error_vec contains the "non-target errors" on each trial, i.e., the (circular) difference between 
% the response and each non-target value (again in the range [-pi, pi]). Say, for example, that we have N=4 on 
% the first trial. Then dist_error_vec{1} contains 3 values, corresponding to the distances between the response 
% and the first, second, and third non-target item, respectively. This information is only required when you 
% fit models that include non-target responses, i.e., all the XXX-NT models.
