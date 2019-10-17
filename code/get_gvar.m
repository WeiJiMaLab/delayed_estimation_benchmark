% function gvar = get_gvar()
% 
% This function returns general settings that are used for all models. 
% Note that to speed things up, the settings differ from those used to 
% generate the results that are reported in the paper.
% 
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2013.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com

function gvar = get_gvar()

% discretization of the error space
error_range = linspace(0,pi,91); 
gvar.error_range = error_range(1:end-1)+diff(error_range(1:2))/2;

% mapping between J and kappa
gvar.kappa_max      = 10000;
gvar.kappa_map      = [linspace(0,10,250) linspace(10.001,gvar.kappa_max,250)];
gvar.J_map          = gvar.kappa_map.*besseli(1,gvar.kappa_map,1)./besseli(0,gvar.kappa_map,1);

% settings of evolutionary algorithm used for optimization
gvar.popSizeStart   = 128;   % size of population at start (Paper: 512)
gvar.popSizeMin     = 32;    % minimum population size 
gvar.popSizeRate    = .98;   % size of new population is popSizeRate times size of current population (unless it is already at minimum)  
gvar.nKids          = 1;     % number of kids per member 
gvar.nGenerations   = 64;    % number of generations to simulate (Paper: 256)
gvar.dls            = -0.01; % difference in ls from generation to generation 
gvar.nMCSamples     = 200;   % number of MC samples to draw when computing model predictions (Paper: 1000)
