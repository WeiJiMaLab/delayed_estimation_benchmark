% This function should do the same as Matlab's built-in randi.m function.
% (some older versions of Matlab don't seem to have this function)
%
% This code accompanies the paper "Conceptualizing and testing working 
% memory models in a three-dimensional model space" by Van den Berg, Awh,
% and Ma, published in Psychological Review, 2013.
%
% For questions/bug reports/etc, please email nronaldvdberg@gmail.com

function Y = randi(X,nr,nc)
if ~exist('nr','var')
    nr=1;
end
if ~exist('nc','var')
    nc=1;
end
Y = ceil(rand(nr,nc)*X);