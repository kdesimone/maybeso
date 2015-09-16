function [ param, pmat, nll ] = fitgddloop( dspace, model, verbose )

% FITGDDLOOP  Fit generalized double detection function or one of its special
%             cases to a proxy decison space; make several fits from random
%             starting points, and choose the best fit
% 
%     usage:  [ param, pmat ] = fitgddloop( dspace, model, verbose )
% 
%     See the help text for fitgdd.m for an explanation of the input and
%     output arguments.
% 
%     The one input argument specific to fitgddloop.m is 'verbose', which
%     flags whether the function should print occasional update messages
%     to the command window so that the user can see progress.

% set default argument
if nargin<3, verbose = 1; end

% choose number of tries
ntries = 20;

% make several tries
errmin = Inf;
for i = 1:ntries
    
    % show progress
    if verbose
        fprintf('    fit %d of %d ...\n',i,ntries);
    end
    
    % make a fit
    [ param_tmp, pmat_tmp, err ] = fitgdd( dspace, model );
    
    % if it's the best fit so far, keep it
    if err<errmin
        param = param_tmp;
        pmat = pmat_tmp;
        errmin = err;
    end
    
end

% return negative log likelihood of the fit
nll = errmin;

end
