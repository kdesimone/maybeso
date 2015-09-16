function [ param, pmat, err ] = fitgdd( dspace, model )

% FITGDD  Fit generalized double detection function or one of its special
%         cases to a proxy decison space
% 
%     usage:  [ param, pmat, err ] = fitgdd( dspace, model )
% 
%     input arguments
%         'dspace' is a struct that describes the proxy decision space; it is the return argument from calcdspace.m
%         'model' is a string describing the model to be fitted to the proxy decision space; the options are
%             'gdd', the GDD function described in Pritchett and Murray (2015)
%             'single1', the single-interval model that only uses stimulus interval 1
%             'single2', the single-interval model that only uses stimulus interval 2
%             'diff',  the difference model
%             'diffg', the difference model with guessing
%             'ddetect', the double detection model
%             'oneline', the difference model with the added flexibility that the decision line can take any orientation
% 
%     return argument
%         'param' is a 1 x 5 matrix of fitted parameters; in terms of Pritchett and Murray's (2015) equation (1), these are
%              param(1) = theta1, the orientation of decision line 1
%              param(2) = theta2, the orientation of decision line 2
%              param(3) = delta1, the position of decision line 1
%              param(4) = delta2, the position of decision line 2
%              param(5) = sigma, the blur parameter explained in the paragraph following Pritchett and Murray's (2015) equation (7)
%         'pmat' is a matrix that shows the fitted proxy decision space; it gives the probability of the observer giving response 2 for each possible value of the proxy decision variable, according to the fitted model
%         'err' is the fit error of param

% seed random number generators
c=clock;
rand('state',round( sum(c) + (1e+8)*(c(6)-floor(c(6))) ));
randn('state',round( sum(c) + (2e+8)*(c(6)-floor(c(6))) ));

% construct a function to fit to the proxy decision space, and also random
% initial values for the parameters; see Pritchett and Murray's (2015) 
% methods section for an explanation of how all these models are special
% cases of the generalized double detection function.
switch model
    case 'gdd'      % full GDD function; fit all parameters
        fitfn = @( dvar, p ) max(min( gddrmat( dvar, p(1), p(2), p(3), p(4), 0.5, p(5) ), 0.99),0.01);
        pinit = [ 2*pi*rand(1,2) -0.5+rand(1,2) 0.1+0.4*rand ];
    case 'single1'  % single decision line, vertical
        fitfn = @( dvar, p ) max(min( gddrmat( dvar, pi, pi, p(1), p(1), 0.5, p(2) ), 0.99),0.01);
        pinit = [ -0.5+rand 0.1+0.4*rand ];
    case 'single2'  % single decision line, horizontal
        fitfn = @( dvar, p ) max(min( gddrmat( dvar, pi/2, pi/2, p(1), p(1), 0.5, p(2) ), 0.99),0.01);
        pinit = [ -0.5+rand 0.1+0.4*rand ];
    case 'diff'     % difference model
        fitfn = @( dvar, p ) max(min( gddrmat( dvar, 3*pi/4, 3*pi/4, p(1), p(1), 0.5, p(2) ), 0.99),0.01);
        pinit = [ -0.5+rand 0.1+0.4*rand ];
    case 'diffg'    % difference model with guessing
        fitfn = @( dvar, p ) max(min( gddrmat( dvar, 3*pi/4, 3*pi/4, p(1), p(2), 0.5, p(3) ), 0.99),0.01);
        pinit = [ -0.5+rand(1,2) 0.1+0.4*rand ];
    case 'ddetect'  % double detection model
        fitfn = @( dvar, p ) max(min( gddrmat( dvar, pi, pi/2, p(1), p(2), 0.5, p(3) ), 0.99),0.01);
        pinit = [ -0.5+rand(1,2) 0.1+0.4*rand ];
    case 'oneline'  % single decision line at any position and orientation
        fitfn = @( dvar, p ) max(min( gddrmat( dvar, p(1), p(1), p(2), p(2), 0.5, p(3) ), 0.99),0.01);
        pinit = [ 2*pi*rand -0.5+rand 0.1+0.4*rand ];
    otherwise
        error('unknown model type ''%s''',model);
end

% shift the decision variable range so that it's approximately centered on
% the origin.  this does not change the fitting problem in principle,
% but in practice it makes the fits converge more reliably.
dmean = mean(dspace.dlist);
dspace.dlist = dspace.dlist - dmean;

% construct the maximum likelihood objective function
errfn = @( p ) -sum(sum(log( binopdf( dspace.kmat, dspace.nmat, fitfn(dspace.dlist,p) ) )));

% fit parameters via simulated annealing
[ pfit, err ] = fminsearch( errfn, pinit );

% get response matrix from fitted parameters
pmat = fitfn( dspace.dlist, pfit );

% create parameter vector
switch model
    case 'gdd'      % full GDD function; fit all parameters
        param = pfit;
    case 'single1'  % single decision line, vertical
        param = [ pi pi pfit(1) pfit(1) pfit(2) ];
    case 'single2'  % single decision line, horizontal
        param = [ pi/2 pi/2 pfit(1) pfit(1) pfit(2) ];
    case 'diff'     % difference model
        param = [ 3*pi/4 3*pi/4 pfit(1) pfit(1) pfit(2) ];
    case 'diffg'    % difference model with guessing
        param = [ 3*pi/4 3*pi/4 pfit(1) pfit(2) pfit(3) ];
    case 'ddetect'  % double detection model
        param = [ pi pi/2 pfit(1) pfit(2) pfit(3) ];
    case 'oneline'  % single decision line at any position and orientation
        param = [ pfit(1) pfit(1) pfit(2) pfit(2) pfit(3) ];
end

% above (lines 62-63) we shifted the decision variable range so that
% it's approximately centered on the origin.  this didn't change the
% decision variable range in the calling function, of course, so we 
% need to adjust the fitted parameters to undo the effect of this shift.
param(3) = param(3) + dmean*( cos(param(1))+sin(param(1)) );
param(4) = param(4) + dmean*( cos(param(2))+sin(param(2)) );
% explanation:  suppose we have coordinate system A with coordinates (x,y)
% and coordinate system B with coordinates (u,v).  suppose further that
% u = x - d, and v = y - d.  consider a decision line in coordinate
% system B:
% 
%     ( cos(theta), sin(theta) )*( u, v ) = c
% 
% the equation of this decision line in coordinate system A is:
% 
%     ( cos(theta), sin(theta) )*( x-d, y-d ) = c
% 
% which is the same as
% 
%     ( cos(theta), sin(theta) )*( x, y ) = c + d*( cos(theta) + sin(theta) )
% 
% thus we switch from coordinate system B to coordinate system A by
% replacing c with c+d*(cos(theta)+sin(theta)).  this is how we transform
% the parameters of each decision line in lines 96-97.

end
