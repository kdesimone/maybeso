% figureS3.m  Find a 95% confidence interval for the width of the guessing
%             region when fitting the difference model with guessing to
%             data from the difference model

clear; clc; clf;

% initialize
addpath(fullfile(pwd,'tools'));  % add tool folder to path
assertstats;                     % check for the statistics toolbox

% to make the simulation as realistic as possible we use the psychometric
% function fitted to the human observer's data in the leftmost panel of Figure 4.
% fitted values:  mu = 0.00, sigma = 0.26
psymet = @( x ) max(min( normcdf( x, 0.00, 0.26 ), 0.99), 0.01);

% to make the simulation as realistic as possible we also use the
% distances from the decision line, and the number of trials at each
% distance from the decision line, shown in the leftmost panel of Figure 4.
ddist = [ -1.33 -1.28 -1.23 -1.18 -1.13 -1.08 -1.03 -0.99 -0.94 ...
          -0.89 -0.84 -0.79 -0.74 -0.69 -0.64 -0.59 -0.54 -0.50 ...
          -0.45 -0.40 -0.35 -0.30 -0.25 -0.20 -0.15 -0.10 -0.05 ...
          -0.01 0.04 0.09 0.14 0.19 0.24 0.29 0.34 0.39 0.44 0.48 ...
          0.53 0.58 0.63 0.68 0.73 0.78 0.83 0.88 0.93 ...
          0.97 1.02 1.07 ];  % stimulus levels, i.e., distances from the fitted decision line
ntrials = [ 5 5 7 5 7 9 10 12 13 22 35 41 46 78 105 160 202 262 ...
            257 327 401 439 445 483 508 538 514 514 535 562 451 ...
            464 455 370 350 237 211 165 118 92 62 69 35 28 15 ...
            14 8 5 4 9 ];  % number of trials at each stimulus level

% find the observer's proportion correct at each stimulus level, according
% to the fitted psychometric function
ptheory = psymet( ddist );

% set number of repetitions
nrepeat = 1000;

% initialize list of guessing region widths
guesswidth = NaN([ nrepeat 1 ]);

% repeat many simulated experiments
fprintf('estimating width of guessing region ...\n');
for t = 1:nrepeat
    
    % get simulated data from the fitted psychometric function, which is
    % an instance of the difference model
    kresp = binornd( ntrials, ptheory );  % number of trials where observer makes response 2
    presp = kresp./ntrials;               % proportion of trials where observer makes response 2

    % define the form of the psychometric function according to the
    % difference model with guessing; see figure4.m, lines 105-136, for an
    % explanation of why the function takes this form
    fitfn_diffg = @( x, p ) max(min( 0.5*normcdf( x, p(1), p(3) ) + 0.5*normcdf( x, p(2), p(3) ), 0.99),0.01);

    % define a maximum likelihood objective function
    errfn_diffg = @( p ) -sum(log( binopdf( kresp, ntrials, fitfn_diffg( ddist, p ) ) ));

    % fit difference model with guessing to simulated data; take best fit
    % from 20 fits with random starting points
    errmin_diffg = Inf;
    for i = 1:20
        [ param_tmp, err ] = fminsearch( errfn_diffg, [ -0.25*0.5*rand(1,2) 0.1+0.3*rand ]);
        if err<errmin_diffg
            param_diffg = param_tmp;
            errmin_diffg = err;
        end
    end
    
    % find guessing interval width
    guesswidth(t) = abs( param_diffg(1)-param_diffg(2) );
    fprintf('  iteration %d of %d, guessing region width = %.4f\n',t,nrepeat,guesswidth(t));
    
    % plot data and fit of difference model with guessing
    plot( ddist, presp, 'ro' ); hold on;
    xlim = get(gca,'XLim');
    fplot( @( x ) fitfn_diffg( x, param_diffg ), xlim, 'r-' );
    hold off; drawnow;
    
end

% report range of guessing region fits
fprintf('guessing region width:  95%% confidence interval ( %.2f, %.2f )\n',quantile(guesswidth,0.025),quantile(guesswidth,0.975));
