% figure4.m  Generate panels of Figure 4 in Pritchett and Murray (2015)
%            i.e., experiment 1, response probability as a function of
%            distance from decision line; also bootstrap an estimate of
%            the width of the guessing region

clear; clc; clf;

% initialize
addpath(fullfile(pwd,'tools'));  % add tool folder to path
assertstats;                     % check for the statistics toolbox

% choose subject (subject = 1, 2, or 3)
subject = 1;
fname = sprintf('data/experiment1_subject%d.txt',subject);
fprintf('processing %s\n',fname);


% part one:  calculate proxy decision space and fit difference model

% load data and select columns
trials = load(fname,'-ascii');
signal = trials(:,3);           % signal order:  1 = white disk first, 2 = black disk first
sigcst = trials(:,6);           % signal contrast
rngseeds = trials(:,[ 4 5 ]);   % random number generator seeds
response = trials(:,7);         % observer response:  1 = judged white disk first, 2 = judged black disk first

% calculate a radially pooled, unit-energy classification image
fprintf('  calculating classification image ...\n');
cimage = calcclass_experiment1( signal, response, rngseeds, subject );

% calculate proxy decision variables
fprintf('  calculating proxy decision variables ...\n');
[ dvar1, dvar2 ] = calcdvar_experiment1( signal, sigcst, rngseeds, cimage );

% calculate proxy decision space
fprintf('  calculating proxy decision space ...\n');
dspace = calcdspace( dvar1, dvar2, response );

% fit the difference model to the decision space
fprintf('  fitting difference model to proxy decision space ...\n');
param = fitgddloop( dspace, 'diff' );


% part two:  fit a normal cdf to response probbility as a function of
% proxy decision variable pair's distance from decision line

% find distance of each ( dvar1(i), dvar2(i) ) pair from decision line
dvec = [ cos(param(1)) sin(param(1)) ]';  % decision vector
crit = param(3);                          % distance of decision line from origin
ddist = [ dvar1 dvar2 ]*dvec - crit;      % distance from decision line

% define normal cdf that saturates at 0.01 and 0.99 for robustness,
% so that occasional random responses from the observer (e.g., keypress
% errors, lapses of attention) don't strongly bias the fit.
fitfn = @( x, p ) max(min( normcdf( x, p(1), p(2) ), 0.99), 0.01);

% make a maximum likelihood objective function
errfn = @( p ) -sum(log( binopdf( response==2, 1, fitfn( ddist, p ) ) ));

% fit the normal cdf; take best of several fits with random starting points
fprintf('  fitting difference model to distance from decision line ...\n');
errmin = Inf;
for i = 1:20
    [ param_tmp, err ] = fminsearch( errfn, [ -0.25+0.5*rand 0.1+0.3*rand ] );
    if err<errmin
        param_norm = param_tmp;
        errmin = err;
    end
end

% make a list of bin centres that spans the range of distances
nbins = 50;                          % number of bins
qmin = quantile(ddist,0.01);         % 0.01 quantile
qmax = quantile(ddist,0.99);         % 0.99 quantile
dlist = linspace(qmin,qmax,nbins);   % list of bin centres
dd = dlist(2)-dlist(1);              % distance between bin centres

% find response probability at each distance from decision line
kresp = zeros(size(dlist));          % number of trials in each bin where observer gave response 2
nresp = zeros(size(dlist));          % number of trials in each bin
presp = zeros(size(dlist));          % proportion of trials in each bin where observer gave response 2
for i = 1:numel(dlist)
    j = abs(ddist-dlist(i))<dd/2;    % find trials in this bin
    kresp(i) = sum(response(j)==2);  % count number of trials with response 2
    nresp(i) = sum(j);               % count number of trials
    presp(i) = kresp(i)/nresp(i);    % find proportion of trials where observer gave response 2
end

% plot the response probabilities and the fitted normal cdf
fplot( @(x) fitfn(x,param_norm), [ -1.3 1.3 ], 'b-' ); hold on;
h = get(gca,'Children');
set(h,'LineWidth',2);
h = plot(dlist(nresp>=20),presp(nresp>=20),'ro'); hold off;
set(h,'MarkerSize',10,'MarkerFaceColor','r');
set(gca,'FontName','helvetica','FontWeight','bold','FontSize',24);
xlabel 'distance from decision line'
ylabel 'response probability'
axis([ -1.3 1.3 0 1 ]);
drawnow;


% part three:  bootstrap the difference model with guessing

% define transition function for difference rule with guessing
fitfn = @( x, p ) max(min( 0.5*normcdf( x, p(1), p(3) ) + 0.5*normcdf( x, p(2), p(3) ), 0.99),0.01);
% explanation:  according to the difference model with guessing, the
% probability of the observer giving response 2 is a stepwise function of
% the distance of the decision variable pair ( dvar1, dvar2 ) from
% the decision line.  if the pair is sufficiently far to one side of the
% decision line, the observer gives response 1.  if the pair is
% sufficiently far to the other side, the observer gives response 2.  in an
% intermediate region, the observer guesses.  here we assume that in the
% intermediate region the observer gives response 2 with probability 0.5.
% thus the response probability is
% 
%     P( response=2 ) = 0.5*step( dvar - criterion1 ) + 0.5*step( dvar - criterion2 )
% 
% here step( x ) = 0 if x<0 and 1 if x>=0.  dvar is the distance of the
% decision variable pair from the decision line.  criterion1 is the
% distance at which the response probability transitions from 0 to 0.5, and
% criterion2 is the distance at which the response probability transitions
% from 0.5 to 1.
% 
% as we show in the supporting information, the proxy decision space is the
% true decision space convolved with a gaussian kernel.  this means that
% if pdvar is the distance of the proxy decision variable pair from the
% decision line, then the response probability is
% 
%     P( response=2 ) = 0.5*normcdf( pdvar, criterion1, sigma ) + 0.5*normcdf( pdvar, criterion2, sigma )
% 
% here normcdf( x, mu, sigma ) is the normal cumulative distribution
% function.  sigma is the standard deviation of the gaussian kernel that
% transforms the decision space into the proxy decision space.  this is
% the function that we fit to the distance of the proxy decision variable
% pair from the decision line in the line above, except that we clip
% the function to the range (0.01,0.99) to provide some robustness
% against occasional lapses of attention and keypress errors.

% initialize
B = 1000;            % number of bootstrap iterations
param = zeros(B,3);  % bootstrapped parameters

% bootstrap
fprintf('  bootstrapping width of guessing region ...\n');
for b = 1:B+1
    
    if b == 1
        % on first pass, use the actual data to get the maximum likelihood estimate
        dstar = ddist;        % distance from decision line
        rstar = response;     % observer's response
    else
        % on subsequent passes, resample the trials with replacement to get
        % bootstrapped estimates
        i = unidrnd( numel(ddist), [ numel(ddist) 1 ] );
        dstar = ddist(i);     % distance from decision line
        rstar = response(i);  % observer's response
    end
    
    % make a maximum likelihood objective function
    errfn = @( p ) -sum(log( binopdf( rstar==2, 1, fitfn( dstar, p ) ) ));
    
    % fit the difference model with guessing; take best of several fits
    % with random starting points
    errmin = Inf;
    for i = 1:20
        [ param_tmp, err ] = fminsearch(errfn,[ -0.25*0.5*rand(1,2) 0.1+0.3*rand ]);
        % if fit error is lower than best fit so far, then keep this fit
        if err<errmin
            param(b,:) = param_tmp;
            errmin = err;
        end
    end
    
    % show progress
    fprintf('    %5d of %5d iterations\n',b,B);
    
end

% report the maximum likelihood estimate and 95% confidence interval for the width of the guessing region
guesswidth = abs(param(:,1)-param(:,2));
fprintf('guessing region width:  maximum likelihood estimate %.2f, 95%% confidence interval ( %.2f, %.2f )\n',guesswidth(1),quantile(guesswidth(2:end),0.025),quantile(guesswidth(2:end),0.975));
