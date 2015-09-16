function [ dvar1, dvar2 ] = calcdvar_experiment2( signal, sigcst, rngseeds, template )

% CALCDVAR_EXPERIMENT2  Calculate decision variables in experiment 2
% 
%     usage:  [ dvar1, dvar2 ] = calcdvar_experiment2( signal, sigcst, rngseeds, template )
% 
%     input arguments
%         'signal' is an n x 1 matrix of 1's and 2's that encode the signal position; 1 = left, 2 = right
%         'sigcst' is an n x 1 matrix of signal contrasts
%         'rngseeds' is an n x 2 matrix of random number generator seeds
%         'template' is the observer's classification image
% 
%     return arguments
%         'dvar1' is an n x 1 matrix of proxy decision variables for stimulus interval 1
%         'dvar2' is an n x 1 matrix of proxy decision variables for stimulus interval 2

% get signal
leftdot = load('signal_experiment2.txt','-ascii');

% set stimulus parameters
stimheightP = 38;  % stimulus height, in pixels
stimwidthP = 76;   % stimulus width, in pixels
basecst = 0.10;    % signal pedestal, in contrast units
noisestd = 0.20;   % noise standard deviation, in contrast units

% initialize proxy decision variables for the two stimulus intervals
dvar1 = NaN([ size(signal,1) 1 ]);
dvar2 = dvar1;

% step through trials
for t=1:size(signal,1)
    
    % seed random number generators
    rand('state',rngseeds(t,1));
    randn('state',rngseeds(t,2));

	% reconstruct signal
    sig = (basecst+sigcst(t))*leftdot + basecst*fliplr(leftdot);  % signal with dot on left
    if signal(t)==2
        sig = fliplr(sig);  % if dot should be on right, flip it
    end
    
    % reconstruct noise
	noise = cnormrnd2(0,noisestd,stimheightP,stimwidthP);
    
    % reconstruct stimulus
    stim = sig + noise;
    
    % get left and right proxy decision variables
    dvar1(t) = sum(sum( stim.*template ));
    dvar2(t) = sum(sum( stim.*fliplr(template) ));

end

end


function r = cnormrnd2( mu, sigma, m, n )

% CNORMRND2  Clipped normal random number generator (version 2)
% 
%     usage:  r = cnormrnd2( mu, sigma, m, n )
% 
%     input arguments
%         'mu' is the mean of the distribution to be sampled
%         'sigma' is the nominal standard deviation of the distribution to be sampled, i.e., we sample from N(mu,sigma^2), but we clip the tails of this distribution
%         [ m n ] is the size of the return argument
% 
%     output argument
%         'r' is a matrix of random numbers sampled from N(m,sigma^2), except that it contains no values outside the interval ( mu-nclip*sigma, mu+nclip*sigma )

r = mu + sigma*randn(m,n);
out = ( (r<=mu-2*sigma) | (r>=mu+2*sigma) );
outf = find(out);
outn = size(outf,1);
okf = find(~out);
r(outf) = r(okf(1:outn));

end
