function [ dvar1, dvar2 ] = calcdvar_experiment1( signal, sigcst, rngseeds, template )

% CALCDVAR_EXPERIMENT1  Calculate proxy decision variables in experiment 1
% 
%     usage:  [ dvar1, dvar2 ] = calcdvar_experiment1( signal, sigcst, rngseeds, template )
% 
%     input arguments
%         'signal' is an n x 1 matrix of 1's and 2's that encode the signal order; 1 = white disk first, 2 = black disk first
%         'sigcst' is an n x 1 matrix of signal contrasts
%         'rngseeds' is an n x 2 matrix of random number generator seeds
%         'template' is the observer's classification image
% 
%     return arguments
%         'dvar1' is an n x 1 matrix of proxy decision variables for stimulus interval 1
%         'dvar2' is an n x 1 matrix of proxy decision variables for stimulus interval 2

% get signals
signalW = load('signal_experiment1.txt','-ascii');  % white bump
signalB = -signalW;                                 % black bump

% set stimulus parameters
stimsizeP = 31;    % stimulus size, in pixels (the stimulus is square)
noisestd = 0.25;   % noise standard deviation, in contrast units

% initialize proxy decision variables for the two stimulus intervals
dvar1 = NaN([ size(signal,1) 1 ]);
dvar2 = dvar1;

% step through trials
for t=1:size(signal,1)
    
    % seed random number generators
    rand('state',rngseeds(t,1));
    randn('state',rngseeds(t,2));
    
    % get noise
    noiseW = cnormrnd1(0,noisestd,2,[ stimsizeP stimsizeP ]);  % noise on white bump
    noiseB = cnormrnd1(0,noisestd,2,[ stimsizeP stimsizeP ]);  % noise on black bump
    
    % reconstruct stimulus
    if signal(t)==1
        stim1 = sigcst(t)*signalW + noiseW;
        stim2 = sigcst(t)*signalB + noiseB;
    else
        stim1 = sigcst(t)*signalB + noiseB;
        stim2 = sigcst(t)*signalW + noiseW;
    end
    
    % get proxy decision variables for stimulus intervals 1 and 2
    dvar1(t) = sum(sum( template.*stim1 ));
    dvar2(t) = sum(sum( template.*stim2 ));
    
end

end


function r = cnormrnd1( mu, sigma, nclip, dim )

% CNORMRND1  Clipped normal random number generator (version 1)
% 
%     usage:  r = cnormrnd1( mu, sigma, nclip, dim )
% 
%     input arguments
%         'mu' is the mean of the distribution to be sampled
%         'sigma' is the nominal standard deviation of the distribution to be sampled, i.e., we sample from N(mu,sigma^2), but we clip the tails of this distribution
%         'nclip' is the number of standard deviations from the mean at which the distribution is clipped
%         'dim' is the dimension of the return argument
% 
%     output argument
%         'r' is a matrix of random numbers sampled from N(m,sigma^2), except that it contains no values outside the interval ( mu-nclip*sigma, mu+nclip*sigma )

% generate random numbers
r = cnorminv( rand(dim), mu, sigma, nclip );

end


function x = cnorminv( p, mu, sigma, nclip )

% CNORMINV  Clipped normal inverse cumulative distribution function
% 
% x = cnorminv( p, mu, sigma, nclip )

% calculate inverse cdf
ptail=0.5*erfc(nclip/sqrt(2));
pstar=ptail+p.*(1-2*ptail);
z=-sqrt(2)*erfcinv(2*pstar);
x=mu+sigma.*z;

end
