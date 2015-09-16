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