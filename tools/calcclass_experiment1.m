function cimage = calcclass_experiment1( signal, response, rngseeds, subject )

% CALCCLASS_EXPERIMENT1  Calculate classification image for experiment 1
% 
%     usage:  cimage = calcclass_experiment1( signal, response, rngseeds, subject )
% 
%     input arguments
%         'signal' is an n x 1 matrix of 1's and 2's that encode the signal order; 1 = white disk first, 2 = black disk first
%         'response' is an n x 1 matrix of 1's and 2's that encode the observer's responses:  1 = judged white disk first, 2 = judged black disk first
%         'rngseeds' is an n x 2 matrix of random number generator seeds
%         'clipradiusP' is the radius, in pixels, of the circle outside which the classification image is set to zero
% 
%     return argument
%         'cimage' is a radially pooled, unit-energy classification image

% set stimulus parameters
stimsizeP = 31;    % stimulus size, in pixels (the stimulus is square)
stimmidP = 16;     % location of centre of signal is (stimmidP, stimmidP)
noisestd = 0.25;   % noise standard deviation, in contrast units

% we set the classification image to zero a certain number of pixels from
% the centre of the signal, in order to to reduce the noise; in
% experiment 1 we find that one subject had a significantly broader
% classification image than the others; here we choose the radius of the
% classification image for the present subject.
clipradiuslist = [ 6 10 6 ];
clipradiusP = clipradiuslist( subject );

% initialize classification image bins
cbin = repmat( { zeros(stimsizeP) }, [ 2 2 ] );  % sum of stimulus noise on trials where stimulus order is i and observer's response is j
ntrials = zeros(2,2);                            % count of number of trials in this bin

% step through trials
for t=1:size(signal,1)
    
    % seed random number generators
    rand('state',rngseeds(t,1));
    randn('state',rngseeds(t,2));
    
    % get noise
    noiseW = cnormrnd1(0,noisestd,2,[ stimsizeP stimsizeP ]);  % noise on white bump
    noiseB = cnormrnd1(0,noisestd,2,[ stimsizeP stimsizeP ]);  % noise on black bump
    
    % combine noise fields using standard method for 2AFC tasks:  noise in
    % stimulus interval 1 minus noise in stimulus interval 2
    if signal(t)==1
        noise = noiseW - noiseB;
    else
        noise = noiseB - noiseW;
    end
    
    % add noise field to appropriate bin
    cbin{signal(t),response(t)} = cbin{signal(t),response(t)} + noise;
    ntrials(signal(t),response(t)) = ntrials(signal(t),response(t)) + 1;
    
end

% calculate classification image
cim = ( cbin{1,1}/ntrials(1,1) + cbin{2,1}/ntrials(2,1) ) - ( cbin{1,2}/ntrials(1,2) + cbin{2,2}/ntrials(2,2) );

% get map of distance in pixels from centre of signal
r = distmap( [ stimsizeP stimsizeP ], [ stimmidP stimmidP ] );

% make radially pooled classification image
rvec = 0:16;              % distances to consider
mvec = NaN(size(rvec));   % mean of classification image at each distance from centre
for i = 1:numel(rvec)
    mvec(i) = mean(cim( r>=rvec(i)-0.5 & r<rvec(i)+0.5 ));
end

% make 2D version of radially pooled image
cimage = zeros(stimsizeP);
cimage(r<=clipradiusP) = linterp( r(r<=clipradiusP), rvec, mvec );

% scale classification image to unit energy
cimage = cimage/sqrt(sum(cimage(:).^2));

end


function r = distmap( dim, originij )

% DISTMAP  Make a matrix whose elements are the distance from a particular
%          element in the matrix
% 
%     usage:  r = distmap( dim, originij )
% 
%     input arguments
%         'dim' gives the size of the matrix being considered, i.e., dim(1) rows and dim(2) columns
%         'originij' gives the (i,j) subscripts of the centre of the matrix
% 
%     output argument
%         'r' is a matrix that gives the distance of each element in the matrix from the element at position 'originij'

% make coordinate matrices
ivec = (1:dim(1)) - originij(1);
jvec = (1:dim(2)) - originij(2);
[j,i] = meshgrid(jvec,ivec);

% calculate distance map
r = sqrt( i.^2 + j.^2 );

end


function y = linterp( x, xref, yref )

% LINTERP  Piecewise linear interpolation function
%
%     usage:  y = linterp( x, xref, yref )
% 
%     input arguments
%         'xref' and 'yref' are n x 1 matrices of reference points
%         'x' is a matrix of x-values for which we want to find piecewise linear interpolations
% 
%     output argument
%         'y' is a matrix of values interpolated from 'xref' and 'yref', in a piecewise linear fashion, at the positions in 'x'

% initialize
y = NaN(size(x));

% step through x values
for i = 1:numel(x)
    
    % find next lower x value, if any
    f = find(xref<=x(i));
    if ~isempty(f)
        [dx1,f1] = min( x(i) - xref(f) );
        i1 = f(f1);
    else
        continue
    end
    
    % find next higher x value, if any
    f = find(xref>=x(i));
    if ~isempty(f)
        [dx2,f2] = min( xref(f) - x(i) );
        i2 = f(f2);
    else
        continue
    end
    
    % interpolate between nearest pair of reference points
    if (dx1+dx2)==0
        y(i) = yref(i1);
    else
        m = (yref(i2)-yref(i1))/(dx1+dx2);
        y(i) = yref(i1) + m*dx1;
    end
    
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
%         'dim' is the size of the return argument
% 
%     output argument
%         'r' is a matrix of random numbers sampled from N(m,sigma^2), except that it contains no values outside the interval ( mu-nclip*sigma, mu+nclip*sigma )

% generate random numbers using the inverse transform method
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
