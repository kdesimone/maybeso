function cimage = calcclass_experiment2( signal, response, rngseeds )

% CALCCLASS_EXPERIMENT2  Calculate classification image for experiment 2
% 
%     usage:  cimage = calcclass_experiment2( signal, response, rngseeds )
% 
%     input arguments
%         'signal' is an n x 1 matrix of 1's and 2's that encode the signal order; 1 = left, 2 = right
%         'response' is an n x 1 matrix of 1's and 2's that encode the observer's responses:  1 = judged increment on left, 2 = judged increment on right
%         'rngseeds' is an n x 2 matrix of random number generator seeds
% 
%     return argument
%         'cimage' is a radially pooled, unit-energy classification image

% set stimulus parameters
stimheightP = 38;  % stimulus height, in pixels
stimwidthP = 76;   % stimulus width, in pixels
noisestd = 0.20;   % noise standard deviation, in contrast units
radiusP = 6;       % radius, in pixels, of the circle outside which we set the classification image to zero

% initialize classification image bins
cbin = repmat( { zeros(stimheightP,stimwidthP) }, [ 2 2 ] );  % sum of stimulus noise on trials where stimulus order is i and observer's response is j
ntrials = zeros(2,2);                                         % count of number of trials in this bin

% step through trials
for t=1:size(signal,1)
    
    % seed random number generators
    rand('state',rngseeds(t,1));
    randn('state',rngseeds(t,2));
    
    % get noise
    noise = cnormrnd2(0,noisestd,stimheightP,stimwidthP);
    
    % add noise field to appropriate bin
    cbin{signal(t),response(t)} = cbin{signal(t),response(t)} + noise;
    ntrials(signal(t),response(t)) = ntrials(signal(t),response(t)) + 1;
    
end

% calculate classification image
cim = ( cbin{1,1}/ntrials(1,1) + cbin{2,1}/ntrials(2,1) ) - ( cbin{1,2}/ntrials(1,2) + cbin{2,2}/ntrials(2,2) );

% get maps of distance in pixels from centres of signals
rleft =  distmap( [ stimheightP stimwidthP ], [ 20 20 ] );  % left dot is centred at matrix element (20,20)
rright = distmap( [ stimheightP stimwidthP ], [ 20 57 ] );  % right dot is centred at matrix element (20,57)

% make radially pooled classification image, combining left and right dots
rvec = 0:10;              % distances to consider
mvec = NaN(size(rvec));   % mean of classification image at each distance from centre
for i = 1:numel(rvec)
    mvec(i) = mean(cim( rleft>=rvec(i)-0.5  & rleft<rvec(i)+0.5  )) - ...  
              mean(cim( rright>=rvec(i)-0.5 & rright<rvec(i)+0.5 ));
end

% make 2D version of radially pooled image
cimage = zeros([ stimheightP stimwidthP ]);
cimage(rleft<=radiusP) = linterp( rleft(rleft<=radiusP), rvec, mvec );

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
