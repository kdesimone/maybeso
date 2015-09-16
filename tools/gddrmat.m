function pmat = gddrmat( dvar, theta1, theta2, delta1, delta2, gamma, sigma )

% GDDRMAT  Generalized double detection (GDD) function response probability matrix
%
%     usage:  pmat = gddrmat( dvar, theta1, theta2, delta1, delta2, gamma, sigma )
%
%     input arguments
%         'dvar' is a n x 1 matrix of evenly spaced decision variable values at which we'll calculate the GDD function
%         'theta1', 'theta2', 'delta1', 'delta2', 'gamma' are the parameters of the GDD function, as defined in Pritchett and Murray's (2015) equation (2)
%         'sigma' is the blur parameter explained in the paragraph following Pritchett and Murray's (2015) equation (7)
%
%     output argument
%         'pmat' is an n x n matrix, where element (i,j) is the probability that the observer gives response 2 when decision variable 1 has value dvar(i) and decision variable 2 has value dvar(j)
% 
%     example:
%         dvar = -1:0.05:1;
%         pmat = gddrmat( dvar, 3*pi/4, 3*pi/4, 0.1, -0.1, 0.5, 0.1 );
%         imagesc( dvar, dvar, pmat );
%         axis xy
%         xlabel 'decision variable 1'
%         ylabel 'decision variable 2'
%         title 'response probability'
%         colorbar

% if sigma is not positive, print a warning message and return a 
% probability matrix that has value 0.5 everywhere
if sigma<=0
%     warning('invalid paramter value:  sigma<=0');
    pmat = repmat( 0.5, [ numel(dvar) numel(dvar) ] );
    return
end

% calculate the size of convolution kernel we'll use on the decision space
% to blur by a 2D gaussian of width sigma
ddvar = dvar(2)-dvar(1);      % distance between decision variable bins
kn = ceil(6*sigma/ddvar)+2;   % six standard deviations, plus a bit more

% extend the list of decision variable values
dvare = [ dvar(1)+ddvar*(-kn:-1) dvar(:)' dvar(end)+ddvar*(1:kn) ]';

% get coordinate matrices for decision variables 1 and 2
[ d1mat, d2mat ] = meshgrid(dvare,dvare);

% assign response probabilities 0 and 1 based on the two decision lines
yes1 = double( ( d1mat*cos(theta1) + d2mat*sin(theta1) ) >= delta1 );  % decision line 1, as defined in Pritchett and Murray's (2015) equation (2)
yes2 = double( ( d1mat*cos(theta2) + d2mat*sin(theta2) ) >= delta2 );  % decision line 2, as defined in Pritchett and Murray's (2015) equation (2)

% in the above two lines we assigned response probabilities 0 and 1 to 
% all elements of the probability matrix.  in fact, matrix elements that
% are split by the decision line have response probabilities in between
% 0 and 1.  from here on we calculate the response probabilities in 
% elements that are split by the decision line

% calculate maximum distance that a square of length ddvar can be from
% decision line 1, and still be split by decision line 1
theta1clip = abs( mod( theta1+pi/4, pi/2 ) - pi/4 );
maxdist1 = ddvar*0.5*abs( sin(theta1clip) + cos(theta1clip) );

% find the distance of the centre of each probability matrix element from
% decision line 1
dist1 = abs( cos(theta1)*d1mat + sin(theta1)*d2mat - delta1 );

% for probability matrix elements that are within the critical distance of
% decision line 1, calculate the response probability
close1 = find(dist1<maxdist1);
for i = close1'
    yes1(i) = pcover( [ d1mat(i) d2mat(i) ], ddvar, theta1, delta1 );
end

% calculate maximum distance that a square of length ddvar can be from
% decision line 2, and still be split by decision line 2
theta2clip = abs( mod( theta2+pi/4, pi/2 ) - pi/4 );
maxdist2 = ddvar*0.5*abs( sin(theta2clip) + cos(theta2clip) );

% find the distance of the centre of each probability matrix element from
% decision line 2
dist2 = abs( cos(theta2)*d1mat + sin(theta2)*d2mat - delta2 );

% for probability matrix elements that are within the critical distance of
% decision line 2, calculate the response probability
close2 = find(dist2<maxdist2);
for i = close2'
    yes2(i) = pcover( [ d1mat(i) d2mat(i) ], ddvar, theta2, delta2 );
end

% at this point we have two matrices, yes1 and yes2, that describe the
% response probabilities according to decision lines 1 and 2, respectively.
% now we construct the probability matrix (pmat) that gives the final
% response probability according to the GDD function.

% *** ok to here

% initialize probability matrix to NaN
pmat = NaN(size(d1mat));

% deal with cases where both decision lines have response probability 0 or 1
pmat( yes1==1 & yes2==1 ) = 1;         % where decision lines agree that the response is 2, assign response 2, i.e., P(response=2) = 1
pmat( yes1==0 & yes2==0 ) = 0;         % where decision lines agree that the response is 1, assign response 1, i.e., P(response=2) = 0
pmat( yes1==1 & yes2==0 ) = gamma;     % where decision lines disagree, guess
pmat( yes1==0 & yes2==1 ) = gamma;     % where decision lines disagree, guess

% deal with cases where one decision rule has response probability 0 or 1
% and the other does not
i = yes1>0 & yes1<1 & yes2==1;  pmat(i) = yes1(i) + (1-yes1(i))*gamma;
i = yes1>0 & yes1<1 & yes2==0;  pmat(i) = yes1(i)*gamma;
i = yes2>0 & yes2<1 & yes1==1;  pmat(i) = yes2(i) + (1-yes2(i))*gamma;
i = yes2>0 & yes2<1 & yes1==0;  pmat(i) = yes2(i)*gamma;

% deal with remaining cases, where neither decision rule has response
% probability 0 or 1
f = find(isnan(pmat))';
for i = f
    
    % make coordinate matrices for the very small region of the decision
    % space represented by this element of the probability matrix
    xvec = linspace(d1mat(i)-ddvar/2,d1mat(i)+ddvar/2,100);
    yvec = linspace(d2mat(i)-ddvar/2,d2mat(i)+ddvar/2,100);
    [x,y] = meshgrid(xvec,yvec);

    % assign response probabilities according to decision lines 1 and 2
    yes1 = double( ( x*cos(theta1) + y*sin(theta1) ) >= delta1 );  % decision line 1
    yes2 = double( ( x*cos(theta2) + y*sin(theta2) ) >= delta2 );  % decision line 2

    % construct probability matrix for this very small region of the
    % decision space
    p = repmat(gamma,size(x));          % where decision lines disagree, assign guessing rate gamma
    p( yes1==1 & yes2==1 ) = 1;         % where decision lines agree, assign response 2, i.e., P(response=2) = 1
    p( yes1==0 & yes2==0 ) = 0;         % where decision lines agree, assign response 1, i.e., P(response=2) = 0
    
    % find average response probability within this small region of
    % decision space
    pmat(i) = mean(p(:));
    
end

% at this point we have the probability matrix pmat that gives the
% response probabilities in the decision space.  now we convert this to
% the proxy decision space by convolving it with a Gaussian kernel
% with space constant sigma.

% construct a 2D gaussian kernel with width sigma and unit volume
kvec = (1:kn) - (floor(kn/2)+1);
kvec = kvec * ddvar;
[kx,ky] = meshgrid(kvec);
kr = sqrt(kx.^2+ky.^2);
k = exp( -0.5*(kr.^2)/sigma^2 );
k = k/sum(k(:));

% convolve the probability matrix pmat with the gaussian kernel
pmat = conv2( pmat, k, 'same' );

% clip the probability matrix pmat to restore the original decision variable range
pmat = pmat(kn+1:end-kn,kn+1:end-kn);

end


function p = pcover( centre, len, theta, delta )

% PCOVER  Proportion of a square that falls on one half of a decision line
% 
%     usage:  p = pcover( centre, len, theta, delta )
% 
%     input arguments
%         'centre' is a 1 x 2 vector that gives the (x,y) coordinates of the centre of the square
%         'len' is the length of the sides of the square
%         'theta' and 'delta' are the parameters of the decision line; theta is in radians;
%             the decision line is x*cos(theta) + y*sin(theta) = delta
% 
%     return arguments
%         'p' is the proportion of the square that falls in the region
%             defined by ( cos(theta), sin(theta) ) * ( x, y ) > delta,
%             where * is the vector dot product

% find vector perpendicular to decision line
v = [ cos(theta) sin(theta) ]';

% find decision line parameters for equivalent problem with unit square centred at the origin
delta = delta - v(:)'*centre(:);
delta = delta/len;

% now we're trying to find the proportion of the unit square, centered
% on the origin, that falls on the side of the decision line that is
% in the direction of the vector v calculated above.

% find point of intersection (if any) with bottom side of square
if cos(theta)~=0
    % if decision line is not horizontal, find the x coordinate of the point
    % where decision line intersects the line y=-0.5
    xb = (delta-(-0.5)*sin(theta))/cos(theta);
else
    % if decision line is horizontal, then there is no point of intersection
    % with the bottom side of the square
    xb = Inf;
end
% flag whether the intersection occurs between x=-0.5 and x=0.5 (i.e.,
% within the width of the square)
bflag = abs(xb)<=0.5;

% find point of intersection (if any) with right side of square
if sin(theta)~=0
    % if decision line is not vertical, find the y coordinate of the point
    % where decision line intersects the line x=0.5
    yr = (delta-0.5*cos(theta))/sin(theta);
else
    % if decision line is vertical, then there is no point of intersection
    % with the right side of the square
    yr = Inf;
end
% flag whether the intersection occurs between y=-0.5 and y=0.5 (i.e.,
% within the height of the square)
rflag = abs(yr)<0.5;

% if decision line intersects both the bottom and the right side of the square,
% then calculate the proportion of the square on one side of decision line.
if bflag && rflag
    % calculate covered area
    p = 0.5*(0.5-xb)*(yr+0.5);
    if v(:)'*[ 0.5 -0.5 ]'<0
        p = 1-p;
    end
    return
end

% if decision line doesn't intersect the bottom side, but does intersect
% the right side, then just rotate the problem -90 degrees.
if ~bflag && rflag
    % rotate -90 degrees and recalculate
    p = pcover( [ 0 0 ], 1, theta-pi/2, delta );
    return
end

% find point of intersection (if any) with top side of square
if cos(theta)~=0
    % if decision line is not horizontal, find the x coordinate of the point
    % where decision line intersects the line y=0.5
    xt = (delta-0.5*sin(theta))/cos(theta);
else
    % if decision line is horizontal, then there is no point of intersection
    % with the top side of the square
    xt = Inf;
end
% flag whether the intersection occurs between x=-0.5 and x=0.5 (i.e.,
% within the width of the square)
tflag = abs(xt)<=0.5;

% if decision line intersects both the bottom and the top side of the
% sdquare, then calculate the proportion of the square on one side of the
% decision line.
if bflag && tflag
    % calculate covered area
    a = 0.5-xt;
    b = 0.5-xb;
    p = 0.5*(a+b)*1;
    if v(:)'*[ 1 0 ]'<0
        p = 1-p;
    end
    return
end

% if decision line doesn't intersect the bottom side, but does intersect
% the top side, then just rotate the problem -180 degrees.
if ~bflag && tflag
    % rotate -180 degrees and recalculate
    p = pcover( [ 0 0 ], 1, theta-pi, delta );
    return
end

% find point of intersection (if any) with left side of square
if sin(theta)~=0
    % if decision line is not vertical, find the y coordinate of the point
    % where decision line intersects the line x=-0.5
    yl = (delta-(-0.5)*cos(theta))/sin(theta);
else
    % if decision line is vertical, then there is no point of intersection
    % with the left side of the square
    yl = Inf;
end
% flag whether the intersection occurs between y=-0.5 and y=0.5 (i.e.,
% within the height of the square)
lflag = abs(yl)<0.5;

% if decision line intersects the bottom side and the left side, then just
% flip the problem left to right.
if bflag && lflag
    % flip left-to-right and recalculate
    p = pcover( [ 0 0 ], 1, pi-theta, delta );
    return
end

% whole square is on one side of decision line, e.g., decision line lies
% exactly along one side of the square.
if delta>0
    p = 0;
else
    p = 1;
end

end
