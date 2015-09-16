function dspace = calcdspace( dvar, response )

% CALCDSPACE  Calculate a proxy decision space
%
%     usage:  dspace = calcdspace( dvar1, dvar2, response )
%
%     input arguments
%         'dvar1' is an n x 1 matrix of proxy decision variables for stimulus interval 1
%         'dvar2' is an n x 1 matrix of proxy decision variables for stimulus interval 2
%         'response' is an n x 1 matrix of 1's and 2's that encode the observer's responses:  1 = judged white disk first, 2 = judged black disk first
%
%     return argument
%         'dspace' is a struct with the following fields
%             'dspace.dlist' is an 1 x n list of evenly spaced values of proxy decision variables
%             'dspace.pmat'  is an n x n matrix; element (i,j) is the proportion of trials on which the observer gave response 2 when 'dvar1' was in a small bin centred on dspace.dlist(i), and 'dvar2' was in a small bin centred on dspace.dlist(j)
%             'dspace.kmat'  is an n x n matrix; element (i,j) is the number     of trials on which the observer gave response 2 when 'dvar1' was in a small bin centred on dspace.dlist(i), and 'dvar2' was in a small bin centred on dspace.dlist(j)
%             'dspace.nmat'  is an n x n matrix; element (i,j) is the number     of trials on which 'dvar1' was in a small bin centred on dspace.dlist(i), and 'dvar2' was in a small bin centred on dspace.dlist(j)

% set number of bins into which we divide the proxy decision variables
nbin = 25;

% find range of proxy decision variables
alld = [ dvar(:) ];     % pool proxy decision variables from the two intervals
qmin = quantile( alld, 0.01 );      % find quantile 0.01
qmax = quantile( alld, 0.99 );      % find quantile 0.99
qmean = (qmin+qmax)/2;              % find the middle of the range of proxy decision variables
qmin = qmean + 1.25*(qmin-qmean);   % choose lowest value of proxy decision variable to consider
qmax = qmean + 1.25*(qmax-qmean);   % choose highest value of proxy decision variable to consider

% make bin boundaries and centres
step = (qmax-qmin)/nbin;            % choose step size between bins
dbins = qmin:step:qmax;             % bin boundaries
dlist = dbins(1:end-1)+(step/2);    % bin centres

% initialize matrices
kmat = NaN(nbin,1);                   % number of trials in each bin where observer gave response 2
nmat = NaN(nbin,1);                   % number of trials in each bin

% find response counts and trial counts in each bin
for j = 1:numel(dbins)-1
    k = (dvar>=dbins(j)) & (dvar<dbins(j+1));
    kmat(j) = sum(response(k)==1);  % count trials where observer gave response 2
    nmat(j) = sum(k);               % count trials
end

% proportion of trials in each bin where observer gave response 2
pmat = kmat./nmat;

% assemble return argument
dspace.dlist = dlist;   % bin centres
dspace.pmat = pmat;     % proportion of trials where observer gave response 2
dspace.kmat = kmat;     % number of trials where observer gave response 2
dspace.nmat = nmat;     % number of trials

return
