% figure2c.m  Generate panels of Figure 2c in Pritchett and Murray (2015)
%             i.e., proxy decision space and GDD function fit for a
%             simulated observer with intrinsic uncertainty

clear; clc; clf;

% initialize
addpath(fullfile(pwd,'tools'));  % add tool folder to path
assertstats;                     % check for the statistics toolbox

% set observer parameters
ulevel = 0;   % number of irrelevant mechanisms being monitored (allowed values are 0 to 8)
sigmae = 1;   % standard deviation of external noise
sigmai = 1;   % standard deviation of internal noise

% get signal contrast that gives 71% correct performance (found empirically)
sigcstlist = [ 1.10 1.35 1.48 1.57 1.64 1.69 1.74 1.78 1.82 ];
sigcst = sigcstlist(ulevel+1);

% set number of simulated trials
ntrials = 10000;

% get proxy decision variable and responses from simulated observer with intrinsic uncertainty
fprintf('simulating model observer ...\n');
[ dvar1, dvar2, response ] = uncertain_sim( ulevel, sigcst, sigmae, sigmai, ntrials );

% calculate proxy decision space
fprintf('calculating proxy decision space ...\n');
dspace = calcdspace( dvar1, dvar2, response );

% fit GDD function to proxy decision space
fprintf('fitting GDD function to proxy decision space ...\n');
param = fitgddloop( dspace, 'gdd' );

% show proxy decision space and fitted GDD function
plotdspace( dspace, param );
