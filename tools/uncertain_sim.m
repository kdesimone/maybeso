function [ dvar1, dvar2, response, pcorrect ] = uncertain_sim( ulevel, sigcst, sigmae, sigmai, ntrials )

% UNCERTAIN_SIM  Simulate an intrinsically uncertain model observer in a 
%                2AFC detection task
% 
%     usage:  [ dvar1, dvar2, response ] = uncertain_sim( ulevel, sigcst, sigmae, sigmai, ntrials )
% 
%     input arguments
%         'ulevel' is the number of irrelevant mechanisms that the observer monitors
%         'sigcst' is the signal contrast
%         'sigmae' is the standard deviation of the external noise
%         'sigmai' is the standard deviation of the internal noise
%         'ntrials' is the number of trials to simulate
% 
%     return arguments
%         'dvar1' is an ntrials x 1 matrix of proxy decision variables for stimulus interval 1
%         'dvar2' is an ntrials x 1 matrix of proxy decision variables for stimulus interval 2
%         'response' is an ntrials x 1 matrix of 1's and 2's that encode the observer's responses:  1 = chose first interval, 2 = chose second interval
%         'pcorrect' is the proportion of correct responses

% here we simulate the observer at the level of the decision variable

% initialize stimulus matrix
% - each row represents a trial
% - each column represents the response of a relevant or irrelevant mechanism
stim1 = zeros([ ntrials ulevel+1 ]);  % stimulus interval 1
stim2 = zeros([ ntrials ulevel+1 ]);  % stimulus interval 2

% set response levels of relevant mechanism (first column)
sigi = ((1:ntrials)'>ntrials/2) + 1;  % mark first and second half of trials
stim1(sigi==1,1) = sigcst;            % signal present in stimulus interval 1 in first half of trials; set first column (i.e., response of relevant mechanism) to sigcst
stim2(sigi==2,1) = sigcst;            % signal present in stimulus interval 2 in second half of trials; set first column (i.e., response of relevant mechanism) to sigcst

% add external noise to relevant mechanism (first column)
stim1(:,1) = stim1(:,1) + sigmae*randn([ ntrials 1 ]);
stim2(:,1) = stim2(:,1) + sigmae*randn([ ntrials 1 ]);

% record proxy decision variables; here we model these as the responses of
% the relevant mechanisms
dvar1 = stim1(:,1);
dvar2 = stim2(:,1);

% add internal noise to responses of all mechanisms
stim1 = stim1 + sigmai*randn([ ntrials ulevel+1 ]);
stim2 = stim2 + sigmai*randn([ ntrials ulevel+1 ]);

% get decision variables:  maximum of all mechanisms
d1 = max(stim1,[],2);
d2 = max(stim2,[],2);

% make responses using difference rule:  choose stimulus interval with largest decision variable
response = (d2>d1) + 1;

% find proportion correct
pcorrect = mean(response==sigi);

end
