function assertstats

% ASSERTSTATS  Halt with an error if the statistics toolbox is not on the path
% 
%     usage:  assertstats

assert(~isempty(which('quantile')),'This package (pnas2015) requires the MATLAB Statistics Toolbox');

end
