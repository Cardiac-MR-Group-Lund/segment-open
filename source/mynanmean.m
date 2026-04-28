function m = mynanmean(x,dim)
% Mean value, ignoring NaNs.

% new implemenation under #2726

% original implemenation can be found under internaltescripts/test_mynanmean_versions

if nargin == 1 || isempty(dim)
  m = mean(x,"omitnan");
else
  m = mean(x,dim,"omitnan");
end