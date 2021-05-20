function catvect=mycat(varargin)
%Same as cat, but removes nans and converts them to handles. Fix bug with
%R2015
%Klas Berggren, comments by Einar Heiberg

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end

for loop = 1:nargin
  varargin{loop} = double(varargin{loop});
end

catvect = cat(varargin{:});
catvect = catvect(~isnan(catvect));