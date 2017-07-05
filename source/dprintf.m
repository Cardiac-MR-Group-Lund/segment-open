function varargout = dprintf(varargin)
% DPRINTF Write formatted data to string. Same as SPRINTF, but uses
% dictionary to translate string into user selected language.
if nargin > 0
  varargin{1} = translation.dictionary(varargin{1});
end
[varargout{1:nargout}] = sprintf(varargin{:});