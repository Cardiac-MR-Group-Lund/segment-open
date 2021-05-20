%------------------------------
function varargout = rprintf(varargin)
%------------------------------
% special function to translate phrases in report
global DATA 
language = DATA.Pref.ReportLanguage;
fromlanguage = 'English';
if nargin > 0
  varargin{1} = translation.dictionary(varargin{1},language,fromlanguage);
end
[varargout{1:nargout}] = sprintf(varargin{:});