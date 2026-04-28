function varargout = rprintf(varargin)
%Special function to translate phrases in report, works as DPRINTF.
%
%See also DPRINTF

%Jelena Bock

global DATA 
language = DATA.Pref.ReportLanguage;
fromlanguage = 'English';
if nargin > 0
  varargin{1} = translation.dictionary(varargin{1},language,fromlanguage);
end
[varargout{1:nargout}] = sprintf(varargin{:});