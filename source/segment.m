function varargout = segment(varargin)
%This is a script that passes calls of segment(fcnname,arg1,arg2...) to the 
%object DATA, if the method fcnname exists. If the method does not exist, 
%the function fcnname is called from segment_main.m, which is the file 
%previously known as segment.m.

global DATA
if nargin < 1
  [varargout{1:nargout}] = segment_main;
else
  if ~isempty(DATA) && ismethod(DATA, varargin{1})
    macro_helper(varargin{:});
    [varargout{1:nargout}] = eval(sprintf('DATA.%s(varargin{2:end});',varargin{1}));
  else
    [varargout{1:nargout}] = segment_main(varargin{:});
  end
end