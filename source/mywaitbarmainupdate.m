function h = mywaitbarmainupdate(varargin)
%H = MYWAITBARDATE(H) 
%
%Update a mywaitbarstructure, and increment 
%counter one step. Graphic is updated depends on
%setting when intializing the structure.
%
%See also MYWAITBARSTART, MYWAITBARCLOSE.

global DATA

if ~isempty(DATA)
  h = DATA.mywaitbarmainupdate(varargin{:}); %Call the potentially overloaded method.
else
  %Well seems not to be overloaded => 
  h = maingui.mywaitbarmainupdate(varargin{:});
end;
