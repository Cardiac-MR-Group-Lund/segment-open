function h = mywaitbarset(varargin)
%H = MYWAITBARSET(frac,h) 
%
%Update a mywaitbarstructure, and set it to fraction
%counter one step. Graphic is updated depends on
%setting when intializing the structure.
%
%See also MYWAITBARSTART, MYWAITBARUPDATE, MYWAITBARCLOSE.

global DATA

if ~isempty(DATA)
  h = DATA.mywaitbarset(varargin{:}); %Call the potentially overloaded method.
else
  %Well seems not to be overloaded => 
  h = maingui.mywaitbarset(varargin{:});
end
