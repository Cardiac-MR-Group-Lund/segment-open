function h = mywaitbarupdate(varargin)
%H = MYWAITBARDATE(H) 
%
%Update a mywaitbarstructure, and increment 
%counter one step. Graphic is updated depends on
%setting when intializing the structure.
%
%NOTE: MYWAITBARSTART is obsoleted and recommended is to replace with
%MYWAITBAR instead.
%
%See also MYWAITBARSTART, MYWAITBARSET, MYWAITBARCLOSE.

global DATA %#ok<*GVMIS> 

if ~isempty(DATA)
  h = DATA.mywaitbarupdate(varargin{:}); %Call the potentially overloaded method.
else
  %Well seems not to be overloaded => 
  h = maingui.mywaitbarupdate(varargin{:});
end
