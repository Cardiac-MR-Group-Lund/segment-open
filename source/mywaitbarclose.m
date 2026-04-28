function mywaitbarclose(varargin)
%MYWAITBARCLOSE
%
%Closes a mywaitbar structure, and removes waitbar.
%
%NOTE: MYWAITBARSTART is obsoleted and recommended is to replace with
%MYWAITBAR instead.
%
%See also MYWAITBARSTART, MYWAITBARUPDATE.

%Einar Heiberg

global DATA %#ok<*GVMIS> 

if ~isempty(DATA)
  DATA.mywaitbarclose(varargin{:}); %Call the potentially overloaded method.
else
  %Well seems not to be overloaded => 
  maingui.mywaitbarclose(varargin{:});
end
