function mywaitbarclose(varargin)
%MYWAITBARCLOSE
%
%Closes a mywaitbar structure, and removes waitbar.
%
%See also MYWAITBARSTART, MYWAITBARUPDATE.

%Einar Heiberg

global DATA

if ~isempty(DATA)
  DATA.mywaitbarclose(varargin{:}); %Call the potentially overloaded method.
else
  %Well seems not to be overloaded => 
  maingui.mywaitbarclose(varargin{:});
end;
