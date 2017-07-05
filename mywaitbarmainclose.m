function mywaitbarmainclose(varargin)
%MYWAITBARCLOSE
%
%Closes a mywaitbar structure, and removes waitbar.
%
%See also MYWAITBARSTART, MYWAITBARUPDATE.

%Einar Heiberg

global DATA

if ~isempty(DATA)
  DATA.mywaitbarmainclose(varargin{:}); %Call the potentially overloaded method.
else
  %Well seems not to be overloaded => 
  maingui.mywaitbarmainclose(varargin{:});
end;
