function h = mywaitbarstart(varargin)
%H = MYWAITBARSTART(NUMITERATIONS,STRING,[IGNORELIMIT],FIGHANDLE)
%
%Similar to waitbar, but first input argument is
%the total number of iterations that will be performed.
%Second input argument is string that will be displayed, 
%and third optional argument is a limit that if fewer
%iterations that this is performed, then there will be 
%no graphical output. Another difference to standard waitbar
%is that no more than 20 graphical steps are displayed to 
%minimize CPU consumption. Function returns a handle to a 
%MYWAITBAROBJECT. FIGHANDLE is optional and indicates alignment
%
%See also MYWAITBARUPDATE, MYWAITBARCLOSE.

%Einar Heiberg

global DATA

if ismethod(DATA,'mywaitbarmainstart')
  h = DATA.mywaitbarstart(varargin{:}); %Call the potentially overloaded method.
else
  %Well seems not to be overloaded => 
  h = maingui.mywaitbarstart(varargin{:});
end;

