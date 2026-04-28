function reportvolume(varargin)
% Volume reporting tools
% Moved out from segment_main by Nisse Lundahl

%#ok<*GVMIS>
%Invoke subfunction
feval(varargin{:}); % FEVAL switchyard

%--------------------------
function curve_Callback(no,path,name)
%--------------------------
%GUI for reporting the volume curve.
global SET NO DATA

if nargin == 0
  no = NO;
end

if nargin <2
  path = getpreferencespath;
end

if nargin <3
  name =  'lvvolumeplot.png';
end

fig=figure('visible','off');
set(fig,'Name',dprintf('Volume curve'),'numbertitle','off');
t = 1000*((0:(SET(no).TSize-1))*SET(no).TIncr);
h = plot(t,SET(no).LVV,'r.-');
set(h,'linewidth',2);
hold on;
xlim([t(1) t(end)]);
% if (SET(no).EPV(1)>SET(no).LVV(1))
%   h = plot(t,SET(no).EPV-SET(no).LVV,'g.-');set(h,'linewidth',2);
%   h = plot([t(1) t(end)],[mean(SET(no).EPV-SET(no).LVV) mean(SET(no).EPV-SET(no).LVV)],'g:');
%   %set(h,'linewidth',2);
% end;
legendstr = makeunitstring(dprintf('LV volume'),'ml');
legend(legendstr,'location','SouthEast'); %'Myocard volume [ml]','Mean myocard volume [ml]',
hold off;
grid on;
% h = title('Volume versus time');set(h,'fontsize',14);
timestr = makeunitstring(dprintf('Time'),'ms');
h = xlabel(timestr);
set(h,'fontsize',14);
volstr = makeunitstring(dprintf('Volume'),'ml');
h = ylabel(volstr);
set(h,'fontsize',14);
set(0, 'currentfigure', fig);
currentFolder = DATA.SegmentFolder;
cd(path)
print(fig,'-dpng', name)
cd(currentFolder)
delete(fig)

%---------------------
function loop_Callback
%---------------------
%Reporter for volume loop %
 
% (not yet implemented)

myfailed('Report volume loop is not yet implemented.');