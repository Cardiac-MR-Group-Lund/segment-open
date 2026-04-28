function [varargout] = createfunctions(varargin)
% Functions for creating graphicsobjects and content in image panels
% Klas

%#ok<*GVMIS>

%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%----------------------------------------------
function addcolorbar(panel,contrast,brightness) 
%----------------------------------------------
%Add colorbar to axes in panel
global DATA SET  

if DATA.GUISettings.ShowColorbar
   no = DATA.ViewPanels(panel);
  if nargin < 2
    contrast = SET(no).IntensityMapping.Contrast;
    brightness = SET(no).IntensityMapping.Brightness;
  end
  %addcolorbar
  h = colorbar(DATA.Handles.imageaxes(panel),'East');
 
  [win, level] = calcfunctions('con2win',contrast, brightness, no);
  if ~isempty(win) && ~isempty(level)
    set(DATA.Handles.imageaxes(panel),'CLim',[(level-win/2) (level+win/2)])
  else
     colorbar(DATA.Handles.imageaxes(panel),'off')
     return;
  end
  if isempty(SET(no).Colormap)
    colormap(h,'gray');
  else
    colormap(h,SET(no).Colormap);
  end
  h.Color = DATA.GUISettings.ForegroundColor;
  %h.Limits =
else
  %removecolorbar
  colorbar(DATA.Handles.imageaxes(panel),'off')
end

%---------------------
function createplotoverlay 
%-----------------------
%Creates graphical objects in volume and flow axes.
global DATA

set([DATA.Handles.volumeaxes,DATA.Handles.flowaxes,DATA.Handles.timebaraxes],'FontSmoothing','off');
volax = DATA.Handles.volumeaxes;
flowax = DATA.Handles.flowaxes;
timeax = DATA.Handles.timebaraxes;

lc = DATA.GUISettings.BarColor;
volaxcolor = DATA.GUISettings.VolumeAxesColor;
esstring = dprintf('ES');
edstring = dprintf('ED');
%time plot
ylim(DATA.Handles.timebaraxes, [0,1]); %always between 0 and 1
DATA.Handles.timebar = line('Parent',timeax,'XData',nan,'YData',nan,'Color',lc,'linewidth',2);
DATA.Handles.timebaraxeshelpers = line('Parent',timeax,'XData',nan,'YData',nan,'Color','k');
DATA.Handles.edtimebartext = text('Parent',timeax,'position',[nan nan],'color',lc,'String',edstring);
DATA.Handles.estimebartext = text('Parent',timeax,'position',[nan nan],'color',lc,'String',esstring);
DATA.Handles.edtimebarline = line('Parent',timeax,'XData',nan,'YData',nan,'Color',lc,'linewidth',2);
DATA.Handles.estimebarline = line('Parent',timeax,'XData',nan,'YData',nan,'Color',lc,'linewidth',2);
set([DATA.Handles.estimebartext DATA.Handles.estimebarline],'ButtonDownFcn','segment(''esedtimebar_Buttondown'',''es'')');
set([DATA.Handles.edtimebartext DATA.Handles.edtimebarline],'ButtonDownFcn','segment(''esedtimebar_Buttondown'',''ed'')');
set([DATA.Handles.timebar,timeax,DATA.Handles.timebaraxeshelpers],'buttondownFcn','segment(''timebar_Buttondown'')');
timestring = makeunitstring(dprintf('Time'),'ms');
xlabel(timeax,timestring,'Color',volaxcolor,'Units','normalized');
ylabel(timeax,timestring,'Color',volaxcolor,'Units','normalized');
set(get(DATA.Handles.timebaraxes,'ylabel'),'Visible','off');
set(timeax,'ytick',[]);
set(volax,'XColor',volaxcolor,'YColor',volaxcolor);

%volume plot
DATA.Handles.timebarlv = line('Parent',volax,'XData',nan,'YData',nan,'Color',lc,'linewidth',2);
DATA.Handles.volumecurve = line('Parent',volax,'XData',nan,'YData',nan,'Color','r','markersize',5,'linestyle','-.');
DATA.Handles.masscurve = line('Parent',volax,'XData',nan,'YData',nan,'Color','b','markersize',5,'linestyle','-.');
DATA.Handles.rvvcurve = line('Parent',volax,'XData',nan,'YData',nan,'Color','m','markersize',5,'linestyle','-.');
DATA.Handles.volumeaxeshelpers = line('Parent',volax,'XData',nan,'YData',nan,'Color','k');
DATA.Handles.estext = text('Parent',volax,'position',[nan nan],'color',lc,'String',esstring);
DATA.Handles.edtext = text('Parent',volax,'position',[nan nan],'color',lc,'String',edstring);
DATA.Handles.esline = line('Parent',volax,'XData',nan,'YData',nan,'Color',lc);
DATA.Handles.edline = line('Parent',volax,'XData',nan,'YData',nan,'Color',lc);
DATA.Handles.lvmtext = text('Parent',volax,'position',[nan nan],'color','b','fontsize',8,'String','LVM','verticalalignment','bottom');
set([DATA.Handles.estext,DATA.Handles.esline],'buttondownFcn','segment(''esed_Buttondown'',''es'')')
set([DATA.Handles.edtext,DATA.Handles.edline],'buttondownFcn','segment(''esed_Buttondown'',''ed'')')
set([volax, DATA.Handles.timebarlv,DATA.Handles.volumecurve,DATA.Handles.masscurve,DATA.Handles.rvvcurve,DATA.Handles.volumeaxeshelpers],...
  'buttondownFcn','segment(''volumeaxes_Buttondown'')')
volstr = makeunitstring(dprintf('Volume'),'ml');
ylabel(volax,volstr,'Color',volaxcolor,'units','normalized');

xlabel(volax,timestring,'Color',volaxcolor,'units','normalized');
set(volax,'XColor',volaxcolor,'YColor',volaxcolor);

%flow plot
if ~isempty(flowax)
  DATA.Handles.timebarflow = line('Parent',flowax,'XData',nan,'YData',nan,'Color',lc,'linewidth',2);
  DATA.Handles.outerbarsflow = line('Parent',flowax,'XData',nan,'YData',nan,'Color',lc,'linewidth',2);
  DATA.Handles.zerolineflow = line('parent', flowax,'xdata',[0 10000000],'ydata',[0 0],'color','k','linestyle',':');
  DATA.Handles.flowcurve = line('Parent',flowax,'XData',nan,'YData',nan,'Color','b','markersize',5,'linestyle','-.','Linewidth',2);
  DATA.Handles.flowtext = text(nan,nan,'[Eddy]', ...
    'Parent',flowax,'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top','Color',DATA.GUISettings.ForegroundColor);
  DATA.Handles.flowaxeshelpers = line('Parent',flowax,'XData',nan,'YData',nan,'Color','k');
  set([flowax,DATA.Handles.zerolineflow,DATA.Handles.flowcurve,DATA.Handles.flowaxeshelpers,DATA.Handles.outerbarsflow],...
    'ButtonDownFcn','segment(''flowaxes_Buttondown'')');
  flowstr = makeunitstring(dprintf('Flow'),'ml/s');
  ylabel(flowax,flowstr,'Color',volaxcolor,'units','normalized');
  xlabel(flowax,timestring,'Color',volaxcolor,'units','normalized');
  set(flowax,'XColor',volaxcolor,'YColor',volaxcolor);
end

%---------------------------------------------
function createaxesandimagehandles(maxnumaxes)
%---------------------------------------------
%Creates the maximum used imageaxes and boxaxes in segment. The field boxaxes is probably unnecessary.
global DATA

% %Sizes in which we place the panels
% left = DATA.GUISettings.LeftGapWidth; %0.12;
% right = 1-DATA.GUISettings.RightGapWidth-0.02;  %Based on that the report panel can never be more than 220.
% bottom = DATA.GUISettings.BottomGapHeight; %0.013;
% top = 1-DATA.GUISettings.TopGapHeight;
% width = right-left;
% height = top-bottom;

%Create image area with grid for boxes.
DATA.Handles.boxaxes = axes('position',...
  [0 0 0 0],'color',[0.1,0.1,0.1],'YTick', [], 'XTick', [],...
  'parent',DATA.fig,'visible','on');

%This is the line in boxaxes which actually makeup the box around all the panels
DATA.Handles.box = line('parent',DATA.Handles.boxaxes,'XData',nan,'YData',nan,'linewidth',2);

%create empty axes
h = [];
for i = 1:maxnumaxes
  h = [h axes('position',...
    [0 0 0 0],...
    'parent',DATA.fig,'visible','off')]; %#ok<AGROW>
end

%Store handles
axis(h,'ij')
set(h,'fontsmoothing','off')
DATA.Handles.imageaxes = h;
set(DATA.Handles.imageaxes,'tag','imageaxes');

if DATA.issegment3dp
  %Create panel handle for viewport (used by myvolshow to host volshow
  DATA.Handles.viewportpanel = uipanel(...
    'Position',h(2).Position,...
    'Parent',DATA.fig,...
    'BackgroundColor',DATA.GUISettings.BackgroundColor,...
    'Visible','off');

  %Create panel handle for viewportgraphics (used by myvolshow to add graphics objects)
  DATA.Handles.viewportgraphicspanel = uipanel(...
    'Position',h(2).Position,...
    'Parent',DATA.fig,...
    'BackgroundColor',DATA.GUISettings.BackgroundColor,...
    'Visible','off');

  %Create sliders (used in 3DPrint for scrolling)
  h = [];
  for i = 1:maxnumaxes
    h = [h uicontrol('Style','slider','Parent',DATA.fig,'Visible','off',...
      'Interruptible','off','BusyAction','cancel')]; %#ok<AGROW>
  end
  DATA.Handles.sliders = h;

end

%Initiate imagehandles in all axes
DATA.Handles.imagehandles = gobjects(1,maxnumaxes);
for i = 1:maxnumaxes
  DATA.Handles.imagehandles(i) = image('parent',DATA.Handles.imageaxes(i),'cdata',[]);
  set(DATA.Handles.imagehandles(i),'tag',sprintf('imagehandle%02d',i));
end

%------------------------
function create3dpoverlay 
%------------------------
global DATA
%Draws and initiates all necessary handles for 3dp viewpanels

%for 3DP the maximum number of panels are
maxnumpanels = 4;
maxsurfaceintersections = 15; %max surface intersections per panel
emptygraphics = gobjects(1,maxnumpanels);
DATA.Handles.rbline = emptygraphics;
DATA.Handles.rgline = emptygraphics;
DATA.Handles.grline = emptygraphics;
DATA.Handles.gbline =    emptygraphics;
DATA.Handles.brline = emptygraphics;
DATA.Handles.bgline = emptygraphics;
DATA.Handles.threedpcontour = emptygraphics;
DATA.Handles.surfaceintersectline = cell(1,maxnumpanels);
DATA.Handles.boxline = emptygraphics;

h = gobjects(1,maxsurfaceintersections); %reserve memory
for i = 1:maxnumpanels
  DATA.Handles.rbline(i) = line('parent',DATA.Handles.imageaxes(i),'XData',[nan nan],'YData',[nan nan],'Color','b','lineStyle','-','tag','rb'); %add lines, later correct in lineupdate;
  DATA.Handles.rgline(i) = line('parent',DATA.Handles.imageaxes(i),'XData',[nan nan],'YData',[nan nan],'Color','g','lineStyle','-','tag','rg');
  
  DATA.Handles.grline(i) = line('parent',DATA.Handles.imageaxes(i),'XData',[nan nan],'YData',[nan nan],'Color','r','lineStyle','-','tag','gr');
  DATA.Handles.gbline(i) = line('parent',DATA.Handles.imageaxes(i),'XData',[nan nan],'YData',[nan nan],'Color','b','lineStyle','-','tag','gb');
  
  DATA.Handles.brline(i) = line('parent',DATA.Handles.imageaxes(i),'XData',[nan nan],'YData',[nan nan],'Color','r','lineStyle','-','tag','br');
  DATA.Handles.bgline(i) = line('parent',DATA.Handles.imageaxes(i),'XData',[nan nan],'YData',[nan nan],'Color','g','lineStyle','-','tag','bg');

  DATA.Handles.boxline(i) = line('parent',DATA.Handles.imageaxes(i),'XData',[nan nan],'YData',[nan nan],'Color','g','lineStyle','--','linewidth',2,'tag','boxline'); %Box around for cropping

  %Create maxsurfaceintersections line handles
  for sloop = 1:maxsurfaceintersections
    h(sloop) = line('parent',DATA.Handles.imageaxes(i),'XData',NaN,'YData',NaN,'Color','g','lineStyle','none','marker','.','tag',sprintf('surfaceintersection%02d',sloop));
  end
  DATA.Handles.surfaceintersectline{i} = h;

  hold(DATA.Handles.imageaxes(i),'on')
  warning('off') %we are fine with the contour not being rendered and only created here
  
  [~,DATA.Handles.threedpcontour(i)] = contour(DATA.Handles.imageaxes(i),nan(2),[127 127]);
  DATA.Handles.threedpcontour(i).LineColor = 'g';
  DATA.Handles.threedpcontour(i).Tag = 'threedpcontour';
  
  warning('on')
  hold(DATA.Handles.imageaxes(i),'off')
end

%--------------------------------
function createimageoverlayall3dp
%--------------------------------
%This function initiates all overlay graphics such as text contours etc for
%Segment 3DPrint

global DATA 

maxnumpanels = 4;

imageaxeshandles = DATA.Handles.imageaxes;
emptygraphics = gobjects(1,maxnumpanels);
DATA.Handles.text = emptygraphics;

%text entries
DATA.Handles.text = (copyobj(text(3,3,'','HorizontalAlignment', 'left', ...
  'VerticalAlignment', 'top','color',[1,1,1],  'Interpreter', 'none','fontsmoothing','off'),imageaxeshandles))';

%----------------------
function createimage3dp 
%----------------------
%This function initiates all overlay graphics such as text contours etc.
% profile on

global DATA 

%maxnumpanels = length(DATA.Handles.imageaxes);
imageaxeshandles = DATA.Handles.imageaxes;

%emptygraphics = gobjects(1,maxnumpanels);

DATA.imagecursor = line('xdata',nan,'ydata',nan,'color','y');
DATA.Handles.text = (copyobj(text(3,3,'','HorizontalAlignment', 'left', ...
  'VerticalAlignment', 'top','color',[1,1,1],  'Interpreter', 'none','fontsmoothing','off'),imageaxeshandles))';

DATA.Handles.selectedslice = (copyobj(line('xdata',nan,'ydata',nan,'color','y','linewidth',1.5),imageaxeshandles))';
DATA.Handles.selectedframe = (copyobj(line('xdata',nan,'ydata',nan,'color','w','linewidth',1.5),imageaxeshandles))';

%-----------------------------
function createimageoverlayall 
%-----------------------------
%This function initiates all overlay graphics such as text contours etc.
% profile on

global DATA 

maxnumpenobjects = DATA.GeneralPenSettings.getmaxnumobjects();
maxnumpanels = length(DATA.Handles.imageaxes);
imageaxeshandles = DATA.Handles.imageaxes;
% the initial test should be with contours measurement lines and measurement text as these
% are the most complex cases.
DATA.imagecursor = gobjects(1,1);

emptygraphics = gobjects(1,maxnumpanels);
DATA.Handles.endocontour = emptygraphics;
DATA.Handles.epicontour = emptygraphics;
DATA.Handles.rvendocontour = emptygraphics;
DATA.Handles.rvepicontour = emptygraphics;
DATA.Handles.endointerp = emptygraphics;
DATA.Handles.epiinterp = emptygraphics;
DATA.Handles.rvendointerp = emptygraphics;
DATA.Handles.rvepiinterp = emptygraphics;
DATA.Handles.lacontour = emptygraphics;
DATA.Handles.racontour = emptygraphics;
DATA.Handles.lainterp = emptygraphics;
DATA.Handles.rainterp = emptygraphics;
DATA.Handles.generalpeninterp = gobjects(maxnumpanels,maxnumpenobjects);
DATA.Handles.generalpencontour = gobjects(maxnumpanels,maxnumpenobjects);
DATA.Handles.scarcontour = emptygraphics;
DATA.Handles.marcontour = emptygraphics;
DATA.Handles.weightedscarcontour = emptygraphics; %Weighted scar contour
DATA.Handles.mocontour = emptygraphics; %microvascular obstruction contour
DATA.Handles.moextentcontour = emptygraphics; %microvascular obstruction contour
DATA.Handles.selectedslice = emptygraphics; %handle to yellow box sourrounding selected slices
DATA.Handles.selectedframe = emptygraphics;
DATA.Handles.cursor = gobjects(1); % this is the temporary line used when drawing with the pen tool

DATA.Handles.centercross = emptygraphics; %handle to the center cross, visible in Segment

DATA.Handles.endocontourintersection = emptygraphics;
DATA.Handles.epicontourintersection = emptygraphics;
DATA.Handles.rvendocontourintersection = emptygraphics;
DATA.Handles.rvepicontourintersection = emptygraphics;
DATA.Handles.text = emptygraphics;

DATA.Handles.planeintersection = gobjects(maxnumpanels);

%DATA.Handles.roi = gobjects(length(DATA.ViewPanels),50); %cannot have more than 20 rois in a panel

DATA.Handles.broi = emptygraphics; %
DATA.Handles.croi = emptygraphics; %
DATA.Handles.groi = emptygraphics; %
DATA.Handles.rroi = emptygraphics; %
DATA.Handles.kroi = emptygraphics; %
DATA.Handles.mroi = emptygraphics; %
DATA.Handles.yroi = emptygraphics; %
DATA.Handles.wroi = emptygraphics; %
DATA.Handles.roicurrent = emptygraphics; %Current selected roi

DATA.Handles.laavplane = emptygraphics;
DATA.Handles.laaxis = emptygraphics;

DATA.Handles.rvavplane = emptygraphics;

DATA.Handles.measurement = emptygraphics;
DATA.Handles.measurementoutsideplane = emptygraphics;
DATA.Handles.measurementtext = gobjects(maxnumpanels,20);%cannot have more than 20 measurement texts in a panel
DATA.Handles.measurementtextline = gobjects(maxnumpanels,20);
if isequal(DATA.ProgramName,'Segment 3DPrint') %Different handling of points in the softwares
  DATA.Handles.point = gobjects(maxnumpanels,20);%cannot have more than 20 pointgroups in a panel
else
  DATA.Handles.point = emptygraphics;
end
DATA.Handles.pointtext = gobjects(maxnumpanels,20);%cannot have more than 20 Point texts in a panel
DATA.Handles.orthoanglehandle = gobjects(1);

%imagecursor is global used in toggleplaceholdermotion to draw custom pen
DATA.imagecursor = gobjects(1);
DATA.imagecursor = line('xdata',nan,'ydata',nan,'color','y');

%then we need to assert it's a line, text etc

%the cursor and orthoanglehandle parents are set on the fly. But so that
%they arent initiated in loading interfaces where
DATA.Handles.cursor = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(1));
DATA.Handles.orthoanglehandle = line('xdata',nan,'ydata',nan,'color','m','marker','.','MarkerSize',20,'parent',imageaxeshandles(1));
orangecolor = helperfunctions('getorangecolor');
%%% fix to copy object handles instead of creating in a loop
lacolor = DATA.GeneralPenSettings.getatriumcolor('la');
racolor = DATA.GeneralPenSettings.getatriumcolor('ra');
DATA.Handles.endocontour = (copyobj(line('xdata',nan,'ydata',nan,'color','r'),imageaxeshandles))';
DATA.Handles.epicontour = (copyobj(line('xdata',nan,'ydata',nan,'color','g'),imageaxeshandles))';
DATA.Handles.rvendocontour = (copyobj(line('xdata',nan,'ydata',nan,'color','m'),imageaxeshandles))';
DATA.Handles.rvepicontour = (copyobj(line('xdata',nan,'ydata',nan,'color','c'),imageaxeshandles))';
DATA.Handles.lacontour = (copyobj(line('xdata',nan,'ydata',nan,'color',lacolor),imageaxeshandles))';
DATA.Handles.racontour = (copyobj(line('xdata',nan,'ydata',nan,'color',racolor),imageaxeshandles))';
DATA.Handles.measurement = (copyobj(line('xdata',nan,'ydata',nan,'color',orangecolor,'marker','+'),imageaxeshandles))';
DATA.Handles.measurementoutsideplane = (copyobj(line('xdata',nan,'ydata',nan,'color','k','marker','+','linestyle','none'),imageaxeshandles))';
DATA.Handles.point = (copyobj(line('xdata',nan,'ydata',nan,'color','y','marker','+','linestyle','none'),imageaxeshandles))';
DATA.Handles.point3D = (copyobj(line('xdata',nan,'ydata',nan,'color','y','marker','o','linestyle','none'),imageaxeshandles))';
DATA.Handles.ortholines = (copyobj(line('xdata',nan,'ydata',nan,'color','w'),imageaxeshandles))';
DATA.Handles.selectedortholine = (copyobj(line('xdata',nan,'ydata',nan,'color','y'),imageaxeshandles))';
DATA.Handles.centercross = (copyobj(line('xdata',nan,'ydata',nan,'color','w','marker','+','MarkerSize',5,'linestyle','none'),imageaxeshandles))';
DATA.Handles.laavplane = (copyobj(line('xdata',nan,'ydata',nan,'color',lacolor,'linestyle',':','marker','o','LineWidth',1.25),imageaxeshandles))';
DATA.Handles.laaxis = (copyobj(line('xdata',nan,'ydata',nan,'color',lacolor,'linestyle',':','marker','.','LineWidth',1.25),imageaxeshandles))';

DATA.Handles.raavplane = (copyobj(line('xdata',nan,'ydata',nan,'color',racolor,'linestyle',':','marker','o','LineWidth',1.25),imageaxeshandles))';

%rois and measurements are special
DATA.Handles.roicurrent = (copyobj(line('xdata',nan,'ydata',nan,'color','b','linewidth',2),imageaxeshandles))';
DATA.Handles.broi = (copyobj(line('xdata',nan,'ydata',nan,'color','b','linewidth',1),imageaxeshandles))';
DATA.Handles.mroi = (copyobj(line('xdata',nan,'ydata',nan,'color','m','linewidth',1),imageaxeshandles))';
DATA.Handles.croi = (copyobj(line('xdata',nan,'ydata',nan,'color','c','linewidth',1),imageaxeshandles))';
DATA.Handles.groi = (copyobj(line('xdata',nan,'ydata',nan,'color','g','linewidth',1),imageaxeshandles))';
DATA.Handles.yroi = (copyobj(line('xdata',nan,'ydata',nan,'color',[1 0.5 0],'linewidth',1),imageaxeshandles))'; %yellow is actually orange
DATA.Handles.wroi = (copyobj(line('xdata',nan,'ydata',nan,'color','w','linewidth',1),imageaxeshandles))';
DATA.Handles.rroi = (copyobj(line('xdata',nan,'ydata',nan,'color','r','linewidth',1),imageaxeshandles))';
DATA.Handles.kroi = (copyobj(line('xdata',nan,'ydata',nan,'color','k','linewidth',1),imageaxeshandles))';

maxnumroitext = 35;
DATA.Handles.roitext = gobjects(maxnumpanels,maxnumroitext);
for j  = 1:maxnumroitext
  DATA.Handles.roitext(:,j) = (copyobj(text(nan,nan,'','HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle','color',[1,1,1],'clipping','on','fontsmoothing','off'),imageaxeshandles))';
end

for j  = 1:20 
  DATA.Handles.measurementtext(:,j) = (copyobj(text(nan,nan,'','FontWeight','bold','interpreter','none','HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle','color',orangecolor,'clipping','on','fontsmoothing','off'),imageaxeshandles))';
  DATA.Handles.measurementtextline(:,j) = (copyobj(line('xdata',nan,'ydata',nan,'color',orangecolor,'linestyle',':','LineWidth',1.25),imageaxeshandles))';
  DATA.Handles.pointtext(:,j) = (copyobj(text(nan,nan,'','FontWeight','bold','interpreter','none','HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle','color','w','clipping','on','fontsmoothing','off'),imageaxeshandles))';
end
for j  = 1:maxnumpanels
    DATA.Handles.planeintersection(:,j) = (copyobj(line('xdata',nan,'ydata',nan,'color','w'),imageaxeshandles))';
end

%text entries
DATA.Handles.text = (copyobj(text(3,3,'','HorizontalAlignment', 'left', ...
  'VerticalAlignment', 'top','color',[1,1,1],  'Interpreter', 'none','fontsmoothing','off'),imageaxeshandles))';

%general pen
color = DATA.GeneralPenSettings.getdefaultcolor();
for loop = 1:maxnumpenobjects
  DATA.Handles.generalpencontour(:,loop) = (copyobj(line('xdata',nan,'ydata',nan,'color',color),imageaxeshandles));
  DATA.Handles.generalpeninterp(:,loop) = (copyobj(line('xdata',nan,'ydata',nan,'marker','.','color',color,'MarkerSize',12,'linestyle','none'),imageaxeshandles));
  DATA.Handles.generalpentext(:,loop) = (copyobj(text(nan,nan,'','HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle','color',[1,1,1],'clipping','on','fontsmoothing','off'),imageaxeshandles));
end

%selected slices handle
DATA.Handles.selectedslice = (copyobj(line('xdata',nan,'ydata',nan,'color','y','linewidth',1.5),imageaxeshandles))';
DATA.Handles.selectedframe = (copyobj(line('xdata',nan,'ydata',nan,'color','w','linewidth',1.5),imageaxeshandles))';
%%% end of fix: copy objects

for i  = 1:maxnumpanels
%   DATA.Handles.generalpencontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y');
%   DATA.Handles.endocontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','r');
%   DATA.Handles.epicontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','g');
%   DATA.Handles.rvendocontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','m');
%   DATA.Handles.rvepicontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','c');
%   DATA.Handles.measurement(i) = line('xdata',nan,'ydata',nan,'color','y','marker','+','parent',DATA.Handles.imageaxes(i));
%   DATA.Handles.point(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y','marker','+','linestyle','none');
%   DATA.Handles.point3D(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y','marker','o','linestyle','none');
%   DATA.Handles.ortholines(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w');
%   DATA.Handles.selectedortholine(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y');
%   DATA.Handles.centercross(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w','marker','+','MarkerSize',5,'linestyle','none');
%   for j  = 1:maxnumpanels
%     DATA.Handles.planeintersection(i,j) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w');
%   end
  
%   %rois and measurements are special
%   DATA.Handles.roicurrent(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','b','linewidth',2);
%   DATA.Handles.broi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','b','linewidth',1);
%   DATA.Handles.mroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','m','linewidth',1);
%   DATA.Handles.croi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','c','linewidth',1);
%   DATA.Handles.groi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','g','linewidth',1);
%   DATA.Handles.yroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y','linewidth',1);
%   DATA.Handles.wroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w','linewidth',1);
%   DATA.Handles.rroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','r','linewidth',1);
%   DATA.Handles.kroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','k','linewidth',1);
  
 
%   for j  = 1:10
%     for c = 'rgbymkwc'
%       DATA.Handles.([c,'roitext'])(i,j) = text(nan,nan,'','HorizontalAlignment', 'right', ...
%         'VerticalAlignment', 'middle','color',[1,1,1], 'parent',DATA.Handles.imageaxes(i),'clipping','on','fontsmoothing','off');
%     end
%   end
%   for j  = 1:20 
%     DATA.Handles.measurementtext(i,j) = text(nan,nan,'','FontWeight','bold','interpreter','none','HorizontalAlignment', 'left', ...
%       'VerticalAlignment', 'middle','color','y', 'parent',DATA.Handles.imageaxes(i),'clipping','on','fontsmoothing','off');
%     DATA.Handles.pointtext(i,j) = text(nan,nan,'','FontWeight','bold','interpreter','none','HorizontalAlignment', 'left', ...
%       'VerticalAlignment', 'middle','color','w', 'parent',DATA.Handles.imageaxes(i),'clipping','on','fontsmoothing','off');
%   end
  
  %plot commands delete axes content if not hold is toggled on
  hold(imageaxeshandles(i),'on')
  
  warning('off') %we are fine with the contour not being rendered and only created here
  [~,DATA.Handles.scarcontour(i)] = contour(imageaxeshandles(i),nan(2),[0.5 0.5]);
  DATA.Handles.scarcontour(i).LineColor = [1 1 0];
  
  [~,DATA.Handles.weightedscarcontour(i)] = contour(imageaxeshandles(i),nan(2),[0.5 0.5]);
  DATA.Handles.weightedscarcontour(i).LineColor = [1 0.5 0.5];
  
  [~,DATA.Handles.mocontour(i)] = contour(imageaxeshandles(i),nan(2),[0.5 0.5]);
  DATA.Handles.mocontour(i).LineColor = [1 0 0];
  DATA.Handles.mocontour(i).LineStyle = ':';
  
  [~,DATA.Handles.moextentcontour(i)] = contour(imageaxeshandles(i),nan(2),[0.99 0.99]);
  DATA.Handles.moextentcontour(i).LineColor = [1 0 0];
  
  [~,DATA.Handles.marcontour(i)] = contour(imageaxeshandles(i),nan(2),[0.5 0.5]);
  DATA.Handles.marcontour(i).LineColor = [1 1 1];
  warning('on')
  
  %interp allocation
  DATA.Handles.endointerp(i) = line('xdata',nan,'ydata',nan,'parent',imageaxeshandles(i),'marker','.','color','r','MarkerSize',12,'linestyle','none');
  DATA.Handles.epiinterp(i) = line('xdata',nan,'ydata',nan,'parent',imageaxeshandles(i),'marker','.','color','g','MarkerSize',12,'linestyle','none');
  DATA.Handles.rvendointerp(i) = line('xdata',nan,'ydata',nan,'parent',imageaxeshandles(i),'marker','.','color','m','MarkerSize',12,'linestyle','none');
  DATA.Handles.rvepiinterp(i) = line('xdata',nan,'ydata',nan,'parent',imageaxeshandles(i),'marker','.','color','c','MarkerSize',12,'linestyle','none');
  DATA.Handles.lainterp(i) = line('xdata',nan,'ydata',nan,'parent',imageaxeshandles(i),'marker','.','color',lacolor,'MarkerSize',12,'linestyle','none');
  DATA.Handles.rainterp(i) = line('xdata',nan,'ydata',nan,'parent',imageaxeshandles(i),'marker','.','color',racolor,'MarkerSize',12,'linestyle','none');
  DATA.Handles.generalpeninterp(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'marker','.','color','y','MarkerSize',12,'linestyle','none');

  %intersections
  DATA.Handles.endocontourintersection(i) = line(imageaxeshandles(i),nan,nan,'marker','.','color','r','MarkerSize',5,'linestyle','none');
  DATA.Handles.epicontourintersection(i) = line(imageaxeshandles(i),nan,nan,'marker','.','color','g','MarkerSize',5,'linestyle','none');
  DATA.Handles.rvendocontourintersection(i) = line(imageaxeshandles(i),nan,nan,'marker','.','color','m','MarkerSize',5,'linestyle','none');
  DATA.Handles.rvepicontourintersection(i) = line(imageaxeshandles(i),nan,nan,'marker','.','color','c','MarkerSize',5,'linestyle','none');
  hold(imageaxeshandles(i),'off')
  
%   %text entries
%   DATA.Handles.text(i) = text(3,3,'','HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top','color',[1,1,1], 'parent',DATA.Handles.imageaxes(i), 'Interpreter', 'none','fontsmoothing','off');
%   
%   %selected slices handle
%   DATA.Handles.selectedslice(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y','linewidth',1.5);
%   DATA.Handles.selectedframe(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w','linewidth',1.5);
  
end
% profile report

%---------------------
function createimages 
%-----------------------
%Initiates all imagehandles
global DATA

%Initiate imagehandles in all axes
DATA.Handles.imagehandles = gobjects(1,length(DATA.ViewPanels));
for i = 1:length(DATA.ViewPanels)
  DATA.Handles.imagehandles(i) = image('parent',DATA.Handles.imageaxes(i),'cdata',[]);
end

%---------------------
function createplot 
%-----------------------
global DATA


DATA.Handles.plot1 = line('xdata',nan,'ydata',nan,'color','r','parent',DATA.Handles.plotaxes);
DATA.Handles.plot2 = line('xdata',nan,'ydata',nan,'color','g','parent',DATA.Handles.plotaxes);
DATA.Handles.plot3 = line('xdata',nan,'ydata',nan,'color','b','parent',DATA.Handles.plotaxes);
DATA.Handles.plot4 = line('xdata',nan,'ydata',nan,'color','c','parent',DATA.Handles.plotaxes);

%---------------------
function createimageoverlay 
%-----------------------
%This function initiates all overlay graphics such as text contours etc.
global DATA

% the initial test should be with contours measurement lines and measurement text as these
% are the most complex cases.
emptygraphics = gobjects(1,length(DATA.ViewPanels));
DATA.Handles.endocontour = emptygraphics;
DATA.Handles.epicontour = emptygraphics;
DATA.Handles.rvendocontour = emptygraphics;
DATA.Handles.rvepicontour = emptygraphics;
DATA.Handles.endointerp = emptygraphics;
DATA.Handles.epiinterp = emptygraphics;
DATA.Handles.rvendointerp = emptygraphics;
DATA.Handles.rvepiinterp = emptygraphics;
DATA.Handles.lacontour = emptygraphics;
DATA.Handles.racontour = emptygraphics;
DATA.Handles.lainterp = emptygraphics;
DATA.Handles.rainterp = emptygraphics;
DATA.Handles.generalpeninterp = emptygraphics;
DATA.Handles.generalpencontour = emptygraphics;
DATA.Handles.scarcontour = emptygraphics;
DATA.Handles.marcontour = emptygraphics;
DATA.Handles.weightedscarcontour = emptygraphics; %Weighted scar contour
DATA.Handles.mocontour = emptygraphics; %microvascular obstruction contour
DATA.Handles.moextentcontour = emptygraphics; %microvascular obstruction contour
DATA.Handles.selectedslice = emptygraphics; %handle to yellow box sourrounding selected slices
DATA.Handles.selectedframe = emptygraphics;
DATA.Handles.cursor = gobjects(1); % this is the temporary line used when drawing with the pen tool

DATA.Handles.endocontourintersection = emptygraphics;
DATA.Handles.epicontourintersection = emptygraphics;
DATA.Handles.rvendocontourintersection = emptygraphics;
DATA.Handles.rvepicontourintersection = emptygraphics;
DATA.Handles.text = emptygraphics;

DATA.Handles.planeintersection = gobjects(length(DATA.ViewPanels));

DATA.Handles.centercross = emptygraphics; %handle to the center cross, visible in Segment

%DATA.Handles.roi = gobjects(length(DATA.ViewPanels),50); %cannot have more than 20 rois in a panel
DATA.Handles.broi = emptygraphics; %
DATA.Handles.croi = emptygraphics; %
DATA.Handles.groi = emptygraphics; %
DATA.Handles.rroi = emptygraphics; %
DATA.Handles.kroi = emptygraphics; %
DATA.Handles.mroi = emptygraphics; %
DATA.Handles.yroi = emptygraphics; %
DATA.Handles.wroi = emptygraphics; %
DATA.Handles.roicurrent = emptygraphics; %Current selected roi
% DATA.Handles.broitext = gobjects(length(DATA.ViewPanels),20);%
% DATA.Handles.croitext = gobjects(length(DATA.ViewPanels),20); %
% DATA.Handles.groitext = gobjects(length(DATA.ViewPanels),20);%
% DATA.Handles.rroitext = gobjects(length(DATA.ViewPanels),20);%
% DATA.Handles.kroitext = gobjects(length(DATA.ViewPanels),20);%
% DATA.Handles.mroitext = gobjects(length(DATA.ViewPanels),20);%
% DATA.Handles.yroitext = gobjects(length(DATA.ViewPanels),20);%
% DATA.Handles.wroitext = gobjects(length(DATA.ViewPanels),20);%

% DATA.Handles.roitext = gobjects(length(DATA.ViewPanels),20);
DATA.Handles.measurement = emptygraphics;
DATA.Handles.measurementoutsideplane = emptygraphics;
DATA.Handles.measurementtext = gobjects(length(DATA.ViewPanels),20);%cannot have more than 20 measurement texts in a panel
DATA.Handles.point = emptygraphics;
DATA.Handles.pointtext = gobjects(length(DATA.ViewPanels),20);%cannot have more than 20 Point texts in a panel
DATA.Handles.orthoanglehandle = gobjects(1);
DATA.Handles.point3D = emptygraphics; %handle to the points in 3D

%then we need to assert it's a line, text etc

%the cursor and orthoanglehandle parents are set on the fly. But so that
%they arent initiated in loading interfaces where
DATA.Handles.cursor = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(1));
DATA.Handles.orthoanglehandle = line('xdata',nan,'ydata',nan,'color','m','marker','.','MarkerSize',20,'parent',DATA.Handles.imageaxes(1));

%Use preferences linewidth if there
if ~isempty(DATA.Pref.LineWidth)
  linewidth = DATA.Pref.LineWidth;
else
  linewidth = 1;
end
orangecolor = helperfunctions('getorangecolor');
for i  = 1:length(DATA.ViewPanels)
  DATA.Handles.generalpencontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y');
  DATA.Handles.endocontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','r', 'LineWidth', linewidth);
  DATA.Handles.epicontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','g', 'LineWidth', linewidth);
  DATA.Handles.rvendocontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','m', 'LineWidth', linewidth);
  DATA.Handles.rvepicontour(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','c', 'LineWidth', linewidth);
  DATA.Handles.measurement(i) = line('xdata',nan,'ydata',nan,'color',orangecolor,'marker','+','parent',DATA.Handles.imageaxes(i));
  DATA.Handles.measurementoutsideplane(i) = line('xdata',nan,'ydata',nan,'color','k','marker','+','parent',DATA.Handles.imageaxes(i)); %used to plot points outside the slice
  DATA.Handles.point(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y','marker','+','linestyle','none');
  DATA.Handles.ortholines(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w');
  DATA.Handles.selectedortholine(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y');
  
  for j  = 1:length(DATA.ViewPanels)
    DATA.Handles.planeintersection(i,j) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w');
  end
  
  %rois and measurements are special
  DATA.Handles.roicurrent(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','b','linewidth',linewidth+1);
  DATA.Handles.broi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','b','linewidth',linewidth);
  DATA.Handles.mroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','m','linewidth',linewidth);
  DATA.Handles.croi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','c','linewidth',linewidth);
  DATA.Handles.groi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','g','linewidth',linewidth);
  DATA.Handles.yroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y','linewidth',linewidth);
  DATA.Handles.wroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w','linewidth',linewidth);
  DATA.Handles.rroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','r','linewidth',linewidth);
  DATA.Handles.kroi(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','k','linewidth',linewidth);
  
  for j  = 1:20
% % %     for c = 'rgbymkwgc'
% % %       DATA.Handles.([c,'roitext'])(i,j) = text(nan,nan,'','HorizontalAlignment', 'right', ...
% % %         'VerticalAlignment', 'middle','color',[1,1,1], 'parent',DATA.Handles.imageaxes(i),'clipping','on');
% % %     end
    DATA.Handles.measurementtext(i,j) = text(nan,nan,'','FontWeight','bold','interpreter','none','HorizontalAlignment', 'left', ...
      'VerticalAlignment', 'middle','color','y', 'parent',DATA.Handles.imageaxes(i),'clipping','on');
    DATA.Handles.pointtext(i,j) = text(nan,nan,'','FontWeight','bold','interpreter','none','HorizontalAlignment', 'left', ...
      'VerticalAlignment', 'middle','color','w', 'parent',DATA.Handles.imageaxes(i),'clipping','on');
  end
  
  %plot commands delete axes content if not hold is toggled on
  hold(DATA.Handles.imageaxes(i),'on')
  
  warning('off') %we are fine with the contour not being rendered and only created here
  [~,DATA.Handles.scarcontour(i)] = contour(DATA.Handles.imageaxes(i),nan(2),[0.5 0.5]);
  DATA.Handles.scarcontour(i).LineColor = [1 1 0];
  DATA.Handles.scarcontour(i).LineWidth = linewidth;
  
  [~,DATA.Handles.weightedscarcontour(i)] = contour(DATA.Handles.imageaxes(i),nan(2),[0.5 0.5]);
  DATA.Handles.weightedscarcontour(i).LineColor = [1 0.5 0.5];
  
  [~,DATA.Handles.mocontour(i)] = contour(DATA.Handles.imageaxes(i),nan(2),[0.5 0.5]);
  DATA.Handles.mocontour(i).LineColor = [1 0 0];
  DATA.Handles.mocontour(i).LineStyle = ':';
  
  [~,DATA.Handles.moextentcontour(i)] = contour(DATA.Handles.imageaxes(i),nan(2),[0.99 0.99]);
  DATA.Handles.moextentcontour(i).LineColor = [1 0 0];
  
  [~,DATA.Handles.marcontour(i)] = contour(DATA.Handles.imageaxes(i),nan(2),[0.5 0.5]);
  DATA.Handles.marcontour(i).LineColor = [1 1 1];
  DATA.Handles.marcontour(i).LineWidth = linewidth;
  warning('on')
  
  %interp allocation
  DATA.Handles.endointerp(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'marker','.','color','r','MarkerSize',12,'linestyle','none');
  DATA.Handles.epiinterp(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'marker','.','color','g','MarkerSize',12,'linestyle','none');
  DATA.Handles.rvendointerp(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'marker','.','color','m','MarkerSize',12,'linestyle','none');
  DATA.Handles.rvepiinterp(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'marker','.','color','c','MarkerSize',12,'linestyle','none');
  DATA.Handles.generalpeninterp(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'marker','.','color','y','MarkerSize',12,'linestyle','none');
  
  %intersections
  DATA.Handles.endocontourintersection(i) = line(DATA.Handles.imageaxes(i),nan,nan,'marker','.','color','r','MarkerSize',4,'linestyle','none');
  DATA.Handles.epicontourintersection(i) = line(DATA.Handles.imageaxes(i),nan,nan,'marker','.','color','g','MarkerSize',4,'linestyle','none');
  DATA.Handles.rvendocontourintersection(i) = line(DATA.Handles.imageaxes(i),nan,nan,'marker','.','color','m','MarkerSize',4,'linestyle','none');
  DATA.Handles.rvepicontourintersection(i) = line(DATA.Handles.imageaxes(i),nan,nan,'marker','.','color','c','MarkerSize',4,'linestyle','none');
  hold(DATA.Handles.imageaxes(i),'off')
  
  %text entries
  DATA.Handles.text(i) = text(3,3,'','HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top','color',[1,1,1], 'parent',DATA.Handles.imageaxes(i), 'Interpreter', 'none');
  
  %selected slices handle
  DATA.Handles.selectedslice(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','y','linewidth',1.5);
  DATA.Handles.selectedframe(i) = line('xdata',nan,'ydata',nan,'parent',DATA.Handles.imageaxes(i),'color','w','linewidth',1.5);
  
end

%--------------------------------
function im = createviewim(panel) 
%--------------------------------
%Initiates all imagehandles and creates DATA.viewim which holds the current slice for all timeframes mapped ready for display.
global DATA SET

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

%Check if calcium score overlay
calciumim = [];
if isfield(SET(no),'CT') && isfield(SET(no).CT,'CalciumMask')
  calciumoverlay = true;
else
  calciumoverlay = false;
end

% check for papillary muscle overlay
stateandicon=viewfunctions('iconson','hidelv');
ispapmusc = false;
if not(stateandicon{1}) && (~isempty(SET(no).PapillaryIM))
 ispapmusc = true;
end

%Check if RGB
isRGB = false;
sz = size(SET(no).IM);
if numel(sz) == 6 && sz(6) == 3
  isRGB = true;
  ispapmusc = false; %disable pap and calcium overlays
  calciumoverlay = false;
end
  
%Later also have a hide function for this :-)

switch DATA.ViewPanelsType{panel}
  case 'trans3DP'
    DATA.ViewIM{panel} = 'cache not used';
    %set zoom and aspectratio for axes
    viewfunctions('updatezoomandaspectratio',panel);
    return

  case 'sag3DP'
    DATA.ViewIM{panel} = 'cache not used';
    %set zoom and aspectratio for axes
    viewfunctions('updatezoomandaspectratio',panel);
    return

  case 'cor3DP'
    DATA.ViewIM{panel} = 'cache not used';
    %set zoom and aspectratio for axes
    viewfunctions('updatezoomandaspectratio',panel);
    return

  case 'speedim'
    DATA.ViewIM{panel} = 'cache not used';
    %set zoom and aspectratio for axes
    viewfunctions('updatezoomandaspectratio',panel);
    return

  case 'viewport'
    DATA.ViewIM{panel} = [];
    return

  case {'one','orth'}
    currentslice = SET(no).CurrentSlice;
    if isempty(currentslice)
      currentslice = 1; %to make sure there is always a current slice
      SET(no).CurrentSlice = 1;
    end
    im = SET(no).IM(:,:,:,currentslice,:,:);
    if calciumoverlay
      calciumim = SET(no).CT.CalciumMask(:,:,:,currentslice);
    end
    if ispapmusc
      papimage = SET(no).PapillaryIM(:,:,:,currentslice);
    end

  case {'montage','montagesegmented','montagerow'}
    if strcmp(DATA.ViewPanelsType{panel},'montagesegmented')%and currentslice
      slicestoinclude = segment('getmontagesegmentedslices',no);
      [rows,cols] = calcfunctions('calcrowscols',no,length(slicestoinclude)); %.cols,.rows
    else
      slicestoinclude = 1:SET(no).ZSize;
      [rows,cols] = calcfunctions('calcrowscols',no); %.cols,.rows
    end
    %Set the correct number or rows and cols
    if strcmp(DATA.ViewPanelsType{panel},'montagerow')
      n = rows*cols;
      rows = min(rows,2);
      cols = ceil(n/rows);
    end
    matrix = [rows cols];
    %Allocate space
    imdim = zoomfunctions.getxysize(no,panel);
    if isRGB
      im = zeros([imdim.XSize*matrix(1) imdim.YSize*matrix(2) SET(no).TSize 1 1 3]);
    else
      if isa(SET(no).IM,'int16')
        im = int16(zeros([imdim.XSize*matrix(1) imdim.YSize*matrix(2) SET(no).TSize]));
      else
        im = zeros([imdim.XSize*matrix(1) imdim.YSize*matrix(2) SET(no).TSize]);
      end
    end
    if ispapmusc
      papimage = false([imdim.XSize*matrix(1) imdim.YSize*matrix(2) SET(no).TSize]);
    end
    for tloop = 1:SET(no).TSize
      for i = 1:numel(slicestoinclude)
        zloop = slicestoinclude(i);
        col = 1 + mod(i-1,matrix(2));
        r = ceil(i/matrix(2));
        if isRGB
          im(...
            (1+(r-1)*imdim.XSize):(r*imdim.XSize),...
            (1+(col-1)*imdim.YSize):(col*imdim.YSize),tloop,:,:,:)...
            = SET(no).IM(imdim.XStart:imdim.XEnd,:,tloop,zloop,:,:,:);
        else
          im(...
            (1+(r-1)*imdim.XSize):(r*imdim.XSize),...
            (1+(col-1)*imdim.YSize):(col*imdim.YSize),tloop)...
            = SET(no).IM(imdim.XStart:imdim.XEnd,imdim.YStart:imdim.YEnd,tloop,zloop);
        end
        if ispapmusc
          papimage(...
            (1+(r-1)*imdim.XSize):(r*imdim.XSize),...
            (1+(col-1)*imdim.YSize):(col*imdim.YSize),tloop)...
            = SET(no).PapillaryIM(imdim.XStart:imdim.XEnd,imdim.YStart:imdim.YEnd,tloop,zloop);
        end
      end
    end
  case 'gla'
    im = viewfunctions('getgla',no);
  case 'vla'
    im = permute(SET(no).IM(:,SET(no).VLA.slice,:,:),[4 1 3 2]);
  case 'hla'
    im = permute(SET(no).IM(SET(no).HLA.slice,:,:,:),[4 2 3 1]);
end

if ispapmusc && ~isempty(papimage)
  im(papimage) = DATA.GUISettings.PapilarColor;
end

%do all resizing based on scale
if (scale==1) %exclude for 1, not sure if saves speed
  outim = im; %Matlab performs copy by reference so should be fine with JIT
else  
  imxsz = size(im,1)*scale;
  imysz = size(im,2)*scale;
  outim = imresize(im,[imxsz imysz],'bilinear');
  if ~isempty(calciumim)
      calciumim = imresize(im,[imxsz imysz],'nearest'); 
  end
end

if calciumoverlay && ~isempty(calciumim)
  im  = calcfunctions('remapuint8viewim',outim,no); %Remap
  rim = im;
  gim = im;
  bim = im;
  
  %colortable for the classes in calciummask
  colortable = [];
  colortable(1).class = uint8(1); %yellow
  colortable(1).rgb = uint8([255 255 0]); 
  colortable(2).class = uint8(2);
  colortable(2).rgb = uint8([0 255 255]); %cyan
  colortable(3).class = uint8(3);
  colortable(3).rgb = uint8([0 255 0]); %green
  colortable(4).class = uint8(4);
  colortable(4).rgb = uint8([255 255 255]); %white
  colortable(5).class = uint8(5);
  colortable(5).rgb = uint8([80 80 80]); %gray
   colortable(6).class = uint8(6);
  colortable(6).rgb = uint8([0 0 255]); %blue
  colortable(7).class = uint8(7);
  colortable(7).rgb = uint8([255 0 0]); %red
   colortable(8).class = uint8(8);
  colortable(8).rgb = uint8([204 0 204]); %lila
  for loop = 1:length(colortable)
    rim(calciumim==colortable(loop).class) = colortable(loop).rgb(1);
    gim(calciumim==colortable(loop).class) = colortable(loop).rgb(2);
    bim(calciumim==colortable(loop).class) = colortable(loop).rgb(3);    
  end
  DATA.ViewIM{panel} = cat(5,rim,gim,bim);
else
  %This is the default
  DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',outim,no);
end

%set zoom and aspectratio for axes
viewfunctions('updatezoomandaspectratio',panel);
viewfunctions('updatetextposition',panel)
if panel==DATA.CurrentPanel
  drawfunctions('drawselectedframe',panel)
end

