function dyssynchrony(varargin)
%Gives dyssynchrony measures such as peaks for segments currently only
%available for SAX as multi long axis doesn't necessarily have the same
%tvector need to create upsample scheme for this.

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

function init(taggroup,type)
global SET DATA

gui = mygui(fullfile('+straintagging','dyssynchrony.fig'));
DATA.GUI.dyssynchrony = gui;
gui.tvec=linspace(0,1,SET(taggroup(1)).TSize);%SET(taggroup(1)).TimeVector;
gui.taggroup=taggroup;
gui.cmap=[1 1 0;...
    1 128/255 128/255;...
    0 1 0;...
    128/255 0 0;...
    0 0 1;...
    0 1 1;...
    1 0 0;...
    1 0 1;...
    .8 .8 .8;...
    240/255 120/255 0;...
    64/255 0 128/255;...
    128/255 64/255 0;...
    0 64/255 0;...
    128/255 128/255 128/255;...
    128/255 128/255 1;...
    0 128/255 128/255;...
    .5 .5 0];

if strcmp(type,'LAX')
 set(gui.handles.typepopupmenu,'String',{'Longitudinal Peak Time';'Longitudinal Strain';'Radial Peak Time';'Radial Strain'})
  include2ch=0;
  include3ch=0;
  include4ch=0;
  
%Get all segments
for no=taggroup
  switch SET(no).ImageViewPlane
    case '2CH'
      segmentalcirc2ch = SET(no).StrainTagging.segmentcirc;
      segmentalrad2ch = SET(no).StrainTagging.segmentrad;
      T=SET(no).TSize;
      include2ch=1;
    case '3CH'
  segmentalcirc3ch = SET(no).StrainTagging.segmentcirc;
      segmentalrad3ch = SET(no).StrainTagging.segmentrad;
      T=SET(no).TSize;
      include3ch=1;
    case '4CH'
      segmentalcirc4ch = SET(no).StrainTagging.segmentcirc;
      segmentalrad4ch = SET(no).StrainTagging.segmentrad;
      T=SET(no).TSize;
      include4ch=1;
  end
end

%calculate global mean strain
if ~include2ch
  globalcirc2ch = NaN*ones(T,1);
  globalrad2ch = NaN*ones(T,1);
  segmentalcirc2ch = NaN*ones(T,7);
  segmentalrad2ch = NaN*ones(T,7);
end

if ~include3ch
  globalcirc3ch = NaN*ones(T,1);
  globalrad3ch = NaN*ones(T,1);
  segmentalcirc3ch = NaN*ones(T,7);
  segmentalrad3ch = NaN*ones(T,7);
end

if ~include4ch
  globalcirc4ch = NaN*ones(T,1);
  globalrad4ch = NaN*ones(T,1);
  segmentalcirc4ch = NaN*ones(T,7);
  segmentalrad4ch = NaN*ones(T,7);
end

bullseyecirc=cell(1,T);%nan(24,4,Tlax);
bullseyeradial=cell(1,T);%nan(24,4,Tlax);

for t=1:T
  bullseyecirc{t} = straintagging.straintagging('getbullseyevalueslax',segmentalcirc2ch,segmentalcirc3ch,segmentalcirc4ch,t,t,t);
  bullseyeradial{t} = straintagging.straintagging('getbullseyevalueslax',segmentalrad2ch,segmentalrad3ch,segmentalrad4ch,t,t,t);
end

else
  set(gui.handles.typepopupmenu,'String',{'Circumferential Peak Time';'Circumferential Strain';'Radial Peak Time';'Radial Strain'})
  no=taggroup;
  T=SET(no).TSize;
  
  for t=1:T
    bullseyecirc{t} = straintagging.straintagging('getbullseyevaluessax',SET(no).StrainTagging.segmentcirc,t,SET(no).StrainTagging.saslices);
    bullseyeradial{t} = straintagging.straintagging('getbullseyevaluessax',SET(no).StrainTagging.segmentrad,t,SET(no).StrainTagging.saslices);
  end
end

for t=1:SET(no).TSize
  bullseyeplot = bullseyeradial{t};
  %basal
  bullseyerad_t(1,t) = mynanmean(bullseyeplot(1:4,4));
  bullseyerad_t(2,t) = mynanmean(bullseyeplot(5:8,4));
  bullseyerad_t(3,t) = mynanmean(bullseyeplot(9:12,4));
  bullseyerad_t(4,t) = mynanmean(bullseyeplot(13:16,4));
  bullseyerad_t(5,t) = mynanmean(bullseyeplot(17:20,4));
  bullseyerad_t(6,t) = mynanmean(bullseyeplot(21:24,4));
  %mid
  bullseyerad_t(7,t) = mynanmean(bullseyeplot(1:4,3));
  bullseyerad_t(8,t) = mynanmean(bullseyeplot(5:8,3));
  bullseyerad_t(9,t) = mynanmean(bullseyeplot(9:12,3));
  bullseyerad_t(10,t) = mynanmean(bullseyeplot(13:16,3));
  bullseyerad_t(11,t) = mynanmean(bullseyeplot(17:20,3));
  bullseyerad_t(12,t) = mynanmean(bullseyeplot(21:24,3));
  %apical
  bullseyerad_t(13,t) = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
  bullseyerad_t(14,t) = mynanmean(bullseyeplot(3:8,2));
  bullseyerad_t(15,t) = mynanmean(bullseyeplot(9:14,2));
  bullseyerad_t(16,t) = mynanmean(bullseyeplot(15:20,2));
  %apex
  bullseyerad_t(17,t) = mynanmean(bullseyeplot(:,1));
  
  bullseyeplot = bullseyecirc{t};
  %basal
  bullseyecirc_t(1,t) = mynanmean(bullseyeplot(1:4,4));
  bullseyecirc_t(2,t) = mynanmean(bullseyeplot(5:8,4));
  bullseyecirc_t(3,t) = mynanmean(bullseyeplot(9:12,4));
  bullseyecirc_t(4,t) = mynanmean(bullseyeplot(13:16,4));
  bullseyecirc_t(5,t) = mynanmean(bullseyeplot(17:20,4));
  bullseyecirc_t(6,t) = mynanmean(bullseyeplot(21:24,4));
  %mid
  bullseyecirc_t(7,t) = mynanmean(bullseyeplot(1:4,3));
  bullseyecirc_t(8,t) = mynanmean(bullseyeplot(5:8,3));
  bullseyecirc_t(9,t) = mynanmean(bullseyeplot(9:12,3));
  bullseyecirc_t(10,t) = mynanmean(bullseyeplot(13:16,3));
  bullseyecirc_t(11,t) = mynanmean(bullseyeplot(17:20,3));
  bullseyecirc_t(12,t) = mynanmean(bullseyeplot(21:24,3));
  %apical
  bullseyecirc_t(13,t) = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
  bullseyecirc_t(14,t) = mynanmean(bullseyeplot(3:8,2));
  bullseyecirc_t(15,t) = mynanmean(bullseyeplot(9:14,2));
  bullseyecirc_t(16,t) = mynanmean(bullseyeplot(15:20,2));
  %apex
  bullseyecirc_t(17,t) = mynanmean(bullseyeplot(:,1));
end

gui.bullseyecirc_t=bullseyecirc_t;
gui.bullseyerad_t=bullseyerad_t;

if isfield(SET(taggroup(1)).StrainTagging,'dyssynchrony')
   
  gui.tfcirc=SET(taggroup(1)).StrainTagging.dyssynchrony.tfcirc;
gui.indcirc=  SET(taggroup(1)).StrainTagging.dyssynchrony.indcirc;
gui.tfrad = SET(taggroup(1)).StrainTagging.dyssynchrony.tfrad;
gui.indrad = SET(taggroup(1)).StrainTagging.dyssynchrony.indrad;
  gui.sections2show = SET(taggroup(1)).StrainTagging.dyssynchrony.include;
  gui.globcircpeaktf=SET(taggroup(1)).StrainTagging.dyssynchrony.globcircpeaktf;
  gui.globradpeaktf=SET(taggroup(1)).StrainTagging.dyssynchrony.globradpeaktf;

else

[~,indcirc]=min(bullseyecirc_t');
[~,indrad]=max(bullseyerad_t');

gui.indcirc=indcirc;
gui.indrad=indrad;

%peaktfs
gui.tfcirc=gui.tvec(indcirc);%SET(no).TimeVector(indcirc);
gui.tfrad=gui.tvec(indrad);%SET(no).TimeVector(indrad);

[~,meancircstrain_ind]=min(nansum(bullseyecirc_t));
[~,meanradstrain_ind]=max(nansum(bullseyerad_t));
gui.globcircpeaktf=meancircstrain_ind./SET(no).TSize;%gui.tvec(meancircstrain_ind);
gui.globradpeaktf=meanradstrain_ind./SET(no).TSize;%gui.tvec(meanradstrain_ind);

gui.sections2show =1:17;
end

%do graphical stuff to gui

for ahaloop=1:17
    [ahastri{ahaloop},pos(ahaloop)] = reportbullseye('aha17nameandpos',ahaloop); %Get name and position of export
end


for i=gui.sections2show
  set(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end

for i=1:17
set(gui.handles.(['radiobutton',num2str(i)]),'String',ahastri{i})
set(gui.handles.(['radiobutton',num2str(i)]),'Callback','straintagging.dyssynchrony(''radio_Callback'')');
set(gui.handles.(['radiobutton',num2str(i)]),'backgroundcolor',gui.cmap(i,:))
end

gui.handles.peakscirc=nan(1,17);
gui.handles.peaksrad=nan(1,17);
gui.handles.graphscirc=nan(1,17);
gui.handles.graphsrad=nan(1,17);
gui.oneatatime=0;

if strcmp(type,'SAX')
  set(gui.handles.unittext,'String', 'circum./rad. [s/T]')
  set(gui.handles.dyssynctext,'String','circum./rad. [std]')
else
  set(gui.handles.unittext,'String', 'longit./rad. [s/T]')
  set(gui.handles.dyssynctext,'String','longit./rad. [std]')
end
  
radio_Callback;

%unclick apex after initiation
if strcmp(type,'SAX')
  set(gui.handles.radiobutton17,'Value',0)
  radio_Callback
end


function setpeak_Callback(type,ind)
global DATA
gui = DATA.GUI.dyssynchrony ;
h=gui.handles.segaxes;
%first find closest cross
switch type
  case 'circ'
    peaks=gui.handles.peakscirc;
    graphs=gui.handles.graphscirc;
  case 'rad'
    peaks=gui.handles.peaksrad;
    graphs=gui.handles.graphsrad;
end

if nargin<2
for i=gui.sections2show
  x(i)=get(peaks(i),'Xdata');
  y(i)=get(peaks(i),'Ydata');
end
[xclick,yclick]=mygetcurrentpoint(h);
[~,ind]=min((xclick-x).^2+(yclick-y).^2);
end

switch type
  case 'circ'
    myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.dyssynchrony(''setpeakmotion'',''circ'',%d)',ind));%d)%sprintf('straintagging.dyssynchrony(''setpeakmotion'',%d)',ind))
    myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.dyssynchrony(''setpeakbuttonup'',''circ'',%d)',ind));%d)%sprintf('straintagging.dyssynchrony(''setpeakmotion'',%d)',ind))
  case 'rad'
    myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.dyssynchrony(''setpeakmotion'',''rad'',%d)',ind));%d)%sprintf('straintagging.dyssynchrony(''setpeakmotion'',%d)',ind))
    myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.dyssynchrony(''setpeakbuttonup'',''rad'',%d)',ind));%d)%sprintf('straintagging.dyssynchrony(''setpeakmotion'',%d)',ind))   
end

function setpeakmotion(type,handleind)%type,ind)
global DATA
gui = DATA.GUI.dyssynchrony ;
h=gui.handles.segaxes;

switch type
  case 'circ'
    peakhandle=gui.handles.peakscirc(handleind);
    graphhandle=gui.handles.graphscirc(handleind);
  case 'rad'
    peakhandle=gui.handles.peaksrad(handleind);
    graphhandle=gui.handles.graphsrad(handleind);
end

[xclicked,~]=mygetcurrentpoint(h);
x=get(graphhandle,'Xdata');
[~,ind]=min((x-xclicked).^2);
y=get(graphhandle,'Ydata');
set(peakhandle,'Xdata',x(ind),'Ydata',y(ind))

for i=gui.sections2show
  %set(gui.handles.(['edit',num2str(i)]),'String',sprintf('%0.2f/%0.2f',gui.tfcirc(i)-gui.globcircpeaktf,gui.tfrad(i)-gui.globradpeaktf))
  set(gui.handles.(['text',num2str(i)]),'String',sprintf('%0.2f/%0.2f',gui.tfcirc(i),gui.tfrad(i)))
end

switch type
  case 'circ'
    gui.tfcirc(handleind)=x(ind);
    gui.indcirc(handleind)=ind;
  case 'rad'
    gui.tfrad(handleind)=x(ind);
    gui.indrad(handleind)=ind;
end
set(gui.handles.dyssyncedit,'String',sprintf('%0.2f/%0.2f',std(gui.tfcirc(gui.sections2show)),std(gui.tfrad(gui.sections2show))));

function setpeakbuttonup(type,handleind)%type,ind)
global DATA
gui = DATA.GUI.dyssynchrony ;
h=gui.handles.segaxes;

switch type
  case 'circ'
    peakhandle=gui.handles.peakscirc(handleind);
    graphhandle=gui.handles.graphscirc(handleind);
  case 'rad'
    peakhandle=gui.handles.peaksrad(handleind);
    graphhandle=gui.handles.graphsrad(handleind);
end

[xclicked,~]=mygetcurrentpoint(h);
x=get(graphhandle,'Xdata');
[~,ind]=min((x-xclicked).^2);
y=get(graphhandle,'Ydata');
set(peakhandle,'Xdata',x(ind),'Ydata',y(ind))

set(gui.handles.fig,'windowbuttonmotionfcn',[]);
set(gui.handles.fig,'windowbuttonupfcn',[])


for i=gui.sections2show
  %set(gui.handles.(['edit',num2str(i)]),'String',sprintf('%0.2f/%0.2f',gui.tfcirc(i)-gui.globcircpeaktf,gui.tfrad(i)-gui.globradpeaktf))
  set(gui.handles.(['text',num2str(i)]),'String',sprintf('%0.2f/%0.2f',gui.tfcirc(i),gui.tfrad(i)))
end
set(gui.handles.dyssyncedit,'String',sprintf('%0.2f/%0.2f',std(gui.tfcirc(gui.sections2show)),std(gui.tfrad(gui.sections2show))));

switch type
  case 'circ'
    gui.tfcirc(handleind)=x(ind);
    gui.indcirc(handleind)=ind;
    myset(gui.handles.fig,'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',%s)','circ'))
  case 'rad'
    gui.tfrad(handleind)=x(ind);
    gui.indrad(handleind)=ind;
    myset(gui.handles.fig,'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',%s)','rad'))
end

function all_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end

radio_Callback

function basal_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
end

for i=1:6
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end
radio_Callback

function mid_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
end

for i=7:12
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end
radio_Callback

function apical_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
end

for i=13:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end
radio_Callback

function anterior_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
end

anteriorinds=[2 8 14];

for i=anteriorinds
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end
radio_Callback

function inferior_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
end

inferiorinds=[5 11 16];

for i=inferiorinds
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end
radio_Callback

function lateral_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
end

lateralinds=[3 4 9 10 15];

for i=lateralinds
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end
radio_Callback

function septal_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end

for i=1:17
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
end

septalinds=[1 6 7 12 13];

for i=septalinds
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end
radio_Callback



function oneatatime_Callback
global DATA
gui = DATA.GUI.dyssynchrony;
oneatatime=get(gui.handles.oneatatimecheckbox,'value');
set(gui.handles.oneatatimecheckbox,'value',oneatatime);
gui.oneatatime=oneatatime;

if oneatatime
  for i=1:17
    val=get(gui.handles.(['radiobutton',num2str(i)]),'Value');
    if val
      ind=i;
      break;
    end
  end
  for i=1:17
    myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
  end
  myset(gui.handles.(['radiobutton',num2str(ind)]),'Value',1);
  gui.lastpressed=ind;
  
  gui.sections2show=ind;
  updateplot;
end

function radio_Callback
global DATA
gui = DATA.GUI.dyssynchrony ;

if gui.oneatatime
  myset(gui.handles.(['radiobutton',num2str(gui.lastpressed)]),'Value',0);
end

tmp=zeros(1,17);
for i=1:17
  tmp(i)=get(gui.handles.(['radiobutton',num2str(i)]),'Value');
end

if gui.oneatatime
  lastpressed=find(tmp);
  if isempty(lastpressed)
    myset(gui.handles.(['radiobutton',num2str(gui.lastpressed)]),'Value',1);
    return;
  else
    myset(gui.handles.(['radiobutton',num2str(lastpressed)]),'Value',1);
    gui.lastpressed=lastpressed;
  end
end

gui.sections2show=find(tmp);
updateplot;


function updateplot
global DATA

gui = DATA.GUI.dyssynchrony ;
myset(gui.handles.fig,'ButtonDownFcn',[])
h=gui.handles.segaxes;
legend(h,'off')

cla(h)
%cla(gui.handles.reportaxes)

set(h,'yticklabelmode','auto')
set(h,'ytickmode','auto')
ylim(h,'auto')
xlim(h,'auto')
% get names of segments

for ahaloop=1:17
    [ahastri{ahaloop},pos(ahaloop)] = reportbullseye('aha17nameandpos',ahaloop); %Get name and position of export
end

for i=gui.sections2show
set(gui.handles.(['text',num2str(i)]),'String',sprintf('%0.2f/%0.2f',gui.tfcirc(i),gui.tfrad(i)))
%set(gui.handles.(['edit',num2str(i)]),'String',sprintf('%0.2f/%0.2f',gui.tfcirc(i)-gui.globcircpeaktf,gui.tfrad(i)-gui.globradpeaktf))
end
set(gui.handles.dyssyncedit,'String',sprintf('%0.2f/%0.2f',std(gui.tfcirc(gui.sections2show)),std(gui.tfrad(gui.sections2show))));

strlist=get(gui.handles.typepopupmenu,'String');
ind=get(gui.handles.typepopupmenu,'Value');
str=strlist{ind};

hold(h,'on')
grid(h,'on');
switch str
  case {'Circumferential Peak Time','Longitudinal Peak Time'}
    set(h,'Ydir','reverse')
    for i=gui.sections2show%1:17
      plot(h,gui.tfcirc(i),i,'k+','markersize',10,'linewidth',3)
    end
    gui.handles.cumplabel=plot(h,[gui.globcircpeaktf,gui.globcircpeaktf],[0,18],'g-','linewidth',2);
    set(h,'ytick',1:17)
    set(h,'yticklabel',ahastri)
    ylim(h,[0,18])
    legend(h,gui.handles.cumplabel,'Cumulative strain peak')
    ylabel(h,'')
    
  case 'Radial Peak Time'
    set(h,'Ydir','reverse')
    for i=gui.sections2show%1:17
      plot(h,gui.tfrad(i),i,'k+','markersize',10,'linewidth',3)
    end
    gui.handles.cumplabel=plot(h,[gui.globradpeaktf,gui.globradpeaktf],[0,18],'g','linewidth',2);
    set(h,'ytick',1:17)
    set(h,'yticklabel',ahastri)
    ylim(h,[0,18])  
    legend(h,gui.handles.cumplabel,'Cumulative strain peak')
    ylabel(h,'')
    
  case {'Circumferential Strain','Longitudinal Strain'}
    %myset(gui.handles.fig,'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',%s)','circ'))
    set(h,'Ydir','Normal')
    for i=gui.sections2show
      gui.handles.graphscirc(i)=plot(h,gui.tvec,gui.bullseyecirc_t(i,:),'.-','linewidth',2,'Color',gui.cmap(i,:));
      gui.handles.peakscirc(i)=plot(h,gui.tfcirc(i),gui.bullseyecirc_t(i,gui.indcirc(i)),'k+','markersize',10,'linewidth',3);
      myset(gui.handles.peakscirc(i),'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',''circ'',%d)',i));
      myset(gui.handles.graphscirc(i),'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',''circ'',%d)',i));
    end
    
    xlim(h,[0,gui.tvec(end)])
    ylabel(h,'%')
  case 'Radial Strain'
    set(h,'Ydir','Normal')
    %myset(gui.handles.fig,'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',%s)','rad'))
    for i=gui.sections2show
      gui.handles.graphsrad(i)=plot(h,gui.tvec,gui.bullseyerad_t(i,:),'.-','linewidth',2,'Color',gui.cmap(i,:));
      gui.handles.peaksrad(i)=plot(h,gui.tfrad(i),gui.bullseyerad_t(i,gui.indrad(i)),'k+','markersize',10,'linewidth',3);
      myset(gui.handles.peaksrad(i),'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',''rad'',%d)',i));
      myset(gui.handles.graphsrad(i),'ButtonDownFcn',sprintf('straintagging.dyssynchrony(''setpeak_Callback'',''rad'',%d)',i));
    end
     
    xlim(h,[0,gui.tvec(end)])
    ylabel(h,'%')
end
xlabel(h,'Heart Cycle [s/T]')
hold(h,'off')



function export_Callback
global DATA SET
%store and close
gui = DATA.GUI.dyssynchrony;
no = gui.taggroup(1);

col=1;
line=1;

outdata{1,1} = 'Patient Name';
outdata{2,1} = 'Patient ID';
outdata{3,1} = 'Heart Rate';
outdata{4,1} = 'Image Type';

outdata{1,2} = SET(no).PatientInfo.Name;
outdata{2,2} = SET(no).PatientInfo.ID;
outdata{3,2} = SET(no).HeartRate;
outdata{4,2} = sprintf('%s %s',SET(no).ImageType,SET(no).ImageViewPlane);

line = 6;

outdata{line,2}='Circ. Peak Time [s/T]';
outdata{line,3}='Rad. Peak Time [s/T]';
outdata{line,4}='Circ. Peak [%]';
outdata{line,5}='Rad. Peak [%]';

for ahaloop=1:17
    [ahastri{ahaloop},pos(ahaloop)] = reportbullseye('aha17nameandpos',ahaloop); %Get name and position of export
end

for i = 1:17
  outdata{line+i,1}=ahastri{i};
end

for i=gui.sections2show;
  outdata{line+i,2}=gui.tfcirc(i);
  outdata{line+i,3}=gui.tfrad(i);
  outdata{line+i,4}=gui.bullseyecirc_t(i,gui.indcirc(i));
  outdata{line+i,5}=gui.bullseyerad_t(i,gui.indrad(i));
end

outdata{line+19,2}='Dyssynchrony [std]';
outdata{line+20,1}='Circ. Peak Time';
outdata{line+21,1}='Rad. Peak Time';
outdata{line+20,2}=std(gui.tfcirc(gui.sections2show));
outdata{line+21,2}=std(gui.tfrad(gui.sections2show));

segment('cell2clipboard',outdata);

function save_Callback
global DATA SET
%store and close
gui = DATA.GUI.dyssynchrony;
for no = gui.taggroup
  SET(no).StrainTagging.dyssynchrony.tfcirc=gui.tfcirc;
  SET(no).StrainTagging.dyssynchrony.indcirc=gui.indcirc;
  SET(no).StrainTagging.dyssynchrony.tfrad=gui.tfrad;
  SET(no).StrainTagging.dyssynchrony.indrad=gui.indrad;
  SET(no).StrainTagging.dyssynchrony.include=gui.sections2show;
  SET(no).StrainTagging.dyssynchrony.globcircpeaktf = gui.globcircpeaktf;
  SET(no).StrainTagging.dyssynchrony.globradpeaktf=gui.globradpeaktf;
end

try
  DATA.GUI.dyssynchrony = close(DATA.GUI.dyssynchrony );
catch   %#ok<CTCH>
  DATA.GUI.dyssynchrony =[];
  delete(gcbf);
end 

function close_Callback
global DATA

try
  DATA.GUI.dyssynchrony = close(DATA.GUI.dyssynchrony );
catch   %#ok<CTCH>
  DATA.GUI.dyssynchrony =[];
  delete(gcbf);
end
