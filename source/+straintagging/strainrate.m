function varargout = strainrate(varargin)
%Gives strainrate measures 

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

 
%---------------------------------------
function init(taggroup,type)
%----------------------------------------
global SET DATA

gui = mygui(fullfile('+straintagging','strainrate.fig'));
DATA.GUI.strainrate = gui;
if strcmp(type,'LAX')
  taggrouptmp=nan(1,3);
  %sort taggroup so we can always use gui.views2show for indexes in LAX
  for tagno=taggroup
    switch SET(tagno).ImageViewPlane
      case '2CH'
        taggrouptmp(1)=tagno;
      case '3CH'
        taggrouptmp(2)=tagno;
      case '4CH'
        taggrouptmp(3)=tagno;
    end
  end
  taggrouptmp(isnan(taggrouptmp))=[];
  taggroup=taggrouptmp;
  gui.taggroup=taggroup;
else
  gui.taggroup=taggroup;
end
gui.oneatatime=0;
gui.lastpressed=[];
gui.tvec = SET(taggroup(1)).StrainTagging.strainrateTvect;
gui.cmap=[1 1 0;...
    0 0 1;...
    0 1 0;...
    0 1 1;...
    1 0 0;...
    1 0 1;...
    .5 .5 0;...
    .8 .8 .8;...
    240/255 120/255 0;...
    64/255 0 128/255;...
    128/255 64/255 0;...
    0 64/255 0;...
    128/255 128/255 128/255;...
    128/255 128/255 1;...
    0 128/255 128/255;...
    128/255 0 0;...
    1 128/255 128/255];

srstack=[];
for no = taggroup
  if isfield(SET(no).StrainTagging,'SR')
    srstack=no;
    break;
  end
end

if isempty(srstack)  %|| isfield(SET(no).StrainTagging.SR,'updated')
  tmpstruct = calcstrainrate(type, taggroup);
else
  tmpstruct=SET(srstack).StrainTagging.SR;
end

%combine structs 
 f = fieldnames(tmpstruct);
 for i = 1:length(f)
    gui.(f{i}) = tmpstruct.(f{i});
 end

switch type
  case 'LAX'
    
    %gui.tvec = %linspace(0,1,SET(taggroup(1)).TSize-2);
    strcell={'2CH','3CH','4CH'};
    gui.radios2show =1:sum(gui.include);%1:length(taggroup);
  set(gui.handles.typepopupmenu ,'String',{'Longitudinal Strain Rate';...
    'Radial Strain Rate';...
    'Mean Longitudinal Strain Rate';...
    'Mean Radial Strain Rate';...
    'Mean Longitudinal Strain with Tangents';...
    'Mean Radial Strain with Tangents';...
    'Longitudinal Strain with Tangents';...
    'Radial Strain with Tangents'});
  used=find(gui.include);
  for i=gui.radios2show
    set(gui.handles.(['radiobutton',num2str(i)]),'String',strcell{used(i)});
    set(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
    set(gui.handles.(['radiobutton',num2str(i)]),'backgroundcolor',gui.cmap(i,:))
  end
%   for i=1:length(gui.include)%gui.radios2show
%     set(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
%     set(gui.handles.(['radiobutton',num2str(i)]),'backgroundcolor',gui.cmap(i,:))
%   end
    set(gui.handles.circtext,'String','Longit. [%/s]');
    
    set(gui.handles.globalcirctext,'String','Longit. [%/(s/T)]');
    set(gui.handles.globalradtext,'String','Rad. [%/(s/T)]');
  case 'SAX'
    
     %graphical
      set(gui.handles.typepopupmenu ,'String',{'Circumferential Strain Rate';...
        'Radial Strain Rate';...
        'Mean Circumferential Strain Rate';...
        'Mean Radial Strain Rate';...
        'Mean Circumferential Strain with Tangents';...
        'Mean Radial Strain with Tangents';...
        'Circumferential Strain with Tangents';...
        'Radial Strain with Tangents'});
      usedslices=find(gui.include);
       gui.radios2show=1:length(usedslices);
       for i=gui.radios2show
         set(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
         set(gui.handles.(['radiobutton',num2str(i)]),'String',sprintf('Slice %d',usedslices(i)))
         set(gui.handles.(['radiobutton',num2str(i)]),'backgroundcolor',gui.cmap(i,:))
       end
        
    set(gui.handles.circtext,'String','Circ. [%/s]');
    
    set(gui.handles.globalcirctext,'String','Circ. [%/s]');
    set(gui.handles.globalradtext,'String','Rad. [%/s]');
    set(gui.handles.autoviewpeakpushbutton, 'String','Auto Detect Slice Peaks')
end

%graphical
for r = gui.radios2show(end)+1:17
  set(gui.handles.(['radiobutton' num2str(r)]),'visible', 'off')
  set(gui.handles.(['text' num2str(r)]),'visible', 'off')
    set(gui.handles.(['text' num2str(17+r)]),'visible', 'off')
    set(gui.handles.(['text' num2str(2*17+r)]),'visible', 'off')
    set(gui.handles.(['text' num2str(3*17+r)]),'visible', 'off')
end

% gui.popupind=get(gui.handles.typepopupmenu,'Value');
  
radio_Callback;

%--------------------------------------------------------------
function resetpeaks_Callback(type)
%--------------------------------------------------------------
global DATA

gui=DATA.GUI.strainrate;

switch type
  case 'global'
    [gui.globalradup, gui.globalradupind] = max(gui.globalsrrad);
    [gui.globalraddown, gui.globalraddownind] = min(gui.globalsrrad);
    [gui.globalcircup, gui.globalcircupind] = max(gui.globalsrcirc);
    [gui.globalcircdown, gui.globalcircdownind] = min(gui.globalsrcirc); 
  case 'views'
    [gui.radup, gui.radupind] = max(gui.srrad,[],2);
    [gui.raddown, gui.raddownind] = min(gui.srrad,[],2);
    [gui.circup, gui.circupind] = max(gui.srcirc,[],2);
    [gui.circdown, gui.circdownind] = min(gui.srcirc,[],2);
end

updateplot


%---------------------------------------------------------------
function [sr] = calcstrainrate(type, taggroup)
%---------------------------------------------------------------
%Calculates strainrate from stack group
global SET


switch type
  case 'LAX'
    sz=[3, SET(taggroup(1)).TSize-2];
    sr.srcirc =nan(sz);
    sr.srrad = nan(sz);
    sr.include=zeros(1,3);
    for i = 1:length(taggroup)
      switch SET(taggroup(i)).ImageViewPlane
        case '2CH'
          if ~all(all(isnan(SET(taggroup(i)).StrainTagging.strainratecircum)))
            sr.srcirc(1,:) = squeeze(nanmean(SET(taggroup(i)).StrainTagging.strainratecircum,2));
            sr.srrad(1,:) = squeeze(nanmean(SET(taggroup(i)).StrainTagging.strainraterad,2));
            sr.include(1)=1;
          else
            sr.include(1)=0;
          end
        case '3CH'
          if ~all(all(isnan(SET(taggroup(i)).StrainTagging.strainratecircum)))
          sr.srcirc(end-1,:) = squeeze(nanmean(SET(taggroup(i)).StrainTagging.strainratecircum,2));
          sr.srrad(end-1,:) = squeeze(nanmean(SET(taggroup(i)).StrainTagging.strainraterad,2));
          sr.include(2)=1;
          else
            sr.include(2)=0;
          end
        case '4CH'
          if ~all(all(isnan(SET(taggroup(i)).StrainTagging.strainratecircum)))
          sr.srcirc(end,:) = squeeze(nanmean(SET(taggroup(i)).StrainTagging.strainratecircum,2));
          sr.srrad(end,:) = squeeze(nanmean(SET(taggroup(i)).StrainTagging.strainraterad,2));       
          sr.include(3)=1;
          else
            sr.include(3)=0;
          end
      end
    end
    sr.srcirc=sr.srcirc(find(sr.include),:);
    sr.srrad=sr.srrad(find(sr.include),:);
    
  case 'SAX'
      sr.srcirc = squeeze(nanmean(SET(taggroup(1)).StrainTagging.strainratecircum,2))';
      sr.srrad = squeeze(nanmean(SET(taggroup(1)).StrainTagging.strainraterad,2))';
      
      tmp1=[];
      tmp2=[];
      sr.include=ones(1,SET(taggroup(1)).ZSize);
      for i = 1:SET(taggroup(1)).ZSize
        if ~all(isnan(sr.srcirc(i,2:end)))
          tmp1=[tmp1;sr.srcirc(i,:)];
          tmp2=[tmp2;sr.srrad(i,:)];
        else
          sr.include(i)=0;
        end
      end
      sr.srcirc=tmp1;
      sr.srrad=tmp2;
end

%Global Strain
  if strcmp(SET(taggroup(1)).ImageViewPlane,'Short-axis')
    sr.globalsrrad = nanmean(sr.srrad,1);
    sr.globalsrcirc = nanmean(sr.srcirc,1);
  else
    %global strainrate normalised if LAX images
    normstrainratecirc=[];
    normstrainraterad =[];
    
    for tagno = taggroup
      normstrainratecirc = [normstrainratecirc, conv2([1 0 -1],1,SET(tagno).StrainTagging.globalcirc,'valid')/(2/(SET(tagno).TSize-1))];
      normstrainraterad = [normstrainraterad, conv2([1 0 -1],1,SET(tagno).StrainTagging.globalrad,'valid')/(2/(SET(tagno).TSize-1))];
    end
    
    sr.globalsrcirc=nanmean(normstrainratecirc,2);
    sr.globalsrrad=nanmean(normstrainraterad,2);
  end
  
  [sr.globalradup, sr.globalradupind] = max(sr.globalsrrad);
  [sr.globalraddown, sr.globalraddownind] = min(sr.globalsrrad);
  [sr.globalcircup, sr.globalcircupind] = max(sr.globalsrcirc);
  [sr.globalcircdown, sr.globalcircdownind] = min(sr.globalsrcirc);

[sr.radup, sr.radupind] = max(sr.srrad,[],2);
[sr.raddown, sr.raddownind] = min(sr.srrad,[],2);
[sr.circup, sr.circupind] = max(sr.srcirc,[],2);
[sr.circdown, sr.circdownind] = min(sr.srcirc,[],2);

%---------------------------------------------------------------
function setpeak_Callback(type,dir)
%--------------------------------------------------------------
global DATA
gui = DATA.GUI.strainrate ;
h=gui.handles.segaxes;

if nargin < 2
  switch type
    case 'circ'
        tmp=[gui.handles.peakscircup,gui.handles.peakscircdown];
    case 'rad'
        tmp=[gui.handles.peaksradup,gui.handles.peaksraddown];
    case 'globalcirc'
       tmp=[gui.handles.globalpeakcircup,gui.handles.globalpeakcircdown];
    case 'globalrad'
       tmp=[gui.handles.globalpeakradup,gui.handles.globalpeakraddown];
  end
  
  if ~any(strcmp(type,{'globalrad','globalcirc'}))
    x=nan(1,gui.radios2show(end)*2);
    y=nan(1,gui.radios2show(end)*2);
    
    for i = [gui.views2show,gui.views2show+gui.radios2show(end)]
      x(i)=get(tmp(i),'Xdata');
      y(i)=get(tmp(i),'Ydata');
    end
    
    [xclick,yclick]=mygetcurrentpoint(h);
    [~,ind]=min((xclick-x).^2+(yclick-y).^2);
    
    if ind<=gui.radios2show(end)
      dir='up';
    else
      dir='down';
    end
  else
    for i = 1:2
      x(i)=get(tmp(i),'Xdata');
      y(i)=get(tmp(i),'Ydata');
    end
    
    [xclick,yclick]=mygetcurrentpoint(h);
    [~,ind]=min((xclick-x).^2+(yclick-y).^2);
    
    if ind==1
      dir = 'up';
    else
      dir = 'down';
    end
  end
end
  
  
  switch type
    case 'circ'
      switch dir
        case 'up'
          peakhandle=gui.handles.peakscircup;
        case 'down'
          peakhandle=gui.handles.peakscircdown;
      end
    case 'rad'
      switch dir
        case 'up'
          peakhandle=gui.handles.peaksradup;
        case 'down'
          peakhandle=gui.handles.peaksraddown;
      end
    case 'globalcirc'
      switch dir
        case 'up'
          peakhandle=gui.handles.globalpeakcircup;
        case 'down'
          peakhandle=gui.handles.globalpeakcircdown;
      end
    case 'globalrad'
      switch dir
        case 'up'
          peakhandle=gui.handles.globalpeakradup;
        case 'down'
          peakhandle=gui.handles.globalpeakraddown;
      end
  end
  
  if ~strcmp(type,{'globalcirc','globalrad'})
    for i = gui.views2show
      x(i)=get(peakhandle(i),'Xdata');
      y(i)=get(peakhandle(i),'Ydata');
    end
  else
    x=get(peakhandle,'Xdata');
    y=get(peakhandle,'Ydata');
  end
  
  [xclick,yclick]=mygetcurrentpoint(h);
  [~,ind]=min((xclick-x).^2+(yclick-y).^2);
  
  switch type
    case 'circ'
      switch dir
        case 'down'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''circ'',''down'',%d)',ind));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''circ'',''down'',%d)',ind));
        case 'up'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''circ'',''up'',%d)',ind));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''circ'',''up'',%d)',ind));
      end
    case 'rad'
      switch dir
        case 'down'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''rad'',''down'',%d)',ind));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''rad'',''down'',%d)',ind));
        case 'up'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''rad'',''up'',%d)',ind));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''rad'',''up'',%d)',ind));
      end
    case 'globalcirc'
      switch dir
        case 'down'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''globalcirc'',''down'')'));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''globalcirc'',''down'')'));
        case 'up'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''globalcirc'',''up'')'));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''globalcirc'',''up'')'));
      end
    case 'globalrad'
      switch dir
        case 'down'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''globalrad'',''down'')'));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''globalrad'',''down'')'));
        case 'up'
          myset(gui.handles.fig,'windowbuttonmotionfcn',sprintf('straintagging.strainrate(''setpeakmotion'',''globalrad'',''up'')'));
          myset(gui.handles.fig,'windowbuttonupfcn',sprintf('straintagging.strainrate(''setpeakbuttonup'',''globalrad'',''up'')'));
      end
  end


%-------------------------------------------------------
function setpeakmotion(type,dir,handleind)%type,ind)
%--------------------------------------------------
global DATA
gui = DATA.GUI.strainrate ;
h=gui.handles.segaxes;

switch type
  case 'circ'
    switch dir
      case 'up'
        peakhandle = gui.handles.peakscircup(handleind);
      case 'down'
        peakhandle = gui.handles.peakscircdown(handleind);
    end
    graphhandle=gui.handles.graphscirc(handleind);
  
  case 'rad'
    switch dir
      case 'up'
        peakhandle = gui.handles.peaksradup(handleind);
      case 'down'
        peakhandle = gui.handles.peaksraddown(handleind);
    end
    graphhandle=gui.handles.graphsrad(handleind);
  
  case 'globalcirc'
    switch dir
      case 'up'
        peakhandle = gui.handles.globalpeakcircup;
      case 'down'
        peakhandle = gui.handles.globalpeakcircdown;
    end
    graphhandle = gui.handles.globalgraphcirc;
  case 'globalrad'
    switch dir
      case 'up'
        peakhandle = gui.handles.globalpeakradup;
      case 'down'
        peakhandle = gui.handles.globalpeakraddown;
    end
    graphhandle = gui.handles.globalgraphrad;
end

[xclicked,~]=mygetcurrentpoint(h);
x=get(graphhandle,'Xdata');
[~,ind]=min((x-xclicked).^2);
y=get(graphhandle,'Ydata');
set(peakhandle,'Xdata',x(ind),'Ydata',y(ind))

if  get(gui.handles.typepopupmenu,'Value')>4 && ind>1
  ind=ind-1;
end

switch type
  case 'circ'
    switch dir
      case 'up'
        gui.circupind(handleind) = ind;
        gui.circup(handleind) = gui.srcirc(handleind,ind);%y(ind);
      case 'down'
        gui.circdownind(handleind) = ind;
        gui.circdown(handleind) = gui.srcirc(handleind,ind);%y(ind);
    end
  case 'rad'
    switch dir
      case 'up'
        gui.radupind(handleind) = ind;
        gui.radup(handleind) = gui.srrad(handleind,ind);%y(ind);
      case 'down'
        gui.raddownind(handleind) = ind;
        gui.raddown(handleind) = gui.srrad(handleind,ind);%y(ind);
    end
  case 'globalrad'
    switch dir
      case 'up'
        gui.globalradupind = ind;
        gui.globalradup = gui.globalsrrad(ind);
      case 'down'
        gui.globalraddownind = ind;
        gui.globalraddown = gui.globalsrrad(ind);
    end
  case 'globalcirc'
    switch dir
      case 'up'
        gui.globalcircupind = ind;
        gui.globalcircup = gui.globalsrcirc(ind);
      case 'down'
        gui.globalcircdownind = ind;
        gui.globalcircdown = gui.globalsrcirc(ind);
    end
end

for i=gui.views2show
  set(gui.handles.(['text',num2str(i)]),'String',sprintf('%0.1f',gui.circup(i)))
  set(gui.handles.(['text',num2str(17+i)]),'String',sprintf('%0.1f',gui.circdown(i)))
  set(gui.handles.(['text',num2str(2*17+i)]),'String',sprintf('%0.1f',gui.radup(i)))
  set(gui.handles.(['text',num2str(3*17+i)]),'String',sprintf('%0.1f',gui.raddown(i)))
end


set(gui.handles.('globalcircuptext'),'String',sprintf('%0.1f',gui.globalcircup))
set(gui.handles.('globalcircdowntext'),'String',sprintf('%0.1f',gui.globalcircdown))
set(gui.handles.('globalraduptext'),'String',sprintf('%0.1f',gui.globalradup))
set(gui.handles.('globalraddowntext'),'String',sprintf('%0.1f',gui.globalraddown))


function setpeakbuttonup(type,dir,handleind)%type,ind)
global DATA
gui = DATA.GUI.strainrate;

if nargin<3
setpeakmotion(type,dir)
else
setpeakmotion(type,dir,handleind)
end
set(gui.handles.fig,'windowbuttonmotionfcn',[]);
set(gui.handles.fig,'windowbuttonupfcn',[])

switch type
  case 'circ'
    myset(gui.handles.fig,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''circ'')'))
  case 'rad'
    myset(gui.handles.fig,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''rad'')'))
  case ''
end

function all_Callback
global DATA
gui = DATA.GUI.strainrate ;

if gui.oneatatime
  val=get(gui.handles.oneatatimecheckbox,'value');
  set(gui.handles.oneatatimecheckbox,'value',~val)
  oneatatime_Callback;
end


for i=gui.radios2show
  myset(gui.handles.(['radiobutton',num2str(i)]),'Value',1);
end

radio_Callback


%-----------------------------------------------
function oneatatime_Callback
%---------------------------------------------
global DATA
gui = DATA.GUI.strainrate;

oneatatime=get(gui.handles.oneatatimecheckbox,'value');
set(gui.handles.oneatatimecheckbox,'value',oneatatime);
gui.oneatatime=oneatatime;

if oneatatime
  ind=[];
  for i=1:17
    val=get(gui.handles.(['radiobutton',num2str(i)]),'Value');
    if val
      ind=i;
      break;
    end
  end
  
  if isempty(ind)
    set(gui.handles.('radiobutton1'),'Value',1)
    ind=1;
  end
  
  for i=1:17
    myset(gui.handles.(['radiobutton',num2str(i)]),'Value',0);
  end
  myset(gui.handles.(['radiobutton',num2str(ind)]),'Value',1);
  gui.lastpressed=ind;
  
  gui.views2show=ind;
  
%   %Auto detect global
% strlist=get(gui.handles.typepopupmenu,'String');
% ind=get(gui.handles.typepopupmenu,'Value');
% str=strlist{ind};

%if any(strcmp(str,{'Mean Longitudinal Strain Rate', 'Mean Radial Strain Rate', 'Mean Circumferential Strain Rate','Mean Longitudinal Strain with Tangents', 'Mean Radial Strain with Tangents', 'Mean Circumferential Strain with Tangents'}))
   setglobalpeaks;
%end

  updateplot;
end
%----------------------------------------------------------------
function radio_Callback
%---------------------------------------------------------------
global DATA
gui = DATA.GUI.strainrate ;

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

gui.views2show=find(tmp);

% %Auto detect global
% strlist=get(gui.handles.typepopupmenu,'String');
% ind=get(gui.handles.typepopupmenu,'Value');
% str=strlist{ind};

%if any(strcmp(str,{'Mean Longitudinal Strain Rate', 'Mean Radial Strain Rate', 'Mean Circumferential Strain Rate','Mean Longitudinal Strain with Tangents', 'Mean Radial Strain with Tangents', 'Mean Circumferential Strain with Tangents'}))
   setglobalpeaks;
%end

updateplot;

%---------------------------------------
function setglobalpeaks
%---------------------------------------
global SET DATA
gui = DATA.GUI.strainrate;
%Auto detect global

if strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')
    gui.globalsrrad = nanmean(gui.srrad(gui.views2show,:),1);
    gui.globalsrcirc = nanmean(gui.srcirc(gui.views2show,:),1);
else
  %global strainrate normalised if LAX images
  normstrainratecirc=[];
  normstrainraterad =[];

  for tagno = gui.taggroup
      normstrainratecirc = [normstrainratecirc, conv2([1 0 -1],1,SET(tagno).StrainTagging.globalcirc,'valid')/(2/(SET(tagno).TSize-1))];
      normstrainraterad = [normstrainraterad, conv2([1 0 -1],1,SET(tagno).StrainTagging.globalrad,'valid')/(2/(SET(tagno).TSize-1))];
  end

  gui.globalsrcirc=nanmean(normstrainratecirc(:,gui.views2show),2);
  gui.globalsrrad=nanmean(normstrainraterad(:,gui.views2show),2);
end

%change value but not time
gui.globalradup=gui.globalsrrad(gui.globalradupind);
gui.globalraddown=gui.globalsrrad(gui.globalraddownind);
gui.globalcircup=gui.globalsrcirc(gui.globalcircupind);
gui.globalcircdown=gui.globalsrcirc(gui.globalcircdownind);


%   [gui.globalradup, gui.globalradupind] = max(gui.globalsrrad);
%   [gui.globalraddown, gui.globalraddownind] = min(gui.globalsrrad);
%   [gui.globalcircup, gui.globalcircupind] = max(gui.globalsrcirc);
%   [gui.globalcircdown, gui.globalcircdownind] = min(gui.globalsrcirc);


%----------------------------------------
function plotNonSegmentedStrainrate(h,tvec,straincurve,type,normalizedtime)
%----------------------------------------
global DATA

if nargin<5
  normalizedtime=0;
end
gui=DATA.GUI.strainrate;

switch type
  case 'globalcirc'
    downslope=gui.globalcircdown;
    inddown=gui.globalcircdownind;
    upslope=gui.globalcircup;
    indup=gui.globalcircupind;
  case 'globalrad'
    downslope=gui.globalraddown;
    inddown=gui.globalraddownind;
    upslope=gui.globalradup;
    indup=gui.globalradupind;
end

plot(h,tvec,straincurve,'.-','linewidth',2,'Color','b')
hold(h,'on')
inddown = inddown+1;
indup=indup+1;
downtime = tvec(inddown);
uptime = tvec(indup);
plot(h,downtime,straincurve(inddown),'r+','markersize',10,'linewidth',3)
plot(h,uptime,straincurve(indup),'k+','markersize',10,'linewidth',3)

tangentplot(h,tvec,straincurve,normalizedtime,inddown,indup,downtime,uptime);

% xlim(h,[tvec(1),tvec(end)])
% ylim(h,'manual')
% xlim(h,'manual')
% x=-SET(gui.taggroup(1)).TSize:SET(gui.taggroup(1)).TSize;
% if ~normalizedtime
%   uptan_t=uptime+x*SET(gui.taggroup(1)).TIncr;
%   downtan_t=downtime+x*SET(gui.taggroup(1)).TIncr;
%   
%   uptan = upslope*SET(gui.taggroup(1)).TIncr*x+straincurve(indup);
%   downtan = downslope*SET(gui.taggroup(1)).TIncr*x+straincurve(inddown);
% else
%   uptan_t=uptime+x/(SET(gui.taggroup(1)).TSize-1);
%   downtan_t=downtime+x/(SET(gui.taggroup(1)).TSize-1);
%   
%   uptan = upslope/(SET(gui.taggroup(1)).TSize-1)*x+straincurve(indup);
%   downtan = downslope/(SET(gui.taggroup(1)).TSize-1)*x+straincurve(inddown);
% end
% 
% gui.handles.tup = plot(h,uptan_t,uptan,'k--');
% gui.handles.tdown = plot(h,downtan_t,downtan,'k--');

% if ~normalizedtime
%  m=straincurve(indup+1)-upslope*SET(gui.taggroup(1)).TIncr*indup;
% else
%   m=straincurve(indup+1)-upslope/(SET(gui.taggroup(1)).TSize-1)*indup;
% end
% 
% xup(1)=(yl(1)-m)/upslope;
% xup(2)=(yl(2)-m)/upslope;
% 
% if ~normalizedtime     
%      m=straincurve(inddown+1)-downslope*SET(gui.taggroup(1)).TIncr*inddown;
% else
%      m=straincurve(inddown+1)-downslope/(SET(gui.taggroup(1)).TSize-1)*inddown;  
% end
% 
% xdown(1)=(yl(1)-m)/downslope;
% % xdown(2)=(yl(2)-m)/downslope;
% gui.handles.tup = plot(h,xup,yl,'k--');
% gui.handles.tdown = plot(h,xdown,yl,'k--');
function tangentplot(h,tvec,straincurve,normalizedtime,inddown,indup,downtime,uptime,downslope,upslope)
global SET DATA

gui=DATA.GUI.strainrate;


if ~strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')
  if indup~=1
    upslope = -(straincurve(indup-1)-straincurve(indup+1))/(2/(SET(gui.taggroup(1)).TSize-1));
  else
    upslope = -(straincurve(indup)-straincurve(indup+1))/(1/(SET(gui.taggroup(1)).TSize-1));
  end
  if inddown~=1
    downslope = -(straincurve(inddown-1)-straincurve(inddown+1))/(2/(SET(gui.taggroup(1)).TSize-1));
  else
    downslope = -(straincurve(inddown)-straincurve(inddown+1))/(1/(SET(gui.taggroup(1)).TSize-1));
  end
else
  if indup~=1
    upslope = -(straincurve(indup-1)-straincurve(indup+1))/(2*SET(gui.taggroup(1)).TIncr);
  else
    upslope = -(straincurve(indup)-straincurve(indup+1))/(1*SET(gui.taggroup(1)).TIncr);
  end
  if inddown~=1
    downslope = -(straincurve(inddown-1)-straincurve(inddown+1))/(2*SET(gui.taggroup(1)).TIncr);
  else
    downslope = -(straincurve(inddown)-straincurve(inddown+1))/(1*SET(gui.taggroup(1)).TIncr);
  end
end

xlim(h,[tvec(1),tvec(end)])
ylim(h,'manual')
xlim(h,'manual')
x=-SET(gui.taggroup(1)).TSize:SET(gui.taggroup(1)).TSize;
if ~normalizedtime
  uptan_t=uptime+x*SET(gui.taggroup(1)).TIncr;
  downtan_t=downtime+x*SET(gui.taggroup(1)).TIncr;
  
  uptan = upslope*SET(gui.taggroup(1)).TIncr*x+straincurve(indup);
  downtan = downslope*SET(gui.taggroup(1)).TIncr*x+straincurve(inddown);
else
  uptan_t=uptime+x/(SET(gui.taggroup(1)).TSize-1);
  downtan_t=downtime+x/(SET(gui.taggroup(1)).TSize-1);
  
  uptan = upslope/(SET(gui.taggroup(1)).TSize-1)*x+straincurve(indup);
  downtan = downslope/(SET(gui.taggroup(1)).TSize-1)*x+straincurve(inddown);
end

gui.handles.tup = plot(h,uptan_t,uptan,'k--');
gui.handles.tdown = plot(h,downtan_t,downtan,'k--');

function updateplot
global DATA SET

gui = DATA.GUI.strainrate ;
myset(gui.handles.fig,'ButtonDownFcn',[])
h=gui.handles.segaxes;
legend(h,'off')

cla(h)
%cla(gui.handles.reportaxes)
%changed curve choice need to redo globalpeaks
% if gui.popupind~=get(gui.handles.typepopupmenu,'Value');
%   setglobalpeaks;
%   gui.popupind=get(gui.handles.typepopupmenu,'Value');
% end

set(h,'yticklabelmode','auto')
set(h,'ytickmode','auto')

ylim(h,'auto')
xlim(h,'auto')
% get names of segments

strlist=get(gui.handles.typepopupmenu,'String');
ind=get(gui.handles.typepopupmenu,'Value');
str=strlist{ind};

hold(h,'on')
grid(h,'on');

used=find(gui.include);

if  ~strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')
  normalizedtime=1;
  tvec=linspace(0,1,SET(gui.taggroup(1)).TSize);
  srtvec=tvec(2:end-1);
else
  srtvec=gui.tvec;
  tvec=SET(gui.taggroup(1)).TimeVector;
  normalizedtime=0;
end


switch str
   case 'Mean Longitudinal Strain with Tangents'
    straincurves=[];
    counter=1;
    
    for tagno=gui.taggroup
        straincurves = [straincurves, SET(tagno).StrainTagging.globalcirc];
    end
    
    straincurve=mynanmean(straincurves(:,gui.views2show),2);
    plotNonSegmentedStrainrate(h,tvec,straincurve,'globalcirc',normalizedtime);
    
  case 'Mean Circumferential Strain with Tangents'
    straincurve=mynanmean(SET(gui.taggroup(1)).StrainTagging.globalcirc(:,used(gui.views2show)),2);
    plotNonSegmentedStrainrate(h,SET(gui.taggroup(1)).TimeVector,straincurve,'globalcirc');
    
  case 'Mean Radial Strain with Tangents'
    straincurves=[];
    
    for tagno=gui.taggroup
        straincurves = [straincurves, SET(tagno).StrainTagging.globalrad];
    end
    if strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')
      straincurve=mynanmean(straincurves(:,used(gui.views2show)),2);
    else
      straincurve=mynanmean(straincurves(:,gui.views2show),2);  
    end
    plotNonSegmentedStrainrate(h,tvec,straincurve,'globalrad',normalizedtime);
    
  case {'Longitudinal Strain with Tangents','Circumferential Strain with Tangents'}
      straincurve = nan(SET(gui.taggroup(1)).TSize,length(gui.taggroup));
      used = find(gui.include);
     for i=gui.views2show
      if strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')
        straincurve(:,i)=SET(gui.taggroup(1)).StrainTagging.globalcirc(:,used(i));
      else
        straincurve(:,i)=SET(gui.taggroup(i)).StrainTagging.globalcirc;
%         switch SET(gui.taggroup(used(i))).ImageViewPlane
%           case '2CH'
%             straincurve(:,1)=SET(gui.taggroup(used(i))).StrainTagging.globalcirc;
%           case '3CH'
%             straincurve(:,end-1)=SET(gui.taggroup(i)).StrainTagging.globalcirc;
%           case '4CH'
%             straincurve(:,end)=SET(gui.taggroup(i)).StrainTagging.globalcirc;
%         end
      end
     end
    for i = gui.views2show
     gui.handles.graphscirc(i)=plot(h,tvec,straincurve(:,i),'.-','linewidth',2,'Color',gui.cmap(i,:));
     gui.handles.peakscircup(i)=plot(h,srtvec(gui.circupind(i)), straincurve(gui.circupind(i)+1,i),'k+','markersize',10,'linewidth',3);
     gui.handles.peakscircdown(i)=plot(h,srtvec(gui.circdownind(i)), straincurve(gui.circdownind(i)+1,i),'r+','markersize',10,'linewidth',3); 
      myset(gui.handles.peakscircup(i),'ButtonDownFcn','');
      myset(gui.handles.peakscircdown(i),'ButtonDownFcn','');
      myset(gui.handles.graphscirc(i),'ButtonDownFcn','');
    end
    
    
    for i = gui.views2show
      tangentplot(h,tvec,straincurve(:,i),normalizedtime,gui.circdownind(i)+1,gui.circupind(i)+1,...
        srtvec(gui.circdownind(i)),srtvec(gui.circupind(i)),gui.circdown(i),gui.circup(i));
    end

    myset(h,'ButtonDownFcn',''); 

  case'Radial Strain with Tangents'
     straincurve = nan(SET(gui.taggroup(1)).TSize,length(gui.views2show));
     used = find(gui.include);
     for i=1:length(used)
       if strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')
         straincurve(:,i)=SET(gui.taggroup(1)).StrainTagging.globalrad(:,used(i));
       else
         straincurve(:,i)=SET(gui.taggroup(i)).StrainTagging.globalrad;
%          switch SET(gui.taggroup(i)).ImageViewPlane
%            case '2CH'
%              straincurve(:,1)=SET(gui.taggroup(i)).StrainTagging.globalrad;
%            case '3CH'
%              straincurve(:,end-1)=SET(gui.taggroup(i)).StrainTagging.globalrad;
%            case '4CH'
%              straincurve(:,end)=SET(gui.taggroup(i)).StrainTagging.globalrad;
%          end
       end
     end
     for i = gui.views2show
     gui.handles.graphsrad(i)=plot(h,tvec,straincurve(:,i),'.-','linewidth',2,'Color',gui.cmap(i,:));
     gui.handles.peaksradup(i)=plot(h,srtvec(gui.radupind(i)), straincurve(gui.radupind(i)+1,i),'k+','markersize',10,'linewidth',3);
     gui.handles.peaksraddown(i)=plot(h,srtvec(gui.raddownind(i)), straincurve(gui.raddownind(i)+1,i),'r+','markersize',10,'linewidth',3); 
     
      myset(gui.handles.peakscircup(i),'ButtonDownFcn','');
      myset(gui.handles.peakscircdown(i),'ButtonDownFcn','');
      myset(gui.handles.graphscirc(i),'ButtonDownFcn','');
    end
    
    for i = gui.views2show
      tangentplot(h,tvec,straincurve(:,i),normalizedtime,gui.raddownind(i)+1,gui.radupind(i)+1,...
        srtvec(gui.raddownind(i)),srtvec(gui.radupind(i)),gui.raddown(i),gui.radup(i));
    end
    myset(h,'ButtonDownFcn',''); 

  
  case {'Circumferential Strain Rate', 'Longitudinal Strain Rate'}
    for i = gui.views2show
     gui.handles.graphscirc(i)=plot(h,srtvec,gui.srcirc(i,:),'.-','linewidth',2,'Color',gui.cmap(i,:));
     gui.handles.peakscircup(i)=plot(h,srtvec(gui.circupind(i)), gui.circup(i),'k+','markersize',10,'linewidth',3);
     gui.handles.peakscircdown(i)=plot(h,srtvec(gui.circdownind(i)), gui.circdown(i),'r+','markersize',10,'linewidth',3); 
     myset(gui.handles.peakscircup(i),'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''circ'',''up'')'));
     myset(gui.handles.peakscircdown(i),'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''circ'',''down'')'));
     myset(gui.handles.graphscirc(i),'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''circ'')'));
    end
     yl=ylim;
     xl=[srtvec(1),srtvec(end)];
     set(h,'XLim',xl,'YLim',yl)
    myset(h,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''circ'')'));
   
     
  case 'Radial Strain Rate'
    for i = gui.views2show
     gui.handles.graphsrad(i)=plot(h,srtvec,gui.srrad(i,:),'.-','linewidth',2,'Color',gui.cmap(i,:));
     gui.handles.peaksradup(i)=plot(h,srtvec(gui.radupind(i)), gui.radup(i),'k+','markersize',10,'linewidth',3);
     gui.handles.peaksraddown(i)=plot(h,srtvec(gui.raddownind(i)), gui.raddown(i),'r+','markersize',10,'linewidth',3); 
     myset(gui.handles.peaksradup(i),'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''rad'',''up'')'));
     myset(gui.handles.peaksraddown(i),'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''rad'',''down'')'));
     myset(gui.handles.graphsrad(i),'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''rad'')'));
    end
     yl=ylim;
     xl=[srtvec(1),srtvec(end)];
     set(h,'XLim',xl,'YLim',yl)
     myset(h,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''rad'')'));
     
   case {'Mean Circumferential Strain Rate', 'Mean Longitudinal Strain Rate'}
     gui.handles.globalgraphcirc=plot(h,srtvec,gui.globalsrcirc,'.-','linewidth',2,'Color','b');
     gui.handles.globalpeakcircup=plot(h,srtvec(gui.globalcircupind), gui.globalcircup,'k+','markersize',10,'linewidth',3);
     gui.handles.globalpeakcircdown=plot(h,srtvec(gui.globalcircdownind), gui.globalcircdown,'r+','markersize',10,'linewidth',3); 
     
     myset(gui.handles.globalpeakcircup,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalcirc'',''up'')'));
     myset(gui.handles.globalpeakcircdown,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalcirc'',''down'')'));
     myset(gui.handles.globalgraphcirc,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalcirc'')'));
    
     yl=ylim;
     xl=[srtvec(1),srtvec(end)];
     set(h,'XLim',xl,'YLim',yl)
    myset(h,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalcirc'')'));
   
     
  case 'Mean Radial Strain Rate'
     gui.handles.globalgraphrad=plot(h,srtvec,gui.globalsrrad,'.-','linewidth',2,'Color','b');
     gui.handles.globalpeakradup=plot(h,srtvec(gui.globalradupind), gui.globalradup,'k+','markersize',10,'linewidth',3);
     gui.handles.globalpeakraddown=plot(h,srtvec(gui.globalraddownind), gui.globalraddown,'r+','markersize',10,'linewidth',3); 
     
     myset(gui.handles.globalpeakradup,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalrad'',''up'')'));
     myset(gui.handles.globalpeakraddown,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalrad'',''down'')'));
     myset(gui.handles.globalgraphrad,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalrad'')'));
    
     yl=ylim;
     xl=[srtvec(1),srtvec(end)];
     set(h,'XLim',xl,'YLim',yl)
    myset(h,'ButtonDownFcn',sprintf('straintagging.strainrate(''setpeak_Callback'',''globalrad'')'));
   
end

if  strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')%~any(strcmp(str,{'Mean Longitudinal Strain with Tangents','Mean Radial Strain with Tangents','Mean Longitudinal Strain Rate','Mean Radial Strain Rate','Mean Longitudinal Strain Rate','Radial Strain Rate'})) || strcmp(SET(gui.taggroup(1)).ImageViewPlane,'Short-axis')
  xlabel(h,'Time [s]')
else
  xlabel(h,'Heart cycle [s/T]')
end

if ind >4
  ylabel(h,'%/s')
else
  ylabel(h,'%')
end
hold(h,'off')

for i=gui.views2show
  set(gui.handles.(['text',num2str(i)]),'String',sprintf('%0.1f',gui.circup(i)))
  set(gui.handles.(['text',num2str(17+i)]),'String',sprintf('%0.1f',gui.circdown(i)))
  set(gui.handles.(['text',num2str(2*17+i)]),'String',sprintf('%0.1f',gui.radup(i)))
  set(gui.handles.(['text',num2str(3*17+i)]),'String',sprintf('%0.1f',gui.raddown(i)))
end

set(gui.handles.('globalcircuptext'),'String',sprintf('%0.1f',gui.globalcircup))
set(gui.handles.('globalcircdowntext'),'String',sprintf('%0.1f',gui.globalcircdown))
set(gui.handles.('globalraduptext'),'String',sprintf('%0.1f',gui.globalradup))
set(gui.handles.('globalraddowntext'),'String',sprintf('%0.1f',gui.globalraddown))

function export_Callback
global DATA SET
gui = DATA.GUI.strainrate ;

col=1;
line=1;
outdata{line,1} = 'Patient Name';
outdata{line,2} = 'Patient ID';
outdata{line,3} = 'Heart Rate';
outdata{line,4} = 'Image Type';

line=2;
outdata{line,1} = SET(no).PatientInfo.Name;
outdata{line,2} = SET(no).PatientInfo.ID;
outdata{line,3} = SET(no).HeartRate;
outdata{line,4} = sprintf('%s %s',SET(no).ImageType,SET(no).ImageViewPlane);
% 
% %Write relevant gui parameters to outdatacell
% out{line,5} = gui.radup;
% sr.radupind = gui.radupind;
% sr.raddown = gui.raddown;
% sr.raddownind = gui.raddownind;
% sr.circup = gui.circup;
% sr.circupind=gui.circupind;
% sr.circdown=gui.circdown;
% sr.circdownind=gui.circdownind;
% sr.srcirc = gui.srcirc;
% sr.srrad = gui.srrad;
% sr.include = gui.include;
% sr.globalcircup=gui.globalcircup;
% sr.globalcircdown=gui.globalcircdown;
% sr.globalcircupind=gui.globalcircupind;
% sr.globalcircdownind=gui.globalcircdownind;
% sr.globalradup=gui.globalradup;
% sr.globalraddown=gui.globalraddown;
% sr.globalradupind=gui.globalradupind;
% sr.globalraddownind=gui.globalraddownind;
% sr.globalsrcirc=gui.globalsrcirc;
% sr.globalsrrad=gui.globalsrrad;

function save_Callback
global DATA SET
gui = DATA.GUI.strainrate ;

%Write relevant gui parameters to SET structure
sr.radup = gui.radup;
sr.radupind = gui.radupind;
sr.raddown = gui.raddown;
sr.raddownind = gui.raddownind;
sr.circup = gui.circup;
sr.circupind=gui.circupind;
sr.circdown=gui.circdown;
sr.circdownind=gui.circdownind;
sr.srcirc = gui.srcirc;
sr.srrad = gui.srrad;
sr.include = gui.include;
sr.globalcircup=gui.globalcircup;
sr.globalcircdown=gui.globalcircdown;
sr.globalcircupind=gui.globalcircupind;
sr.globalcircdownind=gui.globalcircdownind;
sr.globalradup=gui.globalradup;
sr.globalraddown=gui.globalraddown;
sr.globalradupind=gui.globalradupind;
sr.globalraddownind=gui.globalraddownind;
sr.globalsrcirc=gui.globalsrcirc;
sr.globalsrrad=gui.globalsrrad;

for no = gui.taggroup
  SET(no).StrainTagging.SR=sr;
end

close_Callback

function close_Callback
global DATA

try
  DATA.GUI.strainrate = close(DATA.GUI.strainrate );
catch   %#ok<CTCH>
  DATA.GUI.strainrate =[];
  delete(gcbf);
end
