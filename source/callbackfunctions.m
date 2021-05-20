function callbackfunctions(varargin)
% This .m file has the ambition of containing all Callback functions

% Split out by Klas

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
feval(varargin{:}); % FEVAL switchyard

%---------------------
function orthoview_Callback %#ok<DEFNU>
%-----------------------
global SET NO
%First we check so that there are any slices because if not there is no
%point

if SET(NO).ZSize==1
  myfailed('To few slices for orthogonal view.')
  return
end

%this sets all necessary fields
viewfunctions('orthoview',NO);

%This updates the view.
viewfunctions('setview',2,2,ones(1,4)*NO,{'orth','hla','vla','gla'})

%------------------------------
function viewviewport_Callback %#ok<DEFNU>
%------------------------------
%Displays viewport in one panel

global DATA

segment3dp.tools('storetoobject');

viewfunctions('setview',1,1,[],{'speedim'})

DATA.drawlist = [];
DATA.drawlist{1} = {'segment3dp.tools(''helprender'');'};

drawfunctions('setxynan',1);

run = true;
DATA.Handles.configiconholder.indent('view3d',run);

segment3dp.tools('helprender');

DATA.LevelSet.ViewPort.panelhandle.Position = DATA.Handles.imageaxes(1).Position;

%------------------------------
function view4panel3dp_Callback %#ok<DEFNU>
%------------------------------
global DATA

%indent(DATA.Handles.configiconholder,'view4',1);

segment3dp.tools('storetoobject');

viewfunctions('setview',2,2,[],{'trans3DP','speedim','sag3DP','cor3DP'})

if DATA.Handles.configiconholder.findindented('view3d')
  %Update the size
  DATA.LevelSet.ViewPort.panelhandle.Position = DATA.Handles.imageaxes(2).Position;
  DATA.drawlist{2} = {'segment3dp.tools(''helprender'');'};
end

%----------------------------
function selectpanel_Callback(name) %#ok<DEFNU>
%----------------------------
%Select the panel with name

global DATA SET

selpanel = [];
for loop = 1:length(DATA.ViewPanelsType)
  if isequal(DATA.ViewPanelsType{loop},name)
    selpanel = loop;
  end
end

%if ~isempty(selpanel)
%  segment('switchtopanel',selpanel)
%end

if ~isempty(selpanel)
  DATA.CurrentPanel = selpanel;
  no = DATA.ViewPanels(selpanel);
  switch name
    case 'trans3DP'
      SET(no).LevelSet.Pen.Color = 'r';
    case 'sag3DP'
      SET(no).LevelSet.Pen.Color = 'g';
    case 'cor3DP'
      SET(no).LevelSet.Pen.Color = 'b';
  end
  
  for loop = 1:length(DATA.ViewPanelsType)
    if isequal(DATA.ViewPanelsType{loop},name)
      DATA.CurrentPanel = loop;
    end
  end
  
  %Graphical update
  drawfunctions('drawselectedframe',selpanel)
  speedimpanel = find(strcmp(DATA.ViewPanelsType,'speedim'));
  if isempty(speedimpanel)
    return
  end
  drawfunctions('drawpanel',speedimpanel);
  viewfunctions('updatezoomandaspectratio',speedimpanel);
end

%------------------------
function viewhelper(type)
%------------------------
%Set a two panel view with anatomical and speed image
global DATA SET NO

segment3dp.tools('storetoobject');

speedpanel = 2; %currently is always panel number2

switch type
  case 'trans3DP'
    SET(NO).LevelSet.Pen.Color = 'r';
  case 'sag3DP'
    SET(NO).LevelSet.Pen.Color = 'g';
  case 'cor3DP'
    SET(NO).LevelSet.Pen.Color = 'b';
end

if DATA.Handles.configiconholder.findindented('view3d')
  
  %Clear all graphics function from drawlist
  drawfunctions('setxynan',2)
  viewfunctions('setview',1,2,[],{type,'viewport'})
  
  %Update the size
  DATA.LevelSet.ViewPort.panelhandle.Position = DATA.Handles.imageaxes(2).Position;
  
  DATA.drawlist{speedpanel} = {'segment3dp.tools(''helprender'');'};
  
else
  viewfunctions('setview',1,2,[],{type,'speedim'})
  %viewfunctions('updatezoomandaspectratio',1);
  %viewfunctions('updatezoomandaspectratio',2);
end

%----------------------------
function viewtransversal_Callback %#ok<DEFNU>
%----------------------------
viewhelper('trans3DP');

%----------------------------
function viewcoronal_Callback %#ok<DEFNU>
%----------------------------
viewhelper('cor3DP');

%----------------------------
function viewsagittal_Callback %#ok<DEFNU>
%----------------------------
viewhelper('sag3DP');

%----------------------------
function reset3d_Callback %#ok<DEFNU>
%----------------------------
%Reset 3D view

global DATA

if segment3dp.isviewportalive
  %3D display on rotate
  volshowobject = DATA.LevelSet.ViewPort.getvolshowobject;
  
  %Reset to
  newpos = [0 -4 0];
  
  %Set the new camera position
  ax = DATA.LevelSet.ViewPort.getaxeshandle;
  volshowobject.CameraPosition = newpos;
  volshowobject.CameraUpVector = [0 0 -1];
  ax.CameraPosition = newpos;
  ax.CameraUpVector = [0 0 -1];
end

%------------------------------------------
function set3dvieworientation_Callback(view) %#ok<DEFNU>
%------------------------------------------
%Set predefined views

global DATA

newposition = [];
newvector = [];

%Check if viewport is active
if segment3dp.isviewportalive
  
  d = 4;
  
  switch view
    case 'transversal'
      newposition = [0 0 -d];
      newvector = [0 -1 0];
    case 'sagittal'
      newposition = [-d 0 0];
      newvector = [0 0 -1];
    case 'coronal'
      newposition = [0 -d 0];
      newvector = [0 0 -1];
  end
  
  %Make the update
  if ~isempty(newposition)
    DATA.LevelSet.ViewPort.setcameraposition(newposition);
  end
  
  if ~isempty(newvector)
    DATA.LevelSet.ViewPort.setcameraupvector(newvector);
  end
  
  DATA.LevelSet.ViewPort.updatepoints;
end

%-----------------------------------
function set3dviewotherside_Callback %#ok<DEFNU>
%-----------------------------------
%Set predefined views

global DATA

%Check if viewport is active
if segment3dp.isviewportalive
  p = DATA.LevelSet.ViewPort.getcameraposition;
  DATA.LevelSet.ViewPort.setcameraposition(-p);
  DATA.LevelSet.ViewPort.updatepoints;
end

%----------------------------
function zoomin3d_Callback(f) %#ok<DEFNU>
%----------------------------
%Zoom in viewpanel

global DATA

if segment3dp.isviewportalive
  pos = DATA.LevelSet.ViewPort.getvolshowobject.CameraPosition;
  DATA.LevelSet.ViewPort.setcameraposition(pos*(1+f));
  DATA.LevelSet.ViewPort.updatepoints;
end

%--------------------------------------
function rotate3d_Callback(dtheta,dphi) %#ok<DEFNU>
%--------------------------------------
%Rotate 3D view with dtheta, and dphi degrees
%Obs this function is not defined from RLAPFH coordinate system

global DATA

if segment3dp.isviewportalive
  
  %3D display on rotate
  volshowobject = DATA.LevelSet.ViewPort.getvolshowobject;
  pos = volshowobject.CameraPosition;
  
  r = sqrt(sum(pos.*pos));
  theta = acos(pos(3)/r);
  phi = atan2(pos(2),pos(1));
  
  theta = theta+dtheta/180*pi;
  phi = phi+dphi/180*pi;
  
  newpos = [...
    r*sin(theta)*cos(phi) ...
    r*sin(theta)*sin(phi) ...
    r*cos(theta)];
    
  %Set the new camera position
  ax = DATA.LevelSet.ViewPort.getaxeshandle;
  volshowobject.CameraPosition = newpos;
  ax.CameraPosition = newpos;
  
end

%----------------------
function lasso_Callback %#ok<DEFNU>
%----------------------
%Lasso placeholder

myfailed('Lasso not yet available.');

%---------------------
function mmode_Callback %#ok<DEFNU>
%-----------------------
myfailed('This functions is no longer available.')

%---------------------
function play_Callback(panels) %#ok<DEFNU>
%-----------------------
global DATA SET

if nargin == 0
  panels = find(DATA.ViewPanels>0);
end

maxt = 0;
for loop=DATA.ViewPanels(panels)
  if SET(loop).TSize>maxt
    maxt = SET(loop).TSize;
  end
end

if maxt == 1
  myfailed('Need a time resolved image stack')
  undent(DATA.Handles.permanenticonholder,'play',0)
  return
end

%Try different approach where all stacks are played after NO
startframes = SET(DATA.ViewPanels(panels(1))).CurrentTimeFrame;
starttime = now;
beattime = SET(DATA.ViewPanels(panels(1))).BeatTime;

%prior to running the following graphics objects are turned 'off' by
%setting xdata and ydata to nan. This is because these objects need to be
%many in order to e.g draw text at multiple locations.
for p = panels
  % % %   for c = 'cygmkwrb'
  % % %     set(DATA.Handles.([c,'roitext'])(p,:),'Position',[nan,nan])
  % % %   end
  set(DATA.Handles.roitext(p,:),'Position',[nan,nan])
  set(DATA.Handles.measurementtext(p,:),'Position',[nan,nan])
  set(DATA.Handles.pointtext(p,:),'Position',[nan,nan])
  
  %interps not necessary when playing
  type = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
  for t = 1:length(type)
    DATA.Handles.(lower(type{t}))(p).XData = nan;
    DATA.Handles.(lower(type{t}))(p).YData = nan;
  end
end

synchronize = findindented(DATA.Handles.hideiconholder,'synchronize');

DATA.Run=1;
t = maxt;
while DATA.Run%true
  %stopping criteria
  if ~findindented(DATA.Handles.permanenticonholder,'play')
    DATA.Run=0;
  end
  
  for p = panels
    no = DATA.ViewPanels(p);
    if SET(no).TSize>1
      if DATA.Record
        t = 1+mod(t,maxt);
      else
        t = 1+mod(floor(rem(now-starttime,1)*24*3600/(beattime/maxt)+startframes),maxt);
      end
      SET(DATA.ViewPanels(p)).CurrentTimeFrame = max(min(round(SET(no).TSize*(t/maxt)),SET(no).TSize),1);
      drawfunctions('drawpanel',p);
    end
  end
  
  %also update the timebars!
  viewfunctions('updatetimebars',synchronize)
  drawnow %limitrate 
  if DATA.Record  %For Automatic Record Movie functionality
    %drawnow;
    gui = DATA.GUI.Segment;
    DATA.MovieFrame = mygetframe(gui.fig);
    DATA.MovieFrame.colormap=colormap(gui.fig);
    %DATA.MovieFrame = mygetframe(gui.handles.imageaxes);
    export('exportmovierecorder_Callback','newframe');
  end
end

%---------------------
function play_Callback2(panels) %#ok<DEFNU>
%-----------------------
global DATA SET NO

if nargin == 0
  %panels = 1:length(DATA.ViewPanels);
  panels = find(DATA.ViewPanels>0);%1:length(DATA.ViewPanels);
end

maxt = 0;
for loop=DATA.ViewPanels(panels)
  if SET(loop).TSize>maxt
    maxt = SET(loop).TSize;
  end
end

if maxt == 1
  myfailed('Need a time resolved image stack')
  undent(DATA.Handles.permanenticonholder,'play',0)
  return
end

%Try different approach where all stacks are played after NO
%startframes = SET(DATA.ViewPanels(panels(1))).CurrentTimeFrame;
%starttime = now;
%beattime = SET(DATA.ViewPanels(panels(1))).BeatTime;

%prior to running the following graphics objects are turned 'off' by
%setting xdata and ydata to nan. This is because these objects need to be
%many in order to e.g draw text at multiple locations.
for p = panels
  set(DATA.Handles.roitext(p,:),'Position',[nan,nan])
  set(DATA.Handles.measurementtext(p,:),'Position',[nan,nan])
  set(DATA.Handles.pointtext(p,:),'Position',[nan,nan])
  set(DATA.Handles.pointtext(p),'Position',[nan,nan])
  
  %interps not necessary when playing
  type = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
  for t = 1:length(type)
    DATA.Handles.(lower(type{t}))(p).XData = nan;
    DATA.Handles.(lower(type{t}))(p).YData = nan;
  end
end

DATA.Run=1;
while DATA.Run%true
  %stopping criteria
  tic
  if ~findindented(DATA.Handles.permanenticonholder,'play')
    DATA.Run=0;
  end
  
  %Update current timeframe in NO
  if SET(NO).CurrentTimeFrame < SET(NO).TSize
    SET(NO).CurrentTimeFrame = SET(NO).CurrentTimeFrame+1;
  else
    SET(NO).CurrentTimeFrame = 1;
  end
  
  for p = panels
    no = DATA.ViewPanels(p);
    if SET(no).TSize>1
      %t = 1+mod(floor(rem(now-starttime,1)*24*3600/(beattime/maxt)+startframes),maxt);
      %SET(DATA.ViewPanels(p)).CurrentTimeFrame = max(min(round(SET(no).TSize*(t/maxt)),SET(no).TSize),1);
      
      %find the closest matching phase to the currenttimeframe in NO
      [~,t]= min((SET(NO).TimeVector(SET(NO).CurrentTimeFrame)/SET(NO).TimeVector(end)-SET(no).TimeVector/SET(no).TimeVector(end)).^2);
      
      %if t~=SET(no).CurrentTimeFrame || no ~=NO
      SET(no).CurrentTimeFrame = t;
      drawfunctions('drawpanel',p);
      %end
    end
  end
  
  %also update the timebars!
  viewfunctions('updatetimebars')
  drawnow
  
  %then we get the spent time rendering and pause for the residual
  t_time = toc;
  if t_time<SET(NO).TIncr/2
    pause(SET(NO).TIncr/2-t_time)%SET(NO).BeatTime/(SET(NO).TSize-1)-t_time)
  end
  
end

%-------------------------------------------
function segmentclearalllv_Callback %#ok<DEFNU>
%-------------------------------------------
%Clear all LV segmentation, both endo and epi

global DATA

if nargin==1
  if ~yesno('Do you really want to remove all LV segmentation ?',[],DATA.GUI.Segment)
    %     myfailed('Aborted by user.',DATA.GUI.Segment);
    return;
  end
end

segmentation('clearalllv_Callback');

%-------------------------------------------
function segmentclearallrv_Callback %#ok<DEFNU>
%-------------------------------------------
%Clear all RV segmentation, both endo and epi

global DATA

if nargin==1
  if ~yesno('Do you really want to remove all RV segmentation ?',[],DATA.GUI.Segment)
    %     myfailed('Aborted by user.',DATA.GUI.Segment);
    return;
  end
end

segmentation('clearallrv_Callback');

%-------------------------------------------------------
function segmentclearalllvbutsystolediastole_Callback %#ok<DEFNU>
%-------------------------------------------------------
%Clears all LV segmentation except in systole and diastole.

global SET NO DATA

if SET(NO).TSize<2
  myfailed('Not timeresolved data, aborting.',DATA.GUI.Segment);
  return;
end

if isequal(SET(NO).EDT,SET(NO).EST)
  myfailed('Systole and diastole occurs at the same time frame, aborting.',DATA.GUI.Segment);
  return;
end

%Create index structure
ind = true(1,SET(NO).TSize);
ind(SET(NO).EDT) = false;
ind(SET(NO).EST) = false;
arg = struct('endo',true,'epi',true,'rvendo',false,'rvepi',false);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind);
segmentation('removeallinterp_Callback',true,[],arg,indarg);

%Create index structure
ind = true(1,SET(NO).TSize);
ind(SET(NO).EDT) = false;
ind(SET(NO).EST) = false;

if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX(:,ind,:) = NaN;
  SET(NO).EndoY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX(:,ind,:) = NaN;
  SET(NO).EpiY(:,ind,:) = NaN;
end

SET(NO).EndoDraged(ind,:) = false;
SET(NO).EpiDraged(ind,:) = false;

lvsegchanged = true; segment('updatevolume',lvsegchanged);
drawfunctions('drawcontours', DATA.CurrentPanel)
drawfunctions('drawinterp', DATA.CurrentPanel)


%-------------------------------------------------------
function segmentclearallrvbutsystolediastole_Callback %#ok<DEFNU>
%-------------------------------------------------------
%Clears all RV segmentation except in systole and diastole.

global SET NO DATA

if SET(NO).TSize<2
  myfailed('Not timeresolved data, aborting.',DATA.GUI.Segment);
  return;
end

if isequal(SET(NO).EDT,SET(NO).EST)
  myfailed('Systole and diastole occurs at the same time frame, aborting.',DATA.GUI.Segment);
  return;
end

%Create index structure
ind = true(1,SET(NO).TSize);
ind(SET(NO).EDT) = false;
ind(SET(NO).EST) = false;
arg = struct('endo',false,'epi',false,'rvendo',true,'rvepi',true);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind);
segmentation('removeallinterp_Callback',true,[],arg,indarg);

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX(:,ind,:) = NaN;
  SET(NO).RVEndoY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX(:,ind,:) = NaN;
  SET(NO).RVEpiY(:,ind,:) = NaN;
end

SET(NO).EndoDraged(ind,:) = false;
SET(NO).EpiDraged(ind,:) = false;

segment('updatevolume');
drawfunctions('drawcontours',DATA.CurrentPanel)
drawfunctions('drawinterp',DATA.CurrentPanel)

%-----------------------------------------
function contextloadno_Callback(viewpanelstype) %#ok<DEFNU>
%-----------------------------------------
%Loads imagestack no with the desired viewmode
global DATA

no = DATA.LastObject;
%Already out then switch to it and change viewmode
if any(no==DATA.ViewPanels)
  viewfunctions('switchimagestack',no,viewpanelstype)
else %add to first
  viewfunctions('addno2panel',1,no,viewpanelstype)
end

DATA.LastObject = [];


%-----------------------------------------
function pointdeletethis_Callback %#ok<DEFNU>
%-----------------------------------------
%remove point closest to click
global DATA SET
no = DATA.ViewPanels(DATA.CurrentPanel);
pointind = DATA.LastObject(1);
slice = DATA.LastObject(2);
tf = DATA.LastObject(3);

SET(no).Point.X(pointind) = [];
SET(no).Point.Y(pointind) = [];
SET(no).Point.T(pointind) = [];
SET(no).Point.Z(pointind) = [];
SET(no).Point.Label(pointind) = [];

drawfunctions('drawno',no);

%-----------------------------------------
function pointrenamethis_Callback %#ok<DEFNU>
%-----------------------------------------
%rename point closest to click
global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if noobjectshelper
  return
end

pointind = DATA.LastObject(1);
if strcmp(DATA.ProgramName, 'Segment 3DPrint')
  s = myinputdlg({'Enter name'},'Name',1,{sprintf('%s',SET(no).Point.Label{pointind})});
  if isempty(s) || isempty(s{1})
    return;
  else
    SET(no).Point.Label{pointind} = s{1};
  end
else
  % in other programs show a pre-defined set of names
  menuitems = {...
    'Apex',...
    'RV insertion',...
    'RV insertion Anterior',...
    'RV insertion Inferior',...
    'AV plane',...
    'TV plane',...
    'P1',...
    'P2',...
    'General',...
    'Sector start',...
    'User defined ...'};
  
  % s = myinputdlg({'Enter name'},'Name',1,{sprintf('%s',SET(no).Point.Label{pointind})});
  
  s = mymenu(strcat(dprintf('Select a new name for the point')),menuitems,DATA.GUI.Segment);
  if isempty(s)
    myfailed('Invalid name.',DATA.GUI.Segment);
    return;
  elseif s == 0 % cancel was clicked
    return
  elseif s == length(menuitems)
    % user defined name
    name = myinputdlg({'Enter name'},'Name',1,{sprintf('%s',SET(no).Point.Label{pointind})});
    if isempty(name) || isempty(name{1})
      return;
    else
      SET(no).Point.Label{pointind} = name{1};
    end
  else
    SET(no).Point.Label{pointind} = menuitems{s};
  end
end

drawfunctions('drawno',no);

%-----------------------------------------
function pointclearall_Callback %#ok<DEFNU>
%-----------------------------------------
%remove all points

global DATA SET NO

no = NO;%DATA.ViewPanels(DATA.CurrentPanel);

SET(no).Point.X = [];
SET(no).Point.Y = [];
SET(no).Point.T = [];
SET(no).Point.Z = [];
SET(no).Point.Label = {};

drawfunctions('drawno',no);
if any(strcmp(DATA.ProgramName,{'Segment 3DPrint'}))
  if segment3dp.isviewportalive
    DATA.LevelSet.ViewPort.setpoints(NaN,NaN,NaN)
  end
end

%-----------------------------------------
function pointmaketimeresolvedthis_Callback %#ok<DEFNU>
%-----------------------------------------
%make point closest to click timeresolved
global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if noobjectshelper
  return
end

pointind = DATA.LastObject(1);

%nan means timeresolved.
SET(no).Point.T(pointind) = nan;

drawfunctions('drawno',no);

%-----------------------------------------
function pointshowthisframeonly_Callback %#ok<DEFNU>
%-----------------------------------------
%remove point closest to click
global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if noobjectshelper
  return
end

pointind = DATA.LastObject(1);
tf = DATA.LastObject(3);

SET(no).Point.T(pointind) = tf;

drawfunctions('drawno',no);

%-----------------------------------------
function interpdeletepoint %#ok<DEFNU>
%-----------------------------------------
%remove interpolation point closest to click

global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if noobjectshelper
  return
end

pointind = DATA.LastObject(1);
slice = DATA.LastObject(2);
tf = DATA.LastObject(3);
pointtype = DATA.LastObjectType;

SET(no).([pointtype,'X']){tf,slice}(pointind) = [];
SET(no).([pointtype,'Y']){tf,slice}(pointind) = [];

%do resampling of curve when deleting single point.
X = SET(no).([pointtype,'X']){tf,slice};
Y = SET(no).([pointtype,'Y']){tf,slice};
%removes duplicate points and resamples the contour
[x,y] = calcfunctions('resamplecurve',X,Y,DATA.NumPoints-1);
SET(no).([pointtype(1:end-6),'Y'])(:,tf,slice)=[y,y(1)];
SET(no).([pointtype(1:end-6),'X'])(:,tf,slice)=[x,x(1)];

drawfunctions('drawno',no);

%-----------------------------------------
function interpdeletepointthisslicephase %#ok<DEFNU>
%-----------------------------------------
%remove interpolation point closest to click
global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if noobjectshelper
  return
end

pointind = DATA.LastObject(1);
slice = DATA.LastObject(2);
tf = DATA.LastObject(3);
pointtype = DATA.LastObjectType;

SET(no).([pointtype,'X']){tf,slice} = [];
SET(no).([pointtype,'Y']){tf,slice} = [];

drawfunctions('drawno',no);

%-------------------------------------
function measureclearthis_Callback(~)  %#ok<DEFNU>
%-------------------------------------
%Clear this measurement

global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if nargin==0
  if noobjectshelper
    return
  end
  measureind = DATA.LastObject(1);
else
  %Ask
  measureind = measureaskhelper;
end

if isequal(measureind,0)
  %   myfailed('Aborted.');
  return;
end

SET(no).Measure(measureind) = [];
drawfunctions('drawno',no);
segment('updatemeasurement');

%--------------------------------
function measureclearall_Callback %#ok<DEFNU>
%--------------------------------
%Clear all measurements

global DATA SET
no = DATA.ViewPanels(DATA.CurrentPanel);

SET(no).Measure = [];
drawfunctions('drawno',no);
segment('updatemeasurement');

%----------------------------
function m = measureaskhelper
%----------------------------
%Helper function to select a measurement in a menu

global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

m = 0;
if isempty(SET(no).Measure)
  myfailed('No measurements.');
  return;
end

if length(SET(no).Measure)==1
  m = 1;
  return;
end

%Loop over names
namecell = cell(1,length(SET(no).Measure));
for loop = 1:length(namecell)
  namecell{loop} = SET(no).Measure(loop).LongName;
  if isempty(namecell{loop})
    namecell{loop} = 'Empty name';
  end
end

%ask
m = mymenu('Select measurement',namecell);

%-------------------------------------
function measurerenamethis_Callback(~) %#ok<DEFNU>
%-------------------------------------
%Rename measurement

global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if nargin==0
  if noobjectshelper
    return
  end
  measureind = DATA.LastObject(1);
else
  %ask
  measureind = measureaskhelper;
end

if isequal(measureind,0)
  %   myfailed('Aborted.');
  return;
end

tools('enableundo',no);

%Get new name
[stri,lstr] = measureasklabel(measureind);
if ~isempty(stri)
  SET(no).Measure(measureind).Name = stri;
  SET(no).Measure(measureind).LongName = lstr;
else
  myfailed('Invalid name.',DATA.GUI.Segment);
  return;
end

drawfunctions('drawno',no)
segment('updatemeasurement');

%----------------------------------
function doreturn = noobjectshelper
%----------------------------------
%Helper function that checks for LastObject

global DATA

if isempty(DATA.LastObject)
  myfailed('No objects');
  doreturn = true;
else
  doreturn = false;
end

%------------------------------
function [stri,lstr] = measureasklabel(measureind)
%------------------------------
%Asks for a label of a measurement.

global DATA

%Call overloaded method
[stri,lstr] = DATA.measureasklabel(measureind);

segment('updatemeasurement');


%---------------------------------------
function interpdeletepointall(pointtype)
%---------------------------------------
%remove interpolation point closest to click
global DATA SET
no = DATA.ViewPanels(DATA.CurrentPanel);
if nargin < 1
  pointtype = DATA.LastObjectType;
end

SET(no).([pointtype,'X']) = [];
SET(no).([pointtype,'Y']) = [];

drawfunctions('drawno',no);

%-----------------------------------------
function clearallslices_Callback %#ok<DEFNU>
%-----------------------------------------
%Clear all segmentation, both endo and epi, lv and rv mar and scar.

segmentation('clearslices_Callback');

%-----------------------------------------
function segmentclearall_Callback(silent)  %#ok<DEFNU>
%-----------------------------------------
%Clear all segmentation, both endo and epi, lv and rv mar and scar.

global DATA

if nargin == 0 || not(silent)
  msg = dprintf('This removes all existing segmentation of LV, RV, ROI, MaR and scar. Are you sure?');
  if ~yesno(msg,[],DATA.GUI.Segment)
    return;
  end
end

viability('viabilityclear_Callback');
roi('roiclearall_Callback')
segmentation('clearall_Callback');
mar('clearall_Callback');


%-----------------------------------------------------
function segmentclearallbutsystolediastole_Callback
%-----------------------------------------------------
%Clears all segmentation in all timeframes but systole and diastole.
global DATA SET NO

if SET(NO).TSize<2
  myfailed('Not timeresolved data, aborting.',DATA.GUI.Segment);
  return;
end

if isequal(SET(NO).EDT,SET(NO).EST)
  myfailed('Systole and diastole occurs at the same time frame, aborting.',DATA.GUI.Segment);
  return;
end

%Create index structure
ind = true(1,SET(NO).TSize);
ind(SET(NO).EDT) = false;
ind(SET(NO).EST) = false;
arg = struct('endo',true,'epi',true,'rvendo',true,'rvepi',true);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind);
segmentation('removeallinterp_Callback',true,[],arg,indarg);

if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX(:,ind,:) = NaN;
  SET(NO).EndoY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX(:,ind,:) = NaN;
  SET(NO).EpiY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX(:,ind,:) = NaN;
  SET(NO).RVEndoY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX(:,ind,:) = NaN;
  SET(NO).RVEpiY(:,ind,:) = NaN;
end

SET(NO).EndoDraged(ind,:) = false;
SET(NO).EpiDraged(ind,:) = false;

segment('updatevolume');
drawfunctions('drawcontours',DATA.CurrentPanel);
drawfunctions('drawinterp',DATA.CurrentPanel);

%--------------------------------------------
function segmentclearallbutdiastole_Callback
%--------------------------------------------
%Clears all segmentation in all timeframes but diastole.
global DATA SET NO

if SET(NO).TSize<2
  myfailed('Not timeresolved data, aborting.',DATA.GUI.Segment);
  return;
end

%Create index structure
ind = true(1,SET(NO).TSize);
ind(SET(NO).EDT) = false;
arg = struct('endo',true,'epi',true,'rvendo',true,'rvepi',true);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind);
segmentation('removeallinterp_Callback',true,[],arg,indarg);

if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX(:,ind,:) = NaN;
  SET(NO).EndoY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX(:,ind,:) = NaN;
  SET(NO).EpiY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX(:,ind,:) = NaN;
  SET(NO).RVEndoY(:,ind,:) = NaN;
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX(:,ind,:) = NaN;
  SET(NO).RVEpiY(:,ind,:) = NaN;
end

SET(NO).EndoDraged(ind,:) = false;
SET(NO).EpiDraged(ind,:) = false;

segment('updatevolume');
drawfunctions('drawcontours',DATA.CurrentPanel);
drawfunctions('drawinterp',DATA.CurrentPanel);


%------------------------------------
function result = updateintersections_Callback(slices,no,type) %#ok<DEFNU>
%------------------------------------
% Calculates the intersection of a endocardial segmentation
% and SET(no) for
% slices == 'all'
% or
% slices == 'current'
%
% Allowed segmenation 'type' is 'LVEndo' or 'RVEndo'
%
% Calculates intersection for SET(NO) if nargin==1.
% Returns false if calculations fails

% Marten Larsson, June, 2009

global DATA NO SET

if nargin == 1
  no = NO;
  type = 'LV Endocardial';
end

result = true;
segno = [];

switch type
  case 'LV Endocardial'
    % Find sets with LV endocardial segmentation
    for i=1:length(SET)
      if ~isempty(SET(i).EndoX) && any(any(any(~isnan(SET(i).EndoX))))
        segno = [segno i]; %#ok<AGROW>
      end
    end
  case 'RV Endocardial'
    % Find sets with RV endocardial segmentation
    for i=1:length(SET)
      if ~isempty(SET(i).RVEndoX) && any(any(any(~isnan(SET(i).RVEndoX))))
        segno = [segno i]; %#ok<AGROW>
      end
    end
end

if isempty(segno)
  myfailed('No segmentation found',DATA.GUI.Segment);
  result = false;
  return;
end

if length(segno)>1
  segstring = cell(length(segno),1);
  for i = 1: length(segno)
    segstring{i} = dprintf('Image stack %d - %s',segno(i),SET(segno(i)).ImageType);
  end
  
  [selected, ok] = listdlg('ListString',segstring,...
    'SelectionMode','single',...
    'PromptString','Chose image stack to use:',...
    'ListSize',[300 100]);
  
  if ~ok
    result = false;
    return;
  end
  segno = segno(selected);
end

calcsegintersect(segno,no,slices,type);
viewfunctions('setview'); %drawfunctions('drawall');


%----------------------
function selectlvmethod %#ok<DEFNU>
%----------------------
%callback from 3DPrint auto LV icon to select if apply CT or MR auto LV seg

global SET NO

no = NO;
ismr = strfind(SET(no).ImagingTechnique,'MR');
isct = strfind(SET(no).ImagingTechnique,'CT');

if not(isempty(ismr)) && isempty(isct)
  lvsegmentation;
elseif not(isempty(isct)) && isempty(ismr)
  ct.ctlicensecheck('CTLVSegmentation');
else
  %user select if it is CT or MR or cancel
  menuitems{1} = 'MR';
  menuitems{2} = 'CT';
  menuitems{3} = 'Other';
  answer = mymenu('Select imaging technique for the selected image stack.',menuitems);
  if answer
    if answer == 1
      lvsegmentation;
    elseif answer == 2
      ct.ctlicensecheck('CTLVSegmentation');
    elseif answer == 3
      myfailed('No automatic LV segmentation method for other imaging techniques avaliable.');
    end
  else
    disp('LV segmentation canceled.');
  end
end

%---------------------------------------------
function calccalciummask_Callback
%---------------------------------------------
global SET NO
sizesetstruct=size(SET);

temp=[];
tempCaSc=[];

temp=[];
tempCaSc=[];

if sizesetstruct(1,2)==1
  m=1;
else
  m=2;
end

for n=1:sizesetstruct(1,m)
  temp=[temp isempty(SET(n).EndoX)==0];
  tempCaSc=[tempCaSc contains(SET(n).SeriesDescription,'CaSc')];
end
sumtempCaSc= sum(tempCaSc);
if sum(temp)==0
  myfailed('there is no contrastimage')
  return
end
if sum(temp)==1
  seg= find(temp);
elseif sum(temp)==2
  findindex1=find(temp);
  seg = mymenu('Choose segmentation stack',['Image stack no ' num2str(SET(findindex1(1)).Linked)], ['Image stack no ' num2str(SET(findindex1(2)).Linked)]);
end

if sumtempCaSc==0
  myfailed('There is no Ca image')
  return
end
if sumtempCaSc==1
  no=find(tempCaSc);
elseif sumtempCaSc==2
  findindex=find(tempCaSc);
  no = mymenu('In which image stack do you want to detect calcium?',['Image stack no ' num2str(SET(findindex(1)).Linked)], ['Image stack no ' num2str(SET(findindex(2)).Linked)]);
  no=findindex(no);
elseif sumtempCaSc==3
  no = mymenu('In which image stack do you want to detect calcium?',['Image stack no ' num2str(SET(findindex(1)).Linked)], ['Image stack no ' num2str(SET(findindex(2)).Linked)], ['Image stack no ' num2str(SET(findindex(3)).Linked)]);
  no=findindex(no);
end

callbackfunctions('calccalciummask',no,seg,0)

%---------------------------------------------
function splitcalciummask_Callback
%---------------------------------------------?

global SET NO
sizesetstruct=size(SET);

temp=[];
tempCaSc=[];

if sizesetstruct(1,2)==1
  m=1;
else
  m=2;
end

for n=1:sizesetstruct(1,m)
  temp=[temp isempty(SET(n).EndoX)==0];
  tempCaSc=[tempCaSc contains(SET(n).SeriesDescription,'CaSc')];
  %tempCaSc=[tempCaSc contains(SET(n).SeriesDescription,'CorASeq')];
end
sumtempCaSc= sum(tempCaSc);
if sum(temp)==0
 % myfailed('There is no contrast image')
  seg=0;
end
if sum(temp)==1
  seg= find(temp);
elseif sum(temp)>1
  %seg = mymenu('In which image stack do you want to use for the segmentation','e', 'e', 'Image stack no 3','Image stack no 4');
  myfailed('several contrastimages')
  return
end

if sumtempCaSc==0
  myfailed('There is no Ca image')
  return
end
if sumtempCaSc==1
  no=find(tempCaSc);
elseif sumtempCaSc==2
  findindex=find(tempCaSc);
  no = mymenu('In which image stack do you want to detect calcium?',['Image stack no ' num2str(SET(findindex(1)).Linked)], ['Image stack no ' num2str(SET(findindex(2)).Linked)]);
  no=findindex(no);
elseif sumtempCaSc==3
  no = mymenu('In which image stack do you want to detect calcium?',['Image stack no ' num2str(SET(findindex(1)).Linked)], ['Image stack no ' num2str(SET(findindex(2)).Linked)],['Image stack no ' num2str(SET(findindex(3)).Linked)]);
  no=findindex(no);
  
end


callbackfunctions('calccalciummask',no,seg,1)


%---------------------------------------------
function rawca_Callback
%---------------------------------------------?

global SET NO
sizesetstruct=size(SET);

temp=[];
tempCaSc=[];

if sizesetstruct(1,2)==1
  m=1;
else
  m=2;
end

for n=1:sizesetstruct(1,m)
 
  tempCaSc=[tempCaSc contains(SET(n).SeriesDescription,'CaSc')];
  
end
sumtempCaSc= sum(tempCaSc);

seg=0;

if sumtempCaSc==0
  myfailed('There is no Ca image')
  return
end
if sumtempCaSc==1
  no=find(tempCaSc);
elseif sumtempCaSc==2
  findindex=find(tempCaSc);
  no = mymenu('In which image stack do you want to detect calcium?',['Image stack no ' num2str(SET(findindex(1)).Linked)], ['Image stack no ' num2str(SET(findindex(2)).Linked)]);
  no=findindex(no);
elseif sumtempCaSc==3
  no = mymenu('In which image stack do you want to detect calcium?',['Image stack no ' num2str(SET(findindex(1)).Linked)], ['Image stack no ' num2str(SET(findindex(2)).Linked)],['Image stack no ' num2str(SET(findindex(3)).Linked)]);
  no=findindex(no);
  
end


callbackfunctions('calccalciummask',no,seg,1)



%---------------------------
function calccalciummask(no,seg,onORoff) %#ok<DEFNU>
%---------------------------
%Computes calciummask, will be replace by code from Lisa
radius1=[];
global SET NO

if nargin==0
  no = NO;
end
%segim=SET(seg).IM;
im=SET(no).IM;
im = calcfunctions('calctruedata',im,no); %g?r om till Hounsfieldunits
outputmask = uint8((im>=130));

imbin = (im>=130); %allt ?ver 130 Hu r?knas som kalk
imbin = squeeze(imbin);
bw = bwconncomp(imbin);


if bw.NumObjects>0
  for k=1:bw.NumObjects
    thislist=bw.PixelIdxList{1,k};
    [R]=numel(thislist);
    %   Cascore i mm3
    Cascore1 = R*SET(no).ResolutionX*SET(no).ResolutionY*(SET(no).SliceThickness+SET(no).SliceGap);
    if  Cascore1>3
      outputmask(bw.PixelIdxList{k})=uint8(2);
    end
    
  end
end

SET(no).CT.CalciumMask = outputmask;
if seg~=0
lengthsofthisIm=[];
% Pick outsets of 4 or more
sizesetstruct=size(SET);

threedvol=uint8(zeros(SET(seg).XSize,SET(seg).YSize,SET(seg).ZSize));
threedvol2=uint8(zeros(SET(seg).XSize,SET(seg).YSize,SET(seg).ZSize));
threedvolM=uint8(zeros(SET(seg).XSize,SET(seg).YSize,SET(seg).ZSize));
threedvolA=uint8(zeros(SET(seg).XSize,SET(seg).YSize,SET(seg).ZSize));
SET(no).CT.CalciumPenMask=uint8(zeros(SET(seg).XSize,SET(seg).YSize,SET(seg).ZSize));
shortaxisImagelocation=0;
EXlocation= 0;

for n=1:sizesetstruct(1,2)
  
  if strcmp(SET(n).ImageViewPlane,'Short-axis')==1 && contains(SET(n).SeriesDescription,'CaSc')==1
    shortaxisImagelocation=n;
  end
  
  if isempty(SET(n).EndoX)==0
    EXlocation= n;
  else
    EXlocation= seg;
  end
  if strcmp(SET(n).ImageViewPlane,'Transversal')==1 && contains(SET(n).SeriesDescription,'CorASeq')==1
    transversalContrastLocation=n;
  end
end


% kolla vid vilken slice k segmenteringen b?rjar

sizeOfEndoSeg=size(SET(EXlocation).EndoX);
for k=1:sizeOfEndoSeg(1,3)
  if isnan(SET(EXlocation).EndoX(:,:,k))==0
    hereIsNotNan=k;
    break;
  end
end



%  Mask making

meanX=mean(SET(EXlocation).EndoX(:,:,hereIsNotNan));
meanY=mean(SET(EXlocation).EndoY(:,:,hereIsNotNan));

%ta bild 8 snitt upp

imbin1=[];
imbin2=[];
Cascore1=0;
Cascore2=0;

firstZslice=hereIsNotNan-8;
lastZslice=hereIsNotNan+8;
mitsizes=[];
for slice=firstZslice:lastZslice
  
  im1=imadjust(SET(EXlocation).IM(:,:,1,slice));
  
  %ber?kna treshold
  im1=calcfunctions('calctruedata',im1,EXlocation); %g?r om till Hounsfieldunits
  iim=(im1+1024)/ max(max(im1+1024));
  JJ=imhist(iim);
  minloc=178+find(JJ(179:230)==min(JJ(179:230)));
  maxloc=153+find(JJ(154:minloc(1))==max(JJ(154:minloc(1))));
  if length(maxloc(end))>=1 && JJ(maxloc(end))>700
    T1=round((maxloc(end)+minloc(1))/2)/256;
  else
    T1=minloc(1)/256;
  end
  T=T1*max(max(im1+1024))-1024;
  im2=im1>T;
  im2=imfill(im2, 'holes');
  
  imageSizeX =SET(EXlocation).YSize;
  imageSizeY = SET(EXlocation).XSize;
  %om slice>=hereIsNotNan anv?nd befintlig lv segmentering 
  if slice>=hereIsNotNan
      finalmitralmask=zeros(imageSizeY,imageSizeX);
      endoX=round(SET(EXlocation).EndoX(:,:,slice));
      endoY=round(SET(EXlocation).EndoY(:,:,slice));
      for l=1:numel(SET(EXlocation).EndoY(:,:,slice)) 
      finalmitralmask(endoX(l),endoY(l))=1;
      end
      
      se = strel('square',10);
  finalmitralmask=(imdilate(finalmitralmask,se));
  finalmitralmask=imfill(finalmitralmask);
%       figure
%       imshow(finalmitralmask)
  else
  %hitta kantpunkter
  scanwidth=10;
  
  for k=1:100
    sumyXled=0;
    
    sumyXled= sum(im2(round(meanX)-scanwidth:round(meanX)+scanwidth, round(meanY)+k));
    
    if sumyXled/(scanwidth*2)<0.5
      edgeX=round(meanY)+k;
      break
    end
  end
  
  isedgeYcreated=0;
  for k=1:100
    sumyYled=0;
    sumyYled= sum(im2(round(meanX)+k, round(meanY)-scanwidth:round(meanY)+scanwidth));
 
    if sumyYled/(scanwidth*2)<0.5
      edgeY=round(meanX)+k;
       isedgeYcreated=1;
      break
    end
   
  end

 isedgeY2created=0;
  for k=1:100
    sumyYled=0;
    sumyYled= sum(im2(round(meanX)-k, round(meanY)-scanwidth:round(meanY)+scanwidth));
 
    if sumyYled/(scanwidth*2)<0.5
      edgeY2=round(meanX)-k;
       isedgeY2created=1;
      break
    end
  end
  
 if isedgeY2created==0 || isedgeYcreated==0
    newmeanX=meanX;
 else
  newmeanX=(edgeY+edgeY2)/2;
 end
  meanX=newmeanX;
  
  %Create circlemask
  imageSizeX =SET(EXlocation).YSize;
  imageSizeY = SET(EXlocation).XSize;
  [columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
  centerX = meanY;
  centerY = meanX;
  if isedgeYcreated==0
    radius = round(max([edgeX-meanY meanX-edgeY2]));
  else
  radius = round(max([edgeX-meanY edgeY-meanX]));
  end
  radius1=[radius1; radius];
  if radius<30
      radius=30;
  end
  if radius>55
      radius=55;
  end
  circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
  
  immm=imfuse(im2,circlePixels);
  immm=im2.*circlePixels;
  
  linemask= columnsInImage> (round(meanY)-radius);
  
  
  
  L=bwlabel(im2.*linemask);
 
  centerpoints=[];
  centerdists=[];
  sizes=[];
  for k=1:max(max(L))
    [C,R]=find(L==k);
    sizes= [sizes; length(R)];
    centerpoints=[centerpoints; mean(R) mean(C)];
    centerdists= [centerdists; sqrt((mean(R)-meanY)^2+(mean(C)-meanX)^2)];
  end
  maxsize=max(sizes);
  bigobjindex=sizes>(maxsize/4);
  O=find(bigobjindex);
  centerdists2=centerdists(O);
  centerpoints2=centerpoints(O,:);
  EndoObjNbr=find(centerdists==min(centerdists2));
  L= L==EndoObjNbr;
  if maxsize>8000
      L=circlePixels;
  end
  mitsizes=[mitsizes maxsize];
 
%   if sum(sum(L))>15000
%     [C,R]=find(L);
%     xline=mean(R);
%     linemask= columnsInImage> xline-radius;
%     L=L.*linemask;
%     L=bwlabel(L);
%     centerpoints=[];
%     centerdists=[];
%     sizes=[];
%     for k=1:max(max(L))
%       [C,R]=find(L==k);
%       sizes= [sizes; length(R)];
%       centerpoints=[centerpoints; mean(R) mean(C)];
%       centerdists= [centerdists; sqrt((mean(R)-meanY)^2+(mean(C)-meanX)^2)];
%     end
%     maxsize=max(sizes);
%     bigobjindex=sizes>(maxsize/2);
%     O=find(bigobjindex);
%     centerdists2=centerdists(O);
%     centerpoints2=centerpoints(O,:);
%     EndoObjNbr=find(centerdists==min(centerdists2));
%     L= L==EndoObjNbr;
%   end
  finalmitralmask=L;
  
  end
  [yM,xM]=find(finalmitralmask);
 
  %    Aortamask
  
  linemaskA= columnsInImage<min(xM)+5;
  
  pointtolookfor=min(xM);
%    figure
%   imagesc(finalmitralmask)
%   hold on
%   plot(pointtolookfor,meanX,'o')
  iseroded=0;
  L=bwlabel(im2.*linemaskA);
  L2=L;
  for k=1:max(max(L))
    if sum(sum(L==k))>3000 
      
      se = strel('square',7);
      L = imerode(L,se);
      L=bwlabel(L>0);
      iseroded=iseroded+1;
      break
    end
  end
  for k=1:max(max(L))
    if sum(sum(L==k))>3000
      
      se = strel('square',7);
      L = imerode(L,se);
      L=bwlabel(L>0);
      iseroded=iseroded+1;
      break
    end
  end
  centerpoints=[];
  centerdists=[];
  sizes=[];
  isitnegative=[];
  
  for k=1:max(max(L))
    [C,R]=find(L==k);
    sizes= [sizes; length(R)];
    centerpoints=[centerpoints; mean(R) mean(C)];
    centerdists= [centerdists; sqrt((mean(R)-pointtolookfor)^2+(mean(C)-meanX)^2)];
    isitnegative=[isitnegative; (mean(R)-pointtolookfor)<0];
  end

  maxsize=max(sizes);
  bigobjindex=sizes>(maxsize/10);
  O=find(bigobjindex);
  centerdists2=centerdists(O);
  centerpoints2=centerpoints(O,:);
  if isitnegative(find(centerdists==min(centerdists2)))==0
      centerdists2=sort(unique(centerdists2));
      EndoObjNbr=find(centerdists==centerdists2(end-1));
  else
  EndoObjNbr=find(centerdists==min(centerdists2));
  end
  L2=L;
  L= L==EndoObjNbr;
  
  if sum(sum(L))>8000
    
    centerX = pointtolookfor;
    centerY = meanX;
    radius = 50;
    circlePixels = (rowsInImage - centerY).^2 ...
      + (columnsInImage - centerX).^2 <= radius.^2;
    finalaortamask=L.*circlePixels;
  else
    finalaortamask=L ;
  end
  while iseroded>0
    finalaortamask= imdilate(finalaortamask,se);
    iseroded=iseroded-1;
  end
  
  % check if the mask is bad
  
  [RA,CA]=find(finalaortamask);
  [RM,CM]=find(finalmitralmask);
  
  lengthsofthisIm=[lengthsofthisIm sqrt((pointtolookfor-max(CA))^2 + (meanX-mean( RA(find(CA==max(CA)))))^2)];
  
  for i=1:length(CA)
    dist(i) = sqrt((pointtolookfor-CA(i))^2 + (meanX-RA(i))^2) ;
  end
  
  
 
  
  if mean(RA)>mean(RM)
   centerdists2=sort(unique(centerdists2));
   if length(centerdists2)>1
   centerdists2(1)=200;
   end
      EndoObjNbr=find(centerdists==min(centerdists2));
            
if sum(sum(L2==EndoObjNbr))>500 && min(dist(find(CA==max(CA))))>30
      L= L2==EndoObjNbr;

  if sum(sum(L))>8000
    
    centerX = pointtolookfor;
    centerY = meanX;
    radius = 50;
    circlePixels = (rowsInImage - centerY).^2 ...
      + (columnsInImage - centerX).^2 <= radius.^2;
    finalaortamask=L.*circlePixels;
  else
    finalaortamask=L ;
  end
  while iseroded>0
    finalaortamask= imdilate(finalaortamask,se);
    iseroded=iseroded-1;
  end
end
  end
   if min(dist(find(CA==max(CA))))>35 || slice>hereIsNotNan+3
    finalaortamask= finalaortamask*0;
  end
  sizeoA(slice-firstZslice+1)= sum(sum(finalaortamask));
%   if sum(sum(finalaortamask))>5500
%     finalaortamask= finalaortamask*0;
%   end
  
  [yA,xA]=find(finalaortamask);

  
  % Save to mask & Dilate

  outputmask2=uint8(((finalaortamask+finalmitralmask)>0));
  outputmaskMitral=uint8(((finalmitralmask)>0));
  outputmaskAortic=uint8(((finalaortamask)>0));
  

 simm=SET(seg).IM;
simm=squeeze(simm);
% if slice>=(hereIsNotNan-8) && slice<(hereIsNotNan+8)  
%   figure
%   subplot(1,2,1)
%    imagesc(simm(:,:,slice))
%    colormap('gray');
%     hold on
%       
%     [~,c] = contour((finalaortamask)>0);
%     c.LineColor= [0 0 1];
%   
%      [~,c] = contour((finalmitralmask)>0);
%     c.LineColor= [1 0 0];
% slice
%      end
  
 
  se = strel('square',20);
  outputmask3=imfill(imdilate(outputmask2,se));
  outputmask2=edge(imfill(imdilate(outputmask2,se)));
  
  outputmaskMitral=imfill(imdilate(outputmaskMitral,se));
  outputmaskAortic=imfill(imdilate(outputmaskAortic,se));
  
  [yM,xM]=find(outputmaskMitral);
  [yA,xA]=find(outputmaskAortic);
  
  if isempty(max(xA))
      linemaskA= uint8(columnsInImage < min(xM));
  else
  linemaskA= uint8(columnsInImage < max(xA));
  end
  linemaskM= uint8(~linemaskA);
  
  threedvol(:,:,slice) = uint8(outputmask2);
  threedvol2(:,:,slice) = uint8(outputmask3);
  threedvolM(:,:,slice) = uint8(outputmask3.*linemaskM);
  threedvolA(:,:,slice) = uint8(outputmask3.*linemaskA);
  
%  figure
%    imagesc(segim(:,:,1,slice))
%    colormap('gray');
%     hold on
%       
%     [~,c] = contour((outputmask3.*linemaskA)>0);
%     c.LineColor= [0 1 1]; [~,c] = contour((outputmask3.*linemaskM)>0);
%     c.LineColor= [0 0 1]; plot(pointtolookfor,meanX,'o');

  
end
h = waitbar(0,'Please wait.');

threedvol=uint8(threedvol);
outputmask2 = segment3dp.resamplestack(seg,no,threedvol,h); %edge
threedvol2=uint8(threedvol2);
outputmask3 = segment3dp.resamplestack(seg,no,threedvol2,h);
threedvolM=uint8(threedvolM);
outputmaskM = segment3dp.resamplestack(seg,no,threedvolM,h);
threedvolA=uint8(threedvolA);
outputmaskA = segment3dp.resamplestack(seg,no,threedvolA,h);

close(h);

[~,~,zet]=size(outputmask2);
if onORoff==0
  for i=1:zet
    SET(no).CT.CalciumMask(:,:,1,i) =uint8(3)*uint8(uint8(SET(no).CT.CalciumMask(:,:,1,i))==2).*uint8(outputmask3(:,:,i))+  4* uint8((outputmask2(:,:,i)));
  end
elseif onORoff==1
  for i=1:zet
    SET(no).CT.CalciumMask(:,:,1,i) =uint8(7)*uint8(uint8(SET(no).CT.CalciumMask(:,:,1,i))==2).*uint8(outputmaskM(:,:,i)) + uint8(6)*uint8(uint8(SET(no).CT.CalciumMask(:,:,1,i))==2).*uint8(outputmaskA(:,:,i)) +  uint8(4)*uint8(outputmask2(:,:,i));
  end
end
end

    
disp('dsds');

%---------------------------
function removeallcalciumsegmentation_Callback
%---------------------------
%removes all calciumsegmentation

global DATA SET NO
no=NO;
if ~isfield(SET(no).CT, 'CalciumMask')
  myfailed('There is no calcium mask in this image stack');
  return
end

threes=find (SET(no).CT.CalciumMask==3);
sixs=find (SET(no).CT.CalciumMask==6);
sevens=find(SET(no).CT.CalciumMask==7);
eights=find (SET(no).CT.CalciumMask==8);

CalciumMask=SET(no).CT.CalciumMask;

if ~isempty(threes)
  CalciumMask(threes)=5;
end
if ~isempty(sixs)
  CalciumMask(sixs)=5;
end
if ~isempty(sevens)
  CalciumMask(sevens)=5;
end
if ~isempty(eights)
  CalciumMask(eights)=5;
end

SET(no).CT.CalciumMask=CalciumMask;


DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel);

%---------------------------
function removeallcalciumsegmentationonthisslice_Callback
%---------------------------
%removes all calciumsegmentation

global DATA SET NO
no=NO;
if ~isfield(SET(no).CT,'CalciumMask')
  myfailed('There is no calcium mask in this image stack');
  return
end

slice=SET(no).CurrentSlice;
threes=find (SET(no).CT.CalciumMask(:,:,1,slice)==3);
sixs=find (SET(no).CT.CalciumMask(:,:,1,slice)==6);
sevens=find(SET(no).CT.CalciumMask(:,:,1,slice)==7);
eights=find (SET(no).CT.CalciumMask(:,:,1,slice)==8);

CalciumMask=SET(no).CT.CalciumMask(:,:,1,slice);

if ~isempty(threes)
  CalciumMask(threes)=5;
end
if ~isempty(sixs)
  CalciumMask(sixs)=5;
end
if ~isempty(sevens)
  CalciumMask(sevens)=5;
end
if ~isempty(eights)
  CalciumMask(eights)=5;
end

SET(no).CT.CalciumMask(:,:,1,slice)=CalciumMask;


DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel);
