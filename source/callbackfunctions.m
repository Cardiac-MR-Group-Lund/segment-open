function varargout = callbackfunctions(varargin)
% This .m file has the ambition of containing all Callback functions

% Split out by Klas

%#ok<*GVMIS> 

%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%------------------------
function valvetracking_Callback
%------------------------

valvetracking.valvetrackingml('track');

%------------------------
function createcutplaneobject_Callback
%------------------------
%Creates a cutplane object with selected object as parent
global SET NO

O = SET(NO).LevelSet.Object;
O.createcutplaneobject_Callback;

%------------------------
function fdabutton_Callback
%------------------------
%Manage the GUI appareance in FDA theme
global DATA SET

%---Word-around to unfocus the radiobutton
DATA.requestfocus;

%--- Save theme to preferences
DATA.Pref.RunFDAVersion = true;
silent = true;
segpref('save_Callback',silent);

%--- Get all GUI handles
h = DATA.Handles;
handlenames = fieldnames(h);

%--- Verify that current tab is FDA cleared, otherwise switch to LV tab
if ~isempty(h.toggleiconholder.clickedicon)
  if ~h.toggleiconholder.clickedicon.isFDACleared
    executeicon = true;
    h.toggleiconholder.indent('ribbonfunction',executeicon);
  end
end

%--- Make non-FDA cleared menus invisible
for loop = 1:numel(handlenames)
  currenthandle = h.(handlenames{loop});
  try
    if isa(currenthandle,'matlab.ui.container.Menu')
      set(currenthandle,'Visible','off');
    end
  catch me
    mydispexception(me);
  end
end
set(getfdaclearedmenus,'visible','on');

%--- Disable non-FDA cleared icons
fdapreferences = DATA.Pref.RunFDAVersion;
if ~isempty(SET)  %no data loaded, no need to update the icons placeholders
  h.toggleiconholder.enable('',fdapreferences);
  h.permanenticonholder.enable('',fdapreferences);
  h.configiconholder.enable('',fdapreferences);
end

%------------------------
function researchbutton_Callback
%------------------------
%Manage the GUI appareance in Research theme
global DATA SET

%---Word-around to unfocus the radiobutton
DATA.requestfocus;

%--- Save theme to preferences
DATA.Pref.RunFDAVersion = false;
silent = true;
segpref('save_Callback',silent);

%--- Get all GUI handles
h = DATA.Handles;
handlenames = fieldnames(h);

%--- Make all menus invisible
for loop = 1:numel(handlenames)
  currenthandle = h.(handlenames{loop});
  try
    if isa(currenthandle,'matlab.ui.container.Menu')
      if strcmp(currenthandle.Visible,'off')
        set(currenthandle,'Visible','on');
      end
    end
  catch me
    mydispexception(me);
  end
end

%--- Enable all icons
fdapreferences = DATA.Pref.RunFDAVersion;
if ~isempty(SET)  %no data loaded, no need to update the icons placeholders
  h.toggleiconholder.enable('',fdapreferences);
  h.permanenticonholder.enable('',fdapreferences);
  h.configiconholder.enable('',fdapreferences);
end

%------------------------------------------
function handles = getfdaclearedmenus
%------------------------------------------
% Return a vector containing the handles to all the Menu items that FDA
% cleared
global DATA

h = DATA.Handles;

%List of File menu items
filehandles = [ ...
  h.filemenu ...
  h.fileopenfromdiscmenu ...
  h.fileloadnextmenu ...
  h.openpatientdatabasemenu ...
  h.fileloadsegmentationmenu ...
  h.filesaveallmenu ...
  h.filesavetodatabasemenu ...
  h.filesaveallasmenu ...
  h.filesavesubmenu ...
  h.filesavesegdicom ...
  h.filesavecurrentmenu ...
  h.filesavesegmentationmenu ...
  h.filecloseallmenu ...
  h.fileclosemultiplemenu ...
  h.fileclosecurrentimagestack ...
  h.fileresetguimenu ...
  h.filequitmenu ...
  ];

%List of Edit menu items
edithandles = [...
  h.editmenu...
  h.editselectallslices ...
  h.editunselectallslices ...
  h.editautoesedmenu ...
  h.editsetedmenu ...
  h.editsetesmenu ...
  h.editsetfirsttimeframemenu ...
  h.editsetfirsttimeframeforselectedmenu ...
  ];

%List of Image Tools menu items
toolshandles = [...
  h.toolsmenu ...
  h.toolsremovetimeframesmenu ...
  h.editremovecurrenttimeframe ...
  h.editremovenexttimeframe ...
  h.editremoveallbutthistimeframemenu ...
  h.editremoveprevioustimeframe ...
  h.editremoveduplicatetimeframe ...
  h.editremoveallbutedes ...
  h.toolsremovetimeslicesmenu ...
  h.editremovecurrentslicemenu ...
  h.editremoveslicesmenu ...
  h.editremoveunselectedmenu ...
  h.duplicatebasalslicemenu ...
  h.duplicateapicalslicemenu ...
  h.setcolormapmenu ...
  h.originalcolormapmenu ...
  h.colormapgraymenu ...
  h.hsvcolormapmenu ...
  h.jetcolormapmenu ...
  h.hotcolormapmenu ...
  h.spectcolormapmenu ...
  h.gtcolormapmenu ...
  h.toolsfliprotatemenu ...
  h.toolsflipxmenu ...
  h.toolsflipymenu ...
  h.toolsfliptmenu ...
  h.toolsflipzmenu ...
  h.toolsrotate90rightmenu ...
  h.toolsflipxzmenu ...
  h.toolsflipztmenu ...
  h.toolsflipxmirrormenu ...
  h.toolsresamplestackmenu ...
  h.toolsupsamplemenu ...
  h.toolsupsampleslicesmenu ...
  h.toolsupsampletemporalmenu ...
  h.toolscropstackmenu ...
  h.toolsunlinkimagesmenu ...
  h.toolsadjustimagedetailsmenu ...
  ];

%List of ROI menu items
roihandles = [...
  h.roimenu ...
  h.roiimportmenu ...
  h.roicopyall ...
  h.flowaddfixsizeroi ...
  h.roisetroilabelmenu ...
  h.roisetroilcolormenu ...
  h.roideletetemplatemenu ...
  h.roideleteroiintimeframemenu ...
  h.roideleteroinalltimeframesbutthismenu ...
  h.roicopyendomenu ...
  h.roicopyepimenu ...
  h.roicopytolvmenu ...
  h.roisignalintensityanalysis ...
  h.roihistogrammenu ...
  h.roiexportroivaluesmenu ...
  h.roithresholdnumericmenu ...
  h.roithresholdnumericvisual ...
  h.roivisualthresholdingmenu ...
  ];

%List of Annotations menu items
annotationshandles = [...
  h.measuremenu ...
  h.annotationsimportpointsmenu ...
  h.annotationsrenamepointsmenu ...
  h.annotationsclearpointsmenu ...
  h.annotationstemporalfiltermenu ...
  h.annotationsexportpointsmenu ...
  h.exportmeasurements ...
  ];

%List of Segmentation menu items
segmentationhandles = [...
  h.segmentmenu ...
  h.lvmenu ...
  h.lvaisaxmenu ...
  h.lvaisaxedesmenu ...
  h.lvaisaxslicesmenu ...
  h.importsegmenu ...
  h.interchangesegmenu ...
  h.calcendointersectionall ...
  h.calcendointersectioncurrent ...
  h.copyrefinemenu ...
  h.toolscopyendoupwardsmenu ...
  h.toolscopyepiupwardsmenu ...
  h.toolscopyendodownwardsmenu ...
  h.toolscopyepidownwardsmenu ...
  h.toolscopytoalltimeframes ...
  h.clearseginslicesmenu ...
  h.clearsegmentationendo ...
  h.clearsegmentationendothismenu ...
  h.clearsegmentationepi ...
  h.clearsegmentationepithismenu ...
  h.alternativelvmenu ...
  h.alternativelvsegmenu ...
  h.autoendomenu ...
  h.autoepimenu ...
  h.alternativelvrefineendomenu ...
  h.alternativelvrefineepimenu ...
  h.setlongaxismotionmenu ...
  h.autoestimatepapvolumemenu ...
  h.estimatepapvolumefromroimenu ...
  h.adjustpapthresholdmenu ...
  h.segmentresetpapilaryvolumemenu ...
  h.rvmenu ...
  h.rvaisegmenu ...
  ...h.rvinterchangesegmenu ...
  h.copyendolvtorvmenu ...
  h.copyepilvtorvmenu ...
  h.rvclearseginslicesmenu ...
  h.rvclearsegmentationendo ...
  h.rvclearsegmentationendothismenu ...
  h.rvclearsegmentationepi ...
  h.rvclearsegmentationepithismenu ...
  ];

%List of MR menu items
mrhandles = [...
  h.mrmenu ...
  h.flowmenu ...
  h.flowcouplestacks ...
  h.flowcreatevelmag ...
  h.flowdeletevelmag ...
  h.flowcreateangio ...
  h.flowdeleteangio ...
  h.flowconcomittantmenu ...
  h.flowsetvenc ...
  h.flowswitchroisignmenu ...
  h.qpqsmenu ...
  h.pwvmenu ...
  h.transittimetoolmenu ...
  h.flowcurvaturemenu ...
  h.flowalternativetrackmenu ...
  h.manualheartbeats ...
  ];

%List of Analysis menu items
analysishandles = [...
  h.analysismenu ...
  h.reportradvelmenu ...
  h.reportslicemenu ...
  ];

%List of View menu items
viewhandles = [...
  h.viewmenu...
  h.viewallimagestacksmenu ...
  ];

%List of Export menu items
exporthandles = [...
  h.exportmenu ...
  h.exporttoclipboardmenu ...
  h.exporttoclipboardmenunoheader ...
  h.exporttoclipboardthisstackmenu ...
  h.exporttoclipboardthisstackmenunoheader ...
  h.exportvolumecurvemenu ...
  h.exportmmodemenu ...
  h.exportcontoursmenu ...
  h.exportcontourstoclipboardmenu ...
  h.exportallcontourstostlmenu ...
  h.exportcontourtostlmenu ...
  h.exportlvcontourtostlmenu ...
  h.exportvolumeofcontoursmenu ...
  h.exportimagemenu ...
  h.exportscreenshotmenu ...
  h.exportmoviemenu ...
  ];

%List of Utility menu items
utilityhandles = [...
  h.utilitymenu ...
  h.anonymizemenu ...
  h.pseydonymizesubject ...
  h.pseydonymizemat ...
  h.pseydonymizedicom ...
  h.utilitysortdicomstackmenu ...
  h.utilitysortfromcdmenu ...
  h.utilitycreatethumbnailsmenu ...
  h.batchexportmatmenu ...
  h.utilitybatchexportresultsmenu ...
  h.utilitybatchexportsegdicommenu ...
  h.utilitybatchexportroimenu ...
  h.utilitybatchexportinfomenu ...
  h.utilitybatchexport17segmentsmenu ...
  h.utilitybatchexportsegmentwallthicknessmenu ...
  h.utilitybatchexportlvslicevolumemenu ...
  h.batchexportniftimenu ...
  h.batchoperationsonmatmenu ...
  h.batchlvsegmentation ...
  h.batchrvsegmentation ...
  h.batchlvrvsegmentation ...
  h.batchclearsegmentation ...
  ];

%List of Setting menu items
settinghandles = [...
  h.settingsmenu...
  h.setuptoolsmenu...
  h.wizardsetuphelpmenu...
  h.licensetoolsmenu...
  h.generatelicensemenu...
  h.hardwarekeysmenu...
  h.createhardwareinfofilemenu...
  h.installhardwarekeymenu...
  h.checkhardwarekeymenu...
  h.uninstallhardwarekeymenu...
  h.licenseinformationmenu...
  h.removelicensemenu...
  h.preferencesmenu...
  ];

%List of Help menu items
helphandles = [...
  h.helpmenu ...
  h.helpaboutmenu ...
  h.hotkeysmenu ...
  h.supportmenu ...
  h.bugreportmenu ...
  h.helpsupportmenu ...
  h.logfilesmenu ...
  h.openlogfilemenu ...
  h.openlogfoldermenu ...
  h.usermanuals ...
  h.helpvideomenu ...
  h.helpversioncheckmenu ...
  h.citationmenu ...
  h.termsconditionsmenu ...
  h.evalcommandmenu ...
  h.researchmanual ...
  h.fdamanual ...
  ];

%List of context menu items
selectcontextmenuhandles = [...
  h.selectcontextmenu ...
  h.headerselectcontextmenu ...
  h.viewonecontext ...
  h.viewmontagecontext ...
  h.viewmontagerowcontext ...
  h.editunselectallslicesmenu ...
  h.clearsegcontextmenu ...
  h.clearsegthismenu ...
  h.clearallcontextmenu ...
  ];
thumbnailscontextmenuhandles = [...
  h.datasetpreviewmenu ...
  h.headerdatasetpreviewcontextmenu ...
  h.updatemenu ...
  h.thumbnailsviewmontagecontextmenu ...
  h.thumbnailsviewsinglecontextmenu ...
  h.setimagedescriptioncontextmenu ...
  h.setheartratecontextmenu ...
  h.deleteselectedmenu ...
  h.duplicatedataset ...
  ];
roicontextmenuhandles = [...
  h.roicontextmenu ...
  h.headerroicontextmenu ...
  h.plotflowcontextmenu ...
  h.eddycurrentcontextmenu ...
  h.deleteroicontextmenu ...
  h.deleteroithiscontextmenu ...
  h.deleteroiallcontextmenu ...
  h.setroilabelcontextmenu ...
  h.setroicolorcontextmenu ...
  h.copyroiupwardscontextmenu ...
  h.copyroidownwardscontextmenu ...
  h.copyroitoalltimeframescontextmenu ...
  h.refineroicontextmenu ...
  h.switchflowsigncontextmenu ....
  ];
measurecontextmenuhandles = [...
  h.measurecontextmenu ...
  h.headermeasurecontextmenu ...
  h.deletethismeascontextmenu ...
  h.deletallmeascontextmenu ...
  h.renamethismeascontextmenu ... 
  ];
pointcontextmenuhandles = [...
  h.pointcontextmenu ...
  h.headerpointcontextmenu ...
  h.deletethispointcontextmenu ...
  h.deleteallpointscontextmenu ...
  h.renamethispointscontextmenu ...
  h.makepointtimeresolvedcontextmenu ...
  h.keeppointcontextmenu ...
  ];
rvcontextmenuhandles = [...
  h.rvcontextmenu ...
  h.headerrvcontextmenu ...
  h.clearrvcontextmenu ...
  h.clearrvthiscontextmenu ...
  h.clearallrvcontextmenu ...
  h.clearrvexceptedescontextmenu ...
  ];
contrastcontextmenuhandles = [...
  h.contrastcontextmenu ...
  h.headercontrastcontextmenu ...
  h.autocontrastcontextmenu ...
  h.resetcontrastcontextmenu ...
  ];
lvcontextmenuhandles = [...
  h.lvcontextmenu ...
  h.headerlvcontextmenu ...
  h.clearlvcontextmenu ...
  h.clearlvthiscontextmenu ...
  h.clearalllvcontextmenu ...
  h.clearlvexceptedescontextmenu ...
  ];
interpcontextmenuhandles = [...
  h.interppointmenu ...
  h.headerinterppointmenu ...
  h.deleteinterppointcontextmenu ...
  h.deleteallinterppointsthiscontextmenu ...
  h.deleteallinterppointscontextmenu ...
  ];

handles = [filehandles edithandles toolshandles roihandles ...
  annotationshandles segmentationhandles mrhandles analysishandles ...
  viewhandles exporthandles utilityhandles settinghandles helphandles ...
  selectcontextmenuhandles thumbnailscontextmenuhandles ...
  roicontextmenuhandles measurecontextmenuhandles pointcontextmenuhandles ...
  rvcontextmenuhandles contrastcontextmenuhandles lvcontextmenuhandles ...
  interpcontextmenuhandles];

%---------------------
function approve_Callback
%-----------------------
%Callback when pressing checkbox "Approve analysis" in Segment CMR
%Setting reviewed to true for AI AutoMate files

global SET DATA

pressed = findindented(DATA.Handles.approveiconholder,'notapproved');
       
if pressed 
  user = helperfunctions('getuser');
  SET(1).Autoloader.ApprovedTime = datestr(now,'yyyy-mm-dd HH:MM:SS');
  SET(1).Autoloader.ApprovedBy = user;
  SET(1).Autoloader.Approved = true;
else
  SET(1).Autoloader.ApprovedTime ='';
  SET(1).Autoloader.ApprovedBy = '';
  SET(1).Autoloader.Approved = false;
end

%---------------------
function orthoview_Callback
%-----------------------
global SET NO

%First we check so that there are any slices because if not there is no
%point

if SET(NO).ZSize==1
  myfailed('Too few slices for orthogonal view.')
  return
end

%this sets all necessary fields
viewfunctions('orthoview',NO);

%This updates the view.
viewfunctions('setview',2,2,ones(1,4)*NO,{'orth','hla','vla','gla'});

%--------------------------------
function fpsminus_Callback
%--------------------------------
% decrease FPS value
offset = -1;
fpsedit_Callback(offset)

%--------------------------------
function fpsplus_Callback
%--------------------------------
% increase FPS value
offset = 1;
fpsedit_Callback(offset)

%--------------------------------
function fpsedit_Callback(offset)
%--------------------------------
% edit FPS
global DATA
if nargin == 0
  offset = 0;
end
fpsvalue = str2double(get(DATA.Handles.fpsedit,'String'));
if isempty(fpsvalue) || isnan(fpsvalue)
  fpsvalue = 25; % Default value if non-valid value
end
fpsvalue = updatefps(fpsvalue,offset);
set(DATA.Handles.fpsedit,'String', num2str(fpsvalue));

%---------------------------------------------
function fpsvalue = updatefps(fpsvalue,offset)
%---------------------------------------------
% Adjust FPS value and ensure it stays within the range
global SET NO
if nargin == 0
  offset = 0;
end
no = NO;
fpsvalue = fpsvalue + offset;
% Check value to be in range [1 60]
fpsvalue = max(1, min(60, fpsvalue));

SET(no).FPS = fpsvalue;

%---------------------
function play_Callback(panels) 
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
switch DATA.ProgramName
  case 'Segment CT'
    iconholdername = 'permanenticonholder';
  otherwise
    iconholdername = 'playiconholder';
end

if maxt == 1
  myfailed('Need a time resolved image stack')
  undent(DATA.Handles.(iconholdername),'play',0)
  return
end

%Try different approach where all stacks are played after NO
firstno = DATA.ViewPanels(panels(1));
startframes = SET(firstno).CurrentTimeFrame;
starttime = now;
beattime = SET(firstno).BeatTime;

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
  viewfunctions('updatedrawlist',p);
end

synchronize = findindented(DATA.Handles.hideiconholder,'synchronize');

DATA.Run = 1;
t = maxt;

if DATA.issegmentcmr || DATA.issegment
  usefps = true;
  currentno = DATA.ViewPanels(DATA.CurrentPanel);
else
  usefps = false;
end

while DATA.Run
  %stopping criteria
  if ~findindented(DATA.Handles.(iconholdername),'play')
    DATA.Run = 0;
  end
  
  for p = panels
    no = DATA.ViewPanels(p);
    if SET(no).TSize > 1
      if DATA.Record
        t = 1+mod(t,maxt);
        currentf = t;
      else
        if usefps
          fps = SET(currentno).FPS;
          % Calculate time elapsed since starttime in seconds
          elapsedtime = rem(now-starttime,1)*24*3600;
          % Calculate how many frames elapsed
          frameselapsed = floor(elapsedtime * fps);
          % Calculate current frame based on fps
          currentf = 1+mod(frameselapsed+startframes,maxt);
        else
          % original implementation using beattime
          t = 1+mod(floor(rem(now-starttime,1)*24*3600/(beattime/maxt)+startframes),maxt);
          currentf = round(SET(no).TSize*(t/maxt));
        end        
      end
      % Ensure current timeframe is inside [1 TSize] limits
      SET(no).CurrentTimeFrame = max(min(currentf,SET(no).TSize),1);
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
function play_Callback2(panels) 
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
  undent(DATA.Handles.playiconholder,'play',0)
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
  if ~findindented(DATA.Handles.playiconholder,'play')
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
function segmentclearalllv_Callback 
%-------------------------------------------
%Clear all LV segmentation, both endo and epi

global DATA

if nargin==1
  if ~yesno('Do you really want to remove all LV segmentation ?',[],DATA.GUI.Segment)

    return;
  end
end

segmentation('clearalllv_Callback');

%-------------------------------------------
function segmentclearallrv_Callback 
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
function segmentclearalllvbutsystolediastole_Callback 
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
arg = struct('endo',true,'epi',true,'rvendo',false,'rvepi',false, 'la', false, 'ra', false);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind, 'laind', ind, 'raind', ind);
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
function segmentclearallrvbutsystolediastole_Callback 
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
arg = struct('endo',false,'epi',false,'rvendo',true,'rvepi',true, 'la', false, 'ra', false);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind, 'laind', ind, 'raind', ind);
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
function contextloadno_Callback(viewpanelstype) 
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
function pointdeletethis_Callback 
%-----------------------------------------
%remove point closest to click

global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if contains(DATA.CurrentTheme,'3dp')
  pointind = annotationpoint('findpoint_helper');
else
  pointind = DATA.LastObject(1);
  %slice = DATA.LastObject(2);
  %tf = DATA.LastObject(3);
end

SET(no).Point.X(pointind) = [];
SET(no).Point.Y(pointind) = [];
SET(no).Point.T(pointind) = [];
SET(no).Point.Z(pointind) = [];
SET(no).Point.Label(pointind) = [];

drawfunctions('drawno',no);

%-----------------------------------------
function pointrenamethis_Callback 
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
    'MV plane',...
    'TV plane',...
    'P1',...
    'P2',...
    'General',...
    'Sector start',...
    'User defined ...'};
  
  % s = myinputdlg({'Enter name'},'Name',1,{sprintf('%s',SET(no).Point.Label{pointind})});  
%   s = mymenu(strcat(dprintf('Select a new name for the point')),menuitems,DATA.GUI.Segment);
  n = 1;
  fieldlabel = 'Point';
  f(n).Field = fieldlabel;
  labelstr = dprintf('Select a new name for the point');
  labelstr = mysplitstring(labelstr);
  % setup string with newline so the label in myinputstruct is split in 2
  % lines
  f(n).Label = [labelstr{1},newline,labelstr{2}];
  f(n).Default = menuitems;
  [outs,ok] = myinputstruct(f,dprintf('Rename point'),10);
  if ok
    s = outs.(fieldlabel);
  else
    % cancel was clicked
    return
  end
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
function pointclearall_Callback 
%-----------------------------------------
%remove all points

global SET NO

no = NO;%DATA.ViewPanels(DATA.CurrentPanel);

SET(no).Point.X = [];
SET(no).Point.Y = [];
SET(no).Point.T = [];
SET(no).Point.Z = [];
SET(no).Point.Label = {};

drawfunctions('drawno',no);

%-------------------------------
function pointclearallname(name,no) 
%-------------------------------
%Function to delete points of a specific color = name

global DATA SET NO

if nargin < 2
  no = NO;
end

n = length(SET(no).Point.X);
if isequal(n,0)
  return
end

inds = false(1,n);

for loop = 1:n
  if isequal(SET(no).Point.Label{loop},name)
    inds(loop) = true;
  end
end

if contains(DATA.ProgramName,'3DPrint')
  segment3dp.pointtools('deletepointsmenu_helper',[],inds) %[] is unused...
else
  for indloop = 1:length(inds)
    SET(no).Point.X(inds(indloop)) = [];
    SET(no).Point.Y(inds(indloop)) = [];
    SET(no).Point.T(inds(indloop)) = [];
    SET(no).Point.Z(inds(indloop)) = [];
    SET(no).Point.Label(inds(indloop)) = [];
  end
end

%-----------------------------------------
function pointmaketimeresolvedthis_Callback 
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
function pointshowthisframeonly_Callback 
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
function interpdeletepoint 
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
opencontour = false;
datanumpoints = tools('getnumpointsforno',no);
[x,y] = calcfunctions('resamplecurve',X,Y,datanumpoints-1,opencontour);
SET(no).([pointtype(1:end-6),'Y'])(:,tf,slice)=[y,y(1)];
SET(no).([pointtype(1:end-6),'X'])(:,tf,slice)=[x,x(1)];

drawfunctions('drawno',no);
if isequal(tf,SET(no).EDT) && ~isempty(SET(no).StrainMitt)
  updateredo(SET(no).StrainMitt,pointtype);
end


%-----------------------------------------
function interpdeletepointthisslicephase 
%-----------------------------------------
%remove interpolation point closest to click

global DATA SET

no = DATA.ViewPanels(DATA.CurrentPanel);

if noobjectshelper
  return
end

%pointind = DATA.LastObject(1);
slice = DATA.LastObject(2);
tf = DATA.LastObject(3);
pointtype = DATA.LastObjectType;

SET(no).([pointtype,'X']){tf,slice} = [];
SET(no).([pointtype,'Y']){tf,slice} = [];

drawfunctions('drawno',no);

%-------------------------------------
function measureclearthis_Callback(~)  
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
function measureclearall_Callback 
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
function measurerenamethis_Callback(measureind) 
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
%   measureind = measureaskhelper;
end

if isequal(measureind,0)
  return;
end

tools('enableundo',no);

%Get new name
[stri,lstr] = measureasklabel(measureind);
if ~isempty(stri)
  SET(no).Measure(measureind).Name = stri;
  SET(no).Measure(measureind).LongName = lstr;
else
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

%----------------------------------
function imagecomment_Callback
%----------------------------------
%Function to add image comment to SET(NO)
addcommenttostack;
viewfunctions('setview'); %refresh

%----------------------------------
function clearimagecomment_Callback
%----------------------------------
%Function to clear all image comments in SET(NO)
global SET NO
str1 = dprintf('This will remove all image comments in the stack.');
str2 = dprintf('Are you sure?');
str = sprintf('%s %s', str1,str2);
if ~yesno(str)
  return;
end

SET(NO).Comment = [];
viewfunctions('setview'); %refresh

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
function clearallslices_Callback 
%-----------------------------------------
%Clear all segmentation, both endo and epi, lv and rv mar and scar.

segmentation('clearslices_Callback');

%-----------------------------------------
function segmentclearall_Callback(silent)  
%-----------------------------------------
%Clear all segmentation, both endo and epi, lv and rv mar and scar.

global DATA

if nargin == 0 || not(silent)
  msg = dprintf('This removes all existing segmentation of LV, RV, LA, RA, ROI, MaR and scar. Are you sure?');
  if ~yesno(msg,[],DATA.GUI.Segment)
    return;
  end
end

viability('viabilityclear_Callback');
roi('roiclearall_Callback')
segmentation('clearall_Callback');
mar('clearall_Callback');
generalpen.atriumpenfunctions('deleteallobjects');

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
arg = struct('endo',true,'epi',true,'rvendo',true,'rvepi',true, 'la', true, 'ra', true);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind, 'laind', ind, 'raind', ind);
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

types = {'LA','RA'};
for loop = 1:length(types)
  if (~isempty(SET(NO).(types{loop})) && ~isempty(SET(NO).(types{loop}).X))
    SET(NO).(types{loop}).X(:,ind,:) = NaN;
    SET(NO).(types{loop}).Y(:,ind,:) = NaN;
  end
end

segment('updatevolume');
drawfunctions('drawcontours',DATA.CurrentPanel);
drawfunctions('drawcontourslara',DATA.CurrentPanel);
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
arg = struct('endo',true,'epi',true,'rvendo',true,'rvepi',true, 'la', true, 'ra', true);
indarg = struct('endoind',ind,'epiind',ind,'rvendoind',ind,'rvepiind',ind, 'laind', ind, 'raind', ind);
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

types = {'LA','RA'};
for loop = 1:length(types)
  if (~isempty(SET(NO).(types{loop})) && ~isempty(SET(NO).(types{loop}).X))
    SET(NO).(types{loop}).X(:,ind,:) = NaN;
    SET(NO).(types{loop}).Y(:,ind,:) = NaN;
  end
end

segment('updatevolume');
drawfunctions('drawcontours',DATA.CurrentPanel);
drawfunctions('drawcontourslara',DATA.CurrentPanel);
drawfunctions('drawinterp',DATA.CurrentPanel);

%------------------------------------
function result = updateintersections_Callback(slices,no,type) 
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
    segstring{i} = sprintf('%s %d - %s',dprintf('Image stack'),segno(i),SET(segno(i)).ImageType); 
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
