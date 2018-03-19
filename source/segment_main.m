function varargout = segment_main(varargin)
% SEGMENT_MAIN Main file for cardiac image analysis software

% Einar Heiberg

% Revision history
% Written by Einar Heiberg, spring/autumn 2002.
% Continously improved ever since 2002-2016.

% The documentation and version history is found in the file changelog.m

%%%%%%%%%%%%%%%%%%
%%%% Main body %%%
fileinput = false;
if nargin > 0
  %Check if input appears to be a filename
  if regexp(varargin{1}, '^[c-zC-Z]:\\')
    fileinput = true;
  end
end

if nargin == 0 || fileinput  % LAUNCH and initalize GUI
  
  % Fix path
  if isdeployed() && isequal(mexext(), 'mexmaci64')
    cd(getapppath());
  end
  
  %Check if os is supported
  arch = mexext();  
  switch arch
    case {'mexglx','mexmaci','mexmaci64'}
      myfailed('Your platform is not supported. Supported platforms are Windows and Linux 64 bit.');
      return;
  end;
  
  programversion = changelog;
  fig = initializesegment(programversion); %Program version number

  varargout = cell(1,nargout);
  if nargout>0
    varargout{1} = fig;
  end
  
  % Load the mat file
  if fileinput
    openfile('loadfiles', {varargin{1}}, false, []); %#ok<CCAT1>
  end
  
else
  %%%% main clause %%%
  macro_helper(varargin{:}); %future macro recording use
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
end
%%%%%%%%%%%%%%%

%-----------------------------------------------
function fig = initializesegment(programversion)
%-----------------------------------------------
%Initialization of Segment GUI
global DATA 

if isa(DATA,'maingui')
  try
    mydisp('Already running.');
    if ~isempty(DATA.fig)
      figure(DATA.fig);
      fig = DATA.fig;
      return;
    else
      fig = [];
    end;
  catch me
    disp('Program not aborted properly last time, all data will be lost and program restarted.');
    mydispexception(me);
  end;
end;

%--- Find source
if isdeployed()
  %Compiled version 
  disp('Standalone');
else
  %Check if platform is supported.
  ext = mexext;
  all = mexext('all');
  arch = '';
  for loop = 1:length(all)
    if isequal(ext,all(loop).ext)
      arch = all(loop).arch;
    end;
  end;
  
  switch arch
    case {'mac','sol64','glnx86'}
      myfailed('Platform is currently not supported. Mex-files are missing.');
      return;
    otherwise
      disp(['Running from Matlab on platform ' arch '.']);
  end;
           
  %Do nothing
end;

%Make sure fresh start
DATA = []; %#ok<NASGU>
SET = []; %#ok<NASGU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This is where we create object%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA = segmentgui(programversion);%Load standard Segment GUI
DATA.init;
fig = DATA.fig;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------
function resetpreview %#ok<DEFNU>
%---------------------
%Reset preview structure in DATA.Preview

global DATA;

DATA.Preview = genemptypreview(DATA.Pref.datapath);


  %-----------------------------
  function framemode_Callback(type)
  %-----------------------------
  %If type 1 we have pressed the single frame mode button type 2 the multi frame mode button 
  global DATA
  singlebutton=DATA.Handles.singleframemodepushbutton;
  multibutton=DATA.Handles.multiframemodepushbutton;
  switch type
    case 1
      set(singlebutton, 'backgroundcolor','g');
      set(multibutton, 'backgroundcolor',[0.8 0.8 0.8]);
      DATA.thisframeonly_Callback(1)
    case  2
      set(singlebutton, 'backgroundcolor',[0.8 0.8 0.8]);
      set(multibutton, 'backgroundcolor','g');
      DATA.thisframeonly_Callback(0)
  end
  
  %---------------------------------
function preview = genemptypreview(datapath)
%---------------------------------
%Generate an empty preview struct
global DATA

roisize = 150;

try %inside try catch if DATA is not set.
  if isempty(DATA.Preview) || isempty(DATA.Preview.ROISize)
    roisize = 150;
  else
    roisize = DATA.Preview.ROISize;
  end
catch
end;

preview = [];
preview.PathName = datapath;
preview.SelectType = 'AUTO';
preview.SliceThickness = 0;
preview.SliceGap = 0;
preview.ResolutionX = 0;
preview.ResolutionY = 0;
preview.TIncr = 0; %s
preview.TimeVector = 0;
preview.TDelay = 0; %for future use to shift datasets in time
preview.EchoTime = 0;
preview.RepetitionTime = 0;
preview.InversionTime = 0;
preview.FlipAngle = 0;
preview.AccessionNumber = '';
preview.StudyUID = '';
preview.StudyID = '';
preview.NumberOfAverages = 0;
preview.VENC = 0;
preview.GEVENCSCALE = 0;
preview.ImagingTechnique = 'MRSSFP';
preview.ImageType = 'General';
preview.ImageViewPlane = 'Unspecified';
preview.ImagePosition = [0 0 0];
preview.ImageOrientation = [1 0 0 0 1 0];
preview.MultiDataSet = false;
preview.Modality = 'MR';
preview.ItemsSelected = 0;
preview.Cyclic = true;
preview.Rotated = false;
preview.ROISize = roisize;
preview.XMin = 1;
preview.YMin = 1;
preview.XSize = 1;
preview.YSize = 1;
preview.PreviewFile = '';
preview.Stable = false;
preview.Scanner = '';
preview.NoNormalize = false;
preview.NSlices = 1; %Number of Slices per DICOM image.
preview.NFrames = 1; %Number of Frames per DICOM image.
preview.NumSlices = 1; %Number of Slices in final image stack
preview.NumFrames = 1; %Number of Frames in final image stack
preview.VENCDirSkipped = 0;
preview.FileType = '';
preview.Silent = false;
preview.SequenceName = '';         % JU comment:
preview.SeriesDescription = '';    %  new since merge
preview.AcquisitionTime = '';      % added by JS
preview.SeriesNumber = '';         %  new since merge
preview.DICOMImageType = '';       %  new since merge
preview.LoadAll=false;             %  new since merge

%-----------------
function initmenu %#ok<DEFNU>
%-----------------
%Initalize the main menu, for instance adds extra utilities, and plugins.
global DATA

%--- Initialize utilities
utility('init');

%--- Initialize plug-ins
try
  load('plugins.mat');
catch  %#ok<CTCH>
  disp('Could not read plugin file.');
  pluginfiles = {};
end

if isdeployed()
  if ~exist('pluginfiles','var')
    myfailed('Problems reading file containing plugins.');
    return;
  end;
else
  %Use filenames
  f = dir([DATA.SegmentFolder filesep 'plugin_*.m']);
  pluginfiles = cell(1,length(f));
  for loop=1:length(f)
    pluginfiles{loop} = f(loop).name;
  end;
end;

%Initialize menus
for loop=1:length(pluginfiles) 
  stri = pluginfiles{loop};
  stri = stri(1:(end-2)); %Remove .m
  handle = uimenu(DATA.Handles.pluginmenu,...
    'Label','temp',...
    'Callback','');
  namestri = feval(stri,'getname',handle);
  set(handle,'Label',namestri);
end;

%------------------------------------------------------------
function addviewicon_helper(callback,tooltip,cdata,tag,separator)
%------------------------------------------------------------
%Helper fcn to add an icon
global DATA 

if nargin<4
  myfailed('Too few input arguments.', DATA.GUI.Segment);
  return;
end;

if nargin<5
  separator = 'off';
end;

props = [];
props.ClickedCallback = callback;
props.ToolTip = tooltip;
props.Tag = tag;
props.CData = cdata;
props.Separator = separator;
DATA.Handles.(tag) = uitoggletool(DATA.Handles.viewtoolbar,props);

%-------------------------------------
function iconhandles = initviewtoolbar %#ok<DEFNU>
%-------------------------------------
%Initalize view toolbar
global DATA 

if nargout > 0
  %Return names of handles, for use in killhandles function
  %If new handles are added, add those to this list as well!
  iconhandles = {'viewtoolbar'
    'previousallframeicon'
    'previousframeicon'
    'playmovieicon'
    'playallicon'
    'nextframeicon'
    'nextallframeicon'
    'fasterframerateicon'
    'slowerframerateicon'
    'cinetoolicon'
    'viewoneicon'
    'mmodeviewicon'
    'viewallicon'
    'montagerowicon'
    'montagefiticon'
    'orthoviewicon'
    'viewzoominicon'
    'viewzoomouticon'
    'refreshicon'
    'resetlighticon'
    'autocontrasticon'
    'undosegmentationicon'
    'imageinfoicon'
    'patientinfoicon'
    'reportsheeticon'
    'movierecordericon'
    'screenshoticon'
    'generalsegmenticon'
    'mpricon'
    'fusionicon'
    'model3icon'
    'volrendicon'
    'reportpersliceicon'
    'reportbullseyeicon'
    'reportlongaxisicon'
    'reportflowicon'};
  return;
end

%Create menubar
DATA.Handles.viewtoolbar = uitoolbar(DATA.imagefig);

%previous all frame icon
addviewicon_helper('segment(''previousallframe_Callback'')',...
  'Previous frame (all image stacks)',...
  DATA.Icons.prevall,...
  'previousallframeicon');

%previous frame icon
addviewicon_helper('segment(''previousframe_Callback'')',...
  'Previous frame',...
  DATA.Icons.prev,...
  'previousframeicon');

%Play icon
addviewicon_helper('segment(''playmovie_Callback'')',...
  'Play current image stack',...
  DATA.Icons.play,...
  'playmovieicon');

%Playall icon
addviewicon_helper('segment(''playall_Callback'')',...
  'Play all visible image stacks',...
  DATA.Icons.playall,...
'playallicon');

%Next frame icon
addviewicon_helper('segment(''nextframe_Callback'')',...
  'Next frame',...
  DATA.Icons.next,...
  'nextframeicon');

%Next all frame icon
addviewicon_helper('segment(''nextallframe_Callback'')',...
  'Next frame (all image stacks)',...
  DATA.Icons.nextall,...
'nextallframeicon');

%Faster framerate icon
addviewicon_helper('segment(''fasterframerate_Callback'')',...
  'Faster frame rate',...
  DATA.Icons.faster,...
  'fasterframerateicon');

%Slower framerate icon
addviewicon_helper('segment(''slowerframerate_Callback'')',...
  'Slower frame rate',...
  DATA.Icons.slower,...
  'slowerframerateicon');

%Slower framerate icon
addviewicon_helper('segment(''cinetool_Callback'')',...
  'Make thumbnail cine loop (hot key (c) )',...
  DATA.Icons.cinetool,...
  'cinetoolicon');

%view one slice
addviewicon_helper('segment(''viewimage_Callback'',''one'')',...
  'View one slice. Hot key to toogle between one slice and montage is (v)',...
  DATA.Icons.viewone,...
  'viewoneicon',...
  'on');

%mmodeview icon
addviewicon_helper('segment(''viewimage_Callback'',''mmode'')',...
'Mmode view',...
DATA.Icons.mmodeview,...
'mmodeviewicon');

%View all icon
addviewicon_helper('segment(''viewimage_Callback'',''montage'')',...
  'View all slices. Hot key to toogle between one slice and montage is (v)',...
  DATA.Icons.viewall,...
  'viewallicon');

%View montage row
addviewicon_helper('segment(''viewimage_Callback'',''montagerow'')',...
  'View all slices by rows',...
  DATA.Icons.montagerow,...
  'montagerowicon');

%View montage fit
addviewicon_helper('segment(''viewimage_Callback'',''montagefit'')',...
  'View all slices in one row/column',...
  DATA.Icons.montagefit,...
  'montagefiticon');

%View orthogonally
addviewicon_helper('segment(''orthoview'')',...
  'View image in orthogonal view planes',...
  DATA.Icons.orthoview,...
  'orthoviewicon');

%zoom in icon
addviewicon_helper('segment(''viewzoomin_Callback'')',...
'Zoom in',...
DATA.Icons.zoomin,...
'viewzoominicon');

%zoom out icon
addviewicon_helper('segment(''viewzoomout_Callback'')',...
  'Zoom out',...
  DATA.Icons.zoomout,...
  'viewzoomouticon');

%refresh
addviewicon_helper('segment(''viewrefreshall_Callback'')',...
  'Refresh image',...
  DATA.Icons.refresh,...
  'refreshicon');

%resetlight
addviewicon_helper('segment(''resetlight_Callback'')',...
  'Reset light/contrast',...
  DATA.Icons.resetlight,...
  'resetlighticon');

%autocontrast
addviewicon_helper('segment(''autocontrast_Callback'')',...
  'Auto contrast',...
  DATA.Icons.autocontrast,...
  'autocontrasticon');

%undosegmentationicon
addviewicon_helper('tools(''undosegmentation_Callback'')',...
  'Undo latest segmentation, and pin placement(s)',...
  DATA.Icons.undo,...
  'undosegmentationicon');

%image information
addviewicon_helper('tools(''imageinfo_Callback'')',...
  'Show image stack details',...
  DATA.Icons.imageinfo,...
  'imageinfoicon',...
  'on');

%patient info
addviewicon_helper('tools(''viewpatientinfo_Callback'')',...
  'View patient details',...
  DATA.Icons.patientinfo,...
  'patientinfoicon');

%report sheet
addviewicon_helper('reporter.reportsheet',...
  'Create text report sheet',...
  DATA.Icons.reportsheet,...
  'reportsheeticon');

%screenshot
addviewicon_helper('export(''screenshot_Callback'')',...
  'Save screenshot',...
  DATA.Icons.screenshot,...
  'screenshoticon');

%movierecorder
addviewicon_helper('export(''exportmovierecorder_Callback'')',...
  'Movie recorder',...
  DATA.Icons.movie,...
  'movierecordericon');

%general segmentation tool
addviewicon_helper('levelset',...
'General segmentation tool and orthogonal view',...
DATA.Icons.generalsegment,...
'generalsegmenticon');

%mpr
addviewicon_helper('reformater',...
  'MPR (Multiple planar reconstruction)',...
  DATA.Icons.mpr,...
  'mpricon');

%fusion
addviewicon_helper('fusion',...
  'Image fusion',...
  DATA.Icons.fusion,...
  'fusionicon');

%3dmodel icon
addviewicon_helper('report3dmodel',...
  'View 3d model',...
  DATA.Icons.model3,...
  'model3icon');

%volumerendering icon
addviewicon_helper('volumerender',...
  'Volume rendering',...
  DATA.Icons.volrend,...
  'volrendicon');

%reportperslice
addviewicon_helper('slicereport',...
  'Report wall motion per slice',...
  DATA.Icons.reportperslice,...
  'reportpersliceicon');

%bullseyeicon
addviewicon_helper('reportbullseye',...
'Bulls eye plot',...
DATA.Icons.bullseye,...
'reportbullseyeicon');

% %Longaxis report
% addviewicon_helper('longaxisplot',...
% 'Long axis sector report',...
% DATA.Icons.reportlongaxis,...
% 'reportlongaxisicon');

%report flow icon
addviewicon_helper('reportflow',...
'Plot flow (Ctrl-T)',...
DATA.Icons.plotflow,...
'reportflowicon');

%--------------------------
function cinetool_Callback
%--------------------------
%Starts the cinetool that allows simultanues segmentation at the
%same time as it plays.

global DATA SET NO

if SET(NO).TSize==1
  stateandicon=segment('iconson','cineplay');
  stateandicon{2}.isindented=0;
  stateandicon{2}.cdataDisplay=stateandicon{2}.cdata;
  DATA.Handles.configiconholder.render
  return;
end

if isa(DATA.CineTimer,'timer')
  cinewindow('update','kill');
  stateandicon=segment('iconson','cineplay');
  stateandicon{2}.isindented=0;
  stateandicon{2}.cdataDisplay=stateandicon{2}.cdata;
  DATA.Handles.configiconholder.render
  %set(DATA.Handles.cinetoolicon,'state','off');
else
  %set(DATA.Handles.cinetoolicon,'state','on');
  stateandicon=segment('iconson','cineplay');
  stateandicon{2}.isindented=1;
  stateandicon{2}.cdataDisplay=stateandicon{2}.cdataIndent;
  DATA.Handles.configiconholder.render
  cinewindow;
end;


%---------------------------------------------
function killbuttondown=switchtopanel(panel,updateimagestack)
%---------------------------------------------
%Make panel the currentpanel
global DATA SET

killbuttondown=0;

if DATA.Silent
  return;
end;

if panel==DATA.CurrentPanel
  return;
end;

if nargin < 2
  updateimagestack = true;
end

%Make panel current panel
DATA.CurrentPanel = panel;

%Make sure the current image stack, runs updateviewicons
if updateimagestack
  if DATA.CurrentPanel>length(DATA.ViewPanels)
    DATA.CurrentPanel=length(DATA.ViewPanels);
  end
  DATA.switchtoimagestack(DATA.ViewPanels(DATA.CurrentPanel));
end
%Make all black, one highlighted
set(DATA.Handles.imageaxes,...
  'xcolor',[0 0 0],'ycolor',[0 0 0],...
  'linewidth',0.5,...
  'visible','off');
set(DATA.Handles.imageaxes(DATA.CurrentPanel),...
  'xcolor',DATA.GUISettings.AxesColor,'ycolor',DATA.GUISettings.AxesColor,...
  'linewidth',2.5,...
  'visible','on');


updateselectedslices
% for loop=1:prod(DATA.ViewMatrix)
%   no=DATA.ViewPanels(loop);
%   set(DATA.Handles.selectslicehandle{loop},'visible','off');
%   switch DATA.ViewPanelsType{loop}
%     case {'montage','montagerow','montagefit','montagesegmented'}
%       if loop~=DATA.CurrentPanel
%         set(DATA.Handles.selectslicehandle{loop}(SET(no).StartSlice:SET(no).EndSlice),'linestyle','--')%'color',[0.7843  0.7843 0.5]);
%       else
%         set(DATA.Handles.selectslicehandle{loop}(SET(no).StartSlice:SET(no).EndSlice),'color',[1 1 0],'linestyle','-');
%       end
%       set(DATA.Handles.selectslicehandle{loop}(SET(no).StartSlice:SET(no).EndSlice),'visible','on');
%     otherwise
%       %nothing for now.
%   end;
% end;

%Update intersection lines
drawfunctions('drawintersections');

%if any drawing tool or crop tool is selcted, change to pointer and do not
%perform any drawing or cropping

switch DATA.CurrentTool
  case {'drawendo','drawepi','drawrvendo','drawrvepi',...
      'interpendo','interpepi','interprvendo','interprvepi',...
      'drawmarpen','drawmarrubberpen','drawmarrubber',...
      'drawroi','putroi','crop'}
    killbuttondown=1;
end

% switch DATA.CurrentTool
%   case {'drawendo','drawepi','drawrvendo','drawrvepi',...
%       'interpendo','interpepi','interprvendo','interprvepi',...
%       'drawmarpen','drawmarrubberpen','drawmarrubber',...
%       'drawroi','putroi','crop'}
%     
%     DATA.guitoggletool('select',DATA.CurrentTool);
%     updatetool('select');
% end

DATA.updateaxestables('measure');    
DATA.updateaxestables('volume');
DATA.updateaxestables('flow');


%----------------------------
function addtopanels(no,mode)
%----------------------------
%Finds an open space, otherwise increases number of panels

global DATA SET
    
%---Find out what view mode to take
if nargin<2
  if not(isempty(SET(no).Flow))
    mode = 'one';
  elseif not(isempty(SET(no).Scar))
    mode = 'montage';
  elseif (SET(no).ZSize>1)&&(SET(no).ZSize<21)&&(SET(no).XSize*SET(no).YSize<1e4)
    mode = 'montage';
  else
    mode = 'one';    
  end;
end;

%--- if space exists then add it
temp = find(DATA.ViewPanels==0);
if ~isempty(temp)
  DATA.ViewPanels(temp(1)) = no;
  DATA.ViewPanelsType{temp(1)} = mode;
  DATA.ViewIM{temp(1)} = [];
  return;
end;

%--- if space does not exist then add
DATA.ViewPanels = [DATA.ViewPanels no];
DATA.ViewPanelsType = cat(2,DATA.ViewPanelsType,{mode});
[rows,cols] = calcfunctions('calcrowscols',no,SET(no).ZSize);
DATA.ViewPanelsMatrix = [DATA.ViewPanelsMatrix {[rows cols]}];
DATA.ViewIM{length(DATA.ViewPanels)} = [];

%---------------------------------------
function viewimagestackas_Callback(mode) %#ok<DEFNU>
%---------------------------------------
%Select how to view current image stack
global NO

if nargin<1
  mode = 'one';
end;

%Get clicked coordinate
%no = getclickedpreview;%not necessary since when clicking preview
%DATA.switchtoimagestack is called

addtopanels(NO,mode)
drawfunctions('drawall');

%---------------------------
function mainresize_Callback %#ok<DEFNU>
%---------------------------
%This fcn is called when user resizes GUI
global DATA

if isempty(DATA)
  %This prevents validate callbacks from reporting error when opening
  %segment.fig for inspection. 
  return;
end;

try
  if isfield(DATA.GUI,'Segment')
    if not(isempty(DATA.GUI.Segment))
      saveguiposition(DATA.GUI.Segment)
    end
  end
catch me
  disp('Could not do mainresize');
  mydispexception(me);
end
  
try
  figunits=get(DATA.fig,'units');
  set(DATA.fig,'units','pixels');
  pfig = get(DATA.fig,'position');
  set(DATA.fig,'units',figunits);
  
  panelunits = get(DATA.Handles.reportpanel,'units');
  set(DATA.Handles.reportpanel,'units','pixels');
  p = get(DATA.Handles.reportpanel,'position');
  
  rpwidth=round(min(DATA.GUISettings.ReportPanelPixelMax,pfig(3)*0.21)); %DATA.GUISettings.RightGapWidth));
  set(DATA.Handles.reportpanel,'position',[...
    pfig(3)-rpwidth ...
    p(2) ...  
    rpwidth ...
    p(4)]);
  DATA.GUISettings.RightGapWidth = rpwidth/pfig(3);
  set(DATA.Handles.reportpanel,'units',panelunits);  
  
%   %Assert position oficon placeholders.
%   %toggle iconholder
%     pos=plotboxpos(DATA.Handles.ribbonaxes);
%     currentpos=get(DATA.Handles.ribbonaxes,'position');
%     set(DATA.Handles.ribbonaxes,'position',currentpos-[pos(1),0,0,0]);
%   
%   %Configholder
%     pos=plotboxpos(DATA.Handles.configaxes);
%     currentpos=get(DATA.Handles.configaxes,'position');
%     set(DATA.Handles.configaxes,'position',currentpos-[pos(1),0,0,0]);
%     
%     %permanentholder
%     pos=plotboxpos(DATA.Handles.permanentaxes);
%     currentpos=get(DATA.Handles.permanentaxes,'position');
%     set(DATA.Handles.permanentaxes,'position',currentpos-[pos(1),0,0,0]);
%     
    %Render iconplaceholders aswell
    DATA.Handles.toggleiconholder.render
    DATA.Handles.permanenticonholder.render
    if ~isempty(DATA.Handles.configiconholder.cdata) 
      DATA.Handles.configiconholder.render
    end
    
    if ~isempty(DATA.Handles.hideiconholder.cdata) 
      DATA.Handles.hideiconholder.render
    end
    
    try
    rows=DATA.ViewMatrix(1);
    cols=DATA.ViewMatrix(2);
    drawfunctions('drawall',rows,cols);
  catch
    DATA.Handles.toggleiconholder.render
    DATA.Handles.permanenticonholder.render
    if ~isempty(DATA.Handles.configiconholder.cdata) 
      DATA.Handles.configiconholder.render
    end
  end
catch me
  if ~isempty(DATA.fig)
    %Mainresize is called uponloading when .fig is not initialized.
    disp('Could not do mainresize');
    mydispexception(me);
  end;
end;

%---------------------------------
function renderstacksfromdicom(no) %#ok<DEFNU>
%---------------------------------
%Render image stacks in main gui. This function is typically called upon
%loading.

%Do not mess with .Silent here since it will be taken care of by lower
%routine images.

global DATA SET

if ~DATA.Silent
  
  %This is an ugly hack to have PC data load two-panel.
  if isempty(DATA.ViewPanels)
    if (length(SET)==2)
      DATA.ViewMatrix=[1 2];
      DATA.ViewPanels=[1 2];
      if (SET(no).ZSize==1)
        DATA.ViewPanelsType{1} = 'one';
        DATA.ViewPanelsType{2} = 'one';
      else
        DATA.ViewPanelsType{1} = DATA.GUISettings.ViewPanelsTypeDefault;
        DATA.ViewPanelsType{2} = DATA.GUISettings.ViewPanelsTypeDefault;
      end
      [rows1,cols1] = calcfunctions('calcrowscols',1,SET(1).ZSize);
      [rows2,cols2] = calcfunctions('calcrowscols',2,SET(2).ZSize);
      DATA.ViewPanelsMatrix = {[rows1 cols1] [rows2 cols2]};
      
      makeviewim(1,1);
      makeviewim(2,2);
    else
      DATA.ViewPanels=1;
      if (SET(no).ZSize==1)
        DATA.ViewPanelsType{1} = 'one';
      else
        DATA.ViewPanelsType{1} = DATA.GUISettings.ViewPanelsTypeDefault;
      end
      [rows,cols] = calcfunctions('calcrowscols',1,SET(1).ZSize);
      DATA.ViewPanelsMatrix = {[rows cols]}; 
      DATA.ViewMatrix=[1 1];
      makeviewim(DATA.CurrentPanel,no);
    end
  end

  if (~DATA.Preview.Silent)
    %Normal loading
    
    DATA.switchtoimagestack(no,true); %force
    drawfunctions('drawthumbnails',isempty(DATA.DATASETPREVIEW));
    drawfunctions('drawall',DATA.ViewMatrix);
    %Ensure correct light settings.
    %resetlight_Callback;
  elseif (length(SET)==1)
    DATA.switchtoimagestack(no,true); %force
    drawfunctions('drawthumbnails',isempty(DATA.DATASETPREVIEW));
    drawfunctions('drawall',1);
  else
    DATA.switchtoimagestack(no,true); %force
    drawfunctions('drawthumbnails',isempty(DATA.DATASETPREVIEW));
    drawfunctions('drawall',DATA.ViewMatrix);
  end;

end;
%Not great position as it triggers for ever stack
%DATA.dataloadedplaceholders
mydisp('Files loaded.');
%endoffcalculation;

%----------------------------
function update_thumbnail(nos)
%----------------------------
%This fcn updates thumbnail no
global DATA SET

%Check if empty DATASETPREVIEW and generate if so.
if isempty(DATA.DATASETPREVIEW)
  calcfunctions('calcdatasetpreview');
end;

for loop = 1:length(nos)
  no = nos(loop);
  
  %Remap
  if isempty(SET(no).Colormap)
    tempim = calcfunctions('remapuint8',...
      SET(no).IM(:,:,round(SET(no).TSize/2),round(SET(no).ZSize/2)),...
      no,calcfunctions('returnmapping',no,true));
  else
    tempim = calcfunctions('remapuint8',...
      SET(no).IM(:,:,round(SET(no).TSize/2),round(SET(no).ZSize/2)),...
      no);
  end

  % zero padding, elegantly done, no? :) /JU
  sz=size(tempim);
  tempim=padarray(tempim,round((length(tempim)-sz(1:2))/2));

  tempim = imresize(tempim,DATA.GUISettings.ThumbnailSize*[1 1],'bilinear');

  %Store, vertically
  DATA.DATASETPREVIEW((no-1)*DATA.GUISettings.ThumbnailSize+(1:DATA.GUISettings.ThumbnailSize),:,:) = tempim;
end;

set(DATA.Handles.datasetpreviewimage,'cdata',DATA.DATASETPREVIEW);

%---------------------------
function out=thumbnailno(in)
%---------------------------
%Helper fcn to remember what image stack were clicked.
persistent no

if nargin==1
  no = in;
end;
out = no;

%-----------------------------
function thumbnail_Buttondown %#ok<DEFNU>
%-----------------------------
%Buttondown fcn for thumbnails.

global DATA 

thumbsize=DATA.GUISettings.ThumbnailSize;

switch get(DATA.fig,'SelectionType')
  case {'extend','normal'}    
    %--- Prepare to drag the thumbnail
    
    %Set up
    set(DATA.fig,'WindowButtonUpFcn',...
      'segment(''thumbnail_Buttonup'')');
    set(DATA.fig,'WindowButtonMotionFcn',...
      'segment(''thumbnail_Motion'')');
    
    %Get clicked position
    [x,y] = mygetcurrentpoint(DATA.Handles.datasetaxes);
    
    %Find clicked image stack
    no = getclickedpreview(x,y);
    thumbnailno(no); %store
    
    %Create axes
    temp = get(DATA.imagefig,'unit');
    set(DATA.imagefig,'unit','pixels');
    try 
      delete(DATA.Handles.thumbnaildragaxes);
    catch %#ok<CTCH>
    end
    DATA.Handles.thumbnaildragaxes = axes(...
      'unit','pixels',...
      'position',...
      [x-32 y-32 64 64],...
      'parent',DATA.imagefig);
    set(DATA.imagefig,'unit',temp);
    
    %Draw image
    try 
      delete(DATA.Handles.thumbnailimage);
    catch %#ok<CTCH>
    end
    DATA.Handles.thumbnailimage=imagesc(...
       DATA.DATASETPREVIEW(thumbsize*(no-1)+(1:thumbsize),:,:),...
      'parent',DATA.Handles.thumbnaildragaxes);
    axis(DATA.Handles.thumbnaildragaxes,'off');
    
  case 'alt'
    %---Right mouse click
    %Get clicked coordinate
    no=getclickedpreview;

    DATA.switchtoimagestack(no);

    %Bring up popup menu
    [p(1),p(2)] = mygetcurrentpoint(DATA.fig);
    set(DATA.Handles.datasetpreviewmenu,...
      'Position',p,...
      'Visible','on');
end;

%------------------------
function thumbnail_Motion %#ok<DEFNU>
%------------------------
%Motion fcn for thumbnails.
global DATA

%Get coordinate
[x,y] = mygetcurrentpoint(DATA.imagefig);

set(DATA.Handles.thumbnaildragaxes,'position',...
  [x-32 y-32 64 64]);

%--------------------------
function thumbnail_Buttonup %#ok<DEFNU>
%--------------------------
%Buttonup fcn for thumbnails.
global DATA SET NO
%Get coordinate
[x,y] = mygetcurrentpoint(DATA.Handles.boxaxes);

if nargin==1
  x=-1;
end
  
%Restore motion etc
set(DATA.imagefig,'WindowButtonMotionFcn','');
set(DATA.imagefig,'WindowButtonUpFcn','segment(''buttonup_Callback'')');
%Hide the image
try
  delete(DATA.Handles.thumbnaildragaxes);
  delete(DATA.Handles.thumbnailimage);
catch %#ok<CTCH>
end

%Retrieve what image stack chosen
no=thumbnailno;

% Only change panel im if pointer has moved out of the sidebar.
if (size(DATA.ViewPanels)==1)
  ind=1;
elseif (x < 0)
  allreadyout=find(DATA.ViewPanels==no);
  if ~isempty(allreadyout)
    switchtopanel(allreadyout(1));
    return;
  end
  ind=find(DATA.ViewPanels==0,1);
else
  %Find at which image panel we drop it.
  dist = zeros(1,length(DATA.Handles.imageaxes));
  for loop=1:length(DATA.Handles.imageaxes)
    p = get(DATA.Handles.imageaxes(loop),'position');
    p = p(1:2)+0.5*p(3:4); %Center position
    dist(loop) = sqrt(sum(([x y]-p).^2)); %Euklidean distance
  end;
  [~,ind] = min(dist);
  ind=ind(1);   %just in case equal distance..!
end

if ~isempty(ind)
  %Store it
  DATA.ViewPanels(ind) = no;
  if (SET(no).ZSize==1)
    DATA.ViewPanelsType{ind} = 'one';
  else
    DATA.ViewPanelsType{ind} = DATA.GUISettings.ViewPanelsTypeDefault;
  end
  DATA.ViewIM{ind} = [];
  DATA.Overlay(ind) = struct('alphadata', [], 'cdata', []);
  DATA.CurrentPanel = ind;

  oldnos=SET(NO).Linked;
  NO = no;

  % This is to ensure that linkaxes is removed when a non-linked stack is
  % moved into a panel that was previously linked. /JU
  panel=find(ismember(DATA.ViewPanels,NO));
  if ~isempty(panel)&&...
      isappdata(DATA.Handles.imageaxes(panel(1)),'graphics_linkaxes')&&...
      ~ismember(NO,oldnos);
      %isempty(SET(NO).Flow)
    linkaxes(DATA.Handles.imageaxes(panel),'off');
  end
  drawfunctions('drawimageno',NO);
  drawfunctions('drawallslices');
  DATA.switchtoimagestack(NO,true); %force
  
end

%-------------------------
function addnotopanel(no) %#ok<DEFNU>
%-------------------------
    global DATA NO SET
    
    ind=find(DATA.ViewPanels==0,1);
    
    if ~isempty(ind)
      DATA.ViewPanels(ind) = no;
      if (SET(no).ZSize==1)
        DATA.ViewPanelsType{ind} = 'one';
      else
        DATA.ViewPanelsType{ind} = DATA.GUISettings.ViewPanelsTypeDefault;
      end
      DATA.ViewIM{ind} = [];
      DATA.Overlay(ind) = struct('alphadata', [], 'cdata', []);
      DATA.CurrentPanel = ind;
      
      oldnos=SET(NO).Linked;
      NO = no;
      
      % This is to ensure that linkaxes is removed when a non-linked stack is
      % moved into a panel that was previously linked. /JU
      panel=find(ismember(DATA.ViewPanels,NO));
      if ~isempty(panel)&&...
          isappdata(DATA.Handles.imageaxes(panel(1)),'graphics_linkaxes')&&...
          ~ismember(NO,oldnos);
        %isempty(SET(NO).Flow)
        linkaxes(DATA.Handles.imageaxes(panel),'off');
      end
      drawfunctions('drawimageno',NO);
      drawfunctions('drawallslices');
      DATA.switchtoimagestack(NO,true); %force
    end

%----------------------------
function updateselectedslices
%----------------------------
%Graphically updates which slices are selected.
global DATA SET NO

for loop=1:prod(DATA.ViewMatrix)
  no=DATA.ViewPanels(loop);
  set(DATA.Handles.selectslicehandle{loop},'visible','off');
  switch DATA.ViewPanelsType{loop}
    case {'montage','montagerow','montagefit','montagesegmented'}
      if loop~=DATA.CurrentPanel
        set(DATA.Handles.selectslicehandle{loop}(SET(no).StartSlice:SET(no).EndSlice),'linestyle','--')%'color',[0.7843  0.7843 0.5]);
      else
        set(DATA.Handles.selectslicehandle{loop}(SET(no).StartSlice:SET(no).EndSlice),'color',[1 1 0],'linestyle','-');
      end
      set(DATA.Handles.selectslicehandle{loop}(SET(no).StartSlice:SET(no).EndSlice),'visible','on');
    otherwise
      %nothing for now.
  end;
end;

%fix so that linked images are connected
nos = SET(NO).Linked;

panelstodo=find(ismember(DATA.ViewPanels,nos));
nos=DATA.ViewPanels(panelstodo);
for loop=1:length(panelstodo)
  panel=panelstodo(loop);
  no=nos(loop);
  set(DATA.Handles.selectslicehandle{panel},'visible','off');
  switch DATA.ViewPanelsType{panel}
    case {'montage','montagerow','montagefit','montagesegmented'}
      set(DATA.Handles.selectslicehandle{panel}(SET(no).StartSlice:SET(no).EndSlice),'visible','on');
    otherwise
      %nothing for now.
  end;
end;

%--------------------------
function z = remap(im,cmap,c,b) %#ok<DEFNU>
%--------------------------
%Remap data according to cmap

global SET NO

if nargin<2
  cmap = SET(NO).Colormap;
end;
if nargin<4
  c=SET(NO).IntensityMapping.Contrast;
  b=SET(NO).IntensityMapping.Brightness;
end

if isempty(cmap)
  z = im;
else
  map = cmap(:,1);
  switch class(im)
    case 'double'
      map = double(map);
      outsize = size(im);
      im=c*im(:)+(b-0.5);
      z = map(max(min(round(im(:)*256),length(map)),1));
      z = reshape(z,outsize);
    case 'single'
% This is turned off because this is where calls from fusion.m are made,
% and it currently has it's own contrast/brightness settings, that have
% already been applied to colormap. Ideally, this would be reorganized so
% that even remap() accepted a local NO argument. /JU
%      im=c*im+(b-0.5);
      outsize = size(im);
      if nargin > 3
        im=single(c*im(:)+(b-0.5));
      end      
      z = fastremap(im,single(map));      
      z = reshape(z,outsize);
  end;
end;

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
end;

result = true;
segno = [];

switch type
  case 'LV Endocardial'
    % Find sets with LV endocardial segmentation
    for i=1:length(SET)
      if ~isempty(SET(i).EndoX) && any(any(any(~isnan(SET(i).EndoX))))
        segno = [segno i]; %#ok<AGROW>
      end;
    end;
  case 'RV Endocardial'
    % Find sets with RV endocardial segmentation
    for i=1:length(SET)
      if ~isempty(SET(i).RVEndoX) && any(any(any(~isnan(SET(i).RVEndoX))))
        segno = [segno i]; %#ok<AGROW>
      end;
    end;
end;

if isempty(segno)
  myfailed('No segmentation found',DATA.GUI.Segment);
  result = false;
  return;
end;

if length(segno)>1
  segstring = cell(length(segno),1);
  for i = 1: length(segno)
    segstring{i} = sprintf('Image stack %d - %s',segno(i),SET(segno(i)).ImageType);
  end;
  
  [selected, ok] = listdlg('ListString',segstring,...
    'SelectionMode','single',...
    'PromptString','Chose image stack to use:',...
    'ListSize',[300 100]);
  
  if ~ok
    result = false;
    return;
  end;
  segno = segno(selected);
end;

calcsegintersect(segno,no,slices,type);

drawfunctions('drawall');

%-----------------------------------------------
function [intersections, maxintersect] = getendointersection(no) %#ok<DEFNU>
%-----------------------------------------------
% Returns the intersection of the endocardial segmentation and the current
% slice and current time frame of SET(no)

global SET

intersections = [];
maxintersect = [];

slice = SET(no).CurrentSlice;
time = SET(no).CurrentTimeFrame;

if(isfield(SET,'Intersection') && ...
    isstruct(SET(no).Intersection))
    
  isectnum = strncmp('LV Endocardial',{SET(no).Intersection.Type},14); 
  
  if ~isempty(isectnum) && ...
      ~(isempty(SET(no).Intersection(isectnum).Slice(slice).TimeFrame))
    
    intersections = SET(no).Intersection(isectnum)...
      .Slice(slice).TimeFrame(time).Intersection;
    maxintersect = SET(no).Intersection(isectnum).MaxIntersect;
  end;
  
end;

%------------------
function updateflow
%------------------
%calculate flow and graphically update
global DATA SET

if ~isempty(DATA.FlowNO)
  no = DATA.FlowNO;
  if ~isempty(SET(no).Flow) && isfield(SET(no).Flow,'MagnitudeNo') && SET(no).RoiN > 0
    calcfunctions('calcflow',no);
  end
end
DATA.updateaxestables('flow');


%----------------------------------
function updatevolume(lvsegchanged)
%----------------------------------
%Calc volume of segmentation and graphically update.
global DATA SET NO

if nargin < 1
  lvsegchanged = false;
end

rotstring = '';
specstring = sprintf('\n');

%Error check
if isempty(SET(NO).ImageViewPlane)
  SET(NO).ImageViewPlane = 'Unspecified';
end;

calcd = -1;
if ismember(SET(NO).ImageViewPlane,{'2CH','3CH','4CH'})
  [calcd,usednos] = longaxistools('calcbiplanevolume');
  if all(isnan(SET(NO).LVV)) || ~ismember(NO,usednos)
    calcd = 0;
  end
end
 DATA.AxesTables.volume.updateName('INFO','',true);
if calcd > 0
  DATA.AxesTables.volume.updateName('INFO',translation.dictionary(sprintf('Volume from %d LAX stacks',calcd)),true);
  calcfunctions('volume_helper',NO);
elseif calcd == -1
  calcfunctions('calcvolume',NO);
else
  calcfunctions('volume_helper',NO);
end

if DATA.Silent
  return;
end;

DATA.updateaxestables('measure');

% %Create report string
% stri = [];
% 
% %LV
% useml=SET(NO).EDV>1|isnan(SET(NO).EDV);
% lvmed = SET(NO).EPV(SET(NO).EDT)-SET(NO).LVV(SET(NO).EDT)+SET(NO).PV(SET(NO).EDT);
% lvmes = SET(NO).EPV(SET(NO).EST)-SET(NO).LVV(SET(NO).EST)+SET(NO).PV(SET(NO).EST);
% lvm = mynanmean([lvmed lvmes]);
% if SET(NO).Rotated
%   rotstring = dprintf(translation.dictionary('Rotated stack\n'));
% end
% 
% if not(isnan(lvm))&&(lvm>0)
%   if useml
%     stri = [stri ...
%       sprintf('%sED/ES LVM: %3.0f/%3.0f',specstring,lvmed,lvmes)];
%     stri = [stri ...
%       sprintf(' [ml]\nLVM: %3.0f [g]\n',1.05*lvm)];
%   else
%     stri = [stri ...
%       sprintf('%sED/ES LVM: %3.0f/%3.0f',specstring,lvmed*1000,lvmes*1000)];
% 		stri = [stri...
%       sprintf(' [ul]\nLVM: %3.0f [mg]\n',1.05*lvm*1000)];
% 	end
% 	
% else
%   stri = [stri ...
%     sprintf('%sLVM:  - [ml]\n',specstring) ...
%     sprintf('LVM:  - [g]\n')];
% end
% 
% %RV
% rvuseml = SET(NO).RVEDV>1|isnan(SET(NO).RVEDV)|SET(NO).RVESV>1|isnan(SET(NO).RVESV);
% rvmed = SET(NO).RVEPV(SET(NO).EDT)-SET(NO).RVV(SET(NO).EDT);
% rvmes = SET(NO).RVEPV(SET(NO).EST)-SET(NO).RVV(SET(NO).EST);
% rvm = mynanmean([rvmed rvmes]);
% if not(isnan(rvm))&&(rvm>0)
%   if rvuseml
%     stri = [stri ...
%       sprintf('ED/ES RVM: %3.0f/%3.0f',rvmed,rvmes)];
%     stri = [stri ...
%       sprintf(' [ml]\nRVM: %3.0f [g]\n',1.05*rvm)];
%   else
%     stri = [stri ...
%       sprintf('ED/ES RVM: %3.0f/%3.0f',rvmed*1000,rvmes*1000)];
%     stri = [stri...
%       sprintf(' [ul]\nRVM: %3.0f [mg]\n',1.05*rvm*1000)];
%   end
% else
%   stri = [stri ...
%     sprintf('%sRVM:  - [ml]\n',rotstring) ...
%     sprintf('RVM:  - [g]\n')];
% end
% 
% if useml
%   if SET(NO).TSize == 1
%     lvstri = [...
%       sprintf('LV-EDV:%3d [ml]\n',round(SET(NO).EDV)) ...
%       sprintf('HR:%3d [bpm]\n',round(SET(NO).HeartRate))];
%   else
%     lvstri = [...
%       sprintf('LV-EDV:%3d [ml]\n',round(SET(NO).EDV)) ...
%       sprintf('LV-ESV:%3d [ml]\n',round(SET(NO).ESV)) ...
%       sprintf('LV-SV:%3d [ml]\n',round(SET(NO).SV)) ...
%       sprintf('LV-EF:%3d [%%]\n',round(SET(NO).EF*100)) ...
%       sprintf('HR:%3d [bpm]\n',round(SET(NO).HeartRate)) ...
%       sprintf('CO:%1.1f [l/min]\n',(SET(NO).HeartRate*SET(NO).SV)/1000) ...
%       sprintf('LV-PFR:%3d [ml/s]\n',round(SET(NO).PFR)) ...
%       sprintf('LV-PER:%3d [ml/s]\n',round(SET(NO).PER))];
%   end
% else
%   if SET(NO).TSize == 1
%     lvstri = [...
%       sprintf('LV-EDV:%3d [ul]\n',round(1000*SET(NO).EDV)) ...
%       sprintf('HR:%3d [bpm]\n',round(SET(NO).HeartRate))];
%   else
%     lvstri = [...
%       sprintf('LV-EDV:%3d [ul]\n',round(1000*SET(NO).EDV)) ...
%       sprintf('LV-ESV:%3d [ul]\n',round(1000*SET(NO).ESV)) ...
%       sprintf('LV-SV:%3d [ul]\n',round(1000*SET(NO).SV)) ...
%       sprintf('LV-EF:%3d [%%]\n',round(SET(NO).EF*100)) ...
%       sprintf('HR:%3d [bpm]\n',round(SET(NO).HeartRate)) ...
%       sprintf('CO:%1.1f [ml/min]\n',(1000*SET(NO).HeartRate*SET(NO).SV)/1000) ...
%       sprintf('LV-PFR:%3d [ul/s]\n',round(1000*SET(NO).PFR)) ...
%       sprintf('LV-PER:%3d [ul/s]\n',round(1000*SET(NO).PER))];
%   end
% end
% 
% if not(isempty(SET(NO).Scar))
%   %Determine what to type.
%   if SET(NO).Scar.UseWeighting %isequal(SET(NO).Scar.Mode,'undocumented')||isequal(SET(NO).Scar.Mode,'weighted')
%     temp = 'Scar-W';
%   else
%     temp = 'Scar';
%   end;
%   volscale = SET(NO).ResolutionX*SET(NO).ResolutionY*(SET(NO).SliceThickness+SET(NO).SliceGap)/1e3;
%   if isfield(SET(NO).Scar,'NadirIndex')
%     lvstri = [lvstri ...
%       sprintf('%s: %3d [%%] %3d [ml]\n','Scar',round(SET(NO).Scar.Percentage),...
%       round((SET(NO).Scar.Percentage/100)*sum(SET(NO).Scar.MyocardMask(:))*volscale)),...
%       sprintf('Severity Index:%3d [%%]\n',round(SET(NO).Scar.SeverityIndex)) ...
%       sprintf('Nadir Index:%3d [%%]\n',round(SET(NO).Scar.NadirIndex)) ...
%       sprintf('TPD:%3d\n',round(SET(NO).Scar.TPD)) ...
%       sprintf('TPD LAD:%3d\n',round(SET(NO).Scar.TPDLAD)) ...
%       sprintf('TPD LCx:%3d\n',round(SET(NO).Scar.TPDLCx)) ...
%       sprintf('TPD RCA:%3d\n',round(SET(NO).Scar.TPDRCA))];
%   else
%     if useml
%       lvstri = [lvstri ...
%         sprintf('%s:%d%% %d [ml]\nMO:%d%% %d [ml]\nMO-ext:%d%% %d [ml]\n',...
%         temp,...
%         round(SET(NO).Scar.Percentage),...
%         round((SET(NO).Scar.Percentage/100)*sum(SET(NO).Scar.MyocardMask(:))*volscale),...
%         round(SET(NO).Scar.MOPercentage),...
%         round((SET(NO).Scar.MOPercentage/100)*sum(SET(NO).Scar.MyocardMask(:))*volscale),...
%         round(100*sum(SET(NO).Scar.NoReflow(:))/sum(SET(NO).Scar.MyocardMask(:))),...
%         round(sum(SET(NO).Scar.NoReflow(:))*volscale))];
%     else
%       lvstri = [lvstri ...
%         sprintf('%s:%d%% %d [ul]\nMO:%d%% %d [ul]\nMO-ext:%d%% %d [ul]\n',...
%         temp,...
%         round(SET(NO).Scar.Percentage),...
%         round((SET(NO).Scar.Percentage/100)*sum(SET(NO).Scar.MyocardMask(:))*volscale*1000),...
%         round(SET(NO).Scar.MOPercentage),...
%         round((SET(NO).Scar.MOPercentage/100)*sum(SET(NO).Scar.MyocardMask(:))*volscale*1000),...
%         round(100*sum(SET(NO).Scar.NoReflow(:))/sum(SET(NO).Scar.MyocardMask(:))),...
%         round(sum(SET(NO).Scar.NoReflow(:))*volscale*1000))];
%     end
%   end 
% end;
% 
% if not(isempty(SET(NO).AtrialScar))
%   temp = 'Atrial Scar';
%   if useml
%     lvstri = [lvstri ...
%       sprintf('%s: %d%%\n',...
%       temp,...
%       round(SET(NO).AtrialScar.Percentage))];
%   else
%     lvstri = [lvstri ...
%       sprintf('%s: %d%%\n',...
%       temp,...
%       round(SET(NO).AtrialScar.Percentage))];
%   end
% end
% 
% if not(isempty(SET(NO).MaR))
%   volscale = SET(NO).ResolutionX*SET(NO).ResolutionY*(SET(NO).SliceThickness+SET(NO).SliceGap)/1e3;
%   
%   if SET(NO).TSize==1
%     if ~isempty(SET(NO).MaR.MPS.SeverityIndex)
%       lvstri = [lvstri ...
%         sprintf('%s: %3d [%%] %3d [ml]\n','MaR:',round(SET(NO).MaR.Percentage),...
%         round((SET(NO).MaR.Percentage/100)*sum(SET(NO).MaR.MyocardMask(:))*volscale)),...
%         sprintf('Severity Index:%3d [%%]\n',round(SET(NO).MaR.MPS.SeverityIndex)) ...
%         sprintf('Nadir Index:%3d [%%]\n',round(SET(NO).MaR.MPS.NadirIndex)) ...
%         sprintf('TPD:%3d\n',round(SET(NO).MaR.MPS.TPD)) ...
%         sprintf('TPD LAD:%3d\n',round(SET(NO).MaR.MPS.TPDLAD)) ...
%         sprintf('TPD LCx:%3d\n',round(SET(NO).MaR.MPS.TPDLCx)) ...
%         sprintf('TPD RCA:%3d\n',round(SET(NO).MaR.MPS.TPDRCA))];
%     else
%       marpercent = round(SET(NO).MaR.Percentage);
%       marml = round((marpercent/100)*sum(SET(NO).MaR.MyocardMask(:))*volscale);
% 
%       lvstri = [lvstri ...
%         sprintf('MaR: %d%% %d[ml]\n',...
%         marpercent,...
%         marml)];
%     end
% 
%   else
%     maredpercent = round(SET(NO).MaR.Percentage(SET(NO).EDT));
%     marespercent = round(SET(NO).MaR.Percentage(SET(NO).EST));
%     maredml = round((maredpercent/100)*sum(sum(sum(SET(NO).MaR.MyocardMask(:,:,SET(NO).EDT,:))))*volscale);
%     maresml = round((marespercent/100)*sum(sum(sum(SET(NO).MaR.MyocardMask(:,:,SET(NO).EST,:))))*volscale);
%     
%     lvstri = [lvstri ...
%       sprintf('MaR(EDT):%d%%\n%d[ml]\nMaR(EST):%d%% %d[ml]\n',...
%       maredpercent,...
%       maredml,...
%       marespercent,...
%       maresml)];
%   end
% end;
% 
% if rvuseml
%   if SET(NO).TSize == 1
%     rvstri = sprintf('RV-EDV:%3d [ml]\n',round(SET(NO).RVEDV));
%   else
%     rvstri = [
%       sprintf('RV-EDV:%3d [ml]\n',round(SET(NO).RVEDV)) ...
%       sprintf('RV-ESV:%3d [ml]\n',round(SET(NO).RVESV)) ...
%       sprintf('RV-SV:%3d [ml]\n',round(SET(NO).RVSV)) ...
%       sprintf('RV-EF:%3d %%\n',round(SET(NO).RVEF*100))];
%   end
% else
%   if SET(NO).TSize == 1
%     rvstri = sprintf('RV-EDV:%3d [ul]\n',round(1000*SET(NO).RVEDV));
%   else
%     rvstri = [
%       sprintf('RV-EDV:%3d [ul]\n',round(1000*SET(NO).RVEDV)) ...
%       sprintf('RV-ESV:%3d [ul]\n',round(1000*SET(NO).RVESV)) ...
%       sprintf('RV-SV:%3d [ul]\n',round(1000*SET(NO).RVSV)) ...
%       sprintf('RV-EF:%3d %%\n',round(SET(NO).RVEF*100))];
%   end
% end
% %update text
% set(DATA.Handles.volumereporttext,...
%   'String',stri,'visible','on');
% set(DATA.Handles.lvvolumereporttext,...
%   'String',lvstri,'visible','on');
% set(DATA.Handles.rvvolumereporttext,...
%   'String',rvstri,'visible','on');

DATA.updateaxestables('volume');

if lvsegchanged
  if ismember('Strain from tagging',{SET.ImageType})
    taggingno = find(strcmp('Strain from tagging',{SET.ImageType}));
    if ismember(NO,taggingno) && ~isempty(SET(NO).StrainTagging)
      SET(NO).StrainTagging.LVupdated = true;
    end
    for tno = taggingno
      if ~isempty(SET(tno).StrainTagging) && isfield(SET(tno).StrainTagging,'cineno') && NO == SET(tno).StrainTagging.cineno
        if isfield(SET(tno).StrainTagging,'importfromcine') && SET(tno).StrainTagging.importfromcine
          SET(tno).StrainTagging.LVupdated = true;
        end
      end
    end
  end
end

%--------------------------------------
function [x,y,slice] = getclickedcoords
%--------------------------------------
%Find coordinates where the user last clicked. x & y are given in internal
%coordinate system, i.e the functions determines slice in montage view.
global DATA SET NO

%Extract coordinates clicked
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(DATA.CurrentPanel));
panel = DATA.CurrentPanel;
type = DATA.ViewPanelsType{panel};

if any(strcmp(type,{'montage','montagerow','montagefit','sax3','montagesegmented'}))%ismember(type,{'montage','montagerow','montagefit','sax3','montagesegmented'})
  %Find slice
  col = 1+floor((x-0.5)/SET(NO).YSize);
  row = 1+floor((y-0.5)/SET(NO).XSize);
  slice = col+(row-1)*DATA.ViewPanelsMatrix{panel}(2);
  
  %Special case for SAX3 view
  if strcmp(type,'sax3') && slice > 0 && slice <= size(SET(NO).SAX3.slices,1)
    slice = SET(NO).SAX3.slices(slice,SET(NO).CurrentTimeFrame);
  elseif strcmp(type,'montagesegmented')
    slicestoinclude = getmontagesegmentedslices(NO);
    slice = slice + slicestoinclude(1)-1;
  end
  
  %Find coordinates within image
  x = x-(col-1)*SET(NO).YSize;
  y = y-(row-1)*SET(NO).XSize;
elseif strcmp(type,'hla')
  slice = SET(NO).HLA.slice;
elseif strcmp(type,'vla')
  slice = SET(NO).VLA.slice;
elseif strcmp(type,'gla')
  slice = 0;
else
  %set slice
  slice = SET(NO).CurrentSlice;
end;

%--------------------------------
function no=getclickedpreview(~,y)
%--------------------------------
%function which returns the clicked preview image
global DATA

if nargin==0
  [~,y] = mygetcurrentpoint(DATA.Handles.datasetaxes);
end

no = floor(y/DATA.GUISettings.ThumbnailSize)+1;

%----------------------------
function switchtoslice(slice)
%----------------------------
%Called when current slice has changed.
%If slice has changed, make sure montage/one are in sync
%Does nothing if slice = CurrentSlice
global DATA SET NO

if ~isequal(slice,SET(NO).CurrentSlice)||...
    ~isequal(SET(NO).StartSlice,SET(NO).EndSlice)
  %Update linked stacks as well
  nos = SET(NO).Linked;
  for no=unique(nos)
    %copy to linked flows
    SET(no).StartSlice=slice;
    SET(no).CurrentSlice=slice;
    SET(no).EndSlice=slice;
    if ismember(no,DATA.ViewPanels)
      updateoneim(no);
      centeronslice(slice,no,SET(no).ZSize);
      drawfunctions('drawsliceno',no);
    end
  end
  drawfunctions('drawintersections');
  updateselectedslices;
end

%----------------------------------------------------
function switchtolongaxisslice(slice,laxfield,silent)
%----------------------------------------------------
%Called when slice has changed in HLA or VLA view
global SET NO

if nargin < 3
  silent = false;
end
if slice ~= SET(NO).(laxfield).slice
  %Update linked stacks as well
  nos = SET(NO).Linked;
  for no=unique(nos)
    %copy to linked flows
    SET(no).(laxfield).slice=slice;
    if ~silent
      updateoneim(no);
    end
  end
end

if silent
  return
end
centeronslice(slice,NO,SET(NO).ZSize);
drawfunctions('drawsliceno',NO);
drawfunctions('drawintersections');
updateselectedslices;
    
%-----------------------------------
function centeronslice(slice,no,zsz)
%-----------------------------------
%Put slice in center of montagefit view
global DATA
for panel = find(DATA.ViewPanels == no)
  
  if strcmp(DATA.ViewPanelsType{panel},'montagefit')
    viewimsz = size(DATA.ViewIM{panel});
    if viewimsz(1) < viewimsz(2)
      limname = 'XLim';
      limmax = viewimsz(2);
    else
      limname = 'YLim';
      limmax = viewimsz(1);
    end
    slicewidth = limmax/zsz;
    slicemid = slicewidth*(slice-1)+slicewidth/2;
    lims = get(DATA.Handles.imageaxes(panel),limname);
    newlims = slicemid + [-0.5 0.5]*diff(lims);
    set(DATA.Handles.imageaxes(panel),limname,newlims);
    %updatemodeldisplay(no);
  end
end

%-----------------------------
function updatemodeldisplay(~)
%-----------------------------
%Do nothing, introduced to disable excessive calls to updatemodeldisplay

%--------------------
function r = getfieldifcommon(SET, fname) %#ok<DEFNU>
%--------------------
%Helper function  to filesavedicom_Callback
if numel(SET) == 0
  r = [];
  return
end

r = SET(1).(fname);
for n=2:numel(SET)
  if not(isequalwithequalnans(r, SET(n).(fname)))
    r = [];
    return;
  end
end

%--------------------------------
function updateoneim(no)
%--------------------------------
%Updates viewim for 'one view' or 'mmodespatial' 
%Called when changed currentslice
global DATA SET

nos = SET(no).Linked;

%Call make view im if necessary
%This is time consuming, though the alternative will consume memory
for loop=1:length(DATA.ViewPanels)
%  if ...
  if (ismember(DATA.ViewPanels(loop),nos))&&... 
  ismember(DATA.ViewPanelsType{loop},{'one','mmodespatial','ortho','hla','vla','gla'})
    makeviewim(loop,DATA.ViewPanels(loop));
  end;

end;

if ismember('mmodetemporal',DATA.ViewPanelsType)
  updatemmode;
end

%--------------------------------
function movealltowardsbase_Callback
%--------------------------------
% Changes CurrentSlice (and that of parallel image stacks)
% towards base

% Marten Larsson, June, 3, 2009

global SET NO

% Move SET(NO) towards base
if SET(NO).CurrentSlice>1
  slice = SET(NO).CurrentSlice-1;
  SET(NO).CurrentSlice = slice;
  SET(NO).StartSlice = slice;
  SET(NO).EndSlice  = slice;

  updateoneim(NO);
  drawfunctions('drawsliceno',NO);
end;

updateparallelsets;

if length(SET(NO).Linked) > 1
  nos = SET(NO).Linked;

  for no=nos
    SET(no).CurrentSlice = SET(NO).CurrentSlice;
    SET(no).StartSlice = SET(NO).CurrentSlice;
    SET(no).EndSlice  = SET(NO).CurrentSlice;
  end
end

drawfunctions('drawintersections');
updateselectedslices;

%--------------------------------
function movealltowardsapex_Callback
%--------------------------------
% Changes CurrentSlice (and that of parallel image stacks)
% towards apex

% Marten Larsson, June, 3, 2009

global SET NO

% Move SET(NO) towards base
if SET(NO).CurrentSlice < SET(NO).ZSize
  slice = SET(NO).CurrentSlice+1;
  SET(NO).CurrentSlice = slice;
  SET(NO).StartSlice = slice;
  SET(NO).EndSlice  = slice;

  updateoneim(NO);
  drawfunctions('drawsliceno',NO);
end;

updateparallelsets;

if length(SET(NO).Linked) > 1
  nos = SET(NO).Linked;

  for no=nos
    SET(no).CurrentSlice = SET(NO).CurrentSlice;
    SET(no).StartSlice = SET(NO).CurrentSlice;
    SET(no).EndSlice  = SET(NO).CurrentSlice;
  end
end

drawfunctions('drawintersections');
updateselectedslices;

%--------------------------
function updateparallelsets
%--------------------------
% Chages CurrentSlice in SETs that are parallel to SET(NO) to the slice
% closest to SET(NO).CurrentSlice
%
% Help function to movealltowardsbase_Callback and
% movealltowardsapex_Callback

% Marten Larsson, June, 3, 2009

global DATA SET NO

% Get open SETs panels but exclude NO
openpanels = DATA.ViewPanels;
openpanels = unique(openpanels);
openpanels = openpanels(openpanels~=NO&openpanels~=0);

% Find SETs that are parallel to SET(NO), but not a linked flow.
if ~isempty(SET(NO).Flow)
  flownos = [...
    SET(NO).Flow.MagnitudeNo ...
    SET(NO).Flow.PhaseNo ...
    SET(NO).Flow.PhaseX ...
    SET(NO).Flow.PhaseY ...
    SET(NO).Flow.Angio ...
    SET(NO).Flow.VelMag];
else
  flownos = [];
end
parallel  = [];
for i = openpanels
  sameview = orientationcomparison(NO,i);
  linkedflow = any(flownos == i);  % SET(i) is linked to SET(NO) via SET(NO).Flow.
                                   % This case is handled by updateoneim().
  if sameview && ~linkedflow
    parallel = [parallel i]; %#ok<AGROW>
  end;
end;

if not(isempty(parallel))
  % Calculate new CurrentSlide for parallel SETs
  viewdir = cross(SET(NO).ImageOrientation(1:3),SET(NO).ImageOrientation(4:6))';
  currzdist = (SET(NO).CurrentSlice-1)*(SET(NO).SliceThickness+SET(NO).SliceGap);

  % DEBUG
  % disp(sprintf('SET(NO).CurrentSlice z-distance: %d mm', currzdist));

  oldNO = NO;
  oldlink=SET(NO).Linked;
  domontageupdate=0;
  for i = 1:length(parallel)
    % Calculate closest slice in parallel SETs and change their
    % CurrentSlice
    slicethickness = SET(parallel(i)).SliceThickness + SET(parallel(i)).SliceGap;
    zdistances = 0:slicethickness:(slicethickness * (SET(parallel(i)).ZSize-1));
    zdiff = SET(parallel(i)).ImagePosition - SET(NO).ImagePosition;
    zdiff = dot(zdiff,viewdir);
    zdistances =  zdistances - zdiff;
    [~,slice] = min(abs(zdistances - currzdist));
    SET(parallel(i)).CurrentSlice = slice;
    SET(parallel(i)).StartSlice = slice;
    SET(parallel(i)).EndSlice = slice;

    % DEBUG
    %disp(sprintf('SET(%d).CurrentSlice z-distance: %d mm', parallel(i), round(zdistances(slice))));

    % Update images
    
   if  any(strcmp({DATA.ViewPanelsType{DATA.ViewPanels==parallel(i)}},'montage'))
%     updateselectedslices;
    %Trick it claiming that the images are linked for a short wile
   SET(NO).Linked=[SET(NO).Linked parallel(i)];
      domontageupdate=1;
   end
   
   if any(~strcmp({DATA.ViewPanelsType{DATA.ViewPanels==parallel(i)}},'montage'))
    updateoneim(parallel(i)); % this fcn changes NO !!!
    drawfunctions('drawsliceno',parallel(i));
   end
   
    %drawfunctions('drawsliceno',parallel(i));
    NO = oldNO;
  end;
  if domontageupdate
    updateselectedslices
  end
  SET(NO).Linked=oldlink;
    
end;

%----------------------------------------------------
function sameview = orientationcomparison(setindex1,setindex2)
%----------------------------------------------------
% Compares the SET.ImageOrientation between two SETs.
% Help function to updateparallelsets
%
% Return values:
% sameview : true if the orientations are parallel.

% Marten Larsson, June, 3, 2009

global SET

viewdir1 = cross(SET(setindex1).ImageOrientation(1:3),SET(setindex1).ImageOrientation(4:6))';
viewdir2 = cross(SET(setindex2).ImageOrientation(1:3),SET(setindex2).ImageOrientation(4:6))';

if dot(viewdir1, viewdir2) > 0.97 % Threshold correspond to ~15 degrees
  sameview = true;
else
  sameview = false;
  return;
end;

%----------------------------------------
function smoothendowall_Callback(no)
 %----------------------------------------
%Adjust the LV wall so that the thickness is more even. Adjustment is done
%on the endocardial side.

%Einar Heiberg

global SET NO

if nargin<1
  no = NO;
end;

tools('enableundo');

%Find slices
slices = SET(no).StartSlice:SET(no).EndSlice;
t = SET(no).CurrentTimeFrame;

%Clear interp points if there are any
if ~isempty(SET(no).EndoInterpX)
  for zloop=slices
    SET(no).EndoInterpX{t,zloop} = [];
    SET(no).EndoInterpY{t,zloop} = [];
  end
end

%Extract endo & epi
endox = SET(no).EndoX(:,t,slices);
endoy = SET(no).EndoY(:,t,slices);
epix = SET(no).EpiX(:,t,slices);
epiy = SET(no).EpiY(:,t,slices);

%Find mean line
%mx = (endox+epix)/2;
%my = (endoy+epiy)/2;

%Calculate wallthickness
wt = sqrt((endox-epix).^2+(endoy-epiy).^2);

%Create a weight function
weight = wt-repmat(min(wt),[size(wt,1) 1 1]);
weight = weight/max(weight(:));
weight = 0.15*weight.^2;

%Linear interpolation
endox2 = (1-weight).*endox+weight.*epix;
endoy2 = (1-weight).*endoy+weight.*epiy;

%assign new values
SET(no).EndoX(:,t,slices) = endox2;
SET(no).EndoY(:,t,slices) = endoy2;

%Update graphically
%Update
segment('updatemodeldisplay');
segment('updatevolume');
drawfunctions('updatenopanels',no);

%--------------------------------
function movetowardsbase_Callback
%--------------------------------
%Change current slice towards base
global DATA SET NO

laxfield = upper(DATA.ViewPanelsType{DATA.CurrentPanel});
if ismember(laxfield,{'HLA','VLA'})
  if SET(NO).(laxfield).slice>1
    switchtolongaxisslice(SET(NO).(laxfield).slice-1,laxfield);
  end
elseif strcmp(laxfield,'GLA')
  if SET(NO).HLA.slice>1
    silent = true;
    if SET(NO).VLA.slice<=1
      silent = false;
    end
    switchtolongaxisslice(SET(NO).HLA.slice-round(cos(SET(NO).GLA.angle)),'HLA',silent);
  end
  if SET(NO).VLA.slice>1
    switchtolongaxisslice(SET(NO).VLA.slice+round(sin(SET(NO).GLA.angle)),'VLA');
  end
else
  if SET(NO).CurrentSlice>1
    switchtoslice(SET(NO).CurrentSlice-1);
  end;
end

%--------------------------------
function movetowardsapex_Callback
%--------------------------------
%Change current slice towards apex
global DATA SET NO

laxfield = upper(DATA.ViewPanelsType{DATA.CurrentPanel});
if ismember(laxfield,{'HLA','VLA'})
  if SET(NO).(laxfield).slice<SET(NO).(laxfield).maxslice
    switchtolongaxisslice(SET(NO).(laxfield).slice+1,laxfield);
  end
elseif strcmp(laxfield,'GLA')
  if SET(NO).HLA.slice<SET(NO).HLA.maxslice 
    silent = true;
    if SET(NO).VLA.slice>=SET(NO).VLA.maxslice
      silent = false;
    end
    switchtolongaxisslice(SET(NO).HLA.slice+round(cos(SET(NO).GLA.angle)),'HLA',silent);
  end
  if SET(NO).VLA.slice<SET(NO).VLA.maxslice
    switchtolongaxisslice(SET(NO).VLA.slice-round(sin(SET(NO).GLA.angle)),'VLA');
  end
else
  if SET(NO).CurrentSlice<SET(NO).ZSize
    switchtoslice(SET(NO).CurrentSlice+1);
  end;
end

%----------------------------------
function playall_Callback(~)  %#ok<DEFNU>
%----------------------------------
global DATA
DATA.Run = 0;

stateandicon=iconson('play');

icon=stateandicon{2};

if icon.isindented==0
  stopmovie_Callback;
  return;
else
  playall_Helper
end


%----------------------------------
function playall_Helper 
%----------------------------------
%Starts movie display of all visible image stacks.
global DATA SET NO

DATA.Run = 1; %Start the movie :-)

%find what NO's to use
nos = unique(DATA.ViewPanels(:));
nos = nos(nos>0); %Remove zeros
maxt = 0;
for loop=nos'
  if SET(loop).TSize>maxt
    maxt = SET(loop).TSize;
  end;
end;

DATA.StartFrame = SET(NO).CurrentTimeFrame;%DATA.StartFrame = 1;
DATA.StartTime = now;

%Run this loop until stopped by setting DATA.Run==0
t = 1;

if DATA.Record
  DATA.Record = false; %otherwise nextframe will also store...
  for loop=1:SET(NO).TSize
    nextallframe_Callback;
    DATA.MovieFrame = mygetframe(DATA.imagefig);
    export('exportmovierecorder_Callback','newframe');
    
  end;
  return;
end;

%If not record then try to run according to beattime.
if ~DATA.Record
  while true
    ind = find(cat(1,SET(:).TSize)>1);
    if isempty(ind)
      ind = 1;
    else
      ind = ind(1); %Take first
    end;
    beattime = SET(ind).BeatTime;
    
    %Not recording just play along as fast as possible
    nos = unique(DATA.ViewPanels(DATA.ViewPanels > 0));
    for no=nos([SET(nos).TSize] > 1)      
      t = 1+mod(floor(rem(now-DATA.StartTime,1)*24*3600/(beattime/maxt)+DATA.StartFrame),maxt);
      %SET(no).CurrentTimeFrame = max(min(round(SET(no).TSize*(t/maxt)),SET(no).TSize),1);
      
      for linkno = SET(no).Linked
        SET(linkno).CurrentTimeFrame = max(min(round(SET(no).TSize*(t/maxt)),SET(no).TSize),1);
      end
      
      drawfunctions('updatenopanels',no);
      %Some error checking
      if isempty(DATA)
        return; %Someone has exited.
      end;

      if DATA.Run==0
        %Ensure all on the same timeframe
        for loop=1:length(SET)
          SET(loop).CurrentTimeFrame = max(min(round(SET(loop).TSize*(t/maxt)),SET(loop).TSize),1);
        end;
        return;
      end;
           
    end; %Loop over panels
    drawnow%
    pause(0.01)
    
  end; %not recording clause

end;

%------------------------------------
function playmovie_Callback(keypress) %#ok<DEFNU,INUSD>
%------------------------------------
%Starts playing current image stack as a movie.
global DATA SET NO

% if nargin > 0
%   if isequal(get(DATA.Handles.playmovieicon,'enable'),'off')
%     return
%   end
%   if isequal(get(DATA.Handles.playmovieicon,'state'),'off')
%     set(DATA.Handles.playmovieicon,'state','on');
%   elseif isequal(get(DATA.Handles.playmovieicon,'state'),'on')
%     set(DATA.Handles.playmovieicon,'state','off');
%   else %if handle is empty, icon is unavailable
%     playall_Callback('keypressed');
%     return
%   end
% end

% %Opposide since it toggles when pressed
% if isequal(get(DATA.Handles.playmovieicon,'state'),'off')
%   stopmovie_Callback;
%   return;
% end;

%Check if play all is running.
% if isequal(get(DATA.Handles.playallicon,'state'),'on')
%   stopmovie_Callback;
%   return;
% end;

if DATA.Run==1
  stopmovie_Callback;
  return;
end

%figure(DATA.imagefig);
DATA.Run = 1; %Start the movie :-)
%set(DATA.Handles.playmovieicon,'state','on');
DATA.StartFrame = SET(NO).CurrentTimeFrame;
DATA.StartTime = now;

%Run this loop until stopped by setting DATA.Run==0
while true
  if isempty(DATA)
    return; %Someone has exited.
  end;
  drawfunctions('drawsliceno');
  drawnow; %Not expose since want to be able to do callbacks
  if DATA.Record
    drawnow;
    DATA.MovieFrame = mygetframe(DATA.imagefig);
    export('exportmovierecorder_Callback','newframe');
  end;
  if DATA.Run==0
    return;
  end;
  %Do pause to release CPU burden.
  pause(0.5*SET(NO).BeatTime/SET(NO).TSize);
end;

%--------------------------
function stopmovie_Callback
%--------------------------
%End movie display.
global DATA 
DATA.Run = 0;
stateandicon=segment('iconson','play');
stateandicon{2}.undent;
DATA.Handles.permanenticonholder.render;


%--------------------------
  function checkforduplicatehide(name,state) %#ok<DEFNU>
    %--------------------------
    global DATA
    
    icons=[DATA.Handles.permanenticonholder.iconCell{:},DATA.Icons.lviconcell{:},DATA.Icons.rviconcell{:},DATA.Icons.roiflowiconcell{:},DATA.Icons.viabilityiconcell{:},...
      DATA.Icons.analysisiconcell{:}, DATA.Icons.imageiconcell{:}];
    ind=find(strcmp(name,{icons.name}));
    for i = ind
      icons(i).isindented=state;
    end
    
    
%---------------------------------------
function viewhideall_Callback(varargin) %#ok<DEFNU>
%---------------------------------------
%all hide buttons
global DATA
stateandicon=iconson('hideall');
state=stateandicon{1};
hidecell={'hideplus','hidemar','hidescar',... 'hidepins',
  'hideintersections','hideothercontour','hideinterp','hidelv','hiderv','hideroi','hidemeasure','hidepoint','hidetext'};%,'colorbar'};
stateandicon = iconson(hidecell);
availableicons = find(cellfun(@(x) isa(x,'myicon'),stateandicon(:,2)))';
if state
  for i=availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay=icon.cdataIndent;
    icon.isindented=1;
  end
  %   for name=hidecell
%     stateandicon=iconson(name);
%     icon=stateandicon{1,2};
%     icon.cdataDisplay=icon.cdataIndent;
%     icon.isindented=1;
%   end
else
  for i=availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay=icon.cdata;
    icon.isindented=0;
  end
  %   for name=hidecell
%     stateandicon=iconson(name);
%     icon=stateandicon{1,2};
%     icon.cdataDisplay=icon.cdata;
%     icon.isindented=0;
%   end
end
drawfunctions('updatevisibility');
DATA.Handles.configiconholder.render;

%-------------------------
function allhidden
%----------------------
global DATA
stateandicon=iconson('hideall');
hideallstate=stateandicon{1};
hideallicon=stateandicon{2};

if hideallstate
hidecell={'hideplus','hidemar','hidescar',...'hidepins',
  'hideintersections','hideothercontour','hideinterp','hidelv','hiderv','hideroi','hidemeasure','hidepoint','hidetext'};%,'colorbar'};
  stateandicon = iconson(hidecell);
  availableicons = cellfun(@(x) isa(x,'myicon'),stateandicon(:,2));
  state=[stateandicon{availableicons,1}];
%   for name=hidecell
%     stateandicon=iconson(name);
%     state=[state;stateandicon{1}];
%   end
   if any(state==0)
    hideallicon.undent
    DATA.Handles.permanenticonholder.render;
  end
end

%---------------------------------------
function viewhideplus_Callback(varargin) %#ok<DEFNU>
%---------------------------------------
%Toggle visibility of center + 
drawfunctions('updatevisibility');
allhidden

%---------------------------------------
function viewhidepins_Callback(varargin) %#ok<DEFNU>
%---------------------------------------
%Toggle visibility of pins
drawfunctions('updatevisibility');
allhidden
%------------------------------------------------
function viewhideintersections_Callback(varargin) %#ok<DEFNU>
%------------------------------------------------
%Toggle visibility of plane intersections
drawfunctions('updatevisibility');
allhidden
%-----------------------------------------------
function viewhideothercontour_Callback(varargin) %#ok<DEFNU>
%-----------------------------------------------
%Toggle visibility of contours from other image stacks
drawfunctions('updatevisibility');
allhidden
%-----------------------------------------------
function viewhideinterp_Callback(varargin) %#ok<DEFNU>
%-----------------------------------------------
%Toggle visibility of contours from other image stacks
drawfunctions('updatevisibility');
allhidden
%--------------------------------------
function viewhidelv_Callback(varargin) %#ok<DEFNU>
%--------------------------------------
%Toggle visibility of lv segmentation
drawfunctions('updatevisibility');
allhidden
%--------------------------------------
function viewhiderv_Callback(varargin) %#ok<DEFNU>
%--------------------------------------
%Toggle visibility of rv segmentation
drawfunctions('updatevisibility');
allhidden
%--------------------------------------
function viewhideroi_Callback(varargin) %#ok<DEFNU>
%--------------------------------------
%Toggle visibility of roi's
drawfunctions('updatevisibility');
allhidden
%--------------------------------------
function viewhidescar_Callback(varargin) %#ok<DEFNU>
%--------------------------------------
%Toggle visibility of scar contours
%global DATA 
drawfunctions('updatevisibility');
%drawfunctions('drawall',DATA.ViewMatrix);
allhidden
%--------------------------------------
function viewhidescarextent_Callback(varargin) %#ok<DEFNU>
%--------------------------------------
%Toggle visibility of auto scar extent contours
%global DATA 
drawfunctions('updatevisibility');
%drawfunctions('drawall',DATA.ViewMatrix);
allhidden
%--------------------------------------
function viewhidemar_Callback(varargin) %#ok<DEFNU>
%--------------------------------------
%Toggle visibility of MaR contours
%global DATA 
%drawfunctions('drawall',DATA.ViewMatrix);
drawfunctions('updatevisibility');
allhidden
%---------------------------------
function viewhidemeasures_Callback %#ok<DEFNU>
%---------------------------------
%Toggle visibility of measurements
drawfunctions('updatevisibility');
allhidden
%---------------------------------
function viewhidepoints_Callback %#ok<DEFNU>
%---------------------------------
%Toggle visibility of annotation points
drawfunctions('updatevisibility');
allhidden
%---------------------------------------
function viewhidetext_Callback(varargin) %#ok<DEFNU>
%---------------------------------------
%Toggle visibility of text
drawfunctions('updatevisibility');
allhidden
%---------------------------------
function viewhidecolorbar_Callback %#ok<DEFNU>
%---------------------------------
%Toggle visibility of colorbar
global DATA

stateandicon=iconson('colorbar');
state=stateandicon{1};
%stateandicon{2}.isindented=state;
if ~state 
  DATA.GUISettings.ShowColorbar = false;
else
  DATA.GUISettings.ShowColorbar = true;
end
drawfunctions('drawall',DATA.ViewMatrix);
allhidden
%-----------------------------------
function viewhidepap_Callback %#ok<DEFNU>
%---------------------------------
%Toggle visibility of papillary overlay
global DATA NO

makeviewim(DATA.CurrentPanel,NO);
drawfunctions('drawall',DATA.ViewMatrix);
allhidden
%--------------------------
function viewhidemanualinteraction_Callback %#ok<DEFNU>
%--------------------------
%toogle visibility of manual interaction of Scar/MaR
global DATA

drawfunctions('drawall',DATA.ViewMatrix);
allhidden
%-------------------------------
function viewaddtoolbar_Callback %#ok<DEFNU>
%-------------------------------
%Helper function to add a toolbar.
global DATA
h = get(0,'children');
set(h,'menubar','figure');
set(DATA.imagefig,'menubar','none');

%--------------------------------
function viewinterp_Callback(val) %#ok<DEFNU>
%--------------------------------
global DATA
if nargin < 1
  val = ~DATA.Pref.ViewInterpolated;
end
DATA.Pref.ViewInterpolated = val;
% if val
%   set(DATA.Handles.viewpixelsmenu,'Checked','off');
% %   set(DATA.Handles.viewpixelyicon,'state','off');
% else
%   set(DATA.Handles.viewpixelsmenu,'Checked','on');
% %   set(DATA.Handles.viewpixelyicon,'state','on');
% end
DATA.ViewIM = cell(size(DATA.ViewIM));
drawfunctions('drawall');

%----------------------------------
function viewspecial_Callback(mode) %#ok<DEFNU>
%----------------------------------
%Called to switch to predefined layouts of the GUI. Mode is current mode to
%use.

global DATA SET NO

switch mode
  case 'cine'
    cineno = findfunctions('findcineshortaxisno');
    %allcineno = findfunctions('findno'); shall be changed to show all cine
    if isempty(cineno)
      cineno = 1;
    end;
    DATA.ViewMatrix = [1 1];
    DATA.ViewPanels = cineno;
    DATA.ViewPanelsType = {'montage'};
    [rows,cols] = calcfunctions('calcrowscols',cineno,SET(cineno).ZSize);
    DATA.ViewPanelsMatrix = {[rows cols]};
    DATA.ViewIM = {[]}; %Clear viewim, force to recalc
    DATA.Overlay = struct('alphadata',[],'cdata',[]);
    DATA.CurrentPanel = 1;    
    NO=cineno;
    drawfunctions('drawall',DATA.ViewMatrix); 
  case 'lv'
    cineno = DATA.LVNO;
    if ~isempty(cineno)
      DATA.ViewMatrix = [1 1];
      DATA.ViewPanels = cineno;
      DATA.ViewPanelsType = {'montage'};
      [rows,cols] = calcfunctions('calcrowscols',cineno,SET(cineno).ZSize);
      DATA.ViewPanelsMatrix = {[rows cols]};
      DATA.ViewIM = {[]}; %Clear viewim, force to recalc
      DATA.Overlay = struct('alphadata',[],'cdata',[]);
      DATA.CurrentPanel = 1;
      NO=cineno;
      drawfunctions('drawall',DATA.ViewMatrix);
      
    else
      myfailed('No LV stack found.',DATA.GUI.Segment);
      return;
    end
  case 'rv'
    cineno = DATA.RVNO;
    if ~isempty(cineno)
      DATA.ViewMatrix = [1 1];
      DATA.ViewPanels = cineno;
      DATA.ViewPanelsType = {'montage'};
      [rows,cols] = calcfunctions('calcrowscols',cineno,SET(cineno).ZSize);
      DATA.ViewPanelsMatrix = {[rows cols]};
      DATA.ViewIM = {[]}; %Clear viewim, force to recalc
      DATA.Overlay = struct('alphadata',[],'cdata',[]);
      DATA.CurrentPanel = 1;
      NO=cineno;
      drawfunctions('drawall',DATA.ViewMatrix);      
      drawfunctions('drawthumbnails')
    else
      myfailed('No RV stack found.',DATA.GUI.Segment);
      return;
    end
  case 'cinescar'
    cineno = DATA.LVNO;
    if isempty(cineno)
      cineno = findfunctions('findcineshortaxisno');
    end
    scarno = findfunctions('findscarshortaxisno');
    if isempty(scarno)
      myfailed('No scar stack found.',DATA.GUI.Segment);
      return;
    end;
    if ~isempty(cineno) && ~isempty(scarno)
      % plot both scar and cine image
      DATA.ViewMatrix = [2 1];
      DATA.ViewPanels = [cineno(1) scarno(1)];
      DATA.ViewPanelsType = {'montagerow','montagerow'};
      for panel = 1:2
        no = DATA.ViewPanels(panel);
        if SET(no).ZSize>8
          DATA.ViewPanelsMatrix{panel}(1) = 2;
          DATA.ViewPanelsMatrix{panel}(2) = ceil(SET(no).ZSize/2);
        else
          DATA.ViewPanelsMatrix{panel}(1) = 1;
          DATA.ViewPanelsMatrix{panel}(2) = SET(no).ZSize;
        end;
      end
      DATA.ViewIM = {[],[]}; %Clear viewim, force to recalc
      DATA.Overlay = struct('alphadata',[],'cdata',[]);
      DATA.Overlay(2) = struct('alphadata',[],'cdata',[]);
    else
      % only plot scar image
      DATA.ViewMatrix = 1;
      DATA.ViewPanels = scarno(1);
      DATA.ViewPanelsType = {'montagerow'};
      for panel = 1
        no = DATA.ViewPanels(panel);
        if SET(no).ZSize>8
          DATA.ViewPanelsMatrix{panel}(1) = 2;
          DATA.ViewPanelsMatrix{panel}(2) = ceil(SET(no).ZSize/2);
        else
          DATA.ViewPanelsMatrix{panel}(1) = 1;
          DATA.ViewPanelsMatrix{panel}(2) = SET(no).ZSize;
        end;
      end
      DATA.ViewIM = {[]}; %Clear viewim, force to recalc
      DATA.Overlay = struct('alphadata',[],'cdata',[]);
    end
    for i = 1:numel(DATA.ViewPanels)
      if DATA.ViewPanels(i) == scarno(1)
        DATA.CurrentPanel = i;
      end
    end
    NO=DATA.ViewPanels(DATA.CurrentPanel);
    drawfunctions('drawall',DATA.ViewMatrix);
  case 'flow'
    flowno = DATA.FlowNO; 
    if isempty(flowno)
      [~,~,flowno] = findfunctions('findno');
    end
    if isempty(flowno)
      myfailed('No flow stacks found.',DATA.GUI.Segment);
      return;
    end;
    no = flowno(1); %For now take first, maybe late do a toggle...
    nop = SET(no).Flow.PhaseNo;
    if isempty(nop)
      myfailed('Flow view only works for flow image stacks with trough-plane flow.',DATA.GUI.Segment);
      return;
    end;
    DATA.ViewMatrix = [1 2];
    DATA.ViewPanels = [no nop];
    DATA.ViewPanelsType = {'one','one'};    
    DATA.ViewIM = {[],[]};
    DATA.Overlay = struct('alphadata',[],'cdata',[]);
    DATA.Overlay(2) = struct('alphadata',[],'cdata',[]);
    DATA.CurrentPanel = 1;
    NO=no;
    drawfunctions('drawall',DATA.ViewMatrix);
    switchtopanel(1);
  case 'perfusion'
    stressno = findfunctions('findstack','Perfusion Stress');
    restno = findfunctions('findstack','Perfusion Rest');      
    if isempty(stressno) && isempty(restno)
      myfailed('No perfusion stacks found.',DATA.GUI.Segment);
      return;
    end;
    if isempty(stressno)      
      DATA.ViewMatrix = [1 1];
      DATA.ViewPanels = restno;
      DATA.ViewPanelsType = {'montagerow'};
      for panel = 1
        no = DATA.ViewPanels(panel);
        DATA.ViewPanelsMatrix{panel}(1) = 1;
        DATA.ViewPanelsMatrix{panel}(2) = SET(no).ZSize;
      end
      DATA.ViewIM = {[]}; %Clear viewim, force to recalc
      DATA.Overlay = struct('alphadata',[],'cdata',[]);
    elseif isempty(restno)
      DATA.ViewMatrix = [1 1];
      DATA.ViewPanels = stressno;
      DATA.ViewPanelsType = {'montagerow'};
      for panel = 1
        no = DATA.ViewPanels(panel);
        DATA.ViewPanelsMatrix{panel}(1) = 1;
        DATA.ViewPanelsMatrix{panel}(2) = SET(no).ZSize;
      end
      DATA.ViewIM = {[]}; %Clear viewim, force to recalc
      DATA.Overlay = struct('alphadata',[],'cdata',[]);
    else
      DATA.ViewMatrix = [2 1];
      DATA.ViewPanels = [stressno restno];
      DATA.ViewPanelsType = {'montagerow','montagerow'};
      for panel = 1
        no = DATA.ViewPanels(panel);
        DATA.ViewPanelsMatrix{panel}(1) = 1;
        DATA.ViewPanelsMatrix{panel}(2) = SET(no).ZSize;
      end
      DATA.ViewIM = {[],[]}; %Clear viewim, force to recalc
      DATA.Overlay = struct('alphadata',[],'cdata',[]);
      DATA.Overlay(2) = struct('alphadata',[],'cdata',[]);
    end
    NO=DATA.ViewPanels(DATA.CurrentPanel);
    drawfunctions('drawall',DATA.ViewMatrix);
      
  case 'cinescarperf'
  case 'stress'
  otherwise
    myfailed('Unknown option to viewspecial_Callback',DATA.GUI.Segment);
    return;
end;

%------------------------
function zoomhelper(ax,f,no,panel)
%------------------------
%Helper function to zoom image stacks. Zooming is done by changing xlim &
%ylim.

global DATA SET NO
%modified by Jane Sjogren
if nargin<4
  panel=DATA.CurrentPanel;
end
if nargin<3
  no=NO;
end
if nargin<2
  myfailed('Expected two input arguments.',DATA.GUI.Segment);
  return;
end;

if f < 1 && strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},'montagesegmented')
  mfzs = SET(NO).MontageFitZoomState;
  if mfzs(2)-mfzs(1)>=size(DATA.ViewIM{1},1) && ...
      mfzs(4)-mfzs(3)>=size(DATA.ViewIM{1},2)
    viewimage_Callback('montage');
    return
  end
end

%Get old position
temp = [get(ax,'xlim') get(ax,'ylim')];
oldxmid = 0.5*(temp(1)+temp(2));
oldymid = 0.5*(temp(3)+temp(4));
oldxspan = temp(2)-temp(1);
oldyspan = temp(4)-temp(3);

%calculate new xlim and ylim
xlim = [oldxmid - 0.5*oldxspan/f ...
        oldxmid + 0.5*oldxspan/f];
ylim = [oldymid - 0.5*oldyspan/f ...
        oldymid + 0.5*oldyspan/f];
xspan=(xlim(2)-xlim(1));
yspan=(ylim(2)-ylim(1));

% Find view ratio from axes height/width, better than DATA.ViewMatrix
set(ax,'Units','pixels');
pos=get(ax,'Position');
set(ax,'Units','normalized');
viewratio=[pos(4)/pos(3) pos(3)/pos(4)];
%viewratio=[DATA.ViewMatrix(2)/DATA.ViewMatrix(1) DATA.ViewMatrix(1)/DATA.ViewMatrix(2)];

%handle different resolution
xsz = SET(no).XSize;
ysz = SET(no).YSize;
xres = SET(no).ResolutionX;
yres = SET(no).ResolutionY;
switch DATA.ViewPanelsType{panel}    
  case 'hla'
    xsz = SET(no).ZSize;
    xres = SET(no).SliceThickness + SET(no).SliceGap;
  case 'vla'
    xsz = SET(no).ZSize;
    ysz = SET(no).XSize;
    xres = SET(no).SliceThickness + SET(no).SliceGap;
    yres = SET(no).ResolutionX;
  case 'gla'
    xsz = SET(no).ZSize;
    xres = SET(no).SliceThickness + SET(no).SliceGap;
    ysz = size(DATA.ViewIM{panel},2);
    yres = SET(no).ResolutionY*cos(SET(no).GLA.angle)^2+ ...
      SET(no).ResolutionX*sin(SET(no).GLA.angle)^2;
end
viewratio=viewratio.*[yres/xres xres/yres];

%check if xlim and ylim shall be adjusted to handle image with different x and y
%size or view with different x and y size
xsize=size(DATA.ViewIM{panel},2);
ysize=size(DATA.ViewIM{panel},1);
%This is necessary since the size is doubled when using interpolated view.
if DATA.Pref.ViewInterpolated
  xsize = xsize / 2;
  ysize = ysize / 2;
end

if f>1
  if xsize>ysize
    if xspan*viewratio(1)>oldyspan
      ylim=[0.5 ysize+0.5];
    elseif xspan*viewratio(1)<oldyspan
      ylim=[oldymid - 0.5*xspan*viewratio(1)...
            oldymid + 0.5*xspan*viewratio(1)];
    end
  else
    if yspan*viewratio(2)>oldxspan
      xlim=[0.5 xsize+0.5];
    elseif yspan*viewratio(2)<oldxspan
      xlim=[oldxmid - 0.5*yspan*viewratio(2)...
            oldxmid + 0.5*yspan*viewratio(2)];
    end
  end
else
  if xsize>ysize
    if xspan<=xsize
      if yspan*viewratio(2)<oldxspan || yspan>=ysize
        ylim=[0.5 ysize+0.5];
      elseif yspan*viewratio(2)>oldxspan
        ylim=[oldymid - 0.5*xspan*viewratio(1)...
              oldymid + 0.5*xspan*viewratio(1)];
      end
    end
  else
    if yspan<=ysize
      if xspan*viewratio(1)<oldyspan || xspan>=xsize
        xlim=[0.5 xsize+0.5];
      elseif xspan*viewratio(1)>oldyspan*viewratio(2)
        xlim=[oldxmid - 0.5*yspan*viewratio(2)...
              oldxmid + 0.5*yspan*viewratio(2)];
      end
    end
  end
end

set(ax,'xlim',xlim(:),'ylim',ylim(:));

% Get linked axes.
if length(SET(no).Linked) > 1
  nos = SET(no).Linked;
  panelstodo = find(ismember(DATA.ViewPanels,nos));
else
  panelstodo=panel;
%   nos=no;
end;

xlim = get(DATA.Handles.imageaxes(panel),'xlim');
ylim = get(DATA.Handles.imageaxes(panel),'ylim');
for ploop=panelstodo;
  switch DATA.ViewPanelsType{ploop}
    case 'montage'
      SET(DATA.ViewPanels(ploop)).MontageZoomState = [xlim(:);ylim(:)];
    case 'montagerow'
      SET(DATA.ViewPanels(ploop)).MontageRowZoomState = [xlim(:);ylim(:)];
    case {'montagefit','sax3','montagesegmented'}
      SET(DATA.ViewPanels(ploop)).MontageFitZoomState = [xlim(:);ylim(:)];
    case {'one','ortho'}
      SET(DATA.ViewPanels(ploop)).NormalZoomState = [xlim(:);ylim(:)];
    case 'hla'
      SET(DATA.ViewPanels(ploop)).HLA.ZoomState = [xlim(:);ylim(:)];
    case 'vla'
      SET(DATA.ViewPanels(ploop)).VLA.ZoomState = [xlim(:);ylim(:)];
  end;
end

%set dataaspectratio so that adjustment in zooming does not affect the
%aspectratio of the pixels
for ploop=panelstodo
  set(DATA.Handles.imageaxes(ploop),'plotboxaspectratio',[...
    ysz*yres ...
    xsz*xres 1]);
  set(DATA.Handles.imageaxes(ploop),'dataaspectratio',...
    [1/yres ...
    1/xres 1]);
end


%---------------------------
function viewzoomin_Callback %#ok<DEFNU>
%---------------------------
%Zooms in current view panel
global DATA

set(DATA.Handles.viewzoominicon,'state','off');

if strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},'mmodetemporal')
  updatemmode([],'one')
  return
end

if strcmp(get(gcbf,'currentmodifier'),'shift'),
  for loop=1:length(DATA.Handles.imageaxes)
    zoomhelper(DATA.Handles.imageaxes(loop),1.2,DATA.ViewPanels(loop),loop);
  end
else
  zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1.2);
end
drawfunctions('viewupdatetextposition');
drawfunctions('viewupdateannotext');

%----------------------------
function viewzoomout_Callback  %#ok<DEFNU>
%----------------------------
%Zooms out in current view panel
global DATA

set(DATA.Handles.viewzoomouticon,'state','off');

if strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},'mmodetemporal')
  updatemmode([],'three')
  return
end

if strcmp(get(gcbf,'currentmodifier'),'shift'),
  for loop=1:length(DATA.Handles.imageaxes)
    zoomhelper(DATA.Handles.imageaxes(loop),1/1.2,DATA.ViewPanels(loop),loop);
  end
else
  zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1/1.2);
end
drawfunctions('viewupdatetextposition');
drawfunctions('viewupdateannotext');

%--------------------------------------
function viewpandir_Callback(direction) %#ok<DEFNU>
%--------------------------------------
%Pans current view panel
global DATA 

ca = DATA.Handles.imageaxes(DATA.CurrentPanel);
xlim = get(ca,'xlim'); xchange = (xlim(2)-xlim(1))/10;
ylim = get(ca,'ylim'); ychange = (ylim(2)-ylim(1))/10;

switch direction
  case 'up'
    ylim = ylim+ychange;
  case 'down'
    ylim = ylim-ychange;
  case 'left'
    xlim = xlim+xchange;
  case 'right'
    xlim = xlim-xchange;
end;

set(ca,'xlim',xlim,'ylim',ylim);

%------------------------------
function viewhidepanel_Callback %#ok<DEFNU>
%------------------------------
%Hides current panel
global DATA

if length(DATA.ViewPanels)==1
  myfailed('If hiding this panel then no panels will be visible.',DATA.GUI.Segment);
  return;
end;

DATA.ViewPanels(DATA.CurrentPanel) = 0;
DATA.ViewPanelsType{DATA.CurrentPanel} = '';
DATA.ViewIM{DATA.CurrentPanel} = [];
DATA.Overlay(DATA.CurrentPanel) = struct('alphadata',[],'cdata',[]);
drawfunctions('drawall',DATA.ViewMatrix);

%---------------------------------------
function viewhideimagestack_Callback(no) %#ok<DEFNU>
%---------------------------------------
%Remove clicked image stack from view panels
global DATA SET

if nargin==0
  %Clicked position is still valid.
  no = getclickedpreview;%should possibly be removed as in fileclosecurrentimagestack
%   p = get(gca,'CurrentPoint');
%   x=p(1);
%   no = floor(x/DATA.GUISettings.ThumbnailSize)+1;
end;
  
if (no<1)||(no>length(SET))
  return;
end;

if not(any(DATA.ViewPanels==no))
  myfailed('Can not hide a non visible image stack',DATA.GUI.Segment);
  return;
end;

lind = (DATA.ViewPanels~=no)&DATA.ViewPanels; %The last removes 0
if not(sum(lind)>0)
  myfailed('If removing this image stack then no image stacks will be visible.',DATA.GUI.Segment);
  return;
end;

DATA.ViewPanels(lind) = 0;
DATA.ViewPanelsType{lind} = '';
DATA.ViewIM{lind} = [];
DATA.Overlay(lind) = struct('alphadata',[],'cdata',[]);
DATA.CurrentPanel = 1;

drawfunctions('drawall',DATA.ViewMatrix);

%-----------------------------------
function viewimagestack_Callback(no) %#ok<DEFNU>
%-----------------------------------
%Make clicked image stack visible
global DATA SET

if nargin==0
  %Clicked position is still valid.
  no = getclickedpreview;%should possibly be removed as in fileclosecurrentimagestack
end;

if (no<1)||(no>length(SET))
  return;
end;
  
addtopanels(no);
drawfunctions('drawall',length(DATA.ViewPanels));

%-----------------------------------
function viewallimagestacks_Callback
%-----------------------------------
%Displays all image stacks.
global DATA SET

DATA.ViewPanels = [];
DATA.ViewPanelsType = {};
DATA.ViewPanelsMatrix = {};
DATA.ViewIM = {};
DATA.Overlay = [];
DATA.CurrentPanel = 1;

for loop=1:length(SET)
  addtopanels(loop);
end;
drawfunctions('drawall',length(DATA.ViewPanels));

%--------------------------------
function viewimage_Callback(type)
%--------------------------------
%Select type of image to view on screen.
global DATA SET NO
if DATA.Silent
  return;
end;

panel = DATA.CurrentPanel;

oldmode = DATA.ViewPanelsType{panel};

if nargin==0
  type = oldmode;
end;

doall = false;
if ismember(oldmode,{'mmodetemporal','mmodespatial'}) || ...
    ismember(oldmode,{'ortho','hla','vla','gla'})
  %Old mode is mmode or orthoview, clear it.
  ind = true(1,length(DATA.ViewPanelsType));
  for loop=1:length(DATA.ViewPanelsType)
    if ismember(DATA.ViewPanelsType{loop},{'mmodespatial','ortho'})
      panel = loop; %Keep this.
      DATA.CurrentPanel = panel;
      DATA.ViewPanelsType{loop} = 'one'; %Change from mmode to one
    end;
    if ismember(DATA.ViewPanelsType{loop},{'mmodetemporal','hla','vla','gla'})
      ind(loop) = false; %Mark for deletion
    end;    
  end;
  DATA.ViewPanels = DATA.ViewPanels(ind);
  DATA.ViewPanelsType = DATA.ViewPanelsType(ind);
  DATA.ViewPanelsMatrix = DATA.ViewPanelsMatrix(ind);
  DATA.ViewIM = DATA.ViewIM(ind);  
  DATA.Overlay = DATA.Overlay(ind);
  doall = true; %Update complete screen since one panel is removed.
  DATA.CurrentPanel = 1;
  if strcmp(DATA.CurrentTool,'orthoview')
    updatetool('select');
  end
end;
  
switch type
  case {'one','ortho','hla','vla'}

    if length(SET(DATA.ViewPanels(panel)).Linked) > 1
      nos = SET(DATA.ViewPanels(panel)).Linked;
      panelstodo=find(ismember(DATA.ViewPanels,nos));
      nos=DATA.ViewPanels(panelstodo);
    else
      panelstodo=panel;
      nos=NO;
    end;
    
    for loop=1:length(panelstodo)
      panelloop=panelstodo(loop);
      no=nos(loop);

      DATA.ViewPanelsType{panelloop} = type;
      SET(no).StartSlice = SET(no).CurrentSlice;
      SET(no).EndSlice = SET(no).CurrentSlice;
      DATA.ViewIM{panelloop} = [];
      DATA.Overlay(panelloop) = struct('alphadata',[],'cdata',[]);
      makeviewim(panelloop,no);
    end
    
    if doall
      drawfunctions('drawall');
    else
      drawfunctions('drawimageno',DATA.ViewPanels(panel));
    end;

  case 'mmode'
    for loop=1:length(DATA.ViewPanelsType)
      if isequal(DATA.ViewPanelsType{loop},'mmodespatial')
        DATA.ViewPanelsType{loop} = 'one';
      end
      if isequal(DATA.ViewPanelsType{loop},'mmodetemporal')
        %remove existing mmode images
        DATA.ViewPanels(loop) = [];
        DATA.ViewPanelsType(loop) = [];
        DATA.ViewPanelsMatrix(loop) = [];
        DATA.ViewIM(loop) = [];
        DATA.Overlay(loop) = [];
        if panel > loop
          panel = panel-1;
          DATA.CurrentPanel = panel;
        end
        break
        %myfailed('Can not have two mmode images visible at the same time.',DATA.GUI.Segment);
        %return;
      end;
    end;   
    if SET(NO).TSize>1
        
      DATA.ViewPanelsType{panel} = 'mmodespatial';
      SET(NO).StartSlice = SET(NO).CurrentSlice;
      SET(NO).EndSlice = SET(NO).CurrentSlice;

      DATA.ViewIM{panel} = []; %Remove is montage
      empol = struct('alphadata',[],'cdata',[]);
      DATA.Overlay(panel) = empol;
      
      %switch position to make mmodespatial first
      pos = 1;
      temp=DATA.ViewPanels(panel); DATA.ViewPanels(panel)=DATA.ViewPanels(pos); DATA.ViewPanels(pos)=temp;
      temp=DATA.ViewPanelsType{panel}; DATA.ViewPanelsType{panel}=DATA.ViewPanelsType{pos}; DATA.ViewPanelsType{pos}=temp;
      temp=DATA.ViewIM{panel}; DATA.ViewIM{panel}=DATA.ViewIM{pos}; DATA.ViewIM{pos}=temp;      
      DATA.CurrentPanel = 1;
      
      %Add view panel mmodetemporal as number second
      DATA.ViewPanels = cat(2,DATA.ViewPanels(1),DATA.ViewPanels(1),DATA.ViewPanels(2:end));
      DATA.ViewPanelsType = cat(2,DATA.ViewPanelsType(1),{'mmodetemporal'},DATA.ViewPanelsType(2:end));
      DATA.ViewPanelsMatrix = cat(2,DATA.ViewPanelsMatrix(1),{[]},DATA.ViewPanelsMatrix(2:end));
      DATA.ViewIM = cat(2,DATA.ViewIM(1),{[]},DATA.ViewIM(2:end));
      DATA.Overlay = cat(2,DATA.Overlay(1),empol,DATA.Overlay(2:end));
    else
      viewimage_Callback('one');
    end;
    drawfunctions('drawall');
  case 'montage'

    if length(SET(DATA.ViewPanels(panel)).Linked) > 1
      nos = SET(DATA.ViewPanels(panel)).Linked;
      panelstodo=find(ismember(DATA.ViewPanels,nos));
      nos=DATA.ViewPanels(panelstodo);
    else
      panelstodo=panel;
      nos=NO;
    end;

    for loop=1:length(panelstodo)
      panelloop=panelstodo(loop);
      no=nos(loop); 
      DATA.ViewPanelsType{panelloop} = type;
      DATA.ViewIM{panelloop} = [];
      DATA.Overlay(panelloop) = struct('alphadata',[],'cdata',[]);
      makeviewim(panelloop,no); 
    end;

    updatemodeldisplay;
    if doall
      drawfunctions('drawall');
    else
      drawfunctions('drawimageno',DATA.ViewPanels(panel));
    end;
    
  case {'montagerow','montagefit','sax3','montagesegmented'}
    DATA.ViewPanelsType{panel} = type;
    DATA.ViewIM{panel} = [];
    DATA.Overlay(panel) = struct('alphadata',[],'cdata',[]);
    makeviewim(panel,NO);
    updatemodeldisplay;
    if strcmp(type,'montagefit')
      centeronslice(SET(NO).CurrentSlice,NO,SET(NO).ZSize);
    end
    if doall
      drawfunctions('drawall');
    else
      drawfunctions('drawimagepanel',DATA.CurrentPanel);
    end;        
  otherwise
    myfailed('Unknown option to viewimage.',DATA.GUI.Segment);
    return;
end;

drawfunctions('showedits',NO);
drawfunctions('drawintersections');

%--------------------------------------
function viewmanualinteraction_Callback  %#ok<DEFNU>
%--------------------------------------
%Displays a GUI indicating in which slices and timeframes manual
%interaction have been made for LV segmentation.
global DATA SET NO

if isempty(SET(NO).EndoX) || all(isnan(SET(NO).EndoX(:)))
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end;

fig = openfig('manualedit.fig','reuse');
myadjust(fig,DATA.GUI.Segment);
set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

% Generate a structure of handles to pass to callbacks, and store it.
handles = guihandles(fig);

%Fix colors
im = zeros(SET(NO).ZSize,SET(NO).TSize);
for zloop=1:SET(NO).ZSize
  for tloop=1:SET(NO).TSize
    if isnan(SET(NO).EndoX(1,tloop,zloop))
      %Not segmented
      im(zloop,tloop,1) = 0; %r
      im(zloop,tloop,2) = 0; %g
      im(zloop,tloop,3) = 1; %b
    else
      if (~isempty(SET(NO).EndoPinX))&&(not(isempty(SET(NO).EndoPinX{tloop,zloop})))
        %Pin
        im(zloop,tloop,1) = 0; %r
        im(zloop,tloop,2) = 1; %g
        im(zloop,tloop,3) = 1; %b
      else
        if SET(NO).EndoDraged(tloop,zloop)
          %Draged
          im(zloop,tloop,1) = 1; %r
          im(zloop,tloop,2) = 0; %g
          im(zloop,tloop,3) = 1; %b
        else
          %fully automatic
          im(zloop,tloop,1) = 1; %r
          im(zloop,tloop,2) = 0; %g
          im(zloop,tloop,3) = 0; %b
        end;
      end;
    end;
  end;
end;

axes(handles.endoaxes);  %#ok<MAXES>
image([1 SET(NO).TSize]-0.5,[1 SET(NO).ZSize]-0.5,im);
set(gca,'xtick',1:SET(NO).TSize,'ytick',1:SET(NO).ZSize);
ylabel('Apex    =>     Base');
xlabel('Timeframe');
grid on;

%Fix colors epi
im = zeros(SET(NO).ZSize,SET(NO).TSize);
for zloop=1:SET(NO).ZSize
  for tloop=1:SET(NO).TSize
    if isnan(SET(NO).EpiX(1,tloop,zloop))
      %Not segmented
      im(zloop,tloop,1) = 0; %r
      im(zloop,tloop,2) = 0; %g
      im(zloop,tloop,3) = 1; %b
    else
      if ~isempty(SET(NO).EpiPinX)&&~isempty(SET(NO).EpiPinX{tloop,zloop})
        %Pin
        im(zloop,tloop,1) = 0; %r
        im(zloop,tloop,2) = 1; %g
        im(zloop,tloop,3) = 1; %b
      else
        if SET(NO).EpiDraged(tloop,zloop)
          %Draged
          im(zloop,tloop,1) = 1; %r
          im(zloop,tloop,2) = 0; %g
          im(zloop,tloop,3) = 1; %b
        else
          %fully automatic
          im(zloop,tloop,1) = 1; %r
          im(zloop,tloop,2) = 0; %g
          im(zloop,tloop,3) = 0; %b
        end;
      end;
    end;
  end;
end;

axes(handles.epiaxes); %#ok<MAXES>
image([1 SET(NO).TSize]-0.5,[1 SET(NO).ZSize]-0.5,im);
set(gca,'xtick',1:SET(NO).TSize,'ytick',1:SET(NO).ZSize);
ylabel('Apex    =>     Base');
xlabel('Timeframe');
grid on;

%-----------------------------------
function figure_DeleteFcn %#ok<DEFNU>
%-----------------------------------
%Shut down Segment in a controlled manner, and remove
%global DATA SET NO. Note that filequit callback takes
%care of asking user yesno.
global DATA 

delete(DATA.fig);
DATA = [];

%--------------------------
function t = getframenumber
%--------------------------
%Calculates what frame to show when playing a movie if storing a movie 
%then show next frame otherwise user timer info.
global DATA SET NO

if DATA.Record
  t = SET(NO).CurrentTimeFrame+1;
  if t>SET(NO).TSize
    t = 1;
  end;
else
  %User timer info
  t = 1+mod(floor(rem(now-DATA.StartTime,1)*24*3600/(SET(NO).BeatTime/SET(NO).TSize)+DATA.StartFrame),SET(NO).TSize);
end;

%------------------------------------
function selectallslices_Callback
%------------------------------------
%Selects all slices in current image stack
global SET NO

SET(NO).StartSlice = 1;
SET(NO).EndSlice = SET(NO).ZSize;
segment('updateselectedslices');

%--------------------------------------
function unselectallslices_Callback(no)
%--------------------------------------
%Unselects slices in current image stack
global SET NO

if nargin == 0
  no=NO;
end

SET(no).StartSlice = [];
SET(no).EndSlice = [];

%SET(NO).CurrentSlice = [];   
% For various reasons this is a bad idea. Many things rely on it being set,
% and difficult to gaurentee that it is set to something sensible. Best to
% check with StartSlice:EndSlice in cases where 'unselection' is desired.

updateselectedslices;

%---------------------
function unlighttool(h) %#ok<DEFNU>
%---------------------
%Helper function to change color of a tool to represent unselected.
global DATA

set(h,'value',0,'backgroundcolor',DATA.GUISettings.ButtonColor);
%set(h,'value',0,'backgroundcolor',[0.6 0.6 0.6]);

%-----------------------
function highlighttool(h)
%-----------------------
%Helper function to change color of a tool to represent that the tool is
%selected.
global DATA

%set(h,'value',1,'backgroundcolor',[1 1 1]);
% JU: ButtonSelectedColor below was [0.6 0.6 0.6]
set(h,'value',1,'backgroundcolor',DATA.GUISettings.ButtonSelectedColor); 
%temp = get(h(1),'cdata');
%temp(isnan(temp)) = 0.5;
%set(h(1),'cdata',temp);

%-------------------------------
function doputpin_Callback(type)
%-------------------------------
%Put pins. Called when clicked, puts an pin and refines.
global DATA SET NO

if nargin==0
  type = 'endo';
end;

if nargin<2
  dolv=true;
end;

switch type
  case 'endo'
    isendo = true;
    dolv = true;
  case 'epi'
    isendo = false;
    dolv = true;
  case 'rvendo'
    isendo = true;
    dolv = false;    
end;

[x,y,slice] = getclickedcoords;
% If slice has changed, make sure montage/one are in sync
if (slice>SET(NO).ZSize)
  return;
end
switchtoslice(slice);

tools('enableundo');

if dolv
  %LV
  if isendo
    if isempty(SET(NO).EndoPinX)
      SET(NO).EndoPinX = cell(SET(NO).TSize,SET(NO).ZSize);
      SET(NO).EndoPinY = cell(SET(NO).TSize,SET(NO).ZSize);      
    end;
    SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice} = [...
      SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice} ; ...
      y];
    SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice} = [...
      SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice} ; ...
      x];
    if (~DATA.ThisFrameOnly)&&(SET(NO).TSize>1)
      lvpeter('segmentrefineendo_Callback',false,true);      
    end;
    
  else
    if isempty(SET(NO).EpiPinX)
      SET(NO).EpiPinX = cell(SET(NO).TSize,SET(NO).ZSize);
      SET(NO).EpiPinY = cell(SET(NO).TSize,SET(NO).ZSize);      
    end;    
    SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice} = [...
      SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice} ; ...
      y];
    SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice} = [...
      SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice} ; ...
      x];

    if (~DATA.ThisFrameOnly)&&(SET(NO).TSize>1)
      lvpeter('segmentrefineepi_Callback',false);      
    end;
  end;
else
  %RV
  if isendo
    if isempty(SET(NO).RVEndoPinX)
      SET(NO).RVEndoPinX = cell(SET(NO).TSize,SET(NO).ZSize);
      SET(NO).RVEndoPinY = cell(SET(NO).TSize,SET(NO).ZSize);      
    end;
    SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice} = [...
      SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice} ; ...
      y];
    SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice} = [...
      SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice} ; ...
      x];
    
    if (~DATA.ThisFrameOnly) && (SET(NO).TSize>1)
      rv('segmentrefinervendo_Callback',false,true);
    end;
  end;
end;

checkconsistency;
updateselectedslices;
updatemodeldisplay;
drawfunctions('drawsliceno');

%-----------------------------
function volumeaxes_Buttondown %#ok<DEFNU>
%-----------------------------
%Called when user has clicked volume graph, sets current
%timeframe to clicked point.

global DATA SET NO

if NO == DATA.LVNO
  
  [x,y] = mygetcurrentpoint(gca);
%   t = round((x/1000+1.5*SET(NO).TIncr)/(SET(NO).TIncr));
%   t = max(min(t,SET(NO).TSize),1);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
SET(NO).CurrentTimeFrame = t;
  
  %Check if close to ED or ES.
  ylim = get(DATA.Handles.volumeaxes,'ylim');
  if abs(SET(NO).EDT-t)<3
    if abs(ylim(2)-y)<(0.2*(ylim(2)-ylim(1)))
      esed_Buttondown('ed');
      return;
    end;
  end;
  if abs(SET(NO).EST-t)<3
    if abs(y-(0.25*(ylim(2)-ylim(1))+ylim(1)))<(0.2*(ylim(2)-ylim(1))) %ES textg is not really at bottom.
      esed_Buttondown('es');
      return;
    end;
  end;
  DATA.updatevolumeaxes;
  drawfunctions('drawsliceno');
end

%--------------------------
function esed_Buttonup(type) %#ok<DEFNU>
%---------------------------
%Buttonup function when dragging ES or ED markers in volume graph
global DATA SET NO
no = DATA.LVNO;

if NO == no
  %Get position
  [x] = mygetcurrentpoint(gca);
%   t = round((x/1000+1.5*SET(NO).TIncr)/(SET(NO).TIncr));
%   t = max(min(t,SET(NO).TSize),1);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);  
t=max([t,1]);
  %Store position
  switch type
    case 'es'
      SET(NO).EST = t;
    case 'ed'
      SET(NO).EDT = t;
  end;
  
  %disp(sprintf('ed:%d es:%d',SET(NO).EDT,SET(NO).EST));
  
  %Restore
  set(DATA.fig,'WindowButtonUpFcn','segment(''buttonup_Callback'')');
  set(DATA.fig,'WindowButtonMotionFcn','');
  SET(NO).CurrentTimeFrame = t;
  DATA.updatevolumeaxes
  DATA.measurementreportupdate
  updatevolume;
  drawfunctions('drawsliceno');
end

%-------------------------
function esed_Motion(type) %#ok<DEFNU>
%-------------------------
%Motion function when dragging ES or ED markers in volume graph.
global DATA SET NO
no = DATA.LVNO;

if NO == no
  [x] = mygetcurrentpoint(gca);
%   t = round((x/1000+1.5*SET(NO).TIncr)/(SET(NO).TIncr));
%   t = max(min(t,SET(NO).TSize),1);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
  t=max([t,1]);
  switch type
    case 'es'
      p = get(DATA.Handles.estext,'position');
      p(1) = x;
      set(DATA.Handles.estext,'position',p);
      SET(NO).EST = t;
    case 'ed'
      p = get(DATA.Handles.edtext,'position');
      p(1) = x;
      set(DATA.Handles.edtext,'position',p);
      SET(NO).EDT = t;
  end;
  
  SET(NO).CurrentTimeFrame = t;
  drawfunctions('drawsliceno');
  
end

%-----------------------------
function esed_Buttondown(type) 
%-----------------------------
%Buttondown function when dragging ES or ED markers in volume graph.
global DATA SET NO
no = DATA.LVNO;

if NO == no
  switch lower(type)
    case 'ed'
      set(DATA.fig,'WindowButtonUpFcn',...
        'segment(''esed_Buttonup'',''ed'')');
      set(DATA.fig,'WindowButtonMotionFcn',...
        'segment(''esed_Motion'',''ed'')');
      SET(NO).CurrentTimeFrame = SET(NO).EDT;
    case 'es'
      set(DATA.fig,'WindowButtonUpFcn',...
        'segment(''esed_Buttonup'',''es'')');
      set(DATA.fig,'WindowButtonMotionFcn',...
        'segment(''esed_Motion'',''es'')');
      SET(NO).CurrentTimeFrame = SET(NO).EST;
  end;
end

%-----------------------------
function timebaraxes_Buttondown %#ok<DEFNU>
%-----------------------------
%Called when user has clicked time bar graph, sets current
%timeframe to clicked point.

global DATA SET NO

[x,y] = mygetcurrentpoint(gca);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
%t = round(((x-xlim(1))/1000+1.5*SET(NO).TIncr)/(SET(NO).TIncr));
%t = max(min(t,SET(NO).TSize),1);
SET(NO).CurrentTimeFrame = t;

%Check if close to ED or ES.
ylim = get(DATA.Handles.timebaraxes,'ylim');
if abs(SET(NO).EDT-t)<3
  if abs(ylim(2)-y)<(0.2*(ylim(2)-ylim(1)))
    esedtimebar_Buttondown('ed');
    return;
  end;
end;
if abs(SET(NO).EST-t)<3
  if abs(y-(0.25*(ylim(2)-ylim(1))+ylim(1)))<(0.2*(ylim(2)-ylim(1))) %ES textg is not really at bottom.
    esedtimebar_Buttondown('es');
    return;
  end;
end;

nos = SET(NO).Linked;
for loop=setdiff(nos,NO)
  SET(loop).CurrentTimeFrame = SET(NO).CurrentTimeFrame;
  if SET(loop).CurrentTimeFrame > SET(loop).TSize
    SET(loop).CurrentTimeFrame = SET(loop).TSize;
  end;
end;

for nloop=1:length(nos)
  drawfunctions('updatenopanels',nos(nloop));
end;

% drawfunctions('drawsliceno');

%--------------------------
function esedtimebar_Buttonup(type) %#ok<DEFNU>
%---------------------------
%Buttonup function when dragging ES or ED markers.
global DATA SET NO

%Get position
[x] = mygetcurrentpoint(gca);
% t = round((x/1000+1.5*SET(NO).TIncr)/(SET(NO).TIncr));
% t = max(min(t,SET(NO).TSize),1);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
t=max([t,1]);

%Store position
switch type
  case 'es'
    SET(NO).EST = t;
  case 'ed'
    SET(NO).EDT = t;
end;

%Restore
set(DATA.fig,'WindowButtonUpFcn','segment(''buttonup_Callback'')');
set(DATA.fig,'WindowButtonMotionFcn','');
SET(NO).CurrentTimeFrame = t;
updatevolume;
drawfunctions('drawsliceno');

%-------------------------
function esedtimebar_Motion(type) %#ok<DEFNU>
%-------------------------
%Motion function when dragging ES or ED markers in time bar graph.
global DATA SET NO

[x] = mygetcurrentpoint(gca);
% t = round((x/1000+1.5*SET(NO).TIncr)/(SET(NO).TIncr));
% t = max(min(t,SET(NO).TSize),1);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
t=max([t,1]);

switch type
  case 'es'
    p = get(DATA.Handles.estimebartext,'position');
    p(1) = x;
    set(DATA.Handles.estimebartext,'position',p);
    SET(NO).EST = t;    
  case 'ed'
    p = get(DATA.Handles.edtimebartext,'position');
    p(1) = x;
    set(DATA.Handles.edtimebartext,'position',p);
    SET(NO).EDT = t;        
end;

SET(NO).CurrentTimeFrame = t;
drawfunctions('drawsliceno');

%-----------------------------
function esedtimebar_Buttondown(type) 
%-----------------------------
%Buttondown function when dragging ES or ED markers in time bar graph.
global DATA SET NO

switch lower(type)
  case 'ed'
    set(DATA.fig,'WindowButtonUpFcn',...
      'segment(''esedtimebar_Buttonup'',''ed'')');
    set(DATA.fig,'WindowButtonMotionFcn',...
      'segment(''esedtimebar_Motion'',''ed'')');
    SET(NO).CurrentTimeFrame = SET(NO).EDT;
  case 'es'
    set(DATA.fig,'WindowButtonUpFcn',...
      'segment(''esedtimebar_Buttonup'',''es'')');
    set(DATA.fig,'WindowButtonMotionFcn',...
      'segment(''esedtimebar_Motion'',''es'')');
    SET(NO).CurrentTimeFrame = SET(NO).EST;    
end;

%-----------------------------
function flowaxes_Buttondown %#ok<DEFNU>
%-----------------------------
%Called when user has clicked flow graph, sets current
%timeframe to clicked point.

global DATA SET NO

if ismember(NO,SET(DATA.FlowNO).Linked)
  
  [x,y] = mygetcurrentpoint(gca);
%   t = round((x/1000+1.5*SET(NO).TIncr)/(SET(NO).TIncr));
%   t = max(min(t,SET(NO).TSize),1);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
SET(NO).CurrentTimeFrame = t;
  
  nos = SET(NO).Linked;
  for loop=setdiff(nos,NO)
    SET(loop).CurrentTimeFrame = SET(NO).CurrentTimeFrame;
    if SET(loop).CurrentTimeFrame > SET(loop).TSize
      SET(loop).CurrentTimeFrame = SET(loop).TSize;
    end;
  end;
  
  DATA.updateflowaxes;
  
  for nloop=1:length(nos)
    drawfunctions('updatenopanels',nos(nloop));
  end;
%   drawfunctions('drawsliceno');
end

%-------------------------------------
function clickedpin_Callback(type) %#ok<DEFNU>
%-------------------------------------
%Called when user clicks on a pin.
global DATA SET NO

%Check what type of click
switch get(DATA.imagefig,'SelectionType')
  case 'normal'
    [x,y,slice] = getclickedcoords;

    %Find what pin.
    switch type
      case 'endo'
        if isempty(SET(NO).EndoPinX)
          return; %Not sure if can happen, but...
        end;
        xm = SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice};
        ym = SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice};
        set(DATA.imagefig,'WindowButtonMotionFcn',...
          sprintf('%s(''pin_Motion'',''endo'')',mfilename));
      case 'epi'
        if isempty(SET(NO).EpiPinX)
          return; %Not sure if can happen, but...
        end;
        xm = SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice};
        ym = SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice};
        set(DATA.imagefig,'WindowButtonMotionFcn',...
          sprintf('%s(''pin_Motion'',''epi'')',mfilename));
      case 'rvendo'
        if isempty(SET(NO).RVEndoPinX)
          return; %Not sure if can happen, but...
        end;
        xm = SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice};
        ym = SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice};
        set(DATA.imagefig,'WindowButtonMotionFcn',...
          sprintf('%s(''pin_Motion'',''rvendo'')',mfilename));        
      case 'rvepi'
        %EiHMOD: Bug report from Lene R
        if isempty(SET(NO).RVEpiPinX)
          return; %Not sure if can happen, but...
        end;
        xm = SET(NO).RVEpiPinX{SET(NO).CurrentTimeFrame,slice};
        ym = SET(NO).RVEpiPinY{SET(NO).CurrentTimeFrame,slice};
        set(DATA.imagefig,'WindowButtonMotionFcn',...
          sprintf('%s(''pin_Motion'',''rvepi'')',mfilename));        
    end;
    set(DATA.imagefig,'WindowButtonUpFcn','segment(''pin_Buttonup'')');

    %Calculate distance
    [~,tempind] = min(sqrt((xm-y).^2+(ym-x).^2));
    DATA.Pin = tempind; %Store what pin

  case 'alt'
    DATA.contextmenu;
%     p = get(DATA.imagefig,'CurrentPoint');
%     switch type
%       case 'endo'
%         set(DATA.Handles.endopinmenu,...
%           'Position',p,...
%           'Visible','on');
%       case 'epi'
%         set(DATA.Handles.epipinmenu,...
%           'Position',p,...
%           'Visible','on');
%       case 'rvendo'
%         set(DATA.Handles.rvendopinmenu,...
%           'Position',p,...
%           'Visible','on');        
%     end;
end;

%------------------------------
function ok = enablecalculation %#ok<DEFNU>
%------------------------------
%Returns ok, if slices are selected. Sideeffect of
%this function is that it turns on interaction lock, stops movie. Calling
%this function will not allow user to click on a new function before
%endoffcalculation is called. Note that this mechanism needs to be pretty
%safe with try/catch clauses otherwise Segment may hang. Viewrefresh calls
%endoffcalculation.
global DATA SET NO

stopmovie_Callback;

myworkon;
if (not(isempty(SET(NO).StartSlice)))||(SET(NO).ZSize==1)
  if SET(NO).ZSize==1
    SET(NO).StartSlice = 1;
    SET(NO).EndSlice = 1;
  end;
  ok = true;
  DATA.Interactionlock = true;
else
  ok = false;
  myworkoff;
  mydisp('Warning: No slices selected => command ignored.');
end;

%-------------------------
function endoffcalculation
%-------------------------
%Turns off interaction lock, called after a calculation. See above.
global DATA 
DATA.Interactionlock = false;
myworkoff;

%-------------------------------
function viewrefreshall_Callback %#ok<DEFNU>
%-------------------------------
%Main refresh of GUI
global DATA SET NO

mydisp('Full screen refresh');
flushlog;

%To avoid if a button up function is not run.
set(DATA.fig,'WindowButtonMotionFcn','');

if DATA.Silent || (~DATA.DataLoaded)
  return;
end;

try
  oldNO = NO;

  for loop=1:length(SET)
    SET(loop).NormalZoomState = [];
    SET(loop).MontageZoomState = [];
    SET(loop).MontageRowZoomState = [];
  end;

  %---
  try
    delete(DATA.Handles.thumbnailimage);
  catch %#ok<CTCH>
  end

  try
    delete(DATA.Handles.thumbnaildragaxes);
  catch %#ok<CTCH>
  end
  %---
  checkconsistency; %if any problem with LV/RV seg, this should fix it
  updatemodeldisplay;
  drawfunctions('drawall',DATA.ViewMatrix);
  DATA.updatetitle;
  updatetool('select');
  endoffcalculation;
  DATA.switchtoimagestack(oldNO);
catch me
  disp('Could not perform refresh');
  mydispexception(me);
end;

flushlog;

%-----------------------------
function viewrefresh_Callback
%-----------------------------
%Main graphical refresh.
global DATA SET NO

set(DATA.Handles.refreshicon,'State','off');
set([...
  DATA.Handles.hideroiicon ...
  DATA.Handles.hidelvicon ...
  DATA.Handles.hidervicon ...  
  DATA.Handles.hidescaricon ...
  DATA.Handles.hidemaricon],'state','off');

SET(NO).NormalZoomState = [];
SET(NO).MontageZoomState = [];
SET(NO).MontageRowZoomState = [];

DATA.updatetitle;
drawfunctions('drawimageno');
updatevolume;
updatemodeldisplay;
updatetool('select');
endoffcalculation;

%--------------------------------------
function plotmodelrot_Callback(daz,del) %#ok<DEFNU>
%--------------------------------------
%Changes view of 3D model of the segmentation.
global DATA

f = gcf;
figure(4);
[az,el] = view;
view(az+daz,max(min(el+del,90),-90));
figure(f);

if DATA.Record
  drawnow;
  DATA.MovieFrame = mygetframe(4);
  export('exportmovierecorder_Callback','newframe');
end;

%------------------------------
function updatemmodevisibility %#ok<DEFNU>
%------------------------------
%Updates mmode visibility.
global DATA SET NO

rightstate = 'on';
leftstate = 'on';

if (SET(NO).TSize<2)||(isequal(DATA.ViewPanelsType{DATA.CurrentPanel},'flow'))
  rightstate = 'off';
end;

%Apply settings

%left line
set([...
  DATA.Handles.mmodeline ...
  DATA.Handles.mmodempoint1 ...
  DATA.Handles.mmodempoint2 ...
  DATA.Handles.mmodepoint1 ...
  DATA.Handles.mmodepoint2 ...
  DATA.Handles.mmodepointcenter],'visible',leftstate);

%right line
set([...
  DATA.Handles.mmode1line ...
  DATA.Handles.mmode2line ...
  DATA.Handles.mmodetimebar1 ...
  DATA.Handles.mmodetimebar2],'visible',rightstate);
%   DATA.Handles.distancetext

%------------------------
function updatemmodeline
%------------------------
%Updates the mmode line in mmode display.
global DATA SET NO

%obs reversed definiton of xy for mmode line - sorry!

set(DATA.Handles.mmode1line,'ydata',[SET(NO).XSize/2+SET(NO).Mmode.M1 SET(NO).XSize/2+SET(NO).Mmode.M1]);
set(DATA.Handles.mmode2line,'ydata',[SET(NO).XSize/2+SET(NO).Mmode.M2 SET(NO).XSize/2+SET(NO).Mmode.M2]);
[dist,timedist] = calcfunctions('calcmmodedists',NO);
set(DATA.Handles.mmodempoint1,...
  'xdata',SET(NO).Mmode.X-SET(NO).Mmode.Lx*SET(NO).Mmode.M1,...
  'ydata',SET(NO).Mmode.Y-SET(NO).Mmode.Ly*SET(NO).Mmode.M1);
set(DATA.Handles.mmodempoint2,...
  'xdata',SET(NO).Mmode.X-SET(NO).Mmode.Lx*SET(NO).Mmode.M2,...
  'ydata',SET(NO).Mmode.Y-SET(NO).Mmode.Ly*SET(NO).Mmode.M2);
set(DATA.Handles.mmodetimebar1,...
  'xdata',[1 1]*SET(NO).Mmode.T1);
set(DATA.Handles.mmodetimebar2,...
  'xdata',[1 1]*SET(NO).Mmode.T2);
%updates the text
DATA.updateaxestables('measure');
% set(DATA.Handles.distancetext,'string',...
%   sprintf('Distance:%3.2f [mm]\tTime:%3.2f [ms]',dist,timedist));

%------------------------------------
function updatemmode(arg,nbrofcycles) 
%------------------------------------
%Calculate and show mmode image
global DATA SET NO

if SET(NO).TSize<2
  return;
end;

onecycle = true;
if nargin == 2
  dozoom = true;
  switch nbrofcycles
    case 'one'
      onecycle = true;
    case 'three'
      onecycle = false;
  end
else
  dozoom = false;
end
if nargin == 0
  arg = [];
end

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

%Plot the line
set(DATA.Handles.mmodepoint1,'xdata',SET(NO).Mmode.X+SET(NO).Mmode.Lx*20,'ydata',SET(NO).Mmode.Y+SET(NO).Mmode.Ly*20);
set(DATA.Handles.mmodepoint2,'xdata',SET(NO).Mmode.X-SET(NO).Mmode.Lx*20,'ydata',SET(NO).Mmode.Y-SET(NO).Mmode.Ly*20);
set(DATA.Handles.mmodepointcenter,'xdata',SET(NO).Mmode.X,'ydata',SET(NO).Mmode.Y);

nl = SET(NO).XSize/2;
set(DATA.Handles.mmodeline,...
  'xdata',[SET(NO).Mmode.X+SET(NO).Mmode.Lx*nl SET(NO).Mmode.X-SET(NO).Mmode.Lx*nl],...
  'ydata',[SET(NO).Mmode.Y+SET(NO).Mmode.Ly*nl SET(NO).Mmode.Y-SET(NO).Mmode.Ly*nl]);
set(DATA.Handles.mmodempoint1,...
  'xdata',SET(NO).Mmode.X-SET(NO).Mmode.Lx*SET(NO).Mmode.M1,...
  'ydata',SET(NO).Mmode.Y-SET(NO).Mmode.Ly*SET(NO).Mmode.M1);
set(DATA.Handles.mmodempoint2,...
  'xdata',SET(NO).Mmode.X-SET(NO).Mmode.Lx*SET(NO).Mmode.M2,...
  'ydata',SET(NO).Mmode.Y-SET(NO).Mmode.Ly*SET(NO).Mmode.M2);

%Interpolate data
nl = SET(NO).XSize/2;
xi = linspace(SET(NO).Mmode.X+SET(NO).Mmode.Lx*nl,SET(NO).Mmode.X-SET(NO).Mmode.Lx*nl,SET(NO).XSize);
yi = linspace(SET(NO).Mmode.Y+SET(NO).Mmode.Ly*nl,SET(NO).Mmode.Y-SET(NO).Mmode.Ly*nl,SET(NO).XSize);

%Create temporal image
temp = zeros(length(xi),3*SET(NO).TSize);
if isempty(arg)
  for tloop=1:SET(NO).TSize
    temp(:,tloop + SET(NO).TSize*[0 1 2]) = repmat(interp2( ...
      SET(NO).IM(:,:,tloop,SET(NO).CurrentSlice),xi,yi)',1,3);
  end;
else
  for tloop=1:SET(NO).TSize
    temp(:,tloop + SET(NO).TSize*[0 1 2]) = repmat(interp2( ...
      SET(NO).IM(:,:,tloop,SET(NO).CurrentSlice),xi,yi,'nearest')',1,3);
  end;
end;

%Display data
temp = calcfunctions('remapuint8',temp,NO);

panel = [];
for loop=1:length(DATA.ViewPanelsType)
  if isequal(DATA.ViewPanelsType{loop},'mmodetemporal')
    panel = loop;
  end;
end;
if not(isempty(panel))
  set(DATA.Handles.imagehandle(panel),'cdata',temp);
  if dozoom
    imax = [0.5 1 0.5 SET(NO).XSize+0.5];
    if onecycle
      imax(2) = SET(NO).TSize+0.5;
    else
      imax(1) = -SET(NO).TSize+0.5;
      imax(2) = SET(NO).TSize*2+0.5;
    end
    axis(DATA.Handles.imageaxes(panel),imax);
  end
end;

%--------------------------
function mmodecenter_Motion %#ok<DEFNU>
%--------------------------
%Motion function of the mmode center point.
global SET NO

[x,y] = getclickedcoords;
SET(NO).Mmode.X = x;
SET(NO).Mmode.Y = y;
updatemmode('fast');

%---------------------
function mmode1_Motion %#ok<DEFNU>
%---------------------
%Motion function of the first mmode point.
global SET NO

[x,y] = getclickedcoords;
lx = (x-SET(NO).Mmode.X);
ly = (y-SET(NO).Mmode.Y);
len = sqrt(lx.^2+ly.^2);
lx = lx/len;
ly = ly/len;
SET(NO).Mmode.Lx = lx;
SET(NO).Mmode.Ly = ly;
updatemmode('fast');

%---------------------
function mmode2_Motion %#ok<DEFNU>
%---------------------
%Motion function of the second mmode point.
global SET NO

[x,y] = getclickedcoords;
lx = (SET(NO).Mmode.X-x);
ly = (SET(NO).Mmode.Y-y);
len = sqrt(lx.^2+ly.^2);
lx = lx/len;
ly = ly/len;
SET(NO).Mmode.Lx = lx;
SET(NO).Mmode.Ly = ly;
updatemmode('fast');

%-------------------------
function mmode1line_Motion %#ok<DEFNU>
%-------------------------
%Motion function of the first mmode line.
global SET NO

[~,y] = getclickedcoords;
SET(NO).Mmode.M1 = (y-SET(NO).XSize/2);
updatemmodeline

%-------------------------
function mmode2line_Motion %#ok<DEFNU>
%-------------------------
%Motion function of the second mmode line.
global SET NO

[~,y] = getclickedcoords;
SET(NO).Mmode.M2 = (y-SET(NO).XSize/2);
updatemmodeline;

%---------------------------
function mmodempoint1_Motion %#ok<DEFNU>
%---------------------------
%Motion function for mmodepoint one.
global SET NO

[x,y] = getclickedcoords;
SET(NO).Mmode.M1 = (SET(NO).Mmode.X-x)*SET(NO).Mmode.Lx+(SET(NO).Mmode.Y-y)*SET(NO).Mmode.Ly;
updatemmodeline;

%---------------------------
function mmodempoint2_Motion %#ok<DEFNU>
%---------------------------
%Motion function for mmodepoint two.
global SET NO

[x,y] = getclickedcoords;
SET(NO).Mmode.M2 = (SET(NO).Mmode.X-x)*SET(NO).Mmode.Lx+(SET(NO).Mmode.Y-y)*SET(NO).Mmode.Ly;
updatemmodeline;

%-----------------------------
function mmodetimebar1_Motion %#ok<DEFNU>
%-----------------------------
%Motion function for mmodepoint timebar
global SET NO

x = getclickedcoords;
SET(NO).Mmode.T1 = x;
updatemmodeline;

%-----------------------------
function mmodetimebar2_Motion %#ok<DEFNU>
%-----------------------------
%Motion function for mmodepoint timebar
global SET NO

x = getclickedcoords;
SET(NO).Mmode.T2 = x;
updatemmodeline;

%----------------------
function mmode_Buttonup %#ok<DEFNU>
%----------------------
%This function is called when buttonup occurs after draging one
%of the mmode objects.

global DATA 
%Restore so no motion is called
set(DATA.imagefig,'WindowButtonMotionFcn','');

%Restore main buttonup function
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''buttonup_Callback'')','segment'));

%Update mmode display
updatemmode;

%-------------------------------------
function mmodecenter_Buttondown(panel) %#ok<DEFNU>
%-------------------------------------
%Called when mmode center is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodecenter_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%--------------------------------
function mmode1_Buttondown(panel) %#ok<DEFNU>
%--------------------------------
%Called when mmode1 is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode1_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%--------------------------------
function mmode2_Buttondown(panel) %#ok<DEFNU>
%--------------------------------
%Called when mmode1 is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode2_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%------------------------------------
function mmode1line_Buttondown(panel) %#ok<DEFNU>
%------------------------------------
%Called when mmode1line is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode1line_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%------------------------------------
function mmode2line_Buttondown(panel) %#ok<DEFNU>
%------------------------------------
%Called when mmode1line is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode2line_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%-------------------------------------
function mmodepoint1_Buttondown(panel) %#ok<DEFNU>
%-------------------------------------
%Called when mmodepoint1 is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodempoint1_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%-------------------------------------
function mmodepoint2_Buttondown(panel) %#ok<DEFNU>
%-------------------------------------
%Called when mmodepoint1 is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodempoint2_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%---------------------------------------
function mmodetimebar1_Buttondown(panel) %#ok<DEFNU>
%---------------------------------------
%Called when mmodepoint1 is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodetimebar1_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%---------------------------------------
function mmodetimebar2_Buttondown(panel) %#ok<DEFNU>
%---------------------------------------
%Called when mmodepoint1 is pressed down, sets motion and buttonup
%function.
global DATA 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodetimebar2_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''mmode_Buttonup'')',mfilename));

%-------------------------
function setcurrenttimeframe(frame) %#ok<DEFNU>
%-------------------------
%function to set current time frame, input argument is frame to be set to
%current time frame
global SET NO

SET(NO).CurrentTimeFrame=frame;

nos = SET(NO).Linked;
for loop=1:length(nos(nos~=NO))
  SET(nos(loop)).CurrentTimeFrame = SET(NO).CurrentTimeFrame;
end;

for nloop=1:length(nos)
  drawfunctions('updatenopanels',nos(nloop));
end;

%updatevolume;
drawnow('expose'); %Expose does not other callbacks to evaluate.

%--------------------------
function nextframe_Callback %#ok<DEFNU>
%--------------------------
%Displays next timeframe of current image panel. Sideeffect is that the
%movie display is stopped if running.
global DATA SET NO
% 
% set([...
%   DATA.Handles.nextframeicon ...
%   DATA.Handles.playmovieicon ...
%   DATA.Handles.playallicon],'state','off');

stateandicon=segment('iconson','play');
stateandicon{2}.undent;
DATA.Handles.permanenticonholder.render;
DATA.Run = 0;

if SET(NO).TSize<=1
  return;
end;

SET(NO).CurrentTimeFrame = 1+mod(SET(NO).CurrentTimeFrame,SET(NO).TSize);

nos = SET(NO).Linked;
for loop=setdiff(nos,NO)
  SET(loop).CurrentTimeFrame = SET(NO).CurrentTimeFrame;
  if SET(loop).CurrentTimeFrame > SET(loop).TSize
    SET(loop).CurrentTimeFrame = SET(loop).TSize;
  end;
end;

for nloop=1:length(nos)
  drawfunctions('updatenopanels',nos(nloop));
end;

if DATA.Record
  DATA.MovieFrame = mygetframe(DATA.imagefig);
  export('exportmovierecorder_Callback','newframe');
end;

%updatevolume;
drawnow('expose'); %Expose does not other callbacks to evaluate.
  
%-----------------------------
function nextallframe_Callback
%-----------------------------
%Displays next timeframe in currrent image panel and adjust all visiable
%image stacks to the corresponding part of the cardiac cycle.
global DATA SET NO
% 
% set([...
%   DATA.Handles.nextallframeicon ...
%   DATA.Handles.playmovieicon ...
%   DATA.Handles.playallicon],'state','off');

stateandicon=segment('iconson','play');
stateandicon{2}.undent;
DATA.Handles.permanenticonholder.render;
DATA.Run = 0;

%if not timeresolved then return
if SET(NO).TSize<=1
  return;
end;

%Increase one
SET(NO).CurrentTimeFrame = SET(NO).CurrentTimeFrame+1;
if SET(NO).CurrentTimeFrame>SET(NO).TSize
  SET(NO).CurrentTimeFrame = 1;
end;

%Calculate percentage
percent = (SET(NO).CurrentTimeFrame-1)/(SET(NO).TSize-1);

%Loop over all image stacks to update time.
for loop=1:length(SET)
  SET(loop).CurrentTimeFrame = min(max(1,1+round(percent*(SET(loop).TSize-1))),SET(loop).TSize);
end;

%Draw all panels
nos = unique(DATA.ViewPanels(DATA.ViewPanels > 0));
for no=nos([SET(nos).TSize]>1)
  drawfunctions('updatenopanels',no);
end;

updatevolume;
if DATA.Record
  drawnow;
  DATA.MovieFrame = mygetframe(DATA.imagefig);
  export('exportmovierecorder_Callback','newframe');
end;
drawnow('expose'); %Expose does not other callbacks to evaluate.

%------------------------------
function previousframe_Callback %#ok<DEFNU>
%------------------------------
%Displays previous time frame of current panel. 
global DATA SET NO
% 
% set([...
%   DATA.Handles.previousframeicon ...
%   DATA.Handles.playmovieicon ...
%   DATA.Handles.playallicon],'state','off');

stateandicon=segment('iconson','play');
stateandicon{2}.undent;
DATA.Handles.permanenticonholder.render;
DATA.Run = 0;

if SET(NO).TSize<=1
  return;
end;

SET(NO).CurrentTimeFrame = 1+mod(SET(NO).CurrentTimeFrame-2,SET(NO).TSize);

nos = SET(NO).Linked;
for loop=setdiff(nos,NO)
  SET(loop).CurrentTimeFrame = SET(NO).CurrentTimeFrame;
  if SET(loop).CurrentTimeFrame > SET(loop).TSize
    SET(loop).CurrentTimeFrame = 1;
  end;
end;

for nloop=1:length(nos)
  drawfunctions('updatenopanels',nos(nloop));
end;

if DATA.Record
  drawnow;
  DATA.MovieFrame = mygetframe(DATA.imagefig);
  export('exportmovierecorder_Callback','newframe');
end;
drawnow('expose'); %Expose does not other callbacks to evaluate.

%-----------------------------
function previousallframe_Callback
%-----------------------------
%Displays previous time frame in current image panel. For all other visible
%image stacks they are adjusted to show correspondig part of the cardiac
%cycle.

global DATA SET NO
% 
% set([...
%   DATA.Handles.previousallframeicon ...
%   DATA.Handles.playmovieicon ...
%   DATA.Handles.playallicon],'state','off');
stateandicon=segment('iconson','play');
stateandicon{2}.undent;
DATA.Handles.permanenticonholder.render;

DATA.Run = 0;

%if not timeresolved then return
if SET(NO).TSize<=1
  return;
end;

%Decrease one
SET(NO).CurrentTimeFrame = SET(NO).CurrentTimeFrame-1;
if SET(NO).CurrentTimeFrame<1
  SET(NO).CurrentTimeFrame = SET(NO).TSize;
end;

%Calculate percentage
percent = (SET(NO).CurrentTimeFrame-1)/(SET(NO).TSize-1);

%Loop over all image stacks to update time.
for loop=1:length(SET)
  SET(loop).CurrentTimeFrame = min(max(1,1+round(percent*(SET(loop).TSize-1))),SET(loop).TSize);
end;

%Draw all panels
nos = unique(DATA.ViewPanels(DATA.ViewPanels > 0));
for no=nos([SET(nos).TSize]>1)
  drawfunctions('updatenopanels',no);
end;

updatevolume;
if DATA.Record
  drawnow;
  DATA.MovieFrame = mygetframe(DATA.imagefig);
  export('exportmovierecorder_Callback','newframe');
end;

drawnow('expose');%Expose does not other callbacks to evaluate.

%-----------------------------------------------------
function [xout,yout] = checkconsistencyhelper(xin,yin)
%-----------------------------------------------------
%Make sure that the contour is counter clock-wise and that it starts at
%three o clock. Also ensures that the points are evenly distributed.

global DATA

%--- Make sure evenly distributed
diffx = conv2(xin,[1;-1],'valid');
diffy = conv2(yin,[1;-1],'valid');
diffvec = sqrt(diffx.*diffx+diffy.*diffy);

if any((diffvec-mean(diffvec))/mean(diffvec)>2)
  mydisp('Problem with endo => fixed, report to support');
end;
contourlength = cumsum(diffvec);
xout = interp1(contourlength,xin(2:end),linspace(contourlength(1),contourlength(end),DATA.NumPoints-1));
yout = interp1(contourlength,yin(2:end),linspace(contourlength(1),contourlength(end),DATA.NumPoints-1));
xout(end+1) = xout(1);
yout(end+1) = yout(1);

%--- Make sure correct rotation (starting pos)
xr = yout;
yr = xout;
mx = mean(xr);
my = mean(yr);

if sum(unwrap(conv2(angle(complex(xr-mx,yr-my)),[1;-1],'valid')))<0
  disp(sprintf('counterclockwise endo %d slice %d',timeframe,zloop)); %#ok<DSPS>
  xr = fliplr(xr);
  yr = fliplr(yr);
end;

[~,inda] = min(angle(complex(mx-xr,my-yr)));
xout(1:(DATA.NumPoints-inda+1)) = yr(inda:end);
yout(1:(DATA.NumPoints-inda+1)) = xr(inda:end);
xout((DATA.NumPoints+1-inda):end) = yr(1:inda);
yout((DATA.NumPoints+1-inda):end) = xr(1:inda);
        
%------------------------------------------
function checkconsistency(timeframes,slice)
%------------------------------------------
%Check consistency, to prevent earlier
%manual segmentations that have problems with
%direction lef/right
global SET NO DATA

if nargin==0
  timeframes=SET(NO).CurrentTimeFrame; %1;
end;

if nargin<2
  slice=SET(NO).CurrentSlice;
end;

if not(isempty(SET(NO).EndoX)) && not(size(SET(NO).EndoX,1)==DATA.NumPoints)
  tempendox = nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);
  tempendoy = tempendox;
  timeframes = 1:SET(NO).TSize;
  slice = 1:SET(NO).ZSize;  
elseif isempty(SET(NO).EndoX)
  tempendox = []; % nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);
  tempendoy = []; % tempendox;
else
  tempendox = SET(NO).EndoX;
  tempendoy = SET(NO).EndoY;
end
if not(isempty(SET(NO).EpiX)) && not(size(SET(NO).EpiX,1)==DATA.NumPoints)
  tempepix = nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);
  tempepiy = tempepix;
  timeframes = 1:SET(NO).TSize;
  slice = 1:SET(NO).ZSize;  
elseif isempty(SET(NO).EpiX)
  tempepix = []; % nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);;
  tempepiy = []; % tempendox;
else
  tempepix = SET(NO).EpiX;
  tempepiy = SET(NO).EpiY;
end
if not(isempty(SET(NO).RVEndoX)) && not(size(SET(NO).RVEndoX,1)==DATA.NumPoints)
  temprvendox = nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);
  temprvendoy = tempendox;
  timeframes = 1:SET(NO).TSize;
  slice = 1:SET(NO).ZSize;  
elseif isempty(SET(NO).RVEndoX)
  temprvendox = []; % nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);
  temprvendoy = []; % tempendox;
else
  temprvendox = SET(NO).RVEndoX;
  temprvendoy = SET(NO).RVEndoY;
end
if not(isempty(SET(NO).RVEpiX)) && not(size(SET(NO).RVEpiX,1)==DATA.NumPoints)
  temprvepix = nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);
  temprvepiy = tempendox;
  timeframes = 1:SET(NO).TSize;
  slice = 1:SET(NO).ZSize;  
elseif isempty(SET(NO).RVEpiX)
  temprvepix = []; % nan(DATA.NumPoints,SET(NO).TSize,SET(NO).ZSize);
  temprvepiy = []; % tempendox;
else
  temprvepix = SET(NO).RVEpiX;
  temprvepiy = SET(NO).RVEpiY;
end
for timeframe = timeframes
  for zloop=slice;

    %Endo
    if ~isempty(SET(NO).EndoX)
      if not(isnan(SET(NO).EndoX(1,timeframe,zloop)))
        %Call consistenyhelper
        [tempendox(:,timeframe,zloop),tempendoy(:,timeframe,zloop)] = checkconsistencyhelper(...
          SET(NO).EndoX(:,timeframe,zloop),...
          SET(NO).EndoY(:,timeframe,zloop));          
      end;
    end;

    %Epi
    if ~isempty(SET(NO).EpiX)
      if not(isnan(SET(NO).EpiX(1,timeframe,zloop)))
        %Call consistenyhelper
        [tempepix(:,timeframe,zloop),tempepiy(:,timeframe,zloop)] = ...
          checkconsistencyhelper(...
          SET(NO).EpiX(:,timeframe,zloop),...
          SET(NO).EpiY(:,timeframe,zloop));       
      end;
    end;
    
    %RVEndo
    if ~isempty(SET(NO).RVEndoX)
      if not(isnan(SET(NO).RVEndoX(1,timeframe,zloop)))
        %Call consistenyhelper
        [temprvendox(:,timeframe,zloop),temprvendoy(:,timeframe,zloop)] = ...
          checkconsistencyhelper(...
          SET(NO).RVEndoX(:,timeframe,zloop),...
          SET(NO).RVEndoY(:,timeframe,zloop));                
      end;
    end;
    
    %RVEpi
    if ~isempty(SET(NO).RVEpiX)
      if not(isnan(SET(NO).RVEpiX(1,timeframe,zloop)))
        %Call consistenyhelper
        [temprvepix(:,timeframe,zloop),temprvepiy(:,timeframe,zloop)] = ...
          checkconsistencyhelper(...
          SET(NO).RVEpiX(:,timeframe,zloop),...
          SET(NO).RVEpiY(:,timeframe,zloop));                
      end;
    end;
    
  end;
end;
      
if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX = tempendox;
  SET(NO).EndoY = tempendoy;
end
if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX = tempepix;
  SET(NO).EpiY = tempepiy;
end
if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX = temprvendox;
  SET(NO).RVEndoY = temprvendoy;
end
if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX = temprvepix;
  SET(NO).RVEpiY = temprvepiy;
end

%---------------------------------------
function mask = createmask(outsize,y,x)
%---------------------------------------
%Function to generate a mask from a polygon represented with the vectors x
%and y. 
global DATA 

%Check if NaN then return empty
if any(isnan(x)) || any(isnan(y))
  mask = false(outsize);  
  return;
end;

if not(DATA.Pref.IncludeAllPixelsInRoi)
  %mask = roipoly(repmat(uint8(0),outsize),y,x);  
  mask = poly2mask(y,x,outsize(1),outsize(2));
else
  mask = false(outsize);  
  mx = round(mean(x));
  my = round(mean(y));
  x = interp1(x,linspace(1,length(x),1000));
  y = interp1(y,linspace(1,length(y),1000));
  mask(sub2ind(outsize,round(x),round(y))) = true;
  mask = imfill(mask,[mx my],4);
end;

%------------------------------------------------------
function z = reshape2layout(im,no,panel,outsideelement)
%------------------------------------------------------
%Convert a 3D array to an layout:ed image with cols, and rows
global DATA SET NO

if nargin<2
  no = NO;
end;

if nargin<3
  panel = DATA.CurrentPanel;
  if DATA.ViewPanels(panel) ~= no
    panel = find(DATA.ViewPanels == no,1);
  end
end

if nargin<4
  outsideelement = 0;
end;

z = repmat(outsideelement*im(1),DATA.ViewPanelsMatrix{panel}(1)*SET(no).XSize,DATA.ViewPanelsMatrix{panel}(2)*SET(no).YSize);
loop=1;
for slice=1:SET(no).ZSize
  c = 1+mod(loop-1,DATA.ViewPanelsMatrix{panel}(2));
  r = ceil(loop/DATA.ViewPanelsMatrix{panel}(2));
  z(...
    (1+(r-1)*SET(no).XSize):(r*SET(no).XSize),...
    (1+(c-1)*SET(no).YSize):(c*SET(no).YSize)) = im(:,:,slice);
  loop=loop+1;
end;

%----------------------------------------
function helpimagepos(im,sector,ofs,konst) %#ok<DEFNU>
%----------------------------------------
%helper function to display image.
global SET NO

im = im(sector,:,:)-ofs;
im = squeeze(im);
im = im(:,2:end);
h = image(konst*im);
set(h,'xdata',SET(NO).TIncr*[0.5 SET(NO).TSize+0.5]);
set(gca,'yticklabel','','xlim',[0 (SET(NO).TSize-1)*SET(NO).TIncr]);
xlabel('Time');
title(sprintf('Sector %d',sector));

%---------------------------------------
function helpimageposneg(im,sector,konst)
%---------------------------------------
%helper function to display image.
global SET NO

im = im(sector,:,:);
im = squeeze(im);
im = im(:,2:end);
h = image(32+konst*im);
set(h,'xdata',SET(NO).TIncr*[0.5 SET(NO).TSize+0.5]);
set(gca,'yticklabel','','xlim',[0 (SET(NO).TSize-1)*SET(NO).TIncr]);
xlabel('Time');
title(sprintf('Sector %d',sector));

%---------------------------
function resetlight_Callback %#ok<DEFNU>
%---------------------------
%Activated by toolbar icon, different from contrast_Callback (below).

global DATA SET NO

tools('enableundo',NO);
%set(DATA.Handles.resetlighticon,'state','off');
SET(NO).IntensityMapping.Contrast = 1;
SET(NO).IntensityMapping.Brightness = 0.5;
DATA.ViewIM{DATA.CurrentPanel} = [];
DATA.Overlay(DATA.CurrentPanel) = struct('alphadata',[],'cdata',[]);
update_thumbnail(NO);
makeviewim(DATA.CurrentPanel,NO);
drawfunctions('drawcontrastimage',NO); %drawfunctions('drawsliceno',NO);
if DATA.Pref.UseLight 
  DATA.BalloonLevel = -1; %Force update of ballonimage
end;

%---------------------------
function resetlightall_Callback %#ok<DEFNU>
%---------------------------
%Activated by toolbar icon, different from contrast_Callback (below).

global DATA SET NO

for no = 1:length(SET)
  SET(no).IntensityMapping.Contrast = 1;
  SET(no).IntensityMapping.Brightness = 0.5;
  update_thumbnail(no); %drawfunctions('drawsliceno',NO);
  
end
DATA.ViewIM{DATA.CurrentPanel} = [];
DATA.Overlay(DATA.CurrentPanel) = struct('alphadata',[],'cdata',[]);
makeviewim(DATA.CurrentPanel,NO);
drawfunctions('drawcontrastimage',NO);

if DATA.Pref.UseLight
  DATA.BalloonLevel = -1; %Force update of ballonimage
end;

%-----------------------------------
function contrast_Callback(arg,panel) 
%-----------------------------------
%Activated by contrast tool, different from resetlight_Callback (above).

global DATA SET NO
persistent handles

%if nargin==0
%  arg = 'init';
%end;

%switchtopanel(panel);

seltype = get(DATA.imagefig,'SelectionType');
switch seltype
  case 'alt'
    DATA.contextmenu;
  case {'normal','extend'}
    switch arg
%     case 'init'
%       set(DATA.Handles.imagehandle(panel),'ButtonDownFcn',...
%         sprintf('segment(''contrast_Callback'',''down'',%d)',panel));
%       %SET(NO).IntensityMapping.Contrast = 1; %Default values
%       %SET(NO).IntensityMapping.Brightness = 0.5; %Default values
%       %SET(NO).IntensityMapping.Compression = [];
      case 'down'
        switchtopanel(panel);
        tools('enableundo',NO);
        h = DATA.Handles.imageaxes(panel);
        [x,y] = mygetcurrentpoint(h);
        handles.xstart = x;
        handles.ystart = y;
        handles.xsize = get(h,'xlim');
        handles.xsize = handles.xsize(end)-handles.xsize(1);
        handles.ysize = get(h,'ylim');
        handles.ysize = handles.ysize(end)-handles.ysize(1);
        handles.deltacontrast = 0;
        handles.deltabrightness = 0;
        
        switch seltype
          case 'normal'
            %Button down, intialize
            set(gcf,'WindowButtonMotionFcn','segment(''contrast_Callback'',''motion'')');
            set(gcf,'WindowButtonUpFcn','segment(''contrast_Callback'',''up'')');
            handles.im = cell(1,numel(DATA.ViewPanels));
            switch DATA.ViewPanelsType{DATA.CurrentPanel}
              case {'one','ortho'}
                im = SET(NO).IM(:,:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice);
              case {'montage','montagerow','montagefit','sax3'}
                %Create space
                im = reshape2layout(...
                  squeeze(SET(NO).IM(:,:,SET(NO).CurrentTimeFrame,:)));
            end;
            handles.im{DATA.CurrentPanel} = im; %single(im);
%             %hide colormap until button up
%             disp('hej1');
%             handles.stateandiconcolorbar = segment('iconson','colorbar');
%             if handles.stateandiconcolorbar{1}
%               DATA.GUISettings.ShowColorbar = false;
%               drawfunctions('drawcolorbar');
%             end
          case 'extend'
            %Button down, intialize
            set(gcf,'WindowButtonMotionFcn','segment(''contrast_Callback'',''motionall'')');
            set(gcf,'WindowButtonUpFcn','segment(''contrast_Callback'',''upall'')');
            handles.im = cell(1,numel(DATA.ViewPanels));
            for i = 1:numel(DATA.ViewPanels)
              no = DATA.ViewPanels(i);
              if no > 0
                switch DATA.ViewPanelsType{i}
                  case {'one','ortho'}
                    im = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
                  case {'montage','montagerow','montagefit','sax3'}
                    %Create space
                    im = reshape2layout(...
                      squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:)),no,i);
                end;
                handles.im{i} = im; %single(im);
              end
            end
%             %hide colormap until button up
%             disp('hej2');
%             handles.stateandiconcolorbar = segment('iconson','colorbar');
%             if handles.stateandiconcolorbar{1}
%               DATA.GUISettings.ShowColorbar = false;
%               drawfunctions('drawcolorbar');
%             end
        end
      case {'up','upall'}
        panels = DATA.CurrentPanel;
        if strcmp(arg,'upall')
          panels = 1:numel(DATA.ViewPanels);
          panels = panels(DATA.ViewPanels > 0);
        end
        %restore
%         %show colormap again if it was shown before contrast setting
%         if handles.stateandiconcolorbar{1}
%           DATA.GUISettings.ShowColorbar = true;
%         end
        set(gcf,'WindowButtonMotionFcn','');
        set(gcf,'WindowButtonUpFcn','segment(''buttonup_Callback'')');
        precontrast = SET(NO).IntensityMapping.Contrast;
        prebrightness = SET(NO).IntensityMapping.Brightness;
        for panel = panels
          no = DATA.ViewPanels(panel);
          if ~isempty(SET(no).IntensityScaling) && ~isempty(SET(no).IntensityMapping)
            [window,level] = calcfunctions('con2win',...
              precontrast+handles.deltacontrast,...
              prebrightness+handles.deltabrightness,NO);
            [contrast,brightness] = calcfunctions('win2con',window,level,no);
          else
            contrast = precontrast+handles.deltacontrast;
            brightness = prebrightness+handles.deltabrightness;
          end
          SET(no).IntensityMapping.Contrast = contrast;
          SET(no).IntensityMapping.Brightness = brightness;
          DATA.ViewIM{panel} = [];
          DATA.Overlay(panel) = struct('alphadata',[],'cdata',[]);
          makeviewim(panel,no);
        end
        for no = DATA.ViewPanels(panels);
          drawfunctions('updatenopanels',no);
          update_thumbnail(no);
        end
        if DATA.Pref.UseLight
          DATA.BalloonLevel = -1; %Force update of ballonimage
        end;
      case 'wheelUpBrightness' %Brightness only 
        SET(NO).IntensityMapping.Brightness = SET(NO).IntensityMapping.Brightness+0.1;
        DATA.ViewIM{DATA.CurrentPanel} = [];
        DATA.Overlay(DATA.CurrentPanel) = struct('alphadata',[],'cdata',[]);
        makeviewim(DATA.CurrentPanel,NO);
        drawfunctions('updatenopanels',NO);
        update_thumbnail(NO);
        if DATA.Pref.UseLight
          DATA.BalloonLevel = -1; %Force update of ballonimage
        end;
      case 'wheelDownBrightness' %Brightness only
        SET(NO).IntensityMapping.Brightness = SET(NO).IntensityMapping.Brightness-0.1;
        DATA.ViewIM{DATA.CurrentPanel} = [];
        DATA.Overlay(DATA.CurrentPanel) = struct('alphadata',[],'cdata',[]);
        makeviewim(DATA.CurrentPanel,NO);
        drawfunctions('updatenopanels',NO);
        update_thumbnail(NO);
        if DATA.Pref.UseLight
          DATA.BalloonLevel = -1; %Force update of ballonimage
        end;
      case {'motion','motionall'}
        panels = DATA.CurrentPanel;
        if strcmp(arg,'motionall')
          panels = 1:numel(DATA.ViewPanels);
          panels = panels(DATA.ViewPanels > 0);
        end
        h = DATA.Handles.imageaxes(DATA.CurrentPanel);
        [x,y] = mygetcurrentpoint(h);
        handles.deltacontrast = (x-handles.xstart)/handles.xsize;
        handles.deltabrightness = (handles.ystart-y)/handles.ysize;
        for panel = panels
          no = DATA.ViewPanels(panel);
          if ~isempty(SET(no).IntensityScaling) && ~isempty(SET(no).IntensityMapping)
%             [window,level] = calcfunctions('con2win',...
%               SET(NO).IntensityMapping.Contrast+handles.deltacontrast,...
%               SET(NO).IntensityMapping.Brightness+handles.deltabrightness,NO);
            [window,level] = calcfunctions('con2win',...
              SET(no).IntensityMapping.Contrast+handles.deltacontrast,...
              SET(no).IntensityMapping.Brightness+handles.deltabrightness,no);
            [contrast,brightness] = calcfunctions('win2con',window,level,no);
          else
            %contrast = SET(NO).IntensityMapping.Contrast+handles.deltacontrast;
            %brightness = SET(NO).IntensityMapping.Brightness+handles.deltabrightness;
            contrast = SET(no).IntensityMapping.Contrast+handles.deltacontrast;
            brightness = SET(no).IntensityMapping.Brightness+handles.deltabrightness;
         
          end
          im = calcfunctions('remapuint8',handles.im{panel},no,...
            calcfunctions('returnmapping',no),contrast,brightness);
          set(DATA.Handles.imagehandle(panel),'cdata',im);
        end

      otherwise
        myfailed('Unknown option to contrast Callback',DATA.GUI.Segment);
    end;
end;
  

%----------------------------
  function [xlim,ylim] = getbox(no,destno,doindex)
%------------------------

global SET

  x=[];
  y=[];
      
    
  %Source is destination.
  if nargin == 1
    destno=no;
  end
  
  if nargin <3
    doindex=0;
  end
  xsz = SET(destno).XSize;
    ysz = SET(destno).YSize;
    
  diatype=1;
  hasanyseg=0;
  for type={'Endo','Epi','RVEndo','RVEpi'}
    tmp_x = SET(no).([type{1},'X']);
   tmp_y = SET(no).([type{1},'Y']);
    if ~isempty(tmp_x) && ~all(isnan(tmp_x(:)))
      hasanyseg=1;
      x=[x;tmp_x(:)];
      y=[y;tmp_y(:)];
      if ~isempty(regexpi(type{1},'RV'))
        diatype=2;
      end
    end
  end

  if hasanyseg
    %convert all segmentation rlapfh then get min max in destno
    rlapfh=calcfunctions('xyz2rlapfh',no,x,y,ones(length(x),1));
    limits=calcfunctions('rlapfh2xyz',destno,rlapfh(:,1),rlapfh(:,2),rlapfh(:,3));
    dodia=1;
  else
    
    diatype = 3;
    %use auto crop tool
    im= SET(no).IM;
    roisizetime=150;
    roisizenotime=200;
    if SET(no).TSize>1 && SET(no).ZSize~=1
      roisize = roisizetime;
    else
      roisize = roisizenotime;
    end;
    
    nx = roisize/SET(no).ResolutionX;
    ny = roisize/SET(no).ResolutionY;
    [~,~,rlapfh] = autocrop(im,nx,ny,0,80,1,no);
    if ~isempty(rlapfh)
      limits = calcfunctions('rlapfh2xyz',destno,rlapfh(:,1),rlapfh(:,2),rlapfh(:,3));
      dodia=1;
    else
      xlim=[1,xsz];
      ylim=[1,ysz];
      return;
    end
  end
    xmin=min(limits(1,:));
    ymin=min(limits(2,:));
    xmax=max(limits(1,:));
    ymax=max(limits(2,:));
    

    
    if xmax-xmin<20
    xmax=xsz;
    xmin=1;
    dodia=0;
    end
    
    if ymax-ymin<20
    ymax=ysz;
    ymin=1;
    dodia=0;
    end
    
    if dodia
      switch diatype
        case 1
          dia=norm([xmax,ymax]-[xmin,ymin])/2*0.8;
        case 2 %there is RV do very little margin addition
          dia=norm([xmax,ymax]-[xmin,ymin])/2*0.2;
        case 3 %from autocrop do "lagom" margin.
          dia=norm([xmax,ymax]-[xmin,ymin])/2*0.4;
      end
    else
      dia=0;
    end
    xlim=[xmin-dia,xmax+dia];
    ylim=[ymin-dia,ymax+dia];
    
    
    if doindex
      xlim=round(xlim);
      ylim=round(ylim);
      xlim(xlim>xsz)=xsz;
      ylim(ylim>ysz)=ysz;  
      xlim(xlim<1)=1;
      ylim(ylim<1)=1;
    end
   
%     if isempty(xlim)
%       xmin=1;
%       xmax=xsz;
%     else
%       xmin=xlim(1);
%       xmax=xlim(end);
%     end
%     
%     if isempty(ylim)
%       ymin=1;
%       ymax=ysz;
%     else
%       ymin=ylim(1);
%       ymax=ylim(end);
%     end
    

  
% %--------------------
%  function varargout = getbox(no,tono,indexes)
% %-----------------
% global SET
% 
% if nargin<3
%   indexes=0;
% end
%   xsz=SET(no).XSize;
%   ysz=SET(no).YSize;
%   
%   hasanyseg=0;
%   for type={'Endo','Epi','RVEndo','RVEpi'}
%     x=SET(no).([type{1},'X']);
%     if ~isempty(x) && ~all(isnan(x(:)))
%       hasanyseg=1;
%       break
%     end
%   end
%   
%   if hasanyseg
%     xmin=[];
%     ymin=[];
%     xmax=[];
%     ymax=[];
%     for type={'Endo','Epi','RVEndo','RVEpi'}
%       
%       x = SET(no).([type{1},'X']);
%       y = SET(no).([type{1},'Y']);
%       
%       xmin = min([xmin floor(min(x(:)))]);
%       ymin = min([ymin floor(min(y(:)))]);
%       xmax = max([xmax ceil(max(x(:)))]);
%       ymax = max([ymax ceil(max(y(:)))]);
%       
%       dodia = 1;%round(norm([xmax,ymax]-[xmin,ymin])/2*0.8); 
%       
%       if xmin<1
%         xmin=1;
%         dodia=0;
%       end
%       
%       if ymin<1
%         ymin=1;
%         dodia=0;
%       end
%       
%       if xmax>xsz
%         xmax=xsz;
%         dodia=0;
%       end
%       
%       if ymax>ysz
%         ymax=ysz;
%         dodia=0;
%       end
%       
%       %this catches segmentation blow ups
%       if ymax-ymin<20
%         ymin=1;
%         ymax=ysz;
%         dodia=0;
%       end
%       
%       if xmax-xmin<20
%         xmin=1;
%         xmax=xsz;
%         dodia=0;
%       end
%       
%     end
%     
%   else
%   if nargin == 2
%     no=tono;
%     xsz=SET(no).XSize;
%     ysz=SET(no).YSize;
%   end
%     
%    im= SET(no).IM;
%     roisizetime=150;
%      roisizenotime=200;
%      if SET(no).TSize>1 && SET(no).ZSize~=1
%        roisize = roisizetime;
%      else
%        roisize = roisizenotime;
%      end;
%      
%   nx = roisize/SET(no).ResolutionX;
%   ny = roisize/SET(no).ResolutionY;
%     [xlim,ylim,rlapfh] = autocrop(im,nx,ny,0,128,1,no);
%     dodia=0;
%     
%     if isempty(xlim)
%       xmin=1;
%       xmax=xsz;
%     else
%     xmin=xlim(1);
%     xmax=xlim(end);
%     end
%     
%     if isempty(ylim)
%       ymin=1;
%       ymax=ysz;
%     else
%     ymin=ylim(1);
%     ymax=ylim(end);  
%     end
%   end
% %   xlim=[xmin-dia,xmax+dia];
% %   ylim=[ymin-dia,ymax+dia];
%   xlim=[xmin,xmax];
%   ylim=[ymin,ymax];
%   %SET(no).NormalZoomState=[ymin-dia,ymax+dia,xmin-dia,xmax+dia]; 
%   
%    rlapfh = calcfunctions('xyz2rlapfh',no,xlim',ylim',[1,1]');
%   
%   if nargin ==1 && nargout>1
%     if dodia
%       dia=norm([xlim(end),ylim(end)]-[xlim(1),ylim(1)])/2*0.8;
%     else
%       dia=0;
%     end
%     xlim=[xlim(1)-dia,xlim(end)+dia];
%     ylim=[ylim(1)-dia,ylim(end)+dia];
%     
%     if indexes
%       xlim(xlim<1)=1;
%       ylim(ylim<1)=1;
%       
%       xlim(xlim>SET(no).XSize)=SET(no).XSize;
%       ylim(ylim>SET(no).YSize)=SET(no).YSize;
%     end
%   end 
%   
%   if nargin>1 && nargout>1
%     limits = calcfunctions('rlapfh2xyz',tono,rlapfh(:,1),rlapfh(:,2),rlapfh(:,3));
%     
%     %sortlimits
%     xlim=sort(limits(1,:),'ascend');
%     ylim=sort(limits(2,:),'ascend');
%     
%     if dodia
%       dia=norm([xlim(end),ylim(end)]-[xlim(1),ylim(1)])/2*0.8;
%     else
%       dia=0;
%     end
%     xlim=[xlim(1)-dia,xlim(end)+dia];
%     ylim=[ylim(1)-dia,ylim(end)+dia];
%     
%     
%     if indexes     
%     xlim(xlim<1)=1;
%     ylim(ylim<1)=1;
%     
%     xlim(xlim>SET(tono).XSize)=SET(tono).XSize;
%     ylim(ylim>SET(tono).YSize)=SET(tono).YSize;
%     end
%   end
%   
%   if nargout==1
%     varargout={rlapfh};
%   end
%   
%   if nargout==2
%     varargout={round(xlim),round(ylim)};
%   end
%   
%   if nargout==3
%     varargout={rlapfh,round(xlim),round(ylim)};
%   end

  
%-------------------------
  function varargout = autozoom
%-------------------------
% Autozooms everything that is displayed.
global DATA SET

% if isempty(DATA.ZoomState) || DATA.ZoomState == 0
%   DATA.ZoomState = 1;
% elseif DATA.ZoomState==2
%   DATA.ZoomState=0;
% end

no_inds=find(DATA.ViewPanels~=0);
no_inds(~strcmp('one',DATA.ViewPanelsType(no_inds)))=[];

lvno = findfunctions('findcineshortaxisno');

if ~isempty(lvno)
imageorientation = SET(lvno).ImageOrientation;
refaxis = cross(imageorientation(1:3),imageorientation(4:6));
tol=11;
fixednos = false(size(no_inds));

for i=1:length(no_inds)
  no=DATA.ViewPanels(no_inds(i));
  hasanyseg=0;
  for type={'Endo','Epi','RVEndo','RVEpi'}
    x=SET(no).([type{1},'X']);
    if ~isempty(x) && ~all(isnan(x(:)))
      hasanyseg=1;
      break
    end
  end
  
  if hasanyseg && no~=lvno
    break
  end
  
  imageorientation = SET(no).ImageOrientation;
  axis = cross(imageorientation(1:3),imageorientation(4:6));
  score = abs(sum(refaxis.*axis)); %scalar product, and abs
  if acos(score)*180/pi<tol
    [xlim,ylim] = getbox(lvno,no);%,1);
    SET(no).NormalZoomState=[ylim, xlim];
    fixednos(i)=1;
  set(DATA.Handles.imageaxes(no_inds(i)),'xlim',ylim,'ylim',xlim)
  drawfunctions('viewupdatetextposition',no_inds(i));
  drawfunctions('viewupdateannotext',no_inds(i));
  end
end
no_inds(fixednos)=[];
end
%find images with similar imageposition
%group ={};
% tol=0.05;
% notmp=no_inds;
% while  ~isempty(notmp) 
%   tmp=notmp(1);
%   matches=cellfun(@(x) norm(SET(DATA.ViewPanels(tmp)).ImagePosition-x),{SET(DATA.ViewPanels(notmp)).ImagePosition})<tol;
%   if any(matches)
%     group=[group, notmp(matches)];
%     notmp(matches)=[];
%   end
% end

for i=no_inds
  no=DATA.ViewPanels(i);
  [xlim,ylim]=getbox(no);
  %limits = calcfunctions('rlapfh2xyz',no,rlapfh(:,1),rlapfh(:,2),rlapfh(:,3));
  %SET(no).NormalZoomState=[limits(2,:), limits(1,:)];
  SET(no).NormalZoomState=[ylim, xlim];
  set(DATA.Handles.imageaxes(i),'xlim',ylim,'ylim',xlim)
  drawfunctions('viewupdatetextposition',i);
  drawfunctions('viewupdateannotext',i);
end
%   xsz=SET(no).XSize;
%   ysz=SET(no).YSize;
%   
%   hasanyseg=0;
%   for type={'Endo','Epi','RVEndo','RVEpi'}
%     x=SET(no).([type{1},'X']);
%     if ~isempty(x) && ~all(isnan(x(:)))
%       hasanyseg=1;
%       break
%     end
%   end
%   
%    if DATA.Pref.ViewInterpolated
%     scale=2;
%   else
%     scale=1;
%   end
%   
%   if hasanyseg
%     xmin=[];
%     ymin=[];
%     xmax=[];
%     ymax=[];
%     for type={'Endo','Epi','RVEndo','RVEpi'}
%       
%       x = SET(no).([type{1},'X']);
%       y = SET(no).([type{1},'Y']);
%       
%       xmin = min([xmin floor(min(x(:)))]);
%       ymin = min([ymin floor(min(y(:)))]);
%       xmax = max([xmax ceil(max(x(:)))]);
%       ymax = max([ymax ceil(max(y(:)))]);
%       
%       dia = round(norm([xmax,ymax]-[xmin,ymin])/2*0.75); 
%       
%       if xmin<1
%         xmin=1;
%         dia=0;
%       end
%       
%       if ymin<1
%         ymin=1;
%         dia=0;
%       end
%       
%       if xmax>xsz
%         xmax=xsz;
%         dia=0;
%       end
%       
%       if ymax>ysz
%         ymax=ysz;
%         dia=0;
%       end
%       
%       %this catches segmentation blow ups
%       if ymax-ymin<20
%         ymin=1;
%         ymax=ysz;
%         dia=0;
%       end
%       
%       if xmax-xmin<20
%         xmin=1;
%         xmax=xsz;
%         dia=0;
%       end
%       
%     end
%     
%   else
%    im= SET(no).IM;
%     roisizestacks=200;
%      roisizeslices=250;
%      if SET(no).ZSize>1
%        roisize = roisizestacks;
%      else
%        roisize = roisizeslices;
%      end;
%   
%   nx = roisize/SET(no).ResolutionX;
%   ny = roisize/SET(no).ResolutionY;
%     [xlim,ylim,varargout{1}] = autocrop(im,nx,ny,0,128,13,no);
%     dia=0;
%     
%     if isempty(xlim)
%       xmin=1;
%       xmax=xsz;
%     else
%     xmin=xlim(1);
%     xmax=xlim(end);
%     end
%     
%     if isempty(ylim)
%       ymin=1;
%       ymax=ysz;
%     else
%     ymin=ylim(1);
%     ymax=ylim(end);  
%     end
%   end
%   SET(no).NormalZoomState=[ymin-dia,ymax+dia,xmin-dia,xmax+dia];
%  set(DATA.Handles.imageaxes(i),'xlim',[ymin-dia,ymax+dia],'ylim',[xmin-dia,xmax+dia])
%   
% drawfunctions('viewupdatetextposition',i);
% drawfunctions('viewupdateannotext',i);
%end

  

  
% Plan B
%   if ymax-ymin>=xmax-xmin
%     f = 120/(ymax-ymin+1);
%   else
%     f = 120/(xmax-xmin+1);
%   end
%   
%   zoomhelper(no,f);
%end


%------------------------------
function autocontrastall_Callback %#ok<DEFNU>
%------------------------------
%Automatically calculates contrast settings.
global SET DATA
for i = 1:length(SET)
  autocontrast(i,1);
  segment('update_thumbnail',i)
end

for panelloop = 1:length(DATA.ViewPanels)
  panel = DATA.ViewPanels(panelloop);
  drawfunctions('drawcontrastimage',panel);
end

%------------------------------
function autocontrast_Callback %#ok<DEFNU>
%------------------------------
%Automatically calculates contrast settings.
global NO

autocontrast(NO);

%-------------------------
function autocontrast(no,silent)
%-------------------------
%Helper functionk to autocontrast_Callback_
global SET

if nargin <2
  silent=0;
end

%check so that no isn't a phase image
if isfield(SET(no).Flow, 'PhaseNo') && SET(no).Flow.PhaseNo==no
  return
end

%set percentile values
lowerpercentile=0.02;
upperpercentile=0.99;

%calulate histogram for the whole image
myworkon;

sortim = sort(SET(no).IM(:));
lowerthresh = sortim(round(lowerpercentile*length(sortim)));
upperthresh = sortim(round(upperpercentile*length(sortim)));
minInt = sortim(1);
maxInt = sortim(end);

if isempty(lowerthresh)
  lowerthresh=minInt;
end
if isempty(upperthresh)
  upperthresh=maxInt;
end  

%calcualte contrast and brightness so that the intensities below
%lowerthresh and above upperthresh gets saturated
contrast=(maxInt-minInt)/(upperthresh-lowerthresh);
brightness=(maxInt+minInt)/2-(upperthresh+lowerthresh)/2*contrast+0.5;%(maxInt+minInt)/2-(upperthresh-lowerthresh)/2+0.5;%lowerthresh-minInt+0.5;

SET(no).IntensityMapping.Contrast=contrast;
SET(no).IntensityMapping.Brightness=brightness;

%update image
if ~silent 
  drawfunctions('drawcontrastimage',no);
end
myworkoff;

%----------------------------------------
function placetimeresolvedpoints_Callback %#ok<DEFNU>
%----------------------------------------
%When this option is enabled then points are placed timeresolved
%and the name is re-used.

global DATA

%Updating in menu. When checked it changes how points are placed.
if isequal(get(DATA.Handles.placetimeresolvedpoints,'checked'),'off')
  set(DATA.Handles.placetimeresolvedpoints,'checked','on');
else
  set(DATA.Handles.placetimeresolvedpoints,'checked','off');
end;

%-----------------------------------
function measuremove_Callback(dx,dy) %#ok<DEFNU>
%-----------------------------------
%Helper function to move measurements.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

n = DATA.MeasureN;
if (length(SET(no).Measure)<n)||(n<1)
  %Invalid measurement 
  return;
end;
tools('enableundo',no);

SET(no).Measure(n).X = SET(no).Measure(n).X+dx;
SET(no).Measure(n).Y = SET(no).Measure(n).Y+dy;

drawfunctions('drawimageno');

%------------------------------
function measureexport_Callback %#ok<DEFNU>
%------------------------------
%Export measurements.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

if ~isempty(SET(no).Measure)
  stri = [];
  stri = [stri ...
    sprintf('PatientName:\t%s\n',SET(no).PatientInfo.Name) ...
    sprintf('\n') ...
    sprintf('Name\tLength[mm]\n')];
  for loop=1:length(SET(no).Measure)
    stri = [stri ...
      sprintf('%s\t%f\n',...
      SET(no).Measure(loop).Name,...
      SET(no).Measure(loop).Length)]; %#ok<AGROW>
  end;

  clipboard('copy',stri);
  mymsgbox('Results copied to clipboard','Done!',DATA.GUI.Segment);
else
  myfailed('Nothing to export.',DATA.GUI.Segment);
end;

%-----------------------------------
function measureshapeexport_Callback %#ok<DEFNU>
%-----------------------------------
%Export measurement shape
global DATA SET NO
%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

if ~isempty(SET(no).Measure)
  c = cell(6,3*numel(SET(no).Measure));
  c{1,1} = 'PatientName';
  c{1,2} = SET(no).PatientInfo.Name;
  for mnbr = 1:numel(SET(no).Measure)
    c{2,mnbr*3-2} = 'Name';
    c{2,mnbr*3-1} = SET(no).Measure(mnbr).LongName;
    c{3,mnbr*3-2} = 'Length [mm]';
    c{3,mnbr*3-1} = SET(no).Measure(mnbr).Length;
    c(4,mnbr*3+(-2:0)) = {'X [mm]' 'Y [mm]' 'Z [mm]'};
    basepoint = [SET(no).Measure(mnbr).X(1)*SET(no).ResolutionX, ...
      SET(no).Measure(mnbr).Y(1)*SET(no).ResolutionY, ...
      -SET(no).Measure(mnbr).Z(1)*(SET(no).SliceGap + SET(no).SliceThickness)];
    c(5:4+numel(SET(no).Measure(mnbr).X),mnbr*3-2) = num2cell( ...
      [SET(no).Measure(mnbr).X]*SET(no).ResolutionX-basepoint(1));
    c(5:4+numel(SET(no).Measure(mnbr).X),mnbr*3-1) = num2cell( ...
      [SET(no).Measure(mnbr).Y]*SET(no).ResolutionY-basepoint(2));
    c(5:4+numel(SET(no).Measure(mnbr).X),mnbr*3) = num2cell( ...
      -[SET(no).Measure(mnbr).Z]*(SET(no).SliceGap + SET(no).SliceThickness)-basepoint(3));
%     stri = [sprintf('PatientName:\t%s\n',SET(no).PatientInfo.Name) ...
%       sprintf('Name\t%s\nLength[mm]\t%f\nX\tY\tZ\n',SET(no).Measure(mnbr).LongName,...
%       SET(no).Measure(mnbr).Length)];
%     basepoint = [SET(no).Measure(mnbr).X(1)*SET(no).ResolutionX, ...
%       SET(no).Measure(mnbr).Y(1)*SET(no).ResolutionY, ...
%       -SET(no).Measure(mnbr).Z(1)*(SET(no).SliceGap + SET(no).SliceThickness)];
%     for loop=1:length(SET(no).Measure(mnbr).X)
%       stri = [stri sprintf('%f\t%f\t%f\n', ...
%         SET(no).Measure(mnbr).X(loop)*SET(no).ResolutionX-basepoint(1), ...
%         SET(no).Measure(mnbr).Y(loop)*SET(no).ResolutionY-basepoint(2), ...
%         -SET(no).Measure(mnbr).Z(loop)* ...
%         (SET(no).SliceGap + SET(no).SliceThickness)- basepoint(3))]; %#ok<AGROW>
  end
  cell2clipboard(c);
else
  myfailed('Nothing to export.',DATA.GUI.Segment);
end;

%--------------------------------------
function measurepoint_Buttondown(panel) %#ok<DEFNU>
%--------------------------------------
%Buttondown function when clicking on a measurement point/marker.
global DATA SET NO

killbuttondown = switchtopanel(panel);

if killbuttondown
  return
end

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

if not(isfield(SET(no),'Measure'))
  SET(no).Measure = [];
end;

%Get clicked coordinate
[x,y,slice] = getclickedcoords;
% If slice has changed, make sure montage/one are in sync
if (slice>SET(no).ZSize)
  return;
end
if ~ismember(DATA.ViewPanelsType{panel},{'hla','vla','gla'})
  switchtoslice(slice);
end
type = get(DATA.imagefig,'SelectionType');

%Find correct point
ind = NaN;
mindist = 1e10;
[measure, slice] = segment('getmeasurecoords',no,panel);

%Find closest point: 
for loop=1:length(measure)
  if slice<=max(measure(loop).Z) && slice>=min(measure(loop).Z)
    [dist,num] = min(sqrt(...
      (measure(loop).X-y).^2+...
      (measure(loop).Y-x).^2));
    
    if dist<mindist
      ind = loop;
      mindist = dist;
      indnum = num;
    end;
  end;
end;

if isnan(ind)
  myfailed('Could not find measurement point, should not occur.',DATA.GUI.Segment);
  return;
end;

%Prepare DATA.Measure info
DATA.MeasureN = ind;
DATA.MeasureX = SET(no).Measure(ind).X;
DATA.MeasureY = SET(no).Measure(ind).Y;
DATA.MeasureZ = SET(no).Measure(ind).Z;
DATA.MeasureName = SET(no).Measure(ind).Name;
DATA.MeasureT = SET(no).Measure(ind).T;
[DATA.MeasureOffsetY, DATA.MeasureOffsetX] = calcfunctions('calcoffset',slice,[],no,panel);

translateon = strcmp(DATA.CurrentTool,'move'); %SBT20160308

switch type
  case 'alt'
    DATA.contextmenu;
    otherwise
      if translateon %SBT20160308
         %Reset motion and store current mouse position: 
          measure_Motion_translate('reset');
          set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''measure_Motion_translate'');',mfilename)); %SBT20160308
      else %SBT20160308
          %Call measure motion %SBT20160308
          set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''measure_Motion'',%d);',mfilename,indnum)); %SBT20160308
      end %SBT20160308
      set(DATA.imagefig,'WindowButtonUpFcn',... %SBT20160308
              sprintf('%s(''measure_Buttonup'')','segment')); %SBT20160308
end;

%----------------------
function measure_Motion_translate(reset) %#ok<DEFNU> SBT20160308
%----------------------
%Motion function for measurements.
persistent startx starty startslice 
global DATA SET NO

if (nargin==1)
  %Reset coordinates
  [startx,starty,startslice] = getclickedcoords;
else
    no = NO;
    if ~isempty(SET(NO).Parent)
        no = SET(NO).Parent;
    end;
    
    %Compute translational difference:
    [x,y,slice] = getclickedcoords;
    
    %Abort upon slice-change: 
    if (slice~=startslice)
        DATA.measure_Buttonup
        return;
    end
    
    % Only update if marker is inside image: 
    if (x>0.5) && (y>0.5) && (x<SET(NO).YSize+0.5) && (y<SET(NO).XSize+0.5)
        
        DATA.MeasureX = DATA.MeasureX+(y-starty);     %Change to translate
        DATA.MeasureY = DATA.MeasureY+(x-startx);     %Change to translate
%         DATA.MeasureZ = DATA.MeasureZ+(slice-startslice); %Change to translate
        
        %Reset startingpoint for next iteration: 
        starty = y;
        startx = x;
        startslice = slice;
        
        set(DATA.Handles.measureline{DATA.CurrentPanel}{DATA.MeasureN},...
            'xdata',DATA.MeasureOffsetX+DATA.MeasureY,...
            'ydata',DATA.MeasureOffsetY+DATA.MeasureX);
        [ymax,ix] = max(DATA.MeasureOffsetX+DATA.MeasureY);
        set(DATA.Handles.measuretext{DATA.CurrentPanel}(DATA.MeasureN),'position',[...
            ymax+1 ...
            DATA.MeasureX(ix) ...
            0]);
        if isequal(get(DATA.Handles.hidetexticon,'state'),'off') &&...
                (DATA.MeasureN <= length(SET(no).Measure))
            dist = sum(sqrt(...
                (SET(no).ResolutionX*diff(DATA.MeasureX)).^2+...
                (SET(no).ResolutionY*diff(DATA.MeasureY)).^2));
            set(DATA.Handles.measuretext{DATA.CurrentPanel}(DATA.MeasureN),...
                'string',sprintf('%s\n%0.1f [mm]',SET(no).Measure(DATA.MeasureN).Name,dist));
        end
        
    end
end

%----------------------
function measure_Motion(ind) %#ok<DEFNU>
%----------------------
%Motion function for measurements.
global DATA SET NO

if nargin < 1
  ind = numel(DATA.MeasureX);
end

no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

[x,y,slice] = getclickedcoords;

% first line is for montage
% second line is for one-view.
if (x>0.5) && (y>0.5) && (x<SET(NO).YSize+0.5) && (y<SET(NO).XSize+0.5)
 
  DATA.MeasureX(ind) = y;
  DATA.MeasureY(ind) = x;
  DATA.MeasureZ(ind) = slice;
  
  set(DATA.Handles.measureline{DATA.CurrentPanel}{DATA.MeasureN},...
    'xdata',DATA.MeasureOffsetX+DATA.MeasureY,...
    'ydata',DATA.MeasureOffsetY+DATA.MeasureX);
  [ymax,ix] = max(DATA.MeasureOffsetX+DATA.MeasureY);
  set(DATA.Handles.measuretext{DATA.CurrentPanel}(DATA.MeasureN),'position',[...
    ymax+1 ...
    DATA.MeasureX(ix) ...
    0]);
  if isequal(get(DATA.Handles.hidetexticon,'state'),'off') &&...
      (DATA.MeasureN <= length(SET(no).Measure))
    dist = sum(sqrt(...
      (SET(no).ResolutionX*diff(DATA.MeasureX)).^2+...
      (SET(no).ResolutionY*diff(DATA.MeasureY)).^2));
    set(DATA.Handles.measuretext{DATA.CurrentPanel}(DATA.MeasureN),...
      'string',sprintf('%s\n%0.1f [mm]',SET(no).Measure(DATA.MeasureN).Name,dist));
  end

end

%---------------------------------
function measure_Buttondown(panel) %#ok<DEFNU>
%---------------------------------
%Button down function for placing measurements.
global DATA SET NO

killbuttondown = switchtopanel(panel);

if killbuttondown
  return
end

if DATA.Interactionlock
  return;
end;

set(DATA.Handles.hidemeasuresicon,'state','off');

seltype = get(DATA.imagefig,'SelectionType');
switch seltype
  case {'normal','extend'}
    
    %Use to point to mag data set
    no = NO;
    if ~isempty(SET(NO).Parent)
      no = SET(NO).Parent;
    end;
    
    [x,y,slice] = getclickedcoords;
    if ismember(DATA.ViewPanelsType{panel},{'montage','montagerow','montagefit'})
      % If slice has changed, make sure montage/one are in sync
      if (slice>SET(no).ZSize)
        %         ~(strcmp(DATA.ViewPanelsType{panel},'hla') && slice <= SET(no).HLA.maxslice || ...
        %         strcmp(DATA.ViewPanelsType{panel},'vla') && slice <= SET(no).VLA.maxslice)
        return;
      end
      switchtoslice(slice);
    end
    
    DATA.MeasureT = SET(NO).CurrentTimeFrame;
    DATA.MeasureN = length(SET(no).Measure)+1;
    DATA.MeasureName = '';
    
    if strcmp(DATA.ViewPanelsType{panel},'mmodetemporal')
      myfailed('Can not measure in time.',DATA.GUI.Segment);
      return
    end
    [DATA.MeasureOffsetY,DATA.MeasureOffsetX] = calcfunctions('calcoffset',slice,[],no,panel);
    
    set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''measure_Motion'');',mfilename));
    if strcmp(seltype,'normal')
      DATA.MeasureX = [y;y];
      DATA.MeasureY = [x;x];
      DATA.MeasureZ = [slice;slice];
      set(DATA.imagefig,'WindowButtonUpFcn',...
        sprintf('%s(''measure_Buttonup'')','segment'));
    else
      DATA.MeasureX = y;
      DATA.MeasureY = x;
      DATA.MeasureZ = slice;
      set(DATA.imagefig,'WindowButtonUpFcn',...
        sprintf('%s(''measureput_Buttonup'')',mfilename));
    end
    
    hold(DATA.Handles.imageaxes(DATA.CurrentPanel),'on');
    DATA.Handles.measureline{DATA.CurrentPanel}{DATA.MeasureN} = plot(...
      DATA.Handles.imageaxes(DATA.CurrentPanel),DATA.MeasureY,DATA.MeasureX,DATA.GUISettings.MeasureLineSpec);
    set(DATA.Handles.measureline{panel}{DATA.MeasureN},'markersize',DATA.GUISettings.MeasureLineMarkerSize);
    DATA.Handles.measuretext{DATA.CurrentPanel}(DATA.MeasureN) = text(...
      mean(DATA.MeasureY),mean(DATA.MeasureX),'', ...
      'Parent',DATA.Handles.imageaxes(DATA.CurrentPanel));
    hold(DATA.Handles.imageaxes(DATA.CurrentPanel),'off');
  case 'alt'
    DATA.contextmenu;
end

%---------------------------
function measureput_Buttonup %#ok<DEFNU>
%---------------------------
%Called when a measurement point (except the endpoint) is placed
global DATA

switch get(DATA.imagefig,'SelectionType')
  case 'extend'
    DATA.MeasureX(end+1,1) = DATA.MeasureX(end);
    DATA.MeasureY(end+1,1) = DATA.MeasureY(end);
    DATA.MeasureZ(end+1,1) = DATA.MeasureZ(end);
  case 'normal'
    DATA.measure_Buttonup;
end


%--------------------------------
function normal_Buttondown(panel) %#ok<DEFNU>
%--------------------------------
%Called when button down in normal view
global DATA SET 

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

switch get(DATA.imagefig,'SelectionType')
  case 'extend' %shift
    %pan_Buttondown
    contrast_Callback('down',panel) 
  case 'alt'
    DATA.contextmenu;
  case 'open'
    %Double click
    no = DATA.ViewPanels(DATA.CurrentPanel);
    if numel(DATA.ViewPanels) > 1
      DATA.LastView = struct(...
        'Panels',DATA.ViewPanels, ...
        'PanelsType',{DATA.ViewPanelsType}, ...
        'Matrix',DATA.ViewMatrix);
      drawfunctions('drawall',1);
    elseif ~isempty(DATA.LastView) && ismember(no,DATA.LastView.Panels)
      lastview = DATA.LastView;
      DATA.LastView = [];
      panels = lastview.Panels;
      DATA.CurrentPanel = find(panels==no,1);
      drawfunctions('drawimageview',panels, ...
        lastview.Matrix,lastview.PanelsType);
    elseif SET(no).ZSize>1
      viewimage_Callback('montage');
    end;   
end

set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''buttonup_Callback'')','segment'));

%----------------------------
function makeviewim(panel,no)
%----------------------------
%Rearrange data to show all slices
global DATA SET NO

myworkon;

if nargin==0
  panel = DATA.CurrentPanel;
end;

if nargin<2
  %This clause is needed when called from GUI.
  no = NO;
end;

if isempty(SET(no).CurrentSlice)
  SET(no).CurrentSlice=1;
  SET(no).StartSlice=1;
  SET(no).EndSlice=1;
end

im=SET(no).IM;

%Add papillary visualization, connected to LV hide / show
%Weirdest problem ever apparently entering a zero indexes into a ct image
%is extremely timeconsuming. Therefore extra if case
stateandicon=segment('iconson','hidelv');
%if not(stateandicon{1}) && any(SET(no).PapillaryIM==0) % not(isempty(SET(no).PapillaryIM)) && isequal(get(DATA.Handles.hidepapicon,'state'),'off')
if not(stateandicon{1}) && (~isempty(SET(no).PapillaryIM)) % not(isempty(SET(no).PapillaryIM)) && isequal(get(DATA.Handles.hidepapicon,'state'),'off')
  im(SET(no).PapillaryIM)=DATA.GUISettings.PapilarColor;
end

% %Make sure viewpixelyicon is in the right state
% if DATA.Pref.ViewInterpolated
% %   set(DATA.Handles.viewpixelyicon,'state','off');
%   set(DATA.Handles.viewpixelsmenu,'Checked','off');
% else
% %   set(DATA.Handles.viewpixelyicon,'state','on');
%   set(DATA.Handles.viewpixelsmenu,'Checked','on');
% end

scale = 2;
switch DATA.ViewPanelsType{panel}
  case {'one','mmodespatial','ortho'}
    if DATA.Pref.ViewInterpolated
      imxsz = SET(no).XSize*scale;
      imysz = SET(no).YSize*scale;
      imip = cast(zeros(imxsz,imysz,SET(no).TSize),'like',SET(no).IM);
      for t = 1:SET(no).TSize
        imip(:,:,t) = imresize(im(:,:,t,SET(no).CurrentSlice),[imxsz imysz],'bilinear');
      end
      DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',imip,no);
      %SET(no).NormalZoomState = [0.5;imysz+0.5;0.5;imxsz+0.5];
      %SET(no).ResolutionX = SET(no).ResolutionX*imxsz/SET(no).XSize;
      %SET(no).ResolutionY = SET(no).ResolutionY*imxsz/SET(no).YSize;
    else
      %ugly hack to handle orthoview mode in the right cases
      if strcmp(DATA.ViewPanelsType{panel},'ortho')&& ...
          ~isempty(DATA.ViewIM) && ...
          (size(DATA.ViewIM{1},4) ~= SET(no).ZSize || size(DATA.ViewIM{1},5) < size(SET(no).Colormap,2))
        DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im,no);
      elseif ~strcmp(DATA.ViewPanelsType{panel},'ortho') || panel ~= 1
        DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im(:,:,:,SET(no).CurrentSlice),no);
      end
    end
  case 'orthomip'
    if DATA.Pref.ViewInterpolated
      imxsz = SET(no).XSize*scale;
      imysz = SET(no).YSize*scale;
      imip = cast(zeros(imxsz,imysz,SET(no).TSize),'like',SET(no).IM);
      for t = 1:SET(no).TSize
        imip(:,:,t) = imresize(max(im(:,:,t,:),[],4),[imxsz imysz],'bilinear');
      end
      DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',imip,no);
      %SET(no).NormalZoomState = [0.5;imysz+0.5;0.5;imxsz+0.5];
      %SET(no).ResolutionX = SET(no).ResolutionX*imxsz/SET(no).XSize;
      %SET(no).ResolutionY = SET(no).ResolutionY*imxsz/SET(no).YSize;
    else
      %ugly hack to handle orthoview mode in the right cases
      if strcmp(DATA.ViewPanelsType{panel},'ortho')&& ...
          ~isempty(DATA.ViewIM) && ...
          (size(DATA.ViewIM{1},4) ~= SET(no).ZSize || size(DATA.ViewIM{1},5) < size(SET(no).Colormap,2))
        DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im,no);
      elseif strcmp(DATA.ViewPanelsType{panel},'orthomip')&& ...
          ~isempty(DATA.ViewIM) && ...
          (size(DATA.ViewIM{1},4) ~= SET(no).ZSize || size(DATA.ViewIM{1},5) < size(SET(no).Colormap,2))
        DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',max(im,[],4),no);
      elseif ~strncmp(DATA.ViewPanelsType{panel},'ortho',5) || panel ~= 1
        DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im(:,:,:,SET(no).CurrentSlice),no);
      end
    end
  case 'hla'
    if strcmp(DATA.CurrentTool,'orthoview')
      DATA.ViewIM{panel} = permute(DATA.ViewIM{1}(SET(no).HLA.slice,:,:,:,:),[4 2 3 1 5]);
    else
      im = permute(SET(no).IM(SET(no).HLA.slice,:,:,:),[4 2 3 1]);
      DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im,no);
    end
  case 'hlamip'
    im = permute(max(SET(no).IM,[],1),[4 2 3 1]);
    DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im,no);
  case 'vla'
    if strcmp(DATA.CurrentTool,'orthoview')
      DATA.ViewIM{panel} = permute(DATA.ViewIM{1}(:,SET(no).VLA.slice,:,:,:),[4 1 3 2 5]);
    else
      im = permute(SET(no).IM(:,SET(no).VLA.slice,:,:),[4 1 3 2]);
      DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im,no);
    end
  case 'vlamip'
    im = permute(max(SET(no).IM,[],2),[4 1 3 2]);
    DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im,no);
  case 'gla'
    glaangle = SET(no).GLA.angle;
    
    %Define coordinates of midpoint and find image edge
    x0 = SET(no).HLA.slice;
    y0 = SET(no).VLA.slice;
    tlim = zeros(1,4);
    %Here should SET(no).ResolutionX and ResolutionY be introduced somehow
    tlim(1) = SET(no).ResolutionX*(1-SET(no).HLA.slice)/sin(glaangle);
    tlim(2) = SET(no).ResolutionY*(1-SET(no).VLA.slice)/cos(glaangle);
    %tmin = max(tx,ty);
    tlim(3) = SET(no).ResolutionX*(SET(no).XSize-SET(no).HLA.slice)/sin(glaangle);
    tlim(4) = SET(no).ResolutionY*(SET(no).YSize-SET(no).VLA.slice)/cos(glaangle);
    %tmax = min(tx,ty);
    tlim = sort(tlim);
    tmin = tlim(2);
    tmax = tlim(3);
    
    t0 = SET(no).CurrentTimeFrame;
    z0 = SET(no).CurrentSlice;
    
    %Define extension
    sz = round(SET(no).YSize*cos(glaangle)^2+SET(no).XSize*sin(glaangle)^2);
    res = SET(NO).ResolutionY*cos(glaangle)+SET(NO).ResolutionX*abs(sin(glaangle));
    
    if cos(glaangle) >= 0
      x = x0 + sin(glaangle)*(tmin:res:tmax)/SET(no).ResolutionX; %linspace(tmin,tmax,sz);
      y = y0 + cos(glaangle)*(tmin:res:tmax)/SET(no).ResolutionY; %linspace(tmin,tmax,sz);
    else
      x = x0 + sin(glaangle)*(tmax:-res:tmin)/SET(no).ResolutionX; %linspace(tmin,tmax,sz);
      y = y0 + cos(glaangle)*(tmax:-res:tmin)/SET(no).ResolutionY;
    end
    SET(no).GLA.x0 = x(1);
    SET(no).GLA.y0 = y(1);
    
    t = 1:SET(no).TSize;
    z = 1:SET(no).ZSize;
    
    %Define image axes in xy plane and along z and t axis
    xyline = [x;y;t0*ones(size(x));z0*ones(size(x))];
    zline = [x0*ones(size(z));y0*ones(size(z));t0*ones(size(z));z];
    tline = [x0*ones(size(t));y0*ones(size(t));t;z0*ones(size(t))];
    
    X = meshgrid(xyline(1,:),zline(1,:),tline(1,:));
    Y = meshgrid(xyline(2,:),zline(2,:),tline(2,:));
    [~,~,T] = meshgrid(xyline(3,:),zline(3,:),tline(3,:));
    [~,Z] = meshgrid(xyline(4,:),zline(4,:),tline(4,:));
    if strcmp(DATA.CurrentTool,'orthoview')
      if size(DATA.ViewIM{1},3) == 1
        v = interpn(squeeze(DATA.ViewIM{1}),X,Y,Z,'nearest');
      else
        v = interpn(DATA.ViewIM{1},X,Y,T,Z,'nearest');
      end
    else
      if SET(no).TSize == 1
        im = interpn(squeeze(SET(no).IM),X,Y,Z,'nearest');
      else
        im = interpn(SET(no).IM,X,Y,T,Z,'nearest');
      end
      v = calcfunctions('remapuint8viewim',im,no);
    end
    DATA.ViewIM{panel} = v;
  case 'glamip'
    glaangle = SET(no).GLA.angle;
    
    %Define coordinates of midpoint and find image edge
    x0 = SET(no).HLA.slice;
    y0 = SET(no).VLA.slice;
    tlim = zeros(1,4);
    tlim(1) = SET(no).ResolutionX*(1-SET(no).HLA.slice)/sin(glaangle);
    tlim(2) = SET(no).ResolutionY*(1-SET(no).VLA.slice)/cos(glaangle);
    tlim(3) = SET(no).ResolutionX*(SET(no).XSize-SET(no).HLA.slice)/sin(glaangle);
    tlim(4) = SET(no).ResolutionY*(SET(no).YSize-SET(no).VLA.slice)/cos(glaangle);
    tlim = sort(tlim);
    tmin = tlim(2);
    tmax = tlim(3);
    
    t0 = SET(no).CurrentTimeFrame;
    z0 = SET(no).CurrentSlice;
    
    %Define extension
    sz = round(SET(no).YSize*cos(glaangle)^2+SET(no).XSize*sin(glaangle)^2);
    res = SET(NO).ResolutionY*cos(glaangle)+SET(NO).ResolutionX*abs(sin(glaangle));
    
    if cos(glaangle) >= 0
      x = x0 + sin(glaangle)*(tmin:res:tmax)/SET(no).ResolutionX; %linspace(tmin,tmax,sz);
      y = y0 + cos(glaangle)*(tmin:res:tmax)/SET(no).ResolutionY; %linspace(tmin,tmax,sz);
    else
      x = x0 + sin(glaangle)*(tmax:-res:tmin)/SET(no).ResolutionX; %linspace(tmin,tmax,sz);
      y = y0 + cos(glaangle)*(tmax:-res:tmin)/SET(no).ResolutionY;
    end
    SET(no).GLA.x0 = x(1);
    SET(no).GLA.y0 = y(1);
    
    t = 1:SET(no).TSize;
    z = 1:SET(no).ZSize;
    
    %Define image axes in xy plane and along z and t axis
    xyline = [x;y;t0*ones(size(x));z0*ones(size(x))];
    zline = [x0*ones(size(z));y0*ones(size(z));t0*ones(size(z));z];
    tline = [x0*ones(size(t));y0*ones(size(t));t;z0*ones(size(t))];
    
    X = meshgrid(xyline(1,:),zline(1,:),tline(1,:));
    Y = meshgrid(xyline(2,:),zline(2,:),tline(2,:));
    [~,~,T] = meshgrid(xyline(3,:),zline(3,:),tline(3,:));
    [~,Z] = meshgrid(xyline(4,:),zline(4,:),tline(4,:));
    if strcmp(DATA.CurrentTool,'orthoview')
      if size(DATA.ViewIM{1},3) == 1
        v = interpn(squeeze(DATA.ViewIM{1}),X,Y,Z,'nearest');
      else
        v = interpn(DATA.ViewIM{1},X,Y,T,Z,'nearest');
      end
    else
      if SET(no).TSize == 1
        im = interpn(squeeze(SET(no).IM),X,Y,Z,'nearest');
      else
        im = interpn(SET(no).IM,X,Y,T,Z,'nearest');
      end
      v = calcfunctions('remapuint8viewim',im,no);
    end
    DATA.ViewIM{panel} = v;
  case {'realhla','realvla'}
    type = upper(DATA.ViewPanelsType{panel});
    im = SET(no).(type).IM;
    if DATA.Pref.ViewInterpolated
      imxsz = size(im,1)*scale;
      imysz = size(im,2)*scale;
      im = zeros(imxsz,imysz,SET(no).TSize);
      for t = 1:SET(no).TSize
        im(:,:,t) = imresize(SET(no).(type).IM(:,:,t),[imxsz imysz],'bilinear');
      end
      %SET(no).NormalZoomState = [0.5;imysz+0.5;0.5;imxsz+0.5];
      %SET(no).ResolutionX = SET(no).ResolutionX*imxsz/SET(no).XSize;
      %SET(no).ResolutionY = SET(no).ResolutionY*imxsz/SET(no).YSize;
    end
    DATA.ViewIM{panel} = calcfunctions('remapuint8viewim',im(:,:,:),no);
  case 'mmodetemporal'
    DATA.ViewIM{panel} = [];
  case {'montage','montagerow'}
    %Update number of rows and columns
    if isequal(DATA.ViewPanelsType{panel},'montage')
      [rows,cols] = calcfunctions('calcrowscols',no); %.cols,.rows
      DATA.ViewPanelsMatrix{panel} = [rows cols];
    else
      if SET(no).ZSize>8
        DATA.ViewPanelsMatrix{panel}(1) = 2;
        DATA.ViewPanelsMatrix{panel}(2) = ceil(SET(no).ZSize/2);
      else
        DATA.ViewPanelsMatrix{panel}(1) = 1;
        DATA.ViewPanelsMatrix{panel}(2) = SET(no).ZSize;
      end;
    end;
    DATA.ViewIM{panel} = calcfunctions('calcmontageviewim', ...
      no,DATA.ViewPanelsMatrix{panel});
  case 'montagefit'
    %Update number of rows and columns
    if DATA.ViewMatrix(1) < DATA.ViewMatrix(2)
      DATA.ViewPanelsMatrix{panel} = [SET(no).ZSize 1];
    else
      DATA.ViewPanelsMatrix{panel} = [1 SET(no).ZSize];
    end;
    
    %Create space
    if isempty(SET(no).Colormap)
      DATA.ViewIM{panel} = repmat(uint8(0),[SET(no).XSize*DATA.ViewPanelsMatrix{panel}(1) SET(no).YSize*DATA.ViewPanelsMatrix{panel}(2) SET(no).TSize]);      
    else
      DATA.ViewIM{panel} = repmat(uint8(0),[SET(no).XSize*DATA.ViewPanelsMatrix{panel}(1) SET(no).YSize*DATA.ViewPanelsMatrix{panel}(2) SET(no).TSize 3]);
    end

    for tloop=1:SET(no).TSize
      for zloop=1:SET(no).ZSize
        c = 1+mod(zloop-1,DATA.ViewPanelsMatrix{panel}(2));
        r = ceil(zloop/DATA.ViewPanelsMatrix{panel}(2));
        DATA.ViewIM{panel}(...
          (1+(r-1)*SET(no).XSize):(r*SET(no).XSize),...
          (1+(c-1)*SET(no).YSize):(c*SET(no).YSize),tloop,:) = calcfunctions('remapuint8viewim',...
          im(:,:,tloop,zloop),no);
      end;
    end;
  case 'montagesegmented'
    %Montage view of segmentedslices only
    segslices = getmontagesegmentedslices(no);
    [rows,cols] = calcfunctions('calcrowscols',no,numel(segslices)); %.cols,.rows
    DATA.ViewPanelsMatrix{panel} = [rows cols];
    segmentedonly = true;
    DATA.ViewIM{panel} = calcfunctions('calcmontageviewim', ...
      no,DATA.ViewPanelsMatrix{panel},segmentedonly);
  case 'sax3'
    %Update number of rows and columns
    if DATA.ViewMatrix(1) < DATA.ViewMatrix(2)
      DATA.ViewPanelsMatrix{panel} = [3 1];
    else
      DATA.ViewPanelsMatrix{panel} = [1 3];
    end;
       
    %Create space
    if isempty(SET(no).Colormap)
      DATA.ViewIM{panel} = repmat(uint8(0),[SET(no).XSize*DATA.ViewPanelsMatrix{panel}(1) SET(no).YSize*DATA.ViewPanelsMatrix{panel}(2) SET(no).TSize]);      
    else
      DATA.ViewIM{panel} = repmat(uint8(0),[SET(no).XSize*DATA.ViewPanelsMatrix{panel}(1) SET(no).YSize*DATA.ViewPanelsMatrix{panel}(2) SET(no).TSize 3]);
    end

    for tloop=1:SET(no).TSize
      for zloop=1:3
        c = 1+mod(zloop-1,DATA.ViewPanelsMatrix{panel}(2));
        r = ceil(zloop/DATA.ViewPanelsMatrix{panel}(2));
        z = SET(no).SAX3.slices(zloop,tloop);
        DATA.ViewIM{panel}(...
          (1+(r-1)*SET(no).XSize):(r*SET(no).XSize),...
          (1+(c-1)*SET(no).YSize):(c*SET(no).YSize),tloop,:) = calcfunctions('remapuint8viewim',...
          im(:,:,tloop,z),no);
      end;
    end;
end;

if isempty(DATA.Overlay)
  DATA.Overlay = struct('alphadata',[],'cdata',[]);
end
DATA.Overlay(panel) = struct('alphadata',[],'cdata',[]);
myworkoff;

%----------------------
function montage_Motion %#ok<DEFNU>
%----------------------
%Motion function when selection slices in montage view.

global DATA SET NO

if DATA.Interactionlock
  return;
end;

[~,~,slice]=getclickedcoords;
% p = get(gca,'CurrentPoint');
% x=p(1);
% y=p(1,2);
% 
% col = 1+floor(x/SET(NO).YSize);
% row = 1+floor(y/SET(NO).XSize);
% slice = col+(row-1)*DATA.ViewPanelsMatrix{panel}(2);
if slice>=SET(NO).ZSize
  slice = SET(NO).ZSize;
end;
if slice<1
  slice=1;
end;
temp = SET(NO).StartSlice;
SET(NO).StartSlice = min(SET(NO).ZSize,max(SET(NO).StartSlice,1));
SET(NO).EndSlice = min(SET(NO).ZSize,max(SET(NO).EndSlice,1));
SET(NO).StartSlice = min([SET(NO).StartSlice slice SET(NO).EndSlice]);
SET(NO).EndSlice = max([temp slice SET(NO).EndSlice]);

%Update flows as well
if length(SET(NO).Linked) > 1
  nos = SET(NO).Linked;
  for no=unique(nos)
    %copy to linked flows
    SET(no).StartSlice=SET(NO).StartSlice;
    SET(no).CurrentSlice=SET(NO).CurrentSlice;
    SET(no).EndSlice=SET(NO).EndSlice;
  end
end
updateselectedslices;

drawfunctions('drawsliceno');

%--------------------------------
function montage_Buttonup(panel)  %#ok<INUSD,DEFNU>
%--------------------------------
%Button up fucntion for selecting slices in montage view.
global DATA NO

if DATA.Interactionlock
  return;
end;

% updateoneim(NO);

set(DATA.imagefig,'WindowButtonMotionFcn','');
set(DATA.imagefig,'WindowButtonUpFcn','');

drawfunctions('drawsliceno',NO);

%--------------------------------
function montage_Buttondown(panel) %#ok<DEFNU>
%--------------------------------
%Button down function for selecting slices in montage view.
global DATA SET NO

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

%Extract coordinates clicked
[~,~,slice] = getclickedcoords;

if (slice>SET(NO).ZSize)||(slice<1)
  return;
end;

%Check what type of click
switch get(DATA.imagefig,'SelectionType')
  case 'extend'
    %Shift-click to select range.
    %Same selection code as montage_Motion
    temp = SET(NO).StartSlice;
    SET(NO).StartSlice = min(SET(NO).ZSize,max(SET(NO).StartSlice,1));
    SET(NO).EndSlice = min(SET(NO).ZSize,max(SET(NO).EndSlice,1));
    SET(NO).StartSlice = min([SET(NO).StartSlice slice SET(NO).EndSlice]);
    SET(NO).EndSlice = max([temp slice SET(NO).EndSlice]);    

    %Update flows as well
    if length(SET(NO).Linked) > 1
      nos = SET(NO).Linked;
      for no=unique(nos)
        %copy to linked flows
        SET(no).StartSlice=SET(NO).StartSlice;
        SET(no).CurrentSlice=SET(NO).CurrentSlice;
        SET(no).EndSlice=SET(NO).EndSlice;
      end
    end

    updateselectedslices;
    drawfunctions('drawsliceno',NO);
%   pan_Buttondown; %shift-click formerly panned. No longer, though still
%   for one-view.
  case 'normal'
    set(DATA.imagefig,'WindowButtonUpFcn',sprintf('%s(''montage_Buttonup'')',mfilename));
    set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''montage_Motion'');',mfilename));
    switchtoslice(slice);
  case 'alt'
    if ~ismember(slice,SET(NO).StartSlice:SET(NO).EndSlice)
      switchtoslice(slice);
    end
    DATA.contextmenu;
  case 'open'
    no = DATA.ViewPanels(DATA.CurrentPanel);
    if numel(DATA.ViewPanels) > 1
      drawfunctions('drawall',1);
    elseif SET(no).ZSize>1
      viewimage_Callback('one');
    end;    
end;

%----------------------
function pan_Buttondown
%----------------------
%Button down function for panning of current image panel
global DATA 

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''pan_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''pan_Buttonup'')',mfilename));

%------------------------
function pan_Motion(init) %#ok<INUSD>
%------------------------
%Motion function of pan
global DATA
persistent startpos

%If called with argument reset and exit.
if nargin>0
  startpos = [];
  return;
end;

panel = DATA.CurrentPanel;
ax = DATA.Handles.imageaxes(panel);

%Get clicked point
[x,y] = mygetcurrentpoint(ax);

if isempty(startpos)
  startpos = [x y];
end;


xlim = get(ax,'xlim');
ylim = get(ax,'ylim');
if ismember(DATA.ViewPanelsType{panel},{'montagefit','sax3'})
  if all(diff(DATA.ViewPanelsMatrix{panel}) > 0)
    set(ax,...
    'xlim',xlim+startpos(1)-x);
  else
    set(ax,...
    'ylim',ylim+startpos(2)-y);
  end
else
  set(ax,...
    'xlim',xlim+startpos(1)-x,...
    'ylim',ylim+startpos(2)-y);
end

drawfunctions('viewupdatetextposition');
drawfunctions('viewupdateannotext');

%--------------------
function pan_Buttonup %#ok<DEFNU>
%--------------------
%Button up function for panning.
global DATA SET NO

%Reset startpos
pan_Motion(true);

%Restore so no motion is called
set(DATA.imagefig,'WindowButtonMotionFcn','');

%Restore main buttonup function
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''buttonup_Callback'')','segment'));

%Set zoom state for all linked image stacks
nos = SET(NO).Linked;
switch DATA.ViewPanelsType{DATA.CurrentPanel}
  case {'one','ortho'}
    [SET(nos).NormalZoomState] = deal([...
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'xlim')';
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'ylim')']);
  case 'montage'
    [SET(nos).MontageZoomState] = deal([...
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'xlim')';
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'ylim')']);    
  case 'montagerow'
    [SET(nos).MontageRowZoomState] = deal([...
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'xlim')';
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'ylim')']);      
  case {'montagefit','sax3','montagesegmented'}
    [SET(nos).MontageFitZoomState] = deal([...
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'xlim')';
      get(DATA.Handles.imageaxes(DATA.CurrentPanel),'ylim')']);      
end;

drawfunctions('viewupdatetextposition');
drawfunctions('viewupdateannotext');

%--------------------------
function center_Motion %#ok<DEFNU>
%--------------------------
%Motion function of the center point.

global DATA SET NO

[x,y,slice] = getclickedcoords;
panel = DATA.CurrentPanel;
% first line is for montage
% second line is for one-view.
if (slice==SET(NO).CurrentSlice) && ...
   (x>0.5) && (y>0.5) && (x<SET(NO).YSize+0.5) && (y<SET(NO).XSize+0.5)
  SET(NO).CenterX = y;
  SET(NO).CenterY = x;

  switch DATA.ViewPanelsType{panel}
    case {'one','mmodespatial','ortho'}
      set(DATA.Handles.center{panel},...
        'xdata',SET(NO).CenterY,...
        'ydata',SET(NO).CenterX,...
        'color',DATA.centercrossdef);
    case {'montage', 'montagerow','montagefit','sax3'}
      [tempx,tempy] = ndgrid(0:(DATA.ViewPanelsMatrix{panel}(1)-1),0:(DATA.ViewPanelsMatrix{panel}(2)-1));
      set(DATA.Handles.center{panel},...
        'xdata',SET(NO).CenterY+tempy(:)*SET(NO).YSize,...
        'ydata',SET(NO).CenterX+tempx(:)*SET(NO).XSize,...
        'color',DATA.centercrossdef);
  end;
end

%--------------------------------
function center_Buttondown(panel) %#ok<DEFNU>
%--------------------------------
%Called when center '+' is pressed down, sets motion and buttonup fcns.
global DATA SET NO

switchtopanel(panel);

if DATA.Interactionlock
  return;
end;

[~,~,slice] = getclickedcoords;
% If slice has changed, make sure montage/one are in sync
if (slice>SET(NO).ZSize)
  return;
end
switchtoslice(slice);

set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''center_Motion'');',mfilename));
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''center_Buttonup'')',mfilename));

%-----------------------
function center_Buttonup %#ok<DEFNU>
%-----------------------
%This function is called when buttonup occurs after draging center point

global DATA SET NO

%Restore so no motion is called
set(DATA.imagefig,'WindowButtonMotionFcn','');

%Restore main buttonup function
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''buttonup_Callback'')','segment'));

if length(SET(NO).Linked) > 1
  if ~isempty(SET(NO).Parent)
    nos = SET(NO).Linked; 
    for loop=1:length(nos)
      SET(nos(loop)).CenterX = SET(NO).CenterX;
      SET(nos(loop)).CenterY = SET(NO).CenterY;
    end;
  end;
end;
drawfunctions('drawimageno');

%------------------------------
function endo_Buttondown(panel) %#ok<DEFNU>
%------------------------------
%Button down function for manual draw of endocardium.
global DATA 

switchtopanel(panel);

switch get(DATA.imagefig,'SelectionType')
  case 'normal'
    if isequal(get(DATA.Handles.endopenicon,'state'),'on')
      manualdraw_Buttondown('endo',panel,0); %modify old
    else
      manualdraw_Buttondown('endo',panel,1); %add new
    end;
  case 'alt'
    disp('Later add popup menu here.');
end;

%-----------------------------
function epi_Buttondown(panel) %#ok<DEFNU>
%-----------------------------
%Button down function for manual draw of endocardium.
global DATA 

switchtopanel(panel);

switch get(DATA.imagefig,'SelectionType')
  case 'normal'
    if isequal(get(DATA.Handles.epipenicon,'state'),'on')
      manualdraw_Buttondown('epi',panel,0); %modify old
    else
      manualdraw_Buttondown('epi',panel,1); %add new
    end;
  case 'alt'
    disp('Later add popup menu here.');
end;

%-----------------------------
function do = doatrialscar(no)
%-----------------------------
%Helper function to check if user input is for atrial scar or lv scar
%(default)

global SET

do = false;

if isempty(SET(no).RVEndoX)
  return;
end;

if ~isnan(SET(no).RVEndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice))
  if isempty(SET(no).EpiX) || isnan(SET(no).EpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice))
    %RV exist and not LV Epi, then go for atrial scar
    do = true;
  end;
end;
    
%-----------------------------------------
function manualdraw_Buttonup(type,new,obj)  %#ok<DEFNU>
%-----------------------------------------
%Button up function for manual drawing.
global DATA SET NO

if nargin<1
  type = 'endo';
end;

if nargin<2
  new = true;
end;

if nargin==3
  try
    set(obj,'LineWidth',DATA.Pref.LineWidth);
  catch %#ok<CTCH>
  end;
end;

oldpref=DATA.ThisFrameOnly;

%--- Use variable no instead of NO. If flow data
%set then point to magnitude data set.

no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

if ~isempty(SET(no).Flow) && strcmp(type,'roi')
  DATA.ThisFrameOnly = false;
end

%check if resampling is needed
numpointscheck=[size(SET(no).EndoX,1),size(SET(no).RVEndoX,1),size(SET(no).EpiX,1),size(SET(no).RVEpiX,1)];
numpointscheck(isnan(numpointscheck))=[];
numpointscheck(numpointscheck==0)=[];

if any(DATA.NumPoints~=numpointscheck)
  %Do all curve interpolations. this is safety measure only done
  %occasionally when something has happened.
  for loop=1:length(SET)
    segpref('numpointsedithelper', loop);
    segment('updatemodeldisplay',loop);
  end;
end

tools('enableundo',no);

slice=SET(no).CurrentSlice;

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

%Restore
set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
if not(isempty(DATA.Handles.phasecursor))
  set(DATA.Handles.phasecursor,'visible','off');
end;
DATA.Handles.phasecursor = [];

DATA.buttonup_Callback;
set(DATA.imagefig,'WindowButtonUpFcn',sprintf('%s(''buttonup_Callback'')','segment'));%Restore

%Extract data
x = DATA.CursorX(1:DATA.CursorN);
y = DATA.CursorY(1:DATA.CursorN);

%Restore
DATA.CursorN = 0;

%Calculate total length
if length(x)>1  
  len = sqrt(...
    conv2(x,[1 -1],'valid').^2+...
    conv2(y,[1 -1],'valid').^2);
  len = [0;len(:)]; %Add zero first
  len = cumsum(len);
  tempind = find(conv2(len,[1;-1],'valid')~=0); %Remove doublets
  len = [len(1);len(tempind+1)];
  x = [x(1) x(tempind+1)];
  y = [y(1) y(tempind+1)];  
  totallength = len(end);
else
  drawfunctions('drawsliceno');
  % If No Scar data, clear the struct. (Handles empty clicks)
  if ~isempty(SET(NO).Scar)
    if ~any(SET(NO).Scar.Manual(:)) && ~any(SET(NO).Scar.Auto(:))
      viability('viabilityclear_Callback');
    end
  end
  return;
end;

%Check if counterclock wise or clockwise
mx = mean(x);
my = mean(y);

if max(abs(x-mx))>2*SET(no).XSize
  disp('Got strange x-coordinates before interp1.');
  return;
end;

if max(abs(y-my))>2*SET(no).YSize
  disp('Got strange y-coordinates before interp1.');
  return;
end;

if sum(unwrap(conv2(angle(complex(x-mx,y-my)),[1 -1],'valid')))<0
  x = fliplr(x);
  y = fliplr(y);
  %disp('counterclockwise');
else
  %disp('clockwise');
end;

%Resample
xr = interp1(len,x,linspace(0,totallength,DATA.NumPoints),'linear');
yr = interp1(len,y,linspace(0,totallength,DATA.NumPoints),'linear');
xr = xr(:);
yr = yr(:);

if max(abs(xr-mx))>2*SET(no).XSize
  myfailed('Got strange x-coordinates after interp1.',DATA.GUI.Segment);
  return;
end;

if max(abs(yr-my))>2*SET(no).YSize
  myfailed('Got strange y-coordinates after interp1.',DATA.GUI.Segment);
  return;
end;

%This variable used to be dependent on ThisFrameOnly. Check is now done
%later in this fcn ---Correction this is in effect now. /Klas
%if DATA.ThisFrameOnly
  timeframes = SET(no).CurrentTimeFrame;
%else
%  timeframes = 1:SET(no).TSize;%SET(no).CurrentTimeFrame;
%end

%numtimeframes=length(timeframes);

if new
  %--- Add end, resample to equidistant points
  x= [x(:);x(1)];
  y= [y(:);y(1)];
  len = sqrt(...
    conv2(x(:),[1;-1],'valid').^2+...
    conv2(y(:),[1;-1],'valid').^2);
  len = cumsum(len);
  len = [0;len(:)]; %Add zero first
  tempind = find(abs(conv2(len,[1;-1],'valid'))>1e-3);
  len = [len(1);len(tempind+1)]; %remove doublets
  x = [x(1);x(tempind+1)]; 
  y = [y(1);y(tempind+1)];    
  totallength = len(end);
  xr = interp1(len,x,linspace(0,totallength,DATA.NumPoints),'linear');
  yr = interp1(len,y,linspace(0,totallength,DATA.NumPoints),'linear');
  xr = xr(:);
  yr = yr(:);
end;

switch type
  case {'scar','mo'}
    %Check if it is atrial scar drawing, or normal LV scar
    
    %If atria scar then call special functions and exit otherwise continue
    if doatrialscar(no)
      atrialscar('manualdraw_Buttonup',no,type,xr,yr);
      return;
    end;
    
    if isempty(SET(no).Scar)
      viability('viabilityreset_Callback');
      SET(no).Scar.Mode = 'manual';
    end;

    tempmask = createmask([SET(no).XSize SET(no).YSize],xr,yr) & SET(no).Scar.MyocardMask(:,:,slice);
    temp = SET(no).Scar.Manual(:,:,slice);
    
    if isequal(type,'scar')
      temp(tempmask) = int8(1); %Mark manual scar as 1
    else
      temp(tempmask) = int8(2); %Mark manually no reflow as 2
    end;

    SET(no).Scar.Manual(:,:,slice) = temp;
    
    if SET(no).Scar.UpdateDirectly
      viability('viabilitycalc');
    end;
   
    drawfunctions('drawimageno');
    
    return;
  case 'rubberpen'

    %Check if atrial scar then call special code, otherwise continue
    if doatrialscar(no)
      atrialscar('manualdraw_Buttonup',no,type,xr,yr);
      return;
    end;
    
    if isempty(SET(no).Scar)
      myfailed('You need to enter viability view mode before drawing infarct regions.',DATA.GUI.Segment);
      return;
    end;
    
    %Create mask
    tempmask = createmask([SET(no).XSize SET(no).YSize],xr,yr)&SET(no).Scar.MyocardMask(:,:,slice);
    
    %Update scar
    temp = SET(no).Scar.Manual(:,:,slice);
    temp(tempmask) = int8(-1);
    SET(no).Scar.Manual(:,:,slice) = temp;
    
    %Update noreflow
    temp = SET(no).Scar.NoReflow(:,:,slice);
    temp(tempmask) = 0;
    SET(no).Scar.NoReflow(:,:,slice) = temp;
    
    if SET(no).Scar.UpdateDirectly
      viability('viabilitycalc');
    end;
    drawfunctions('drawimageno');

    return;
  case 'rubber'
    if isempty(SET(no).Scar)
      myfailed('You need to enter viability view mode before drawing infarct regions.',DATA.GUI.Segment);
      return;
    end;
    tempmask = createmask([SET(no).XSize SET(no).YSize],xr,yr)&SET(no).Scar.MyocardMask(:,:,slice);
    temp = SET(no).Scar.Manual(:,:,slice);
    temp(tempmask) = int8(0);
    SET(no).Scar.Manual(:,:,slice) = temp;
    if SET(no).Scar.UpdateDirectly
      viability('viabilitycalc');
    end;
    
    drawfunctions('drawimageno');    
    %if not(specialgui)
    %  viability('viabilityshowedits_Callback','on');
    %end;

    return;
  case 'drawmarpen'
    mar('createmyocardmask');
    SET(no).MaR.Mode = 'manual';
    tempmask = createmask([SET(no).XSize SET(no).YSize],xr,yr)&SET(no).MaR.MyocardMask(:,:,timeframes,slice);
    temp = SET(no).MaR.Manual(:,:,timeframes,slice);

    temp(tempmask) = int8(1); %Mark manual scar as 1
    SET(no).MaR.Manual(:,:,timeframes,slice) = temp;
    mar('update');

    drawfunctions('drawimageno');
    return;
  case 'drawmarrubberpen'
    mar('createmyocardmask');
    tempmask = createmask([SET(no).XSize SET(no).YSize],xr,yr)&SET(no).MaR.MyocardMask(:,:,timeframes,slice);
    temp = SET(no).MaR.Manual(:,:,timeframes,slice);
    temp(tempmask) = int8(-1);
    SET(no).MaR.Manual(:,:,timeframes,slice) = temp;
    mar('update');

    drawfunctions('drawimageno');
    return;
  case 'drawmarrubber'
    mar('createmyocardmask');
    tempmask = createmask([SET(no).XSize SET(no).YSize],xr,yr)&SET(no).MaR.MyocardMask(:,:,timeframes,slice);
    temp = SET(no).MaR.Manual(:,:,timeframes,slice);
    temp(tempmask) = int8(0);
    SET(no).MaR.Manual(:,:,timeframes,slice) = temp;
    mar('update');

    drawfunctions('drawimageno');
    return;
end;

%--- Type is either roi, endo or epi => more advanced...
pindensity = 3;

if new
  %--- New
  pindensity = 3;
  
  %Find startpoint

  [~,inda] = min(angle(complex(mx-xr,my-yr)));

  %Reshape
  xr = xr';
  yr = yr';
  %Store it
  switch type
    case 'endo'
      %New endo
      if isempty(SET(no).EndoX)
        SET(no).EndoX = nan([DATA.NumPoints SET(no).TSize SET(no).ZSize]);
        SET(no).EndoY = SET(no).EndoX;
      end;
      SET(no).EndoX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).EndoY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).EndoX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).EndoY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';
%       SET(no).EndoX(1:(end-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end)',1,numtimeframes);
%       SET(no).EndoY(1:(end-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end)',1,numtimeframes);
%       SET(no).EndoX((end+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda)',1,numtimeframes);
%       SET(no).EndoY((end+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda)',1,numtimeframes);
    case 'epi'
      %New epi
      if isempty(SET(no).EpiX)
        SET(no).EpiX = nan([DATA.NumPoints SET(no).TSize SET(no).ZSize]);
        SET(no).EpiY = SET(no).EpiX;
      end;      
      SET(no).EpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).EpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).EpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).EpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';
%       SET(no).EpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end)',1,numtimeframes);
%       SET(no).EpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end)',1,numtimeframes);
%       SET(no).EpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda)',1,numtimeframes);
%       SET(no).EpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda)',1,numtimeframes);
    case 'rvendo'
      %New rv-endo
      if isempty(SET(no).RVEndoX)
        SET(no).RVEndoX = nan([DATA.NumPoints SET(no).TSize SET(no).ZSize]);
        SET(no).RVEndoY = SET(no).RVEndoX;
      end;
            SET(no).RVEndoX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).RVEndoY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).RVEndoX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).RVEndoY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';
%       SET(no).RVEndoX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end)',1,numtimeframes);
%       SET(no).RVEndoY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end)',1,numtimeframes);
%       SET(no).RVEndoX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda)',1,numtimeframes);
%       SET(no).RVEndoY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda)',1,numtimeframes);
    case 'rvepi'
      %New rv-epi
      if isempty(SET(no).RVEpiX)
        SET(no).RVEpiX = nan([DATA.NumPoints SET(no).TSize SET(no).ZSize]);
        SET(no).RVEpiY = SET(no).RVEpiX;
      end;      
      SET(no).RVEpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).RVEpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).RVEpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).RVEpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';
%       SET(no).RVEpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end)',1,numtimeframes);
%       SET(no).RVEpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end)',1,numtimeframes);
%       SET(no).RVEpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda)',1,numtimeframes);
%       SET(no).RVEpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda)',1,numtimeframes);
    case 'roi'
      %New roi
      roiok = DATA.manualdraw_Buttonup_roi(no,xr,yr,slice);
      if ~roiok
        return
      end
    otherwise
      myfailed(dprintf('Unknown draw type:%s should not occur.',type),DATA.GUI.Segment);
      return;
  end;
end;

if not(new)
  %--- Not closed try to fit in

  %Extract
  switch type
    case 'endo'
      contoury = SET(no).EndoX(:,timeframes,SET(no).CurrentSlice)';
      contourx = SET(no).EndoY(:,timeframes,SET(no).CurrentSlice)';
    case 'epi'
      contoury = SET(no).EpiX(:,timeframes,SET(no).CurrentSlice)';
      contourx = SET(no).EpiY(:,timeframes,SET(no).CurrentSlice)';
    case 'rvendo'
      contoury = SET(no).RVEndoX(:,timeframes,SET(no).CurrentSlice)';
      contourx = SET(no).RVEndoY(:,timeframes,SET(no).CurrentSlice)';      
    case 'rvepi'
      contoury = SET(no).RVEpiX(:,timeframes,SET(no).CurrentSlice)';
      contourx = SET(no).RVEpiY(:,timeframes,SET(no).CurrentSlice)';      
    case 'roi'
      no = roi('roifindmag');
      %--- Try to find correct roi
      oldmin=1e9;
      for rloop=1:SET(no).RoiN
        if SET(no).CurrentSlice==SET(no).Roi(rloop).Z
          temp = min(sqrt(...
            (SET(no).Roi(rloop).X(:,SET(no).CurrentTimeFrame)-y(1)).^2+...
            (SET(no).Roi(rloop).Y(:,SET(no).CurrentTimeFrame)-x(1)).^2));
          if temp<oldmin
            SET(no).RoiCurrent=rloop;
            oldmin = temp;
          end;
        end;
      end;
      %Extract data
      contoury = SET(no).Roi(SET(no).RoiCurrent).X(:,timeframes)';
      contourx = SET(no).Roi(SET(no).RoiCurrent).Y(:,timeframes)';
  end;

  %Find closest point to start/end point
  [~,startind] = min((contourx-x(1)).^2+(contoury-y(1)).^2);
  [~,endind] = min((contourx-x(end)).^2+(contoury-y(end)).^2);

  if startind>endind
    temp = startind;
    startind = endind;
    endind = temp;
  end;

  %--- Find the two parts where the contour is cut.
  %part1
  part1x = contourx(startind:endind);
  part1y = contoury(startind:endind);

  %part2
  part2x = [contourx(endind:end) contourx(2:startind)];
  part2y = [contoury(endind:end) contoury(2:startind)];

  %Test join to fit
  if ...
      ( (part1x(end)-x(1))^2 + (part1y(end)-y(1))^2 ) < ...
      ( (part1x(end)-x(end))^2 + (part1y(end)-y(end))^2 )
    part1x = [part1x x];
    part1y = [part1y y];
  else
    part1x = [part1x fliplr(x)];
    part1y = [part1y fliplr(y)];
  end;
  area1 = stablepolyarea(part1x,part1y);

  if ...
      ( (part2x(end)-x(1))^2 + (part2y(end)-y(1))^2 ) < ...
      ( (part2x(end)-x(end))^2 + (part2y(end)-y(end))^2 )
    part2x = [part2x x];
    part2y = [part2y y];
  else
    part2x = [part2x fliplr(x)];
    part2y = [part2y fliplr(y)];
  end;
  area2 = stablepolyarea(part2x,part2y);

  %Take the largest area
  if area1>area2
    newcontourx = part1x;
    newcontoury = part1y;
  else
    newcontourx = part2x;
    newcontoury = part2y;
  end;

  %Resample to equidistant points
  newlen = sqrt(...
    conv2(newcontourx,[1 -1],'valid').^2+...
    conv2(newcontoury,[1 -1],'valid').^2);
  newlen = [0;newlen(:)]; %Add zero first
  newlen = cumsum(newlen);
  tempind = find(abs(conv2(newlen,[1;-1],'valid'))>1e-3);
  newlen = [newlen(1);newlen(tempind+1)];
  newcontourx = [newcontourx(1) newcontourx(tempind+1)];
  newcontoury = [newcontoury(1) newcontoury(tempind+1)]; 
  newtotallength = newlen(end);

  xr = interp1(newlen,newcontourx,linspace(0,newtotallength,DATA.NumPoints),'linear');
  yr = interp1(newlen,newcontoury,linspace(0,newtotallength,DATA.NumPoints),'linear');

  %Find startpoint
  %[tempres,inda] = min(abs(angle(complex(mx-xr,my-yr))));
  [~,inda] = min(angle(complex(mx-xr,my-yr)));  
  
  %Reshape
  xr = xr';
  yr = yr';

  %Store it
  switch type
    case 'endo'
      SET(no).EndoX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).EndoY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).EndoX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).EndoY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';
      
      %       SET(no).EndoX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end),1,numtimeframes);
%       SET(no).EndoY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end),1,numtimeframes);
%       SET(no).EndoX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda),1,numtimeframes);
%       SET(no).EndoY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda),1,numtimeframes);
      %Remove interpolation points if existing
      if ~isempty(SET(no).EndoInterpX)
        [SET(no).EndoInterpX{timeframes,SET(no).CurrentSlice}] = deal([]);
        [SET(no).EndoInterpY{timeframes,SET(no).CurrentSlice}] = deal([]);
      end
    case 'epi'
      
      SET(no).EpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).EpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).EpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).EpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';
%       SET(no).EpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end),1,numtimeframes);
%       SET(no).EpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end),1,numtimeframes);
%       SET(no).EpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda),1,numtimeframes);
%       SET(no).EpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda),1,numtimeframes);
      %Remove interpolation points if existing
      if ~isempty(SET(no).EpiInterpX)
        [SET(no).EpiInterpX{timeframes,SET(no).CurrentSlice}] = deal([]);
        [SET(no).EpiInterpY{timeframes,SET(no).CurrentSlice}] = deal([]);
      end
    case 'rvendo'
      SET(no).RVEndoX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).RVEndoY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).RVEndoX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).RVEndoY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';%       SET(no).RVEndoX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end),1,numtimeframes);
%       SET(no).RVEndoY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end),1,numtimeframes);
%       SET(no).RVEndoX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda),1,numtimeframes);
%       SET(no).RVEndoY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda),1,numtimeframes);      
%       %Remove interpolation points if existing
      if ~isempty(SET(no).RVEndoInterpX)
        [SET(no).RVEndoInterpX{timeframes,SET(no).CurrentSlice}] = deal([]);
        [SET(no).RVEndoInterpY{timeframes,SET(no).CurrentSlice}] = deal([]);
      end
     case 'rvepi'
      SET(no).RVEpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = yr(inda:end)';
      SET(no).RVEpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = xr(inda:end)';
      SET(no).RVEpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = yr(1:inda)';
      SET(no).RVEpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = xr(1:inda)';
       %       SET(no).RVEpiX(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(yr(inda:end),1,numtimeframes);
%       SET(no).RVEpiY(1:(DATA.NumPoints-inda+1),timeframes,SET(no).CurrentSlice) = repmat(xr(inda:end),1,numtimeframes);
%       SET(no).RVEpiX((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(yr(1:inda),1,numtimeframes);
%       SET(no).RVEpiY((DATA.NumPoints+1-inda):end,timeframes,SET(no).CurrentSlice) = repmat(xr(1:inda),1,numtimeframes);
      %Remove interpolation points if existing
      if ~isempty(SET(no).RVEpiInterpX)
        [SET(no).RVEpiInterpX{timeframes,SET(no).CurrentSlice}] = deal([]);
        [SET(no).RVEpiInterpY{timeframes,SET(no).CurrentSlice}] = deal([]);
      end
    case 'roi'
      SET(no).Roi(SET(no).RoiCurrent).X(1:(DATA.NumPoints-inda+1),timeframes) = yr(inda:end)';
      SET(no).Roi(SET(no).RoiCurrent).Y(1:(DATA.NumPoints-inda+1),timeframes) = xr(inda:end)';
      SET(no).Roi(SET(no).RoiCurrent).X((DATA.NumPoints+1-inda):end,timeframes) = yr(1:inda)';
      SET(no).Roi(SET(no).RoiCurrent).Y((DATA.NumPoints+1-inda):end,timeframes) = xr(1:inda)';
      roi('roiforceapply');
      
      %Calculate area and intensity of ROI
      thisframeonly = true;
      [~,SET(no).Roi(SET(no).RoiCurrent).Area(SET(no).CurrentTimeFrame)] = ...
        calcfunctions('calcroiarea',no,SET(no).RoiCurrent,thisframeonly);
      [m,sd]=calcfunctions('calcroiintensity',no,SET(no).RoiCurrent,false,thisframeonly);
      SET(no).Roi(SET(no).RoiCurrent).Mean(SET(no).CurrentTimeFrame) = m;
      SET(no).Roi(SET(no).RoiCurrent).StD(SET(no).CurrentTimeFrame) = sd;
  end;
end; %Done not closed contour

%--- Double check rotation direction
switch type
  case 'endo'
    xr = SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    yr = SET(no).EndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
  case 'epi'
    xr = SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    yr = SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
  case 'rvendo'
    xr = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    yr = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
  case 'rvepi'
    xr = SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    yr = SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);  
  case 'roi'
    xr = SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame);
    yr = SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame);
end;

mx = mean(xr);
my = mean(yr);

if sum(unwrap(conv2(angle(complex(xr-mx,yr-my)),[1;-1],'valid')))<0
  %disp('counterclockwise');
else
  disp('Wrong direction detected. Fixed but please report to author');
  switch type
    case 'endo'
      for tloop=timeframes
        SET(no).EndoX(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).EndoX(:,tloop,SET(no).CurrentSlice));
        SET(no).EndoY(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).EndoY(:,tloop,SET(no).CurrentSlice));
      end;
    case 'epi'
      for tloop=timeframes
        SET(no).EpiX(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).EpiX(:,tloop,SET(no).CurrentSlice));
        SET(no).EpiY(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).EpiY(:,tloop,SET(no).CurrentSlice));
      end;
    case 'rvendo'
      for tloop=timeframes
        SET(no).RVEndoX(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).RVEndoX(:,tloop,SET(no).CurrentSlice));
        SET(no).RVEndoY(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).RVEndoY(:,tloop,SET(no).CurrentSlice));
      end;      
    case 'rvepi'
      for tloop=timeframes
        SET(no).RVEpiX(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).RVEpiX(:,tloop,SET(no).CurrentSlice));
        SET(no).RVEpiY(:,tloop,SET(no).CurrentSlice) = flipud(SET(no).RVEpiY(:,tloop,SET(no).CurrentSlice));
      end;      
    case 'roi'
      for tloop=timeframes
        SET(no).Roi(SET(no).RoiCurrent).X(:,tloop) = flipud(SET(no).Roi(SET(no).RoiCurrent).X(:,tloop));
        SET(no).Roi(SET(no).RoiCurrent).Y(:,tloop) = flipud(SET(no).Roi(SET(no).RoiCurrent).Y(:,tloop));
      end;
      roi('roiforceapply');
  end;
end;

%Check if should copy to all timeframes
switch type
  case 'endo'
    if not(DATA.ThisFrameOnly)
      %Copy segmentation to all timeframes
      for tloop=1:SET(no).TSize
        SET(no).EndoX(:,tloop,SET(no).CurrentSlice) = SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        SET(no).EndoY(:,tloop,SET(no).CurrentSlice) = SET(no).EndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      end;
    end;
  case 'epi'
    if not(DATA.ThisFrameOnly)
      %Copy segmentation to all timeframes
      for tloop=1:SET(no).TSize
        SET(no).EpiX(:,tloop,SET(no).CurrentSlice) = SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        SET(no).EpiY(:,tloop,SET(no).CurrentSlice) = SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      end;
    end;
  case 'rvendo'
    if not(DATA.ThisFrameOnly)
      %Copy segmentation to all timeframes
      for tloop=1:SET(no).TSize
        SET(no).RVEndoX(:,tloop,SET(no).CurrentSlice) = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        SET(no).RVEndoY(:,tloop,SET(no).CurrentSlice) = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      end;
    end;    
  case 'rvepi'
    if not(DATA.ThisFrameOnly)
      %Copy segmentation to all timeframes
      for tloop=1:SET(no).TSize
        SET(no).RVEpiX(:,tloop,SET(no).CurrentSlice) = SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        SET(no).RVEpiY(:,tloop,SET(no).CurrentSlice) = SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      end;
    end;    
  case 'roi'
    if DATA.ThisFrameOnly && new
      %Remove segmentation from all timeframes other than current
      SET(no).Roi(SET(no).RoiCurrent).T = SET(no).CurrentTimeFrame;
      roix = SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame);
      roiy = SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame);
      SET(no).Roi(SET(no).RoiCurrent).X = nan(size(SET(no).Roi(SET(no).RoiCurrent).X));
      SET(no).Roi(SET(no).RoiCurrent).Y = nan(size(SET(no).Roi(SET(no).RoiCurrent).Y));
      SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame) = roix;
      SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame) = roiy;
    end;    
end;

%--- Add pins
if (isequal(type,'endo')||isequal(type,'epi'))&&DATA.Pref.AddPoints
  %Resample to equidistant points
  x=x(:);
  y=y(:);
  len = sqrt(...
    conv2(x,[1;-1],'valid').^2+...
    conv2(y,[1;-1],'valid').^2);
  len = [0;len(:)]; %Add zero first
  len = cumsum(len);
  tempind = find(conv2(len,[1;-1],'valid')~=0);
  len = [len(1);len(tempind+1)]; %Remove doublets
  x = [x(1);x(tempind+1)];
  y = [y(1);y(tempind+1)];
  totallength = len(end);
  xi = interp1(len,x,linspace(0,totallength,round(totallength/pindensity)),'linear');
  yi = interp1(len,y,linspace(0,totallength,round(totallength/pindensity)),'linear');
  xi = xi(:);
  yi = yi(:);

  %Put pins
  switch type
    case 'endo'
      if isempty(SET(no).EndoPinX)
        SET(no).EndoPinX = cell(SET(no).TSize,SET(no).ZSize);
        SET(no).EndoPinY = cell(SET(no).TSize,SET(no).ZSize);
      end;
      SET(no).EndoPinX{SET(no).CurrentTimeFrame,slice} = ...
        [SET(no).EndoPinX{SET(no).CurrentTimeFrame,slice};yi];
      SET(no).EndoPinY{SET(no).CurrentTimeFrame,slice} = ...
        [SET(no).EndoPinY{SET(no).CurrentTimeFrame,slice};xi];
    case 'epi'
      if isempty(SET(no).EpiPinX)
        SET(no).EpiPinX = cell(SET(no).TSize,SET(no).ZSize);
        SET(no).EpiPinY = cell(SET(no).TSize,SET(no).ZSize);
      end;      
      SET(no).EpiPinX{SET(no).CurrentTimeFrame,slice} = ...
        [SET(no).EpiPinX{SET(no).CurrentTimeFrame,slice};yi];
      SET(no).EpiPinY{SET(no).CurrentTimeFrame,slice} = ...
        [SET(no).EpiPinY{SET(no).CurrentTimeFrame,slice};xi];
  end;
end;

%Double check!
if isequal(type,'endo')||isequal(type,'epi')
  checkconsistency(SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
end;

%Finalize & update
updatemodeldisplay;
lvsegchanged = true; updatevolume(lvsegchanged);
if isequal(type,'roi') 
  docalc = false; %Intensity and area already updated
  %This is needed to initiate ROI handles. Slice update won't do
  for panel = find(ismember(DATA.ViewPanels,[no SET(no).Children]))
    if strcmp(DATA.ViewPanelsType{panel},'montage')
      DATA.plotrois(panel,no,docalc);
    else
      DATA.drawroiinpanel(panel,docalc);
    end
  end
  updatetool; %to set buttondown fcn for new contours
else
  drawfunctions('drawallslices');
end

if ~isempty(SET(no).StrainTagging)
  SET(no).StrainTagging.LVupdated=1;
end

if oldpref~=DATA.ThisFrameOnly
  DATA.ThisFrameOnly = oldpref;
end

%-------------------------
function manualdraw_Motion %#ok<DEFNU>
%-------------------------
%Motion function for manual drawings.
global DATA SET NO

[x,y,slice] = getclickedcoords;

xofs = 0;
yofs = 0;
if any(strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},{'montage','montagerow','montagefit','sax3','montagesegmented'}))
  if ~isequal(slice,SET(NO).CurrentSlice)
    if abs(slice-SET(NO).CurrentSlice)==1
      xofs = slice-SET(NO).CurrentSlice;
    else
      yofs = sign(slice-SET(NO).CurrentSlice);
    end;
  end;
end;

% first line is for montage
% second line is for one-view.
if (slice==SET(NO).CurrentSlice) && ...
   (x>0.5) && (y>0.5) && (x<SET(NO).YSize+0.5) && (y<SET(NO).XSize+0.5)

  DATA.CursorN = DATA.CursorN+1;
  DATA.CursorX(DATA.CursorN) = x+xofs*SET(NO).XSize;
  DATA.CursorY(DATA.CursorN) = y+yofs*SET(NO).YSize;
  
  set(cat(2,DATA.Handles.cursor(DATA.CurrentPanel), DATA.Handles.phasecursor),...
    'XData',DATA.CursorXOfs+DATA.CursorX(1:DATA.CursorN),...
    'YData',DATA.CursorYOfs+DATA.CursorY(1:DATA.CursorN));

end
%------------------------------------------------
  function buttondowntoggler(caller,panel)
%-----------------------------------------------
%when clicking in a image toggles to the correct handle buttondown
global SET NO DATA

[x,y,slice] = getclickedcoords;  

distlimit = DATA.Pref.ContourAdjustDistance;

dist=2*distlimit*ones(1,5);

if ~isempty(SET(NO).EndoX)
  dist(1) = min(sqrt(...
    (SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
    (SET(NO).EndoY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
end;

if ~isempty(SET(NO).EpiX)
  dist(2) = min(sqrt(...
    (SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
    (SET(NO).EpiY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
end;

if ~isempty(SET(NO).RVEndoX)
  dist(3) = min(sqrt(...
    (SET(NO).RVEndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
    (SET(NO).RVEndoY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
end;

if ~isempty(SET(NO).RVEpiX)
  dist(4) = min(sqrt(...
    (SET(NO).RVEpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
    (SET(NO).RVEpiY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
end

%Calc distance to contour
no = roi('roifindmag');
for rloop=1:SET(no).RoiN
  if ~isempty(find(SET(no).Roi(rloop).Z==SET(no).CurrentSlice,1)) && ~isempty(find(SET(no).Roi(rloop).T==SET(no).CurrentTimeFrame,1))
    temp = min(sqrt(...
      (SET(no).Roi(rloop).X(:,SET(no).CurrentTimeFrame)-y).^2+...
      (SET(no).Roi(rloop).Y(:,SET(no).CurrentTimeFrame)-x).^2));
    if temp<dist(end)
      dist(end) = temp;
    end;
  end;
end;

typecell={'endo','epi','rvendo','rvepi','roi'};
[~,ind]=min(dist);

if any(dist<distlimit)
  type = typecell{ind};
else
  type = 'image';
end

if strcmp('move',caller)
  move_Buttondown(type,panel);
else
  scale_Buttondown(type,panel);
end

%---------------------------------------------
function manualdraw_Buttondown(panel,type,new)
%---------------------------------------------
%Button down function for manual drawings.
global DATA SET NO

killbuttondown = switchtopanel(panel);

if killbuttondown
  return
end

if isequal(DATA.CurrentTool,'select')
  return;
end

if DATA.Interactionlock
  return;
end;

if nargin<2
  myfailed('Expected two or more input arguments.',DATA.GUI.Segment);
  return;
end;

if nargin<3
  new = true;
end;

%Check what type of click
switch get(DATA.imagefig,'SelectionType')
  case 'normal'

    %Setup
    [x,y,slice] = getclickedcoords;    
    % If slice has changed, make sure montage/one are in sync
    if (slice>SET(NO).ZSize)
      return;
    end
    switchtoslice(slice);
    
    %Moved this before setting of buttonup and button 
     if any(strcmp({'scar','mo','rubberpen'},type))
        %if isempty(SET(NO).Scar)
        hasLVseg = (not(isempty(SET(NO).EndoX)) && not(all(all(squeeze(isnan(SET(NO).EndoX(1,:,:))))))) && (not(isempty(SET(NO).EpiX)) && not(all(all(squeeze(isnan(SET(NO).EpiX(1,:,:)))))));
        hasRVseg = not(isempty(SET(NO).RVEndoX)) && not(all(all(squeeze(isnan(SET(NO).RVEndoX(1,:,:))))));
        if SET(NO).TSize>1 || hasRVseg+hasLVseg < 1
          myfailed('Does not seem to be viability image stack.',DATA.GUI.Segment);
          updatetool('select')
          stateandicon_scar=iconson('scarpen');
          stateandicon_mo=iconson('mopen');
          stateandicon_drp=iconson('rubberscar');
          stateandicon_scar{2}.undent;
          stateandicon_mo{2}.undent;
          stateandicon_drp{2}.undent;
          DATA.Handles.configiconholder.iconCell{1}.isindented=1;
          DATA.Handles.configiconholder.iconCell{1}.cdataDisplay=DATA.Handles.configiconholder.iconCell{1}.cdataIndent;
          DATA.Handles.configiconholder.render;
          return;
        end;
     end
     
    %Set up buttonup. If 'obj' is set during the execution below, then
    %WindowButtonUpFcn is remapped with obj passed on at the bottom of
    %this function. Obj is used to bold and then unbold an item.
    obj = [];
    set(DATA.imagefig,'WindowButtonUpFcn',...
      sprintf('%s(''manualdraw_Buttonup'',''%s'',%d )',mfilename,type,new));

    set(DATA.Handles.cursor(panel),'visible','on',...
      'xdata',NaN,'ydata',NaN);
    
    set(DATA.imagefig,'WindowButtonMotionFcn',...
      sprintf('%s(''manualdraw_Motion'')',mfilename))
    

    DATA.CursorX = nan(1,256);
    DATA.CursorY = nan(1,256);
    DATA.CursorN = 0;
    
    [DATA.CursorYOfs,DATA.CursorXOfs] = calcfunctions('calcoffset',slice,[],NO,panel);

    distlimit = DATA.Pref.ContourAdjustDistance;
    
    set(DATA.Handles.cursor(panel),'linestyle','-');

    switch type
      case 'endo'
      %Remove interpolation points if existing
      if ~isempty(SET(NO).EndoInterpX)
        [SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
        [SET(NO).EndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
      end
        set(DATA.Handles.cursor(panel),'color','r');
        if new
          %Calc distance to contour
          if isempty(SET(NO).EndoX)
            dist = 2*distlimit;
          else
            dist = min(sqrt(...
              (SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
              (SET(NO).EndoY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
          end;
          if dist<distlimit
            new = false;
          end;
        end;
        if not(new)
          obj=DATA.Handles.endocontour(panel);
        end;
      case 'epi'
      %Remove interpolation points if existing
      if ~isempty(SET(NO).EpiInterpX)
        [SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
        [SET(NO).EpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
      end
        set(DATA.Handles.cursor(panel),'color','g');
        if new
          %Calc distance to contour
          if isempty(SET(NO).EpiX)
            dist = 2*distlimit;
          else
            dist = min(sqrt(...
              (SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
              (SET(NO).EpiY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
          end;
          if dist<distlimit
            new = false;
          end;          
        end;
        if not(new)
          obj=DATA.Handles.epicontour(panel);
        end;
      case 'rvendo'
      if ~isempty(SET(NO).RVEndoInterpX)
        [SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
        [SET(NO).RVEndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
      end
        set(DATA.Handles.cursor(panel),'color','m');
        if new
          %Calc distance to contour
          if isempty(SET(NO).RVEndoX)
            dist = 2*distlimit;
          else
            dist = min(sqrt(...
              (SET(NO).RVEndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
              (SET(NO).RVEndoY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
          end;
          if dist<distlimit
            new = false;
          end;
        end;
        if not(new)
          obj=DATA.Handles.rvendocontour(panel);
        end;
      case 'rvepi'
      if ~isempty(SET(NO).RVEpiInterpX)
        [SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
        [SET(NO).RVEpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}] = deal([]);
      end
        set(DATA.Handles.cursor(panel),'color','c');
        if new
          %Calc distance to contour
          if isempty(SET(NO).RVEpiX)
            dist = 2*distlimit;
          else
            dist = min(sqrt(...
              (SET(NO).RVEpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-y).^2+...
              (SET(NO).RVEpiY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)-x).^2));
          end;
          if dist<distlimit
            new = false;
          end;
        end;
        if not(new)
          obj=DATA.Handles.rvepicontour(panel);
        end;
      case 'roi'
        
        %toggle hide roi off
        stateandicon=segment('iconson','hideroi');
        stateandicon{2}.isindented=0;
        stateandicon{2}.cdataDisplay=stateandicon{2}.cdata;
        DATA.Handles.configiconholder.render
        feval(stateandicon{2}.execute);
        
        set(DATA.Handles.cursor(panel),'color',DATA.GUISettings.DefaultROIDrawColor);
        DATA.Handles.phasecursor = []; %reset
        no = roi('roifindmag');
        if ~isempty(SET(NO).Flow)
          %Find panel for phase image
          if isempty(SET(NO).Flow.PhaseNo)
            phasepanel = [];
          else
            phasepanel = find(DATA.ViewPanels==SET(NO).Flow.PhaseNo);
          end;
          if ~isempty(phasepanel)
            phasepanel = phasepanel(1); %take first if many
            hold(DATA.Handles.imageaxes(phasepanel),'on');
            DATA.Handles.phasecursor = plot(DATA.Handles.imageaxes(phasepanel),...
              NaN,NaN,'b-');
            hold(DATA.Handles.imageaxes(phasepanel),'off');            
          end;
        end;
        %Calc distance to contour
        oldmin = 1e10;
        temproicurrent = SET(no).RoiN;
        for rloop=1:SET(no).RoiN
          if ~isempty(find(SET(no).Roi(rloop).Z==SET(no).CurrentSlice,1)) && ~isempty(find(SET(no).Roi(rloop).T==SET(no).CurrentTimeFrame,1))
            temp = min(sqrt(...
              (SET(no).Roi(rloop).X(:,SET(no).CurrentTimeFrame)-y).^2+...
              (SET(no).Roi(rloop).Y(:,SET(no).CurrentTimeFrame)-x).^2));
            if temp<oldmin
              temproicurrent=rloop;
              oldmin = temp;
            end;
          end;
        end;
        if new && (oldmin<distlimit)
          new = false;
        end;
        if not(new)
          SET(no).RoiCurrent = temproicurrent;
          drawfunctions('updatenopanels',no);
          obj=DATA.Handles.roicontour{panel}(SET(no).RoiCurrent);
        end;                
      case {'scar','mo','rubberpen','rubber'}
        if isempty(SET(NO).Scar)
%           if SET(NO).TSize>1
%             myfailed('Does not seem to be viability image stack.',DATA.GUI.Segment);
%             return;
%           end;
          
          %Check if doing atrial scar or normal scar
          if doatrialscar(NO)
            %Need to do nothing. Initialization is taken care of in
            %buttonup
          else
            viability('viabilityreset_Callback','manual');
            SET(NO).Scar.Mode='manual';
            viability('viabilitymenu');
          end;
        end;
        switch type
          case 'scar'
            set(DATA.Handles.cursor(panel),'color','y');
          case 'mo'
            set(DATA.Handles.cursor(panel),'color','r');            
          case 'rubberpen'
            set(DATA.Handles.cursor(panel),'color','y','linestyle',':');
          case 'rubber'
            set(DATA.Handles.cursor(panel),'color','w','linestyle','-');
        end;
      case 'drawmarpen'
        set(DATA.Handles.cursor(panel),'color','w','linestyle','-');
      case 'drawmarrubberpen'
        set(DATA.Handles.cursor(panel),'color','w','linestyle',':');
    end;

    if ~isempty(obj)
      set(obj,'LineWidth',2);
      set(DATA.imagefig,'WindowButtonUpFcn',...
        sprintf('%s(''manualdraw_Buttonup'',''%s'',%d,%16.16f )',mfilename,type,new,obj));
    end
  case 'extend'
    pan_Buttondown;
  case 'alt'
    if strcmp(type,'roi')
      roi('roiget','nomenu','coord');
      drawfunctions('drawallslices');
    end
    DATA.contextmenu;
end;

%--------------------------------------------------
function [xout,yout] = interphelper(pinx,piny)
%--------------------------------------------------
%Interpolates points to create a contour without loops.
global DATA 

%Since RV addition we need to manage nan entries here 
pinx=pinx(~isnan(pinx));
piny=piny(~isnan(piny));

xout = []; %#ok<NASGU>
yout = []; %#ok<NASGU>

if (length(pinx)<4)
  xout = NaN*ones(DATA.NumPoints,1);
  yout = xout;
  return;
end
%--- Generate output

reppinx = [pinx;pinx;pinx];
reppiny = [piny;piny;piny];

%Find distance between them
d = [0;cumsum(sqrt(diff(reppinx).^2+diff(reppiny).^2))];

%Remve duplicate points
ind = [1;diff(d)]>1e-3;
reppinx = reppinx(ind);
reppiny = reppiny(ind);
d = d(ind);

%Interpolate, periodic boundaries
xout = interp1(d,reppinx,linspace(d(length(reppinx)/3+1),d(2*length(reppinx)/3+1),DATA.NumPoints),'pchip');

yout = interp1(d,reppiny,linspace(d(length(reppinx)/3+1),d(2*length(reppinx)/3+1),DATA.NumPoints),'pchip');

[xout,yout]=mypoly2cw(xout,yout);
  
% if isopengui('strain.fig')
%   straintagging.straintagging('updateinterp',1)
%   straintagging.straintagging('updateimages')
% end

%-------------------------------------
function tool=interptoolfromcoords(x,y,slice)
%-------------------------------------
%get tool from point
%valid for 'interpendo' 'interpepi' 'interprvendo' 'interprvepi'

global SET NO;

dist_min=ones(1,4)*1e6;
ret_tool={'interpendo' 'interpepi' 'interprvendo' 'interprvepi'};

if ~isempty(SET(NO).EndoInterpX),
  [dist] = min(sqrt(...
      (SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,slice}-y).^2+...
      (SET(NO).EndoInterpY{SET(NO).CurrentTimeFrame,slice}-x).^2));
    if ~isempty(dist),
      dist_min(1)=dist;
    end
end

if ~isempty(SET(NO).EpiInterpX)
    [dist] = min(sqrt(...
      (SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,slice}-y).^2+...
      (SET(NO).EpiInterpY{SET(NO).CurrentTimeFrame,slice}-x).^2));
     if ~isempty(dist),
      dist_min(2)=dist;
     end
end

if ~isempty(SET(NO).RVEndoInterpX),
    [dist] = min(sqrt(...
      (SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,slice}-y).^2+...
      (SET(NO).RVEndoInterpY{SET(NO).CurrentTimeFrame,slice}-x).^2));
     if ~isempty(dist),
      dist_min(3)=dist;
     end
end

if ~isempty(SET(NO).RVEpiInterpX),
    [dist] = min(sqrt(...
      (SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,slice}-y).^2+...
      (SET(NO).RVEpiInterpY{SET(NO).CurrentTimeFrame,slice}-x).^2));
     if ~isempty(dist),
      dist_min(4)=dist;
     end
end

[~,ind]=min(dist_min);
tool=ret_tool{ind};

%---------------------------------------------
function interpdeletepointthisslicephase %#ok<DEFNU>
%---------------------------------------------
%Delete interp points this slice and phase.
global DATA SET NO

tools('enableundo')
tool=DATA.CurrentTool;
if ~strcmp(tool,{'interpendo' 'interpepi' 'interprvendo' 'interprvepi'}),
  [x,y,slice] = getclickedcoords;   
   tool=interptoolfromcoords(x,y,slice);
end

switch tool
  case 'interpendo'
    SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
    SET(NO).EndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
  case 'interpepi'
    SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
    SET(NO).EpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
  case 'interprvendo'
    SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
    SET(NO).RVEndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
  case 'interprvepi'
    SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
    SET(NO).RVEpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=[];
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(NO).StrainTagging) && isfield(SET(NO).StrainTagging, 'LVupdated')
  SET(NO).StrainTagging.LVupdated = 1;
end

updatemodeldisplay;
drawfunctions('drawsliceno');

%---------------------------------------------
function interpdeletepointthisslice %#ok<DEFNU>
%---------------------------------------------
%Delete interp points this slice.
global DATA SET NO

tools('enableundo')
tool=DATA.CurrentTool;
if ~strcmp(tool,{'interpendo' 'interpepi' 'interprvendo' 'interprvepi'}),
  [x,y,slice] = getclickedcoords;   
   tool=interptoolfromcoords(x,y,slice);
end

for tloop=1:SET(NO).TSize
  switch tool
    case 'interpendo'
      SET(NO).EndoInterpX{tloop,SET(NO).CurrentSlice}=[];
      SET(NO).EndoInterpY{tloop,SET(NO).CurrentSlice}=[];
    case 'interpepi'
      SET(NO).EpiInterpX{tloop,SET(NO).CurrentSlice}=[];
      SET(NO).EpiInterpY{tloop,SET(NO).CurrentSlice}=[];
    case 'interprvendo'
      SET(NO).RVEndoInterpX{tloop,SET(NO).CurrentSlice}=[];
      SET(NO).RVEndoInterpY{tloop,SET(NO).CurrentSlice}=[];
    case 'interprvepi'
      SET(NO).RVEpiInterpX{tloop,SET(NO).CurrentSlice}=[];
      SET(NO).RVEpiInterpY{tloop,SET(NO).CurrentSlice}=[];
  end
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(NO).StrainTagging) && isfield(SET(NO).StrainTagging, 'LVupdated')
  SET(NO).StrainTagging.LVupdated = 1;
end

updatemodeldisplay;
drawfunctions('drawsliceno');

%---------------------------------------------
function interpdeletepointall %#ok<DEFNU>
%---------------------------------------------
%Delete interp points for all slices timeframes
global DATA SET NO

tools('enableundo')
tool=DATA.CurrentTool;
if ~strcmp(tool,{'interpendo' 'interpepi' 'interprvendo' 'interprvepi'}),
  [x,y,slice] = getclickedcoords;   
   tool=interptoolfromcoords(x,y,slice);
end

for tloop=1:SET(NO).TSize
  for zloop=1:SET(NO).ZSize
    switch tool
      case 'interpendo'
        SET(NO).EndoInterpX{tloop,zloop}=[];
        SET(NO).EndoInterpY{tloop,zloop}=[];
      case 'interpepi'
        SET(NO).EpiInterpX{tloop,zloop}=[];
        SET(NO).EpiInterpY{tloop,zloop}=[];
      case 'interprvendo'
        SET(NO).RVEndoInterpX{tloop,zloop}=[];
        SET(NO).RVEndoInterpY{tloop,zloop}=[];
      case 'interprvepi'
        SET(NO).RVEpiInterpX{tloop,zloop}=[];
        SET(NO).RVEpiInterpY{tloop,zloop}=[];
    end
  end
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(NO).StrainTagging) && isfield(SET(NO).StrainTagging, 'LVupdated')
  SET(NO).StrainTagging.LVupdated = 1;
end

updatemodeldisplay;
drawfunctions('drawsliceno');

%---------------------------------------------
function interpdeletepoint %#ok<DEFNU>
%---------------------------------------------
%Delete interp point.
global DATA SET NO

[x,y,slice] = getclickedcoords;    
distlimit = DATA.Pref.ContourAdjustDistance;

tool=DATA.CurrentTool;
if ~strcmp(tool,{'interpendo' 'interpepi' 'interprvendo' 'interprvepi'}),
   tool=interptoolfromcoords(x,y,slice);
end

switch tool
  case 'interpendo'
    [dist,delind] = min(sqrt(...
      (SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-y).^2+...
      (SET(NO).EndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-x).^2));
    if ~(dist<distlimit)
      return;
    end
    tools('enableundo')
    numofpins=length(SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice});
    ind=true(1,numofpins);
    ind(delind)=false;
    SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
    SET(NO).EndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).EndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
  case 'interpepi'
    [dist,delind] = min(sqrt(...
      (SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-y).^2+...
      (SET(NO).EpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-x).^2));
    if ~(dist<distlimit)
      return;
    end
    tools('enableundo')
    numofpins=length(SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice});
    ind=true(1,numofpins);
    ind(delind)=false;
    SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
    SET(NO).EpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).EpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
  case 'interprvendo'
    [dist,delind] = min(sqrt(...
      (SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-y).^2+...
      (SET(NO).RVEndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-x).^2));
    if ~(dist<distlimit)
      return;
    end
    tools('enableundo')
    numofpins=length(SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice});
    ind=true(1,numofpins);
    ind(delind)=false;
    SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
    SET(NO).RVEndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).RVEndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
  case 'interprvepi'
    [dist,delind] = min(sqrt(...
      (SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-y).^2+...
      (SET(NO).RVEpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}-x).^2));
    if ~(dist<distlimit)
      return;
    end
    tools('enableundo')
    numofpins=length(SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice});
    ind=true(1,numofpins);
    ind(delind)=false;
    SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
    SET(NO).RVEpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=...
      SET(NO).RVEpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}(ind);
end

switch tool
  case 'interpendo'
    contourexists=...
         ~isempty(SET(NO).EndoX)&&...
         ~any(isnan(SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,slice)));
  case 'interpepi'
    contourexists=...
         ~isempty(SET(NO).EpiX)&&...
         ~any(isnan(SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,slice)));
  case 'interprvendo'
    contourexists=...
         ~isempty(SET(NO).RVEndoX)&&...
         ~any(isnan(SET(NO).RVEndoX(:,SET(NO).CurrentTimeFrame,slice)));
  case 'interprvepi'
    contourexists=...
         ~isempty(SET(NO).RVEpiX)&&...
         ~any(isnan(SET(NO).RVEpiX(:,SET(NO).CurrentTimeFrame,slice)));
  otherwise
    contourexists=0;
end
if contourexists
  interpdraw_Buttondown([],[],true);
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(NO).StrainTagging) && isfield(SET(NO).StrainTagging, 'LVupdated')
  SET(NO).StrainTagging.LVupdated = 1;
end

updatemodeldisplay;
drawfunctions('drawsliceno');

%-----------------------------------------------
function [x,y]=pointsfromcontour(X,Y,nbr_points)
%-----------------------------------------------
%select points from X,Y with probability higher
%when derivative magnitude in respect to distance to mass centre
%is large.

%  x_m=mean(X);
%  y_m=mean(Y);
%  r=sqrt((X-x_m).^2+(Y-y_m).^2);
% % 
% p=[1;abs(diff(r))];%./sum(abs(diff(r)));
% % P=[0 cumsum(p)'];
% 
% 
 ind=false(length(X),1);
ind(1:floor(length(X)/nbr_points):end)=true;
% ind(p(ind)<0.3)=false;
% s = rand('twister');
% rand('twister',5489);%<-- semi random for now
% ran=rand(1,nbr_points);
% rand('twister',s); %<-random seed again
% 
% nbr_points=floor(nbr_points/2);
% for i=1:nbr_points,
%   ind(sum(ran(i)<P))=true;
% end

y=Y(ind);
x=X(ind);

%--------------------------------------------------
function interpdrawGuessPoints_Callback(type,panel) %#ok<DEFNU>
%--------------------------------------------------
%Guess points in current selection then moves to corrosponding tool.
global SET NO DATA

if (nargin<2)||isempty(panel)
  panel=DATA.CurrentPanel; 
else
  switchtopanel(panel);
end

if (nargin<1)||isempty(type)
  %Annoying but easy way to get type from context menu calls. 
  % Means that you can't change another type of contour from one tool. 
  switch DATA.CurrentTool
    case 'interpendo'
      type='endo';
    case 'interpepi'
      type='epi';
    case 'interprvendo'
      type='rvendo';
    case 'interprvepi'
      type='rvepi';
    otherwise
      disp('Count not find tool. Should not happen if handles are properly assigned/reassigned.');
      return;
  end
end;

[~,~,slice] = getclickedcoords; %only care about slice
if (slice>SET(NO).ZSize)
  return;
end
switchtoslice(slice);

switch type,
  case 'endo'
    if isempty(SET(NO).EndoX)
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
    else
      if any(isnan(SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)))
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
      end
    end
    if ~isempty(SET(NO).EndoInterpX)
      if ~isempty(SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice})
        if ~yesno('There is pin content in current selection. Do you want to remove it?',[],DATA.GUI.Segment)
          return;
        end
      end
    end

    [x,y]=pointsfromcontour(...
      SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),...
      SET(NO).EndoY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),16);

    SET(NO).EndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=x;
    SET(NO).EndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=y;
    
    if size(SET(NO).EndoInterpX,1)<SET(NO).TSize||...
       size(SET(NO).EndoInterpX,2)<SET(NO).ZSize,
         SET(NO).EndoInterpX = cell(SET(NO).TSize,SET(NO).ZSize);
         SET(NO).EndoInterpY = cell(SET(NO).TSize,SET(NO).ZSize);      
    end;
  case 'epi'
    if isempty(SET(NO).EpiX)
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
    else
      if any(isnan(SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)))
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
      end
    end
    if ~isempty(SET(NO).EpiInterpX)
      if ~isempty(SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice})
        if ~yesno('There is pin content in current selection. Do you want to remove it?',[],DATA.GUI.Segment)
          return;
        end
      end
    end

    [x,y]=pointsfromcontour(...
      SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),...
      SET(NO).EpiY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),16);

    SET(NO).EpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=x;
    SET(NO).EpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=y;
    
    if size(SET(NO).EpiInterpX,1)<SET(NO).TSize||...
       size(SET(NO).EpiInterpX,2)<SET(NO).ZSize,
         SET(NO).EpiInterpX = cell(SET(NO).TSize,SET(NO).ZSize);
         SET(NO).EpiInterpY = cell(SET(NO).TSize,SET(NO).ZSize);      
    end;
  case 'rvendo'
    if isempty(SET(NO).RVEndoX)
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
    else
      if any(isnan(SET(NO).RVEndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)))
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
      end
    end
    if ~isempty(SET(NO).RVEndoInterpX)
      if ~isempty(SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice})
        if ~yesno('There is pin content in current selection. Do you want to remove it?',[],DATA.GUI.Segment)
          return;
        end
      end
    end

    [x,y]=pointsfromcontour(...
      SET(NO).RVEndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),...
      SET(NO).RVEndoY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),16);

    SET(NO).RVEndoInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=x;
    SET(NO).RVEndoInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=y;
    
    if size(SET(NO).RVEndoInterpX,1)<SET(NO).TSize||...
       size(SET(NO).RVEndoInterpX,2)<SET(NO).ZSize,
         SET(NO).RVEndoInterpX = cell(SET(NO).TSize,SET(NO).ZSize);
         SET(NO).RVEndoInterpY = cell(SET(NO).TSize,SET(NO).ZSize);      
    end;
  case 'rvepi'
   if isempty(SET(NO).RVEpiX)
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
    else
      if any(isnan(SET(NO).RVEpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)))
        mywarning('Could not find contour to base guess on. Aborting!',DATA.GUI.Segment);
        return;
      end
    end
    if ~isempty(SET(NO).RVEpiInterpX)
      if ~isempty(SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice})
        if ~yesno('There is pin content in current selection. Do you want to remove it?',[],DATA.GUI.Segment)
          return;
        end
      end
    end

    [x,y]=pointsfromcontour(...
      SET(NO).RVEpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),...
      SET(NO).RVEpiY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice),16);

    SET(NO).RVEpiInterpX{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=x;
    SET(NO).RVEpiInterpY{SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice}=y;
    
    if size(SET(NO).RVEpiInterpX,1)<SET(NO).TSize||...
       size(SET(NO).RVEpiInterpX,2)<SET(NO).ZSize,
         SET(NO).RVEpiInterpX = cell(SET(NO).TSize,SET(NO).ZSize);
         SET(NO).RVEpiInterpY = cell(SET(NO).TSize,SET(NO).ZSize);      
    end;
end

interpdraw_Buttondown(panel,type,true);

%-----------------------------------------------------------------------------------
function [interpx,interpy] = interpdrawhelper(interpx,interpy,contourx,contoury,x,y)
%-----------------------------------------------------------------------------------
%Helper fcn to interpdrawbuttonup
%Side effects is update of DATA.Pin.
global DATA SET NO

new=true;
tf = SET(NO).CurrentTimeFrame;
slice = SET(NO).CurrentSlice;
distlimit = DATA.Pref.ContourAdjustDistance;

if ~isempty(interpx)&&~isempty(interpx{tf,slice})
  [dist,tempind] = min(sqrt(...
    (interpx{tf,slice}-y).^2+...
    (interpy{tf,slice}-x).^2));
  if dist<distlimit
    new = false;
    DATA.Pin = tempind; %Store what pin
  end;
end;

if ~new
  return
end

%New interpolation point. User has not clicked close to old point.

%Prepare to store point
if isempty(interpx)||...
    size(interpx,1)<SET(NO).TSize||...
    size(interpx,2)<SET(NO).ZSize,
  interpx = cell(SET(NO).TSize,SET(NO).ZSize);
  interpy = cell(SET(NO).TSize,SET(NO).ZSize);
end;

newcontourx = [];
newcontoury = [];
origx=x;
origy=y;
if (~isempty(contourx)) && (~isnan(contourx(1,tf,slice)))
  %Segmentation exist.
  if isempty(interpx{tf,slice})
    %Contour exist but no previous points>=add many points.
    %Interpolate points
    newcontourx = contourx(:,tf,slice);
    newcontoury = contoury(:,tf,slice);
    n = size(contourx,1);
    nstep = n/DATA.Pref.NumInterpPoints; %Later take from preferences.
    newcontourx = newcontourx(round(1:nstep:n));
    newcontoury = newcontoury(round(1:nstep:n));
    newcontourx = newcontourx(:);
    newcontoury = newcontoury(:);
    y = []; %Do not add the point
    x = []; %Do not add the point
  else
    %Contour and previous points exist => find closest gap
    pinx = [interpx{tf,slice}; y];
    piny = [interpy{tf,slice}; x];
    contx = contourx(:,tf,slice);
    conty = contoury(:,tf,slice);
    
    pinxrep=repmat(pinx',[length(contx) 1]);
    contxrep=repmat(contx,[1 length(pinx)]);
    pinyrep=repmat(piny',[length(conty) 1]);
    contyrep=repmat(conty,[1 length(piny)]);
    pindist2cont = (pinxrep-contxrep).^2+(pinyrep-contyrep).^2;
    [~,mindistindex] =min(pindist2cont);
    [~,sortindex] =sort(mindistindex);
    pinx=pinx(sortindex);
    piny=piny(sortindex);
    
    interpx{tf,slice}=pinx;
    interpy{tf,slice}=piny;
    
    y = []; %Do not add the point
    x = []; %Do not add the point
  end;
end;

%Add points
interpx{tf,slice} = [...
  interpx{tf,slice} ; ...
  newcontourx;y];
interpy{tf,slice} = [...
  interpy{tf,slice} ; ...
  newcontoury;x];

%Store what pin
pinx=interpx{tf,slice};
piny=interpy{tf,slice};
[~,index]=min((pinx-origy).^2+(piny-origx).^2);
DATA.Pin = index;

%---------------------------------------------
function interpdraw_Buttondown(panel,type,forcedraw) 
%---------------------------------------------
%Button down function to draw interp points.
global DATA SET NO
%profile on

if DATA.Interactionlock
  return;
end;

if (nargin<1)||isempty(panel)
  panel=DATA.CurrentPanel; %#ok<NASGU>
  killbuttondown = 0;
else
  killbuttondown = switchtopanel(panel);
end

if killbuttondown
  return
end

if isequal(DATA.CurrentTool,'select')
  return;
end

no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = true;
% end

if (nargin<2)||isempty(type)
  %Annoying but easy way to get type from context menu calls. 
  % Means that you can't change another type of contour from one tool. 
  switch DATA.CurrentTool
    case 'interpendo'
      type='endo';
    case 'interpepi'
      type='epi';
    case 'interprvendo'
      type='rvendo';
    case 'interprvepi'
      type='rvepi';
    otherwise
      %shouldn't happen, if handles are properly assigned/reassigned.
      return;
  end
end;

if nargin<3
  forcedraw=false;
end

[x,y,slice] = getclickedcoords;    
% If slice has changed, make sure montage/one are in sync
if (slice>SET(no).ZSize)
  return;
end
switchtoslice(slice);

%Check what type of click
if ~(forcedraw)
  switcher=get(DATA.imagefig,'SelectionType');
  tools('enableundo') %only push undo stack when not forcedrawn
else
  switcher = 'extend';
end

switch switcher
  case 'normal'
    
    % Distance check. Check to see if we haven't clicked close to something

    nomotion=false; 

    %Add new points depending if contour exists etc.
    switch type
      case 'endo'
        [SET(no).EndoInterpX,SET(no).EndoInterpY] = interpdrawhelper(...
          SET(no).EndoInterpX,...
          SET(no).EndoInterpY,...
          SET(no).EndoX,...
          SET(no).EndoY,x,y);
      case 'epi'
        [SET(no).EpiInterpX,SET(no).EpiInterpY] = interpdrawhelper(...
          SET(no).EpiInterpX,...
          SET(no).EpiInterpY,...
          SET(no).EpiX,...
          SET(no).EpiY,x,y);
      case 'rvendo'
        [SET(no).RVEndoInterpX,SET(no).RVEndoInterpY] = interpdrawhelper(...
          SET(no).RVEndoInterpX,...
          SET(no).RVEndoInterpY,...
          SET(no).RVEndoX,...
          SET(no).RVEndoY,x,y);           
      case 'rvepi'
        [SET(no).RVEpiInterpX,SET(no).RVEpiInterpY] = interpdrawhelper(...
          SET(no).RVEpiInterpX,...
          SET(no).RVEpiInterpY,...
          SET(no).RVEpiX,...
          SET(no).RVEpiY,x,y);        
    end;
    
    if ~(nomotion)
      set(DATA.imagefig,'WindowButtonMotionFcn',...
        sprintf('%s(''interppointMotion'',''%s'')',mfilename,type));
      set(DATA.imagefig,'WindowButtonUpFcn',...
        sprintf('%s(''interppointButtonup'',''%s'')',mfilename,type));
    end
  case 'extend' %shift click
    switch type
      case 'endo'
        contourexists=...
          ~isempty(SET(no).EndoX)&&...
          ~any(isnan(SET(no).EndoX(:,SET(no).CurrentTimeFrame,slice)));
        interpcellexist = ~isempty(SET(no).EndoInterpX);
      case 'epi'
        contourexists=...
          ~isempty(SET(no).EpiX)&&...
          ~any(isnan(SET(no).EpiX(:,SET(no).CurrentTimeFrame,slice)));
        interpcellexist = ~isempty(SET(no).EpiInterpX);
      case 'rvendo'
        contourexists=...
          ~isempty(SET(no).RVEndoX)&&...
          ~any(isnan(SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,slice)));
        interpcellexist = ~isempty(SET(no).RVEndoInterpX);
      case 'rvepi'
        contourexists=...
          ~isempty(SET(no).RVEpiX)&&...
          ~any(isnan(SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,slice)));
        interpcellexist = ~isempty(SET(no).RVEpiInterpX);
    end
    
    if ~interpcellexist && ~contourexists
      return;
    end

    switch type
      case 'endo'
        if isempty(SET(no).EndoX)
          SET(no).EndoX = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
          SET(no).EndoY = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
        end;
        [xout,yout] = interphelper(...
          SET(no).EndoInterpX{SET(no).CurrentTimeFrame,slice},...
          SET(no).EndoInterpY{SET(no).CurrentTimeFrame,slice});
        if ~isempty(xout)
          SET(no).EndoX(:,SET(no).CurrentTimeFrame,slice) = xout;
          SET(no).EndoY(:,SET(no).CurrentTimeFrame,slice) = yout;
        end
      case 'epi'
        if isempty(SET(no).EpiX)
          SET(no).EpiX = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
          SET(no).EpiY = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
        end;
        [xout,yout] = interphelper(...
          SET(no).EpiInterpX{SET(no).CurrentTimeFrame,slice},...
          SET(no).EpiInterpY{SET(no).CurrentTimeFrame,slice});
        if ~isempty(xout)
         SET(no).EpiX(:,SET(no).CurrentTimeFrame,slice) = xout;
         SET(no).EpiY(:,SET(no).CurrentTimeFrame,slice) = yout;
        end
      case 'rvendo'
        if isempty(SET(no).RVEndoX)
          SET(no).RVEndoX = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
          SET(no).RVEndoY = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
        end;
        [xout,yout] = interphelper(...
          SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,slice},...
          SET(no).RVEndoInterpY{SET(no).CurrentTimeFrame,slice});
        if ~isempty(xout) 
          SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,slice) = xout;
          SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,slice) = yout;
        end
      case 'rvepi'
        if isempty(SET(no).RVEpiX)
          SET(no).RVEpiX = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
          SET(no).RVEpiY = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
        end;
        [xout,yout] = interphelper(...
          SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,slice},...
          SET(no).RVEpiInterpY{SET(no).CurrentTimeFrame,slice});
        if ~isempty(xout)
          SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,slice) = xout;
          SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,slice) = yout;
        end
      otherwise
        %nothing
    end;
    %Check if should copy to all timeframes. Contour is copied if:
    %  - SingleFrameOnly is actived
    %  - There are no interpolation points in any other other timeframes.
    switch type
      case 'endo'
        isemptymatrix=cellfun(@isempty,SET(no).EndoInterpX);
        if not(DATA.ThisFrameOnly)&&...
           (sum(~isemptymatrix(:,SET(no).CurrentSlice))<2)
          %Copy segmentation to all slices
          for tloop=1:SET(no).TSize
            SET(no).EndoX(:,tloop,SET(no).CurrentSlice) = SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
            SET(no).EndoY(:,tloop,SET(no).CurrentSlice) = SET(no).EndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
          end;
        end;
      case 'epi'
        isemptymatrix=cellfun(@isempty,SET(no).EpiInterpX);
        if not(DATA.ThisFrameOnly)&&...
           (sum(~isemptymatrix(:,SET(no).CurrentSlice))<2)
          %Copy segmentation to all slices
          for tloop=1:SET(no).TSize
            SET(no).EpiX(:,tloop,SET(no).CurrentSlice) = SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
            SET(no).EpiY(:,tloop,SET(no).CurrentSlice) = SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
          end;
        end;
      case 'rvendo'
        isemptymatrix=cellfun(@isempty,SET(no).RVEndoInterpX);
        if not(DATA.ThisFrameOnly)&&...
           (sum(~isemptymatrix(:,SET(no).CurrentSlice))<2)
          %Copy segmentation to all slices
          for tloop=1:SET(no).TSize
            SET(no).RVEndoX(:,tloop,SET(no).CurrentSlice) = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
            SET(no).RVEndoY(:,tloop,SET(no).CurrentSlice) = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
          end;
        end;    
      case 'rvepi'
        isemptymatrix=cellfun(@isempty,SET(no).RVEpiInterpX);
        if not(DATA.ThisFrameOnly)&&...
           (sum(~isemptymatrix(:,SET(no).CurrentSlice))<2)
          %Copy segmentation to all slices
          for tloop=1:SET(no).TSize
            SET(no).RVEpiX(:,tloop,SET(no).CurrentSlice) = SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
            SET(no).RVEpiY(:,tloop,SET(no).CurrentSlice) = SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
          end;
        end;    
    end;
    lvsegchanged = true; 
    updatevolume(lvsegchanged);
  case 'alt'
    DATA.contextmenu;
end;

%Test
if ~isempty(SET(no).StrainTagging)
  SET(no).StrainTagging.LVupdated=1;
end

checkconsistency;
updateselectedslices;
updatemodeldisplay;
drawfunctions('drawsliceno');

%--------------------
function interppointButtonup(type) %#ok<DEFNU>
%--------------------
%Button up function for interp points.
global DATA SET NO

no = NO;
if ~isempty(SET(NO).Parent)
 no = SET(NO).Parent;
end;

switch type
  case 'endo'
    contourexists=...
         ~isempty(SET(no).EndoX)&&...
         ~any(isnan(SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)));
  case 'epi'
    contourexists=...
         ~isempty(SET(no).EpiX)&&...
         ~any(isnan(SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)));
  case 'rvendo'
    contourexists=...
         ~isempty(SET(no).RVEndoX)&&...
         ~any(isnan(SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)));
  case 'rvepi'
    contourexists=...
         ~isempty(SET(no).RVEpiX)&&...
         ~any(isnan(SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)));
end

checkconsistency;
set(DATA.imagefig,'WindowButtonMotionFcn','');
if contourexists
  interpdraw_Buttondown([],[],true);
end
%set(DATA.fig,'WindowButtonMotionFcn',@DATA.toggleplaceholdermotion)
updatemodeldisplay;
updatevolume;
drawfunctions('drawallslices');
%profile report

%---------------------------------------------
function interppointMotion(type) %#ok<DEFNU>
%---------------------------------------------
%Motion function for interp points.
global DATA SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

%Setup
[x,y,slice] = getclickedcoords;
%extra safety check
if slice>SET(NO).ZSize
  slice=SET(NO).CurrentSlice;
end

switch type
  case 'endo'
    if isempty(SET(no).EndoInterpX)
      return;
    end;
    temp = SET(no).EndoInterpX{SET(no).CurrentTimeFrame,slice};
    if isempty(temp)
      return
    end
    temp(DATA.Pin,1) = y; %Make sure column vector
    SET(no).EndoInterpX{SET(no).CurrentTimeFrame,slice} = temp;
    temp = SET(no).EndoInterpY{SET(no).CurrentTimeFrame,slice};
    temp(DATA.Pin,1) = x;
    SET(no).EndoInterpY{SET(no).CurrentTimeFrame,slice} = temp;
    if ~isempty(SET(no).EndoX)&&...
        ~all(isnan(SET(no).EndoX(:,SET(no).CurrentTimeFrame,slice)))
      [xout,yout] = interphelper(...
        SET(no).EndoInterpX{SET(no).CurrentTimeFrame,slice},...
        SET(no).EndoInterpY{SET(no).CurrentTimeFrame,slice});
      if ~isempty(xout)
        SET(no).EndoX(:,SET(no).CurrentTimeFrame,slice) = xout;
        SET(no).EndoY(:,SET(no).CurrentTimeFrame,slice) = yout;
      end
    end
  case 'epi'
    if isempty(SET(no).EpiInterpX)
      return;
    end;    
    temp = SET(no).EpiInterpX{SET(no).CurrentTimeFrame,slice};
    if isempty(temp)
      return
    end
    temp(DATA.Pin,1) = y;
    SET(no).EpiInterpX{SET(no).CurrentTimeFrame,slice} = temp;
    temp = SET(no).EpiInterpY{SET(no).CurrentTimeFrame,slice};
    temp(DATA.Pin,1) = x;
    SET(no).EpiInterpY{SET(no).CurrentTimeFrame,slice} = temp;
    if ~isempty(SET(no).EpiX)&&...
        ~all(isnan(SET(no).EpiX(:,SET(no).CurrentTimeFrame,slice)))
      [xout,yout] = interphelper(...
        SET(no).EpiInterpX{SET(no).CurrentTimeFrame,slice},...
        SET(no).EpiInterpY{SET(no).CurrentTimeFrame,slice});
      if ~isempty(xout)
        SET(no).EpiX(:,SET(no).CurrentTimeFrame,slice) = xout;
        SET(no).EpiY(:,SET(no).CurrentTimeFrame,slice) = yout;
      end
    end
  case 'rvendo'
    if isempty(SET(no).RVEndoInterpX)
      return;
    end;        
    temp = SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,slice};
    if isempty(temp)
      return
    end
    temp(DATA.Pin,1) = y;
    SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,slice} = temp;
    temp = SET(no).RVEndoInterpY{SET(no).CurrentTimeFrame,slice};
    temp(DATA.Pin,1) = x;
    SET(no).RVEndoInterpY{SET(no).CurrentTimeFrame,slice} = temp;
    if ~isempty(SET(no).RVEndoX)&&...
        ~all(isnan(SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      [xout,yout] = interphelper(...
        SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,slice},...
        SET(no).RVEndoInterpY{SET(no).CurrentTimeFrame,slice});
      if ~isempty(xout)
        SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,slice) = xout;
        SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,slice) = yout;
      end
    end
  case 'rvepi'
    if isempty(SET(no).RVEpiInterpX)
      return;
    end;      
    temp = SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,slice};
    if isempty(temp)
      return
    end
    temp(DATA.Pin,1) = y;
    SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,slice} = temp;
    temp = SET(no).RVEpiInterpY{SET(no).CurrentTimeFrame,slice};
    temp(DATA.Pin,1) = x;
    SET(no).RVEpiInterpY{SET(no).CurrentTimeFrame,slice} = temp;
    if ~isempty(SET(no).RVEpiX)&&...
        ~all(isnan(SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,slice)))
      [xout,yout] = interphelper(...
        SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,slice},...
        SET(no).RVEpiInterpY{SET(no).CurrentTimeFrame,slice});
      if ~isempty(xout)
        SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,slice) = xout;
        SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,slice) = yout;
      end
    end
end;

if ~all(strcmp(DATA.ViewPanelsType(DATA.ViewPanels == no),'one'))
  updatemodeldisplay;
end
drawfunctions('drawsliceno');

%------------------------------------
function scale_Buttondown(type,panel) %#ok<DEFNU>
%------------------------------------
%Button down function for scaling of objects / contours.
global DATA SET NO

if nargin<1
  return;
end;

if nargin<2
  panel = DATA.CurrentPanel;
end;

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

[~,~,slice] = getclickedcoords;
% If slice has changed, make sure montage/one are in sync
if (slice>SET(no).ZSize)
  return;
end
switchtoslice(slice);

tools('enableundo')
switchtopanel(panel);

[DATA.CursorYOfs,DATA.CursorXOfs] = calcfunctions('calcoffset',slice,[],no,panel);

switch get(DATA.imagefig,'SelectionType')
  case 'alt'
    return;
  case 'normal'
    switch type
      case 'image'
        return;
      case 'endo'
        if isempty(SET(no).EndoX)
          myfailed('No LV endocardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        y = SET(no).EndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);    
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','r');            
      case 'epi'
        if isempty(SET(no).EpiX)
          myfailed('No LV epicardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        y = SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);    
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','g');
      case 'rvendo'
        if isempty(SET(no).RVEndoX)
          myfailed('No RV endocardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        y = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);    
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','m');        
      case 'rvepi'
        if isempty(SET(no).RVEpiX)
          myfailed('No RV epicardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        y = SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);    
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','c');    
      case 'roi'
        %Find closest ROI

        if SET(no).RoiN<1
          myfailed('No ROIs available.',DATA.GUI.Segment);
          return;
        end;

        [ycl,xcl] = getclickedcoords;
        mindist = 1e10;
        n = NaN;
        for loop=1:SET(no).RoiN 
          if ~isempty(find(SET(no).Roi(loop).Z==SET(no).CurrentSlice,1))&&...
              ~isempty(find(SET(no).Roi(loop).T==SET(no).CurrentTimeFrame,1))
            xr = SET(no).Roi(loop).X(:,SET(no).CurrentTimeFrame);
            yr = SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame);
            thisdist = min(sqrt((xr-xcl).^2+(yr-ycl).^2));
            if thisdist<mindist
              n = loop;
              mindist = thisdist;
            end;
          end
        end;

        if isnan(n)
          myfailed('No ROIs available.',DATA.GUI.Segment);
          return;
        end;

        SET(no).RoiCurrent = n;
        %Extract data
        x = SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame);
        y = SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame);

        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','b');
    end;
  otherwise
    return; %'open', 'extended'
end;

DATA.CursorX = x;
DATA.CursorY = y;
set(DATA.Handles.cursor(DATA.CurrentPanel),...
  'xdata',DATA.CursorY+DATA.CursorXOfs,...
  'ydata',DATA.CursorX+DATA.CursorYOfs,...
  'visible','on');

set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('segment(''scale_Buttonup'',''%s'')',type));
motionfcn = 'segment(''scale_Motion'');'; %get(DATA.imagefig,'WindowButtonMotionFcn')];
set(DATA.imagefig,'WindowButtonMotionFcn',motionfcn);

%---------------------------
function scale_Motion(reset) %#ok<INUSD>
%---------------------------
%Motion function for scaling.
persistent startslice startrad
global DATA SET NO

if nargin==1
  startrad = [];
  return;
end;

if isempty(startrad)

  %Get coordinate 
  [starty,startx,startslice] = getclickedcoords;
  
  mx = mean(DATA.CursorX);
  my = mean(DATA.CursorY);
   
  startrad = sqrt((mx-startx).^2+(my-starty).^2);
  
end

%Get coordinate 
[y,x,slice] = getclickedcoords;

mx = mean(DATA.CursorX);
my = mean(DATA.CursorY);

rad = sqrt((mx-x).^2+(my-y).^2);
f = rad/startrad;

%If different slice then buttonup
if (slice~=startslice) ...
    || (~isempty(DATA.CursorX) && ...
     (any((mx+(DATA.CursorX-mx)*f)<0.5)||any((mx+(DATA.CursorX-mx)*f)>SET(NO).YSize))) ...
    || (~isempty(DATA.CursorY) && ...
     (any((my+(DATA.CursorY-my)*f)<0.5)||any((my+(DATA.CursorY-my)*f)>SET(NO).XSize)))
  scale_Buttonup('abort');
end;

set(DATA.Handles.cursor(DATA.CurrentPanel),...
  'ydata',mx+DATA.CursorYOfs+(DATA.CursorX-mx)*f,...
  'xdata',my+DATA.CursorXOfs+(DATA.CursorY-my)*f);    

%----------------------------
function scale_Buttonup(type)
%----------------------------
%Button up function for scaling
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

DATA.buttonup_Callback;
set(DATA.imagefig,'WindowButtonUpFcn',sprintf('%s(''buttonup_Callback'')','segment'));
set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = true;
% end
if DATA.ThisFrameOnly && ~strcmp(type,'roi')
  tf = SET(no).CurrentTimeFrame;
else
  tf = 1:SET(no).TSize;
end;

slice = SET(no).StartSlice:SET(no).EndSlice;

%Extract the current data
x = get(DATA.Handles.cursor(DATA.CurrentPanel),'ydata')-DATA.CursorYOfs;
y = get(DATA.Handles.cursor(DATA.CurrentPanel),'xdata')-DATA.CursorXOfs;

x = x(1); %one point is enough
y = y(1);

%Find mean (same for both since scaling)
mx = mean(DATA.CursorX);
my = mean(DATA.CursorY);

startrad = sqrt((mx-DATA.CursorX(1)).^2+(my-DATA.CursorY(1)).^2);
rad = sqrt((mx-x).^2+(my-y).^2);

%Calculate factor
f = rad/startrad;
pinwarn=false;
switch type
  case 'endo'
    SET(no).EndoX(:,tf,slice) = mx+f*(SET(no).EndoX(:,tf,slice)-mx);
    SET(no).EndoY(:,tf,slice) = my+f*(SET(no).EndoY(:,tf,slice)-my);
    if ~isempty(SET(no).EndoPinX)
      [SET(no).EndoPinX,SET(no).EndoPinY,pinwarn]=pinresize(...
        SET(no).EndoPinX,SET(no).EndoPinY,tf,slice,[mx my],f);
    end
    if ~isempty(SET(no).EndoInterpX)
      [SET(no).EndoInterpX,SET(no).EndoInterpY,pinwarn]=pinresize(...
        SET(no).EndoInterpX,SET(no).EndoInterpY,tf,slice,[mx my],f);
    end
  case 'epi'
    SET(no).EpiX(:,tf,slice) = mx+f*(SET(no).EpiX(:,tf,slice)-mx);
    SET(no).EpiY(:,tf,slice) = my+f*(SET(no).EpiY(:,tf,slice)-my);        
    if ~isempty(SET(no).EpiPinX)
      [SET(no).EpiPinX,SET(no).EpiPinY,pinwarn]=pinresize(...
        SET(no).EpiPinX,SET(no).EpiPinY,tf,slice,[mx my],f);
    end
    if ~isempty(SET(no).EpiInterpX)
      [SET(no).EpiInterpX,SET(no).EpiInterpY,pinwarn]=pinresize(...
        SET(no).EpiInterpX,SET(no).EpiInterpY,tf,slice,[mx my],f);
    end
  case 'rvendo'
    SET(no).RVEndoX(:,tf,slice) = mx+f*(SET(no).RVEndoX(:,tf,slice)-mx);
    SET(no).RVEndoY(:,tf,slice) = my+f*(SET(no).RVEndoY(:,tf,slice)-my);        
    if ~isempty(SET(no).RVEndoPinX)
      [SET(no).RVEndoPinX,SET(no).RVEndoPinY,pinwarn]=pinresize(...
        SET(no).RVEndoPinX,SET(no).RVEndoPinY,tf,slice,[mx my],f);
    end
    if ~isempty(SET(no).RVEndoInterpX)
      [SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,pinwarn]=pinresize(...
        SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,tf,slice,[mx my],f);
    end
  case 'rvepi'
    SET(no).RVEpiX(:,tf,slice) = mx+f*(SET(no).RVEpiX(:,tf,slice)-mx);
    SET(no).RVEpiY(:,tf,slice) = my+f*(SET(no).RVEpiY(:,tf,slice)-my);            
    if ~isempty(SET(no).RVEpiPinX)
      [SET(no).RVEpiPinX,SET(no).RVEpiPinY,pinwarn]=pinresize(...
        SET(no).RVEpiPinX,SET(no).RVEpiPinY,tf,slice,[mx my],f);
    end
    if ~isempty(SET(no).RVEpiInterpX)
      [SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,pinwarn]=pinresize(...
        SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,tf,slice,[mx my],f);
    end
  case 'roi'
    SET(no).Roi(SET(no).RoiCurrent).X(:,tf) = mx+f*(SET(no).Roi(SET(no).RoiCurrent).X(:,tf)-mx);
    SET(no).Roi(SET(no).RoiCurrent).Y(:,tf) = my+f*(SET(no).Roi(SET(no).RoiCurrent).Y(:,tf)-my);
    roi('roiforceapply'); %ensure that all rois are equal size if this option is checked.
    [~,SET(no).Roi(SET(no).RoiCurrent).Area] = ...
      calcfunctions('calcroiarea',no,SET(no).RoiCurrent);
    [m,sd]=calcfunctions('calcroiintensity',no,SET(no).RoiCurrent);
    SET(no).Roi(SET(no).RoiCurrent).Mean = m;
    SET(no).Roi(SET(no).RoiCurrent).StD = sd;
    
    if ~isempty(DATA.FlowNO) && ismember(DATA.FlowNO,SET(no).Linked) && ~isempty(DATA.FlowROI)
      if SET(no).RoiN <1
        DATA.FlowROI = [];
      else
        calcfunctions('calcflow',no);
      end
    end
    
    DATA.updateaxestables('area',no);
    
  case 'abort'
    beep;
    DATA.CursorX = [];
    DATA.CursorY = [];
    scale_Motion('reset');
    set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
    return;
end;

if (pinwarn)
  mywarning('Pins moved outside frame deleted.',DATA.GUI.Segment);
end

%updatemodeldisplay(no);
%drawimagepanel;
%set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');

%Reset
DATA.CursorX = [];
DATA.CursorY = [];
scale_Motion('reset');

updatemodeldisplay;
drawfunctions('drawallslices');

updatevolume;
DATA.updateaxestables('t2star');
%-------------------------------
function hide(names,states)
%--------------------------
global DATA
stateandicon = segment('iconson',names);

for i =1:size(stateandicon,1)
  stateandicon{i,2}.isindented=states(i);
end

DATA.Handles.configiconholder.render;
drawfunctions('updatevisibility');

%-------------------------------------------
    function state = iconson(name)
    %------------------------------------------
    %if given name of button return state if run with no input returns
    %states of all buttons. if given cell with multiple button return icons
    %in order of request.
    global DATA
    
    state=0;
    icons=[DATA.Handles.permanenticonholder.iconCell{:},DATA.Icons.lviconcell{:},DATA.Icons.rviconcell{:},DATA.Icons.roiflowiconcell{:},DATA.Icons.viabilityiconcell{:},...
      DATA.Icons.analysisiconcell{:}, DATA.Icons.imageiconcell{:},DATA.Icons.hidecell{:}];
    N=length(icons);
    if nargin==1
      %return icon state
    if ~iscell(name)
      ind=find(strcmp(name,{icons.name}));
      if length(ind)>1
      ind=ind(1);
      end
      
      if isempty(ind)
        state={0,nan};
      else
        state={icons(ind).isindented,icons(ind)};
      end
      
    else
      n=length(name);
      state=cell(n,2);
      %Get all icons and states for names in cell
      counter=1;
      indlist=zeros(1,n);
      for i=1:n
        ind=find(strcmp(name(i),{icons.name}));
        
        if isempty(ind)
          ind=nan;
        end
        
        if length(ind)>1
          ind=ind(1);
        end
        
        indlist(i)=ind;
      end
      
      for i=1:n
        if isnan(indlist(i))
          state{i,1} = 0;
          state{i,2} = nan;
        else
          state{i,1} = icons(indlist(i)).isindented;
          state{i,2} = icons(indlist(i));
        end
       end
     end
    else
      state=cell(N,2);
       for i=1:N
         state{i,1} = icons{i}.name;
         state{i,2} = icons{i}.isindented;
       end
    end


%------------------------------------
function dragepi_Buttondown(type,panel) %#ok<DEFNU>
%------------------------------------
%Button down function for scaling of objects / contours.
global DATA SET NO

if nargin<1
  return;
end;

if nargin<2
  panel = DATA.CurrentPanel;
end;

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

[~,~,slice] = getclickedcoords;
% If slice has changed, make sure montage/one are in sync
if (slice>SET(no).ZSize)
  return;
end
switchtoslice(slice);

tools('enableundo')
switchtopanel(panel);

[DATA.CursorYOfs,DATA.CursorXOfs] = calcfunctions('calcoffset',slice,[],no,panel);

switch get(DATA.imagefig,'SelectionType')
  case 'alt'
    return;
  case 'normal'
    switch type
      case 'image'
        return;
      case {'one','allslices'}
        if isempty(SET(no).EpiX)
          myfailed('No LV epicardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
        y = SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);    
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','g');
    end
  otherwise
    return;
end;

DATA.CursorX = x;
DATA.CursorY = y;
set(DATA.Handles.cursor(DATA.CurrentPanel),...
  'xdata',DATA.CursorY+DATA.CursorXOfs,...
  'ydata',DATA.CursorX+DATA.CursorYOfs,...
  'visible','on');

set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('segment(''dragepi_Buttonup'',''%s'')',type));
motionfcn = ['segment(''dragepi_Motion'');' get(DATA.imagefig,'WindowButtonMotionFcn')];
set(DATA.imagefig,'WindowButtonMotionFcn',motionfcn);

%---------------------------
function dragepi_Motion(reset) %#ok<INUSD>
%---------------------------
%Motion function for scaling.
persistent startslice startrad
global DATA SET NO

if nargin==1
  startrad = [];
  return;
end;

if isempty(startrad)
  %Get coordinate 
  [starty,startx,startslice] = getclickedcoords;  
  mx = mean(DATA.CursorX);
  my = mean(DATA.CursorY);   
  startrad = sqrt((mx-startx).^2+(my-starty).^2);  
end

%Get coordinate 
[y,x,slice] = getclickedcoords;

mx = mean(DATA.CursorX);
my = mean(DATA.CursorY);

rad = sqrt((mx-x).^2+(my-y).^2);
f = rad/startrad;

%If different slice then buttonup
if (slice~=startslice) ...
    || (~isempty(DATA.CursorX) && ...
     (any((mx+(DATA.CursorX-mx)*f)<0.5)||any((mx+(DATA.CursorX-mx)*f)>SET(NO).YSize))) ...
    || (~isempty(DATA.CursorY) && ...
     (any((my+(DATA.CursorY-my)*f)<0.5)||any((my+(DATA.CursorY-my)*f)>SET(NO).XSize)))
  dragepi_Buttonup('abort');
end;

set(DATA.Handles.cursor(DATA.CurrentPanel),...
  'ydata',mx+DATA.CursorYOfs+(DATA.CursorX-mx)*f,...
  'xdata',my+DATA.CursorXOfs+(DATA.CursorY-my)*f);    

%----------------------------
function dragepi_Buttonup(type)
%----------------------------
%Button up function for scaling
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

DATA.buttonup_Callback;
set(DATA.imagefig,'WindowButtonUpFcn',sprintf('%s(''buttonup_Callback'')','segment'));
set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');

tf = SET(no).CurrentTimeFrame;

%Extract the current data
x = get(DATA.Handles.cursor(DATA.CurrentPanel),'ydata')-DATA.CursorYOfs;
y = get(DATA.Handles.cursor(DATA.CurrentPanel),'xdata')-DATA.CursorXOfs;

x = x(1); %one point is enough
y = y(1);

%Find mean (same for both since scaling)
mx = mean(DATA.CursorX);
my = mean(DATA.CursorY);

startrad = sqrt((mx-DATA.CursorX(1)).^2+(my-DATA.CursorY(1)).^2);
rad = sqrt((mx-x).^2+(my-y).^2);

%Calculate factor
f = rad/startrad;
meanres = mean([SET(no).ResolutionX, SET(no).ResolutionY]);
%expand the epicardium (not in the outflow tract (wallthickness = 0 in the 30% most basal LV slices))
switch type
  case {'one','allslices'}
    switch type
      case 'one'
        slices = SET(no).CurrentSlice;
      case 'allslices'
        startslice = find(~isnan(SET(no).EpiX(1,tf,:)),1,'first');
        endslice = find(~isnan(SET(no).EpiX(1,tf,:)),1,'last');
        slices = startslice:endslice;
        if length(slices) > 2
          mostbasalslice = startslice+0.3*length(slices);
        else
          mostbasalslice = endslice;
        end
    end
    for slice = slices
      %calculate wallthickness
      switch type
        case 'one'
          dowallthickness = existfunctions('existendoinselected',no,tf,slice);
        case 'allslices'
          dowallthickness = (existfunctions('existendoinselected',no,tf,slice) && slice < mostbasalslice);
      end
      if dowallthickness
        wallthickness = nan(1,length(SET(no).EpiX(:,tf,slice)));
        for zloop = 1:length(SET(no).EpiX(:,tf,slice))
          tempwallthickness = sqrt(...
            (SET(no).EpiX(zloop,tf,slice)-SET(no).EndoX(:,tf,slice)).^2+ ...
            (SET(no).EpiY(zloop,tf,slice)-SET(no).EndoY(:,tf,slice)).^2);
          wallthickness(zloop) = min(tempwallthickness);
        end
        isnotoutflowtract = find(wallthickness>(1/meanres));
      else
        isnotoutflowtract = 1:length(SET(no).EpiX(:,tf,slice));
      end
      
      %shrink/expand epi
      mx = mean(SET(no).EpiX(:,tf,slice));
      my = mean(SET(no).EpiY(:,tf,slice));
      SET(no).EpiX(isnotoutflowtract,tf,slice) = mx+f*(SET(no).EpiX(isnotoutflowtract,tf,slice)-mx);
      SET(no).EpiY(isnotoutflowtract,tf,slice) = my+f*(SET(no).EpiY(isnotoutflowtract,tf,slice)-my);
    end
  case 'abort'
    beep;
    DATA.CursorX = [];
    DATA.CursorY = [];
    dragepi_Motion('reset');
    set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
    return;
end;

%Reset
DATA.CursorX = [];
DATA.CursorY = [];
dragepi_Motion('reset');

updatemodeldisplay;
drawfunctions('drawallslices');%was drawsliceno but changed in order to update intersectionpoints
updatevolume;
DATA.updateaxestables('t2star');


%------------------
function rotatetemp(alpha) %#ok<DEFNU>
%------------------
%This function should late be replaced with qa tool that rotates objects,
%just as scale and move does. 

global NO SET

if nargin==0
  alpha = 10;
end;

x = SET(NO).EndoX(:,:,SET(NO).CurrentSlice);
y = SET(NO).EndoY(:,:,SET(NO).CurrentSlice);
mx = mean(x(:));
my = mean(y(:));

%Translate to origo
x = x - mx;
y = y - my;

%Define rotation
alpha = alpha/180*pi; %from degrees
sinalpha = sin(alpha);
cosalpha = cos(alpha);

%Rotate
xnew = x*cosalpha+y*sinalpha;
ynew = -x*sinalpha+y*cosalpha;

%Translage back
xnew = xnew + mx;
ynew = ynew + my;

SET(NO).EndoX(:,:,SET(NO).CurrentSlice) = xnew;
SET(NO).EndoY(:,:,SET(NO).CurrentSlice) = ynew;

updatemodeldisplay;

drawfunctions('drawallslices');%was drawsliceno but changed in order to update intersectionpoints

%-----------------------------------
function move_Buttondown(type,panel) %#ok<DEFNU>
%-----------------------------------
%Button down function for moving / translating objects / contours.
global DATA SET NO

if nargin<1
  type = 'image';
end;

if nargin<2
  panel=DATA.CurrentPanel;
end;

switchtopanel(panel);

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

selectiontype = get(DATA.imagefig,'SelectionType');
switch selectiontype
  case 'alt'
    DATA.contextmenu;
  case 'normal'
    %--- Normal  click
    [y,x,slice] = getclickedcoords;
    % If slice has changed, make sure montage/one are in sync
    if (slice>SET(no).ZSize)
      return;
    end
    if ~isempty(SET(no).StartSlice)
      if slice < SET(no).StartSlice || slice > SET(no).EndSlice
        switchtoslice(slice);
      end
    end;
    tools('enableundo')
    
    [DATA.CursorYOfs,DATA.CursorXOfs] = calcfunctions('calcoffset',slice,[],no,panel);
    
    %Called when clicked on an object and move tool.
    switch type
      case 'endo'
        if isempty(SET(no).EndoX)
          myfailed('No LV endocardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).EndoX(:,SET(no).CurrentTimeFrame,slice);
        y = SET(no).EndoY(:,SET(no).CurrentTimeFrame,slice);
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','r');
      case 'epi'
        if isempty(SET(no).EpiX)
          myfailed('No LV epicardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).EpiX(:,SET(no).CurrentTimeFrame,slice);
        y = SET(no).EpiY(:,SET(no).CurrentTimeFrame,slice);
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','g');
      case 'rvendo'
        if isempty(SET(no).RVEndoX)
          myfailed('No RV endocardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,slice);
        y = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,slice);
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','m');
      case 'rvepi'
        if isempty(SET(no).RVEpiX)
          myfailed('No RV epicardium available.',DATA.GUI.Segment);
          return;
        end;
        x = SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,slice);
        y = SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,slice);
        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','c');
      case 'allseg'
        
      case 'roi'
        %Find closest ROI
        if SET(no).RoiN<1
          myfailed('No ROIs available.',DATA.GUI.Segment);
          return;
        end;
        
        mindist = 1e10;
        n = NaN;
        for loop=1:SET(no).RoiN
          if ~isempty(find(SET(no).Roi(loop).Z==slice,1))&& ...
            ~isempty(find(SET(no).Roi(loop).T==SET(no).CurrentTimeFrame,1))
            xr = SET(no).Roi(loop).X(:,SET(no).CurrentTimeFrame);
            yr = SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame);
            thisdist = min(sqrt((xr-x).^2+(yr-y).^2));
            if thisdist<mindist
              n = loop;
              mindist = thisdist;
            end;
          end
        end;

        if isnan(n)
          myfailed('No ROIs available.',DATA.GUI.Segment);
          return;
        end;

        SET(no).RoiCurrent = n;
        %Extract data
        x = SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame);
        y = SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame);

        set(DATA.Handles.cursor(DATA.CurrentPanel),'color','b');
      case 'image'
%        if ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'one','mmodespatial'})
          pan_Buttondown;
%        end
        return;
    end;

    %generellt
    DATA.CursorX = x;
    DATA.CursorY = y;    
    move_Motion('reset');
    set(DATA.imagefig,'WindowButtonUpFcn',...
      sprintf('segment(''move_Buttonup'',''%s'')',type));
    motionfcn = 'segment(''move_Motion'');'; %get(DATA.imagefig,'WindowButtonMotionFcn')];
    set(DATA.imagefig,'WindowButtonMotionFcn',motionfcn);

    set(DATA.Handles.cursor(DATA.CurrentPanel),...
      'xdata',DATA.CursorY+DATA.CursorXOfs,...
      'ydata',DATA.CursorX+DATA.CursorYOfs,...
      'visible','on');
  otherwise
    return; %case 'open', 'extended'
end;
    

%--------------------------
function move_Motion(reset) %#ok<INUSD>
%--------------------------
%Motion function for moving / translating objects / contours.
persistent startx starty startslice 
global DATA SET NO

if (nargin==1)
  startx = [];
  starty = [];
end;

if isempty(startx)
  %Get coordinate 
  [starty,startx,startslice] = getclickedcoords;
else
  %Get coordinate
  [y,x,slice] = getclickedcoords;
  %If different slice then buttonup
  if (slice~=startslice) ...
      || (isempty(DATA.CursorX) && ...
      (any((DATA.CursorX-startx+x)<0.5)||any((DATA.CursorX-startx+x)>SET(NO).YSize))) ...
      || (isempty(DATA.CursorY) && ...
      (any((DATA.CursorY-starty+y)<0.5)||any((DATA.CursorY-starty+y)>SET(NO).XSize)))   
    move_Buttonup('abort');
  end;
  if nargin==0
    if (slice~=startslice) 
      move_Buttonup('abort');
    end;
  end;
  
set(DATA.Handles.cursor(DATA.CurrentPanel),...
    'ydata',DATA.CursorYOfs+DATA.CursorX-startx+x,...
    'xdata',DATA.CursorXOfs+DATA.CursorY-starty+y);
  %set(DATA.Handles.cursor(DATA.CurrentPanel),...
   % 'ydata',DATA.CursorYOfs+DATA.CursorX-startx+x,...
   % 'xdata',DATA.CursorXOfs+DATA.CursorY-starty+y);
end;

%---------------------------
function move_Buttonup(type)
%---------------------------
%Button up function for moving functions.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

DATA.buttonup_Callback;
set(DATA.imagefig,'WindowButtonUpFcn',sprintf('%s(''buttonup_Callback'')','segment'));
set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');

%Extract the data
x = get(DATA.Handles.cursor(DATA.CurrentPanel),'ydata')-DATA.CursorYOfs;
y = get(DATA.Handles.cursor(DATA.CurrentPanel),'xdata')-DATA.CursorXOfs;

%Calculate displacement
dx = mean(x)-mean(DATA.CursorX);
dy = mean(y)-mean(DATA.CursorY);

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = true;
% end
if DATA.ThisFrameOnly %&& ~strcmp(type, 'roi')
  tf = SET(no).CurrentTimeFrame;
else
  tf = 1:SET(no).TSize;
end;

slice = SET(no).StartSlice:SET(no).EndSlice;
pinwarn=false;
switch type
  case 'endo'
      SET(no).EndoX(:,tf,slice) = SET(no).EndoX(:,tf,slice)+dx;
      SET(no).EndoY(:,tf,slice) = SET(no).EndoY(:,tf,slice)+dy;
      if ~isempty(SET(no).EndoPinX)
        [SET(no).EndoPinX,SET(no).EndoPinY,pinwarn]=pinresize(...
          SET(no).EndoPinX,SET(no).EndoPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).EndoInterpX)
        [SET(no).EndoInterpX,SET(no).EndoInterpY,pinwarn]=pinresize(...
          SET(no).EndoInterpX,SET(no).EndoInterpY,tf,slice,[dx dy]);
      end
  case 'epi'
      SET(no).EpiX(:,tf,slice) = SET(no).EpiX(:,tf,slice)+dx;
      SET(no).EpiY(:,tf,slice) = SET(no).EpiY(:,tf,slice)+dy;
      if ~isempty(SET(no).EpiPinX)
        [SET(no).EpiPinX,SET(no).EpiPinY,pinwarn]=pinresize(...
          SET(no).EpiPinX,SET(no).EpiPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).EpiInterpX)
        [SET(no).EpiInterpX,SET(no).EpiInterpY,pinwarn]=pinresize(...
          SET(no).EpiInterpX,SET(no).EpiInterpY,tf,slice,[dx dy]);
      end
  case 'rvendo'
      SET(no).RVEndoX(:,tf,slice) = SET(no).RVEndoX(:,tf,slice)+dx;
      SET(no).RVEndoY(:,tf,slice) = SET(no).RVEndoY(:,tf,slice)+dy;
      if ~isempty(SET(no).RVEndoPinX)
        [SET(no).RVEndoPinX,SET(no).RVEndoPinY,pinwarn]=pinresize(...
          SET(no).RVEndoPinX,SET(no).RVEndoPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).RVEndoInterpX)
        [SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,pinwarn]=pinresize(...
          SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,tf,slice,[dx dy]);
      end
  case 'rvepi'
      SET(no).RVEpiX(:,tf,slice) = SET(no).RVEpiX(:,tf,slice)+dx;
      SET(no).RVEpiY(:,tf,slice) = SET(no).RVEpiY(:,tf,slice)+dy;    
      if ~isempty(SET(no).RVEpiPinX)
        [SET(no).RVEpiPinX,SET(no).RVEpiPinY,pinwarn]=pinresize(...
          SET(no).RVEpiPinX,SET(no).RVEpiPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).RVEpiInterpX)
        [SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,pinwarn]=pinresize(...
          SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,tf,slice,[dx dy]);
      end
  case 'roi'
      SET(no).Roi(SET(no).RoiCurrent).X(:,tf) = SET(no).Roi(SET(no).RoiCurrent).X(:,tf)+dx;
      SET(no).Roi(SET(no).RoiCurrent).Y(:,tf) = SET(no).Roi(SET(no).RoiCurrent).Y(:,tf)+dy;
      roi('roiforceapply'); %ensure that all rois are equal size if this option is checked.
      [m,sd]=calcfunctions('calcroiintensity',no,SET(no).RoiCurrent);
      SET(no).Roi(SET(no).RoiCurrent).Mean = m;
      SET(no).Roi(SET(no).RoiCurrent).StD = sd;
       %update flow result panel
       if ~isempty(DATA.FlowNO) && ismember(DATA.FlowNO,SET(no).Linked) && ~isempty(DATA.FlowROI)
         if SET(no).RoiN <1
           DATA.FlowROI = [];
         else
           calcfunctions('calcflow',no);
         end
       end
       DATA.updateaxestables('area',no);
  case 'abort'
    beep;
    DATA.CursorX = [];
    DATA.CursorY = [];
    set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
    return;
    %myfailed('You can not translate between slices. Aborted.',DATA.GUI.Segment);
end;

if (pinwarn)
  mywarning('Pins moved outside frame deleted.',DATA.GUI.Segment);
end

set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
updatemodeldisplay;
drawfunctions('drawallslices');% was drawsliceno but changed in order to update intersection points
updatevolume;
DATA.updateaxestables('t2star');

%Reset
DATA.CursorX = [];
DATA.CursorY = [];


%-----------------------------------
function moveall_Buttondown(type,panel) %#ok<DEFNU>
%-----------------------------------
%Button down function for moving / translating all objects / contours.
global DATA SET NO

if nargin<1
  type = 'image';
end;

if nargin<2
  panel=DATA.CurrentPanel;
end;

switchtopanel(panel);

% statesandicons=segment('iconson',{'hidelv','hiderv'});
% lvicon=statesandicons{1,2};
% rvicon=statesandicons{2,2};
% 
% %trick that buttons are pressed
% lvicon.isindented=1;
% rvicon.isindented=1;
% 
% %updatevisibility
% viewhiderv_Callback;
% viewhidelv_Callback;

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

selectiontype = get(DATA.imagefig,'SelectionType');
switch selectiontype
  case 'alt'
    DATA.contextmenu;
  case 'normal'
    %--- Normal  click
    [y,x,slice] = getclickedcoords;
    % If slice has changed, make sure montage/one are in sync
    if (slice>SET(no).ZSize)
      return;
    end
    if ~isempty(SET(no).StartSlice)
      if slice < SET(no).StartSlice || slice > SET(no).EndSlice
        switchtoslice(slice);
      end
    end;
    tools('enableundo')
    
    [DATA.CursorYOfs,DATA.CursorXOfs] = calcfunctions('calcoffset',slice,[],no,panel);
    x=[];
    y=[];
    if ~isempty(SET(no).EpiX)
      x = [x;nan;SET(no).EpiX(:,SET(no).CurrentTimeFrame,slice)];
      y = [y;nan;SET(no).EpiY(:,SET(no).CurrentTimeFrame,slice)];
    end
    
    if ~isempty(SET(no).EndoX)
      x = [x;nan;SET(no).EndoX(:,SET(no).CurrentTimeFrame,slice)];
      y = [y;nan;SET(no).EndoY(:,SET(no).CurrentTimeFrame,slice)];
    end
    
    if ~isempty(SET(no).RVEndoX)
      x = [x;nan;SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,slice)];
      y = [y;nan;SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,slice)];
    end
    
    if ~isempty(SET(no).RVEpiX)
      x = [x;nan;SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,slice)];
      y = [y;nan;SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,slice)];
    end
    set(DATA.Handles.cursor(DATA.CurrentPanel),'color','w');
      %Called when clicked on an object and move tool.
%     switch type
%       case 'endo'
%         if isempty(SET(no).EndoX)
%           myfailed('No LV endocardium available.',DATA.GUI.Segment);
%           return;
%         end;
%         x = SET(no).EndoX(:,SET(no).CurrentTimeFrame,slice);
%         y = SET(no).EndoY(:,SET(no).CurrentTimeFrame,slice);
%         set(DATA.Handles.cursor(DATA.CurrentPanel),'color','r');
%       case 'epi'
%         if isempty(SET(no).EpiX)
%           myfailed('No LV epicardium available.',DATA.GUI.Segment);
%           return;
%         end;
%         x = SET(no).EpiX(:,SET(no).CurrentTimeFrame,slice);
%         y = SET(no).EpiY(:,SET(no).CurrentTimeFrame,slice);
%         set(DATA.Handles.cursor(DATA.CurrentPanel),'color','g');
%       case 'rvendo'
%         if isempty(SET(no).RVEndoX)
%           myfailed('No RV endocardium available.',DATA.GUI.Segment);
%           return;
%         end;
%         x = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,slice);
%         y = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,slice);
%         set(DATA.Handles.cursor(DATA.CurrentPanel),'color','m');
%       case 'rvepi'
%         if isempty(SET(no).RVEpiX)
%           myfailed('No RV epicardium available.',DATA.GUI.Segment);
%           return;
%         end;
%         x = SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,slice);
%         y = SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,slice);
%         set(DATA.Handles.cursor(DATA.CurrentPanel),'color','c');
%       case 'roi'
%         %Find closest ROI
%         if SET(no).RoiN<1
%           myfailed('No ROIs available.',DATA.GUI.Segment);
%           return;
%         end;
%         
%         mindist = 1e10;
%         n = NaN;
%         for loop=1:SET(no).RoiN
%           if ~isempty(find(SET(no).Roi(loop).Z==slice,1))&& ...
%             ~isempty(find(SET(no).Roi(loop).T==SET(no).CurrentTimeFrame,1))
%             xr = SET(no).Roi(loop).X(:,SET(no).CurrentTimeFrame);
%             yr = SET(no).Roi(loop).Y(:,SET(no).CurrentTimeFrame);
%             thisdist = min(sqrt((xr-x).^2+(yr-y).^2));
%             if thisdist<mindist
%               n = loop;
%               mindist = thisdist;
%             end;
%           end
%         end;
% 
%         if isnan(n)
%           myfailed('No ROIs available.',DATA.GUI.Segment);
%           return;
%         end;
% 
%         SET(no).RoiCurrent = n;
%         %Extract data
%         x = SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame);
%         y = SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame);
% 
%         set(DATA.Handles.cursor(DATA.CurrentPanel),'color','b');
%       case 'image'
% %        if ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'one','mmodespatial'})
%           pan_Buttondown;
% %        end
%         return;
%     end;

    %generellt
    DATA.CursorX = x;
    DATA.CursorY = y;    
    move_Motion('reset');
    set(DATA.imagefig,'WindowButtonUpFcn',...
      sprintf('segment(''moveall_Buttonup'',''%s'')',type));
    motionfcn = 'segment(''move_Motion'');'; %get(DATA.imagefig,'WindowButtonMotionFcn')];
    set(DATA.imagefig,'WindowButtonMotionFcn',motionfcn);

    set(DATA.Handles.cursor(DATA.CurrentPanel),...
      'xdata',DATA.CursorY+DATA.CursorXOfs,...
      'ydata',DATA.CursorX+DATA.CursorYOfs,...
      'visible','on');
  otherwise
    return; %case 'open', 'extended'
end;
    

%--------------------------
function moveall_Motion(reset) %#ok<INUSD>
%--------------------------
%Motion function for moving / translating all objects / contours.
persistent startx starty startslice 
global DATA SET NO

no = NO;
if (nargin==1)
  startx = [];
  starty = [];
end;

if isempty(startx)
  %Get coordinate 
  [starty,startx,startslice] = getclickedcoords;
else
  %Get coordinate
  [y,x,slice] = getclickedcoords;

  %If different slice then buttonup
  if (slice~=startslice) ...
      || (isempty(DATA.CursorX) && ...
      (any((DATA.CursorX-startx+x)<0.5)||any((DATA.CursorX-startx+x)>SET(NO).YSize))) ...
      || (isempty(DATA.CursorY) && ...
      (any((DATA.CursorY-starty+y)<0.5)||any((DATA.CursorY-starty+y)>SET(NO).XSize)))   
    moveall_Buttonup('abort');
  end;
  if nargin==0
    if (slice~=startslice) 
      moveall_Buttonup('abort');
    end;
  end;
  
 dx =x-startx;% mean(x)-mean(DATA.CursorX);
 dy =y-starty;%mean(y)-mean(DATA.CursorY);
% 
 tools('translatecontours',dx,dy);
  %set(DATA.Handles.cursor(DATA.CurrentPanel),...
  %  'ydata',DATA.CursorYOfs+DATA.CursorX-startx+x,...
  %  'xdata',DATA.CursorXOfs+DATA.CursorY-starty+y);
end;
% x = get(DATA.Handles.cursor(DATA.CurrentPanel),'ydata')-DATA.CursorYOfs;
% y = get(DATA.Handles.cursor(DATA.CurrentPanel),'xdata')-DATA.CursorXOfs;
% % 
%  dx =x-startx;% mean(x)-mean(DATA.CursorX);
%  dy =y-starty;%mean(y)-mean(DATA.CursorY);
% % 
%  tools('translatecontours',dx,dy);

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = true;
% end
% if DATA.ThisFrameOnly
%   tf = SET(no).CurrentTimeFrame;
% else
%   tf = 1:SET(no).TSize;
% end;
% slice = SET(no).StartSlice:SET(no).EndSlice;
% pinwarn=false;
% %Extract the data
% x = get(DATA.Handles.cursor(DATA.CurrentPanel),'ydata')-DATA.CursorYOfs;
% y = get(DATA.Handles.cursor(DATA.CurrentPanel),'xdata')-DATA.CursorXOfs;
% %Calculate displacement
% dx = mean(x)-mean(DATA.CursorX);
% dy = mean(y)-mean(DATA.CursorY);
% disp(mean(x)), disp(mean(DATA.CursorX))
% if ~isempty(SET(no).EndoX) && ~isnan(SET(no).EndoX(1,tf,slice))
%   SET(no).EndoX(:,tf,slice) = SET(no).EndoX(:,tf,slice)+dx;
%   SET(no).EndoY(:,tf,slice) = SET(no).EndoY(:,tf,slice)+dy;
%   if ~isempty(SET(no).EndoPinX)
%     [SET(no).EndoPinX,SET(no).EndoPinY,pinwarn]=pinresize(...
%       SET(no).EndoPinX,SET(no).EndoPinY,tf,slice,[dx dy]);
%   end
%   if ~isempty(SET(no).EndoInterpX)
%     [SET(no).EndoInterpX,SET(no).EndoInterpY,pinwarn]=pinresize(...
%       SET(no).EndoInterpX,SET(no).EndoInterpY,tf,slice,[dx dy]);
%   end
% end
% if ~isempty(SET(no).EpiX) && ~isnan(SET(no).EpiX(1,tf,slice))
%   SET(no).EpiX(:,tf,slice) = SET(no).EpiX(:,tf,slice)+dx;
%   SET(no).EpiY(:,tf,slice) = SET(no).EpiY(:,tf,slice)+dy;
%   if ~isempty(SET(no).EpiPinX)
%     [SET(no).EpiPinX,SET(no).EpiPinY,pinwarn]=pinresize(...
%       SET(no).EpiPinX,SET(no).EpiPinY,tf,slice,[dx dy]);
%   end
%   if ~isempty(SET(no).EpiInterpX)
%     [SET(no).EpiInterpX,SET(no).EpiInterpY,pinwarn]=pinresize(...
%       SET(no).EpiInterpX,SET(no).EpiInterpY,tf,slice,[dx dy]);
%   end
% end
% if ~isempty(SET(no).RVEndoX) && ~isnan(SET(no).RVEndoX(1,tf,slice))
%   SET(no).RVEndoX(:,tf,slice) = SET(no).RVEndoX(:,tf,slice)+dx;
%   SET(no).RVEndoY(:,tf,slice) = SET(no).RVEndoY(:,tf,slice)+dy;
%   if ~isempty(SET(no).RVEndoPinX)
%     [SET(no).RVEndoPinX,SET(no).RVEndoPinY,pinwarn]=pinresize(...
%       SET(no).RVEndoPinX,SET(no).RVEndoPinY,tf,slice,[dx dy]);
%   end
%   if ~isempty(SET(no).RVEndoInterpX)
%     [SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,pinwarn]=pinresize(...
%       SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,tf,slice,[dx dy]);
%   end
% end
% if ~isempty(SET(no).RVEpiX) && ~isnan(SET(no).RVEpiX(1,tf,slice))
%   SET(no).RVEpiX(:,tf,slice) = SET(no).RVEpiX(:,tf,slice)+dx;
%   SET(no).RVEpiY(:,tf,slice) = SET(no).RVEpiY(:,tf,slice)+dy;
%   if ~isempty(SET(no).RVEpiPinX)
%     [SET(no).RVEpiPinX,SET(no).RVEpiPinY,pinwarn]=pinresize(...
%       SET(no).RVEpiPinX,SET(no).RVEpiPinY,tf,slice,[dx dy]);
%   end
%   if ~isempty(SET(no).RVEpiInterpX)
%     [SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,pinwarn]=pinresize(...
%       SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,tf,slice,[dx dy]);
%   end
% end
% if SET(no).RoiN > 0
%   for roin = 1:SET(no).RoiN
%     if ismember(SET(no).Roi(roin).Z,slice)
%       SET(no).Roi(roin).X(:,tf) = SET(no).Roi(roin).X(:,tf)+dx;
%       SET(no).Roi(roin).Y(:,tf) = SET(no).Roi(roin).Y(:,tf)+dy;
%       roi('roiforceapply'); %ensure that all rois are equal size if this option is checked.
%       [m,sd]=calcfunctions('calcroiintensity',no,SET(no).RoiCurrent);
%       SET(no).Roi(roin).Mean = m;
%       SET(no).Roi(roin).StD = sd;
%     end
%   end
% end
% 
% if (pinwarn)
%   mywarning('Pins moved outside frame deleted.',DATA.GUI.Segment);
% end
    

%set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
%updatemodeldisplay;
%drawfunctions('drawallslices');% was drawsliceno but changed in order to update intersection points

%---------------------------
function moveall_Buttonup(type)
%---------------------------
%Button up function for translating all objects / contours
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

% statesandicons=segment('iconson',{'hidelv','hiderv'});
% lvicon=statesandicons{1,2};
% rvicon=statesandicons{2,2};
% 
% %trick that buttons aren't pressed
% lvicon.isindented=0;
% rvicon.isindented=0;
% 
% %updatevisibility
% viewhiderv_Callback;
% viewhidelv_Callback;

DATA.buttonup_Callback;
set(DATA.imagefig,'WindowButtonUpFcn',sprintf('%s(''buttonup_Callback'')','segment'));
set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');

%Extract the data
x = get(DATA.Handles.cursor(DATA.CurrentPanel),'ydata')-DATA.CursorYOfs;
y = get(DATA.Handles.cursor(DATA.CurrentPanel),'xdata')-DATA.CursorXOfs;

%Calculate displacement
dx = nanmean(x)-nanmean(DATA.CursorX);
dy = nanmean(y)-nanmean(DATA.CursorY);

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = true;
% end
if DATA.ThisFrameOnly
  tf = SET(no).CurrentTimeFrame;
else
  tf = 1:SET(no).TSize;
end;

slice = SET(no).StartSlice:SET(no).EndSlice;
pinwarn=false;

switch type
  case 'abort'
    beep;
    DATA.CursorX = [];
    DATA.CursorY = [];
    set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
    return;
    %myfailed('You can not translate between slices. Aborted.',DATA.GUI.Segment);
  otherwise
    if ~isempty(SET(no).EndoX) %&& all(~isnan(SET(no).EndoX(1,tf,slice)))
      SET(no).EndoX(:,tf,slice) = SET(no).EndoX(:,tf,slice)+dx;
      SET(no).EndoY(:,tf,slice) = SET(no).EndoY(:,tf,slice)+dy;
      if ~isempty(SET(no).EndoPinX)
        [SET(no).EndoPinX,SET(no).EndoPinY,pinwarn]=pinresize(...
          SET(no).EndoPinX,SET(no).EndoPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).EndoInterpX)
        [SET(no).EndoInterpX,SET(no).EndoInterpY,pinwarn]=pinresize(...
          SET(no).EndoInterpX,SET(no).EndoInterpY,tf,slice,[dx dy]);
      end
    end
    if ~isempty(SET(no).EpiX) %&& all(~isnan(SET(no).EpiX(1,tf,slice)))
      SET(no).EpiX(:,tf,slice) = SET(no).EpiX(:,tf,slice)+dx;
      SET(no).EpiY(:,tf,slice) = SET(no).EpiY(:,tf,slice)+dy;
      if ~isempty(SET(no).EpiPinX)
        [SET(no).EpiPinX,SET(no).EpiPinY,pinwarn]=pinresize(...
          SET(no).EpiPinX,SET(no).EpiPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).EpiInterpX)
        [SET(no).EpiInterpX,SET(no).EpiInterpY,pinwarn]=pinresize(...
          SET(no).EpiInterpX,SET(no).EpiInterpY,tf,slice,[dx dy]);
      end
    end
    if ~isempty(SET(no).RVEndoX) %&& all(~isnan(SET(no).RVEndoX(1,tf,slice)))
      SET(no).RVEndoX(:,tf,slice) = SET(no).RVEndoX(:,tf,slice)+dx;
      SET(no).RVEndoY(:,tf,slice) = SET(no).RVEndoY(:,tf,slice)+dy;
      if ~isempty(SET(no).RVEndoPinX)
        [SET(no).RVEndoPinX,SET(no).RVEndoPinY,pinwarn]=pinresize(...
          SET(no).RVEndoPinX,SET(no).RVEndoPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).RVEndoInterpX)
        [SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,pinwarn]=pinresize(...
          SET(no).RVEndoInterpX,SET(no).RVEndoInterpY,tf,slice,[dx dy]);
      end
    end
    if ~isempty(SET(no).RVEpiX) %&& all(~isnan(SET(no).RVEpiX(1,tf,slice)))
      SET(no).RVEpiX(:,tf,slice) = SET(no).RVEpiX(:,tf,slice)+dx;
      SET(no).RVEpiY(:,tf,slice) = SET(no).RVEpiY(:,tf,slice)+dy;
      if ~isempty(SET(no).RVEpiPinX)
        [SET(no).RVEpiPinX,SET(no).RVEpiPinY,pinwarn]=pinresize(...
          SET(no).RVEpiPinX,SET(no).RVEpiPinY,tf,slice,[dx dy]);
      end
      if ~isempty(SET(no).RVEpiInterpX)
        [SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,pinwarn]=pinresize(...
          SET(no).RVEpiInterpX,SET(no).RVEpiInterpY,tf,slice,[dx dy]);
      end
    end
    if SET(no).RoiN > 0
      for roin = 1:SET(no).RoiN
        if ismember(SET(no).Roi(roin).Z,slice)
          SET(no).Roi(roin).X(:,tf) = SET(no).Roi(roin).X(:,tf)+dx;
          SET(no).Roi(roin).Y(:,tf) = SET(no).Roi(roin).Y(:,tf)+dy;
          roi('roiforceapply'); %ensure that all rois are equal size if this option is checked.
          [m,sd]=calcfunctions('calcroiintensity',no,SET(no).RoiCurrent);
          SET(no).Roi(roin).Mean = m;
          SET(no).Roi(roin).StD = sd;
        end
      end
      if ~isempty(DATA.FlowNO) && ismember(DATA.FlowNO,SET(no).Linked) && ~isempty(DATA.FlowROI)
        if SET(no).RoiN <1
          DATA.FlowROI = [];
        else
          calcfunctions('calcflow',no);
        end
      end
    end
end;

if (pinwarn)
  mywarning('Pins moved outside frame deleted.',DATA.GUI.Segment);
end

set(DATA.Handles.cursor(DATA.CurrentPanel),'visible','off');
updatemodeldisplay;
drawfunctions('drawallslices');% was drawsliceno but changed in order to update intersection points
updatevolume;
DATA.updateaxestables('t2star');

%Reset
DATA.CursorX = [];
DATA.CursorY = [];


%-----------------------------------
function putpin_Callback(panel,type) %#ok<DEFNU>
%-----------------------------------
%Callback to put pins.
global DATA 

switch type
  case 'endo'
    isendo = true;
    dolv = true;
  case 'epi'
    isendo = false;
    dolv = true;    
  case 'rvendo'
    isendo = true;    
    dolv = false;    
end;

if dolv
  %--- LV
  if isendo
    highlighttool(DATA.Tools.putendopin);
    set(DATA.Handles.imagehandle(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''endo'')',mfilename,panel));
    set(DATA.Handles.endocontour(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''endo'')',mfilename,panel));
    set(DATA.Handles.imagehandle(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''endo'')',mfilename,panel));
  else
    highlighttool(DATA.Tools.putepipin);
    set(DATA.Handles.imagehandle(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''epi'')',mfilename,panel));
    set(DATA.Handles.endocontour(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''epi'')',mfilename,panel));
    set(DATA.Handles.imagehandle(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''epi'')',mfilename,panel));
  end;
else
  %--- RV
  if isendo
    highlighttool(DATA.Tools.putrvendopin);
    set(DATA.Handles.imagehandle(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''rvendo'')',mfilename,panel));
    set(DATA.Handles.rvendocontour(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''rvendo'')',mfilename,panel));
    set(DATA.Handles.imagehandle(panel),'ButtonDownFcn',...
      sprintf('%s(''putpin_Buttondown'',%d,''rvendo'')',mfilename,panel));
  end;
  
end;

%-------------------------------------
function putpin_Buttondown(panel,type) %#ok<DEFNU>
%-------------------------------------
%Button down function for pins.
global DATA 

switchtopanel(panel);

switch get(DATA.imagefig,'SelectionType')
  case 'normal'
    %tools('enableundo') %called in doputpin_Callback();
    doputpin_Callback(type);
  case 'alt'
    %show menu
    DATA.contextmenu;
end;

%--------------------
function pin_Buttonup %#ok<DEFNU>
%--------------------
%Button up function for pins.
global DATA 

checkconsistency;
set(DATA.imagefig,'WindowButtonMotionFcn','');
updatemodeldisplay;
updatevolume;

%------------------------
function pin_Motion(type) %#ok<DEFNU>
%------------------------
%Motion function for pins.
global DATA SET NO

[x,y,slice] = getclickedcoords;

switch type
  case 'endo'
    temp = SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice};
    temp(DATA.Pin) = y;
    SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice} = temp;
    temp = SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice};
    temp(DATA.Pin) = x;
    SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice} = temp;
  case 'epi'    
    temp = SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice};
    temp(DATA.Pin) = y;
    SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice} = temp;
    temp = SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice};
    temp(DATA.Pin) = x;
    SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice} = temp;
  case 'rvendo'
    temp = SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice};
    temp(DATA.Pin) = y;
    SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice} = temp;
    temp = SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice};
    temp(DATA.Pin) = x;
    SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice} = temp;    
end

updatemodeldisplay;
drawfunctions('drawsliceno');

%------------------------------
function [PinX,PinY,pinwarn]=pinresize(PinX,PinY,tf,slice,mean,f) 
%------------------------------
%Helper function to move pins when resizing contours.
global SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

mx=mean(1);
my=mean(2);
if (nargin<6)
  f=1;
  dx=0;
  dy=0;
else
  dx=mx;
  dy=my;
end
pinwarn=false;

for tloop=tf
  for zloop=slice
    if ~isempty(PinX{tloop,zloop})
      PinX{tloop,zloop} = mx+f*(PinX{tloop,zloop}-dx);
      PinY{tloop,zloop} = my+f*(PinY{tloop,zloop}-dy);
      ind= (PinX{tloop,zloop}>0.5)...
        &(PinY{tloop,zloop}>0.5)...
        &(PinX{tloop,zloop}<SET(no).XSize+0.5)...
        &(PinY{tloop,zloop}<SET(no).YSize+0.5);
      if any(~ind(:))
        PinX{tloop,zloop} = PinX{tloop,zloop}(ind);
        PinY{tloop,zloop} = PinY{tloop,zloop}(ind);
        pinwarn=true;
      end
    end
  end
end

%------------------------------
function timebar_Motion %#ok<DEFNU>
%------------------------------
%Update current timeframe when user clicked in time bar graph.
global DATA SET NO

[x] = mygetcurrentpoint(DATA.Handles.timebaraxes);
x = x/1000;

%Find timeframe
SET(NO).CurrentTimeFrame = 1+round(x/SET(NO).TIncr);
SET(NO).CurrentTimeFrame = min(max(1,SET(NO).CurrentTimeFrame),SET(NO).TSize);

%update time bar line
t = SET(NO).TimeVector*1000;
set(DATA.Handles.timebar,'xdata',[t(SET(NO).CurrentTimeFrame) t(SET(NO).CurrentTimeFrame)],...
  'ydata',get(DATA.Handles.timebaraxes,'ylim'));

%Update
% profile on;
% drawfunctions('drawsliceno');
% profile report;

nos = SET(NO).Linked;
%updatenopanels(NO) is performed in drawsliceno.

%additional speed enhancement is that stateandicon is supplied from
%buttondown. However not sure if possible.
stateandicon=segment('iconson',{'hidescar','hidemar','hideall','play'});
for nloop=1:length(nos)
  drawfunctions('updatenopanels',nos(nloop),stateandicon);
end;

drawnow;

%----------------------------------
function timebar_Buttondown %#ok<DEFNU>
%----------------------------------
%Button down function for the time bar in the time bar graph
global DATA

if DATA.Interactionlock
  return;
end;

%stateandicon=segment('iconson',{'hidescar','hidemar','hideall','play'});
%set(DATA.fig,'WindowButtonMotionFcn',@ segment_main.timebar_Motion(stateandicon));
set(DATA.fig,'WindowButtonMotionFcn',...
  sprintf('%s(''timebar_Motion'')',mfilename));
set(DATA.fig,'WindowButtonUpFcn',...
  sprintf('%s(''buttonup_Callback'')','segment'));

%------------------------------
function timebarlv_Motion %#ok<DEFNU>
%------------------------------
%Update current timeframe when user clicked in volume graph.
global DATA SET NO

[x] = mygetcurrentpoint(DATA.Handles.volumeaxes);
x = x/1000;

%Find timeframe
SET(NO).CurrentTimeFrame = 1+round(x/SET(NO).TIncr);
SET(NO).CurrentTimeFrame = min(max(1,SET(NO).CurrentTimeFrame),SET(NO).TSize);

%update time bar line
t = SET(NO).TimeVector*1000;
set(DATA.Handles.timebarlv,'xdata',[t(SET(NO).CurrentTimeFrame) t(SET(NO).CurrentTimeFrame)],...
  'ydata',get(DATA.Handles.volumeaxes,'ylim'));
%Update
drawfunctions('drawsliceno');

%----------------------------------
function timebarlv_Buttondown %#ok<DEFNU>
%----------------------------------
%Button down function for the time bar in the volume graph.
global DATA NO

if NO == DATA.LVNO
  if DATA.Interactionlock
    return;
  end;
  
  set(DATA.fig,'WindowButtonMotionFcn',...
    sprintf('%s(''timebarlv_Motion'')',mfilename));
  set(DATA.fig,'WindowButtonUpFcn',...
    sprintf('%s(''buttonup_Callback'')','segment'));
end

%------------------------------
function timebarflow_Motion %#ok<DEFNU>
%------------------------------
%Update current timeframe when user clicked in volume graph.
global DATA SET NO

[x] = mygetcurrentpoint(DATA.Handles.flowaxes);
x = x/1000;

%Find timeframe
SET(NO).CurrentTimeFrame = 1+round(x/SET(NO).TIncr);
SET(NO).CurrentTimeFrame = min(max(1,SET(NO).CurrentTimeFrame),SET(NO).TSize);
%update time bar line
t = SET(NO).TimeVector*1000;
set(DATA.Handles.timebarflow,'xdata',[t(SET(NO).CurrentTimeFrame) t(SET(NO).CurrentTimeFrame)],...
  'ydata',get(DATA.Handles.flowaxes,'ylim'));
%Update
nos = SET(NO).Linked;
for nloop=1:length(nos)
  drawfunctions('updatenopanels',nos(nloop));
end;
% drawfunctions('drawsliceno');

%----------------------------------
function timebarflow_Buttondown %#ok<DEFNU>
%----------------------------------
%Button down function for the time bar in the volume graph.
global DATA SET NO

if NO == DATA.FlowNO || ismember(NO,SET(DATA.FlowNO).Linked)
  if DATA.Interactionlock
    return;
  end;
  
  set(DATA.fig,'WindowButtonMotionFcn',...
    sprintf('%s(''timebarflow_Motion'')',mfilename));
  set(DATA.fig,'WindowButtonUpFcn',...
    sprintf('%s(''buttonup_Callback'')','segment'));
end

%--------------------------------
function fasterframerate_Callback %#ok<DEFNU>
%--------------------------------
%Makes movie play faster
global SET%DATA SET

%set(DATA.Handles.fasterframerateicon,'state','off');
ind = find(cat(1,SET(:).TSize)>1);
if isempty(ind)
  ind = 1;
end;

for no=ind'
  SET(no).BeatTime = SET(no).BeatTime*0.9;
end

%DATA.StartTime = now;
%if ~DATA.Run 
 % playmovie_Callback;
%end;

%--------------------------------
function slowerframerate_Callback %#ok<DEFNU>
%--------------------------------
%Makes movie play slower
global SET%DATA SET

%set(DATA.Handles.slowerframerateicon,'state','off');
%If running make sure always update the first one.
ind = find(cat(1,SET(:).TSize)>1);
if isempty(ind)
  ind = 1;
end;

for no=ind'
  SET(no).BeatTime = SET(no).BeatTime/0.9;
end

%DATA.StartTime = now;
%if ~DATA.Run
%  playmovie_Callback;
%end;

%--------------
function ctrlc %#ok<DEFNU>
%--------------
%function to handle ctrl-c keypress, which is disabled
global DATA
mywarning('Ctrl-C is disabled',DATA.GUI.Segment);

%--------------------------------------------------------
function [varargout] = cell2clipboard(outdata,writetofile) %#ok<DEFNU>
%--------------------------------------------------------
%Converts a cell to a string that is output to clipboard.
%If more than 8000 cells are written then an .xls file is
%written instead. Note that this used active-X on Windows and requires
%Excel to be installed on the computer.

global DATA

if nargout>0
  varargout = cell(1,nargout);
end;

if nargin<2
  writetofile = false;
end;

if not(writetofile) && ( (numel(outdata)<8000)||(nargout>0) )
  stri = [];
  for rloop=1:size(outdata,1)
    for cloop=1:size(outdata,2)
      temp = outdata{rloop,cloop};
      switch class(temp)
        case 'logical'
          if temp
            stri = [stri sprintf('true\t')]; %#ok<AGROW>
          else
            stri = [stri sprintf('false\t')]; %#ok<AGROW>
          end;
        case {'uint8','int8','uint32'}
          stri = [stri sprintf('%d\t',temp)]; %#ok<AGROW>
        case {'double','single'}
          if length(temp)>1
            myfailed('Vector not allowed in one cell.',DATA.GUI.Segment);
            return;
          end;
          if isempty(temp)||isnan(temp)
            stri = [stri sprintf('\t')]; %#ok<AGROW>
          else
            if isempty(temp)
              stri = [stri sprintf('\t')]; %#ok<AGROW>
            else
              stri = [stri sprintf('%f\t',temp)]; %#ok<AGROW>
            end;
          end;
        case 'char'
          %If a char has ascii 0 its a null char this might cause trouble
          %for export.
          nullremove=double(temp)==0;
          if any(nullremove)
            temp(nullremove)=[];
          end
          
          if isempty(temp)
            stri = [stri sprintf('\t')]; %#ok<AGROW>
          else
            stri = [stri sprintf('%s\t',strtrim(temp))]; %#ok<AGROW>
          end;
        otherwise
          myfailed('Unknown object type when exporting.',DATA.GUI.Segment);
          disp(class(temp))
          return;
      end;
    end;
    stri = [stri sprintf('\n')]; %#ok<AGROW>
  end;
  if nargout==0
    clipboard('copy',stri);
    mymsgbox('Data copied to clipboard.','Done!',DATA.GUI.Segment);
  else
    varargout{1} = stri;
  end;
else
  [filename, pathname] = myuiputfile('*.xls', 'Save as Excel file.');
  if isequal(filename,0) || isequal(pathname,0)
    disp('Aborted');
    return;
  else
    [sucess,message] = xlswrite(fullfile(pathname,filename),outdata);
    if sucess
      mymsgbox('File successfully written.','Done!',DATA.GUI.Segment);
    else
      myfailed(dprintf('Error:%s',message.message),DATA.GUI.Segment);
      return;
    end
  end;
end;
flushlog;

%--------------------------
function mccheckconsistency %#ok<DEFNU>
%--------------------------
%Check consistency, to prevent earlier
%manual segmentations that have problems with
%direction left/right
global DATA SET NO


timeframes=1:SET(NO).TSize;
slice=1:SET(NO).ZSize;

for timeframe = timeframes
  for zloop=slice;

    %Epi
    if ~isempty(SET(NO).EpiX)
      if not(isnan(SET(NO).EpiX(1,timeframe,zloop)))

        xr = SET(NO).EpiY(:,timeframe,zloop);
        yr = SET(NO).EpiX(:,timeframe,zloop);
        mx = SET(NO).CenterY;
        my = SET(NO).CenterX;
        
        %Find point closest to cursor
        [~,inda] = min(abs(complex(mx-xr,my-yr)));
        
        SET(NO).EpiX(1:(DATA.NumPoints-inda+1),timeframe,zloop) = yr(inda:end);
        SET(NO).EpiY(1:(DATA.NumPoints-inda+1),timeframe,zloop) = xr(inda:end);
        SET(NO).EpiX((DATA.NumPoints+1-inda):end,timeframe,zloop) = yr(1:inda);
        SET(NO).EpiY((DATA.NumPoints+1-inda):end,timeframe,zloop) = xr(1:inda);
      end;
    end;
  end;
end;

%------------------------------------------
function [x,y,name] = askcontour(queststri) %#ok<DEFNU>
%------------------------------------------
%Show menu so that user can indicate what contour to use.
%Used by levelset to import contours, and by export function
%to export contours as ascii file.
%
%returns contour in x, and y, and a name of the contour.
global DATA SET NO

if nargin==0
  queststri = 'Choose which contour to import';
end;

x = [];
y = [];
name = '';

m=mymenu(queststri,...
  {'LV Endocardium (red)',...
  'LV Epicardium (green)',...
  'RV Endocardium (magenta)',...
  'RV Epicardium (cyan)'},DATA.GUI.Segment);

switch m
  case 1
    %--- Import endocardium
    if isempty(SET(NO).EndoX)
      myfailed('No LV endocardium available.',DATA.GUI.Segment);
      return;
    end;
    y=SET(NO).EndoY;
    x=SET(NO).EndoX;
    name='endocardium';
  case 2
    %--- Import epicardium
    if isempty(SET(NO).EpiX)
      myfailed('No LV epicardium available.',DATA.GUI.Segment);
      return;
    end;
    y=SET(NO).EpiY;
    x=SET(NO).EpiX;
    name='epicardium';
  case 3
    %--- Import RVendocardium
    if isempty(SET(NO).RVEndoX)
      myfailed('No RV endocardium available.',DATA.GUI.Segment);
      return;
    end;
    y=SET(NO).RVEndoY;
    x=SET(NO).RVEndoX;
    name='RVendocardium';
  case 4
    %--- Import RVepicardium
    if isempty(SET(NO).RVEpiX)
      myfailed('No RV epicardium available.',DATA.GUI.Segment);
      return;
    end;
    y=SET(NO).RVEpiY;
    x=SET(NO).RVEpiX;
    name='RVepicardium';
  otherwise
    myfailed('User aborted.',DATA.GUI.Segment);
end

%----------------------------
function evalcommand_Callback %#ok<DEFNU>
%----------------------------
%Function to evaluate matlab commands from compiled version.
global DATA SET NO %#ok<NUSED>

eval_str=inputdlg('command:','Execute matlab command');
if ~isempty(eval_str)
  try
    eval(cell2mat(eval_str));
  catch me
    mydispexception(me);
  end
end

segmenthelp('openthislogfile_Callback');

%---------------------------------
function changewheel_Callback(h,e) %#ok<DEFNU,INUSL>
%---------------------------------
%scrollwheel with modifer.
%Tab are not included but you can if you like
global DATA
if not(DATA.DataLoaded)
  return;
end

wheelup=e.VerticalScrollCount>0;

modifier=char(get(gcbf,'currentmodifier'));

switch size(modifier,1)
  case 2
    mod_str=[deblank(modifier(1,:)) '-' deblank(modifier(2,:))];
  case 3
    mod_str=[deblank(modifier(1,:)) '-' ...
             deblank(modifier(2,:)) '-' ...
             deblank(modifier(3,:))];
  otherwise
    mod_str=modifier;
end

switch mod_str
  case 'shift'  
    if wheelup,
      previousallframe_Callback;
    else
     nextallframe_Callback;
    end
  case 'alt'
    if wheelup,
      zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1.2);
    else
      zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1/1.2);
    end
    drawfunctions('viewupdatetextposition');
    drawfunctions('viewupdateannotext');
  case 'control'
    if wheelup,
      thumbnailslider_Callback('up');
    else
      thumbnailslider_Callback('down');      
    end
  case 'shift-alt'
  case 'control-alt'
  case 'shift-control'
    if wheelup,
      contrast_Callback('wheelUpBrightness',DATA.CurrentPanel);
    else
      contrast_Callback('wheelDownBrightness',DATA.CurrentPanel); 
    end
  case 'shift-control-alt'
  otherwise %no modifier
    if e.VerticalScrollCount<0,
      if DATA.Synchronize
        movealltowardsbase_Callback;
      else
        movetowardsbase_Callback;
      end
    else
      if DATA.Synchronize
        movealltowardsapex_Callback;
      else
        movetowardsapex_Callback;
      end
    end  
end

%---------------------------------
function recursekeypressfcn(h,fcn)
%---------------------------------
%Helper function to create callbacks to keypressed function.
if nargin<2
  fcn = @(x,y)segment('keypressed',x,y);
end;

%Start with the current handle.
try
  set(h,'keypressfcn',fcn);
catch %#ok<CTCH>
end;

children = get(h,'children');
for loop=1:length(children)
  recursekeypressfcn(children(loop),fcn);
end;

%------------------------------------------
function thumbnailslider_Callback(varargin)
%------------------------------------------
%Slider callback for thumbnail slider.
global DATA

temp=DATA.VisibleThumbnails;
updateslider(varargin{:})
recalculatepreview=isempty(DATA.DATASETPREVIEW) | length(DATA.VisibleThumbnails)~=length(temp);
sliderupdated=1;
drawfunctions('drawthumbnails',recalculatepreview,sliderupdated);

%------------------------------
function updateslider(whattodo)
%------------------------------
%Update when user changes in thumbnail slider

global DATA SET

if isempty(DATA.Handles.thumbnailslider) || (~DATA.DataLoaded) || DATA.Silent
  DATA.VisibleThumbnails=1:length(SET);
  return;
end;

%Get range
slidermin=1;
slidermax=min(max(length(SET)-DATA.Pref.NumberVisibleThumbnails+1,1),length(SET));

%get slider value
slidervalue=slidermax-round(mygetvalue(DATA.Handles.thumbnailslider))+1;
slidervalue=min(max(slidervalue,slidermin),slidermax);

if nargin==0
  %No argument simply user changed slider
  if isempty(DATA.VisibleThumbnails)
    DATA.VisibleThumbnails=1:min(DATA.Pref.NumberVisibleThumbnails,length(SET));
  end;

  firstthumbnail=slidervalue;
  lastthumbnail=max(min(firstthumbnail+DATA.Pref.NumberVisibleThumbnails-1,length(SET)),1);
else
  switch whattodo
    case 'down'
      firstthumbnail=max(DATA.VisibleThumbnails(1)-1,1);
      lastthumbnail=min(firstthumbnail+DATA.Pref.NumberVisibleThumbnails-1,length(SET));
      slidervalue = firstthumbnail;
    case 'up'
      lastthumbnail=min(DATA.VisibleThumbnails(end)+1,length(SET));
      firstthumbnail=max(lastthumbnail-DATA.Pref.NumberVisibleThumbnails+1,1);
      slidervalue = firstthumbnail;      
  end;
end;

%Do the update
if slidermin==slidermax
  set(DATA.Handles.thumbnailslider,...
  'min', 0,...
  'max',1,...
  'value',1,...
  'visible','off',...
  'enable','off');
else
  %sliderstep=[1/(length(SET)-DATA.Pref.NumberVisibleThumbnails+1),(slidermax-slidermin+1)];%[0.25/(slidermax-slidermin) 0.5/(slidermax-slidermin)];
  sliderstep = [1/(length(SET)-DATA.Pref.NumberVisibleThumbnails),0.1];
  set(DATA.Handles.thumbnailslider,...
    'min',slidermin,...
    'max',slidermax,...
    'sliderstep',sliderstep,...
    'value',slidermax-slidervalue+1,...
    'visible','on',...
    'enable','on');
  if length(SET) <DATA.Pref.NumberVisibleThumbnails
      set(DATA.Handles.thumbnailslider,'visible','off','enable','off');
  end    
end

DATA.VisibleThumbnails=firstthumbnail:lastthumbnail;

%---------------------
function saveguipositiontodisk %#ok<DEFNU>
%---------------------
%Save guipositions to disk
global DATA

pathname = getpreferencespath;
GUIPositions = DATA.GUIPositions; %#ok<NASGU> %Saved to file
try
  save([pathname filesep '.segment_guipositions.mat'],'GUIPositions', DATA.Pref.SaveVersion);
catch %#ok<CTCH>
  myfailed('Could not save GUI positions. Write permission? Disk full?');
  return;
end;

disp('GUI Positions saved.')

%--------------------------------
function numericversion=getnumericversion
%-----------------------------------------
%function called to get versionnumber after R. For example in '1.9 R4040' this function returns 4040
global DATA

stringversion=DATA.ProgramVersion;
rpos=find(upper(stringversion)=='R');
numericversion=str2double(stringversion(rpos+1:end));

%--------------------------
function resetguipositions %#ok<DEFNU>
%--------------------------
%Resets GUI positions of Segment.
global DATA

%GUI positions default values
GUIPositions(1).FileName='segment.fig';
GUIPositions(1).Position=[0.05 0.05 0.9 0.85];%normalized position

DATA.GUIPositions=GUIPositions;

guinames=fieldnames(DATA.GUI);
for loop=1:length(guinames)
  gui=getfield(DATA.GUI,guinames{loop}); %#ok<GFLD>
  if ~isempty(gui)
    try
      setguiposition(gui);
    catch me
      disp(sprintf('Could not reset GUI %s',guinames{loop})); %#ok<DSPS>
      mydispexception(me);
    end
  end;
end;

%---------------------------------------------------
function corrupted=checkcorrupteddataforautomaticsave(setstruct) %#ok<DEFNU>
%---------------------------------------------------
%This function checks if the data (SET) is corrupted due to corrupted
%loading when loading files which has been saved with older saveversion
% see ticket 502 in wush for more details on the bug
global SET

corrupted=false;
if nargin==0%check SET(:).IM
  for loop=1:length(SET)
    if all(SET(loop).IM(:)==0)
      corrupted=true;
    end
  end
else %check setstruct(:).IM
  for loop=1:length(setstruct)
    if all(setstruct(loop).IM(:)==0)
      corrupted=true;
    end
  end
end

%-------------------------------------------------------
function slicestoinclude = getmontagesegmentedslices(no)
%-------------------------------------------------------
%Get slices to include in montage segmented view.
global SET

slicestoinclude = find(findfunctions('findslicewithendo',no))';
if min(slicestoinclude) > 1
  slicestoinclude = [min(slicestoinclude)-1 slicestoinclude];
end
if max(slicestoinclude < SET(no).ZSize)
  slicestoinclude = [slicestoinclude max(slicestoinclude)+1];
end