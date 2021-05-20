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
  
  %Check if os is supported
  arch = mexext();  
  switch arch
    case {'mexglx','mexmaci','mexmaci64'}
      myfailed('Your platform is not supported. Supported platforms are Windows and Linux 64 bit.');
      return;
  end
  
  if not(isdeployed())
    %source code version, check matlab version
    try
      matlabversion = ver;
      [toolboxindex] = find(ismember({matlabversion.Name},'MATLAB'));
      if isempty(toolboxindex), toolboxindex=1; end
      if not(strcmp(matlabversion(toolboxindex).Release,'(R2019a)'))
        myfailed(sprintf('Recommended Matlab version for Segment is R2019a.\nPlease use this version for best performance of Segment.\nVersions prior to Matlab 2015a will NOT work for this version of Segment.'));
      end
    catch
    end
  end
  
  programversion = changelog;
  fig = initializesegment(programversion); %Program version number
  try
    compilerpragmas;
  catch
  end
  
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
    disp('Already running.');
    if ~isempty(DATA.fig)
      figure(DATA.fig);
      fig = DATA.fig;
      return;
    else
      fig = [];
    end
  catch me
    disp('Program not aborted properly last time, all data will be lost and program restarted.');
    mydispexception(me);
  end
end

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
    end
  end
  
  switch arch
    case {'mac','sol64','glnx86'}
      myfailed('Platform is currently not supported. Mex-files are missing.');
      return;
    otherwise
      disp(['Running from Matlab on platform ' arch '.']);
  end
           
  %Do nothing
end

%Make sure fresh start
DATA = []; %#ok<NASGU>
SET = []; %#ok<NASGU>

%This is where we create object
DATA = segmentgui(programversion);%Load standard Segment GUI
DATA.GUISettings.ShowColorbar = false;  %until the colorbar is impelemented

% Register where segment.m is located and if we are running from source
if isdeployed()
  [status, result] = system('set PATH');
  pathname = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
  DATA.SegmentFolder=pathname;
  cd(pathname);
else
  [segmentfolder_, ~, ~] = fileparts(which('segment'));
  DATA.SegmentFolder = segmentfolder_;
  clear segmentfolder_
end

DATA.init;

try
  fig = DATA.fig;
catch
  disp('Software initialization aborted')
  fig = [];
  return
end

checkpath(DATA.SegmentFolder); %ensure running on correct path


%---------------------
function resetpreview %#ok<DEFNU>
%---------------------
%Reset preview structure in DATA.Preview

global DATA;

DATA.Preview = genemptypreview(DATA.Pref.datapath);

%---------------------------------
function preview = genemptypreview(datapath)
%---------------------------------
%Generate an empty preview struct
global DATA

roisize = 150;

try %inside try catch if DATA is not set.
  if isempty(DATA.Preview) || isempty(DATA.Preview.ROISize)
    roisize = 'full'; %150;
  else
    roisize = DATA.Preview.ROISize;
  end
catch
end

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
  end
else
  %Use filenames
  f = dir([DATA.SegmentFolder filesep 'plugin_*.m']);
  pluginfiles = cell(1,length(f));
  for loop=1:length(f)
    pluginfiles{loop} = f(loop).name;
  end
end

%Initialize menus
for loop=1:length(pluginfiles) 
  stri = pluginfiles{loop};
  stri = stri(1:(end-2)); %Remove .m
  handle = uimenu(DATA.Handles.pluginmenu,...
    'Label','temp',...
    'Callback','');
  namestri = feval(stri,'getname',handle);
  set(handle,'Label',namestri);
end

%--------------------------------
function singleframemode_Callback %#ok<DEFNU>
%--------------------------------
%define if single (one) or all frames mode
global DATA
DATA.ThisFrameOnly = not(DATA.Handles.configiconholder.findindented('selectoneall'));

%--------------------------
function cinetool_Callback %#ok<DEFNU>
%--------------------------
%Starts the cinetool that allows simultanues segmentation at the
%same time as it plays.

global DATA SET NO

if SET(NO).TSize==1
  try %icon might be there and might not
  stateandicon=viewfunctions('iconson','cineplay');
  stateandicon{2}.isindented=0;
  stateandicon{2}.cdataDisplay=stateandicon{2}.cdata;
  DATA.Handles.configiconholder.render
  catch
  end
  myfailed('Need a timeresolved image stack.')
  return;
end

if isa(DATA.CineTimer,'timer')
  cinewindow('update','kill');
  
  try %icon might be there and might not
  stateandicon=viewfunctions('iconson','cineplay');
  stateandicon{2}.isindented=0;
  stateandicon{2}.cdataDisplay=stateandicon{2}.cdata;
  DATA.Handles.configiconholder.render
  catch 
  end
else
  try %icon might be there and might not
  stateandicon=viewfunctions('iconson','cineplay');
  stateandicon{2}.isindented=1;
  stateandicon{2}.cdataDisplay=stateandicon{2}.cdataIndent;
  DATA.Handles.configiconholder.render
  catch
  end
  cinewindow;
end

%----------------------------
function addtopanels(no,mode) %#ok<DEFNU>
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
  end
end

%--- if space exists then add it
temp = find(DATA.ViewPanels==0);
if ~isempty(temp)
  DATA.ViewPanels(temp(1)) = no;
  DATA.ViewPanelsType{temp(1)} = mode;
  DATA.ViewIM{temp(1)} = [];
  return;
end

%--- if space does not exist then add
DATA.ViewPanels = [DATA.ViewPanels no];
DATA.ViewPanelsType = cat(2,DATA.ViewPanelsType,{mode});
[rows,cols] = calcfunctions('calcrowscols',no,SET(no).ZSize);
DATA.ViewPanelsMatrix = [DATA.ViewPanelsMatrix {[rows cols]}];
DATA.ViewIM{length(DATA.ViewPanels)} = [];

%---------------------------
function mainresize_Callback %#ok<DEFNU>
%---------------------------
%This fcn is called when user resizes GUI
global DATA

if isempty(DATA)
  %This prevents validate callbacks from reporting error when opening
  %segment.fig for inspection. 
  return;
end

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
  
  rpwidth=round(min(DATA.GUISettings.ReportPanelPixelMax,pfig(3)*DATA.GUISettings.RightGapWidth)); %0.21DATA.GUISettings.RightGapWidth));
  set(DATA.Handles.reportpanel,'position',[...
    pfig(3)-rpwidth ...
    p(2) ...  
    rpwidth ...
    p(4)]);
  set(DATA.Handles.reportpanel,'units',panelunits);  
  
  if any(strcmp(DATA.ProgramName,{'Segment 3DPrint'}))
      panelunits = get(DATA.Handles.printuipanel,'units');
      set(DATA.Handles.printuipanel,'units','pixels');
       rpwidth=round(pfig(3)*DATA.GUISettings.RightGapWidth); %0.21DATA.GUISettings.RightGapWidth));
      set(DATA.Handles.printuipanel,'position',[...
      pfig(3)-rpwidth ...
      p(2) ...
      rpwidth ...
      p(4)]);
    set(DATA.Handles.printuipanel,'units',panelunits)

    if ~isempty(DATA.Handles.iconholder2.cdata)
      DATA.Handles.iconholder2.render
    end
    
    if ~isempty(DATA.Handles.iconholder3.cdata)
      DATA.Handles.iconholder3.render
    end
    
    if ~isempty(DATA.Handles.iconholder4.cdata)
      DATA.Handles.iconholder4.render
    end
  end
  
  
  DATA.GUISettings.RightGapWidth = rpwidth/pfig(3);
  
  %Render iconplaceholders aswell
  DATA.Handles.toggleiconholder.render
  if not(contains(DATA.ProgramName,'3D')) 
    DATA.Handles.permanenticonholder.render
    if ~isempty(DATA.Handles.hideiconholder.cdata)
      DATA.Handles.hideiconholder.render
    end
  end
  
  if ~isempty(DATA.Handles.configiconholder.cdata)
    DATA.Handles.configiconholder.render
  end    
    
  try
    if ~isempty(DATA.ViewMatrix)
      rows=DATA.ViewMatrix(1);
      cols=DATA.ViewMatrix(2);
      if length(DATA.ViewPanelsType) == 4 && all(strcmp(DATA.ViewPanelsType, {'orth', 'hla', 'vla', 'gla'}))
        viewfunctions('setview',rows,cols,DATA.ViewPanels,DATA.ViewPanelsType);
      else
        viewfunctions('setview',rows,cols); %drawfunctions('drawall',rows,cols);
      end
    end
  catch me
    mydispexception(me)
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
  end
end

%---------------------------------
function renderstacksfromdicom(no) %#ok<DEFNU>
%---------------------------------
%Render image stacks in main gui. This function is typically called upon
%loading.

%Do not mess with .Silent here since it will be taken care of by lower
%routine images.

%New routine doesnt need this function perhaps

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
      DATA.ViewPanelsMatrix = {[cols1 rows1] [cols2 rows2]};
    else
      DATA.ViewPanels=1;
      if (SET(no).ZSize==1)
        DATA.ViewPanelsType{1} = 'one';
      else
        DATA.ViewPanelsType{1} = DATA.GUISettings.ViewPanelsTypeDefault;
      end
      if strcmp(DATA.ViewPanelsType,'montage')
        [rows,cols] = calcfunctions('calcrowscols',1,SET(1).ZSize);
      else
          rows = 1;
          cols = 1;
      end
      DATA.ViewPanelsMatrix = {[cols rows]}; 
      DATA.ViewMatrix=[1 1];
    end
  end

%   if (~DATA.Preview.Silent)
%     %Normal loading
%     
% %DATA.switchtoimagestack(no,true); %force
%     %drawfunctions('drawthumbnails',isempty(DATA.DATASETPREVIEW));
%     
%     %The refresh starts everything
%     %viewfunctions('setview',DATA.ViewMatrix(1),DATA.ViewMatrix(2))
%   end;

end

disp('Files loaded.');

%endoffcalculation;

%----------------------------
function update_thumbnail(nos)
%----------------------------
%This fcn updates thumbnail no
global DATA SET

%if DATA.Silent
%  return
%end

%Check if empty DATASETPREVIEW and generate if so.
if isempty(DATA.DATASETPREVIEW)
  calcfunctions('calcdatasetpreview');
end

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
end

if ~DATA.Silent
  set(DATA.Handles.datasetpreviewimage,'cdata',DATA.DATASETPREVIEW);
end

%---------------------------
function out=thumbnailno(in)
%---------------------------
%Helper fcn to remember what image stack were clicked.
persistent no

if nargin==1
  no = in;
end
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
    
    %then we want the fig position
    [x,y] = mygetcurrentpoint(DATA.fig);
    
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
    %Set up
    set(DATA.fig,'WindowButtonUpFcn',...
      'segment(''thumbnail_Buttonup'')');
    %Get clicked position
    [x,y] = mygetcurrentpoint(DATA.Handles.datasetaxes);
    
    %Find clicked image stack
    no = getclickedpreview(x,y);
    thumbnailno(no); %store
    
    %add no in DATA.lastobject field so that it can be loaded into panels
    %if user clicks this
    DATA.LastObject = no;

    %Bring up popup menu
    [p(1),p(2)] = mygetcurrentpoint(DATA.fig);
    set(DATA.Handles.datasetpreviewmenu,...
      'Position',p,...
      'Visible','on');
end

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
global DATA NO
%Get coordinate
[x,y] = mygetcurrentpoint(DATA.Handles.boxaxes);

if nargin==1
  x=-1;
end
  
%Restore motion etc
set(DATA.fig,'WindowButtonMotionFcn',@DATA.toggleplaceholdermotion);
set(DATA.fig,'WindowButtonUpFcn','buttonupfunctions(''buttonup_Callback'')');
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
    viewfunctions('switchpanel',allreadyout(1));
    return;
  end
  ind = find(DATA.ViewPanels==0,1);
else
  if strcmp(DATA.CurrentTheme,'3dp') %If in 3dp mode disperse the image according to current view.
    
    DATA.LevelSet = [];
    
    segment3dp.tools('check3dpfields');
    segment3dp.tools('storetoobject');
    oldno = NO;
    NO = no;
    ok = segment3dp.tools('init3DP',0,1);
    if ~ok
      NO = oldno;
    end
    
    segment3dp.tools('update3DP');

    drawfunctions('drawthumbnailframes')
    return
  elseif strcmp(DATA.ViewPanelsType{1},'orth')
    %then we switch to single view mode
    viewfunctions('setview',1, 1, no,{'one'});
    drawfunctions('drawthumbnailframes')
    return
  else
    %Find at which image panel we drop it.
    %dist = zeros(1,length(DATA.Handles.imageaxes));
    dist=zeros(1,length(DATA.ViewPanels));
   % for loop=1:length(DATA.Handles.imageaxes)
   for loop=1:length(DATA.ViewPanels)
      p = get(DATA.Handles.imageaxes(loop),'position');
      p = p(1:2)+0.5*p(3:4); %Center position
      dist(loop) = sqrt(sum(([x y]-p).^2)); %Euklidean distance
    end
    [~,ind] = min(dist);
    ind=ind(1);   %just in case equal distance..!
  end
end

if ~isempty(ind)
    viewfunctions('addno2panel',ind,no)  
end

%--------------------------
function z = remap(im,cmap,c,b) %#ok<DEFNU>
%--------------------------
%Remap data according to cmap

global SET NO

if nargin<2
  cmap = SET(NO).Colormap;
end
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
  end
end


%-------------------------
function updatemeasurement %#ok<DEFNU>
%-------------------------
%calculate measurement and graphically update
global DATA

DATA.updatemeasurementreport;

%updateplot
DATA.updatemeasurementaxes;

%----------------------------------
function updatevolume(lvsegchanged)
%----------------------------------
%Calc volume of segmentation and graphically update.
global DATA SET NO

if nargin < 1
  lvsegchanged = false;
end

rotstring = '';
specstring = newline;

%Error check
if isempty(SET(NO).ImageViewPlane)
  SET(NO).ImageViewPlane = 'Unspecified';
end

calcd = -1;
if ismember(SET(NO).ImageViewPlane,{'2CH','3CH','4CH'})
  
  [calcd,usednos] = longaxistools('calcbiplanevolume');
%   if all(isnan(SET(NO).LVV)) || ~ismember(NO,usednos)
%     calcd = 0;
%   end
  
  for no = DATA.LVNO
    calcfunctions('volume_helper',no);
  end

  %if current stack is a new contribution to a lax volume calculation in
  %DATA.LVNO
  %we need to add it to that
  if any(ismember(usednos,DATA.LVNO))
    LAX_group = findfunctions('findlaxset',1);
    LAX_group = LAX_group(LAX_group~=0);
    if ~isempty(LAX_group)
      DATA.LVNO = LAX_group;
      str = [];
      for i = 1:length(LAX_group)
        str = [str, num2str(LAX_group(i)),','];
      end
      str(end) = [];
      if isfield(DATA.Handles,'lvstackpushbutton')
        set(DATA.Handles.lvstackpushbutton,'String',sprintf('Stack #%s',str));
      end
    else
      if isfield(DATA.Handles,'lvstackpushbutton')
        set(DATA.Handles.lvstackpushbutton,'String',sprintf('Stack #%d',DATA.LVNO(1)));
      end
    end
  end
end

if calcd == -1
  calcfunctions('calcvolume',NO);
else
  calcfunctions('volume_helper',NO);
end

% if calcd > 0
%   for no = DATA.LVNO
%     calcfunctions('volume_helper',no);
%   end
% elseif calcd == -1
%   calcfunctions('calcvolume',NO);
% else
%   calcfunctions('volume_helper',NO);
% end

if DATA.Silent
  return;
end

%update all reports
DATA.updatelvreport
DATA.updatervreport
% DATA.updateflowreport

%updateplot
DATA.updatevolumeaxes
%DATA.updatetimebaraxes


% if strcmp(DATA.ProgramName,'Segment')
%     %update all reports
%     DATA.updatelvreport
%     DATA.updatervreport
%     %DATA.updatemeasurementreport
% 
%     %updateplot
%     DATA.updatevolumeaxes
% elseif DATA.LVNO == NO || DATA.RVNO == NO
%     %update all reports
%     DATA.updatelvreport
%     DATA.updatervreport
%     %DATA.updatemeasurementreport
%     
%     %updateplot
%     DATA.updatevolumeaxes
% end

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
end

%--------------------------------
function no=getclickedpreview(~,y)
%--------------------------------
%function which returns the clicked preview image
global DATA

if nargin==0
  [~,y] = mygetcurrentpoint(DATA.Handles.datasetaxes);
end

no = floor(y/DATA.GUISettings.ThumbnailSize)+1;
    

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
end

%----------------------------------------
function smoothendowall_Callback(no)
 %----------------------------------------
%Adjust the LV wall so that the thickness is more even. Adjustment is done
%on the endocardial side.

%Einar Heiberg

global SET NO DATA

if nargin<1
  no = NO;
end
tools('connectinterpolation',no,{'EndoInterp','EpiInterp'});
%We need to have both endo and epi to do this
endotfs = findfunctions('findframeswithsegmentation','Endo',no,SET(no).CurrentSlice);
epitfs = findfunctions('findframeswithsegmentation','Epi',no,SET(no).CurrentSlice);

if ~(endotfs(SET(no).CurrentTimeFrame) && epitfs(SET(no).CurrentTimeFrame))
    myfailed('Need both Endo and Epi segmentation to even out wall')
    return
end

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
segment('updatevolume');
drawfunctions('drawno',no)

%--------------------------
function stopmovie_Callback
%--------------------------
%End movie display.
global DATA 

DATA.Run = 0;

try
  stateandicon=viewfunctions('iconson','play');
  stateandicon{2}.undent;
  DATA.Handles.permanenticonholder.render;
catch
end

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
    end
    myworkon;
    viewfunctions('setview',1,1,cineno,{'montage'});
    myworkoff;
  case 'lv'
    if ~isempty(DATA.LVNO)
      cineno = DATA.LVNO(1);
    else
      lvrvstacks = findfunctions('findcineshortaxisno',true);
      if not(isempty(lvrvstacks)), cineno = lvrvstacks(1); else cineno = []; end
      if not(isempty(cineno))
        DATA.LVNO = cineno;
      end
    end
    if ~isempty(cineno)
      myworkon;
      viewfunctions('setview',1,1,cineno,{'montage'});
      myworkoff;
    else
      myfailed('No LV stack found.',DATA.GUI.Segment);
      return;
    end
  case 'rv'    
    if ~isempty(DATA.RVNO)
      cineno = DATA.RVNO(1);
    else
      lvrvstacks = findfunctions('findcineshortaxisno',true);
      if not(isempty(lvrvstacks)), cineno = lvrvstacks(2); else cineno = []; end
      if not(isempty(cineno))
        DATA.RVNO = cineno;
      end
    end
    if ~isempty(cineno)
      myworkon;
      viewfunctions('setview',1,1,cineno,{'montage'});
      myworkoff;
    else
      myfailed('No RV stack found.',DATA.GUI.Segment);
      return;
    end
  case 'cinescar'
    if isempty(DATA.LVNO)
      cineno = findfunctions('findcineshortaxisno',false, 28);
    else
      cineno = DATA.LVNO(1);
    end
    scarno = findfunctions('findscarshortaxisno');
    if isempty(scarno)
        myfailed('No scar stack found.',DATA.GUI.Segment);
      return;
    end
    myworkon;
    if ~isempty(cineno) && ~isempty(scarno)
      % plot both scar and cine image
      viewfunctions('setview',2,1,[cineno,scarno],{'montagerow','montagerow'},[scarno 2]);
    else
      % only plot scar image
      viewfunctions('setview',1,1,[scarno],{'montage'});
    end
    myworkoff;
  case 'flow'
    flowno = DATA.FlowNO; 
    if isempty(flowno)
      [~,~,flowno] = findfunctions('findno');
    end
    if isempty(flowno)
      myfailed('No flow stacks found.',DATA.GUI.Segment);
      return;
    end
    haveflow = zeros(length(flowno),1); 
    for flownoloop = 1:length(flowno)
      haveflow(flownoloop) = not(isempty(SET(flowno(flownoloop)).Flow)); 
    end
    if sum(haveflow) == 0
      myfailed('No flow stacks found.',DATA.GUI.Segment);
      return;    
    end
    flowno = flowno(haveflow==1);
    no = flowno(1); %For now take first, maybe late do a toggle...    
%     DATA.FlowNO = no; %set also flowno to be found
    nom = SET(no).Flow.MagnitudeNo;
    nop = SET(no).Flow.PhaseNo;
    if isempty(nop)
      myfailed('Flow view only works for flow image stacks with trough-plane flow.',DATA.GUI.Segment);
      return;
    end
    myworkon;
    viewfunctions('setview',1,2,[nom,nop],{'one','one'});
    myworkoff;
  case 'perfusion'
    stressno = findfunctions('findstack','Perfusion Stress');
    restno = findfunctions('findstack','Perfusion Rest');      
    if isempty(stressno) && isempty(restno)
      myfailed('No perfusion stacks found.',DATA.GUI.Segment);
      return;
    end
    myworkon;
    if isempty(stressno)      
      viewfunctions('setview',1,1,restno,{'montagerow'});
    elseif isempty(restno)
      viewfunctions('setview',1,1,restno,{'montagerow'});
    else
      viewfunctions('setview',2,1,[stressno,restno],{'montagerow','montagerow'});
    end
    myworkoff;
  case 'strainsax'
    ftno = strcmp('Feature tracking',{SET.ImageType});
    saxno = strcmp('Short-axis',{SET.ImageViewPlane});
    [points,strainnos] = sort(ftno+saxno,'descend');
    strainno = strainnos(1);
    if points(1)<1 || isempty(strainno)
      myfailed('No strain stack found',DATA.GUI.Segment);
      return;
    end
    myworkon;
    viewfunctions('setview',1,1,strainno,{'montagerow'});
    myworkoff;
  case 'strainlax'
    twochno = strcmp('2CH',{SET.ImageViewPlane});
    threechno = strcmp('3CH',{SET.ImageViewPlane});
    fourchno = strcmp('4CH',{SET.ImageViewPlane});
    for noloop = 1:length(SET)
      havestrain(noloop) = not(isempty(SET(noloop).StrainTagging));
    end
    [points,strain2chnos] = sort(twochno+havestrain,'descend');
    [points,strain3chnos] = sort(threechno+havestrain,'descend');
    [points,strain4chnos] = sort(fourchno+havestrain,'descend');
    strain2chno = strain2chnos(1);
    strain3chno = strain3chnos(1);
    strain4chno = strain4chnos(1);
    if isempty(strain2chno) && isempty(strain3chno) && isempty(strain4chno)
      myfailed('No strain stack found',DATA.GUI.Segment);
      return;
    end
    myworkon;
    if not(isempty(strain2chno)) && not(isempty(strain3chno)) && not(isempty(strain4chno))
      viewfunctions('setview',1,3,[strain2chno,strain3chno,strain4chno],{'one','one','one'});
    elseif not(isempty(strain2chno)) && not(isempty(strain3chno))
      viewfunctions('setview',1,2,[strain2chno,strain3chno],{'one','one'});
    elseif not(isempty(strain2chno)) && not(isempty(strain4chno))
      viewfunctions('setview',1,2,[strain2chno,strain4chno],{'one','one'});
    elseif not(isempty(strain3chno)) && not(isempty(strain4chno))
      viewfunctions('setview',1,2,[strain3chno,strain4chno],{'one','one'});
    elseif not(isempty(strain2chno))
      viewfunctions('setview',1,1,strain2chno,{'montagerow'});
    elseif not(isempty(strain3chno))
      viewfunctions('setview',1,1,strain3chno,{'montagerow'});
    elseif not(isempty(strain4chno))
      viewfunctions('setview',1,1,strain4chno,{'montagerow'});
    end
    myworkoff;
  case 'cinescarperf'
  case 'stress'
  otherwise
    myfailed('Unknown option to viewspecial_Callback',DATA.GUI.Segment);
    return;
end
%adjust visible thumbnails so that NO is visible
nosafter = intersect(NO:NO+DATA.Pref.NumberVisibleThumbnails-1,1:length(SET));
nosbefore = intersect(NO-DATA.Pref.NumberVisibleThumbnails+1:NO,1:length(SET));
newnos = union(nosbefore,nosafter);
if length(SET) > DATA.Pref.NumberVisibleThumbnails
  newnos = newnos(1:DATA.Pref.NumberVisibleThumbnails);
  if min(DATA.VisibleThumbnails)-NO > 0
    whattodo = 'up';
  else
    whattodo = 'down';
  end
  DATA.VisibleThumbnails = newnos;
  thumbnailslider_Callback(whattodo);
end
segment('updatemeasurement');
DATA.updatetimebaraxes

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
end

if f < 1 && strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},'montagesegmented')
  mfzs = SET(NO).MontageFitZoomState;
  if mfzs(2)-mfzs(1)>=size(DATA.ViewIM{1},1) && ...
      mfzs(4)-mfzs(3)>=size(DATA.ViewIM{1},2)
    %viewimage_Callback('montage');
    viewfunctions('setviewtype','montage');
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
end

xlim = get(DATA.Handles.imageaxes(panel),'xlim');
ylim = get(DATA.Handles.imageaxes(panel),'ylim');
for ploop=panelstodo
  switch DATA.ViewPanelsType{ploop}
    case 'montage'
      SET(DATA.ViewPanels(ploop)).MontageZoomState = [xlim(:);ylim(:)];
    case 'montagerow'
      SET(DATA.ViewPanels(ploop)).MontageRowZoomState = [xlim(:);ylim(:)];
    case {'montagefit','sax3','montagesegmented'}
      SET(DATA.ViewPanels(ploop)).MontageFitZoomState = [xlim(:);ylim(:)];
    case {'one','orth'}
      SET(DATA.ViewPanels(ploop)).NormalZoomState = [xlim(:);ylim(:)];
    case 'hla'
      SET(DATA.ViewPanels(ploop)).HLA.ZoomState = [xlim(:);ylim(:)];
    case 'vla'
      SET(DATA.ViewPanels(ploop)).VLA.ZoomState = [xlim(:);ylim(:)];
  end
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

if strcmp(DATA.CurrentTheme,'3dp')
  segment3dp.tools('zoomin_Callback')
  return
end

if strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},'mmodetemporal')
  updatemmode([],'one')
  return
end

if strcmp(get(gcbf,'currentmodifier'),'shift')
  for loop=1:length(DATA.Handles.imageaxes)
    zoomhelper(DATA.Handles.imageaxes(loop),1.2,DATA.ViewPanels(loop),loop);
  end
else
  zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1.2);
end
%drawfunctions('viewupdatetextposition');
%drawfunctions('viewupdateannotext');

%----------------------------
function viewzoomout_Callback  %#ok<DEFNU>
%----------------------------
%Zooms out in current view panel
global DATA

if strcmp(DATA.CurrentTheme,'3dp')
  segment3dp.tools('zoomout_Callback')
  return
end

if strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},'mmodetemporal')
  updatemmode([],'three')
  return
end

if strcmp(get(gcbf,'currentmodifier'),'shift')
  for loop=1:length(DATA.Handles.imageaxes)
    zoomhelper(DATA.Handles.imageaxes(loop),1/1.2,DATA.ViewPanels(loop),loop);
  end
else
  zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1/1.2);
end
%drawfunctions('viewupdatetextposition');
%drawfunctions('viewupdateannotext');

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
end

set(ca,'xlim',xlim,'ylim',ylim);

%--------------------------------------
function viewmanualinteraction_Callback  %#ok<DEFNU>
%--------------------------------------
%Displays a GUI indicating in which slices and timeframes manual
%interaction have been made for LV segmentation.
global DATA SET NO

if isempty(SET(NO).EndoX) || all(isnan(SET(NO).EndoX(:)))
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end

fig = mygui('manualedit.fig');
% Generate a structure of handles to pass to callbacks, and store it.
handles = fig.handles;
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
        end
      end
    end
  end
end

axes(handles.endoaxes);  
image([1 SET(NO).TSize]-0.5,[1 SET(NO).ZSize]-0.5,im);
set(gca,'xtick',1:SET(NO).TSize,'ytick',1:SET(NO).ZSize);
ylabel(dprintf('Apex    =>     Base'));
xlabel(dprintf('Timeframe'));
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
        end
      end
    end
  end
end

axes(handles.epiaxes); 
image([1 SET(NO).TSize]-0.5,[1 SET(NO).ZSize]-0.5,im);
set(gca,'xtick',1:SET(NO).TSize,'ytick',1:SET(NO).ZSize);
ylabel(dprintf('Apex    =>     Base'), 'Color', DATA.GUISettings.ForegroundColor);
xlabel(dprintf('Timeframe'), 'Color', DATA.GUISettings.ForegroundColor);
grid on;

%Make plot axis color match Segment's theme
set([handles.endoaxes handles.epiaxes],...
  'XColor', DATA.GUISettings.ForegroundColor,...
  'YColor', DATA.GUISettings.ForegroundColor);

%Set legend color
set(handles.frame1,'Backgroundcolor', [0 0 1]); %blue
set(handles.frame2,'Backgroundcolor', [1 0 1]); %pink
set(handles.frame3,'Backgroundcolor', [1 0 0]); %red
set(handles.frame4,'Backgroundcolor', [0 1 1]); %cyan

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
  end
else
  %User timer info
  t = 1+mod(floor(rem(now-DATA.StartTime,1)*24*3600/(SET(NO).BeatTime/SET(NO).TSize)+DATA.StartFrame),SET(NO).TSize);
end

%-----------------------------
function volumeaxes_Buttondown %#ok<DEFNU>
%-----------------------------
%Called when user has clicked volume graph, sets current
%timeframe to clicked point.

global DATA SET NO

if any(NO == DATA.LVNO) || strcmp(DATA.ProgramName,'Segment')
  
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
    end
  end
  if abs(SET(NO).EST-t)<3
    if abs(y-(0.25*(ylim(2)-ylim(1))+ylim(1)))<(0.2*(ylim(2)-ylim(1))) %ES textg is not really at bottom.
      esed_Buttondown('es');
      return;
    end
  end
  DATA.updatevolumeaxes;
  panels = find(DATA.ViewPanels == NO);
  for p = panels
      DATA.ViewIM{p} = [];
      drawfunctions('drawpanel',p);
  end
  viewfunctions('updatetimebars')
end

%--------------------------
function esed_Buttonup(type) %#ok<DEFNU>
%---------------------------
%Buttonup function when dragging ES or ED markers in volume graph
global DATA SET NO


if any(NO ==  DATA.LVNO) || strcmp(DATA.ProgramName,'Segment')
  %Get position
  [x] = mygetcurrentpoint(gca);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);  
t=max([t,1]);
  %Store position
  switch type
    case 'es'
      SET(NO).EST = t;
    case 'ed'
      SET(NO).EDT = t;
  end
  
  %disp(sprintf('ed:%d es:%d',SET(NO).EDT,SET(NO).EST));
  
  %Restore
  set(DATA.fig,'WindowButtonUpFcn','buttonupfunctions(''buttonup_Callback'')');
  set(DATA.fig,'WindowButtonMotionFcn','');
  SET(NO).CurrentTimeFrame = t;
  DATA.updatevolumeaxes
  %we also need to update the timebar ed es
  switch type
      case 'es'
          p = get(DATA.Handles.estimebartext,'position');
          p(1) = SET(NO).TimeVector(t)*1000;
          set(DATA.Handles.estimebartext,'position',p);
          set(DATA.Handles.estimebarline,'XData',[SET(NO).TimeVector(t),SET(NO).TimeVector(t)]*1000);
      case 'ed'
          p = get(DATA.Handles.edtimebartext,'position');
          p(1) = SET(NO).TimeVector(t)*1000;
          set(DATA.Handles.edtimebartext,'position',p);
          set(DATA.Handles.edtimebarline,'XData',[SET(NO).TimeVector(t),SET(NO).TimeVector(t)]*1000);
  end
  
  updatevolume;
  for p = find(DATA.ViewPanels==NO)
      drawfunctions('drawpanel',p);
  end
  
  viewfunctions('updatetimebars')
end

%-------------------------
function esed_Motion(type) %#ok<DEFNU>
%-------------------------
%Motion function when dragging ES or ED markers in volume graph.
global DATA SET NO

if any(NO == DATA.LVNO) || strcmp(DATA.ProgramName,'Segment')
  [x] = mygetcurrentpoint(gca);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
  t=max([t,1]);
  switch type
    case 'es'
      p = get(DATA.Handles.estext,'position');
      p(1) = SET(NO).TimeVector(t)*1000;
      set(DATA.Handles.estext,'position',p);
      SET(NO).EST = t;
      set(DATA.Handles.esline,'xdata',[p(1) p(1)]);
    case 'ed'
      p = get(DATA.Handles.edtext,'position');
      p(1) = SET(NO).TimeVector(t)*1000;
      set(DATA.Handles.edtext,'position',p);
      SET(NO).EDT = t;
      set(DATA.Handles.edline,'xdata',[p(1) p(1)]);
  end
  
  SET(NO).CurrentTimeFrame = t;
  panels = find(DATA.ViewPanels == NO);
  for p = panels
      DATA.ViewIM{p} = [];
      drawfunctions('drawpanel',p);
  end

end

%-----------------------------
function esed_Buttondown(type) 
%-----------------------------
%Buttondown function when dragging ES or ED markers in volume graph.
global DATA SET NO

if any(NO == DATA.LVNO) || strcmp(DATA.ProgramName,'Segment')
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
  end
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
  end
end
if abs(SET(NO).EST-t)<3
  if abs(y-(0.25*(ylim(2)-ylim(1))+ylim(1)))<(0.2*(ylim(2)-ylim(1))) %ES textg is not really at bottom.
    esedtimebar_Buttondown('es');
    return;
  end
end

nos = SET(NO).Linked;
for loop=setdiff(nos,NO)
  SET(loop).CurrentTimeFrame = SET(NO).CurrentTimeFrame;
  if SET(loop).CurrentTimeFrame > SET(loop).TSize
    SET(loop).CurrentTimeFrame = SET(loop).TSize;
  end
end

for p = find(ismember(DATA.ViewPanels,nos))
    DATA.ViewIM{p} = [];
  drawfunctions('drawpanel',p);
end

%update timebars
viewfunctions('updatetimebars')
set(DATA.fig,'WindowButtonUpFcn','buttonupfunctions(''buttonup_Callback'')');
set(DATA.fig,'WindowButtonMotionFcn','segment(''timebar_Motion'')');

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
end

%Restore
set(DATA.fig,'WindowButtonUpFcn','buttonupfunctions(''buttonup_Callback'')');
set(DATA.fig,'WindowButtonMotionFcn',@DATA.toggleplaceholdermotion);
SET(NO).CurrentTimeFrame = t;
updatevolume;
panels = find(DATA.ViewPanels == NO);
for p = panels
    DATA.ViewIM{p} = [];
    drawfunctions('drawpanel',p);
end

%update timebars
viewfunctions('updatetimebars')

%-------------------------
function esedtimebar_Motion(type) %#ok<DEFNU>
%-------------------------
%Motion function when dragging ES or ED markers in time bar graph.
global DATA SET NO

[x] = mygetcurrentpoint(gca);
xlim=get(DATA.Handles.timebaraxes,'xlim');
[~,t]=min((SET(NO).TimeVector/SET(NO).TimeVector(end)-x/xlim(2)).^2);
t=max([t,1]);

switch type
  case 'es'
    p = get(DATA.Handles.estimebartext,'position');
    p(1) = SET(NO).TimeVector(t)*1000;
    set(DATA.Handles.estimebartext,'position',p);
    SET(NO).EST = t; 
    set(DATA.Handles.estimebarline,'XData',[SET(NO).TimeVector(t),SET(NO).TimeVector(t)]*1000);    
  case 'ed'
    p = get(DATA.Handles.edtimebartext,'position');
    p(1) = SET(NO).TimeVector(t)*1000;
    set(DATA.Handles.edtimebartext,'position',p);
    SET(NO).EDT = t;        
    set(DATA.Handles.edtimebarline,'XData',[SET(NO).TimeVector(t),SET(NO).TimeVector(t)]*1000);    
end

SET(NO).CurrentTimeFrame = t;
panels = find(DATA.ViewPanels == NO);
for p = panels
    DATA.ViewIM{p} = [];
    drawfunctions('drawpanel',p);
end

%update timebars
viewfunctions('updatetimebars')

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
end

%update timebars
viewfunctions('updatetimebars')

%-----------------------------
function flowaxes_Buttondown %#ok<DEFNU>
%-----------------------------
%Called when user has clicked flow graph, sets current
%timeframe to clicked point.

global DATA SET NO

if ismember(NO,SET(DATA.FlowNO).Linked) || strcmp(DATA.ProgramName,'Segment')
  
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
    end
  end
  
  panels = find(ismember(DATA.ViewPanels ,nos));
  for p = panels
      DATA.ViewIM{p} = [];
      drawfunctions('drawpanel',p);
  end
  
  viewfunctions('updatetimebars')
end

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
  end
  ok = true;
  DATA.Interactionlock = true;
else
  ok = false;
  myworkoff;
  disp('Warning: No slices selected => command ignored.');
end

%-------------------------
function endoffcalculation %#ok<DEFNU>
%-------------------------
%Turns off interaction lock, called after a calculation. See above.
global DATA 
DATA.Interactionlock = false;
myworkoff;

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
end

%------------------------------
function updatemmodevisibility %#ok<DEFNU>
%------------------------------
%Updates mmode visibility.
global DATA SET NO

rightstate = 'on';
leftstate = 'on';

if (SET(NO).TSize<2)||(isequal(DATA.ViewPanelsType{DATA.CurrentPanel},'flow'))
  rightstate = 'off';
end

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
% DATA.updateaxestables('measure');
% set(DATA.Handles.distancetext,'string',...
%   sprintf('Distance:%3.2f [mm]\tTime:%3.2f [ms]',dist,timedist));

%------------------------------------
function updatemmode(arg,nbrofcycles) 
%------------------------------------
%Calculate and show mmode image
global DATA SET NO

if SET(NO).TSize<2
  return;
end

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
  end
else
  for tloop=1:SET(NO).TSize
    temp(:,tloop + SET(NO).TSize*[0 1 2]) = repmat(interp2( ...
      SET(NO).IM(:,:,tloop,SET(NO).CurrentSlice),xi,yi,'nearest')',1,3);
  end
end

%Display data
temp = calcfunctions('remapuint8',temp,NO);

panel = [];
for loop=1:length(DATA.ViewPanelsType)
  if isequal(DATA.ViewPanelsType{loop},'mmodetemporal')
    panel = loop;
  end
end
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
end

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

% %-------------------------------------
% function mmodecenter_Buttondown(panel) %#ok<DEFNU>
% %-------------------------------------
% %Called when mmode center is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodecenter_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %--------------------------------
% function mmode1_Buttondown(panel) %#ok<DEFNU>
% %--------------------------------
% %Called when mmode1 is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode1_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %--------------------------------
% function mmode2_Buttondown(panel) %#ok<DEFNU>
% %--------------------------------
% %Called when mmode1 is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode2_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %------------------------------------
% function mmode1line_Buttondown(panel) %#ok<DEFNU>
% %------------------------------------
% %Called when mmode1line is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode1line_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %------------------------------------
% function mmode2line_Buttondown(panel) %#ok<DEFNU>
% %------------------------------------
% %Called when mmode1line is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmode2line_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %-------------------------------------
% function mmodepoint1_Buttondown(panel) %#ok<DEFNU>
% %-------------------------------------
% %Called when mmodepoint1 is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodempoint1_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %-------------------------------------
% function mmodepoint2_Buttondown(panel) %#ok<DEFNU>
% %-------------------------------------
% %Called when mmodepoint1 is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodempoint2_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %---------------------------------------
% function mmodetimebar1_Buttondown(panel) %#ok<DEFNU>
% %---------------------------------------
% %Called when mmodepoint1 is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodetimebar1_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));
% 
% %---------------------------------------
% function mmodetimebar2_Buttondown(panel) %#ok<DEFNU>
% %---------------------------------------
% %Called when mmodepoint1 is pressed down, sets motion and buttonup
% %function.
% global DATA 
% 
% switchtopanel(panel);
% 
% if DATA.Interactionlock
%   return;
% end
% 
% set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''mmodetimebar2_Motion'');',mfilename));
% set(DATA.imagefig,'WindowButtonUpFcn',...
%   sprintf('%s(''mmode_Buttonup'')',mfilename));

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
end

% for nloop=1:length(nos)
%   drawfunctions('updatenopanels',nos(nloop));
% end
viewfunctions('setview');

%updatevolume;
drawnow('expose'); %Expose does not other callbacks to evaluate.


%-----------------------------------------------------
function [xout,yout] = checkconsistencyhelper(xin,yin,numpoints)
%-----------------------------------------------------
%Make sure that the contour is counter clock-wise and that it starts at
%three o clock. Also ensures that the points are evenly distributed.

global DATA

if nargin == 2
    numpoints = DATA.NumPoints;
end

%--- Make sure evenly distributed
diffx = conv2(xin,[1;-1],'valid');
diffy = conv2(yin,[1;-1],'valid');
diffvec = sqrt(diffx.*diffx+diffy.*diffy);

if any((diffvec-mean(diffvec))/mean(diffvec)>2)
  %disp('Problem with endo => fixed, report to support');
end
contourlength = cumsum(diffvec);
xout = interp1(contourlength,xin(2:end),linspace(contourlength(1),contourlength(end),numpoints-1));
yout = interp1(contourlength,yin(2:end),linspace(contourlength(1),contourlength(end),numpoints-1));
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
end

[~,inda] = min(angle(complex(mx-xr,my-yr)));
xout(1:(numpoints-inda+1)) = yr(inda:end);
yout(1:(numpoints-inda+1)) = xr(inda:end);
xout((numpoints+1-inda):end) = yr(1:inda);
yout((numpoints+1-inda):end) = xr(1:inda);
        
%--------------------------------------------
function checkconsistency(timeframes,slice,no)
%---------------------------------------------
%Check consistency, to prevent earlier
%manual segmentations that have problems with
%direction lef/right
global SET NO DATA

if nargin==0
  timeframes=SET(NO).CurrentTimeFrame; %1;
end

if nargin<2
  slice=SET(NO).CurrentSlice;
end

if nargin<3
  no = NO;
end

if not(isempty(SET(no).EndoX)) && not(size(SET(no).EndoX,1)==DATA.NumPoints)
  tempendox = nan(DATA.NumPoints-1,SET(no).TSize,SET(no).ZSize);
  tempendoy = tempendox;
  timeframes = 1:SET(no).TSize;
  slice = 1:SET(no).ZSize;  
elseif isempty(SET(no).EndoX)
  tempendox = []; % nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
  tempendoy = []; % tempendox;
else
  tempendox = SET(no).EndoX(1:end-1,:,:);
  tempendoy = SET(no).EndoY(1:end-1,:,:);
end
if not(isempty(SET(no).EpiX)) && not(size(SET(no).EpiX,1)==DATA.NumPoints)
  tempepix = nan(DATA.NumPoints-1,SET(no).TSize,SET(no).ZSize);
  tempepiy = tempepix;
  timeframes = 1:SET(no).TSize;
  slice = 1:SET(no).ZSize;  
elseif isempty(SET(no).EpiX)
  tempepix = []; % nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);;
  tempepiy = []; % tempendox;
else
  tempepix = SET(no).EpiX(1:end-1,:,:);
  tempepiy = SET(no).EpiY(1:end-1,:,:);
end
if not(isempty(SET(no).RVEndoX)) && not(size(SET(no).RVEndoX,1)==DATA.NumPoints)
  temprvendox = nan(DATA.NumPoints-1,SET(no).TSize,SET(no).ZSize);
  temprvendoy = tempendox;
  timeframes = 1:SET(no).TSize;
  slice = 1:SET(no).ZSize;  
elseif isempty(SET(no).RVEndoX)
  temprvendox = []; % nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
  temprvendoy = []; % tempendox;
else
  temprvendox = SET(no).RVEndoX(1:end-1,:,:);
  temprvendoy = SET(no).RVEndoY(1:end-1,:,:);
end
if not(isempty(SET(no).RVEpiX)) && not(size(SET(no).RVEpiX,1)==DATA.NumPoints)
  temprvepix = nan(DATA.NumPoints-1,SET(no).TSize,SET(no).ZSize);
  temprvepiy = tempendox;
  timeframes = 1:SET(no).TSize;
  slice = 1:SET(no).ZSize;  
elseif isempty(SET(no).RVEpiX)
  temprvepix = []; % nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
  temprvepiy = []; % tempendox;
else
  temprvepix = SET(no).RVEpiX(1:end-1,:,:);
  temprvepiy = SET(no).RVEpiY(1:end-1,:,:);
end
for timeframe = timeframes
  for zloop=slice

    %Endo
    if ~isempty(SET(no).EndoX)
      if not(isnan(SET(no).EndoX(1,timeframe,zloop)))
        %Call consistenyhelper
        [tempendox(:,timeframe,zloop),tempendoy(:,timeframe,zloop)]=calcfunctions('resamplecurve',SET(no).EndoX(:,timeframe,zloop),SET(no).EndoY(:,timeframe,zloop),DATA.NumPoints-1);
%           [tempendox(:,timeframe,zloop),tempendoy(:,timeframe,zloop)] = checkconsistencyhelper(...
%             SET(no).EndoX(:,timeframe,zloop),...
%             SET(no).EndoY(:,timeframe,zloop));
      end
    end

    %Epi
    if ~isempty(SET(no).EpiX)
      if not(isnan(SET(no).EpiX(1,timeframe,zloop)))
        [tempepix(:,timeframe,zloop),tempepiy(:,timeframe,zloop)] = ...
          calcfunctions('resamplecurve',SET(no).EpiX(:,timeframe,zloop),SET(no).EpiY(:,timeframe,zloop),DATA.NumPoints-1);
%            [tempepix(:,timeframe,zloop),tempepiy(:,timeframe,zloop)] = ...
%           checkconsistencyhelper(...
%           SET(no).EpiX(:,timeframe,zloop),...
%           SET(no).EpiY(:,timeframe,zloop));
      end
    end
    
    %RVEndo
    if ~isempty(SET(no).RVEndoX)
      if not(isnan(SET(no).RVEndoX(1,timeframe,zloop)))
        [temprvendox(:,timeframe,zloop),temprvendoy(:,timeframe,zloop)] = ...
          calcfunctions('resamplecurve',SET(no).RVEndoX(:,timeframe,zloop),SET(no).RVEndoY(:,timeframe,zloop),DATA.NumPoints-1);
%         [temprvendox(:,timeframe,zloop),temprvendoy(:,timeframe,zloop)] = ...
%           checkconsistencyhelper(...
%           SET(no).RVEndoX(:,timeframe,zloop),...
%           SET(no).RVEndoY(:,timeframe,zloop)); 
      end
    end
    
    %RVEpi
    if ~isempty(SET(no).RVEpiX)
      if not(isnan(SET(no).RVEpiX(1,timeframe,zloop)))
        %Call consistenyhelper
        [temprvepix(:,timeframe,zloop),temprvepiy(:,timeframe,zloop)] = ...
          calcfunctions('resamplecurve', SET(no).RVEpiX(:,timeframe,zloop),SET(no).RVEpiY(:,timeframe,zloop),DATA.NumPoints-1);
%           checkconsistencyhelper(...
%           SET(no).RVEpiX(:,timeframe,zloop),...
%           SET(no).RVEpiY(:,timeframe,zloop));          
      end
    end
    
  end
end
      
if ~isempty(SET(no).EndoX)
  SET(no).EndoX = cat(1,tempendox,tempendox(1,:,:));
  SET(no).EndoY = cat(1,tempendoy,tempendoy(1,:,:));
end
if ~isempty(SET(no).EpiX)
  SET(no).EpiX = cat(1,tempepix,tempepix(1,:,:));
  SET(no).EpiY = cat(1,tempepiy,tempepiy(1,:,:));
end
if ~isempty(SET(no).RVEndoX)
  SET(no).RVEndoX = cat(1,temprvendox,temprvendox(1,:,:));
  SET(no).RVEndoY = cat(1,temprvendoy,temprvendoy(1,:,:));
end
if ~isempty(SET(no).RVEpiX)
  SET(no).RVEpiX = cat(1,temprvepix,temprvepix(1,:,:));
  SET(no).RVEpiY = cat(1,temprvepiy,temprvepiy(1,:,:));
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
end

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
end


%---------------------------
function resetlight_Callback %#ok<DEFNU>
%---------------------------
%Activated by toolbar icon, different from contrast_Callback (below).

global DATA SET NO

tools('enableundo',NO);
SET(NO).IntensityMapping.Contrast = 1;
SET(NO).IntensityMapping.Brightness = 0.5;
DATA.ViewIM{DATA.CurrentPanel} = [];
if ~isempty(DATA.Overlay)
  DATA.Overlay(DATA.CurrentPanel) = struct('alphadata',[],'cdata',[]);
end
update_thumbnail(NO);
drawfunctions('drawno',NO)
createfunctions('addcolorbar',DATA.CurrentPanel)
if DATA.Pref.UseLight 
  DATA.BalloonLevel = -1; %Force update of ballonimage
end

%---------------------------
function resetlightall_Callback %#ok<DEFNU>
%---------------------------
%Activated by toolbar icon, different from contrast_Callback (below).

global DATA SET

for no = 1:length(SET)
  SET(no).IntensityMapping.Contrast = 1;
  SET(no).IntensityMapping.Brightness = 0.5;
  update_thumbnail(no); %drawfunctions('drawsliceno',NO);
end

for p = find(DATA.ViewPanels)
    DATA.ViewIM{p} = [];
    drawfunctions('drawimages',p);
    createfunctions('addcolorbar',p)
end

if DATA.Pref.UseLight
  DATA.BalloonLevel = -1; %Force update of ballonimage
end
  
%-----------------------------------------------
function [xlim,ylim] = getbox(no,destno,doindex)
%-----------------------------------------------
%find x-limits and y-limits for a zoom box
    
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
    x = [x;tmp_x(:)];
    y = [y;tmp_y(:)];
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
  im = SET(no).IM;
  roisizetime = 150;
  roisizenotime = 200;
  if SET(no).TSize>1 && SET(no).ZSize~=1
    roisize = roisizetime;
  else
    roisize = roisizenotime;
  end
  
  nx = roisize/SET(no).ResolutionX;
  ny = roisize/SET(no).ResolutionY;
  [~,~,rlapfh] = autocrop(im,nx,ny,0,80,1,no);
  if ~isempty(rlapfh)
    limits = calcfunctions('rlapfh2xyz',destno,rlapfh(:,1),rlapfh(:,2),rlapfh(:,3));
    dodia=1;
  else
    xlim = [1,xsz];
    ylim = [1,ysz];
    return;
  end
end

xmin = min(limits(1,:));
ymin = min(limits(2,:));
xmax = max(limits(1,:));
ymax = max(limits(2,:));

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
      dia = norm([xmax,ymax]-[xmin,ymin])/2*0.8;
    case 2 %there is RV do very little margin addition
      dia = norm([xmax,ymax]-[xmin,ymin])/2*0.2;
    case 3 %from autocrop do "lagom" margin.
      dia = norm([xmax,ymax]-[xmin,ymin])/2*0.4;
  end
else
  dia=0;
end
xlim = [xmin-dia,xmax+dia];
ylim = [ymin-dia,ymax+dia];

if doindex
  xlim = round(xlim);
  ylim = round(ylim);
  xlim(xlim>xsz) = xsz;
  ylim(ylim>ysz) = ysz;
  xlim(xlim<1) = 1;
  ylim(ylim<1) = 1;
end

%----------------
function autozoom %#ok<DEFNU>
%----------------
%New autozoom autozooms everything that is displayed.

global DATA SET

panel = DATA.CurrentPanel;
panels = DATA.ViewPanels;
no = panels(panel);
origno = no;
zdirno = cross(SET(no).ImageOrientation(1:3),SET(no).ImageOrientation(4:6));

%Check which nos to do
nos = no;
for loop = 1:length(SET)

  loopno = loop;
  if ~isequal(loopno,0)
    zdir = cross(SET(loopno).ImageOrientation(1:3),SET(loopno).ImageOrientation(4:6));
    a = acos(abs(sum(zdir.*zdirno)))/pi*180;
    if a<20
      nos = [nos loopno];
    end
  end
end

%remove duplicates
nos = union(nos,[]);

%Loop over to find contours  
hasanyseg = false;
hasrv = false;
x = [];
y = [];

for loop = 1:length(nos)
  no = nos(loop);
  for type={'Endo','Epi','RVEndo','RVEpi'}
    tmp_x = SET(no).([type{1},'X']);
    tmp_y = SET(no).([type{1},'Y']);
    if ~isempty(tmp_x) && ~all(isnan(tmp_x(:)))
      hasanyseg = true;
      tmp_x = tmp_x(~isnan(tmp_x(:)));      
      tmp_y = tmp_y(~isnan(tmp_y(:)));
      
      %Convert to origno's coordinate system
      pos = calcfunctions('xyz2rlapfh',no,tmp_x,tmp_y,repmat(SET(no).ZSize/2,size(tmp_x)));
      pos = calcfunctions('rlapfh2xyz',origno,pos(:,1),pos(:,2),pos(:,3));
      
      x = [x;pos(1,:)'];
      y = [y;pos(2,:)'];
      if ~isempty(regexpi(type{1},'RV'))
        hasrv = true;
      end
    end
  end
end

if hasrv
  fmargin = 1.15;
else
  fmargin = 1.5;
end

if ~hasanyseg
  return;
end

%Compute limit & center (in origno coordinate system)
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);
xc = mean(x);
yc = mean(y);

%Loop over panels and check if update
for loop = 1:length(panels)
  
  panel = loop;
  no = panels(loop);
  
  if ismember(no,nos) && isequal(DATA.ViewPanelsType{panel},'one')

    %get previous zoomstate
    zoomstate = viewfunctions('getnewzoomstate',panel,no);
    
    %Convert xc,yc
    pos = calcfunctions('xyz2rlapfh',origno,[xc;xmin;xmax],[yc;ymin;ymax],repmat(SET(origno).ZSize/2,3,1));
    pos = calcfunctions('rlapfh2xyz',no,pos(:,1),pos(:,2),pos(:,3));
    xcno = pos(1,1); ycno = pos(2,1);
    xminno = pos(1,2); yminno = pos(2,2);
    xmaxno = pos(1,3); ymaxno = pos(2,3);
    
    oldxrange = zoomstate(4)-zoomstate(3);
    oldyrange = zoomstate(2)-zoomstate(1);
    oldxc = mean(zoomstate(3:4));
    oldyc = mean(zoomstate(1:2));
    translatex = oldxc-xcno; 
    translatey = oldyc-ycno;     
    
    %Translate
    oldxc = oldxc-translatex;
    oldyc = oldyc-translatey;
    
    %Find out zoom
    fx = max((oldxc-xminno)/(oldxrange/2),(xmaxno-oldxc)/(oldxrange/2));
    fy = max((oldyc-yminno)/(oldyrange/2),(ymaxno-oldyc)/(oldyrange/2));
    f = max(fx,fy)*fmargin;
    
    newzoomstate = [oldyc-oldyrange/2*f oldyc+oldyrange/2*f oldxc-oldxrange/2*f oldxc+oldxrange/2*f];
    
    SET(no).NormalZoomState = newzoomstate;
    
    %update it graphically
    viewfunctions('updatezoomandaspectratio',panel);
    
    %we want to update the text position and the frame
    viewfunctions('updatetextposition',panel)
        
  end
  
end

drawfunctions('drawselectedframe',DATA.CurrentPanel);

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
  %panel = DATA.ViewPanels(panelloop);
  DATA.ViewIM{panelloop} = [];
  drawfunctions('drawimages',panelloop);
end

%------------------------------
function autocontrast_Callback %#ok<DEFNU>
%------------------------------
%Automatically calculates contrast settings.
global NO
no = NO;
autocontrast(no);
segment('update_thumbnail',no)

%-------------------------
function autocontrast(no,silent)
%-------------------------
%Helper functionk to autocontrast_Callback_
global SET DATA

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
for p = find(DATA.ViewPanels==no)
    DATA.ViewIM{p} = [];
    drawfunctions('drawimages',p)
    createfunctions('addcolorbar',p)
end
% if ~silent 
%   drawfunctions('drawcontrastimage',no);
% end
myworkoff;

%-----------------------------------
function measuremove_Callback(dx,dy) %#ok<DEFNU>
%-----------------------------------
%Helper function to move measurements.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

n = DATA.MeasureN;
if (length(SET(no).Measure)<n)||(n<1)
  %Invalid measurement 
  return;
end
tools('enableundo',no);

SET(no).Measure(n).X = SET(no).Measure(n).X+dx;
SET(no).Measure(n).Y = SET(no).Measure(n).Y+dy;

viewfunctions('setview');  %drawfunctions('drawimageno');

%------------------------------
function measureexport_Callback %#ok<DEFNU>
%------------------------------
%Export measurements.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

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
  end

  clipboard('copy',stri);
  mymsgbox('Results copied to clipboard','Done!',DATA.GUI.Segment);
else
  myfailed('Nothing to export.',DATA.GUI.Segment);
end

%-----------------------------------
function measureshapeexport_Callback %#ok<DEFNU>
%-----------------------------------
%Export measurement shape
global DATA SET NO
%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

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
  end
  cell2clipboard(c);
else
  myfailed('Nothing to export.',DATA.GUI.Segment);
end


%--------------------------
function parallelout = updateparallelsets(no)
%--------------------------
% Chages CurrentSlice in SETs that are parallel to SET(NO) to the slice
% closest to SET(NO).CurrentSlice
%
% Marten Larsson, June, 3, 2009

global DATA SET 

% Get open SETs panels but exclude no
openpanels = DATA.ViewPanels;
openpanels = unique(openpanels);
openpanels = openpanels(openpanels~=no&openpanels~=0);

% Find SETs that are parallel to SET(no), but not a linked flow.
if ~isempty(SET(no).Flow)
  flownos = [...
    SET(no).Flow.MagnitudeNo ...
    SET(no).Flow.PhaseNo ...
    SET(no).Flow.PhaseX ...
    SET(no).Flow.PhaseY ...
    SET(no).Flow.Angio ...
    SET(no).Flow.VelMag];
else
  flownos = [];
end
parallel  = [];
for i = openpanels
  sameview = orientationcomparison(no,i);
  linkedflow = any(flownos == i);  % SET(i) is linked to SET(no) via SET(no).Flow.
                                   % This case is handled by updateoneim().
  if sameview && ~linkedflow
    parallel = [parallel i]; %#ok<AGROW>
  end
end

if not(isempty(parallel))
  % Calculate new CurrentSlide for parallel SETs
  %viewdir = cross(SET(no).ImageOrientation(1:3),SET(no).ImageOrientation(4:6))';
  %currzdist = (SET(no).CurrentSlice-1)*(SET(no).SliceThickness+SET(no).SliceGap);
  pos = calcfunctions('xyz2rlapfh',no,SET(no).XSize/2,SET(no).YSize/2,SET(no).CurrentSlice); %Center of image and current slice
  for i = 1:length(parallel)
    
    %--- Old code
    % Calculate closest slice in parallel SETs and change their
    % CurrentSlice
    %slicethickness = SET(parallel(i)).SliceThickness + SET(parallel(i)).SliceGap;
    %zdistances = 0:slicethickness:(slicethickness * (SET(parallel(i)).ZSize-1));
    %zdiff = SET(parallel(i)).ImagePosition - SET(no).ImagePosition;
    %zdiff = dot(zdiff,viewdir);
    %zdistances =  zdistances - zdiff;
    %[~,slice] = min(abs(zdistances - currzdist));
   
    
    %--- New code
    [outpos] = calcfunctions('rlapfh2xyz',parallel(i),pos(1),pos(2),pos(3));
    slice = min(max(1,round(outpos(3))),SET(parallel(i)).ZSize);
    
    %--- Store 
    SET(parallel(i)).CurrentSlice = slice;
    SET(parallel(i)).StartSlice = slice;
    SET(parallel(i)).EndSlice = slice;    
    
  end
end

if nargout == 1
    parallelout = parallel;
end

%-----------------------------
  function do = doatrialscar(no)
%-----------------------------
%Helper function to check if user input is for atrial scar or lv scar
%(default)

global SET

do = false;

if isempty(SET(no).RVEndoX)
  return;
end

if ~isnan(SET(no).RVEndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice))
  if isempty(SET(no).EpiX) || isnan(SET(no).EpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice))
    %RV exist and not LV Epi, then go for atrial scar
    do = true;
  end
end
    
%---------------------------------------------
  function updatervstack
%---------------------------------------------
%Updated the definition of RV stack number

global DATA NO

if isempty(DATA.RVNO)
  DATA.RVNO = NO;
  if isfield(DATA.Handles,'rvstackpushbutton')
    set(DATA.Handles.rvstackpushbutton,'String', sprintf('Stack #%d',DATA.RVNO))
  end
end
segment('updatevolume');
% DATA.updateaxestables('volume')

%---------------------------------------------
  function updatelvstack
%---------------------------------------------
%Updated the definition of LV stack number

global DATA NO

if isempty(DATA.LVNO) && isfield(DATA.Handles,'lvstackpushbutton')
  DATA.LVNO = NO;
  set(DATA.Handles.lvstackpushbutton,'String', sprintf('Stack #%d',DATA.LVNO))
end
segment('updatevolume');
% DATA.updateaxestables('volume')  

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




%------------------
function rotatetemp(alpha) %#ok<DEFNU>
%------------------
%This function should late be replaced with qa tool that rotates objects,
%just as scale and move does. 

global NO SET

if nargin==0
  alpha = 10;
end

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

viewfunctions('setview'); 
% drawfunctions('drawallslices');%was drawsliceno but changed in order to update intersectionpoints

%------------------------------
function timebar_Motion 
%------------------------------
%Update current timeframe when user clicked in time bar graph.
global DATA SET NO

[x] = mygetcurrentpoint(DATA.Handles.timebaraxes);
x = x/1000;

%Find timeframe
SET(NO).CurrentTimeFrame = 1+round(x/SET(NO).TIncr);
SET(NO).CurrentTimeFrame = min(max(1,SET(NO).CurrentTimeFrame),SET(NO).TSize);

panels = find(ismember(DATA.ViewPanels,SET(NO).Linked));
for p = panels
    DATA.ViewIM{p} = [];
    drawfunctions('drawpanel',p)
end

%update timebars
viewfunctions('updatetimebars')
drawnow;

%----------------------------------
function timebar_Buttondown %#ok<DEFNU>
%----------------------------------
%Button down function for the time bar in the time bar graph
global DATA

if DATA.Interactionlock
  return;
end

%stateandicon=segment('iconson',{'hidescar','hidemar','hideall','play'});
%set(DATA.fig,'WindowButtonMotionFcn',@ segment_main.timebar_Motion(stateandicon));
set(DATA.fig,'WindowButtonMotionFcn',...
  sprintf('%s(''timebar_Motion'')',mfilename));
set(DATA.fig,'WindowButtonUpFcn',...
  sprintf('%s(''buttonup_Callback'')','buttonupfunctions'));

timebar_Motion


%--------------------------------
function fasterframerate_Callback %#ok<DEFNU>
%--------------------------------
%Makes movie play faster
global DATA SET

%First we need to pause the running by setting DATA.Run = 0; Then do a
%pause in the software so playing is taken down safely
wason = 0;
if DATA.Run
    wason = 1;
    undent(DATA.Handles.permanenticonholder,'play',0)
    %pause(0.1)
end

ind = find(cat(1,SET(:).TSize)>1);
if isempty(ind)
  ind = 1;
end

for no=ind'
  SET(no).BeatTime = SET(no).BeatTime*0.9;
end

% if wason
%     indent(DATA.Handles.permanenticonholder,'play',1)    
% end
%--------------------------------
function slowerframerate_Callback %#ok<DEFNU>
%--------------------------------
%Makes movie play slower
global DATA SET

wason = 0;
if DATA.Run
    wason = 1;
    undent(DATA.Handles.permanenticonholder,'play',0)
    %pause(0.1)
end

%If running make sure always update the first one.
ind = find(cat(1,SET(:).TSize)>1);
if isempty(ind)
  ind = 1;
end

for no=ind'
  SET(no).BeatTime = SET(no).BeatTime/0.9;
end

% if wason
%     indent(DATA.Handles.permanenticonholder,'play',1)    
% end

%--------------
function ctrlc %#ok<DEFNU>
%--------------
%function to handle ctrl-c keypress, which is disabled
global DATA
mywarning('Ctrl-C is disabled',DATA.GUI.Segment);

%--------------------------------------------------------
function [varargout] = cell2clipboard(outdata,writetofile) 
%--------------------------------------------------------
%Converts a cell to a string that is output to clipboard.
%If more than 8000 cells are written then an .xls file is
%written instead. Note that this used active-X on Windows and requires
%Excel to be installed on the computer.

if nargout>0
  varargout = cell(1,nargout);
end

if nargin<2
  writetofile = false;
end

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
          end
        case {'uint8','int8','uint32'}
          stri = [stri sprintf('%d\t',temp)]; %#ok<AGROW>
        case {'double','single'}
          if length(temp)>1
            myfailed('Vector not allowed in one cell.');
            return;
          end
          if isempty(temp)||isnan(temp)
            stri = [stri sprintf('\t')]; %#ok<AGROW>
          else
            if isempty(temp)
              stri = [stri sprintf('\t')]; %#ok<AGROW>
            else
              stri = [stri sprintf('%f\t',temp)]; %#ok<AGROW>
            end
          end
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
          end
        otherwise
          myfailed('Unknown object type when exporting.');
          disp(class(temp))
          return;
      end
    end
    stri = [stri newline]; %#ok<AGROW>
  end
  if nargout==0
    clipboard('copy',stri);
    mymsgbox('Data copied to clipboard.','Done!');
  else
    varargout{1} = stri;
  end
else
  [filename, pathname] = myuiputfile('*.xlsx',dprintf('Save as Excel file.'));
  if isequal(filename,0) || isequal(pathname,0)
    disp('Aborted');
    return;
  else
    [sucess,message] = xlswrite(fullfile(pathname,filename),outdata);
    if sucess
      mymsgbox('File successfully written.','Done!');
    else
      myfailed(dprintf('Error:%s',message.message));
      return;
    end
  end
end
flushlog;

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
end

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
    end
    y=SET(NO).EndoY;
    x=SET(NO).EndoX;
    name='endocardium';
  case 2
    %--- Import epicardium
    if isempty(SET(NO).EpiX)
      myfailed('No LV epicardium available.',DATA.GUI.Segment);
      return;
    end
    y=SET(NO).EpiY;
    x=SET(NO).EpiX;
    name='epicardium';
  case 3
    %--- Import RVendocardium
    if isempty(SET(NO).RVEndoX)
      myfailed('No RV endocardium available.',DATA.GUI.Segment);
      return;
    end
    y=SET(NO).RVEndoY;
    x=SET(NO).RVEndoX;
    name='RVendocardium';
  case 4
    %--- Import RVepicardium
    if isempty(SET(NO).RVEpiX)
      myfailed('No RV epicardium available.',DATA.GUI.Segment);
      return;
    end
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

eval_str=myinputdlg('command:','Execute matlab command');
if ~isempty(eval_str)
  try
    eval(cell2mat(eval_str));
  catch me
    mydispexception(me);
  end
end

segmenthelp('openthislogfile_Callback');

%---------------------------------
function changewheel_Callback(h,e) %#ok<DEFNU>
%---------------------------------
%scrollwheel with modifer.
%Tab are not included but you can if you like
global SET DATA

if not(DATA.DataLoaded)
  return;
end

if strcmp(DATA.CurrentTheme,'3dp') && any(strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},{'trans3DP','sag3DP','cor3DP','speedim','viewport'}))
  segment3dp.tools('scrollwheelfcn3dp',h,e)
  return;
end

wheelup = e.VerticalScrollCount>0;

modifier = char(get(gcbf,'currentmodifier'));

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

%Calculate no
no = DATA.ViewPanels(DATA.CurrentPanel);

switch mod_str
  case 'shift'
    if SET(no).TSize>1
      if wheelup
        viewfunctions('switchtimeframe',1);
      else
        viewfunctions('switchtimeframe',-1);
      end
    end
%   case 'alt'
%     if wheelup,
%       zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1.2);
%     else
%       zoomhelper(DATA.Handles.imageaxes(DATA.CurrentPanel),1/1.2);
%     end
  case 'control'
    if wheelup
      thumbnailslider_Callback('up');
    else
      thumbnailslider_Callback('down');      
    end
  case 'shift-alt'
  case 'control-alt'
  case 'shift-control'
    if wheelup
      %contrast_Callback('wheelUpBrightness',DATA.CurrentPanel);
      SET(no).IntensityMapping.Brightness = SET(no).IntensityMapping.Brightness+0.1;
      DATA.ViewIM{DATA.CurrentPanel} = [];
      drawfunctions('drawpanel',DATA.CurrentPanel);
      drawfunctions('drawthumbnails',1);
    else
      SET(no).IntensityMapping.Brightness = SET(no).IntensityMapping.Brightness-0.1;
      DATA.ViewIM{DATA.CurrentPanel} = [];
      drawfunctions('drawpanel',DATA.CurrentPanel);
      drawfunctions('drawthumbnails',1); 
    end
  case 'shift-control-alt'
  otherwise %no modifier
    if SET(no).ZSize>1
      
      if wheelup
        viewfunctions('switchslice',1);
      else
        viewfunctions('switchslice',-1);
      end
      
    end
end

%---------------------------------
function recursekeyreleasefcn(h,fcn)
%---------------------------------
%Helper function to create callbacks to keyrelease function.
if nargin<2
  fcn = @(x,y)segment('keyreleased',x,y);
end

allowedtypes = {'figure', 'uicontrol', 'uipushtool', 'uitable', 'uitoolbar'};

if any(contains(allowedtypes,h.Type))
  set(h,'KeyReleaseFcn',fcn);
end
if any(contains({'figure','uipanel','uibuttongroup'},h.Type))
  children = get(h,'children');
  for loop=1:length(children)
    recursekeyreleasefcn(children(loop),fcn);
  end
end

%---------------------------------
function recursekeypressfcn(h,fcn)
%---------------------------------
%Helper function to create callbacks to keypressed function.
if nargin<2
  fcn = @(x,y)segment('keypressed',x,y);
end

allowedtypes = {'figure', 'uicontrol', 'uipushtool', 'uitable', 'uitoolbar'};

%Start with the current handle.
if any(contains(allowedtypes,h.Type))
  set(h,'keypressfcn',fcn);
end

if any(contains({'figure','uipanel','uibuttongroup'},h.Type))
  children = get(h,'children');
  for loop=1:length(children)
    recursekeypressfcn(children(loop),fcn);
  end
end

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
end

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
  end

  firstthumbnail=slidervalue;
  lastthumbnail=max(min(firstthumbnail+DATA.Pref.NumberVisibleThumbnails-1,length(SET)),1);
else
  switch whattodo
    case 'down'
      firstthumbnail=max(DATA.VisibleThumbnails(1),1);
      lastthumbnail=min(firstthumbnail+DATA.Pref.NumberVisibleThumbnails-1,length(SET));
      slidervalue = firstthumbnail;
    case 'up'
      lastthumbnail=min(DATA.VisibleThumbnails(end),length(SET));
      firstthumbnail=max(lastthumbnail-DATA.Pref.NumberVisibleThumbnails+1,1);
      slidervalue = firstthumbnail;      
  end
end

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
  sliderstep = [1/(length(SET)-DATA.Pref.NumberVisibleThumbnails),2/(length(SET)-DATA.Pref.NumberVisibleThumbnails)];
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
GUIPositions = DATA.GUIPositions; %Saved to file
try
  save([pathname filesep '.segment_guipositions.mat'],'GUIPositions', DATA.Pref.SaveVersion);
catch %#ok<CTCH>
  myfailed('Could not save GUI positions. Write permission? Disk full?');
  return;
end

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
GUIPositions(1).FileName = 'segment.fig';
GUIPositions(1).Position = [0.05 0.05 0.9 0.85];%normalized position

DATA.GUIPositions = GUIPositions;

guinames = fieldnames(DATA.GUI);
for loop = 1:length(guinames)
  gui = getfield(DATA.GUI,guinames{loop}); %#ok<GFLD>
  if ~isempty(gui)
    try
      setguiposition(gui);
      mainresize_Callback;
    catch me
      disp(sprintf('Could not reset GUI %s',guinames{loop})); %#ok<DSPS>
      mydispexception(me);
    end
  end
end

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

slicestoinclude = find(findfunctions('findslicewithendo',no)+findfunctions('findslicewithepi',no)+findfunctions('findslicewithrvendo',no))';
if isempty(slicestoinclude)
  slicestoinclude = SET(no).CurrentSlice;
end
slicestoinclude=slicestoinclude(1):slicestoinclude(end);
%this is a test
if SET(no).StartSlice<slicestoinclude(1)
  slicestoinclude=[SET(no).StartSlice:slicestoinclude(1)-1,slicestoinclude];
else
  if min(slicestoinclude) > 1
    slicestoinclude = [min(slicestoinclude)-1 slicestoinclude];
  end
end

if SET(no).EndSlice>slicestoinclude(end)
  slicestoinclude=[slicestoinclude,slicestoinclude(end)+1:SET(no).EndSlice];
else
  if max(slicestoinclude) < SET(no).ZSize
    slicestoinclude = [slicestoinclude max(slicestoinclude)+1];
  end
end

%-------------------------------------------------------
 function maxlvdiameter_Callback(type)
%-------------------------------------------------------
     %Get maximum LV diameter.
     global SET DATA
     maxtmp=0;
     maxdia=0;
     Points=[];
     maxtf=0;

     Zslice=[];
     ind=0;
     csax=[];
     csax=findfunctions('findctsaxwithsegmentation','Endo');
     
     %Do check if possible to get LV diameter
     if isempty(csax)
         mywarning(dprintf('Unable to get max diameter since no short axis with Endo segmentation is available.'))
        return
     end
 
     for i=1:length(SET(csax).Measure)
         if strcmp(SET(csax).Measure(i).LongName,'Automatic LV SAX max diameter')
             ind=i;
         end
     end
     
     if ind == 0
         ind=length(SET(csax).Measure)+1;
     end
     for tf = 1:SET(csax).TSize
         [maxtmp, Points, Zslice] = calcfunctions('maxsaxdiameter',csax,tf,'LV');

         if maxtmp>maxdia
             maxdia=maxtmp;
             maxtf=tf;
             SET(csax).Measure(ind).X=Points(:,1);
             SET(csax).Measure(ind).Y=Points(:,2);
             SET(csax).Measure(ind).Z=[Zslice;Zslice];
             SET(csax).Measure(ind).Length=maxdia;
             SET(csax).Measure(ind).Name='LV SAX max';
             SET(csax).Measure(ind).LongName='Automatic LV SAX max diameter';
             SET(csax).Measure(ind).T=tf;
         end
     end

 %NO=csax;
 %SET(csax).CurrentSlice = SET(csax).Measure(end).Z;
 SET(csax).CurrentTimeFrame = maxtf;
%drawfunctions('drawimageview', csax, [1,1],{'one'})
 %drawfunctions('drawmeasures',csax,DATA.CurrentPanel)
 %DATA.Handles.measureline{DATA.CurrentPanel}{length(SET(csax).Measure)+1} = plot(...
 %      DATA.Handles.imageaxes(DATA.CurrentPanel),[],[],DATA.GUISettings.MeasureLineSpec);
 viewfunctions('setview',1,1,csax,{'one'}); %drawfunctions('drawimageview', csax, [1,1],{'one'})
 %switchtoslice(SET(csax).Measure(ind).Z(1))%drawfunctions('drawall',1,1)
 sliceincr = SET(csax).Measure(ind).Z(1) - SET(csax).CurrentSlice;
 viewfunctions('switchslice',sliceincr,1);
 
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
  end

end



%-------------------------------------------------------
 function maxrvdiameter_Callback(type)
%-------------------------------------------------------
 %Get maximum RV diameter on SAX images.
     
 global SET DATA
 maxtmp=0;
 maxdia=0;
 Points=[];
 maxtf=0;
 %Zslicetmp=0;
 Zslice=[];
 ind=0;
 csax=[];
 LVsax=findfunctions('findctsaxwithsegmentation','Endo');
 csax=findfunctions('findctsaxwithsegmentation','RVEndo');

 %Do check if possible to get RV diameter

  %   RV/LV check
if isempty(csax)&&isempty(LVsax)
     mywarning(dprintf('Unable to get max RV diameter since no short axis with RV and LV segmentation is available.'))
     return
end
 
 %   RV/LV Epi check
if isempty(csax)&& (isempty(SET(LVsax).EpiX))
     mywarning(dprintf('Unable to get max RV diameter since no short axis with RV and LV Epi segmentation is available.'))
     return
 end
  %   RV check
 if isempty(csax)
     mywarning(dprintf('Unable to get max RV diameter since no short axis with RV Endo segmentation is available.'))
     return
 end
 
  %   LV Epi check
 if (isempty(SET(csax).EpiX))
     mywarning(dprintf('Unable to get max RV diameter since LV Epi segmentation is not available.'))
     return
 end

 
 for i=1:length(SET(csax).Measure)
     if strcmp(SET(csax).Measure(i).LongName,'Automatic RV SAX max diameter')
         ind=i;
     end
 end
 
 if ind == 0
     ind=length(SET(csax).Measure)+1;
 end
 for tf = 1:SET(csax).TSize
     [maxtmp, Points, Zslice] = calcfunctions('maxsaxdiameter',csax,tf,'RV');
     if maxtmp>maxdia
         maxdia=maxtmp;
         maxtf=tf;
         SET(csax).Measure(ind).X=Points(:,1);
         SET(csax).Measure(ind).Y=Points(:,2);
         SET(csax).Measure(ind).Z=[Zslice;Zslice];
         SET(csax).Measure(ind).Length=maxdia;
         SET(csax).Measure(ind).Name='RV SAX max';
         SET(csax).Measure(ind).LongName='Automatic RV SAX max diameter';
         SET(csax).Measure(ind).T=tf;
     end
 end
 % RV and LV slice check
 if maxdia==0
      mywarning(dprintf('Only including slices with both RV Endo and LV Epi segmentation'))
     return
 end

 SET(csax).CurrentTimeFrame = maxtf;
 viewfunctions('setview',1,1,csax,{'one'}); %drawfunctions('drawimageview', csax, [1,1],{'one'})
 %switchtoslice(SET(csax).Measure(ind).Z(1))%drawfunctions('drawall',1,1)
 sliceincr = SET(csax).Measure(ind).Z(1) - SET(csax).CurrentSlice;
 viewfunctions('switchslice',sliceincr,1)