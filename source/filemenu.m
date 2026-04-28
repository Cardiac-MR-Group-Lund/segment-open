function varargout = filemenu(varargin)
% FILEMENU
% File menu callbacks

% Nisse Lundahl

%#ok<*GVMIS> 
%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%--------------------------------------
function loadednext = loadnext_Callback %#ok<*DEFNU>
%--------------------------------------
%Load next .mat file in the current data folder.
global DATA SET NO 

if ~DATA.DataLoaded
  loadednext = 0;
  return;
end

if length(SET)<1
  loadednext = 0;
  return;
end

%Get filename from first image stack.
filename = SET(1).FileName;

[pathname,filename,~] = fileparts(filename);

f = dir([pathname filesep '*.mat']);

pos = [];
for loop=1:length(f)
  if isequal(f(loop).name,[filename '.mat'])
    pos = loop;
  end
end

if isempty(pos)
  myfailed('Could not find the current file, removed?');
  loadednext = 0;
  return;
end

if (pos+1)<=length(f)
  
  %Ask to quite nicely
  DATA.filecloseall_Callback;
  
  if ~isempty(SET)
    myfailed('Loading next aborted.');
    loadednext = 0;
    return;
  end
  
  try
    %Load
    SET = []; %Help Matlab memory garbage collector
    load([pathname filesep f(pos+1).name],'-mat');

    %Assign
    SET = setstruct;
    clear setstruct;
    
    NO = 1; 
    
    %Call to intialize all variables correcly after loaded data.
    DATA.Preview.PreviewFile = [pathname filesep f(pos+1).name];
    
    openfile('setupstacksfrommat',NO);
    segment('renderstacksfrommat');
 
     if SET(1).Autoloader.Autoloaded && SET(1).Autoloader.Autoanalysed %file was created with autoloader and analysed
      DATA.CurrentTheme = 'review';
     end

    %Does all graphical updates
    if not(isempty(DATA.ViewMatrix))
      DATA.init_graphics;
    end
    
    DATA.NeedToSave = false;
   
    %Show warnings from Autoloader
    if SET(1).Autoloader.Autoloaded %file was created with autoloader
      set(DATA.Handles.autoloaderwarningsmenu, 'Enable', 'on');
      if ~SET(1).Autoloader.Approved
         if ~isempty(SET(1).Autoloader.Warnings)
           warningfunctions('showautoloaderwarnings', 'AI AutoMate', SET(1).Autoloader.Autoanalysed);
         end
      end
    else
      set(DATA.Handles.autoloaderwarningsmenu, 'Enable', 'off');
    end
  catch me
    mydispexception(me);
  end
  loadednext = 1;
else
  myfailed('Last file');
  loadednext = 0;
end

%---------------------------------------------
function savetopatientdatabaseandpacs_Callback
%---------------------------------------------
%Callback to save image stacks to patientdatabase and PACS in one go.

savetopatientdatabase_Callback;
savetopacs_Callback;

%--------------------------------------
function savetopatientdatabase_Callback(saveinbackground,autocomment)
%--------------------------------------
%Callback to save image stacks to patientdatabase. Uses functions in
%patientdatabase.

%If the option 'saveinbackground' is set to true, then the comment
%'autocomment' is automatically added to the study in the patient database.
%This is used when exporting to REDCap database for example.

global DATA SET NO

if nargin < 1
  saveinbackground = false;
end

if nargin < 2
  autocomment = [];
end

tools('closeallongoinginterpolations');
[pathname,filename,ext] = fileparts(SET(NO).FileName);

if ~saveinbackground
  if DATA.Pref.AskAddComment
    try
      comment = patientdatabase('addcomment_init');
    catch me
      mydispexception(me);
    end
  else
    comment = [];
  end
else
  comment = autocomment;
end

try
  patientdatabase('savetodatabase',pathname,filename,ext,comment,saveinbackground);
catch me
  mydispexception(me);
end

%------------------------
function saveall_Callback
%------------------------
%Saves all image stacks to one .mat file. Calls GUI method 
%filesaveallas_Callback which is the workhorse when saving image stacks.

global DATA SET

tools('closeallongoinginterpolations');
filenames = unique({SET.FileName});
filenames = filenames(~cellfun('isempty', filenames));
if numel(filenames) > 1
  filenames = strrep(filenames,'\','\\');
  i = mymenu('Multiple file names were found in data. Select where to save', ...
    filenames);
  if i > 0
    fname = filenames{i};
  else
    myfailed('Saving aborted');
    return
  end
else
  if ~isempty(filenames)
    fname = filenames{1};
  else
    fname = filenames;
  end
end
[pathname,filename,ext] = fileparts(fname);

%--- Normal research mode
if isequal(ext,'.mat')
  DATA.filesaveallas_Callback(pathname,[filename ext]);
else
  %Call Save as instead
  DATA.filesaveallas_Callback;
end

DATA.LastSaved = now;
flushlog;

%--------------------------------------------------------------------------------
function fail = saveallas_helper(pathname,filename,topatientdatabase,silent,fast)
%--------------------------------------------------------------------------------
%Save all image stacks to the file specified. It also stores current view
%and modes etc.

global DATA SET NO
    
if nargin < 3
  topatientdatabase = false; %Default save to disk not patientdatabase
  %Setting this to true only changes the final saved display message.
end

if nargin < 4
  silent = DATA.Silent;
end

if nargin < 5
  fast = false;
end

fail = 0; %optimistic

%Check if there are different patientnames, only performed if not silent
if ~silent
  differentnames = false;
  differencestring = '';
  patientname = SET(1).PatientInfo.Name;
  namestri = dprintf('Name');
  stackstri = dprintf('Stack');
  for loop = 2:length(SET)
    thisname = SET(loop).PatientInfo.Name;
    if ~isequal(patientname,thisname)
      differentnames = true;
      differencestring = sprintf('%s%s: %d, %s: "%s"\n',differencestring,stackstri,loop,namestri,thisname);
    end
  end
  if differentnames
    originalname = sprintf('Stack: %d, %s:"%s"\n',1,namestri,patientname);
    messagestri = dprintf('Different patient names detected.');
    messagestri = [messagestri newline originalname differencestring];
    messagestri = [messagestri newline dprintf('It is recommended to abort and update patient info.') ' '];
    messagestri = [messagestri dprintf('Continue anyway') '?'];
    if ~yesno(messagestri)
      fail = true;
      return
    end
  end
end

%--- Get filename and pathname
if nargin==0
  temp = DATA.SegmentFolder;
  if exist(DATA.Pref.exportpath,'dir')
    cd(DATA.Pref.exportpath);
  else
    mydisp('Warning: Export path does not exist, please check preferences.');
  end
  
  filename = DATA.generatesavefilename;
  
  [filename,pathname] = myuiputfile(filename,dprintf('Save all image stacks to one file'));
  cd(temp);
else
  %Take from input two argument
  if nargin<2
    myfailed('Expected two input arguments.',DATA.GUI.Segment);
    return;
  end
end

if isequal(filename,0)
%   myfailed('Operation cancelled.',DATA.GUI.Segment);
  fail = 1;
  return;
end

if length(filename)>4
  if ~isequal(filename(end-3:end),'.mat')
    filename = [filename '.mat'];
  end
else
  filename = [filename '.mat'];
end

myworkon(DATA.fig);

%Store info about current user (e.g. username, computer name, date)
%according to user preferences
if DATA.Pref.SaveUserInfo
  SET(NO).LastUserInfo.Username = getenv('username');
  SET(NO).LastUserInfo.ComputerName = getenv('computername');
  SET(NO).LastUserInfo.ModificationDate = datestr(now,'yyyy-mm-dd HH:MM');
else
  SET(NO).LastUserInfo = '';
end

%Store filename
for loop=1:length(SET)
  SET(loop).FileName = [pathname filesep filename];
end

%Store current view
SET(1).View = [];
SET(1).View.ViewPanels = DATA.ViewPanels;
SET(1).View.ViewPanelsType = DATA.ViewPanelsType;
SET(1).View.ViewPanelsMatrix = DATA.ViewPanelsMatrix;
SET(1).View.ViewMatrix = DATA.ViewMatrix;
SET(1).View.ThisFrameOnly = DATA.ThisFrameOnly;
SET(1).View.CurrentPanel = DATA.CurrentPanel;
SET(1).View.CurrentTheme = DATA.CurrentTheme;
SET(1).View.RelevantMode = DATA.ShowRelevantStacksOnly;
SET(1).View.CurrentTool = DATA.CurrentTool;
SET(1).View.LVNO = DATA.LVNO;
SET(1).View.RVNO = DATA.RVNO;
SET(1).View.FlowNO = DATA.FlowNO;
SET(1).View.FlowROI = DATA.FlowROI;

%get preview
if DATA.ShowRelevantStacksOnly
  % create preview of all stacks
  preview = calcfunctions('calcdatasetpreview_helper');
else
  preview = DATA.DATASETPREVIEW; %saved to file
end

%get info 
info = SET(NO).PatientInfo;
%correct study date 
info.AcquisitionDate = helperfunctions('getmostrecentstudydate');
%fake data
info.NFrames = 0;
info.NumSlices = 0;
info.ResolutionX = 0;
info.ResolutionY = 0;
info.SliceThickness = 0;
info.SliceGap = 0;
info.TIncr = 0;
info.EchoTime = 0;
info.FlipAngle = 0;
info.AccessionNumber = '';
info.StudyUID = '';
info.StudyID = '';
info.NumberOfAverages = 0;
info.RepetitionTime = 0;
info.InversionTime = 0;
info.TDelay = 0;
info.VENC = 0;
info.Scanner = '';
info.ImagingTechnique = '';
info.ImageType = '';
info.ImageViewPlane = '';
info.IntensityScaling = 1;
info.IntensityOffset = 0;
info.MultiDataSet = true;
info.Modality = SET(1).Modality;
if ~isempty(SET(1).Autoloader)
  %only AutoMate format in patient database if it is autoloaded and not
  %reviewed
  info.AutoMateformat = SET(1).Autoloader.Autoloaded && ~SET(1).Autoloader.Approved;
else
  info.AutoMateformat = false;
end

im = [];  %saved to file

if DATA.issegment3dp()
  myworkon(DATA.fig)
  for loop = 1:length(SET)

    %Compress infinite number of objects and force them to be compressed
    SET(loop).LevelSet.Object.docompress(inf,true); 

    %Restore image enhancement
    if isfield(SET(loop),'IMBackup') && (~isempty(SET(loop).IMBackup))
      SET(loop).IM = SET(loop).IMBackup;
      SET(loop).IMBackup = [];
    end

  end
end

setstruct = SET;  %saved to file

lastwarn('');
warnchk1 = 'not saved';
warnchk2 = 'corrupt';

retried = false;
try

  %Check size
  w = whos('SET');

  myworkon(DATA.fig);
  drawnow;
  threshold = 1.5; %smaller than 1.5 GB

  if ((w.bytes/1e9) < threshold) && (~fast) 
    %File is smaller than threshold => default and it is not fast
    logdisp('Normal save')
    fail = true; %#ok<NASGU> 
    tic
    save(fullfile(pathname,filename),'preview','info','im','setstruct', DATA.Pref.SaveVersion);
    t = toc;
    logdisp(sprintf('Saving took %0.5g s',t));
    fail = false;
  else
    %Save with nocompression
    logdisp('Save with compression')
    retried = true; %It will be no use re-saving it
    fail = true; %#ok<NASGU> %pesimistic
    tic
    save(fullfile(pathname,filename),'preview','info','im','setstruct', '-v7.3','-nocompression');
    t = toc;
    logdisp(sprintf('Saving took %0.5g s',t));
    fail = false;
  end

  warnmsg = lastwarn;
  lastwarn('');
  if contains(warnmsg,warnchk1) || contains(warnmsg,warnchk2)

    if ~retried
      mywarning('First save attempt failed, retrying using -v7.3 version.');
      logdisp('String contained corrupt or not saved')
      fail = true; %#ok<NASGU> 
      save(fullfile(pathname,filename),'preview','info','im','setstruct','-v7.3','-nocompression'); %added nocompression 2024-09-26 /EiH
      fail = false;
      retried = true;
      warnmsg = lastwarn;
      lastwarn('');
      if ~isempty(warnmsg)
        myfailed('Dataset cannot be saved');
        fail = 1;
        myworkoff(DATA.fig);
        return;
      end

    else
      mywarning('Problem saving');
    end

  end
  myworkoff(DATA.fig);
catch me %#ok<CTCH>
  mydispexception(me);

  if ~retried
    mywarning('First save attempt failed, retrying using -v7.3 version.');
    lastwarn('');
    fail = false;
    try
      save(fullfile(pathname,filename),'preview','info','im','setstruct','-v7.3','-nocompression'); %added nocompression 2024-09-26 /EiH
    catch me
      mydispexception(me)
      fail = true;
    end
    retried = true;
    warnmsg = lastwarn;
    lastwarn('');
    if ~isempty(warnmsg)
      fail = true;
    end
  else
    fail = 1;
  end

  if fail || not(retried)
    if contains(me.message,warnchk1) || contains(me.message,warnchk2)
      stri = dprintf('Save failed. Data may be corrupted.');
      if DATA.issegment3dp
        stri = [stri ' ' dprintf('Try to delete unecessary objects or crop image volume.')];
      end
      stri = [stri ' ' dprintf('Retry saving under a different name.')];
      myfailed(stri);
    else
      myfailed('Could not save data. Write permission? Disk full?',DATA.GUI.Segment);
    end
  end
  myworkoff(DATA.fig);
  return;
end

if not(silent)
  if topatientdatabase
    mymsgbox('Image stacks saved to patient database.');
    logdisp('Image stacks saved to patient database.');
  else
    filename = SET(NO).FileName;
    if length(SET)>1
      stri = dprintf('All image stacks and segmentation stored to %s.',filename);
    else
      stri = dprintf('Image stack and segmentation stored to %s.',filename);
    end
    mymsgbox(stri,'Save successful.',DATA.GUI.Segment);
    logdisp('Save successful');
  end
end
if ~isempty(SET(1).FileName)
  str = sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),DATA.LoggedUser,sprintf('Saved file %s',SET(1).FileName),SET(1).PatientInfo.ID);
else
  str = sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),DATA.LoggedUser,sprintf('Saved file %s',SET(1).PathName),SET(1).PatientInfo.ID);
end

DATA.adduserevent(str);
%DATA.adduserevent([' Saved file: ', SET(1).FileName])
%DATA.adduserevent(['Time:' datestr(now,'yyyymmddHHMMSS')])
%DATA.adduserevent([' Patient name: ', SET(1).PatientInfo.Name])
%DATA.adduserevent([' Patient ID: ', SET(1).PatientInfo.ID])
DATA.LastSaved = now;
DATA.NeedToSave = false;
% set(DATA.Handles.filesaveicon,'enable','off');
DATA.updatetitle;
flushlog;

%---------------------------
function savetoopenapps_Callback
%---------------------------
%Send image stacks to the folder provide by command line parameter
% save into the siemens catalogue
  try
    applysiemensoptions('savetosiemens_Callback')
  catch e
    myfailed(e.message);
  end

%---------------------------
function savetopacs_Callback
%---------------------------
%Send image stacks to PACS. This function should display a list of
%available PACS (.con files) and when user has selected store files on disk
%temporarily and then send the files to the PACS.

global SET

% save to pacs
try
 successful = pacsaccess('savetopacs');
 if successful == 1
    mymsgbox('Image stacks saved to PACS.')
 end
catch e
  myfailed(e.message);
end

%Remove temp file used to export to PACS
pacsaccess('removetempsegdicom');

%----------------------------
function savecurrent_Callback 
%----------------------------
%Save current image set to file. Note that this is the old Segment file
%format and this fcn may soon be depreciated.
global DATA SET NO

if not(DATA.DataLoaded)
  myfailed('No data to save.',DATA.GUI.Segment);
  return;
end

temp = DATA.SegmentFolder;
if exist(DATA.Pref.exportpath,'dir')
  cd(DATA.Pref.exportpath);
else
  mydisp('Warning: Export path does not exist, please check preferences.');
end

[filename,pathname] = myuiputfile('*.mat',dprintf('Save current image stack'));
cd(temp);

if isequal(filename,0)
%   myfailed('Operation cancelled.',DATA.GUI.Segment);
  return;
end

%IM and preview are saved to disk
preview = SET(NO).IM(:,:,round(SET(NO).TSize/2),round(SET(NO).ZSize/2)); 
info = SET(NO).PatientInfo;
info.NFrames = SET(NO).TSize;
info.NumSlices = SET(NO).ZSize;
info.ResolutionX = SET(NO).ResolutionX;
info.ResolutionY = SET(NO).ResolutionY;
info.SliceThickness = SET(NO).SliceThickness;
info.SliceGap = SET(NO).SliceGap;
info.TIncr = SET(NO).TIncr;
info.EchoTime = SET(NO).EchoTime;
info.RepetitionTime = SET(NO).RepetitionTime;
info.InversionTime = SET(NO).InversionTime;
info.FlipAngle = SET(NO).FlipAngle;
info.AccessionNumber = SET(NO).AccessionNumber;
info.StudyUID = SET(NO).StudyUID;
info.StudyID = SET(NO).StudyID;
info.NumberOfAverages = SET(NO).NumberOfAverages;
info.TDelay = SET(NO).TDelay;
info.VENC = SET(NO).VENC;
info.Scanner = SET(NO).Scanner;
info.ImagingTechnique = SET(NO).ImagingTechnique;
info.ImageType = SET(NO).ImageType;
info.ImageViewPlane = SET(NO).ImageViewPlane;
info.IntensityScaling = SET(NO).IntensityScaling;
info.IntensityOffset = SET(NO).IntensityOffset;
info.Modality = SET(NO).Modality;
im = calcfunctions('calctruedata',SET(NO).IM,NO); 
setstruct = SET(SET(NO).Linked);
linkies = setstruct(1).Linked;
newlinkies = 1:length(linkies);
dadno = find(cellfun(@isempty,{setstruct.Parent}));
flowstruct = setstruct(dadno).Flow;
if ~isempty(flowstruct)
  fnames = fieldnames(flowstruct);
  for i = 1:length(fnames)
    floval = flowstruct.(fnames{i});
    if ~isempty(floval)
      flowstruct.(fnames{i}) = find(linkies == floval);
    end
  end
end
kidsno = setdiff(newlinkies,dadno);
setstruct(dadno).Children = kidsno;
setstruct(dadno).Linked = newlinkies;
setstruct(dadno).Flow = flowstruct;
for i = kidsno
  setstruct(i).Parent = dadno;
  setstruct(i).Linked = newlinkies;
  setstruct(i).Flow = flowstruct;
end

try
  myworkon(DATA.fig);
  save(fullfile(pathname,filename),'preview','info','im','setstruct');
  myworkoff(DATA.fig);
catch %#ok<CTCH>
  myfailed('Could not save data. Write permission? Disk full?',DATA.GUI.Segment);
  myworkoff(DATA.fig);
  return;
end
mymsgbox('Current image stack stored.','',DATA.GUI.Segment);

DATA.LastSaved = now;
flushlog;

%----------------------------------------------------
function savesegmentation_Callback(pathname,filename)
%----------------------------------------------------
%Saves segmentation as a .seg file. This way of saving contours is not
%recommended and may be depreceiated.
global DATA SET NO

segment('stopmovie_Callback');

%---Find what directory
if nargin==0
  pathname = DATA.Preview.PathName;
  %pathname = DATA.Pref.exportpath;
  [filename,pathname] = myuiputfile([pathname filesep '*.seg'],dprintf('Save segmentation as'));
  if isequal(filename,0)||isequal(pathname,0)
    mydisp('Save segmentation aborted by user');
    return;
  end
else
  if nargin<2
    myfailed('Too few input arguments.',DATA.GUI.Segment);
    return;
  end
end

[~,~,ext] = fileparts(fullfile(pathname,filename));
if isempty(ext)
  filename = [filename '.seg'];
end

SEG.ProgramVersion = DATA.ProgramVersion;
SEG.EndoX = SET(NO).EndoX;
SEG.EndoY = SET(NO).EndoY;
SEG.EndoPinX = SET(NO).EndoPinX;
SEG.EndoPinY = SET(NO).EndoPinY;
SEG.EndoInterpX = SET(NO).EndoInterpX;
SEG.EndoInterpY = SET(NO).EndoInterpY;
SEG.EndoDraged = SET(NO).EndoDraged;
SEG.EpiDraged = SET(NO).EpiDraged;
SEG.EpiX = SET(NO).EpiX;
SEG.EpiY = SET(NO).EpiY;
SEG.EpiPinX = SET(NO).EpiPinX;
SEG.EpiPinY = SET(NO).EpiPinY;
SEG.EpiInterpX = SET(NO).EpiInterpX;
SEG.EpiInterpY = SET(NO).EpiInterpY;
SEG.RVEndoX = SET(NO).RVEndoX;
SEG.RVEndoY = SET(NO).RVEndoY;
SEG.RVEpiX = SET(NO).RVEpiX;
SEG.RVEpiY = SET(NO).RVEpiY;
SEG.RVEndoPinX = SET(NO).RVEndoPinX;
SEG.RVEndoPinY = SET(NO).RVEndoPinY;
SEG.RVEndoInterpX = SET(NO).RVEndoInterpX;
SEG.RVEndoInterpY = SET(NO).RVEndoInterpY;
SEG.RVEpiPinX = SET(NO).RVEpiPinX;
SEG.RVEpiPinY = SET(NO).RVEpiPinY;
SEG.RVEpiInterpX = SET(NO).RVEpiInterpX;
SEG.RVEpiInterpY = SET(NO).RVEpiInterpY;
SEG.Name = SET(NO).PatientInfo.Name;
SEG.PatientInfo = SET(NO).PatientInfo;
SEG.HeartRate = SET(NO).HeartRate;
SEG.LVV = SET(NO).LVV;
SEG.PV = SET(NO).PV;
SEG.EPV = SET(NO).EPV;
SEG.EDV = SET(NO).EDV;
SEG.ESV = SET(NO).ESV;
SEG.SV = SET(NO).SV;
SEG.EDT = SET(NO).EDT;
SEG.EST = SET(NO).EST;
SEG.PFR = SET(NO).PFR;
SEG.PFRT = SET(NO).PFRT;
SEG.PER = SET(NO).PER;
SEG.PERT = SET(NO).PERT;
SEG.RVV = SET(NO).RVV;
SEG.RVM = SET(NO).RVM;
SEG.RVEPV = SET(NO).RVEPV;
SEG.RVEDV = SET(NO).RVEDV;
SEG.RVESV = SET(NO).RVESV;
SEG.RVEF = SET(NO).RVEF;
SEG.XMin = SET(NO).XMin;
SEG.YMin = SET(NO).YMin;
SEG.ImagingTechnique = SET(NO).ImagingTechnique;
SEG.SectorRotation = SET(NO).SectorRotation;
SEG.Longaxis = SET(NO).Longaxis;
SEG.AutoLongaxis = SET(NO).AutoLongaxis;
SEG.UseLight = DATA.Pref.UseLight;
SEG.SliceThickness = SET(NO).SliceThickness;
SEG.SliceGap = SET(NO).SliceGap;
SEG.ResolutionX = SET(NO).ResolutionX;
SEG.ResolutionY = SET(NO).ResolutionY;
SEG.TIncr = SET(NO).TIncr;
SEG.TDelay = SET(NO).TDelay;
SEG.VENC = SET(NO).VENC;
SEG.GEVENCSCALE = SET(NO).GEVENCSCALE;
SEG.Rotated = SET(NO).Rotated;
SEG.Roi = SET(NO).Roi;
SEG.RoiCurrent = SET(NO).RoiCurrent;
SEG.RoiN = SET(NO).RoiN;
SEG.OrgXSize = SET(NO).OrgXSize;
SEG.OrgYSize = SET(NO).OrgYSize;
SEG.OrgZSize = SET(NO).OrgZSize;
SEG.OrgTSize = SET(NO).OrgTSize;
if not(isempty(SET(NO).Scar))
  SEG.Scar = SET(NO).Scar;
end
try
  myworkon(DATA.fig);
  save(fullfile(pathname,filename),'SEG', DATA.Pref.SaveVersion);
  myworkoff(DATA.fig);
catch %#ok<CTCH>
  myfailed('Could not save data. Write permission? Disk full?',DATA.GUI.Segment);
  myworkoff(DATA.fig);
  return;
end

DATA.LastSaved = now;

%---------------------------------
function success = savesegdicom_Callback(filename) 
%---------------------------------
%Save image stack as DICOM file
global DATA SET NO

success = false;
if not(DATA.DataLoaded)
  myfailed('Nothing to save.',DATA.GUI.Segment);
  return;
end

%--- Get filename and pathname
if(nargin < 1)
  temp = DATA.SegmentFolder;
  if exist(DATA.Pref.exportpath,'dir')
    cd(DATA.Pref.exportpath);
  else
    mydisp('Warning: Export path does not exist, please check preferences.');
  end
  filename = [SET(NO).PatientInfo.Name '.segdicom'];
  [filename,pathname] = myuiputfile(filename,dprintf('Save all image stacks to one file'));
  cd(temp);
  if isequal(filename,0)
%     myfailed('Operation cancelled.',DATA.GUI.Segment);
    return;
  end
  if length(filename)>=9
    if ~isequal(filename(end-8:end),'.segdicom')
      filename = [filename '.segdicom'];
    end
  else
    filename = [filename '.segdicom'];
  end
  filename = fullfile(pathname, filename);
end

%Store current view
SET(1).View = [];
SET(1).View.ViewPanels = DATA.ViewPanels;
SET(1).View.ViewPanelsType = DATA.ViewPanelsType;
SET(1).View.ViewPanelsMatrix = DATA.ViewPanelsMatrix;
SET(1).View.ViewMatrix = DATA.ViewMatrix;
SET(1).View.ThisFrameOnly = DATA.ThisFrameOnly;
SET(1).View.CurrentPanel = DATA.CurrentPanel;
SET(1).View.CurrentTheme = DATA.CurrentTheme;
SET(1).View.RelevantMode = DATA.ShowRelevantStacksOnly;
SET(1).View.CurrentTool = DATA.CurrentTool;
SET(1).View.LVNO = DATA.LVNO;
SET(1).View.RVNO = DATA.RVNO;
SET(1).View.FlowNO = DATA.FlowNO;
SET(1).View.FlowROI = DATA.FlowROI;

% Save the data
study_uid = segment('getfieldifcommon',SET, 'StudyUID');
%study_id = segment('getfieldifcommon',SET, 'StudyID');
id = SET(1).PatientInfo.ID;
for no = 2:numel(SET)
  if ~strcmp(SET(no).PatientInfo.ID,id)
    id = [];
    break
  end
end

if ~isempty(id)
  PatientInfo = SET(1).PatientInfo;
else
  PatientInfo = [];
end

try
  if isempty(PatientInfo)
    patname = 'Mixed patients';
    patid = '';
    patbd = '';
    patsex = '';
  else
    patname = PatientInfo.Name;
    patid = PatientInfo.ID;
    patbd = PatientInfo.BirthDate;
    patsex = PatientInfo.Sex;
  end
  devicemodelname = helperfunctions('getdevicemodelname');
  seriesdescription = helperfunctions('getseriesdescriptionfordevice');
  segdicomfile.create(filename, SET, study_uid, ...
    patname, patid, ...
    patbd, patsex,DATA.Pref.Pacs.SwitchTags,...
    devicemodelname,seriesdescription);

catch me
  if strcmp(me.identifier, 'SEGMENT:ERROR')
    mydispexception(me)
    myfailed(me.message);
  else
    rethrow(me);
  end
end
success = true;

%--------------------------------------------------------
function closecurrentimagestack_Callback(frompreviewmenu) 
%--------------------------------------------------------
%Close current image stack, i.e the current image stack is deleted. It
%also takes care of eventual cross couplings between image stacks.
global DATA SET NO

flushlog;

if nargin>0 && frompreviewmenu
%don't need to get clicked coordinates since when bringing up previemenu switchimagestack is called
% checking clicked coordinates also results in errors since CurrentPoint
% not alwasy is inside the previewaxes but rather in the previewmenu when
% clicking close 

%code was previously:
%   %Get clicked coordinate
%   p = get(gca,'CurrentPoint');
%   y=p(1,2); % vertical
%   no = floor(y/DATA.GUISettings.ThumbnailSize)+1;
%   DATA.switchtoimagestack(no);
  no = segment_main('thumbnailno');
%   if no ~= NO
%     % found mismatch between stacks 
%     ind = mymenu('Mismatched stacks numbers. Choose stack to delete',['Image stack no ' num2str(no)], ['Image stack no ' num2str(NO)]);
%     if ind == 0
%       %cancel option
%       return
%     elseif ind == 1 
%     elseif ind == 2
%       no = NO;
%     end
%   end
else
  no = NO;
end

if ~DATA.Silent
  if (DATA.NeedToSave)
    qstr = dprintf('Unsaved changes: are you sure you want to delete image stack %d?',no);
    if ~yesno(qstr,[],DATA.GUI.Segment)
      return;
    end
  end
end
logdisp(['Deleted stack ',num2str(no)])
if length(SET) < 2
  DATA.filecloseall_Callback(true);
  return;
end
closeimagestack(no);

%If Segment 3DPrint then take extra care when deleting image stack
if DATA.issegment3dp
  try
    delete(DATA.LevelSet.ViewPort);
    DATA.LevelSet.ViewPort = [];
  catch
    %Do nothing if could not be deleted, likely already deleted or empty
  end
  viewfunctions('setview'); %Same as clicking refresh
end

flushlog;

%----------------------------------------------
function closemultipleimagestacks_Callback(arg)
%----------------------------------------------
%Close multiple user selected image stacks. 
%Useful when a lot of stacks are open.

global DATA SET NO

if nargin < 1
  %Launch GUI
  fig = mygui('closemultiple.fig');
  handles = fig.handles;
  stackc = cell(1,numel(SET));
  for no = 1:numel(SET)
    stackc{no} = sprintf('%d. %s, %s',no,SET(no).ImageType,SET(no).ImageViewPlane);
  end
  set(handles.imagestackslistbox,'String',stackc);
else
  handles = guihandles(gcbf);
  %to enable multiple selection, set max = 2.0
  nos = mygetvalue(handles.imagestackslistbox);
  if strcmp(arg,'update')
    kids = [SET(nos).Children];
    set(handles.imagestackslistbox,'Value',union(nos,kids));
  elseif strcmp(arg,'doclose')
    myworkon(DATA.fig)
    if isequal(length(nos), length(SET))
      %All imagestacks are selected=> Call closeallimagestack function:
      segment('filecloseall_Callback');
    else
      for no = nos(end:-1:1)
        closeimagestack(no)
      end
      viewfunctions('setview',1,1,NO,{'one'})%segment('switchtoimagestack',NO,true); %NO set in subfcn, force
      if DATA.ShowRelevantStacksOnly
        viewfunctions('relevantmode_Callback');
      end
    end
    myworkoff(DATA.fig);
    flushlog;
    close(handles.figure1);
  end
end

%---------------------------
function closeimagestack(no)
%---------------------------
%Close image stack no. Takes care of possible cross couplings 
%between image stacks

global DATA SET NO

%Delete autoloader warning for the stack if there is one
if ~isempty(SET(1).Autoloader)
  if SET(1).Autoloader.Autoloaded
    warnings = SET(1).Autoloader.Warnings;
    if ~isempty(warnings)
      autoloaderstackid = SET(no).StackID;
      ind = find(autoloaderstackid == [warnings.stackID]);
      warnings(ind) = []; %#ok<FNDSB>
      SET(1).Autoloader.Warnings = warnings;
    end
  end
end

%Make sure report is preserved, if removing stack number 1
if no == 1
  SET(2).Report = SET(1).Report;
  SET(2).Autoloader = SET(1).Autoloader;
end


%Close all associated GUIs
DATA.closeallnoguis(no);
resetFlowNO = 0;

if ~isempty(SET(no).Parent) 
  %Remove as overlay from parent stack
  parentno = SET(no).Parent;
  if isequal(SET(parentno).Overlay,no)
    SET(parentno).Overlay = [];
  end
  % remove flow calculations from parent stack  
   SET(parentno).Flow = [];
   % go over ROI in parent stack
   numrois = length(SET(parentno).Roi);
   for r = 1:numrois
     if isfield(SET(parentno).Roi(r),'Flow')
       SET(parentno).Roi(r).Flow = [];
       resetFlowNO = 1;
     end
     if isfield(SET(parentno).Roi(r),'FlowSnake')
       SET(parentno).Roi(r).FlowSnake =[];
     end
   end
end

%If parent, remove all child stacks
ind = true(1,length(SET));
family = [no SET(no).Children];
newlinkies = setdiff(SET(no).Linked,family);
for noloop = newlinkies
  SET(noloop).Linked = newlinkies;
  SET(noloop).Children = setdiff(SET(noloop).Children,family);
end
ind(family) = false;
%correct global LVNO, RVNO and FlowNO 

resetLVNO = 0;
resetRVNO = 0;

if not(isempty(DATA.LVNO))
  if any(family == DATA.LVNO)
    resetLVNO = 1;
  elseif any(family < DATA.LVNO)
    decreseby = length(find(family<DATA.LVNO));
    DATA.LVNO = DATA.LVNO-decreseby;
  end
end
if not(isempty(DATA.RVNO))
  if any(family == DATA.RVNO)
    resetRVNO = 1;
  elseif any(family < DATA.RVNO)
    decreseby = length(find(family < DATA.RVNO));
    DATA.RVNO = DATA.RVNO-decreseby;
  end
end
if not(isempty(DATA.FlowNO))
  if any(family == DATA.FlowNO)
    resetFlowNO = 1;
  elseif any(family < DATA.FlowNO)
    decreseby = length(find(family < DATA.FlowNO));
    DATA.FlowNO = DATA.FlowNO-decreseby;
  end
end

%Here it is actually deleted
SET = SET(ind);

numstacks = length(SET);
if not(isempty(DATA.LVNO)) && (any(DATA.LVNO < 1) || any(DATA.LVNO > numstacks))
  DATA.LVNO = [];
  resetLVNO = 1;
end
if not(isempty(DATA.RVNO)) && (any(DATA.RVNO < 1) || any(DATA.RVNO > numstacks))
  DATA.RVNO = [];
  resetRVNO = 1;
end
if not(isempty(DATA.FlowNO)) && (any(DATA.FlowNO < 1) || any(DATA.FlowNO > numstacks))
  DATA.FlowNO = [];
  resetFlowNO = 1;
end

%fix references for linked images
linkfields = {'Parent' 'Children' 'Linked'};
taggroup = zeros(1,numstacks);
laxgroup = taggroup;

for setloop = 1:numstacks
  for i = 1:length(linkfields)
      linkfield = linkfields{i};
      lks = SET(setloop).(linkfield);
      for stackloop = 1:length(lks)
        %decrease by 1 for each removed stack index < this stack
        lks(stackloop) = lks(stackloop) - sum(family < lks(stackloop));
      end
      SET(setloop).(linkfield) = lks;
  end
  if isstruct(SET(setloop).Flow)
    flowfields = {'MagnitudeNo','PhaseNo','PhaseX','PhaseY','Angio','VelMag'};
    for i = 1:length(flowfields)
      flowfield = flowfields{i};
      ffs = SET(setloop).Flow.(flowfield);
      for stackloop = 1:length(ffs)
        %decrease by 1 for each removed stack index < this stack
        ffs(stackloop) = ffs(stackloop) - sum(family < ffs(stackloop));
      end
      SET(setloop).Flow.(flowfield) = ffs;
    end
  end
  if isfield(SET(setloop),'StrainTagging') && ~isempty(SET(setloop).StrainTagging)
    if isfield(SET(setloop).StrainTagging,'cineno')
      if ismember(SET(setloop).StrainTagging.cineno,family)
        SET(setloop).StrainTagging.cineno = setloop;
        SET(setloop).StrainTagging.cineupdated = true;        
      else
        increment = sum(family <= setloop);
        if length(family) > 1 && increment == 1
          % increase increment to 2, because we are deleting a linked stack
          increment = increment + 1;
        end
        SET(setloop).StrainTagging.cineno = SET(setloop).StrainTagging.cineno - increment;
      end
      % write value of cineno into the corresponding taggroup array
      taggroup(1,setloop) = SET(setloop).StrainTagging.cineno;
    end
  end
  if isfield(SET(setloop),'StrainMitt') && ~isempty(SET(setloop).StrainMitt)
    if strcmpi(SET(setloop).StrainMitt.ImageViewPlane,'lax')
      if ismember(SET(setloop).StrainMitt.OriginalNO,family)
        SET(setloop).StrainMitt.OriginalNO = setloop;
      else
        increment = sum(family <= SET(setloop).StrainMitt.OriginalNO);
          if length(family) > 1 && increment == 1
            % increase increment to 2, because we are deleting a linked stack
            increment = increment + 1;
          end
          SET(setloop).StrainMitt.OriginalNO = SET(setloop).StrainMitt.OriginalNO - increment;
          
          if ~isempty(SET(setloop).StrainMitt.LAXGroup)
            laxgroup(1,setloop) = SET(setloop).StrainMitt.OriginalNO;
          end
      end
    else
      if ismember(SET(setloop).StrainMitt.OriginalNO,family)
        SET(setloop).StrainMitt.OriginalNO = setloop;
      else
        increment = sum(family <= SET(setloop).StrainMitt.OriginalNO);
          if length(family) > 1 && increment == 1
            % increase increment to 2, because we are deleting a linked stack
            increment = increment + 1;
          end
          SET(setloop).StrainMitt.OriginalNO = SET(setloop).StrainMitt.OriginalNO - increment;
      end
    end
  end
end

% cleanup taggroup 
taggroup = (nonzeros(taggroup))';
if ~isempty(taggroup)
  for tagloop = taggroup
    % write new taggroup into SET struct
    SET(tagloop).StrainTagging.taggroup = taggroup;
  end
end

% cleanup laxgroup 
laxgroup = (nonzeros(laxgroup))';
if ~isempty(laxgroup)
  for tagloop = laxgroup
    % write new taggroup into SET struct
    SET(tagloop).StrainMitt.LAXGroup = laxgroup;
  end
end

%Things for Segment 3DPrint
if DATA.issegment3dp

  %Loop over sets and ensure that O.no is correcct
  for loop = 1:length(SET)
    O = SET(loop).LevelSet.Object;
    O.NO = loop;
  end

  %Update 
  NO = 1;
  O = SET(NO).LevelSet.Object;
  if ~isempty(O)
    O.updateobjectlist();
    if isfield(DATA.LevelSet,'ViewPort')
      v = DATA.LevelSet.ViewPort;
      if (~isempty(v)) && (~isdeleted(v))
        delete(v); %Close  the viewport
      end
    end
  end

end

if isempty(SET)
  DATA.filecloseall_Callback(true);
  return;
end

NO = 1;

%Reset NOS if deleted 
cinestacks = findfunctions('findcineshortaxisno','true');

if resetLVNO 
  if ~isempty(cinestacks)
    DATA.LVNO = cinestacks(1);
  else
    DATA.LVNO = [];
  end
end
if resetRVNO
  if ~isempty(cinestacks)
    DATA.RVNO = cinestacks(2);
  else
    DATA.RVNO = [];
  end
end
if resetFlowNO
  [DATA.FlowNO, DATA.FlowROI] = findfunctions('findflowaxisno');
end

if DATA.issegment3dp
  segment3dp.callbackfunctions('view4panel_Callback');
  newno = NO;
else  
  % check if it was deleted in "ShowRelevantMode"
  if DATA.ShowRelevantStacksOnly
    DATA.RelevantStacks = setdiff(DATA.RelevantStacks, family, 'stable');
    if ~isempty(DATA.RelevantStacks)
      % Calculate decrements for all elements
      decrements = arrayfun(@(x) sum(family < x), DATA.RelevantStacks);
      DATA.RelevantStacks = DATA.RelevantStacks - decrements;
      NO = DATA.RelevantStacks(1);
    else
      viewfunctions('relevantmode_Callback')
      return
    end
  end
  newno = NO;
  DATA.ViewMatrix = [1 1];
  viewfunctions('setview',1,1,newno,{'one'}); %updates result tables
end
drawfunctions('drawthumbnails',newno);

%----------------------------------------------------
function loadsegmentation_Callback(pathname,filename) 
%----------------------------------------------------
%Loads segmentation to current image stack from a .seg file.
global DATA SET NO

tools('enableundo');

if nargin==0
  %Find what directory
  pathname = DATA.Preview.PathName;
  
  [filename,pathname] = myuigetfile([pathname filesep '*.seg'],'Load segmentation from');
  if isequal(filename,0)
    mydisp('Load segmentation aborted by user');
    return;
  end
else
  if nargin<2
    myfailed('Expected two input arguments.',DATA.GUI.Segment);
    return;
  end
end

SEG = []; %#ok<NASGU> 
try
  myworkon(DATA.fig);
  load(fullfile(pathname,filename),'SEG','-mat');
  myworkoff(DATA.fig);
catch %#ok<CTCH>
  myfailed('Could not load data. Data corrupted?',DATA.GUI.Segment);
  myworkoff(DATA.fig);
  return;
end

if isempty(SEG)
  myfailed('The file did not contain any valid segmentation.',DATA.GUI.Segment);
  return;
end

%Change this check later
if ~isempty(SEG.EndoX)
  if size(SEG.EndoX,3)~=SET(NO).ZSize
    myfailed('Number of slices must be the same as current number of slices.',DATA.GUI.Segment);
    return;
  end

  if size(SEG.EndoX,2)~=SET(NO).TSize
    myfailed('Number of timeframes must be the same as current number timeframes.',DATA.GUI.Segment);
    return;
  end
end

if str2num(removechars(SEG.ProgramVersion)) > str2num(removechars(DATA.ProgramVersion)) %#ok<ST2NM>
  mydisp('Warning, data has been generated with a later version that this program.');
end

if str2num(removechars(SEG.ProgramVersion))<0.98 %#ok<ST2NM>
  mydisp('File is older than v0.98');
  is098 = true;
else
  is098 = false;
end

try
  %--- Ok should be trustworthy now
  [~,~,ext] = fileparts(SET(NO).FileName);
  if not(DATA.Silent)
    if not(isequal(ext,'.mat'))
      if yesno('Loading segmentation to cropped image data?',[],DATA.GUI.Segment)
        xofs = SEG.XMin-SET(NO).XMin;
        yofs = SEG.YMin-SET(NO).YMin;
      else
        xofs = 0;
        yofs = 0;
      end
    else
      xofs = 0;
      yofs = 0;
    end
  else
    %Silent
    xofs = 0;
    yofs = 0;
  end
  
  SET(NO).EndoX = SEG.EndoX+xofs;
  SET(NO).EndoY = SEG.EndoY+yofs;
  SET(NO).EpiX = SEG.EpiX+xofs;
  SET(NO).EpiY = SEG.EpiY+yofs;
  if is098
    datanumpoints = tools('getnumpointsforno',NO);
    [SET(NO).EpiX,SET(NO).EpiY] = calcfunctions('resamplemodel',SET(NO).EpiX,SET(NO).EpiY,datanumpoints);
  end

  SET(NO).EndoPinX = SEG.EndoPinX; %They are translated below.
  SET(NO).EndoPinY = SEG.EndoPinY;
  SET(NO).EpiPinX = SEG.EpiPinX;
  SET(NO).EpiPinY = SEG.EpiPinY;
  
  SET(NO).RVEndoPinX = SEG.RVEndoPinX;
  SET(NO).RVEndoPinY = SEG.RVEndoPinY;  
  SET(NO).RVEpiPinX = SEG.RVEpiPinX;
  SET(NO).RVEpiPinY = SEG.RVEpiPinY;
  
  if isfield(SEG,'EndoInterpX')
    SET(NO).EndoInterpX = SEG.EndoInterpX;
    SET(NO).EndoInterpY = SEG.EndoInterpY;
    SET(NO).EpiInterpX = SEG.EpiInterpX;
    SET(NO).EpiInterpY = SEG.EpiInterpY;
    if isfield(SEG,'RVEndoInterpX')
      SET(NO).RVEndoInterpX = SEG.RVEndoInterpX;
      SET(NO).RVEndoInterpY = SEG.RVEndoInterpY;
    elseif isfield(SEG,'EndoRVInterpX')%backwardscompability
      SET(NO).RVEndoInterpX = SEG.EndoRVInterpX;
      SET(NO).RVEndoInterpY = SEG.EndoRVInterpY;   
    end
    SET(NO).RVEpiInterpX = SEG.RVEpiInterpX;
    SET(NO).RVEpiInterpY = SEG.RVEpiInterpY;
  else
    SET(NO).EndoInterpX = [];
    SET(NO).EndoInterpY = [];
    SET(NO).EpiInterpX = [];
    SET(NO).EpiInterpY = [];
    SET(NO).RVEndoInterpX = [];
    SET(NO).RVEndoInterpY = [];  
    SET(NO).RVEpiInterpX = [];
    SET(NO).RVEpiInterpY = [];
  end
  
  SET(NO).EndoDraged = SEG.EndoDraged;
  SET(NO).EpiDraged = SEG.EpiDraged;

  SET(NO).SectorRotation = SEG.SectorRotation;
  SET(NO).EDT = SEG.EDT;
  SET(NO).EST = SEG.EST;
  
  if isfield(SEG,'Colormap')
    SET(NO).Colormap=SEG.Colormap;
  end
  
  if not(isfield(SEG,'ImagingTechnique'))
    %Old file format
    SET(NO).ImagingTechnique = SEG.ImageType;
  else
    SET(NO).ImagingTechnique = SEG.ImagingTechnique;
    if not(isfield(SEG,'ImageType'))
      SET(NO).ImageType = 'General';
    else
      SET(NO).ImageType = SEG.ImageType;
    end
  end

  if not(isfield(SEG,'OrgXSize'))
    if not(DATA.Silent)
      mywarning('Old file format did not contain original DICOM image size, guess on 256.',DATA.GUI.Segment);
    end
    SET(NO).OrgXSize = 256;
    SET(NO).OrgYSize = 256;
    SET(NO).OrgTSize = SET(NO).TSize;
    SET(NO).OrgZSize = SET(NO).ZSize;
  end
  
  if isfield(SEG,'ExcludePapilars')
      mydisp('Ignoring ExcludePapilars on imported segmentation');
   % DATA.ExcludePapilars = SEG.ExcludePapilars;
    %set(DATA.Handles.excludepapilarscheckbox,'value',DATA.ExcludePapilars);
  end
  
  if isfield(SEG,'AutoLongAxis')
    SET(NO).AutoLongaxis = SEG.AutoLongaxis;
  end
  
  if isfield(SEG,'UseLight')
      mydisp('Ignoring UseLight on imported segmentation');
    %DATA.UseLight = SEG.UseLight;
    %set(DATA.Handles.uselightcheckbox,'value',SEG.UseLight);
  end
  
  if isfield(SEG,'Scar')
    if not(DATA.Silent)
      mywarning('Contains scar data, loading not 100% supported yet.',DATA.GUI.Segment);
    end
    SET(NO).Scar = SEG.Scar;
    SET(NO).Scar.Mode = 'manualthreshold';
    if not(isfield(SET(NO).Scar,'NoReflow'))
      SET(NO).Scar.NoReflow = repmat(uint8(0),size(SET(NO).Scar.Manual));
    end
  end
  
  if isfield(SEG,'RVEndoX')
    SET(NO).RVEndoX = SEG.RVEndoX;
    SET(NO).RVEndoY = SEG.RVEndoY;
    SET(NO).RVEpiX = SEG.RVEpiX;
    SET(NO).RVEpiY = SEG.RVEpiY;
  end
  
  SET(NO).Longaxis = SEG.Longaxis;
    
  if ~isempty(SET(NO).EndoPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).EndoPinX{tloop,zloop} = SET(NO).EndoPinX{tloop,zloop}+xofs;
        SET(NO).EndoPinY{tloop,zloop} = SET(NO).EndoPinY{tloop,zloop}+yofs;
      end
    end
  end
  
  if ~isempty(SET(NO).EpiPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).EpiPinX{tloop,zloop} = SET(NO).EpiPinX{tloop,zloop}+xofs;
        SET(NO).EpiPinY{tloop,zloop} = SET(NO).EpiPinY{tloop,zloop}+yofs;
      end
    end
  end

  if ~isempty(SET(NO).RVEndoPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).RVEndoPinX{tloop,zloop} = SET(NO).RVEndoPinX{tloop,zloop}+xofs;
        SET(NO).RVEndoPinY{tloop,zloop} = SET(NO).RVEndoPinY{tloop,zloop}+yofs;
      end
    end
  end

  if ~isempty(SET(NO).RVEpiPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).RVEpiPinX{tloop,zloop} = SET(NO).RVEpiPinX{tloop,zloop}+xofs;
        SET(NO).RVEpiPinY{tloop,zloop} = SET(NO).RVEpiPinY{tloop,zloop}+yofs;
      end
    end
  end
  
  if isfield(SEG,'Rotated')
    SET(NO).Rotated = SEG.Rotated;
  end
  
  if isfield(SEG,'RoiN')
    SET(NO).RoiCurrent = SEG.RoiCurrent;
    SET(NO).RoiN = SEG.RoiN;
    if isfield(SEG,'Roi')
      SET(NO).Roi=SEG.Roi;
      for loop=1:SET(NO).RoiN
        SET(NO).Roi(loop).X=SET(NO).Roi(loop).X+xofs;
        SET(NO).Roi(loop).Y=SET(NO).Roi(loop).Y+yofs;
      end
    else
      for loop=1:SET(NO).RoiN
        SET(NO).Roi(loop).X = SEG.RoiX(:,:,loop)+xofs;
        SET(NO).Roi(loop).Y = SEG.RoiY(:,:,loop)+yofs;
        SET(NO).Roi(loop).T = 1:SET(NO).TSize;
        SET(NO).Roi(loop).Z = SEG.RoiZ(loop);
        SET(NO).Roi(loop).Name = SEG.RoiName{loop};
        SET(NO).Roi(loop).Sign = SEG.RoiSign(loop);
        SET(NO).Roi(loop).LineSpec = SEG.RoiLineSpec{loop};
      end
    end
  end
  
  if isfield(SEG,'Measure')
    SET(NO).Measure = SEG.Measure;
    for loop=1:length(SET(NO).Measure)
      SET(NO).Measure(loop).X = SET(NO).Measure(loop).X+xofs;
      SET(NO).Measure(loop).Y = SET(NO).Measure(loop).Y+xofs;
    end
  end
  
  if isfield(SEG,'Point')
    SET(NO).Point = SEG.Point;
    SET(NO).Point.X = SET(NO).Point.X+xofs;
    SET(NO).Point.Y = SET(NO).Point.Y+yofs;    
  end
  
  if isfield(SEG,'GEVENCSCALE')
    SET(NO).GEVENCSCALE = SEG.GEVENCSCALE;
  end
  
catch me
  myfailed('Something went wrong during loading. Undoing...',DATA.GUI.Segment);
  mydispexception(me);
  tools('undosegmentation_Callback');
end

if size(SET(NO).EndoDraged,1)~=SET(NO).TSize
  SET(NO).EndoDraged = repmat(SET(NO).EndoDraged(:),1,SET(NO).TSize);
  SET(NO).EpiDraged = repmat(SET(NO).EpiDraged(:),1,SET(NO).TSize);
end

%Prevent problems from old .seg files
segment('checkconsistency',1:SET(NO).TSize,1:SET(NO).ZSize);

segment('updatevolume');
%segment('viewrefresh_Callback');
viewfunctions('setview')

%---------------------
function quit_Callback(varargin)
%---------------------
%Quit Segment
global DATA SET NO

if ~isempty(varargin)
  str = varargin{1};
else
  str = 'Closing software.';
end

if ~isempty(DATA)
  if DATA.isSiemensVersion
    % save results first
    success = applysiemensoptions('savetosiemens_Callback',true);
    if success
      mydisp('Data was successfully saved');
    else
      mydisp('Could not save data');
    end
    applysiemensoptions('copyfinalresultstosyngovia');
    % move final results from temporal folder to the result folder
    % delete possible report pdf 
    folderpath = getpreferencespath;
    delete([folderpath filesep '*.pdf'])
    fig = DATA.fig;
    DATA = [];
    SET = [];
    NO = [];
    clear DATA
    clear SET
    clear NO
    delete(fig); %Call destroyer to take it down nice, see code =>
    %close('all'); %Added do not know if it works.
    close all hidden % closes all figures, even the hidden ones
    commandlinehelper('reset');
  else
    if DATA.quit
      DATA.adduserevent(sprintf('%s\t%s\t%s\t%s', datestr(now,'yyyy-mm-dd HH:MM'),DATA.LoggedUser,str,'-'));
      %DATA.adduserevent(['Time:' datestr(now,'yyyymmddHHMMSS')])
      saveguiposition(DATA.GUI.Segment);
      segment('saveguipositiontodisk');
      if DATA.Pref.UseCache && DATA.issegmentcmr
        cachefunctions('deletecache');
      end
      fig = DATA.fig;
      DATA = [];
      SET = [];
      NO = [];
      clear DATA
      clear SET
      clear NO
      delete(fig); %Call destroyer to take it down nice, see code =>
      %close('all'); %Added do not know if it works.
      %close all hidden % closes all figures, even the hidden ones
      close all force % closes all figures
      commandlinehelper('reset');
    end
  end
else
  delete(gcbo);
end