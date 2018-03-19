function varargout = filemenu(varargin)
% FILEMENU
% File menu callbacks

% Nisse Lundahl

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
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
end;

if length(SET)<1
  loadednext = 0;
  return;
end;

%Get filename from first image stack.
filename = SET(1).FileName;

[pathname,filename,ext] = fileparts(filename);

f = dir([pathname filesep '*.mat']);

pos = [];
for loop=1:length(f)
  if isequal(f(loop).name,[filename '.mat']);
    pos = loop;
  end;
end;

if isempty(pos)
  myfailed('Could not find the current file, removed?');
  loadednext = 0;
  return;
end;

if (pos+1)<=length(f)
  
  %Ask to quite nicely
  DATA.filecloseall_Callback;
  
  if ~isempty(SET)
    myfailed('Loading next aborted.');
    loadednext = 0;
    return;
  end;
  
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
    
    DATA.NeedToSave = false;
    
  catch me
    mydispexception(me);
  end;
  loadednext = 1;
else
  myfailed('Last file');
  loadednext = 0;
end;

%--------------------------------------
function savetopatientdatabase_Callback
%--------------------------------------
%Callback to save image stacks to patientdatabase. Uses functions in
%patientdatabase.
global DATA SET NO

set(DATA.Handles.databaseaddicon,'state','off');

[pathname,filename,ext] = fileparts(SET(NO).FileName);

try
  patientdatabase('savetodatabase',pathname,filename,ext);
catch me
  mydispexception(me);
end;

%------------------------
function saveall_Callback
%------------------------
%Saves all image stacks to one .mat file. Calls GUI method 
%filesaveallas_Callback which is the workhorse when saving image stacks.
global DATA SET

filenames = unique({SET.FileName});
if numel(filenames) > 1
  i = mymenu('Multiple file names were found in data. Select where to save', ...
    filenames);
  if i > 0
    fname = filenames{i};
  else
    myfailed('Saving aborted');
    return
  end
else
  fname = filenames{1};
end
[pathname,filename,ext] = fileparts(fname);

%--- Normal research mode
if isequal(ext,'.mat')
  DATA.filesaveallas_Callback(pathname,[filename ext]);
else
  %Call Save as instead
  DATA.filesaveallas_Callback;
end;

DATA.LastSaved = now;
flushlog;

%--------------------------------------------------------------------
function fail = saveallas_helper(pathname,filename,topatientdatabase)
%--------------------------------------------------------------------
%Save all image stacks to the file specified. It also stores current view
%and modes etc.

global DATA SET NO
    
if nargin<3
  topatientdatabase = false; %Default save to disk not patientdatabase
  %Setting this to true only changes the final saved display message.
end;

fail=0;

set(DATA.Handles.filesaveicon,'state','off');

% if not(DATA.DataLoaded)
%   myfailed('Nothing to save.',DATA.GUI.Segment);
%   return;
% end;

%--- Get filename and pathname
if nargin==0
  temp = pwd;
  if exist(DATA.Pref.exportpath,'dir')
    cd(DATA.Pref.exportpath);
  else
    mydisp('Warning: Export path does not exist, please check preferences.');
  end;
  
  filename = DATA.generatesavefilename;
  
  [filename,pathname] = myuiputfile(filename,'Save all image stacks to one file');
  cd(temp);
else
  %Take from input two argument
  if nargin<2
    myfailed('Expected two input arguments.',DATA.GUI.Segment);
    return;
  end;
end;

if isequal(filename,0)
  myfailed('Operation cancelled.',DATA.GUI.Segment);
  fail=1;
  return;
end;

if length(filename)>4
  if ~isequal(filename(end-3:end),'.mat')
    filename = [filename '.mat'];
  end;
else
  filename = [filename '.mat'];
end;

myworkon;

%Store filename
for loop=1:length(SET)
  SET(loop).FileName = [pathname filesep filename];
end;

%Store current view
SET(1).View = [];
SET(1).View.ViewPanels = DATA.ViewPanels;
SET(1).View.ViewPanelsType = DATA.ViewPanelsType;
SET(1).View.ViewPanelsMatrix = DATA.ViewPanelsMatrix;
SET(1).View.ViewMatrix = DATA.ViewMatrix;
SET(1).View.ThisFrameOnly = DATA.ThisFrameOnly;
SET(1).View.CurrentPanel = DATA.CurrentPanel;
SET(1).View.CurrentTheme = DATA.CurrentTheme;
SET(1).View.CurrentTool = DATA.CurrentTool;
preview = DATA.DATASETPREVIEW; %#ok<NASGU> %saved to file
info = SET(NO).PatientInfo;

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
im = []; %#ok<NASGU> %saved to file
setstruct = SET; %#ok<NASGU> %saved to file

lastwarn('');
warnchk = 'cannot be saved';
try
  save(fullfile(pathname,filename),'preview','info','im','setstruct', DATA.Pref.SaveVersion);
  warnmsg = lastwarn;
  lastwarn('');
  if ~isempty(strfind(warnmsg,warnchk))
    save(fullfile(pathname,filename),'preview','info','im','setstruct','-v7.3')
    warnmsg = lastwarn;
    lastwarn('');
    if ~isempty(strfind(warnmsg,warnchk))
      error('Dataset cannot be saved');
    end
  end
  myworkoff;
catch %#ok<CTCH>
  if ~DATA.Silent
    myfailed('Could not save data. Write permission? Disk full?',DATA.GUI.Segment);
  end
  fail = 1;
  myworkoff;
  return;
end;

if not(DATA.Silent)
  if topatientdatabase
    mymsgbox('Image Stacks Saved to Patient Database.');
  else
    if length(SET)>1
      stri = dprintf('All image stacks and segmentation stored to %s.',...
        SET(NO).FileName);
    else
      stri = dprintf('Image stack and segmentation stored to %s.',...
        SET(NO).FileName);
    end;
    mymsgbox(stri,'Save successful.',DATA.GUI.Segment);
  end;
end;

DATA.LastSaved = now;
DATA.NeedToSave = false;
set(DATA.Handles.filesaveicon,'enable','off');
DATA.updatetitle;
flushlog;

%---------------------------
function savetopacs_Callback
%---------------------------
%Send image stacks to PACS. This function should display a list of
%available PACS (.con files) and when user has selected store files on disk
%temporarily and then send the files to the PACS.
global DATA

set(DATA.Handles.pacsaddicon,'state','off');

try
  pacs('savetopacs');
  mymsgbox('Image stacks saved to PACS.')
catch e
  myfailed(e.message);
end

%----------------------------
function savecurrent_Callback 
%----------------------------
%Save current image set to file. Note that this is the old Segment file
%format and this fcn may soon be depreciated.
global DATA SET NO

if not(DATA.DataLoaded)
  myfailed('No data to save.',DATA.GUI.Segment);
  return;
end;

temp = pwd;
if exist(DATA.Pref.exportpath,'dir')
  cd(DATA.Pref.exportpath);
else
  mydisp('Warning: Export path does not exist, please check preferences.');
end;

[filename,pathname] = myuiputfile('*.mat','Save current image stack');
cd(temp);

if isequal(filename,0)
  myfailed('Operation cancelled.',DATA.GUI.Segment);
  return;
end;

%IM and preview are saved to disk
preview = SET(NO).IM(:,:,round(SET(NO).TSize/2),round(SET(NO).ZSize/2)); %#ok<NASGU>
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
im = calcfunctions('calctruedata',SET(NO).IM,NO); %#ok<NASGU>
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
  myworkon;
  save(fullfile(pathname,filename),'preview','info','im','setstruct');
  myworkoff;
catch %#ok<CTCH>
  myfailed('Could not save data. Write permission? Disk full?',DATA.GUI.Segment);
  myworkoff;
  return;
end;
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
  [filename,pathname] = myuiputfile([pathname filesep '*.seg'],'Save segmentation as');
  if isequal(filename,0)||isequal(pathname,0)
    mydisp('Save segmentation aborted by user');
    return;
  end;
else
  if nargin<2
    myfailed('Too few input arguments.',DATA.GUI.Segment);
    return;
  end;
end;

[temppath,tempname,ext] = fileparts(fullfile(pathname,filename));
if isempty(ext)
  filename = [filename '.seg'];
end;

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
end;
try
  myworkon;
  save(fullfile(pathname,filename),'SEG', DATA.Pref.SaveVersion);
  myworkoff;
catch %#ok<CTCH>
  myfailed('Could not save data. Write permission? Disk full?',DATA.GUI.Segment);
  myworkoff;
  return;
end;

DATA.LastSaved = now;

%---------------------------------
function success = savesegdicom_Callback(filename) %#ok<DEFNU>
%---------------------------------
%Save image stack as DICOM file
global DATA SET NO

success = false;
if not(DATA.DataLoaded)
  myfailed('Nothing to save.',DATA.GUI.Segment);
  return;
end;

%--- Get filename and pathname
if(nargin < 1)
  temp = pwd;
  if exist(DATA.Pref.exportpath,'dir')
      cd(DATA.Pref.exportpath);
  else
      mydisp('Warning: Export path does not exist, please check preferences.');
  end;
  filename = [SET(NO).PatientInfo.Name '.segdicom'];
  [filename,pathname] = myuiputfile(filename,'Save all image stacks to one file');
  cd(temp);
  if isequal(filename,0)
    myfailed('Operation cancelled.',DATA.GUI.Segment);
    return;
  end;
  if length(filename)>=9
    if ~isequal(filename(end-8:end),'.segdicom')
      filename = [filename '.segdicom'];
    end;
  else
    filename = [filename '.segdicom'];
  end;
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
SET(1).View.CurrentTool = DATA.CurrentTool;

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
    segdicomfile.create(filename, SET, study_uid, ...
      'Mixed patients', '', '', '',DATA.Pref.Pacs.SwitchTags);
  else
    segdicomfile.create(filename, SET, study_uid, ...
      PatientInfo.Name, PatientInfo.ID, ...
      PatientInfo.BirthDate, PatientInfo.Sex,DATA.Pref.Pacs.SwitchTags);
  end
catch e
  if strcmp(e.identifier, 'SEGMENT:ERROR')
    myfailed(e.message);
  else
    rethrow(e);
  end
end
success = true;

%--------------------------------------------------------
function closecurrentimagestack_Callback(frompreviewmenu) %#ok<INUSD>
%--------------------------------------------------------
%Close current image stack, i.e the current image stack is deleted. It
%also takes care of eventual cross couplings between image stacks.
global DATA SET NO

flushlog;

if nargin>0
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
end;

if ~DATA.Silent
  if (DATA.NeedToSave)
    if ~yesno('Unsaved changes: are you sure you want to delete the current image stack?',[],DATA.GUI.Segment);
      return;
    end;
  end
end;

if length(SET) < 2
  DATA.filecloseall_Callback(true);
  return;
end
closeimagestack(NO);
DATA.switchtoimagestack(NO,true); %force
flushlog;

%----------------------------------------------
function closemultipleimagestacks_Callback(arg)
%----------------------------------------------
%Close multiple user selected image stacks. 
%Useful when a lot of stacks are open.
global SET NO

if nargin < 1
  %Launch GUI
  fig = openfig('closemultiple.fig');
  handles = guihandles(fig);
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
    if isequal(length(nos), length(SET))
        %All imagestacks are selected=> Call closeallimagestack function: 
        segment('filecloseall_Callback');
    else
        for no = nos(end:-1:1)
            closeimagestack(no)
        end
        segment('switchtoimagestack',NO,true); %NO set in subfcn, force
    end
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
resetLVNO=0;
resetRVNO=0;
resetFlowNO=0;

if no<=DATA.LVNO
resetLVNO=1;
end

if no<=DATA.RVNO
resetRVNO=1;
end

if no<=DATA.FlowNO
resetFlowNO=1;
end

%Make sure report is preserved, if removing stack number 1
if no == 1
  SET(2).Report = SET(1).Report;
end

if ~isempty(SET(no).Flow)
  SET(no).Flow=[];
end

%Close all associated GUIs
DATA.closeallnoguis(no);

%Remove as overlay from parent stack
if ~isempty(SET(no).Parent) && isequal(SET(SET(no).Parent).Overlay,no)
  SET(SET(no).Parent).Overlay = [];
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
if ~isempty(SET(no).Flow)
  DATA.updateaxestables('flowremovestack',family);
  DATA.updateaxestables('areaclearall');
end
DATA.updateaxestables('t2star');
SET = SET(ind);

%fix references for linked images
linkfields = {'Parent' 'Children' 'Linked'};

for setloop = 1:length(SET)
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
        SET(setloop).StrainTagging.cineno = SET(setloop).StrainTagging.cineno - sum(family <= setloop);
      end
    end
  end
end

if isempty(SET)
  DATA.filecloseall_Callback(true);
  return;
end;

%Take away from viewpanels
%zer = zeros(1,length(DATA.ViewPanels));
tempind = not(DATA.ViewPanels==no);
DATA.ViewPanels = DATA.ViewPanels(tempind);
DATA.ViewPanelsType = DATA.ViewPanelsType(tempind);
DATA.ViewPanelsMatrix = DATA.ViewPanelsMatrix(tempind);
DATA.ViewIM = DATA.ViewIM(tempind);
DATA.Overlay = DATA.Overlay(tempind);
DATA.ViewMatrix = [];

% %DATA.ViewPanels have internal references to old image stack ordering.
% for loop=1:length(DATA.ViewPanels)
%   temp = DATA.ViewPanels(loop);
%   if temp>0
%     newpos = zer; %zer frome above;
%     newpos(temp) = 1;
%     newpos = newpos(tempind); %tempind from above
%     newpos = find(newpos);
%     if ~isempty(newpos)
%       DATA.ViewPanels(loop) = newpos;
%     else
%       DATA.ViewPanels(loop) = 0;
%     end;
%   end;
% end;

%update NO
if isempty(DATA.ViewPanels) || all(DATA.ViewPanels==0)
  NO=1;
else
  tmpno = DATA.ViewPanels(find(DATA.ViewPanels,1));
  if tmpno<no
    NO=tmpno;
  else
    NO=tmpno-1;
  end
end
% NO = NO-length(find(ind==0));
% if NO < 1
%   NO = 1;
% end

%Update datasetpreview
ind = repmat(ind,DATA.GUISettings.ThumbnailSize,1);
% vertical
DATA.DATASETPREVIEW = DATA.DATASETPREVIEW(ind(:)',:);

DATA.CurrentPanel = DATA.CurrentPanel-1;
if DATA.CurrentPanel<1
  DATA.CurrentPanel = 1;
end;

% if ~isempty(DATA.FlowNO) && no<DATA.FlowNO
%   DATA.FlowNO=[];
% end
% 
% if ~isempty(DATA.LVNO) && no<DATA.LVNO
%   DATA.LVNO=[];
% end
% 
% if ~isempty(DATA.RVNO) && no<DATA.RVNO
%   DATA.LVNO=[];
% end

% %Switch data does some backup before...
% SET(NO).StartSlice = SET(NO).StartSlice;
% SET(NO).EndSlice = SET(NO).EndSlice;
% SET(NO).CurrentTimeFrame = SET(NO).CurrentTimeFrame;

%if isempty(SET(NO).Scar) && strcmp(DATA.CurrentTheme,'scar')
%  DATA.updateicons('lv');
%end
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

if resetLVNO
  DATA.updateaxestables('volume')
end


if resetRVNO
  DATA.updateaxestables('volume')
end


DATA.ViewMatrix=[1 1];
drawfunctions('drawall',1);
drawfunctions('drawthumbnails');



%----------------------------------------------------
function loadsegmentation_Callback(pathname,filename) %#ok<DEFNU>
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
  end;
else
  if nargin<2
    myfailed('Expected two input arguments.',DATA.GUI.Segment);
    return;
  end;
end;

SEG = [];
try
  myworkon;
  load(fullfile(pathname,filename),'SEG','-mat');
  myworkoff;
catch %#ok<CTCH>
  myfailed('Could not load data. Data corrupted?',DATA.GUI.Segment);
  myworkoff;
  return;
end;

if isempty(SEG)
  myfailed('The file did not contain any valid segmentation.',DATA.GUI.Segment);
  return;
end;

%Change this check later
if ~isempty(SEG.EndoX)
  if size(SEG.EndoX,3)~=SET(NO).ZSize
    myfailed('Number of slices must be the same as current number of slices.',DATA.GUI.Segment);
    return;
  end;

  if size(SEG.EndoX,2)~=SET(NO).TSize
    myfailed('Number of timeframes must be the same as current number timeframes.',DATA.GUI.Segment);
    return;
  end;
end;

if str2num(removechars(SEG.ProgramVersion))>str2num(removechars(DATA.ProgramVersion)) %#ok<ST2NM>
  mydisp('Warning, data has been generated with a later version that this program.');
end;

if str2num(removechars(SEG.ProgramVersion))<0.98 %#ok<ST2NM>
  mydisp('File is older than v0.98');
  is098 = true;
else
  is098 = false;
end;

try
  %--- Ok should be trustworthy now
  [~,~,ext] = fileparts(SET(NO).FileName);
  if not(DATA.Silent)
    if not(isequal(ext,'.mat'))
      if yesno('Loading segmentation to cropped image data?',[],DATA.GUI.Segment);
        xofs = SEG.XMin-SET(NO).XMin;
        yofs = SEG.YMin-SET(NO).YMin;
      else
        xofs = 0;
        yofs = 0;
      end;
    else
      xofs = 0;
      yofs = 0;
    end;
  else
    %Silent
    xofs = 0;
    yofs = 0;
  end;
  
  SET(NO).EndoX = SEG.EndoX+xofs;
  SET(NO).EndoY = SEG.EndoY+yofs;
  SET(NO).EpiX = SEG.EpiX+xofs;
  SET(NO).EpiY = SEG.EpiY+yofs;
  if is098
    [SET(NO).EpiX,SET(NO).EpiY] = calcfunctions('resamplemodel',SET(NO).EpiX,SET(NO).EpiY,DATA.NumPoints);
  end;

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
    end;
  end;

  if not(isfield(SEG,'OrgXSize'))
    if not(DATA.Silent)
      mywarning('Old file format did not contain original DICOM image size, guess on 256.',DATA.GUI.Segment);
    end;
    SET(NO).OrgXSize = 256;
    SET(NO).OrgYSize = 256;
    SET(NO).OrgTSize = SET(NO).TSize;
    SET(NO).OrgZSize = SET(NO).ZSize;
  end;
  
  if isfield(SEG,'ExcludePapilars')
      mydisp('Ignoring ExcludePapilars on imported segmentation');
   % DATA.ExcludePapilars = SEG.ExcludePapilars;
    %set(DATA.Handles.excludepapilarscheckbox,'value',DATA.ExcludePapilars);
  end;
  
  if isfield(SEG,'AutoLongAxis')
    SET(NO).AutoLongaxis = SEG.AutoLongaxis;
  end;
  
  if isfield(SEG,'UseLight')
      mydisp('Ignoring UseLight on imported segmentation');
    %DATA.UseLight = SEG.UseLight;
    %set(DATA.Handles.uselightcheckbox,'value',SEG.UseLight);
  end;
  
  if isfield(SEG,'Scar')
    if not(DATA.Silent)
      mywarning('Contains scar data, loading not 100% supported yet.',DATA.GUI.Segment);
    end;
    SET(NO).Scar = SEG.Scar;
    SET(NO).Scar.Mode = 'manualthreshold';
    if not(isfield(SET(NO).Scar,'NoReflow'))
      SET(NO).Scar.NoReflow = repmat(uint8(0),size(SET(NO).Scar.Manual));
    end;
  end;
  
  if isfield(SEG,'RVEndoX')
    SET(NO).RVEndoX = SEG.RVEndoX;
    SET(NO).RVEndoY = SEG.RVEndoY;
    SET(NO).RVEpiX = SEG.RVEpiX;
    SET(NO).RVEpiY = SEG.RVEpiY;
  end;
  
  SET(NO).Longaxis = SEG.Longaxis;
    
  if ~isempty(SET(NO).EndoPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).EndoPinX{tloop,zloop} = SET(NO).EndoPinX{tloop,zloop}+xofs;
        SET(NO).EndoPinY{tloop,zloop} = SET(NO).EndoPinY{tloop,zloop}+yofs;
      end;
    end;
  end;
  
  if ~isempty(SET(NO).EpiPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).EpiPinX{tloop,zloop} = SET(NO).EpiPinX{tloop,zloop}+xofs;
        SET(NO).EpiPinY{tloop,zloop} = SET(NO).EpiPinY{tloop,zloop}+yofs;
      end;
    end;
  end;

  if ~isempty(SET(NO).RVEndoPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).RVEndoPinX{tloop,zloop} = SET(NO).RVEndoPinX{tloop,zloop}+xofs;
        SET(NO).RVEndoPinY{tloop,zloop} = SET(NO).RVEndoPinY{tloop,zloop}+yofs;
      end;
    end;
  end;

  if ~isempty(SET(NO).RVEpiPinX)
    for tloop=1:SET(NO).TSize
      for zloop=1:SET(NO).ZSize
        SET(NO).RVEpiPinX{tloop,zloop} = SET(NO).RVEpiPinX{tloop,zloop}+xofs;
        SET(NO).RVEpiPinY{tloop,zloop} = SET(NO).RVEpiPinY{tloop,zloop}+yofs;
      end;
    end;
  end;
  
  if isfield(SEG,'Rotated')
    SET(NO).Rotated = SEG.Rotated;
  end;
  
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
  end;
  
  if isfield(SEG,'Measure')
    SET(NO).Measure = SEG.Measure;
    for loop=1:length(SET(NO).Measure)
      SET(NO).Measure(loop).X = SET(NO).Measure(loop).X+xofs;
      SET(NO).Measure(loop).Y = SET(NO).Measure(loop).Y+xofs;
    end;
  end;
  
  if isfield(SEG,'Point')
    SET(NO).Point = SEG.Point;
    SET(NO).Point.X = SET(NO).Point.X+xofs;
    SET(NO).Point.Y = SET(NO).Point.Y+yofs;    
  end;
  
  if isfield(SEG,'GEVENCSCALE')
    SET(NO).GEVENCSCALE = SEG.GEVENCSCALE;
  end;
  
catch me
  myfailed('Something went wrong during loading. Undoing...',DATA.GUI.Segment);
  mydispexception(me);
  tools('undosegmentation_Callback');
end;

if size(SET(NO).EndoDraged,1)~=SET(NO).TSize
  SET(NO).EndoDraged = repmat(SET(NO).EndoDraged(:),1,SET(NO).TSize);
  SET(NO).EpiDraged = repmat(SET(NO).EpiDraged(:),1,SET(NO).TSize);
end;

%Prevent problems from old .seg files
segment('checkconsistency',1:SET(NO).TSize,1:SET(NO).ZSize);

segment('updatemodeldisplay');
segment('updatevolume');
segment('viewrefresh_Callback');

%---------------------
function quit_Callback
%---------------------
%Quit Segment
global DATA SET NO

if ~isempty(DATA)
  if DATA.quit
    saveguiposition(DATA.GUI.Segment);
    segment('saveguipositiontodisk');

    fig=DATA.fig;
    DATA = [];
    SET = [];
    NO = [];
    clear DATA
    clear SET
    clear NO
    delete(fig); %Call destroyer to take it down nice, see code =>
    close('all'); %Added do not know if it works.
  end
else
  delete(gcbo);
end