function varargout = segpref(varargin)
%SEGPREF Helper function for SEGMENT, GUI for preferences

%Einar Heiberg

global DATA 

if nargin == 0  % LAUNCH GUI
  
  DATA.setprefhandles;
  update;
  fig = DATA.PrefHandles.fig;
  translation.translatealllabels(fig);
  if nargout > 0
    varargout{1} = fig;
  end
  
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
  
  try
    macro_helper(varargin{:});
    if (nargout)
      [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    else
      feval(varargin{:}); % FEVAL switchyard
    end
  catch me
    mydispexception(me);
  end
  
end

%--------------
function update
%--------------
global DATA

DATA.updateprefhandles;

%---------------------
function undo_Callback %#ok<DEFNU>
%---------------------
%Loads from disk and graphical update
loadpreferences;
if isopengui('segpref.fig')
  update;
end
if isopengui('segprefadvanced.fig')
  updateadvanced;
end

% %----------------------------
% function init_userlogging
% %----------------------------
% global DATA
% val = mygetvalue(DATA.PrefHandles.userloggingcheckbox)
% 
% if val == 1
%       %Initiate log file. Overloaded in most other GUI's
%     pathname =  DATA.Pref.userpath;
%     DATA.UserFile = [pathname filesep sprintf('segmentlog_%s.log',datestr(now,'yyyymmddHHMMSS'))];
%     fid = fopen(g.LogFile,'w');
%     if isequal(fid,-1)
%       myfailed(dprintf('Could not create .log file %s.',g.LogFile));
%     else
%       fclose(fid); %Close the file, to make it empty
%       DATA.startlog(DATA.UserFile); %Start diary.
%     end;
% end


%-----------------------------
function userlogging_Callback %#ok<DEFNU>
%-----------------------------
%started the user logging

global DATA
val = mygetvalue(DATA.PrefHandlesAdvanced.userloggingcheckbox);

% [~,streng] = system('net localgroup Admin');
% [~,strswe] = system('net localgroup Administratörer');
usrname = getenv('USERNAME');
% if isempty(regexp(streng,usrname,'once')) && isempty(regexp(strswe,usrname,'once'))
%   myfailed('You do not have administrator priviliges')
%   return
% end

%Display message about potentially require admin rights.
adminrequirement;

if val == 1
  answer = yesno(sprintf('You are currently logged in as %s. Do you wish to proceed?',usrname));
else
  str = sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),usrname,'User logging stopped.','-');
  DATA.adduserevent(str);
  set(DATA.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','off');
  DATA.Pref.UserLogging = 0;
  return
end

if answer~=1
  DATA.Pref.UserLogging = 0;
  set(DATA.PrefHandlesAdvanced.userloggingcheckbox,'Value',0);
  set(DATA.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','off');
  return
end

DATA.Pref.UserLogging = 1;
if val == 1
  %set path for storing current user.
  if isempty(DATA.Pref.UserLogPath)
    abort = setuserlogpath;
    if abort
      DATA.Pref.UserLogging = 0;
      return
    end
  end

  DATA.init_userlogging
  if DATA.Pref.UserLogging
    set(DATA.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','on');
  else
    set(DATA.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','off');
  end
  
else
  DATA.Pref.UserLogging = 0;
  set(DATA.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','off');
end

%------------------------------
function abort = setuserlogpath 
%------------------------------
% set the path for storing user logging

global DATA
abort=0;
%set path for storing current user.
temp = myuigetdir(pwd,'Select folder for storing user information');
if temp == 0
  mydisp('aborted by user')
  abort = 1;
  return
end
DATA.Pref.UserLogPath=[temp filesep sprintf('segmentuserlog.txt')];

fail = DATA.adduserevent(sprintf('%s\t%s\t%s\t%s','DATE','USER','EVENT','PATIENT ID'));

if fail
  set(DATA.PrefHandlesAdvanced.userloggingcheckbox,'Value',DATA.Pref.UserLogging);
  set(DATA.PrefHandlesAdvanced.userlogpathtext,'String',DATA.Pref.UserLogPath);
  set(DATA.PrefHandlesAdvanced.setuserlogpathpushbutton,'Enable','off');
  return
end

DATA.adduserevent(sprintf('%s\t%s\t%s\t%s',datestr(now,'yyyy-mm-dd HH:MM'),getenv('USERNAME'),sprintf('Changed log path to %s',DATA.Pref.UserLogPath),'-'));
set(DATA.PrefHandlesAdvanced.userlogpathtext,'String',DATA.Pref.UserLogPath);
 %set(DATA.PrefHandles.setuserlogpathstring,DATA.Pref.userpath)

%----------------------------
function setbackgroundcolor(backgroundcolor) %#ok<DEFNU>
%----------------------------
%set the background color for the main interface

global DATA SET

if nargin == 0
  mode = mygetvalue(DATA.PrefHandles.backgroundcolorpopupmenu);
  switch mode
    case 1 %Light gray
      backgroundcolor=[0.94,0.94,0.94];
    case 2
      backgroundcolor=[0.34,0.34,0.34];
    case 3
      backgroundcolor=[0 0 0];
  end
else
  %this is default if not used new addon.
  mode = 1;
  
  if all(backgroundcolor == [0.94,0.94,0.94])
    mode=1;
  end
  
  if all(backgroundcolor == [0.34,0.34,0.34])
    mode=2;
  end
  
  if all(backgroundcolor == [0,0,0])
    mode=3;
  end
end
DATA.Pref.GUIBackgroundColor = backgroundcolor;
DATA.GUISettings.BackgroundColor=DATA.Pref.GUIBackgroundColor;
DATA.GUISettings.BoxAxesColor=[0.1 0.1 0.1];
if mode == 3
  slidercolor=[0.5 0.5 0.5];
else
  slidercolor=[0.86 0.86 0.86];
end

set(DATA.fig,'Color',backgroundcolor);
set(DATA.Handles.barpanel,'BackgroundColor',backgroundcolor)
set(DATA.Handles.distancetext,'BackgroundColor',backgroundcolor)

if mode==3
  set(DATA.Handles.distancetext,'Foreground',[1,1,1])
  set(get(DATA.Handles.timebaraxes,'Xlabel'),'Color','white');
  set(DATA.Handles.timebaraxes,'XColor','white');
  DATA.GUISettings.TimebarAxesColor='white';
else
  DATA.GUISettings.TimebarAxesColor='black';
  set(get(DATA.Handles.timebaraxes,'Xlabel'),'Color','black');
  set(DATA.Handles.timebaraxes,'XColor','black');
end
set(DATA.Handles.flowuipanel,'BackgroundColor',backgroundcolor)
set(DATA.Handles.lvuipanel,'BackgroundColor',backgroundcolor)
set(DATA.Handles.measurementuipanel,'BackgroundColor',backgroundcolor)
set(DATA.Handles.reportpanel,'BackgroundColor',backgroundcolor)
set(DATA.Handles.slider1text,'BackgroundColor',backgroundcolor)
set(DATA.Handles.slider2text,'BackgroundColor',backgroundcolor)

if mode == 3
  set(DATA.Handles.volumeaxes,'Color',[0.34,0.34,0.34])
  set(DATA.Handles.flowaxes,'Color',[0.34,0.34,0.34])
  set(DATA.Handles.lvuipanel,'HighlightColor',[1,1,1])
  set(DATA.Handles.flowuipanel,'HighlightColor',[1,1,1])
  set(DATA.Handles.measurementuipanel,'HighlightColor',[1,1,1])
  set(DATA.Handles.timebaraxes,'Color',[0.34,0.34,0.34])
  set(DATA.Handles.timebaraxes,...
    'XColor',[1,1,1],...
    'YColor',[1,1,1]);
  DATA.GUISettings.VolumeColorGraph=[0.34,0.34,0.34];
  DATA.GUISettings.TimebarAxesColor=[1,1,1];
  DATA.GUISettings.VolumeAxesColor=[1,1,1];
  DATA.AxesTables.volume.fontcolor = [1 1 1]; % [1 1 1];%
  DATA.AxesTables.flow.fontcolor = [1 1 1];
  DATA.AxesTables.measurement.fontcolor = [1 1 1];
  DATA.GUISettings.BarColor=[1,0.6,0.2]; %brandgul
  foregroundcolor=[1,1,1];
  buttonbackgroundcolor=[0.34,0.34,0.34];
else
  set(DATA.Handles.volumeaxes,'Color',[0,0,0])
  set(DATA.Handles.flowaxes,'Color',[0,0,0])
  set(DATA.Handles.lvuipanel,'HighlightColor',[0,0,0])
  set(DATA.Handles.flowuipanel,'HighlightColor',[0,0,0])
  set(DATA.Handles.measurementuipanel,'HighlightColor',[0,0,0])
  DATA.GUISettings.VolumeAxesColor=[0,0,0];
  set(DATA.Handles.timebaraxes,'Color',[1,1,1])
  set(DATA.Handles.timebaraxes,...
    'XColor',[0,0,0],...
    'YColor',[0,0,0]);
  DATA.GUISettings.TimebarAxesColor=[0,0,0];
  DATA.GUISettings.VolumeColorGraph=[1,1,1];
  DATA.GUISettings.BarColor=DATA.GUISettings.BarColorDefault;
  DATA.AxesTables.volume.fontcolor = [0 0 0]; % [1 1 1];%
  DATA.AxesTables.flow.fontcolor = [0 0 0];
  DATA.AxesTables.measurement.fontcolor = [0 0 0];
  foregroundcolor=[0,0,0];
  buttonbackgroundcolor=backgroundcolor;
end

DATA.updateflowaxes;
DATA.updatevolumeaxes;
if ~isempty(SET)
  DATA.updatetimebaraxes;
end

if ~strcmp(DATA.ProgramName,'Segment')
  set(DATA.Handles.lvstackpushbutton,'BackgroundColor',buttonbackgroundcolor,'foregroundcolor',foregroundcolor)
  set(DATA.Handles.rvstackpushbutton,'BackgroundColor',buttonbackgroundcolor,'foregroundcolor',foregroundcolor)
  set(DATA.Handles.flowstackpushbutton,'BackgroundColor',buttonbackgroundcolor,'foregroundcolor',foregroundcolor)
end

set(DATA.Handles.slider1text,'ForegroundColor',foregroundcolor)
set(DATA.Handles.slider2text,'ForegroundColor',foregroundcolor)
 set(DATA.Handles.thumbnailslider,'backgroundcolor',slidercolor);
 
%LV and RV report table
DATA.AxesTables.volume.backgroundcolor = backgroundcolor;%[0.94 0.94 0.94]; %[0 0 0];%
DATA.AxesTables.flow.backgroundcolor = backgroundcolor;%[0.94 0.94 0.94]; %[0 0 0];%
DATA.AxesTables.measurement.backgroundcolor = backgroundcolor;%[0.94 0.94 0.94]; %[0 0 0];%

DATA.AxesTables.volume.draw
DATA.AxesTables.flow.draw
DATA.AxesTables.measurement.draw

%------------------------
function default_Callback %#ok<DEFNU>
%------------------------
%Sets default preferences called from GUI
global DATA

DATA.defaultpref;
update;

%---------------------
function save_Callback(silent) 
%---------------------
%Save preferences to disk
%
%When changing in this function also make the same changes in
%loadpreferences/savetodisk function.

if nargin < 1
  silent = false;
end

global DATA

pathname = getpreferencespath;
pathnamesaveall = pwd; %Segment folder
Pref = DATA.Pref; %#ok<NASGU> %Saved to file
%check if preferences for all users already exist
if exist([pathnamesaveall filesep 'default_preferences.mat']) %,'Pref', DATA.Pref.SaveVersion);
  if not(silent)
    myfailed('Default preferences for all users exists and thereby used prior to local preferences.',DATA.GUI.Segment);

    %lägg till att läsa default preferences och skriva över lokala.
    %Uppdatera gui.
  end
end

try
  save([pathname filesep '.segment_preferences.mat'],'Pref', DATA.Pref.SaveVersion);
catch %#ok<CTCH>
  myfailed('Could not save preferences. Write permission? Disk full?',DATA.PrefHandles.fig);
  return;
end;

disp('Preferences saved.');

%---------------------------
function savetoall_Callback %#ok<DEFNU>
%---------------------------
%Save preferences for all users. Saves to the file defaul_preferences.mat.
%If this file exists, then the PACS and Server settings are copied (and
%overwritten) for the user upon loading.

global DATA

%Display message about potentially require admin rights.
adminrequirement

Pref = DATA.Pref; %#ok<NASGU> %Saved to file

pathname = pwd;
try
  save([pathname filesep 'default_preferences.mat'],'Pref', DATA.Pref.SaveVersion);
  disp('Preferences saved to all users (file default preferences).');
catch
  myfailed('Could not save preferences. Write permission? Disk full?',DATA.PrefHandles.fig);
end;

%--------------------------------
function closepushbutton_Callback %#ok<DEFNU>
%--------------------------------
global DATA

close(DATA.PrefHandles.fig);
DATA.PrefHandles = [];

if DATA.DataLoaded
  segment('viewinterp_Callback',DATA.Pref.ViewInterpolated);
  segment('thumbnailslider_Callback');
  drawfunctions('drawthumbnails');
  segment('viewrefreshall_Callback');

end

%---------------------------
function datapath_Callback %#ok<DEFNU>
%---------------------------
global DATA

if exist(DATA.Pref.datapath,'dir')
  temp = myuigetdir(DATA.Pref.datapath,'Select folder for data');
else
  temp = myuigetdir(pwd,'Select folder for data');
end;
if ~(isempty(temp)||isequal(temp,0))
  DATA.Pref.datapath = temp;
  DATA.Preview.PathName = temp;
end;

update;
  
%---------------------------
function exportpath_Callback %#ok<DEFNU>
%---------------------------
global DATA

if exist(DATA.Pref.exportpath,'dir')
  temp = myuigetdir(DATA.Pref.exportpath,'Select folder for exporting');
else
  temp = myuigetdir(pwd,'Select folder for exporting');
end;

if ~(isempty(temp)||isequal(temp,0))
  DATA.Pref.exportpath = temp;
end;

update;

%-----------------------
function cdpath_Callback %#ok<DEFNU>
%-----------------------
global DATA

temp = myuigetdir(pwd,'Select drive letter for CD or mount location:');

if ~(isempty(temp)||isequal(temp,0))
  DATA.Pref.CDPath = temp;
else
  myfailed('Aborted.',DATA.PrefHandles.fig);
end
  
update;

%-------------------------------------
function openGL_Callback
%-------------------------------------

global DATA

val = get(DATA.PrefHandles.openGLcheckbox,'Value');
if val == 0
  disp('Selecting hardware openGL.');
  opengl hardware
else
  disp('Selecting software openGL.');
  opengl software
end

%Write to log-file
opengl info

DATA.Pref.OpenGLSoftware = val;

update;

%---------------------------
function tempfolder_Callback %#ok<DEFNU>
%---------------------------
%Sets temporaryfolder used by Segment Server and PACS connection (subfolder tempsearch) 

global DATA

if ~isempty(DATA.Pref.Pacs.ImageBasePath)
  if isdir([DATA.Pref.Pacs.ImageBasePath filesep 'TEMP'])
    suggestedfolder = [DATA.Pref.Pacs.ImageBasePath filesep 'TEMP'];    
  else
    suggestedfolder = DATA.Pref.Pacs.ImageBasePath;
  end
else
  suggestedfolder = pwd;
end
folder = myuigetdir(suggestedfolder,'Select folder for temporary storage for patient database and PACS connection. Recommended to keep on SSD disc');
  
DATA.Pref.Pacs.TempStoragePath = folder;
  
updateadvanced;

%-------------------------------------
function imagebasepath_Callback(folder) %#ok<DEFNU>
%-------------------------------------
%Set folder for patientdatabase. If called with folder then do not ask and
%no update.

global DATA

if nargin==0
  if ~isempty(DATA.Pref.Pacs.ImageBasePath)
    suggestedfolder = DATA.Pref.Pacs.ImageBasePath;
  else
    suggestedfolder = pwd;
  end
    folder = myuigetdir(suggestedfolder,'Select folder for patient database and PACS connection. Recommended to keep on SSD disc');
end;

if ~(isempty(folder)||isequal(folder,0))
  
  DATA.Pref.Pacs.ImageBasePath = folder;
  
  %Check if ends with filesep, then remove. This happens when for instance
  %selecting root folder, i.e M:\
  if isequal(folder(end),filesep)
    folder = folder(1:(end-1));
  end;
  
  %Assign other folders
  DATA.Pref.Pacs.AnalysedPath = [folder filesep 'Analysed'];  
  DATA.Pref.Pacs.DicomPath = [folder filesep 'DICOM'];
  
end;

if nargin==0
  updateadvanced;
end;

reportpath_Callback
tempfolder_Callback;

%------------------------------
function reportpath_Callback 
%------------------------------
global DATA

if ~isempty(DATA.Pref.Pacs.ImageBasePath)
  suggestedfolder = DATA.Pref.Pacs.ReportsheetPath;
else
  suggestedfolder = pwd;
end
temp = myuigetdir(suggestedfolder,'Select base folder for patient reports');

if ~(isempty(temp)||isequal(temp,0))
  DATA.Pref.Pacs.ReportsheetPath = temp;
end;

updateadvanced;

%------------------------------
function pafreportpath_Callback %#ok<DEFNU>
%------------------------------
global DATA

temp = myuigetdir([pwd filesep 'PAF'],'Select base folder for PAF report');

if ~(isempty(temp)||isequal(temp,0))
  DATA.Pref.Pacs.PAFPathname = temp;
end;

updateadvanced;

%-------------------------------
function copyexportpath_Callback %#ok<DEFNU>
%-------------------------------
global DATA

DATA.Pref.exportpath = DATA.Pref.datapath;

update;

%---------------------------------
function autosavecheckbox_Callback %#ok<DEFNU>
%---------------------------------
global DATA

DATA.Pref.AutoSave = get(DATA.PrefHandles.autosavecheckbox,'Value');

%-------------------------------
function anonymcheckbox_Callback %#ok<DEFNU>
%-------------------------------
global DATA

DATA.Pref.AnonymMode = get(DATA.PrefHandles.anonymcheckbox,'Value');

if DATA.Pref.AnonymMode
  mymsgbox('This does only blind data when displayed. No changes are made to files.');  
end;

%--------------------------------------
function endocenterradiobutton_Callback %#ok<DEFNU>
%--------------------------------------
global DATA

DATA.Pref.EndoCenter = 1;
set(DATA.PrefHandles.endocenterradiobutton,'Value',1);
set(DATA.PrefHandles.epicenterradiobutton,'Value',0);

%-------------------------------------
function epicenterradiobutton_Callback %#ok<DEFNU>
%-------------------------------------
global DATA

DATA.Pref.EndoCenter = 0;
set(DATA.PrefHandles.endocenterradiobutton,'Value',0);
set(DATA.PrefHandles.epicenterradiobutton,'Value',1);

%-----------------------------------
function blackwhitecheckbox_Callback %#ok<DEFNU>
%-----------------------------------
global DATA

DATA.Pref.BlackWhite = get(DATA.PrefHandles.blackwhitecheckbox,'Value');

%-------------------------------------
function addpointscheckbox_Callback %#ok<DEFNU>
%-------------------------------------
global DATA

DATA.Pref.AddPoints = get(DATA.PrefHandles.addpointscheckbox,'Value');

%------------------------------
function linewidthedit_Callback %#ok<DEFNU>
%------------------------------
global DATA

stri = mygetedit(DATA.PrefHandles.linewidthedit);
[linewidth,ok] = str2num(stri);

if not(ok)
  myfailed('Invalid line width.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.linewidthedit,'String',sprintf('%0.5g',DATA.Pref.LineWidth));
  return;
end;

if linewidth<0.25
  myfailed('Too small line width, minimum 0.25.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.linewidthedit,'String',sprintf('%0.5g',DATA.Pref.LineWidth));
  return;  
end;

if linewidth>5
  myfailed('Too large line width, maximum 5.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.linewidthedit,'String',sprintf('%0.5g',DATA.Pref.LineWidth));
  return;  
end;

%Store & update
DATA.Pref.LineWidth = linewidth;
set(DATA.PrefHandles.linewidthedit,'String',sprintf('%0.5g',DATA.Pref.LineWidth));

%------------------------------
function markersizeedit_Callback %#ok<DEFNU>
%------------------------------
global DATA

stri = mygetedit(DATA.PrefHandles.markersizeedit);
[markersize,ok] = str2num(stri);

if not(ok)
  myfailed('Invalid marker size.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.markersizeedit,'String',sprintf('%0.5g',DATA.Pref.MarkerSize));
  return;
end;

if markersize<1
  myfailed('Too small marker size, minimum 1.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.markersizeedit,'String',sprintf('%0.5g',DATA.Pref.MarkerSize));
  return;  
end;

if markersize>20
  myfailed('Too large marker size, maximum 20.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.markersizeedit,'String',sprintf('%0.5g',DATA.Pref.MarkerSize));
  return;  
end;

%Store & update
DATA.Pref.MarkerSize = markersize;
set(DATA.PrefHandles.markersizeedit,'String',sprintf('%0.5g',DATA.Pref.MarkerSize));

%------------------------------
function numberthumbnailsedit_Callback %#ok<DEFNU>
%------------------------------
global DATA

stri = mygetedit(DATA.PrefHandles.numberthumbnailsedit);
[numberthumbnails,ok] = str2num(stri);

if not(ok)
  myfailed('Invalid number of thumbnails.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.numberthumbnailsedit,'String',sprintf('%0.5g',DATA.Pref.NumberVisibleThumbnails));
  return;
end;

if numberthumbnails<1
  myfailed('At least one thumbnail must be visible.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.numberthumbnailsedit,'String',sprintf('%0.5g',DATA.Pref.NumberVisibleThumbnails));
  return;  
end;

if numberthumbnails>20
  myfailed('A maximum of 20 thumbnails can be visible.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.numberthumbnailsedit,'String',sprintf('%0.5g',DATA.Pref.NumberVisibleThumbnails));
  return;  
end;

%Store & update
DATA.Pref.NumberVisibleThumbnails = numberthumbnails;
set(DATA.PrefHandles.numberthumbnailsedit,'String',sprintf('%0.5g',DATA.Pref.NumberVisibleThumbnails));

%-----------------------------------
function radialprofilesedit_Callback %#ok<DEFNU>
%-----------------------------------
global DATA

stri = mygetedit(DATA.PrefHandles.radialprofilesedit);

[radialprofiles,ok] = str2num(stri);
if not(ok)
  myfailed('Invalid number of points.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.radialprofilesedit,'String',sprintf('%d',DATA.Pref.RadialProfiles));
  return;
end;

radialprofiles = round(radialprofiles);

if radialprofiles<80
  myfailed('Can not be smaller than 80.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.radialprofilesedit,'String',sprintf('%d',DATA.Pref.RadialProfiles));
  return;
end;

if radialprofiles>1000
  myfailed('Can not be larger than 1000.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.radialprofilesedit,'String',sprintf('%d',DATA.Pref.RadialProfiles));
  return;
end;

%Store & update
DATA.Pref.RadialProfiles = radialprofiles;
set(DATA.PrefHandles.radialprofilesedit,'String',sprintf('%d',DATA.Pref.RadialProfiles));

%---------------------------------
function interppointsedit_Callback %#ok<DEFNU>
%---------------------------------
global DATA 

stri = mygetedit(DATA.PrefHandles.interppointsedit);

[interppoints,ok] = str2num(stri);
if not(ok)
  myfailed('Invalid number of points.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.interppointsedit,'String',sprintf('%d',DATA.Pref.NumInterpPoints));
  return;
end;

interppoints = round(interppoints);
DATA.Pref.NumInterpPoints = interppoints;
set(DATA.PrefHandles.interppointsedit,'String',sprintf('%d',DATA.Pref.NumInterpPoints));

%------------------------------
function numpointsedit_Callback  %#ok<DEFNU>
%------------------------------
global DATA SET

stri = mygetedit(DATA.PrefHandles.numpointsedit);

[numpoints,ok] = str2num(stri);
if not(ok)
  myfailed('Invalid number of points.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.numpointsedit,'String',sprintf('%d',DATA.Pref.NumPoints));
  return;
end;

numpoints = round(numpoints);

if numpoints<80
  myfailed('Can not be smaller than 80.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.numpointsedit,'String',sprintf('%d',DATA.Pref.NumPoints));
  return;
end;

if numpoints>1000
  myfailed('Can not be larger than 1000.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.numpointsedit,'String',sprintf('%d',DATA.Pref.NumPoints));
  return;
end;
%Store & update
DATA.Pref.NumPoints = numpoints;
set(DATA.PrefHandles.numpointsedit,'String',sprintf('%d',DATA.Pref.NumPoints));
DATA.NumPoints = DATA.Pref.NumPoints;
set(DATA.PrefHandles.numpointstext,'String',sprintf('%d',DATA.NumPoints));
h = mywaitbarstart(length(SET),'Please wait, resampling curves.');

for loop=1:length(SET)
  numpointsedithelper(loop);
  segment('updatemodeldisplay',loop);
  h = mywaitbarupdate(h);
end;
mywaitbarclose(h);

%---------------------------
function numpointsedithelper(loop)
%---------------------------
global DATA SET

if ~isempty(SET(loop).EndoX)
  [SET(loop).EndoX,SET(loop).EndoY] = calcfunctions('resamplemodel',...
    SET(loop).EndoX,SET(loop).EndoY,DATA.NumPoints);
end;
if ~isempty(SET(loop).EpiX)
  [SET(loop).EpiX,SET(loop).EpiY] = calcfunctions('resamplemodel',...
    SET(loop).EpiX,SET(loop).EpiY,DATA.NumPoints);
end;
if ~isempty(SET(loop).RVEndoX)
  [SET(loop).RVEndoX,SET(loop).RVEndoY] = calcfunctions('resamplemodel',...
    SET(loop).RVEndoX,SET(loop).RVEndoY,DATA.NumPoints);
end;
if ~isempty(SET(loop).RVEpiX)
  [SET(loop).RVEpiX,SET(loop).RVEpiY] = calcfunctions('resamplemodel',...
    SET(loop).RVEpiX,SET(loop).RVEpiY,DATA.NumPoints);
end;
for rloop=1:SET(loop).RoiN
  [SET(loop).Roi(rloop).X,SET(loop).Roi(rloop).Y] = calcfunctions('resamplemodel',...
    SET(loop).Roi(rloop).X,SET(loop).Roi(rloop).Y,DATA.NumPoints);
end

% if isopengui('strain.fig')
%   mywarning('Strain analysis needs to be redone after changing Number of points along contour. Closing Strain GUI.')
%   straintagging.straintagging('close_Callback')
% end

%------------------------------------------
function contouradjustdistanceedit_Callback %#ok<DEFNU>
%------------------------------------------
global DATA

stri = mygetedit(DATA.PrefHandles.contouradjustdistanceedit);

[dist,ok] = str2num(stri);
if not(ok)
  myfailed('Invalid distance.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.contouradjustdistanceedit,'String',...
    sprintf('%0.5g',DATA.Pref.ContourAdjustDistance));
  return;
end;

if dist<0
  myfailed('Distance may not be negative.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.contouradjustdistanceedit,'String',...
    sprintf('%0.5g',DATA.Pref.ContourAdjustDistance));  
  return;
end;

if dist>40
  myfailed('Too large distance, max 40.',DATA.PrefHandles.fig);
  set(DATA.PrefHandles.contouradjustdistanceedit,'String',...
    sprintf('%0.5g',DATA.Pref.ContourAdjustDistance));  
  return;
end
  
DATA.Pref.ContourAdjustDistance = dist;
set(DATA.PrefHandles.contouradjustdistanceedit,'String',...
    sprintf('%0.5g',DATA.Pref.ContourAdjustDistance));  

%---------------------------------------
function includeallpixelsinroi_Callback %#ok<DEFNU>
%---------------------------------------
global DATA

DATA.Pref.IncludeAllPixelsInRoi = get(DATA.PrefHandles.includeallpixelsinroicheckbox,'Value');

%-------------------------
function donotask_Callback %#ok<DEFNU>
%-------------------------
global DATA

DATA.Pref.DoNotAsk = get(DATA.PrefHandles.donotaskcheckbox,'value');

%------------------------------
function dicomportedit_Callback %#ok<DEFNU>
%------------------------------
global DATA

temp = mygetedit(DATA.PrefHandlesAdvanced.dicomportedit);
temp = removechars(temp);

DATA.Pref.Server.DICOMPort = temp;
set(DATA.PrefHandlesAdvanced.dicomportedit,'string',temp);

%----------------------------
function aetitleedit_Callback %#ok<DEFNU>
%----------------------------
global DATA

temp = mygetedit(DATA.PrefHandlesAdvanced.aetitleedit);
temp = upper(temp);

DATA.Pref.Server.AETitle = temp;
set(DATA.PrefHandlesAdvanced.aetitleedit,'string',temp);

%-----------------------------------
function receiveoptionsedit_Callback %#ok<DEFNU>
%-----------------------------------
global DATA
stri = mygetedit(DATA.PrefHandlesAdvanced.receiveoptionsedit);
DATA.Pref.Server.ReceiveOptions = stri;

%--------------------------------
function sendoptionsedit_Callback %#ok<DEFNU>
%--------------------------------
global DATA
stri = mygetedit(DATA.PrefHandlesAdvanced.sendoptionsedit);
DATA.Pref.Pacs.SendOptions = stri;

%---------------------------
function switchtags_Callback
%---------------------------
global DATA
DATA.Pref.Pacs.SwitchTags = mygetvalue(DATA.PrefHandlesAdvanced.switchtagscheckbox);

%----------------------------------
function normalizephaseask_Callback %#ok<DEFNU>
%----------------------------------
global DATA
DATA.Pref.Dicom.NormalizePhase=[];
DATA.normalizephaseupdate;

%----------------------------------
function normalizephaseyes_Callback %#ok<DEFNU>
%----------------------------------
global DATA
DATA.Pref.Dicom.NormalizePhase=1;
DATA.normalizephaseupdate;

%---------------------------------
function normalizephaseno_Callback %#ok<DEFNU>
%---------------------------------
global DATA
DATA.Pref.Dicom.NormalizePhase=0;
DATA.normalizephaseupdate;    

%--------------------------------------
function uselightcontrast_Callback %#ok<DEFNU>
%--------------------------------------
global DATA
DATA.BalloonLevel = -1; %will force recalc of balloon image.
%DATA.UseLight = get(DATA.PrefHandles.uselightcheckbox,'value');
DATA.Pref.UseLight = get(DATA.PrefHandles.uselightcheckbox,'value');

%----------------------------------
function viewinterpolated_Callback %#ok<DEFNU>
%----------------------------------
global DATA

DATA.Pref.ViewInterpolated = get(DATA.PrefHandles.viewinterpolatedcheckbox,'value');

%----------------------------------
function bgcolor_Callback %#ok<DEFNU>
%----------------------------------
global DATA

DATA.Pref.BackgroundColor = get(DATA.PrefHandles.bgcolorcheckbox,'value');

%----------------------------------
function force16bit_Callback %#ok<DEFNU>
%----------------------------------
global DATA
DATA.Pref.Dicom.Force16Bit = get(DATA.PrefHandlesAdvanced.force16bitcheckbox,'value');

%-------------------------------
function fasterpreview_Callback %#ok<DEFNU> Callback from segpref.fig
%-------------------------------
global DATA
DATA.Pref.Dicom.FasterPreview = get(DATA.PrefHandlesAdvanced.fasterpreviewcheckbox,'value');

%-------------------------------
function showseriesdescription_Callback %#ok<DEFNU> Callback from segpref.fig
%-------------------------------
global DATA
DATA.Pref.ShowSeriesDescription = ...
            get(DATA.PrefHandles.showseriesdescriptioncheckbox, 'value');
          
%-------------------------------
function hidefilesunix_Callback %#ok<DEFNU> Callback from segpref.fig
%-------------------------------
global DATA
DATA.Pref.HideFilesUnix =get(DATA.PrefHandles.hidefilesunixcheckbox, 'value');

%----------------------------------------
function fastpreviewloadcheckbox_Callback %#ok<DEFNU> Callback from segpref.fig
%----------------------------------------
global DATA
DATA.Pref.FastPreviewLoad = get(DATA.PrefHandles.fastpreviewloadcheckbox,'value');

%-------------------------------------
function checkversioncheckbox_Callback %#ok<DEFNU>
%-------------------------------------
global DATA
DATA.Pref.CheckVersion = get(DATA.PrefHandles.checkversioncheckbox,'value');

%-------------------------------------
function useproxyserver_Callback(value) 
%-------------------------------------
global DATA

if nargin==0
  value = get(DATA.PrefHandles.useproxyservercheckbox,'value');  
end;

DATA.Pref.UseProxyServer = value;
set(DATA.PrefHandles.useproxyservercheckbox,'value',value);

%Check if not properly defined
if DATA.Pref.UseProxyServer  
  if isequal(DATA.Pref.ProxySettings.HostName,'')
    myfailed('No proxy server configured. Aborted.');
    defineproxyserver_Callback;
    DATA.Pref.UseProxyServer = false;
    set(DATA.PrefHandles.useproxyservercheckbox,'value',0);
  end;
end;

if ~DATA.Pref.UseProxyServer
  disp('No proxy server used');
end;

useproxyserver(DATA.Pref.UseProxyServer);

%----------------------------------
function defineproxyserver_Callback 
%----------------------------------
global DATA

prompt = {...
  'HostName',...
  'Port',...
  'UserName',...
  'Password'};

def = {...
  DATA.Pref.ProxySettings.HostName ...
  DATA.Pref.ProxySettings.Port ...
  DATA.Pref.ProxySettings.UserName ...
  DATA.Pref.ProxySettings.Password };

res = inputdlg(prompt,'Enter data',1,def);

if not(length(res)==length(prompt))
  myfailed('Aborted.');
  return;
end;

DATA.Pref.ProxySettings.HostName = res{1};
DATA.Pref.ProxySettings.Port = res{2};
DATA.Pref.ProxySettings.UserName = res{3};
DATA.Pref.ProxySettings.Password  = res{4};

if isequal(DATA.Pref.ProxySettings.HostName,'') || isempty(DATA.Pref.ProxySettings.HostName)
  disp('No proxy server used');
  useproxyserver_Callback(false);  
else
  useproxyserver_Callback(true);  
end;

%-----------------------------------
function language_Callback(language) %#ok<DEFNU>
%-----------------------------------
global DATA

if isfield(DATA.Handles,'englishmenu')
  set(DATA.Handles.englishmenu,'Checked','off');
end
if isfield(DATA.Handles,'svenskamenu')
  set(DATA.Handles.svenskamenu,'Checked','off');
end
% if isfield(DATA.Handles,'italianomenu')
%   set(DATA.Handles.italianomenu,'Checked','off');
% end
if isfield(DATA.Handles,'deutschmenu')
  set(DATA.Handles.deutschmenu,'Checked','off');
end
% set([DATA.Handles.englishmenu DATA.Handles.svenskamenu DATA.Handles.italianomenu DATA.Handles.deutschmenu],'Checked','off');
try
  eval(sprintf('set(DATA.Handles.%smenu,''Checked'',''on'')',lower(language)));
catch % if the current language does not exist in the current software
  DATA.Pref.Language = 'English';
  language = 'English';
  eval(sprintf('set(DATA.Handles.%smenu,''Checked'',''on'')',lower(language)));
end
prevlanguage = DATA.Pref.Language;
DATA.Pref.Language = language;
if strcmp(prevlanguage,language)
  prevlanguage = 'English';
end

if isfield(DATA.Handles,'toggleiconholder');
  DATA.setribbonimages(language);
end
set(get(DATA.Handles.volumeaxes,'xlabel'),'string',translation.dictionary(...
  'Time [ms]'))%get(get(DATA.Handles.volumeaxes,'xlabel'),'string')))
set(get(DATA.Handles.volumeaxes,'ylabel'),'string',translation.dictionary(...
  'Volume [ml]'))%get(get(DATA.Handles.volumeaxes,'ylabel'),'string')))
set(get(DATA.Handles.flowaxes,'xlabel'),'string',translation.dictionary(...
  'Time [ms]'))%get(get(DATA.Handles.flowaxes,'xlabel'),'string')))
set(get(DATA.Handles.flowaxes,'ylabel'),'string',translation.dictionary(...
  'Flow [ml/s]'))%get(get(DATA.Handles.flowaxes,'ylabel'),'string')))
set(get(DATA.Handles.timebaraxes,'xlabel'),'string',translation.dictionary(...
  'Time [ms]'))%get(get(DATA.Handles.timebaraxes,'xlabel'),'string')))
% set(DATA.Handles.lvstackpushbutton,'String',dprintf('Stack #%d',no));
if isfield(DATA.Handles,'lvstackpushbutton')
  if ~isempty(DATA.LVNO)    
    set(DATA.Handles.lvstackpushbutton,'String',dprintf(translation.dictionary('Stack #%d'),DATA.LVNO));
  else
    set(DATA.Handles.lvstackpushbutton,'String',dprintf(translation.dictionary('Set stack')));
  end
end
if isfield(DATA.Handles,'rvstackpushbutton')
  if ~isempty(DATA.RVNO)    
    set(DATA.Handles.rvstackpushbutton,'String',dprintf(translation.dictionary('Stack #%d'),DATA.RVNO));
  else
    set(DATA.Handles.rvstackpushbutton,'String',dprintf(translation.dictionary('Set stack')));
  end
end
if isfield(DATA.Handles,'flowstackpushbutton')
  if ~isempty(DATA.FlowNO)    
    set(DATA.Handles.flowstackpushbutton,'String',dprintf(translation.dictionary('Stack #%d'),DATA.FlowNO));
  else
    set(DATA.Handles.flowstackpushbutton,'String',dprintf(translation.dictionary('Set stack')));
  end
end


%DATA.flowreportupdate
if isfield(DATA.AxesTables,'flow')
DATA.AxesTables.flow.draw();
%DATA.AxesTables.flow.updateName('Backward','Backward',true);
%DATA.AxesTables.flow.updateName('Forward','Forward',true);
%DATA.AxesTables.flow.updateTitle('Flow',true);
end

guis = fieldnames(DATA.GUI);
for i = 1:numel(guis)
  if isa(DATA.GUI.(guis{i}),'mygui')
    translation.translatealllabels(DATA.GUI.(guis{i}).fig,prevlanguage);
  end
end
silent = true;
save_Callback(silent); %Always save


%---------------------------
function webbrowser_Callback %#ok<DEFNU> Callback from segpref.fig
%---------------------------
global DATA

chosenwebbrowser=get(DATA.PrefHandles.webbrowserpopupmenu,'value');

browserexist=false;
browser='explorer';
browserpath='internal command';
switch chosenwebbrowser
  case 1
    if ispc
      browserexist=true;
      DATA.Pref.WebBrowser='explorer';
    end
  case 2
    browser = 'Mozilla Firefox\firefox.exe';
    browserpath = getenv('ProgramFiles');
    browserexist=exist([browserpath filesep 'Mozilla Firefox'],'dir');
    if browserexist
      browserexist=exist([browserpath filesep 'Mozilla Firefox' filesep 'firefox.exe'],'file');
      if browserexist
        DATA.Pref.WebBrowser=[browserpath filesep browser];
      end
    end    
  case 3
    browserpath = getenv('ProgramFiles');
    browser = 'other';
    [filename,pathname] = myuigetfile('*.*','Select web browser to use',browserpath);
    if (pathname~=0)
      browserexist=true;
      DATA.Pref.WebBrowser=[pathname filename];
    end    
end
if not(browserexist)
  myfailed(dprintf('Web browser %s does not exist in %s. Choose another web browser.',browser, browserpath));
end


%----------------------------------
function advancedsettings_Callback %#ok<DEFNU> called from segpref.fig
%----------------------------------

global DATA

fig = openfig('segprefadvanced.fig','reuse');
set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
translation.translatealllabels(fig);

% Generate a structure of handles to pass to callbacks, and store it.
DATA.PrefHandlesAdvanced = guihandles(fig);
DATA.PrefHandlesAdvanced.fig = fig;

updateadvanced;

%----------------------
function updateadvanced
%----------------------
global DATA

DATA.updateprefhandlesadvanced;

%--------------------------------
function defaultadvanced_Callback %#ok<DEFNU> called from segprefadvanced.fig
%--------------------------------
global DATA

%DICOM Communication Settings
DATA.Pref.Pacs.ImageBasePath = [pwd filesep 'Database'];
DATA.Pref.Pacs.TempStoragePath = [pwd filesep 'Database' filesep 'TempStorage'];
DATA.Pref.Pacs.DicomPath = [pwd filesep 'Database' filesep 'DICOM'];
DATA.Pref.Pacs.AnalysedPath = [pwd filesep 'Database' filesep 'Analysed'];
DATA.Pref.Pacs.ReportsheetPath = [pwd filesep 'Database' filesep 'Report'];
DATA.Pref.Server.DICOMPort = '104';
aetitle = getenv('COMPUTERNAME');
if isempty(aetitle)
  aetitle = 'SEGMENT';
end;
DATA.Pref.Server.AETitle = aetitle;

%DICOM interpretation
DATA.Pref.Dicom.FasterPreview = true;
DATA.Pref.Dicom.Force16Bit = false;
DATA.Pref.Dicom.NormalizePhase = []; 

DATA.Pref.Server.ReceiveOptions = '';
DATA.Pref.Pacs.SendOptions = '';

updateadvanced;

%--------------------------------
function closeadvancedpushbutton_Callback %#ok<DEFNU> %called from segprefadvanced.fig
%--------------------------------
global DATA

close(DATA.PrefHandlesAdvanced.fig);
DATA.PrefHandlesAdvanced = [];

%--------------------------------
function generatebatfile_Callback %#ok<DEFNU>
%--------------------------------
%Function to generate .bat files for server and sorter

global DATA

if ~checkdatabaselicense
  return;
end;

adminrequirement; %message that it may require admin rights.

segmentpath = pwd;

%Create path for logfile if it does not already 
logfilepath = [DATA.Pref.Pacs.TempStoragePath filesep 'logfiles'];
sucess = mymkdir(logfilepath);
if ~sucess
  myfailed('Could not create folder for logfile. Check write permission.');
  return;
end;

%--- Prepare data for bat file
if ispc
  callname = 'storescp.exe';
else
  callname = 'storescp';
end;
outpath = DATA.Pref.Pacs.TempStoragePath;
dicomport = DATA.Pref.Server.DICOMPort;
aetitle = DATA.Pref.Server.AETitle;
receiveoptions = DATA.Pref.Server.ReceiveOptions;
  
%--- Create the segment storage server file
[fid, message] = fopen('segmentstorageserver.bat','wt');

%Check if could create.
if fid<0
  myfailed(sprintf('Could not create file. Detailed message:%s',message))
  return;
end;

%Create string & write
fprintf(fid,'"%s%s%s" -od "%s" %s -aet "%s" %s \n',... 
  segmentpath,... %path to executable, this is segmentpath
  filesep,... %filesep
  callname,... %name of executable
  outpath,... %pathname for output directory (temporary directory)
  receiveoptions,... %receiveoptions
  aetitle,... %aetitle
  dicomport); %port name  

%Close the file
fclose(fid);

%--- Create the segmentsorterbatfile

%Create the segment storage server file
[fid, message] = fopen('segmentsorterserver.bat','wt');

%Check if could create.
if fid<0
  myfailed(sprintf('Could not create file. Detailed message:%s',message))
  return;
end;

%Create string & write
if ispc
  callname = 'segmentserversorter.exe';
else
  callname = 'segmentserversorter';
end;
databasepath = DATA.Pref.Pacs.ImageBasePath;

%segmentserversorter(DATA.Pref.Pacs.TempStoragePath, DATA.Pref.Pacs.ImageBasePath);

fprintf(fid,'"%s%s%s" "%s" "%s"  \n',... 
  segmentpath,... %path to executable, this is segmentpath
  filesep,... %filesep
  callname,... %name of executable
  outpath,... %pathname for output directory (temporary directory)
  databasepath);

%Close the file
fclose(fid);

mymsgbox('Created storage and sorter .bat files.');

%-------------------------------
function installservices_Callback %#ok<DEFNU>
%-------------------------------
%Install the segmentstorageserver & segmentsorterservice as services.

if ~checkdatabaselicense
  return;
end;
forbidden={'segmentserversorter.exe', 'nssm.exe', 'storescu.exe'};

isrunning=zeros(size(forbidden));

for i=length(forbidden)
  [~,result] = system(sprintf('tasklist /FI "imagename eq %s" /fo table /nh',forbidden{i}));
  if ~isempty(regexp(lower(result),forbidden{i}, 'once'))%strcmp(lower(result(1:length(forbidden{i}))),forbidden{i})
    isrunning(i)=1;
  end
end

if any(isrunning)
   answer=yesno(sprintf('To proceed with install, the following services need to shutdown:\n\n%s\n%s\n\n Do this now?','SegmentStorageServer' , 'SegmentSorterService'));
   if answer
     segpref('editservices_Callback')
   else
     return
   end
end

segmentpath = pwd;

%--- StorageServer
storagebatfilename = [pwd filesep 'segmentstorageserver.bat']; 

if ~exist(storagebatfilename,'file')
  myfailed('No storage .bat file exist, please generate.');
  return;
end;

stri = sprintf('"%s%snssm.exe" install SegmentStorageServer "%s"',...
  segmentpath,...
  filesep,...
  storagebatfilename);

[status,result] = system(stri);
statushelper('install service',status,result);
if isequal(status,0)
  mymsgbox('SegmentStorageServer service have been installed. It is not yet running. Either start manually (under Edit), or restart computer.');
end;

%--- SorterServer
sorterbatfilename = [pwd filesep 'segmentsorterserver.bat']; 

if ~exist(sorterbatfilename,'file')
  myfailed('No sorter .bat file exist, please generate.');
  return;
end;

stri = sprintf('"%s%snssm.exe" install SegmentSorterServer "%s"',...
  segmentpath,...
  filesep,...
  sorterbatfilename);

[status,result] = system(stri);
statushelper('install service',status,result);
if isequal(status,0)
  mymsgbox('SegmentSorterServer service have been installed. It is not yet running. Either start manually (under Edit), or restart computer.');
end;

% %isrunning=wasrunning
% if any(isrunning)
%     try
%       tostart=find(isrunning);
%       %To start and stop you need administrator priviliges, is pause sufficient?
%       for i=tostart
%         [~,result] = system(sprintf('net start %s',forbidden{isrunning(i)}));
%       end
%     catch
%       mywarning('Failed to start previously ongoing process, are administrator priviliges on?')
%       return;
%     end
% end

%--------------------------------
function deletestorageservice_Callback %#ok<DEFNU>
%--------------------------------
%Delete SegmentStorageServer as service.

if ~checkdatabaselicense
  return;
end;

mymsgbox('Make sure that the SegmentStorageServer is stopped before deleting.');

segmentpath = pwd;

%StorageService
stri = sprintf('"%s%snssm.exe" remove SegmentStorageServer',...
  segmentpath,...
  filesep);
[status,result] = system(stri);
statushelper('remove service',status,result);

%--------------------------------
function deletesorterservice_Callback %#ok<DEFNU>
%--------------------------------
%Delete SegmentSorterServer as service.

if ~checkdatabaselicense
  return;
end;

mymsgbox('Make sure that the SegmentSorterServer is stopped before deleting.');

segmentpath = pwd;

%SorterService
stri = sprintf('"%s%snssm.exe" remove SegmentSorterServer',...
  segmentpath,...
  filesep);
[status,result] = system(stri);
statushelper('remove service',status,result);

%-----------------------------
function editservices_Callback %#ok<DEFNU>
%-----------------------------
%Starts system editor for services.

if ~checkdatabaselicense
  return;
end;

mymsgbox('Stop/Start SegmentStorageServer and SegmentSorterService. Segment will not be responsive before you close services editor. Now starting services.msc.');
system('services.msc');
mymsgbox('Services editor closed.');

%-----------------------------------
function statushelper(stri,status,result)
%-----------------------------------
%Helper function to generate error messages

if ~isequal(status,0)
  myfailed(sprintf('Could not %s. Message:%s',stri,result));
end;

%------------------------------------
function ok = checkdatabaselicense
%------------------------------------
global DATA

ok = false;

if not(isequal(3,getmodule(2,'D',[],true))) %'' => make license check silent.
  myfailed('Your license does not include this module.',DATA.GUI.Segment);
  return;
end;

ok = true;