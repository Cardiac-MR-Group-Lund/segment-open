function loadpreferences
%Load Segment preferences. If none exist in preferences path, try to move
%preferences.mat from main Segment folder to preferences path. If that does
%not exist, try to copy default_preferences.mat to preferences path. This
%allows the admin to specify preferences to be used every time a new user
%starts using Segment.

global DATA 

%Einar Heiberg

pathname = getpreferencespath;

if ispc
  if not(exist([pathname filesep '.segment_preferences.mat'],'file'))
    %mywarning('Moving system preferences file to new user location.');
    mydisp('Moving system preferences file to new user location.');
    
    %Old old location
    if exist('preferences.mat','file')
      mymovefile([pwd filesep 'preferences.mat'],[pathname filesep '.segment_preferences.mat']);
      %     dos(sprintf('move "%s" "%s"',...
      %       [pwd filesep 'preferences.mat'],...
      %       [pathname filesep '.segment_preferences.mat']));
    elseif exist('default_preferences.mat','file')
      mycopyfile([pwd filesep 'default_preferences.mat'],[pathname filesep '.segment_preferences.mat']);
    end;
    
    %old location just try it...
    mymovefile([pwd filesep '.segment_preferences.mat'],[pathname filesep '.segment_preferences.mat']);
    %   dos(sprintf('move "%s%s..%s%s" "%s%s%s"',...
    %     pathname,filesep,filesep,'.segment_preferences.mat',...
    %     pathname,filesep,'.segment_preferences.mat'));%Det gar inte att hitta filen.
    
    %Remove old log file
    %dos(sprintf('delete "%s%s..%s%s"',pathname,filesep,filesep,'.segment.log'));%delete ar inte ett internt kommando, externt kommando, program eller kommandofil.
  end
end;

%Try to find a more suitable place for the 
try
  load([pathname filesep '.segment_preferences.mat'],'Pref');
catch %#ok<CTCH>
  mydisp('Could not read preferences. Probably a new installation.');
  Pref.datapath = [pwd filesep 'Examples'];
  Pref.exportpath = Pref.datapath;
  Pref.AnonymMode = false;
  Pref.AddPoints = false;
  DATA.Pref = Pref;
  preferencesbackward;
  DATA.defaultpref;
  Pref = DATA.Pref;
  savetodisk; %Save the preferences to disk
end;
DATA.Pref = Pref;
DATA.Pref.AutoSave = false;

%Run backwards compability test
preferencesbackward;

if isa(DATA,'maingui')
  DATA.standardpref;
  DATA.updatepacslabels;
end

% %this part is overwritten if we find a default preferences
% if not(isfield(DATA.Pref,'UserLogging'))
%   DATA.Pref.UserLogging = false;
% end;
% 
% if not(isfield(DATA.Pref,'UserLogPath'))
%   DATA.Pref.UserLogPath = [];
% end;

%Check if we should read in default preferences for PACS, Server etc.
filename = [pwd filesep 'default_preferences.mat'];
if exist(filename,'file')  
  load(filename,'Pref','-mat');
  mydisp('Found default preferences, overriding PACS, DICOM, Server settings, and User logging.');
  mydisp(' ');
  
  %Copy data
  DATA.Pref.Pacs = Pref.Pacs;
  DATA.Pref.Dicom = Pref.Dicom;
  DATA.Pref.Server = Pref.Server;
  DATA.Pref.Gensvar = Pref.Gensvar;
  DATA.Pref.AETitle = Pref.AETitle;
  DATA.Pref.UseProxyServer = Pref.UseProxyServer;
  DATA.Pref.ProxySettings = Pref.ProxySettings;
    
  if isfield(Pref,'UserLogging')
    DATA.Pref.UserLogging = Pref.UserLogging;
  end
  if isfield(Pref,'UserLogPath')
    DATA.Pref.UserLogPath = Pref.UserLogPath;
  end
  
  %Run extra backwards compability test
  preferencesbackward;

end;

%Select OpenGL according to preferences
if DATA.Pref.OpenGLSoftware
  opengl software;
else
  opengl hardware;
end;
disp('OpenGL info:');
opengl info;

%---------------------------
function preferencesbackward
%---------------------------
global DATA

%Moved this up from Pacs clause. It is needed.
if not(isfield(DATA.Pref,'AETitle'))
  aetitle = getenv('COMPUTERNAME');
  if isempty(aetitle)
    aetitle = 'SEGMENT';
  end;
  DATA.Pref.AETitle = aetitle;
end;

%Backwards compability
if not(isfield(DATA.Pref,'Pacs'))
  %These fields can (and have to) be skipped if restruct has been done
  if not(isfield(DATA.Pref,'AnonymMode'))
    DATA.Pref.AnonymMode = false;
  end;
  if not(isfield(DATA.Pref,'EndoCenter'))
    DATA.Pref.EndoCenter = true;
  end;
  if not(isfield(DATA.Pref,'BlackWhite'))
    DATA.Pref.BlackWhite = false;
  end;
  if not(isfield(DATA.Pref,'LineWidth'))
    DATA.Pref.LineWidth = 1;
  end;
  if not(isfield(DATA.Pref,'NumPoints'))
    DATA.Pref.NumPoints = 80;
  end;
  if not(isfield(DATA.Pref,'RadialProfiles'))
    DATA.Pref.RadialProfiles = DATA.Pref.NumPoints;
    DATA.Pref.NumPoints = 80;
  end;
  if not(isfield(DATA.Pref,'LearnMode'))
    DATA.Pref.LearnMode = true;
  end;
  if not(isfield(DATA.Pref,'UndoHistory'))
    DATA.Pref.UndoHistory = 10;
  end;
  if not(isfield(DATA.Pref,'IncludeAllPixelsInRoi'))
    DATA.Pref.IncludeAllPixelsInRoi = false;
  end;
  if not(isfield(DATA.Pref,'AutoSave'))
    DATA.Pref.AutoSave = true;
  end;
  if not(isfield(DATA.Pref,'ContourAdjustDistance'))
    DATA.Pref.ContourAdjustDistance = 2.5;
  end;
  if not(isfield(DATA.Pref,'AllowDicomCache'))
    DATA.Pref.AllowDicomCache = false; %Set to false since obsoleted
  end;
  if not(isfield(DATA.Pref,'DoNotAsk'))
    DATA.Pref.DoNotAsk = false;
  end;
  if not(isfield(DATA.Pref,'DICOMPort'))
    DATA.Pref.DICOMPort = '104';
  end;
  if not(isfield(DATA.Pref,'NormalizePhase'))
    DATA.Pref.NormalizePhase = []; %EH bugfix: was without .Pref.
  end
  if not(isfield(DATA.Pref,'NumberThumbnails'))
    DATA.Pref.NumberVisibleThumbnails = 7;
  end
  if not(isfield(DATA.Pref,'MarkerSize'))
    DATA.Pref.MarkerSize = 5;
  end

  if not(isfield(DATA.Pref,'AETitle'))
    aetitle = getenv('COMPUTERNAME');
    if isempty(aetitle)
      aetitle = 'SEGMENT';
    end;
    DATA.Pref.AETitle = aetitle;
  end;
  
  %added by eriks
  if not(isfield(DATA.Pref,'UseLight'))
    DATA.Pref.UseLight = false;
  end
  
  if not(isfield(DATA.Pref,'ViewInterpolated'))
    DATA.Pref.ViewInterpolated = true;
  end;
  
   if not(isfield(DATA.Pref,'UserLogging'))
    DATA.Pref.UserLogging = false;
  end;
  
   if not(isfield(DATA.Pref,'UserLogPath'))
    DATA.Pref.UserLogPath = [];
  end;
  
  %added by EH:
  if not(isfield(DATA.Pref,'Force16Bit'))
    DATA.Pref.Force16Bit = false;
  end
  
  %added by EH:
  %When true segment uses imread(filename,1) to load only first frame
  %for multiframe DICOM images
  if not(isfield(DATA.Pref,'FasterPreview'))
    DATA.Pref.FasterPreview = true;
  end;
  
  %--- Paths  
  if not(isfield(DATA.Pref,'PacsTempStoragePath'))
    DATA.Pref.PacsTempStoragePath = [pwd filesep 'Database' filesep 'TEMP'];
  end;
  if not(isfield(DATA.Pref,'DicomPath'))
    DATA.Pref.DicomPath = [pwd filesep 'Database' filesep 'DICOM'];
  end;
  if not(isfield(DATA.Pref,'AnalysedPath'))
    DATA.Pref.AnalysedPath = [pwd filesep 'Database' filesep 'Analysed'];
  end;
  if not(isfield(DATA.Pref,'reportsheetpath'))
    DATA.Pref.reportsheetpath = [pwd filesep 'Database' filesep 'Report'];
  end
  if not(isfield(DATA.Pref,'CDPath'))
    DATA.Pref.CDPath = 'D:';
  end;
  if not(isfield(DATA.Pref,'ImageBasePath'))
    DATA.Pref.ImageBasePath = [pwd filesep 'Database'];
  end;
  if not(isfield(DATA.Pref,'NumInterpPoints'))
    DATA.Pref.NumInterpPoints = 15;
  end;
  
  % JT 090929 #411
  if not(isfield(DATA.Pref, 'ShowSeriesDescription'))
    DATA.Pref.ShowSeriesDescription = false;
  end
  
  % JT 091027 #412
  if not(isfield(DATA.Pref, 'HideFilesUnix'))
    DATA.Pref.HideFilesUnix = true;
  end
  
  %JS 100407
  if not(isfield(DATA.Pref, 'SaveVersion'))
    DATA.Pref.SaveVersion = '-v7';
  end
  
  %JS 100511
  if not(isfield(DATA.Pref, 'WebBrowser'))
    DATA.Pref.WebBrowser = 'explorer';
  end
  
  %JS 100511
  if not(isfield(DATA.Pref, 'RetrieveMode'))
    DATA.Pref.RetrieveMode = 2;
  end
  %JS 100511
  if not(isfield(DATA.Pref, 'DebugLog'))
    DATA.Pref.DebugLog = 0;
  end
  %JS 100511
  if not(isfield(DATA.Pref, 'VerboseLog'))
    DATA.Pref.VerboseLog = 0;
  end
  
  %EH 110217
  if not(isfield(DATA.Pref,'RetrieveOptions'))
    DATA.Pref.RetrieveOptions = ''; %'--prefer-little --propose-little --bit-preserving'; %This is used for movescu to configure.
  end;
  if not(isfield(DATA.Pref,'QueryOptions'))
    DATA.Pref.QueryOptions = ''; %'--propose-little'; %This is used for findscu to configure.
  end;
  if not(isfield(DATA.Pref,'ReceiveOptions'))
    DATA.Pref.ReceiveOptions = ''; %'--prefer-little'; %This is used for storescp to configure.
  end;
  if not(isfield(DATA.Pref,'SendOptions'))
    DATA.Pref.SendOptions = ''; %This is used for storescu to configure.
  end;
  
  %NL 110426
  if not(isfield(DATA.Pref,'FastPreviewLoad'))
    DATA.Pref.FastPreviewLoad = true;
  end
  
  %NL 110901
  if not(isfield(DATA.Pref, 'GensvarUrl'))
    DATA.Pref.GensvarUrl = 'http://';
    DATA.Pref.GensvarUsername = '';
    DATA.Pref.GensvarPassword = '';
  end
  
  %NL 111011
  prefrestruct;
end

if not(isfield(DATA.Pref,'UserLogging'))
  DATA.Pref.UserLogging = false;
end;

if not(isfield(DATA.Pref,'UserLogPath'))
  DATA.Pref.UserLogPath = [];
end;

%NL111012
if not(isfield(DATA.Pref.Pacs,'patientconfig'))
  DATA.Pref.Pacs.patientconfig = [1 0 1 1 0 0];
  DATA.Pref.Pacs.studyconfig = [0 0 1 1 1 0 1 1 1 1 0 0 1];
  DATA.Pref.Pacs.seriesconfig = [0 0 1 0 1 1 1 1 1 1 0 0 0];
  DATA.Pref.Pacs.retrievemodel = [1 0 0];
end

%NL111017
if not(isfield(DATA.Pref,'CheckVersion'))
  DATA.Pref.CheckVersion = true;
end

%NL130311
if not(isfield(DATA.Pref.Pacs,'PAFPathname'))
  DATA.Pref.Pacs.PAFPathname = [pwd filesep 'PAF'];
end

%NL130925
if not(isfield(DATA.Pref.Pacs,'SwitchTags'))
  DATA.Pref.Pacs.SwitchTags = 1;
end

%NL140402
if not(isfield(DATA.Pref,'Language'))
  DATA.Pref.Language = 'English';
end

%NL140429
if not(isfield(DATA.Pref,'BackgroundColor'))
  DATA.Pref.BackgroundColor = false;
end

%NL150310
if isfield(DATA.Pref,'ExcludePapilars')
  DATA.Pref = rmfield(DATA.Pref,'ExcludePapilars');
end

%EH160521
if not(isfield(DATA.Pref,'UseProxyServer'))
  DATA.Pref.UseProxyServer = false;
  DATA.Pref.ProxySettings = [];
  DATA.Pref.ProxySettings.HostName = '';
  DATA.Pref.ProxySettings.Port = '';
  DATA.Pref.ProxySettings.UserName = '';
  DATA.Pref.ProxySettings.Password = '';
  useproxyserver(false); %turn off;
else  
  useproxyserver(DATA.Pref.UseProxyServer); %turn on/off;
end;
  
%Obsoleted => disable
DATA.Pref.AllowDicomCache = false; %EH: Set to false since obsoleted

if not(isfield(DATA.Pref,'OpenGLSoftware'))
  DATA.Pref.OpenGLSoftware = false; %Klas: before hardware was default.
end;

if not(isfield(DATA.Pref,'GUIBackgroundColor'))
    DATA.Pref.GUIBackgroundColor = [0.94 0.94 0.94]; 
end;

%--------------------
function prefrestruct
%--------------------
%Restructures DATA.Pref according to PACS refactory project
global DATA
pref = DATA.Pref;

if ~isfield(pref,'Pacs')
  pacs = struct;
  [pacs,pref] = mvfield(pacs,pref,'reportsheetpath','ReportsheetPath');
  [pacs,pref] = mvfield(pacs,pref,'PacsTempStoragePath','TempStoragePath');
  [pacs,pref] = mvfield(pacs,pref,'DicomPath');
  [pacs,pref] = mvfield(pacs,pref,'AnalysedPath');
  [pacs,pref] = mvfield(pacs,pref,'ImageBasePath');
  [pacs,pref] = mvfield(pacs,pref,'DebugLog');
  [pacs,pref] = mvfield(pacs,pref,'VerboseLog');
  [pacs,pref] = mvfield(pacs,pref,'RetrieveMode');
  [pacs,pref] = mvfield(pacs,pref,'RetrieveOptions');
  [pacs,pref] = mvfield(pacs,pref,'QueryOptions');
  [pacs,pref] = mvfield(pacs,pref,'SendOptions');
  pref.Pacs = pacs;
end

if ~isfield(pref,'Dicom')
  dicom = struct;
  [dicom,pref] = mvfield(dicom,pref,'Force16Bit');
  [dicom,pref] = mvfield(dicom,pref,'FasterPreview');
  [dicom,pref] = mvfield(dicom,pref,'NormalizePhase');
  pref.Dicom = dicom;
end

if ~isfield(pref,'Gensvar')
  gensvar = struct;
  [gensvar,pref] = mvfield(gensvar,pref,'GensvarUrl','Url');
  [gensvar,pref] = mvfield(gensvar,pref,'GensvarUsername','Username');
  [gensvar,pref] = mvfield(gensvar,pref,'GensvarPassword','Password');
  pref.Gensvar = gensvar;
end

if ~isfield(pref,'Server')
  server = struct;
  [server,pref] = mvfield(server,pref,'DICOMPort');
  [server,pref] = mvfield(server,pref,'AETitle');
  [server,pref] = mvfield(server,pref,'ReceiveOptions');
  pref.Server = server;
end

DATA.Pref = pref;

%--------------------------------------------------------------------------
function [tostruct, fmstruct] = mvfield(tostruct, fmstruct, fname, newname)
%--------------------------------------------------------------------------
%Moves field fname from fmstruct to tostruct
if nargin < 4
  newname = fname;
end
tostruct.(newname) = fmstruct.(fname);
fmstruct = rmfield(fmstruct,fname);

%-------------------
function savetodisk
%-------------------
global DATA
%This is a copy of save_Callback from segpref.m. In this way we avoid
%to compile in segment when compiling segmentserver.

pathname = getpreferencespath;
Pref = DATA.Pref; %#ok<NASGU> %Saved to file
try
  save([pathname filesep '.segment_preferences.mat'],'Pref');
catch %#ok<CTCH>
  myfailed('Could not save preferences. Write permission? Disk full?',DATA.PrefHandles.fig);
  return;
end;

disp('Preferences saved.');