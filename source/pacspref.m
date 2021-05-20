function varargout = pacspref(varargin)
%PACSPREF Helper function for SEGMENT, GUI for PACS preferences

%Nils Lundahl

if nargin == 0  % LAUNCH GUI
  
  init;
  
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

%------------
function init
%------------
%Initiate Pacs Preferences GUI
global DATA

fig = mygui('pacspref.fig');

% Generate a structure of handles to pass to callbacks
handles = fig.handles;
DATA.PrefHandlesPacs = handles;
DATA.PrefHandlesPacs.fig = handles.figure1;

patientkids = get(handles.patientlevelpanel,'Children');
studykids = get(handles.studylevelpanel,'Children');
serieskids = get(handles.serieslevelpanel,'Children');
set(patientkids(2:end-1),'Callback','pacspref(''queryconfig_Callback'',''patient'')');
set(studykids(2:end-1),'Callback','pacspref(''queryconfig_Callback'',''study'')');
set(serieskids(2:end-1),'Callback','pacspref(''queryconfig_Callback'',''series'')');
set(handles.patientradio,'SelectionChangeFcn','pacspref(''queryconfig_Callback'',''patient'')')
set(handles.studyradio,'SelectionChangeFcn','pacspref(''queryconfig_Callback'',''study'')');
set(handles.seriesradio,'SelectionChangeFcn','pacspref(''queryconfig_Callback'',''series'')');
set(handles.retrieveradio,'SelectionChangeFcn','pacspref(''retrieveradio_Callback'')');

update;

%--------------
function update
%--------------
%Update handles
global DATA

%General PACS and Database Settings
set(DATA.PrefHandlesPacs.verboselogcheckbox,'Value',DATA.Pref.Pacs.VerboseLog);
set(DATA.PrefHandlesPacs.debuglogcheckbox,'Value',DATA.Pref.Pacs.DebugLog);

%PACS Query Configurations
%Patient Level
cvec = DATA.Pref.Pacs.patientconfig;
set(DATA.PrefHandlesPacs.piprradiobutton,'value',cvec(1));
set(DATA.PrefHandlesPacs.pipsrradiobutton,'value', cvec(2));
set(DATA.PrefHandlesPacs.ppnmcheckbox,'value', cvec(3));
set(DATA.PrefHandlesPacs.ppidcheckbox,'value', cvec(4));
set(DATA.PrefHandlesPacs.ppbdcheckbox,'value', cvec(5));
set(DATA.PrefHandlesPacs.ppsexcheckbox,'value', cvec(6));

%Study Level
cvec = DATA.Pref.Pacs.studyconfig;
set(DATA.PrefHandlesPacs.stiprradiobutton,'value', cvec(1));
set(DATA.PrefHandlesPacs.stisrradiobutton,'value', cvec(2));
set(DATA.PrefHandlesPacs.stipsrradiobutton,'value', cvec(3));
set(DATA.PrefHandlesPacs.stpnmcheckbox,'value', cvec(4));
set(DATA.PrefHandlesPacs.stpidcheckbox,'value', cvec(5));
set(DATA.PrefHandlesPacs.stpsexcheckbox,'value', cvec(6));
set(DATA.PrefHandlesPacs.ststidcheckbox,'value', cvec(7));
set(DATA.PrefHandlesPacs.ststdtcheckbox,'value', cvec(8));
set(DATA.PrefHandlesPacs.ststtcheckbox,'value', cvec(9));
set(DATA.PrefHandlesPacs.ststdscheckbox,'value', cvec(10));
set(DATA.PrefHandlesPacs.stancheckbox,'value', cvec(11));
set(DATA.PrefHandlesPacs.stmodscheckbox,'value', cvec(12));
set(DATA.PrefHandlesPacs.stmodcheckbox,'value', cvec(13));

%Series Level
cvec = DATA.Pref.Pacs.seriesconfig;
set(DATA.PrefHandlesPacs.seiprradiobutton,'value', cvec(1));
set(DATA.PrefHandlesPacs.seisrradiobutton,'value', cvec(2));
set(DATA.PrefHandlesPacs.seipsrradiobutton,'value', cvec(3));
set(DATA.PrefHandlesPacs.sepnmcheckbox,'value', cvec(4));
set(DATA.PrefHandlesPacs.sepidcheckbox,'value', cvec(5));
set(DATA.PrefHandlesPacs.sestdtcheckbox,'value', cvec(6));
set(DATA.PrefHandlesPacs.sestidcheckbox,'value', cvec(7));
set(DATA.PrefHandlesPacs.seseidcheckbox,'value', cvec(8));
set(DATA.PrefHandlesPacs.sesetcheckbox,'value', cvec(9));
set(DATA.PrefHandlesPacs.sesedscheckbox,'value', cvec(10));
set(DATA.PrefHandlesPacs.seancheckbox,'value', cvec(11));
set(DATA.PrefHandlesPacs.semodscheckbox,'value', cvec(12));
set(DATA.PrefHandlesPacs.semodcheckbox,'value', cvec(13));

set(DATA.PrefHandlesPacs.queryoptionsedit,'String',DATA.Pref.Pacs.QueryOptions);

%Retrieve Configurations
set(DATA.PrefHandlesPacs.riprradiobutton,'Value',DATA.Pref.Pacs.retrievemodel(1));
set(DATA.PrefHandlesPacs.risrradiobutton,'Value',DATA.Pref.Pacs.retrievemodel(2));
set(DATA.PrefHandlesPacs.ripsrradiobutton,'Value',DATA.Pref.Pacs.retrievemodel(3));

set(DATA.PrefHandlesPacs.retrievemodepopupmenu,'Value',DATA.Pref.Pacs.RetrieveMode);
set(DATA.PrefHandlesPacs.retrieveoptionsedit,'String',DATA.Pref.Pacs.RetrieveOptions);

% PACS connections
updatepacscon;

%---------------------
function updatepacscon
%---------------------
% Updates the connection listbox

global DATA

% Update the list of possible connections.
f = dir('*.con');
li = cell(1,length(f));
for loop=1:length(f)
  load(f(loop).name,'res','-mat');
  li{loop} = res.DescriptiveName;
end

% update v
v = mygetlistbox(DATA.PrefHandlesPacs.connecttolistbox);
v = min([length(li) v]); % v shouldn't be bigger then li
if (v == 0) && (~isempty(li)) % v shouldn't be zero if li not empty
  v = 1;
end

% update listbox, first set value 1 to avoid warnings
set(DATA.PrefHandlesPacs.connecttolistbox, 'value', 0);
set(DATA.PrefHandlesPacs.connecttolistbox, 'String', li);
set(DATA.PrefHandlesPacs.connecttolistbox, 'value', v);

%update labels if in GUI that wants you to.
DATA.updatepacslabels;

%-----------------------------
function verboselog_Callback 
%----------------------------
global DATA

DATA.Pref.Pacs.VerboseLog=get(DATA.PrefHandlesPacs.verboselogcheckbox,'value');

%-----------------------------
function debuglog_Callback 
%----------------------------
global DATA

DATA.Pref.Pacs.DebugLog=get(DATA.PrefHandlesPacs.debuglogcheckbox,'value');

%-----------------------------------
function queryconfig_Callback(level)
%-----------------------------------
%Callback from query configuration radiobuttons and checkboxes
global DATA

switch level
  case 'patient'
    cvec = zeros(1,6);
    cvec(1) = get(DATA.PrefHandlesPacs.piprradiobutton,'value');
    cvec(2) = get(DATA.PrefHandlesPacs.pipsrradiobutton,'value');
    cvec(3) = get(DATA.PrefHandlesPacs.ppnmcheckbox,'value');
    cvec(4) = get(DATA.PrefHandlesPacs.ppidcheckbox,'value');
    cvec(5) = get(DATA.PrefHandlesPacs.ppbdcheckbox,'value');
    cvec(6) = get(DATA.PrefHandlesPacs.ppsexcheckbox,'value');
    DATA.Pref.Pacs.patientconfig = cvec;
  case 'study'
    cvec = zeros(1,13);
    cvec(1) = get(DATA.PrefHandlesPacs.stiprradiobutton,'value');
    cvec(2) = get(DATA.PrefHandlesPacs.stisrradiobutton,'value');
    cvec(3) = get(DATA.PrefHandlesPacs.stipsrradiobutton,'value');
    cvec(4) = get(DATA.PrefHandlesPacs.stpnmcheckbox,'value');
    cvec(5) = get(DATA.PrefHandlesPacs.stpidcheckbox,'value');
    cvec(6) = get(DATA.PrefHandlesPacs.stpsexcheckbox,'value');
    cvec(7) = get(DATA.PrefHandlesPacs.ststidcheckbox,'value');
    cvec(8) = get(DATA.PrefHandlesPacs.ststdtcheckbox,'value');
    cvec(9) = get(DATA.PrefHandlesPacs.ststtcheckbox,'value');
    cvec(10) = get(DATA.PrefHandlesPacs.ststdscheckbox,'value');
    cvec(11) = get(DATA.PrefHandlesPacs.stancheckbox,'value');
    cvec(12) = get(DATA.PrefHandlesPacs.stmodscheckbox,'value');
    cvec(13) = get(DATA.PrefHandlesPacs.stmodcheckbox,'value');
    DATA.Pref.Pacs.studyconfig = cvec;
  case 'series'
    cvec = zeros(1,13);
    cvec(1) = get(DATA.PrefHandlesPacs.seiprradiobutton,'value');
    cvec(2) = get(DATA.PrefHandlesPacs.seisrradiobutton,'value');
    cvec(3) = get(DATA.PrefHandlesPacs.seipsrradiobutton,'value');
    cvec(4) = get(DATA.PrefHandlesPacs.sepnmcheckbox,'value');
    cvec(5) = get(DATA.PrefHandlesPacs.sepidcheckbox,'value');
    cvec(6) = get(DATA.PrefHandlesPacs.sestdtcheckbox,'value');
    cvec(7) = get(DATA.PrefHandlesPacs.sestidcheckbox,'value');
    cvec(8) = get(DATA.PrefHandlesPacs.seseidcheckbox,'value');
    cvec(9) = get(DATA.PrefHandlesPacs.sesetcheckbox,'value');
    cvec(10) = get(DATA.PrefHandlesPacs.sesedscheckbox,'value');
    cvec(11) = get(DATA.PrefHandlesPacs.seancheckbox,'value');
    cvec(12) = get(DATA.PrefHandlesPacs.semodscheckbox,'value');
    cvec(13) = get(DATA.PrefHandlesPacs.semodcheckbox,'value');
    DATA.Pref.Pacs.seriesconfig = cvec;
end

%-----------------------------
function retrievemode_Callback % %called from segprefadvanced.fig
%-----------------------------
global DATA

DATA.Pref.Pacs.RetrieveMode=get(DATA.PrefHandlesPacs.retrievemodepopupmenu,'value');

%------------------------------------
function retrieveoptionsedit_Callback 
%------------------------------------
global DATA

stri = mygetedit(DATA.PrefHandlesPacs.retrieveoptionsedit);
DATA.Pref.Pacs.RetrieveOptions = stri;

%------------------------------
function retrieveradio_Callback 
%------------------------------
global DATA

cvec = zeros(1,3);
cvec(1) = get(DATA.PrefHandlesPacs.riprradiobutton,'value');
cvec(2) = get(DATA.PrefHandlesPacs.risrradiobutton,'value');
cvec(3) = get(DATA.PrefHandlesPacs.ripsrradiobutton,'value');
DATA.Pref.Pacs.retrievemodel = cvec;

%---------------------------------
function queryoptionsedit_Callback 
%---------------------------------
global DATA

stri = mygetedit(DATA.PrefHandlesPacs.queryoptionsedit);
DATA.Pref.Pacs.QueryOptions = stri;

%---------------------------
function addpacscon_Callback %#ok<*DEFNU>
%---------------------------
% Executes when user press add in advanced dicom pref.
% Ask user for input, stores con file and call updatepacscon.

global DATA

% Ask user for input
s = [];
s(1).Field = 'Descriptive_Name';
s(1).Label = dprintf('Descriptive Name');
s(1).Default = 'Name               ';
n = 2;
s(n).Field = 'PACS_AE_Title';
s(n).Label = dprintf('PACS AE Title');
s(n).Default = 'HOSPITAL_PACS'; %Called_AE
n = n+1;
s(n).Field = 'PACS_IP_Number';
s(n).Label = dprintf('PACS IP Number');
s(n).Default = '127.0.0.1'; %Peer_IP
n = n+1;
s(n).Field = 'PACS_Port_Number';
s(n).Label = dprintf('PACS Port Number');
s(n).Default = '4006'; %Peer_Port
n = n+1;
s(n).Field = 'Segment_Port_Number';
s(n).Label = dprintf('Segment Port Number');
s(n).Default = '104'; %Port
n = n+1;
s(n).Field = 'Retrieve_AE_Title';
s(n).Label = dprintf('Retrieve AE Title');
s(n).Default = DATA.Pref.Server.AETitle; %Take from Segment Server AE Title
n = n+1;
s(n).Field = 'Query_AE_Title';
s(n).Label = dprintf('Query AE Title');
s(n).Default = DATA.Pref.Server.AETitle; %Take from Segment Server AE Title
n = n+1;
s(n).Field = 'StoreSCU_AE_Title';
s(n).Label = dprintf('StoreSCU AE Title');
s(n).Default = DATA.Pref.Server.AETitle; %Take from Segment Server AE Title
[res,ok] = myinputstruct(s,'',20); %20 = width of editbox
if not(ok)
  return;
end

temp = res; 
s = [];
s.DescriptiveName = temp.Descriptive_Name;
s.Called_AE = temp.PACS_AE_Title;
s.Peer_IP = temp.PACS_IP_Number;
s.Peer_Port = temp.PACS_Port_Number;
s.Port = temp.Segment_Port_Number;
s.Retrieve_AE_Title = temp.Retrieve_AE_Title;
s.Query_AE_Title = temp.Query_AE_Title;
s.StoreSCU_AE_Title = temp.StoreSCU_AE_Title;
clear temp;

% Store the con file
% [filename,pathname] = myuiputfile('*.con','Save Connection As: (place in Segment folder)');
% if isequal(filename,0)||isequal(pathname,0)
%   return;
% end

%The .con file is saved in main software folder as DescriptiveName.con
filename=[s.DescriptiveName '.con'];
pathname=DATA.SegmentFolder;
res = s;
res.Calling_AE = 'DEFAULT_CALLING_AE'; % Added for legacy reasons
save(fullfile(pathname,filename),'res', DATA.Pref.SaveVersion);

% calls updatepacscon
updatepacscon;

%----------------------------
function editpacscon_Callback
%----------------------------
% Executed when user press edit in advanced pref.
% Displays and let user edit content of a .con file.

global DATA

% Load the file
v = mygetlistbox(DATA.PrefHandlesPacs.connecttolistbox);
f = dir('*.con');
try
  load(f(v).name, 'res', '-mat');
catch e
  myfailed(dprintf('Could not open the selected connection. Error is ''%s''',e.message),DATA.GUI.PacsCon);
  updatepacscon;
  return
end

%Remove non used field from display
res = rmfield(res, 'Calling_AE'); 

%Prepare variable with new fieldnames
s = [];
n = 1;
s(n).Field = 'Descriptive_Name';
s(n).Label = dprintf('Descriptive Name');
s(n).Default = res.DescriptiveName;
n = n+1;
s(n).Field = 'PACS_AE_Title';
s(n).Label = dprintf('PACS AE Title');
s(n).Default = res.Called_AE;
n = n+1;
s(n).Field = 'PACS_IP_Number';
s(n).Label = dprintf('PACS IP Number');
s(n).Default = res.Peer_IP;
n = n+1;
s(n).Field = 'PACS_Port_Number';
s(n).Label = dprintf('PACS Port Number');
s(n).Default = res.Peer_Port;
n = n+1;
s(n).Field = 'Segment_Port_Number';
s(n).Label = dprintf('Segment Port Number');
s(n).Default = res.Port;
n = n+1;
if isfield(res,'Retrieve_AE_Title')
  s(n).Field = 'Retrieve_AE_Title';
  s(n).Label = dprintf('Retrieve AE Title');
  s(n).Default = res.Retrieve_AE_Title;
  n = n+1;
else  
  s(n).Field = 'Retrieve_AE_Title';
  s(n).Label = dprintf('Retrieve AE Title');  
  s(n).Default = DATA.Pref.Server.AETitle; %Take from Segment Server AETitle
  n = n+1;
end
if isfield(res,'Query_AE_Title')
  s(n).Field = 'Query_AE_Title';
  s(n).Label = dprintf('Query AE Title');
  s(n).Default = res.Query_AE_Title;
  n = n+1;
else
  s(n).Field = 'Query_AE_Title';
  s(n).Label = dprintf('Query AE Title');
  s(n).Default = DATA.Pref.Server.AETitle; %Take from Segment Server AETitle
  n = n+1;
end
if isfield(res,'StoreSCU_AE_Title')
  s(n).Field = 'StoreSCU_AE_Title';
  s(n).Label = dprintf('StoreSCU AE Title');
  s(n).Default = res.StoreSCU_AE_Title;
  n = n+1; %#ok<NASGU>
else  
  s(n).Field = 'StoreSCU_AE_Title';
  s(n).Label = dprintf('StoreSCU AE Title');
  s(n).Default = DATA.Pref.Server.AETitle; %Take from Segment Server AETitle
  n = n+1; %#ok<NASGU>
end

% Let user edit the file and save it
[s,ok] = myinputstruct(s,'',20); %20 is width of editboxes

if not(ok)
  return;
end

%Convert to internal fieldname representation
temp = s; 
res = [];
res.DescriptiveName = temp.Descriptive_Name;
res.Called_AE = temp.PACS_AE_Title;
res.Peer_IP = temp.PACS_IP_Number;
res.Peer_Port = temp.PACS_Port_Number;
res.Port = temp.Segment_Port_Number;
res.Retrieve_AE_Title = temp.Retrieve_AE_Title;
res.Query_AE_Title = temp.Query_AE_Title;
res.StoreSCU_AE_Title = temp.StoreSCU_AE_Title;
clear temp;

res.Calling_AE = 'DEFAULT_CALLING_AE';


try
  save(f(v).name, 'res', DATA.Pref.SaveVersion);
catch
  myfailed('Could not save connection file. Disk write protected?');
end

%Update the list
updatepacscon;

%------------------------------
function deletepacscon_Callback
%------------------------------
% Executes when user press delete in advanced pref.
% Deletes selected pacscon.

global DATA

% Delete connection file
v = mygetlistbox(DATA.PrefHandlesPacs.connecttolistbox);
f = dir('*.con');
if v>length(f)
  errormsg='Number of files and number of connections does not match.';
  myfailed(errormsg,DATA.GUI.PacsCon);
  updatepacscon;
  return;
end
ok = mydel(f(v).name);
if not(ok)
  errormsg=sprintf(...
    'Could not remove the connection file %s.\n\n\Error message:\n\n%s\n',...
    f(v).name, w);
  myfailed(errormsg,DATA.GUI.PacsCon);
  return;
end

% Update listbox value and updates pacscon
if length(f)>1
  set(DATA.PrefHandlesPacs.connecttolistbox, 'value', 1);
else
  set(DATA.PrefHandlesPacs.connecttolistbox, 'value', 0);
end
updatepacscon;

%------------------------
function default_Callback
%------------------------
default;
update;

%---------------
function default
%---------------
global DATA

DATA.Pref.Pacs.VerboseLog=0;
DATA.Pref.Pacs.DebugLog=0;

DATA.Pref.Pacs.patientconfig = [0 1 1 1 0 0];
DATA.Pref.Pacs.studyconfig = [0 1 0 1 1 0 1 1 1 1 1 1 0];
DATA.Pref.Pacs.seriesconfig = [0 1 0 1 1 1 1 1 1 1 1 1 0];
DATA.Pref.Pacs.QueryOptions = ''; %'--propose-little';

DATA.Pref.Pacs.retrievemodel = [0 1 0];
DATA.Pref.Pacs.RetrieveMode=4;
DATA.Pref.Pacs.RetrieveOptions = ''; %'--prefer-little --propose-little --bit-preserving';

%---------------------
function undo_Callback
%---------------------
loadpreferences;
update;

%----------------------
function close_Callback
%----------------------
global DATA

closereq;
DATA.PrefHandlesPacs = [];
