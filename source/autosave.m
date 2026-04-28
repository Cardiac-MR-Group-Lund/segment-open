function [varargout] = autosave(varargin)
%Tools for autosaving
%
%Autosave was included many years ago in Segment but is now re-introduced
%in 2020-10-17. Previous implementation was problematic due to overwriting,
%and also wrote to working folder that could have been a slow network disc.
%New implementation writes series of files to getpreferencespath.
%
%Autosave clears files such as:
% 1) there no more than 500 MB in old files (=> max space 1 GB)
% 2) older files than 24 hours.
% 3) no more than 10 files stored.
%
%Potential problems with new implementation
%1) Potential to get non anonymized data in getpreferencespath
%2) If the files are very big then a lot of data will be in the getpreferencespath (max 1GB)
%3) There is no guarantee that the files are save in a bad state du to
%ongoing computations, some precautions are performed.

%Einar Heiberg

%#ok<*GVMIS> 

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%--------------
function init 
%--------------
%called from openfile.m

global DATA 

timerobj = timerfind('Name','AutoSaveTimer');
if isempty(timerobj) 
  if DATA.Pref.AutoSave %check preferences 
    startit; 
  end
else
  % timer exist
  if DATA.Pref.AutoSave
    % restart timer if AutoSave is on
    startit(timerobj)
  else
    % stop timer gracefully because it shouldn't be on based on preferences
    removeit(timerobj);
  end    
end

%---------------
function startit(timerobj)
%---------------
%Starts the timer
if nargin == 0 || ~isa(timerobj,'Timer')
  % no timer object was provide or the provided object is not a timer
  timerobj = timerfind('Name','AutoSaveTimer');
end

%This code is inspired from sectra.m
istimeron = false;
if ~isempty(timerobj)
  %it is a timer, ensure it is stopped properly
  stop(timerobj);
  delete(timerobj);
  istimeron = true;
end
clear timerobj
t = 15*60; %15 minutes in seconds
timerobj = timer('StartDelay',t,...
                 'Period',t,...
                 'ExecutionMode','fixedSpacing',...
                 'TimerFcn','autosave(''doautosave'')',...
                 'Name','AutoSaveTimer');
start(timerobj);
if istimeron
  logdisp(sprintf('Autosave timer restarted with %0.5g minutes interval.',t/60));
else
  logdisp(sprintf('Autosave timer started with %0.5g minutes interval.',t/60));
end

%--------------
function stopit 
%--------------
%Stop autosave timer, called from maingui.m

timerobj = timerfind('Name','AutoSaveTimer');
%Try to cleanup
try
  cleanup;
catch me
  logdisp('Could not perform cleanup before stopping autosave');
  mydispexception(me);
end

if ~isempty(timerobj)
  try
    stop(timerobj);
  catch me
    mydispexception(me)
    logdisp('Could not stop autosave timer')
  end
end

%--------------------------
function removeit(timerobj)
%--------------------------
%Stop autosave timer

%Try to cleanup
try
  cleanup;
catch me
  logdisp('Could not perform cleanup before stopping autosave');
  mydispexception(me);
end
if nargin == 0 || ~isa(timerobj,'Timer')
  % no timer object was provide or the provided object is not a timer
  timerobj = timerfind('Name','AutoSaveTimer');
end

%Try to stop
if ~isempty(timerobj)
  try
    stop(timerobj);
    delete(timerobj);
    logdisp('Autosave timer deleted.');
  catch me
    logdisp('Something went wrong stopping autosave timer.');
    mydispexception(me);  
  end
end

%------------------
function doautosave
%------------------
%Performs the autosave, first performs cleanup

global DATA SET

%Try to cleanup prior to save
try
  cleanup;
catch me
  logdisp('Could not perform autosave cleanup.');
  mydispexception(me);
end

%Do not save if no data loaded.
if ~DATA.DataLoaded
  return;
end

%Do not need to save if NeedToSave is not updated
if ~DATA.NeedToSave
  return;
end

%check if there is 'myworkon' running
pointer = get(DATA.imagefig,'pointer');

while isequal(pointer,'watch')
  pause(5); %pause 5 seconds to wait for ongoing computations to finish, user interface is not locked while we perform pause
  pointer = get(DATA.imagefig,'pointer');
end

%Now the current pointer is no longer a waitbar, this indicates to use that
%it is relatively safe to do autosave

%Save it
topatientdatabase = false;
silent = true;
fast = true;
pathname = getpreferencespath;
oldfilename = SET(1).FileName;

filename = sprintf('autosave-%s.mat',datestr(now,'yyyymmddTHHMMSS',now));
h = msgbox(dprintf('Autosaving %s%s%s',pathname,filesep,filename)); %Reason for msgbox instead of mymsgbox is that msgbox is non modular (i.e code does not stop and wait for ok).
logdisp(sprintf('Autosaving %s%s%s',pathname,filesep,filename));
try
  filemenu('saveallas_helper',pathname,filename,topatientdatabase,silent,fast);
  setfilenamehelper(oldfilename);
catch
  setfilenamehelper(oldfilename);
end
logdisp('Autosave done.');

try
  delete(h); %Try to delete the message box
catch
end

%-----------------------------------
function setfilenamehelper(filename)
%-----------------------------------
%Sets filename for all stacks.

global SET

for loop = 1:length(SET)
  SET(loop).FileName = filename;
end

%---------------
function cleanup
%---------------
%Removes old autosaves

%Removes:
%1) autosaves older than 24 hours
%2) only keeps 4 files (in practise there will be 5 files as cleanup is
%called prior to doautosave.

pathname = getpreferencespath;

%Get potential autosave files to delete
f = dir([pathname filesep 'autosave*.mat']);

if isempty(f)
  logdisp('No autosave files to delete.',true);
  return;
end

logdisp('Checking for autosave files to clean.');

%Files are returned alphanumerically, we should delete "the first ones"
numfiles = length(f);

%find cumulative size of files, counting backwards
cumbytes = zeros(numfiles,1);
totbytes = 0;
for loop = numfiles:-1:1
  totbytes = totbytes + f(loop).bytes;
  cumbytes(loop) = totbytes;
end

%Delete so that there is no more than 500 MB stored
bytelimit = 500*1024*1024;
dodelete = cumbytes>bytelimit;

%Also delete old files
fileage = now-cat(1,f(:).datenum); %unit is days
dodelete = dodelete | (fileage>1) ; %files older than 24 hours

if sum(~dodelete)>10
  ind = find(~dodelete);
  dodelete(ind(1)) = true; %delete the first one
end

%loop over to delete
if sum(dodelete)>0
  for loop = 1:length(dodelete)
    if dodelete(loop)
      logdisp(sprintf('Deleting autosave file %s',f(loop).name));
      delete([pathname filesep f(loop).name]);
    end
  end
else
  logdisp('No autosave files to delete.');
end


