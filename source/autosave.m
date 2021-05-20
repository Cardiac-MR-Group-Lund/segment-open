function [varargout] = autosave(varargin)
%Tools for autosaving

%Einar Heiberg

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

%macro_helper(varargin{:}); %No need for macro recorder
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init %#ok<DEFNU>
%------------
%Init function check if should start timer

global DATA

if DATA.Pref.AutoSave
  startit;
else  
  DATA.AutoSaveTimer = [];
end

%---------------
function startit
%---------------
%Starts the timer

global DATA

if ~isfield(DATA,'AutoSaveTimer')
  DATA.AutoSaveTimer = [];
end

%This code is inspired from sectra.m
if ~isempty(DATA.AutoSaveTimer) && isa(DATA.AutoSaveTimer,'timer')
  %it is a timer, ensure it is stopped properly
  stop(DATA.AutoSaveTimer);
  delete(DATA.AutoSaveTimer);
end

t = 15*60; %15 minutes in seconds
DATA.AutoSaveTimer = timer('StartDelay',t,'Period',t,'ExecutionMode','fixedSpacing');
DATA.AutoSaveTimer.TimerFcn = 'autosave(''doautosave'')';
start(DATA.AutoSaveTimer);
disp('Autosave timer started.');

%--------------
function stopit
%--------------
%Stop autosave timer

global DATA

%Try to cleanup
try
  cleanup;
catch me
  disp('Could not perform cleanup before stopping autosave');
  mydispexception(me);
end

%Try to stop
try
  stop(DATA.AutoSaveTimer);
  delete(DATA.AutoSaveTimer);
  DATA.AutoSaveTimer = [];
  disp('Autosave time deleted.');
catch me
  disp('Something went wrong stopping autosave timer.');
  mydispexception(me);  
end

%------------------
function doautosave %#ok<DEFNU>
%------------------
%Performs the autosave, first performs cleanup

global DATA

%Do not save if no data loaded.
if ~DATA.DataLoaded
  return;
end

%Try to cleanup prior to save
try
  cleanup;
catch me
  disp('Could not perform cleanup.');
  mydispexception(me);
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
pathname = getpreferencespath;
filename = sprintf('autosave-%s.mat',datestr(now,'yyyymmddTHHMMSS',now));
h = msgbox(dprintf('Autosaving %s%s%s',pathname,filesep,filename)); %Reason for msgbox instead of mymsgbox is that msgbox is non modular (i.e code does not stop and wait for ok).
disp(sprintf('Autosaving %s%s%s',pathname,filesep,filename)); %#ok<DSPS>
filemenu('saveallas_helper',pathname,filename,topatientdatabase,silent)
disp('Autosave done.');

try
  delete(h); %Try to delete the message box
catch
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
  disp('No autosave files to delete.');
  return;
end

disp('Checking for autosave files to clean.');

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
      disp(sprintf('Deleting autosave file %s',f(loop).name)); %#ok<DSPS>
      delete([pathname filesep f(loop).name]);
    end
  end
else
  disp('No files to delete.');
end


