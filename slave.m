function slave(name)
%SLAVE Slave to do work from the myserver class.
%
%Should be compiled to standalone file.
%
%See also MYCLIENTSERVER.

%Todo later implement optional argument for log-file

%Communication protocol is described in myclientserver

global DATA 

checkpath();

DATA.RecordMacro = false; %may be used by subfunctions

%Say hello
disp(sprintf('This is slave %s, helper program to Segment.',name));

if nargin<1
  error('Not enough input arguments, expected one input arguments.');
end;

%Set up struct to define server.
status = [];
status.Name = removeforbiddenchars(name);
status.Verbose = true;
status.Path = getpreferencespath;
status.Pause = 0.1; %s between checks
status.LastMessageID = 0;

disp('Starting...');
disp('');
disp('********************************************');
disp('**** You should not close this window *****');
disp('********************************************');

%Load preferences from Segment preferences. Used for thumbnail creation
%etc.
loadpreferences;

%create output file
outfid = fopen(getoutboxname(status),'w');
if outfid<1
  %Could not open file for writing
  pause(1);
  outfid = fopen([status.Path filesep '.' name '.outbox'],'w');
  if outfid<1
    error('Could not open outbox at second attempt.');
  end;
end;

%Close the file again, it is now empty
fclose(outfid);

arequitting = false;

while ~arequitting;
  
  [stri,status] = read(status); %Read from the inbox
  
  if ~isempty(stri)
    if status.Verbose
      disp(sprintf('Received string: %s',stri));
    end;
    arequitting = receive(status,stri);
  else 
    pause(status.Pause);
  end;

end;

disp('****************************************');
disp('**** You may now close this window *****');
disp('****************************************');

%------------------------------------
function name = getoutboxname(status)
%------------------------------------
name = [status.Path filesep '.' status.Name '.outbox'];

%------------------------------------
function name = getinboxname(status)
%------------------------------------
name = [status.Path filesep '.' status.Name '.inbox'];

%----------------------------------------
function [fullstri,status] = read(status)
%----------------------------------------
%Reads message from inbox. If no inbox or message exist,

infid = fopen(getinboxname(status),'r');
if infid<1
  fullstri = '';
else
  %Try to read from the file
  fullstri = '';
  stri = '';
  try
    while ~feof(infid) && ~isnumeric(stri);
      stri = fgetl(infid);
      if ~isnumeric(stri) && ~isempty(stri)
        fullstri = [fullstri stri sprintf('\n')]; %#ok<AGROW>
      end;
    end;
    fclose(infid);
    
    %Decode the message to get the id number
    if length(fullstri)>1
      try
        [~,~,id] = decodemessage(fullstri);
      catch me
        mydispexception(me);
        return;
      end;
    
      if id<=status.LastMessageID
        %Already read this message, ignore it.
        fullstri = '';
      else
        status.LastMessageID = id;
      end;
    end;
    
    %--Remove message
    %Deleted this section since it is dangerous
    %infid = fopen(getinboxname(status),'w');
    %fclose(infid); 
    
    %avoid []
    if ~ischar(fullstri)
      fullstri = '';
    end;
    
  catch %#ok<CTCH>
    fullstri = '';
  end;
end;

%---------------------------------------
function [are_quitting] = receive(status,stri)
%---------------------------------------
%Called when a string is received.
global DATA  %#ok<NUSED>

% per default we are not quitting
are_quitting = false;

if isempty(stri)
  disp('Empty string. Timeout?');
  return;
end;

try
  [~,~,id,cmd,args] = decodemessage(stri);
catch me
  mydispexception(me);
  return;
end;

switch lower(cmd)
  case 'ping'
    send(status,0,id,'ping','pong');
  case 'kill'
    send(status,0,id,'kill','dying');
    pause(0.5);
    are_quitting = true;
  case 'cd'
    try
      cd(args);
      send(status,0,id,'cd','successful');
    catch me
      mydispexception(me);
      send(status,-1,id,'cd',args);
    end;
  case 'cache'
    if status.Verbose
      disp(sprintf('Caching DICOM info is obsoleted %s',args));
    end;
    send(status,-1,id,'cache','failed');    
  case 'thumbnails'
    if status.Verbose    
      disp(sprintf('Caching thumbnails in %s',args));
    end;
    try
      thumbnails('createthumbnails_helper',args,true) %creates in given folder, true is silent
      send(status,0,id,'thumbnail','successful');
    catch me
      mydispexception(me);      
      send(status,-1,id,'thumbnail','failed');
    end;
  case 'pause'
    %Implemented for testing.
    try
      args = str2double(args);
      if status.Verbose
        disp(sprintf('Pause for %d seconds',args));
      end;
      try
        pause(args);
        send(status,0,id,'pause','successful');
      catch me
        mydispexception(me);
        send(status,-1,id,'pause','successful');
      end;
    catch me
      disp(sprintf('Could not convert %s to number.',args));
      mydispexception(me);
      send(status,-1,id,'pause','failed, could not convert argument to seconds.');
    end;
  case 'system'
    [stat,w] = system(args);
    if (stat~=0) %||any(findstr(lower(w),'failed'))
      send(status,-1,id,'system',validstring(w));
    else
      send(status,0,id,'system',validstring(w));
    end;
  otherwise
    disp(sprintf('Unknown command %s',cmd));
    send(status,-1,id,'unknown','unknown command');
end;

%----------------------------------
function send(status,s,id,cmd,args)
%----------------------------------
%Send something back to slave.

if nargin<1
  error('Expected at least one input arguments.');
end;

if nargin<2
  s = 0;
end;

if nargin<3
  cmd = '';
end;

if nargin<4
  args = '';
end;

stri = encodemessage(status.Name,s,id,cmd,validstring(args));

if status.Verbose
  disp(sprintf('Sending string: %s',stri));
end;

%Read to ensure it is not already containing stuff.
outfid = fopen(getoutboxname(status),'r');
if outfid>0
  readstri = fgetl(outfid);
  fclose(outfid);
  if ~isempty(readstri) && ischar(readstri)
    disp(sprintf('Data already present %s',readstri));
  end;
end;

outfid = fopen(getoutboxname(status),'w');
if outfid<1
  pause(status.Pause);
else
  %Okey to open => write and close
  fprintf(outfid,'%s\n',stri);
  fclose(outfid);
  return; %to avoid second attempt
end;

%Second attempt
outfid = fopen(getoutboxname(status),'w');
if outfid<1
  error('Could not send message back to client server');
else
  %Okey to open => write and close
  fprintf(outfid,'%s',stri);
  fclose(outfid);
  return;
end;

%--------------------------------
function stri = validstring(stri)
%--------------------------------
%Ensures that result is a valid string without carrage returns or |

stri(stri=='|') = ':';
%stri(stri==sprintf('\n')) = ' ';