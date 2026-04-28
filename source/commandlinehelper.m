function [varargout] = commandlinehelper(arg,varargin)
%Helper function to handle commandline options

%Einar Heiberg

persistent s %store options in persistent variable

%--- Define commandline options
c = [];
n = 1;
c(n).switch = 'StudiesFolder'; c(n).nargin = 1; c(n).default = ''; c(n).description = '<dir>: folder software reads input files (dicoms), each subfolder is one series';
n = n + 1;
c(n).switch = 'DICOMStorageFolder'; c(n).nargin = 1; c(n).default = ''; c(n).description = '<dir>: folder where output dicom files are stored';
n = n + 1;
c(n).switch = 'LogFolder'; c(n).nargin = 1; c(n).default = ''; c(n).description = '<dir>: folder where log-files are stored';
n = n + 1;
c(n).switch = 'NoQuit'; c(n).nargin = 0; c(n).default = false; c(n).description = ': Disable quit options';
n = n + 1;
c(n).switch = 'NoPopUp'; c(n).nargin = 0; c(n).default = false; c(n).description = ': Disable pop-up when starting';
n = n + 1;
c(n).switch = 'License'; c(n).nargin = 1; c(n).default = ''; c(n).description = '<file>: License file';
n = n + 1;
c(n).switch = 'TempStorageFolder'; c(n).nargin = 1; c(n).default = ''; c(n).description = '<dir>: folder for temporal file storage';
n = n + 1;
c(n).switch = 'UseDefaultStudiesFolder'; c(n).nargin = 0; c(n).default = false; c(n).description = '<dir>: Create studies folder if it not existing according to default preference settings';
n = n + 1;
c(n).switch = 'UseComputerAETitle'; c(n).nargin = 0; c(n).default = false; c(n).description = '<dir>: Use computer AETitle in case if in segpref was different one assigned';
n = n + 1;
c(n).switch = 'NoGUI'; c(n).nargin = 0; c(n).default = false; c(n).description = ': Starts and quits the software in the background. No GUI available.';
n = n + 1;
c(n).switch = 'AutoMate'; c(n).nargin = 0; c(n).default = false; c(n).description = ': Starts AI AutoMate';
n = n + 1;
c(n).switch = 'MatfileStorageFolder'; c(n).nargin = 1; c(n).default = ''; c(n).description = '<dir>: folder where output mat files are stored in the autoloader';
n = n + 1;
c(n).switch = 'Autoanalyse'; c(n).nargin = 1; c(n).default = '00000'; c(n).description = 'Perform automatic analysis in the autoloader';
n = n + 1;
c(n).switch = 'SaveToDatabase'; c(n).nargin = 0; c(n).default = false; c(n).description = 'Save to patient database';
n = n + 1;
c(n).switch = 'DatabasePath'; c(n).nargin = 1; c(n).default = ''; c(n).description = '<dir>: folder for patient database';
n = n + 1;
c(n).switch = 'SkipSecCapImages'; c(n).nargin = 0; c(n).default = false; c(n).description = 'Skip to load Secondary Captures in AI AutoMate';
n = n + 1;
c(n).switch = 'SkipScoutImages'; c(n).nargin = 0; c(n).default = false; c(n).description = 'Skip to load Scout images in AI AutoMate';
n = n + 1;
c(n).switch = 'Segment3DP'; c(n).nargin = 0; c(n).default = false; c(n).description = 'Run AI AutoMate for Segment 3DPrint';
n = n + 1;
c(n).switch = 'SendtoPACS'; c(n).nargin = 0; c(n).default = false; c(n).description = 'Send to PACS in AI AutoMate';

switch arg
  case 'parse'
    %do nothing here, do later after switch clause
  case 'getparameters'
    %return parameters
    varargout = cell(1,1);
    varargout{1} = s;
    return
  case 'storetopref'
    %store settings to prefstruct
    storetopref(s);
    return;
  case 'reset'
    %reset persistent variable
    s = [];
    for loop = 1:length(c)
      s.(c(loop).switch) = c(loop).default;
    end
  otherwise
    error('Invalid options to commandlinehelper');
end

%--- parse clause

if isempty(varargin)
  return
end

disp('------------------------------');    
disp('Parsing commandline options');

%Initialize s
s = [];
for loop = 1:length(c)
  s.(c(loop).switch) = c(loop).default;
end

%Add -help
n = length(c)+1;
c(n).switch = 'help'; c(n).nargin = 0; c(n).default = '';

%--- parse options
parameters = varargin;

while ~isempty(parameters)
  found = false;
  thisparameter = parameters{1}; %first of parameters
  if isequal(thisparameter,'-help')
    displayhelp(c);
    parameters = parameters(2:end);
  else
    %loop over to find
    for loop = 1:length(c)
      if isequal(['-' lower(c(loop).switch)],lower(thisparameter))
        found = true;
        numargs = c(loop).nargin;
        option = c(loop).switch;
        if numargs==1
          arg = cleanstring(parameters{2}); %more than one argument not yet implemented
          parameters = parameters(3:end);
          disp(sprintf(' * %s=%s',c(loop).switch,arg)); %#ok<DSPS>
        else
          arg = true; %to indicate that it is found
          parameters = parameters(2:end);
          disp(sprintf(' * %s=true',c(loop).switch)); %#ok<DSPS>
        end
        
      end
    end
    
    if found
      s.(option) = arg; %store it
    else
      %warn user
      disp(sprintf(' * option %s not recognized, ignored.',thisparameter)); %#ok<DSPS>
      parameters = parameters(2:end); %Remove parameters
    end

    
  end
  
end

disp('------------------------------');

%----------------------
function displayhelp(c)
%----------------------
%Displays helptext

disp(sprintf('Help information, version %s',changelog));  %#ok<DSPS>
disp(' ');
for loop = 1:(length(c)-1)
  disp(sprintf('-%s %s',lower(c(loop).switch),c(loop).description));  %#ok<DSPS>
  disp(' ');
end

%----------------------
function storetopref(s)
%----------------------
%Store settings to pref struct (DATA.Pref)

global DATA

if isempty(s)
  return
end

f = fieldnames(s);

for loop = 1:length(f)
  DATA.Pref.(f{loop}) = s.(f{loop});
end

%--------------------------------
function stri = cleanstring(stri)
%--------------------------------
%Cleans string from " signs

%Remove "
stri = stri(stri~='"');

%Remove '
c = '''';
stri = stri(stri~=c);
