function [varargout] = commandlinehelper(arg,varargin)
%Helper function to handle commandline options

%Einar Heiberg

persistent s %store options in persistent variable

%--- Define commandline options
c = [];
c(1).switch = 'StudiesFolder'; c(1).nargin = 1; c(1).default = ''; c(1).description = '<dir>: folder software reads input files (dicoms), each subfolder is one series';
c(2).switch = 'DICOMStorageFolder'; c(2).nargin = 1; c(2).default = ''; c(2).description = '<dir>: folder where output dicom files are stored';
c(3).switch = 'LogFolder'; c(3).nargin = 1; c(3).default = ''; c(3).description = '<dir>: folder where log-files are stored';
c(4).switch = 'NoQuit'; c(4).nargin = 0; c(4).default = false; c(4).description = ': Disable quit options';
c(5).switch = 'NoPopUp'; c(5).nargin = 0; c(5).default = false; c(5).description = ': Disable pop-up when starting';
c(6).switch = 'License'; c(6).nargin = 1; c(6).default = ''; c(6).description = '<file>: License file';
c(7).switch = 'TempStorageFolder'; c(7).nargin = 1; c(7).default = ''; c(7).description = '<dir>: folder for temporal file storage';
c(8).switch = 'UseDefaultStudiesFolder'; c(8).nargin = 0; c(8).default = false; c(8).description = '<dir>: Create studies folder if it not existing according to default preference settings';
c(9).switch = 'UseComputerAETitle'; c(9).nargin = 0; c(9).default = false; c(9).description = '<dir>: Use computer AETitle in case if in segpref was different one assigned';
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
