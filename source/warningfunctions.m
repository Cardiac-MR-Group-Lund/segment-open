function [varargout] = warningfunctions(varargin)
% Functions for collecting warnings

%Fanny Månefjord, Medviso, 2023

%#ok<*GVMIS>

%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%----------------------------
function addwarning(s)
%----------------------------
% Add warning to DATA.Warnings with corresponding StackID
%Used in AI AutoMate

global DATA

if contains(s, 'For more details open log') %don't include line about log file in warning list
  s = s(1 : end - 76);
end

s = erase(s, char(10)); %delete newline
troubleloading = contains(s, 'Could not load');

if ~isempty(DATA.StackID) && ~troubleloading
  warningstruct.stackID = DATA.StackID;
elseif troubleloading && length(s) == 37 %no info about stack
  return;
else
  warningstruct.stackID = 0;
end
warningstruct.warnmsg = s;

if isempty(DATA.Warnings)
  DATA.Warnings = warningstruct; %first time
else
  DATA.Warnings(end + 1) = warningstruct;
end

%----------------------------
function clearwarning(s)
%----------------------------
% Clear warning s if it exists in DATA.Warnings

global DATA

warnings = DATA.Warnings;
if isempty(warnings)
  return
end

if any(contains({warnings.warnmsg}, s))
  index = contains({warnings.warnmsg}, s);
  warnings(index) = [];
end
DATA.Warnings = warnings;

%----------------------------
function clearallwarnings
%----------------------------
% Clear all warnings in DATA.Warnings

global DATA

DATA.Warnings = [];

%----------------------------
function showpatientdbwarnings(titlestr)
%----------------------------
% Show database warnings in a mymsgbox

global DATA;

title = dprintf(titlestr);
str = [];

for w = 1 : length(DATA.Warnings)
  if ~isempty(DATA.Warnings(w))
    warning = DATA.Warnings(w).warnmsg;
    if isempty(str)
      str = warning;
    else
      if ~contains(str, warning) %not many of the same warning
        str = [str newline  '- ' warning];
      end
    end
  end
end
if isempty(str)
  str = dprintf('No warnings to display');
end

warningstring = convertCharsToStrings(str); %string to get scroll box in myinfostruct
s(1).Text = warningstring;
myinfostruct(s, title, 85);

%----------------------------
function str = makewarninglist(warnings)
%----------------------------
% Make pretty list of warnings for AI AutoMate
global SET

setIDs = {SET.StackID};
str = [];
len = length(warnings);
for w = 1 : len
  id = warnings(w).stackID;
  stackno = find(cellfun(@(x) isequal(x, id), setIDs));
  if ~isempty(stackno)
    stackstr = num2str(stackno); %to do: prettier list of stacks, 1-6, 8, 11
    warningstr = sprintf('%s %s: \n%s \n\n', dprintf('Stack'), stackstr, warnings(w).warnmsg);
  else
    warningstr = sprintf('%s \n \n', warnings(w).warnmsg);
  end

  if isempty(str) %first time
    str = warningstr;
  end
  if ~contains(str, warningstr) %don't add the same warning
    str = [str warningstr];
  end
end

%----------------------------
function movewarningstoDATA
%----------------------------
% Move autoloader warnings from SET to DATA

global DATA SET

if SET(1).Autoloader.Autoloaded
  if ~isempty(SET(1).Autoloader.Warnings) %warnings exists
    DATA.Warnings = SET(1).Autoloader.Warnings;
  end
end

%----------------------------
function showautoloaderwarnings(titlestr, autoanalysed)
%----------------------------
% Show autoloader warnings in a mymsgbox

global SET DATA;

if nargin<1
  title = dprintf('AI AutoMate notifications');
else
  title = dprintf(titlestr);
end
warnings = makewarninglist(SET(1).Autoloader.Warnings);
warningstring = convertCharsToStrings(warnings);
if nargin < 2
  autoanalysed = SET(1).Autoloader.Autoanalysed;
end
autoloaded = SET(1).Autoloader.Autoloaded;
n = 1;
if autoloaded
  if autoanalysed
    analysedstr = dprintf('This file was automatically loaded and segmented by AI AutoMate.');
    s(n).Text = analysedstr;
    n = n + 1;
    if contains(DATA.ProgramName, 'CMR')
      approvestr = [dprintf('Please verify the segmentations by approving them.')];
      s(n).Text = approvestr;
      n = n + 1;
    end

    if ~contains(DATA.ProgramName, '3DP')
      %review info
      if SET(1).Autoloader.Approved
        str1 = dprintf('Approved by:');
        str2 = dprintf('at');
        reviewstr = sprintf('%s %s %s %s',str1,SET(1).Autoloader.ApprovedBy,str2,SET(1).Autoloader.ApprovedTime);
      else
        reviewstr = [dprintf('Not approved.')];
      end
      s(n).Text = reviewstr;
      n = n + 1;
    end

  else
    analysedstr = dprintf('This file was automatically loaded by AI AutoMate.');
    s(n).Text = analysedstr;
    n = n + 1;
  end
else
  mymsgbox(dprintf('Not an AI AutoMate file'));
  return
end

if ~isempty(warningstring)
  s(n).Text = dprintf('DICOM file loading:');
  n = n + 1;
  s(n).Text = warningstring;
  amountwarnings = count(warningstring, ':');
  height = min(amountwarnings*3, 16);
  s(n).Height = height;
end

pref = DATA.Pref.DoNotAsk;
DATA.Pref.DoNotAsk = false;

if n > 1
  myinfostruct(s, title, 82); %82 is the length of analysedstr
else
  s(n).Text = dprintf('This file was automatically loaded by AI AutoMate.');
  myinfostruct(s, title, 35);
end

DATA.Pref.DoNotAsk = pref;
