function varargout = segmentopen(varargin)
%#ok<*GVMIS>

%%%% Main body %%%
fileinput = false;
if nargin > 0
  %Check if input appears to be a filename
  if regexp(varargin{1}, '^[c-zC-Z]:\\')
    fileinput = true;
  end
end
% LAUNCH and initalize GUI
programversion = changelog;
fig = initializesegmentopen(programversion); %Program version number

varargout = cell(1,nargout);
if nargout>0
  varargout{1} = fig;
end

% Load the mat file
if fileinput && ~isempty(fig)
  openfile('loadfiles', {varargin{1}}, false, []);
end

%-----------------------------------------------
function fig = initializesegmentopen(programversion)
%-----------------------------------------------
%Initialization of Segment CMR GUI
global DATA

if ~isempty(DATA)
  try
    mydisp('Already running.');
    if ~isempty(DATA.fig)
      figure(DATA.fig);
      fig = DATA.fig;
      return;
    else
      fig = []; %#ok<NASGU>
    end
  catch me
    disp('Program not aborted properly last time, all data will be lost and program restarted.');
    mydispexception(me);
  end
end

%Make sure fresh start
DATA = [];
SET = []; %#ok<NASGU>

%Get arguments from command line
s = commandlinehelper('getparameters');

if isfield(s,'NoGUI') && (s.NoGUI)
  %Start and quit software without displaying the GUI
  fig = [];
  return
end

%This is where we create object
DATA = opengui(programversion);%Load custom GUI



[segmentfolder_, ~, ~] = fileparts(which('segment'));
DATA.SegmentFolder = segmentfolder_;
clear segmentfolder_


DATA.init;
try
  fig = DATA.fig;
catch
  fig = [];
  logdisp('Software initialization aborted');
  return
end

checkpath(DATA.SegmentFolder); %ensure running on correct path
