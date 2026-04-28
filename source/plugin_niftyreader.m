function varargout = plugin_niftyreader(fcn,varargin)
% Plugin to load NIFTY files into Segment.

%Einar Heiberg

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end

switch fcn
  case 'getname'
    varargout = cell(1,1);

    %Register callbacks
    uimenu(varargin{1},'Label','Load nifti file','Callback','plugin_niftyreader(''loadnifty_Callback'')');
    set(varargin{1},'Callback','');

    %Return title
    varargout{1} = 'Nifti Image Loader';
  case 'getdependencies'
    %Here: List all depending files. This is required if your plugin should
    %be possible to compile to the stand-alone version of Segment.
    varargout = cell(1,4);

    %M-files
    varargout{1} = {};

    %Fig-files
    varargout{2} = {};

    %Mat-files
    varargout{3} = {};

    %Mex-files
    varargout{4} = {};

  otherwise
    macro_helper(fcn,varargin{:});
    [varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard
end

%---------------------
function resetpreview
%---------------------
%Resets The DATA.Preview struct. Fills in some general information.

global DATA

segment('resetpreview');

DATA.Preview.PatientInfo.Name = 'Image loaded via imageloader plugin';
DATA.Preview.PatientInfo.ID = '';
DATA.Preview.PatientInfo.BirthDate = '';
DATA.Preview.PatientInfo.Sex = '';
DATA.Preview.PatientInfo.Age = '';
DATA.Preview.HeartRate = 60;
DATA.Preview.PatientInfo.AcquisitionDate = datestr(now);
DATA.Preview.PatientInfo.BSA = 0;
DATA.Preview.PatientInfo.BMI = 0;
DATA.Preview.PatientInfo.Weight = 0;
DATA.Preview.PatientInfo.Length = 0;
DATA.Preview.AcquisitionTime = '00.00';
DATA.Preview.ResolutionX = 1;
DATA.Preview.ResolutionY = 1;
DATA.Preview.SliceThickness = 5;
DATA.Preview.SliceGap = 0;
DATA.Preview.TIncr = 1;
DATA.Preview.TriggerTime = 0;
DATA.Preview.Bitstored = 16;

%-------------------------
function loadnifty_Callback %#ok<DEFNU>
%-------------------------
% Callback functions when user press the loadfile button in
% the menu. Ask the user for a filename and load that file.

global DATA SET NO

%selecting the folder
pathname = DATA.Pref.datapath;
%pathname = uigetdir(pathname,'Select a folder');

[filename,filepath] = myuigetfile(sprintf('%s%s*.nii',pathname,filesep));

if filename == 0
  return; % User aborted
end

%Construct path 
f = fullfile(filepath, filename); %Mod by Julius

%Get info on the nifty file %J
info = niftiinfo(f); 

%Read the volme
vol = double(niftiread(f)); %Mod by Julius
% if ndims(vol) == 4 
%   vol = permute(vol, [2, 1, 3, 4]); 
% else
%   vol = permute(vol, [2, 1, 3]); %Mod by Julius
% end

%Normalize
vol = uint8(255*(vol - min(vol, [], 'all'))./max(vol, [], 'all')); %Mod by Julius

%Assign into Segment
NO = length(SET)+1;

% load SET(NO).IM
SET(NO).IM = single(permute(vol,[1 2 4 3])); % As in Segment's data format X*Y*T*Z; NIFTII format is X*Y*Z*T

% Set default values to DATA.Preview
resetpreview();
DATA.Preview.XSize = size(SET(NO).IM,1);
DATA.Preview.YSize = size(SET(NO).IM,2);
DATA.Preview.TSize = size(SET(NO).IM,3);
DATA.Preview.ZSize = size(SET(NO).IM,4);
DATA.Preview.ResolutionX = info.PixelDimensions(2);
DATA.Preview.ResolutionY = info.PixelDimensions(1);
DATA.Preview.SliceThickness = info.PixelDimensions(3);
DATA.Preview.TimeVector = linspace(0,1,DATA.Preview.TSize);

% Load it
openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%Does all graphical updates
if not(isempty(DATA.ViewMatrix))
  DATA.init_graphics;
end
 
