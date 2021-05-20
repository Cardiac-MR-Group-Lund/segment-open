function varargout = plugin_imageloader(fcn,varargin)
% Plugin to load non DICOM images into Segment.
%
% Main plugin functions. Works as switchyard and also
% reports dependencies.

%Written by Jonatan Wulcan. Slightly commented by Einar Heiberg
%Enhanced by Einar Heiberg with tools to load multiple ones into separate
%image stacks.

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end

switch fcn
  case 'getname'
    varargout = cell(1,1);

    %Register callbacks
    uimenu(varargin{1},'Label','Load single file','Callback','plugin_imageloader(''loadfile_Callback'')');
    uimenu(varargin{1},'Label','Load files from directory to one stack','Callback','plugin_imageloader(''loaddir_Callback'')');
    uimenu(varargin{1},'Label','Load files from directory to separate stacks','Callback','plugin_imageloader(''loadmany_Callback'')');
    set(varargin{1},'Callback','');

    %Return title
    varargout{1} = 'Image Loader';
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
function loadfile_Callback %#ok<DEFNU>
%-------------------------
% Callback functions when user press the loadfile button in
% the menu. Ask the user for a filename and load that file.

global DATA SET NO

%selecting the folder
pathname = DATA.Pref.datapath;
%pathname = uigetdir(pathname,'Select a folder');

[filename,pathname] = myuigetfile(sprintf('%s%s*.*',pathname,filesep));
if filename == 0
  return; % User aborted
end

NO = length(SET)+1;

% load SET(NO).IM
SET(NO).IM = loadsinglegray(fullfile(pathname,filename));

% Set default values to DATA.Preview
resetpreview();
DATA.Preview.XSize = size(SET(NO).IM,1);
DATA.Preview.YSize = size(SET(NO).IM,2);

% Load it
openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%Does all graphical updates
if not(isempty(DATA.ViewMatrix))
  DATA.init_graphics;
end
 
%-------------------------
function loadmany_Callback %#ok<DEFNU>
%-------------------------
% Callback routine to load multiple files from a folder to separare image
% stacks.

global DATA SET NO

pathname = myuigetdir(DATA.Pref.datapath,'Select folder to load jpg files from');
if pathname == 0
  return; % User aborted
end
files = dir([pathname filesep '*.jpg']);
if isempty(files)
  files = dir([pathname filesep '*.JPG']);
end
if isempty(files)
  myfailed('Could not find any files.');
end

%Get factor
[f,ok] = mygetnumber('Enter upsampling/downsampling factor (<1 => downsample)','Resampling Factor',0.1,0,10);
if ~ok
  myfailed('Invalid upsampling/downsampling factor or aborted.');
  return;
end
  
h = waitbar(0,'Please wait.');
for fileloop = 1:length(files)
  NO = length(SET)+1;

  % load SET(NO).IM
  tempim = loadsinglegray([pathname filesep files(fileloop).name]);  
  
  %Downsample
  tempim = single(imresize(tempim,f,'bicubic'));
  
  SET(NO).IM = tempim;
  
  % Set default values to DATA.Preview
  resetpreview();
  DATA.Preview.XSize = size(SET(NO).IM,1);
  DATA.Preview.YSize = size(SET(NO).IM,2);
  
  % Load it
  openfile('setupstacksfromdicom',NO);
  segment('renderstacksfromdicom',NO);
  
  waitbar(fileloop/length(files),h);
end
close(h);

%Does all graphical updates
if not(isempty(DATA.ViewMatrix))
  DATA.init_graphics;
end
 
%------------------------
function loaddir_Callback %#ok<DEFNU>
%------------------------
% Load directory callback. Executes when user press the
% load direcory button in the menu. Asks user for a directory
% then loads all files in that dir as a single stack.

global DATA SET NO

pathname = myuigetdir(DATA.Pref.datapath,'Select folder to load');
if pathname == 0
  return; % User aborted
end
files = dir(pathname);

NO = length(SET)+1;

% load SET(NO).IM
im = single(zeros([0 0 0 0]));
z = 1;
for n=1:length(files)
  if isimage(files(n))
    % Load pixel_data
    pixel_data = loadsinglegray(fullfile(pathname, files(n).name));

    % Check image size
    if z ~= 1 && any(size(im(:,:, 1, 1)) ~= size(pixel_data))
      myfailed('Images must be of the same size');
      return
    end

    % Save pixel_data to SET(NO).IM
    im(:, :, 1, z) = pixel_data;
    z = z+1;
  end
end
SET(NO).IM = im;

% Set default values to DATA.Preview
resetpreview();
DATA.Preview.XSize = size(SET(NO).IM,1);
DATA.Preview.YSize = size(SET(NO).IM,2);
DATA.Preview.ZSize = size(SET(NO).IM,4);

% Load it
openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%Does all graphical updates
if not(isempty(DATA.ViewMatrix))
  DATA.init_graphics;
end
 
%-----------------------
function r = isimage(f)
%-----------------------
% Returns true if file seems to be a supported image
% else returns false.

r = false;
if f.isdir
  return;
end
if regexpi(f.name, '\.jpg$');
  r = true;
  return;
end
if regexpi(f.name, '\.tif$');
  r = true;
  return;
end
if regexpi(f.name, '\.png$');
  r = true;
  return;
end

%---------------------------------------
function pixel_data = loadsinglegray(f)
%---------------------------------------
% Loads an image and convert it to grayscale
% and class single.

pixel_data = imread(f);
if size(pixel_data, 3) == 3
  pixel_data = rgb2gray(pixel_data);
end

pixel_data = single(pixel_data);

