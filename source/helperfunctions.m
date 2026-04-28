function varargout = helperfunctions(varargin)
% Useful functions to use in the code

%#ok<*GVMIS>

if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-------------------------------------
function tf = gettimeframe(str,no)
%-------------------------------------
arguments
  str {mustBeMember(str,{'ed','es',''})}
  no = []
end
global SET

if isempty(no)
  global NO  %#ok<TLEV> 
  no = NO;
end

switch str
  case 'ed'
    tf = SET(no).EDT;
  case 'es'
    tf = SET(no).EST;
  otherwise
    tf = SET(no).CurrentTimeFrame;
end

%-------------------------------------
function [lvrvongoing,laraongoing,generalpenongoing] = isanyinterpongoing(no)
%-------------------------------------
%Helper function to check if any interpolation contour is ongoing.
global DATA SET

lvrvongoing = false;
laraongoing = false;
generalpenongoing = false;

%Check any other contour
if SET(no).EndoInterpOngoing || SET(no).EpiInterpOngoing || ...
    SET(no).RVEndoInterpOngoing || SET(no).RVEpiInterpOngoing
  lvrvongoing = true;
end

%Check LA/RA interpolation
if (~isempty(SET(no).LA) && SET(no).LA.InterpOngoing) || ...
    (~isempty(SET(no).RA) && SET(no).RA.InterpOngoing)
  laraongoing = true;
end

%Check General Pen interpolation. If any ongoing interpolation is detected,
%the loop can break.
for ind = 1:DATA.GeneralPenSettings.getnumobjects(no)
  if SET(no).GeneralPenObjects(ind).InterpOngoing
    generalpenongoing = true;
    break
  end
end

%-------------------------------------
function storeinterptocontour(no,panel,type,timeframe,slice,objectind,threshold)
%-------------------------------------
%Helper function to store the delineated interpolated contour to the normal
%contour. The contour is resampled, stored and drawn.
global SET

%Get interpolated contour
numpoints = tools('getnuminterppointsforno',no,type,slice,objectind);
x = helperfunctions('parsesetfield',SET(no),type,'X',objectind,timeframe,slice); %get interp contour
y = helperfunctions('parsesetfield',SET(no),type,'Y',objectind,timeframe,slice); %get interp contour

%Resample contour
% resamples the contour
isinterpongoing = helperfunctions('isinterpongoing',SET(no),type,objectind);
if ~isinterpongoing && length(x) > 1 && ~(any(contains(type,{'LA','RA'})))
  opencontour = false;
  closethecurve = true;
  numpointstoresample = numpoints - 1;
else
  opencontour = true;
  closethecurve = false;
  numpointstoresample = numpoints;
end
[x,y] = calcfunctions('resamplecurve',x,y,numpointstoresample,opencontour,closethecurve);

%write the results back to the contour field in the SET struct if we've
%placed more than N points (threshold)
if length(x) > threshold
  switch type
    case 'Roi'
      %would be nice to introduce RoiInterpX...
    otherwise
      if isinterpongoing
        typeX = helperfunctions('parsesetfield',SET(no),type(1:end-6),'X',objectind); %get contour
        expectedlength = length(typeX(:,timeframe,slice));
        if numpoints < expectedlength
          % fill up with NaNs
          x = cat(2,x,nan(1,expectedlength-numpoints));
          y = cat(2,y,nan(1,expectedlength-numpoints));
        end
      end

      if ~opencontour
        xvec = [x,x(1)]';
        yvec = [y,y(1)]';
      else
        xvec = x';
        yvec = y';
      end
      helperfunctions('assignsetfield',no,type(1:end-6),'X',xvec,objectind,timeframe,slice);
      helperfunctions('assignsetfield',no,type(1:end-6),'Y',yvec,objectind,timeframe,slice);

      %Draw contour
      switch type
        case 'GeneralPenInterp'
          drawfunctions('drawcontoursgeneralpen',panel);
        case {'LAInterp','RAInterp'}
          drawfunctions('drawcontourslara',panel,type(1:end-6));
        otherwise
          drawfunctions('drawcontours',panel,type(1:end-6));
      end
  end
elseif length(x) == threshold
  switch type
    case 'GeneralPenInterp'
      SET(no).GeneralPenObjects(objectind).InterpOngoing = true;
    case {'LAInterp','RAInterp'}
      SET(no).(type(1:end-6)).InterpOngoing = true;
    otherwise
      SET(no).([type(1:end-6),'InterpOngoing']) = true;
  end
  drawfunctions('updateinterpolationsettings',panel,type);
end

%----------------------------------------------
function [no,slice,timeframe,x,y] = getinterpparameters(panel)
%----------------------------------------------
%Helper function to find general parameters for interp functions
global DATA SET

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

%get closest interpolation point
[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',panel,y,x);
slices = viewfunctions('slicesinpanel',panel);
timeframe = SET(no).CurrentTimeFrame;

%normalize clicked position to contour domain
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
imdim = zoomfunctions.getxysize(no,panel);
xt = scale*((xl-1)*imdim.XSize- imdim.XStart +1);
yt = scale*((yl-1)*imdim.YSize- imdim.YStart +1);

x = (x-xt)/scale;
y = (y-yt)/scale;


%------------------------------------------------------------------------
function bool = isinterpongoing(setstruct,type,objectind)
%------------------------------------------------------------------------
%Function to parse the value of the field "InterpOngoing", independently
%from field depth. Useful for functions that work both on LV and RV
%contours (legacy) and on contours using the class clpen.
arguments
  setstruct %struct of the specified stack, i.e. SET(NO)
  type {mustBeMember(type,{'GeneralPenInterp','LAInterp','RAInterp',...
    'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})}
  objectind {mustBeInteger} = []; %used for general pen objects
end

%Parse the specified field
fullfield = getfullfield_helper(type,'Ongoing');

%Find field names
separator = strfind(fullfield,'.');
temp = setstruct;
while ~isempty(separator)
  temp = temp.(fullfield(1:separator(1)-1))(objectind);
  fullfield = fullfield(separator(1)+1:end);
  separator = strfind(fullfield,'.');
end

%Get field value
bool = temp.(fullfield);

%------------------------------------------------------------------------
function fieldvalue = parsesetfield(setstruct,type,xy,objectind,timeframe,slice)
%------------------------------------------------------------------------
%Function to parse field value, independently from field depth. Useful for
%functions that work both on LV and RV contours (legacy) and on contours
%using the class clpen.
arguments
  setstruct %struct of the specified stack, i.e. SET(NO)
  type {mustBeMember(type,{'GeneralPen','LA','RA','Endo','Epi','RVEndo','RVEpi', ...
    'GeneralPenInterp','LAInterp','RAInterp','EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})}
  xy {mustBeMember(xy,{'X','Y'})} %get X or Y coordinates
  objectind {mustBeInteger} = 1; %used for general pen objects
  timeframe {mustBeInteger} = [];
  slice {mustBeInteger} = [];
end

%Parse the specified field
fullfield = getfullfield_helper(type,xy);

%Find field names
separator = strfind(fullfield,'.');
temp = setstruct;
while ~isempty(separator)
  if ~isempty(temp.(fullfield(1:separator(1)-1)))
    temp = temp.(fullfield(1:separator(1)-1))(objectind);
    fullfield = fullfield(separator(1)+1:end);
    separator = strfind(fullfield,'.');
  else
    fieldvalue = [];
    return
  end
end

%Get field value
if ~isempty(timeframe) && ~isempty(slice)
  fieldvalue = temp.(fullfield){timeframe,slice};
else
  fieldvalue = temp.(fullfield);
end

%------------------------------------------------------------------------
function assignsetfield(no,type,xy,valuetoassign,objectind,timeframe,slice,ind)
%------------------------------------------------------------------------
%Function to assign field value, independently from field depth. Useful for
%functions that work both on LV and RV contours (legacy) and on contours
%using the class clpen.
arguments
  no
  type {mustBeMember(type,{'GeneralPen','LA','RA','Endo','Epi','RVEndo','RVEpi', ...
    'GeneralPenInterp','LAInterp','RAInterp','EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})}
  xy {mustBeMember(xy,{'X','Y'})} %assign to X or Y coordinates
  valuetoassign
  objectind {mustBeInteger} = 1; %used for general pen objects
  timeframe {mustBeInteger} = [];
  slice {mustBeInteger} = [];
  ind {mustBeInteger} = []; %specific point
end
global SET

%Parse the specified field
fullfield = getfullfield_helper(type,xy);
fields = split(fullfield,'.');

%Find field names
separator = strfind(fullfield,'.');
temp = SET(no);
while ~isempty(separator)
  temp = temp.(fullfield(1:separator(1)-1))(objectind);
  fullfield = fullfield(separator(1)+1:end);
  separator = strfind(fullfield,'.');
end

%Assign value to specified field
if numel(fields) > 1
  if ~isempty(timeframe) && ~isempty(slice)
    if iscell(temp.(fullfield))
      if ~isempty(objectind)
        if ~isempty(ind)
          SET(no).(fields{1})(objectind).(fields{2}){timeframe,slice}(ind,1) = valuetoassign;
        else
          SET(no).(fields{1})(objectind).(fields{2}){timeframe,slice} = valuetoassign;
        end
      else
        if ~isempty(ind)
          SET(no).(fields{1}).(fields{2}){timeframe,slice}(ind,1) = valuetoassign;
        else
          SET(no).(fields{1}).(fields{2}){timeframe,slice} = valuetoassign;
        end
      end
    else
      if ~isempty(objectind)
        SET(no).(fields{1})(objectind).(fields{2})(:,timeframe,slice) = valuetoassign;
      else
        SET(no).(fields{1}).(fields{2})(:,timeframe,slice) = valuetoassign;
      end
    end
  else
    if ~isempty(objectind)
      SET(no).(fields{1})(objectind).(fields{2}) = valuetoassign;
    else
      SET(no).(fields{1}).(fields{2}) = valuetoassign;
    end
  end
else
  if ~isempty(timeframe) && ~isempty(slice)
    if iscell(temp.(fullfield))
      if ~isempty(ind)
        SET(no).(fields{1}){timeframe,slice}(ind,1) = valuetoassign;
      else
        SET(no).(fields{1}){timeframe,slice} = valuetoassign;
      end
    else
      SET(no).(fields{1})(:,timeframe,slice) = valuetoassign;
    end
  else
    SET(no).(fields{1}) = valuetoassign;
  end
end

%------------------------------------------------------------------------
function fullfield = getfullfield_helper(type,xy)
%------------------------------------------------------------------------
%Helper function for parsesetfield and assignsetfield.
arguments
  type {mustBeMember(type,{'GeneralPen','LA','RA','Endo','Epi','RVEndo','RVEpi', ...
    'GeneralPenInterp','LAInterp','RAInterp','EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})}
  xy {mustBeMember(xy,{'X','Y','Ongoing'})} %get X or Y coordinates or InterpOngoing flag
end

switch type
  %need to parse property under field
  case 'GeneralPen'
    fullfield = ['GeneralPenObjects.' xy];
  case 'LA'
    fullfield = ['LA.' xy];
  case 'RA'
    fullfield = ['RA.' xy];
  case 'GeneralPenInterp'
    fullfield = ['GeneralPenObjects.Interp' xy];
  case 'LAInterp'
    fullfield = ['LA.Interp' xy];
  case 'RAInterp'
    fullfield = ['RA.Interp' xy];
  otherwise
    %EndoInterp,EpiInterp,RVendoInterp,RVepiInterp
    %no need to parse property under field
    fullfield = [type xy];
end

%------------------------------------------------------------------------
function [mostrecentstudydate, studytime] = getmostrecentstudydate(setstruct)
%------------------------------------------------------------------------
%Helper function to find the most recent study date in a dataset. This is
% useful in case the datasat contains several studies.
arguments
  setstruct = [];
end

global SET

if isempty(setstruct)
  setstruct = SET;
end

%if no dates are provided, then look on the whole dataset
if numel(setstruct) > 1
  datecell = arrayfun(@(x) x.PatientInfo.AcquisitionDate, setstruct, 'UniformOutput', false);
else
  %Only one stack
  mostrecentstudydate = setstruct(1).PatientInfo.AcquisitionDate;
  studytime = setstruct(1).AcquisitionTime;
  return
end

% Find indices of empty cells
emptyIndices = cellfun(@isempty, datecell);
datecell = datecell(~emptyIndices); %delete empty cells
if isempty(datecell)
  mostrecentstudydate = '';
  studytime = '';
  return
end
%sort the dates
studydatesarray = datetime(datecell,'InputFormat','yyyyMMdd','Format','yyyyMMdd');
[~, studyind] = sort(studydatesarray,'descend');
ind = studyind==1;

%get the most recent one and the corresponding study time
mostrecentstudydate = removeforbiddenchars(char(studydatesarray(ind)));
studytime = setstruct(ind).AcquisitionTime;

%------------------------
function changelinecolorfromyellowtoorange(parentaxes)
%------------------------
%Helper function to find all children in an axes that are ploted in yellow
%and change the color to orange. This is a work-around for plotting Flow
%results based on the current implementation of ROIs.
global DATA

yellowlines = findobj('Parent',parentaxes,'Color','y');
set(yellowlines,'Color',DATA.Handles.yroi(1).Color);

%------------------------
function setfdabuttongroupcolor(backgroundcolor,foregroundcolor)
%------------------------
%Set background and foreground colors for the radiobuttons that controls
%the FDA cleared version of Segment. Called in segpref('setbackgroundcolor').
global DATA
h = DATA.Handles;
set([h.fdabuttongroup, h.fdaradiobutton, h.researchradiobutton], ...
  'BackgroundColor',backgroundcolor,'ForegroundColor',foregroundcolor);

%--------------
function getcitation
%--------------
%write correct citation for Segment into a .txt. file

% open file to write in
pathname = getpreferencespath;
txtfile = [pathname filesep 'segmentcitation.txt'];
citefid = fopen(txtfile,'wt');
% write program version
[~,~,~,version] = changelog;
txt = sprintf('All image analysis was done using the freely available software Segment %s (Medviso, Lund, Sweden) [1]',version);
fprintf(citefid,sprintf('%s\n',txt));

citetxt = ['[1] E. Heiberg, J. Sjögren, M. Ugander, M. Carlsson, H. Engblom, and H. Arheden, ',newline,...
  'Design and Validation of Segment – ',...
  'a Freely Available Software for Cardiovascular Image Analysis, ',newline,...
  'BMC Medical Imaging, 10:1, 2010.'];
fprintf(citefid,sprintf('%s\n',citetxt));

homepage = 'https://medviso.com/how-to-refer/';
txt = sprintf('For more details see also %s',homepage);
fprintf(citefid,sprintf('\n\n%s\n',txt));

fclose(citefid);
myopenfile(txtfile);

%--------------
function checkforundeletedfiles(filename)
%--------------
%check if delete failed
if exist(filename, 'file')
  warnstr = dprintf('Could not delete temporary file!\n%s\nIMPORTANT: Contains patient data.',filename);
  myfailed(warnstr);
  stri = sprintf('explorer.exe /select,%s', filename);
  system(stri);
end

%--------------
function devicemodelname = getdevicemodelname
%--------------
% function to get full device name
global DATA

devicemodelname = [DATA.ProgramName,' ',DATA.ProgramVersion];

%--------------
function seriesdescription = getseriesdescriptionfordevice
%--------------
% get series descirption for device
global DATA

seriesdescription = [DATA.ProgramName,' data'];

%------------------------------------
function romnum = getromannumber(num)
%------------------------------------
% convert number to roman number in string format
romansymbols = {'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};
values = [100, 90, 50, 40, 10, 9, 5, 4, 1];

romnum = '';
for loop = 1:length(values)
  while num >= values(loop)
    romnum = [romnum, romansymbols{loop}]; %#ok<AGROW>
    num = num - values(loop);
  end
end

%------------------------------
function [px, py, figwidth, figheight] = getfigposition(figwidth,figheight,parenthandle)
%------------------------------
% get figure position based on requested height and width, that does
% positioning in the center of parent handle and does not go outside the screen
arguments
  figwidth {mustBeNumeric}
  figheight {mustBeNumeric}
  parenthandle = [];
end

global DATA
if ~isempty(DATA)
  if isempty(parenthandle)
    screensize = DATA.GUI.Segment.fig.Position;
  elseif isa(parenthandle,'mygui')
    screensize = parenthandle.fig.Position;
  else
    screensize = get(groot,"ScreenSize");
  end
else
  screensize = get(groot,"ScreenSize");
end
taskbar = 40; % windows taskbar size
titlebar = 25; % size of gray title bar in matlab figures
screenwidth = screensize(3);
screenheight = screensize(4)-taskbar-titlebar; %height-taskbar-titlebar
if figheight > screenheight
  figheight = screenheight;
end
if figwidth > screenwidth
  figwidth = screenwidth;
end
% get the starting point of the figure
px = floor((screenwidth-figwidth)/2)+screensize(1);
py = floor((screenheight-figheight)/2)+screensize(2)+taskbar;


%--------------------------------------------
function deleteallfilesfromdirectory(dirname)
%--------------------------------------------
% fucntion to remove all files from provided directory
% can be used, when removing the whole directory fails

files = dir(fullfile(dirname, '*'));

% loop over each file in the directory
for loop = 1:numel(files)
  % Check if the file is not a directory
  if ~files(loop).isdir
    filepath = fullfile(dirname, files(loop).name);
    % Delete the file
    delete(filepath);
  end
end

%------------------------
function user = getuser()
%------------------------
%Gets current user

persistent storeduser %Use persistent to avoid to many calls to getenv

if isempty(storeduser)
  user = getenv('username');
  if ~isempty(user)
    user(1) = upper(user(1));
  end
  storeduser = user;
else
  user = storeduser;
end

%--------------------------------------------------------
function childpos = getfigposbasedonparent(parentfig,childfig)
%--------------------------------------------------------
% get position of the child figure in pixels, that would place child figure
% in the middle above the parent figure
childpos = getpixelposition(childfig);
if isempty(parentfig)
  % take original child position, since parent figure is empty
  return
end
parentpos = getpixelposition(parentfig);


posx = parentpos(1) + (parentpos(3) - childpos(3))/2;
posy = parentpos(2) + (parentpos(4) - childpos(4))/2;

childpos(1) = posx;
childpos(2) = posy;


%---------------------------------
function tf = iswindows11
%---------------------------------
%Function to find if the OS running Segment is Win11
osstr = system_dependent('getos');
inddigit = isstrprop(osstr, 'digit');
verstr = osstr(inddigit);
tf = false;
if ~isempty(verstr)
  if strcmp(verstr(1:2),'11')
    tf = true;
  end
end

%--------------------------------------------------------
function [ver,subver,releaseid] = parseversionstr(verstr)
%--------------------------------------------------------
%Function to parse the software version number
tt = sscanf(verstr, '%d.%d.%d.%d R%d%s');
if ~(numel(tt) == 5 || numel(tt) == 6)
  error('SEGMENT:ERROR', 'Couldn''t parse version string');
end
ver = tt(5);
if numel(tt) == 6
  subver = char(tt(6));
else
  subver = 'a';
end
releaseid = regexp(verstr,'\d\.\d\.\d\.\d','match','once');

%------------------------------------
function answ = checkifgpuuptodate
%------------------------------------
answ = true;

infostruct = rendererinfo;
if isfield(infostruct,'Vendor') && contains(infostruct.Vendor,'NVIDIA',IgnoreCase=true)
  if isfield(infostruct,'Details')
    infodetails = infostruct.Details;
    if isfield(infodetails,'RendererDriverReleaseDate')
      mydisp(['Driver Release Date: ',infodetails.RendererDriverReleaseDate])
      driverdate = datetime(infodetails.RendererDriverReleaseDate, 'InputFormat', 'yyyy-MM-dd');
      ageindays = days(datetime('today')-driverdate);
      if ageindays > 365 % more than year
        answ = false;
      end
    end
  end
end

%--------------------------------------
function fps = getdefaultfps(imagetype)
%--------------------------------------
if nargin == 0 || isempty(imagetype)
  imagetype = 'cine';
end
if contains(imagetype,'perfusion',IgnoreCase=true)
  fps = 10;
else
  fps = 25;
end

%------------------------------------
function orangecolor = getorangecolor
%------------------------------------
% get RGB values for ornage color in measurements
orangecolor = [1 0.5 0];

%-----------------------
function resetrandomseed
%-----------------------
% Reset seed for MATLAB's random number generator using current time for non-repetitive,
% non-predictable random number creation, used for example field StackUID
% in SET struct
reset(RandStream.getGlobalStream,sum(100*clock));

%----------------------------------
function ok = checknvfile(nvtofind)
%----------------------------------
% Checks if the NV file for a given Name exists on path

global DATA
structArray = DATA.Networks;
ok = false;

% Find index of struct with matching Name
idx = find(strcmp({structArray.Name}, nvtofind), 1);

if isempty(idx)
  logdisp('"%s" not found', nvtofind);
  return;
end
dirpath = getsegmentdirpath;
% Construct NV path
nvfilename = [dirpath,filesep,'nv',filesep,structArray(idx).NV{1}];
% Check if path exist
ok = isfile(nvfilename);

%-----------------------------------
function dirpath = getsegmentdirpath
%----------------------------------
% Function to get directory path to file storage
global DATA
if isdeployed
  foldername = DATA.getsoftwarenamenospace;
  if length(foldername) > 12
    foldername = foldername(1:12);
  end
  dirpath = [ctfroot,filesep,foldername];
else
  dirpath = DATA.SegmentFolder;
end


