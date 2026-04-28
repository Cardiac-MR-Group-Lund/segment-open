function [varargout] = slicelocationtool(varargin)
%GUI to fix slice location problems

%Einar Heiberg

%#ok<*GVMIS> 

if nargin==0
  varargin = {'init_Callback'};
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%---------------------
function init_Callback 
%---------------------
%Init GUI
global DATA 

logdisp('Slice location tool');

%Check if File Loader is open, otherwise throw error
try
  files2load = openfile('getfiles2load');
catch
  disp('Slice location tool: Could not find files to load. File Loader might not be open.');
  return;
end

%Analyse data
filenames = cell(1,length(files2load));
slicelocations = zeros(1,length(files2load));
imageorientations = filenames;
imagepositions = filenames;
triggertimes = zeros(size(filenames));
dinfocell = filenames;
gantrytilt = zeros(size(filenames));
ignored = false(size(filenames));

h = mywaitbarstart(length(files2load),'Please wait',10,DATA.GUI.OpenFile); %length(files2load));
for loop = 1:length(files2load)
  try
    filenumber = openfile('getfilenumber',files2load{loop});
    filenames{loop} = sprintf('%d',filenumber);
    if length(filenames{loop})>10
      filenames{loop} = filenames{loop}((end-4):end);
    end
    dinfo = dicominfo(files2load{loop});
    dinfocell{loop} = dinfo;
    try
      gantrytilt(loop) = dinfo.GantryDetectorTilt;
    catch
    end

    try
      imageorientations{loop} = dinfo.ImageOrientation;
    catch
      imageorientations{loop} = dinfo.ImageOrientationPatient;
    end

    try
      imagepositions{loop} = dinfo.ImagePosition;
    catch
      imagepositions{loop} = dinfo.ImagePositionPatient;
    end

    try
      triggertimes(loop) = dinfo.TriggerTime;
    catch
      triggertimes(loop) = 0;
    end

    if isequal(rem(loop,10),0)
      h = mywaitbarset(loop/length(files2load),h);
    end
  catch
    %Something went wrong, ignore this file
    ignored(loop) = true;
    logdisp(sprintf('Problems parsing file, ignoring %s',files2load{loop}));
  end

end

mywaitbarclose(h);

%--- Remove ignored files
files2load = files2load(~ignored);
slicelocations = slicelocations(~ignored);
imageorientations = imageorientations(~ignored);
imagepositions = imagepositions(~ignored);
triggertimes = triggertimes(~ignored); %Seems like we do not use it..
gantrytilt = gantrytilt(~ignored);
dinfocell = dinfocell(~ignored);

sortedtriggertimes = sort(triggertimes);

if any(diff(sortedtriggertimes) > 0)
  mywaitbarclose(h);
  failstr = sprintf('%s.\n%s!',dprintf('Time resolved images are not supported'),dprintf('Delete all timeframes except this'));
  myfailed(failstr);
  return
end

if any(abs(gantrytilt)>0)
  y = yesno('Gantry tilt is not supported, use "Gantry Tilt Tool" to convert images. Use anyway?');
  if ~y
    mywaitbarclose(h);
    myworkoff;
    return
  end
end

%--- Compute slicelocations
for loop = 1:length(files2load)
  zdir = cross(imageorientations{loop}(1:3),imageorientations{loop}(4:6));
  slicelocations(loop) = sum(imagepositions{loop}.*zdir);
end

%Create figure
DATA.GUI.SliceLocationTool = mygui('slicelocationtool.fig');
gui = DATA.GUI.SliceLocationTool;

%Store to gui
gui.page = 1;
gui.files2load = files2load;
gui.filenames = filenames;
gui.filenumber = filenumber;
gui.slicelocations = slicelocations;
gui.imagepositions = imagepositions;
gui.imageorientations = imageorientations;
gui.dinfocell = dinfocell;
gui.included = true(1,length(files2load));
gui.timeframes = guesstimeframes;
gui.nfiles = length(files2load);

%updated later
gui.minincluded = min(slicelocations);
gui.maxincluded = max(slicelocations);
gui.minincludedind = 1;
gui.maxincludedind = length(files2load);
update;

%------------------------------------
function timeframes = guesstimeframes
%------------------------------------
%Guess on number of timeframes
global DATA

gui = DATA.GUI.SliceLocationTool;

discretelocations = union(gui.slicelocations,[]);

numatlocation = nan(1,length(discretelocations));
for loop = 1:length(discretelocations)
  numatlocation(loop) = sum(gui.slicelocations==discretelocations(loop));
end

timeframes = median(numatlocation); 

%check
if isequal(rem(length(gui.files2load),timeframes),0)
else
  timeframes = 1;
end

%--------------
function update
%--------------
%Update each page

global DATA

gui = DATA.GUI.SliceLocationTool;

set(gui.handles.plotaxes,'visible','on');

switch gui.page
  case 1
    updatepage1;
  case 2
    updatepage2;
  case 3
    updatepage3;
  case 4
    updatepage4;    
end

%-------------------
function updatepage1
%-------------------
%View page 1

global DATA

gui = DATA.GUI.SliceLocationTool;
gui.ind = 1:length(gui.files2load); %not sorted

%graphical
headlinetext = sprintf('%s %d/%d',dprintf('Step'),1,4);
infotext = [...
  dprintf('Plot shows slice location versus files (not sorted).') ' ' ...
  dprintf('Zoom to see more details.')];
set(gui.handles.previouspushbutton,'Enable','on');

set(gui.handles.headlinetext,'String',headlinetext);
set(gui.handles.infotext,'String',infotext);
set(gui.handles.timeframestext,'String',dprintf('Time frames'));
set(gui.handles.timeframesedit,'String',sprintf('%d',gui.timeframes));

%get local variables
bc = DATA.GUISettings.BackgroundColor;
fc = DATA.GUISettings.ForegroundColor;
ax = gui.handles.plotaxes;

slicelocations = gui.slicelocations;

plot(ax,NaN,NaN); %empty
hold(ax,'on');
for loop = 1:length(gui.files2load)
  h = plot(ax,loop,slicelocations(loop),'o'); %Plot them one by one
  if gui.included(loop)
    set(h,'Color',fc,'UserData',loop);
  else
    set(h,'Color',[0.3 0.3 0.3],'UserData',loop);
  end
end
hold(ax,'off');

set(ax,'Color',bc,'XColor',fc,'YColor',fc);
ylabel(ax,dprintf('SliceLocation'),'Color',fc);
xlabel(ax,dprintf('Filenumber in list'),'Color',fc);
zoom(gui.fig,'on');

%-------------------
function updatepage2
%-------------------
%View page 1, sorted on slice location

global DATA

gui = DATA.GUI.SliceLocationTool;

slicelocations = gui.slicelocations;
[sortedslicelocations,gui.ind] = sort(slicelocations); %sorted

gui.minincluded = min(slicelocations(gui.included));
gui.maxincluded = max(slicelocations(gui.included));
gui.minincludedind = find(sortedslicelocations==gui.minincluded,1,'first');
gui.maxincludedind = find(sortedslicelocations==gui.maxincluded,1,'first');

%graphical
headlinetext = sprintf('%s %d/%d',dprintf('Step'),2,4);
infotext = [...
  dprintf('Plot shows slice location versus sorted files.') ' ' ...
  dprintf('Color line show estimated slice location.')];
set(gui.handles.headlinetext,'String',headlinetext);
set(gui.handles.infotext,'String',infotext);
set(gui.handles.previouspushbutton,'Enable','on');

%get local variables
bc = DATA.GUISettings.BackgroundColor;
fc = DATA.GUISettings.ForegroundColor;
ax = gui.handles.plotaxes;

plot(ax,NaN,NaN); %empty
hold(ax,'on');
if gui.timeframes>1
  for loop = 1:round(gui.nfiles/gui.timeframes)
    h = plot(ax,(loop-1)*gui.timeframes+[1 gui.timeframes],[sortedslicelocations(1+(loop-1)*gui.timeframes) sortedslicelocations(loop*gui.timeframes)],'-');
    switch rem(loop,3)
      case 0
        set(h,'Color',[1 1 0],'Linewidth',2);
      case 1
        set(h,'Color',[1 0 0],'Linewidth',2);
      case 2
        set(h,'Color',[0 1 0],'Linewidth',2);
    end
  end
else
  h = plot(ax,[gui.minincludedind gui.maxincludedind],[gui.minincluded gui.maxincluded],'-'); set(h,'Color',[1 1 0],'Linewidth',2);
end

%Plot all points over it
for loop = 1:gui.nfiles
  h = plot(ax,loop,slicelocations(gui.ind(loop)),'o'); 
  if gui.included(gui.ind(loop))
    set(h,'Color',fc,'UserData',gui.ind(loop));
  else
    set(h,'Color',[0.3 0.3 0.3],'UserData',gui.ind(loop));
  end
end
hold(ax,'off');

set(ax,'Color',bc,'XColor',fc,'YColor',fc);
ylabel(dprintf('SliceLocation'));
xlabel(dprintf('File'));

%-------------------
function updatepage3
%-------------------
%Page 3, difference to regression line

global DATA

gui = DATA.GUI.SliceLocationTool;
sortedslicelocations = sort(gui.slicelocations);

%graphical
headlinetext = sprintf('%s %d/%d',dprintf('Step'),3,4);
infotext = [...
  dprintf('Plot shows difference to estimated slice location.') ' ' ...
  ];

set(gui.handles.headlinetext,'String',headlinetext);
set(gui.handles.infotext,'String',infotext);

%y = kx+m, x=loop-1
m = gui.minincluded;
k = (gui.maxincluded-gui.minincluded)/(gui.maxincludedind-gui.minincludedind);

%get local variables
bc = DATA.GUISettings.BackgroundColor;
fc = DATA.GUISettings.ForegroundColor;
ax = gui.handles.plotaxes;

%Loop to plot
cla(ax,'reset');
hold(ax,'on');
ylim = 0.1;
for loop = 1:length(gui.files2load)

  %Compute slice
  if gui.timeframes>1
    slice = floor((loop-1)/gui.timeframes);
  else
    slice = loop-1;
  end
  
  %Compute distance
  y = k*slice*gui.timeframes+m;
  d = y-sortedslicelocations(loop);
  ylim = max(ylim,d);
  
  h = plot(gui.handles.plotaxes,loop,d,'ko');
  if gui.included(gui.ind(loop))
    set(h,'Color',fc,'UserData',gui.ind(loop));
  else
    set(h,'Color',[0.3 0.3 0.3],'UserData',gui.ind(loop));
  end
end
hold(ax,'off')

set(ax,'ylim',[-ylim,ylim]);
set(ax,'Color',bc,'XColor',fc,'YColor',fc);
ylabel(ax,'SliceLocation','Color',fc);
xlabel(ax,'Filenumber in list','Color',fc);

paramstr = dprintf('Difference');
labelstr = makeunitstring(paramstr,'mm');
ylabel(labelstr);

xlabel(dprintf('File'));

%-------------------
function updatepage4
%-------------------
%This actually fix the files

global DATA

gui = DATA.GUI.SliceLocationTool;

headlinetext = sprintf('%s %d/%d',dprintf('Step'),4,4);
infotext = '';
set(gui.handles.headlinetext,'String',headlinetext);
set(gui.handles.infotext,'String',infotext);
plot(gui.handles.plotaxes,0,0,'ko');
set(gui.handles.plotaxes,'visible','off');

if ~yesno('This step will create modified files, are you sure?')
  %ues clicked no
  gui.page = gui.page-1;
  update;
  return
end

%find new folder name
filename = gui.files2load{1};
pos = find(filename==filesep);
if length(pos)<2
  myfailed('Failed.');
  return
end
basepathname = filename(1:(pos(end-1)-1));
foldername = filename((pos(end-1)+1):(pos(end)-1));
newfoldername = [foldername '-modified-' datestr(now,'HHMMSS')];

mymkdir([basepathname filesep newfoldername]);

%Sort on slicelocation
[~,ind] = sort(gui.slicelocations);

%find positon for min/max slice location and details for computing new
%imagepositions
startposition = gui.imagepositions{ind(1)};

imageorientation = gui.imageorientations{ind(1)};
zdir = cross(imageorientation(1:3),imageorientation(4:6));

%Loop over files
createmode = 'copy';
h = mywaitbarstart(length(gui.files2load),'Please wait',1,[]); %length(files2load));
for loop = 1:gui.nfiles
  
  if gui.included(gui.ind(loop))
    %This is the file we are working with
    newfilename = sprintf('out%06d.dcm',loop-1);
    outfilename = [basepathname filesep newfoldername filesep newfilename];
    
    %y = kx+m, x=loop-1
    m = gui.minincluded;
    k = (gui.maxincluded-gui.minincluded)/(gui.maxincludedind-gui.minincludedind);

    %Compute slice
    if gui.timeframes>1
      slice = floor((loop-1)/gui.timeframes);
    else
      slice = loop-1;
    end
  
    %Compute distance
    y = k*slice*gui.timeframes+m;
    
    %Compute new position
    newposition = startposition+y*zdir;
    newslicelocation = sum(newposition.*zdir);
    
    %Get file info
    dinfo = gui.dinfocell{gui.ind(loop)};
    
    %Read image
    im = dicomread(dinfo);
    
    %Assign new position
    if isfield(dinfo,'ImagePosition')
      dinfo.ImagePosition = newposition(:);
      dinfo.ImageOrientation = imageorientation; %ensure all are the same, take first
    else
      dinfo.ImagePositionPatient = newposition(:);
      dinfo.ImageOrientationPatient = imageorientation;
    end
    
    %Other tags we fix
    dinfo.SliceLocation = newslicelocation;
    dinfo.InstanceNumber = loop;
    dinfo.AcquisitionNumber = loop;
    dinfo.DerivationDescription = 'Sorted';
    
    %write
    try
      dicomwrite(im,outfilename,dinfo,'WritePrivate',true,'CreateMode',createmode);
    catch
      if isequal(createmode,'copy')
        createmode = 'create';
      end
      try
        dicomwrite(im,outfilename,dinfo,'WritePrivate',true,'CreateMode',createmode);
      catch
        myfailed('Failed to write file.')
        mywaitbarclose(h);
        return
      end

    end
    
    h = mywaitbarset(loop/gui.nfiles,h);
        
  end
  
end

mywaitbarclose(h);

mymsgbox('Done');

close_Callback;
openfile('refresh_Callback');

%---------------------
function next_Callback 
%---------------------
%Next page

global DATA

gui = DATA.GUI.SliceLocationTool;

gui.page = gui.page+1;
if gui.page>4
  gui.page = 4;
end
update;

%-------------------------
function previous_Callback 
%-------------------------
%Previous page

global DATA

gui = DATA.GUI.SliceLocationTool;

gui.page = gui.page-1;
if gui.page<1
  gui.page = 1;
end
update;

%------------------------------- 
function timeframesedit_Callback 
%-------------------------------
%Edit number of timeframes

global DATA

gui = DATA.GUI.SliceLocationTool;

s = mygetedit(gui.handles.timeframesedit);
v = str2double(s);
if isnan(v)
  v = gui.timeframes;
end

if ~isequal(rem(gui.nfiles,v),0)
  myfailed('Number of files and timeframes does not match.');
  v = 1;
end

%Assign
gui.timeframes = v;
set(gui.handles.timeframesedit,'String',sprintf('%d',gui.timeframes));

update;

%------------------------------
function insidehelper(newvalue)
%------------------------------
global DATA

gui = DATA.GUI.SliceLocationTool;

ax = gui.handles.plotaxes;
h = get(ax,'Children');
xlim = get(ax,'XLim');
ylim = get(ax,'YLim');

%Loop over all objects and check if inside
for loop = 1:length(h)
  ind = get(h(loop),'UserData');
  if ~isempty(ind)
    %Check if inside
    x = get(h(loop),'XData');
    y = get(h(loop),'YData');
    if (x>=xlim(1)) && (x<=xlim(2)) && (y>=ylim(1)) && (y<=ylim(2))
      gui.included(ind) = newvalue;
    end
  end
end

%------------------------
function exclude_Callback 
%------------------------
%Exclude visible points

insidehelper(false);
update;
%------------------------
function include_Callback 
%------------------------
%Exclude visible points

insidehelper(true);
update;

%----------------------
function close_Callback 
%----------------------
%Close the tool

global DATA

try
  gui = DATA.GUI.SliceLocationTool;
  delete(gui.fig)
catch
  delete(gfc);
  DATA.GUI.SliceLocationTool = [];
end

%---------------------
function scramblefiles 
%---------------------
%Function to scramble files. Used for debugging

files2load = openfile('getfiles2load');

%find new folder name
filename = files2load{1};
pos = find(filename==filesep);
if length(pos)<2
  myfailed('Failed.');
  return
end
foldername = filename(1:(pos(end)-1));

%Get files
f = dir([foldername filesep '*']);
f = f(~cat(1,f(:).isdir)); %remove . & ..

ind = randperm(length(f));

for loop = 1:length(f)
  movefile([foldername filesep f(loop).name],[foldername filesep sprintf('scrambled%06d.dcm',ind(loop))]);
end

