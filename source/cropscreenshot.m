function varargout = cropscreenshot(varargin)
% MATLAB code for cropscreenshot.fig
%#ok<*GVMIS> 

if nargin < 1 || isempty(varargin{1})
  varargin{1} = 'init';
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init 
%------------
%Initialize the GUI

global DATA  

gui = mygui('cropscreenshot.fig');
gui.SelectedGUI = 'Segment'; %default GUI to display in the crop region
DATA.GUI.CropScreenshot = gui;
gui.handles.screenshotimgh = image([],'Parent',gui.handles.screenshotaxes);
axis(gui.handles.screenshotaxes,'image','off');
gui.handles.screenshotaxes.Toolbar.Visible = 'off';
gui.handles.cropbox = line(gui.handles.screenshotaxes,nan, nan,'Color','y','ButtonDownFcn', ...
'cropscreenshot(''cropbox_Buttondown'')', 'LineWidth', 0.75);
gui.handles.framebox = line(gui.handles.screenshotaxes,nan, nan,'Color','k','ButtonDownFcn', ...
'cropscreenshot(''cropbox_Buttondown'')', 'LineWidth', 2);
gui.handles.resizeu = line(gui.handles.screenshotaxes,nan, nan,'Color','y','LineStyle','None','Marker','o','ButtonDownFcn', ...
'cropscreenshot(''cropbox_Buttondown'')','LineWidth', 0.75);
gui.handles.resizel = line(gui.handles.screenshotaxes,nan, nan,'Color','y','LineStyle','None','Marker','o','ButtonDownFcn', ...
'cropscreenshot(''cropbox_Buttondown'')','LineWidth', 0.75);

selectedbtn = gui.handles.selectcropregioncuibuttongroup.SelectedObject;
switch selectedbtn.Tag
  case 'imagepanelradiobutton'
    imagepanelradiobutton_Callback
  case 'userselectionradiobutton'
    setimage;  
    resetcropboundaries;
    updatecropbox;
end

%-------------------------------------
function imagepanelradiobutton_Callback 
%---------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;
setimage;
gui.handles.resetpushbutton.Enable = 'off';

%-------------------------------------
function userselectionradiobutton_Callback  
%---------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;
try
  gui.handles.resetpushbutton.Enable = 'on';
  setimage;
  updatecropbox;
catch me
  mydispexception(me)
end

%-------------------------------------
function resetcropboundaries 
%-------------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;
[sizeX, sizeY,~] = size(gui.handles.screenshotimgh.CData);
gui.startx  = 1;
gui.starty = 1;
gui.endx = sizeY;
gui.endy = sizeX;
%-------------------------------------
function resetframe
%-------------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;
[sizeX, sizeY,~] = size(gui.handles.screenshotimgh.CData);
selectedbtn = gui.handles.selectcropregioncuibuttongroup.SelectedObject;
switch selectedbtn.Tag
  case 'imagepanelradiobutton'
    offset = 2;
  case 'userselectionradiobutton'
    offset = 5;
end

startx  = 1-offset;
starty = 1-offset;
endx = sizeY+offset;
endy = sizeX+offset;
boxX = [startx,endx,endx,startx,startx];
boxY = [starty,starty,endy,endy,starty];

gui.handles.framebox.XData = boxX;
gui.handles.framebox.YData = boxY;

%-------------------------------------
function cropbox_Buttondown
%---------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;

[y,x] = mygetcurrentpoint(gui.handles.screenshotaxes);

%distance to startx/y
dstart = norm([gui.startx,gui.starty]-[y,x]);
dstop = norm([gui.endx,gui.endy]-[y,x]);

if dstop <= 25 && dstop < dstart
  %scale box lower point
  pointtoscale = 'scalel'; 
elseif dstart <= 25
  %scale box upper point
  pointtoscale = 'scaleu';
else
  return  
end

gui.fig.WindowButtonMotionFcn = sprintf('cropscreenshot(''cropbox_Motion'',''%s'')',pointtoscale);
gui.fig.WindowButtonUpFcn = 'cropscreenshot(''cropbox_Buttonup'')';

%-------------------------------------
function cropbox_Motion(type)
%---------------------------------
%motion function to move cropping box points.
global DATA 
gui = DATA.GUI.CropScreenshot;
[x,y] = mygetcurrentpoint(gui.handles.screenshotaxes);
[sizeX, sizeY,~] = size(gui.handles.screenshotimgh.CData);
switch type
  case 'scaleu'
    startx = min(sizeY,max(1,x));
    starty = min(sizeX,max(1,y));
    endx = gui.endx;
    endy = gui.endy;
  case 'scalel'
    startx = gui.startx;
    starty = gui.starty;
    endx = min(sizeY,max(1,x));
    endy = min(sizeX,max(1,y));
end

%update params
gui.startx = startx;
gui.starty = starty;
gui.endx = endx;
gui.endy = endy;

updatecropbox

%-------------------------------------
function cropbox_Buttonup
%---------------------------------
%motion function for mid slice axes. Move box points.
global DATA 
gui = DATA.GUI.CropScreenshot;
gui.fig.WindowButtonMotionFcn = '';
gui.fig.WindowButtonUpFcn = '';

%-------------------------------------
function updatecropbox
%---------------------------------
%updates cropbox
global DATA 
gui = DATA.GUI.CropScreenshot;
boxX = [gui.startx,gui.endx,gui.endx,gui.startx,gui.startx];
boxY = [gui.starty,gui.starty,gui.endy,gui.endy,gui.starty];

gui.handles.cropbox.XData = boxX;
gui.handles.cropbox.YData = boxY;

gui.handles.resizeu.XData = gui.startx;
gui.handles.resizeu.YData = gui.starty;
gui.handles.resizel.XData = gui.endx;
gui.handles.resizel.YData = gui.endy;

%-------------------------------------
function reset_Callback
%---------------------------------
resetcropboundaries;
userselectionradiobutton_Callback;

%-------------------------------------
function close_Callback 
%---------------------------------
global DATA

try
  DATA.GUI.CropScreenshot = close(DATA.GUI.CropScreenshot);
catch me
  mydispexception(me)%#ok<CTCH>
  DATA.GUI.CropScreenshot =[];
  delete(gcbf);
end

%-------------------------------------
function save_Callback
%---------------------------------
[fname, extindex] = setfilename;
im = getimage;
successful =  writeimage(im,fname,extindex);
if successful
  close_Callback;
end

%-------------------------------------
function [fname, extindex] = setfilename
%-------------------------------------
global DATA SET
gui = DATA.GUI.CropScreenshot;
selectedbtn = gui.handles.savefileuibuttongroup.SelectedObject;
switch selectedbtn.Tag
  case 'fileradiobutton'
    [filename, pathname,extindex] = myuiputfile(...
      { '*.png','PNG image (*.png)';...
      '*.jpg','JPEG image (*.jpg)';...
      '*.bmp','BMP image (*.bmp)';...
      '*.tif','TIFF image (*.tif)'},...
      dprintf('Save file as'),'screenshot');
  case 'pacsradiobutton'
    defname = 'screenshot';
    filename = char(myinputdlg({'Enter File Name'},'',1,{defname}));
    isnotok = 1;
    while isnotok
      if isempty(filename)
        return;
      elseif ~isvarname(filename)
        myfailed(dprintf('File name %s is not valid.\nFile name can contain letters, digits, and underscores',filename),DATA.GUI.Segment);
        isnotok = 1;
        filename = char(myinputdlg({'Enter File Name'},'',1,{defname}));
      elseif isvarname(filename)
        break;
      end
    end
    pathname = getpreferencespath;
    extindex = 5;
  case 'reportradiobutton'
    name = removeforbiddenchars(SET(1).PatientInfo.Name);
    if isempty(name)
      name = 'Hidden';
    end
    pathname = fullfile(DATA.Pref.Pacs.ReportsheetPath, ...
      name,'Screenshots');
    if ~exist(pathname,'dir')
      sucess = mkdir(pathname);
      if ~sucess
        myfailed('Could not create screenshot directory. Aborted');
        return
      end
    end
    reportdir = dir(fullfile(pathname,'screenshot*.png'));
    nbrs = zeros(1,numel(reportdir));
    for i = 1:numel(reportdir)
      nbrs(i) = sscanf(reportdir(i).name(12:15),'%f');
    end
    newnbr = min(setdiff(1:9999,nbrs));
    filename = sprintf('screenshot_%04.0f.png',newnbr);
    extindex = 1;
    if ~isfield(DATA.Pref.Pacs, 'ScreenshotPath')
      DATA.Pref.Pacs.ScreenshotPath = pathname;
    elseif ~strcmp(DATA.Pref.Pacs.ScreenshotPath, pathname)
      DATA.Pref.Pacs.ScreenshotPath = pathname;
    end
  case 0
    filename = 0;
end

fname = fullfile(pathname,filename);
%Add extension if necessary
[~,~,ext] = fileparts(fname);
if isempty(ext)
  switch extindex
    case 1
      fname = [fname '.png'];
    case 2
      fname = [fname '.jpg'];
    case 3
      fname = [fname '.bmp'];
    case 4
      fname = [fname '.tif'];
    case 5
      fname = [fname '.dcm'];
  end
end

%-------------------------------------
function successful = writeimage(im,fname,extindex)
%-------------------------------------
global DATA
successful = false;
saveoptions = {
  'png',...
  'jpg',...
  'bmp',...
  'tif',...
  'pacs'
  };
if isnumeric(extindex)
  formatstr = saveoptions{extindex};
else
  formatstr = extindex;
end
switch formatstr
  case 'png'
    try
      imwrite(im,fname,'png','bitdepth',8,'software','Segment',...
        'creationtime',datestr(now));
      disp('Export successful');
      successful = true;
    catch me
      mydispexception(me);
      disp('Export failed');
    end
  case {'jpeg','jpg'}
    try
      imwrite(im,fname,'jpg','quality',100);
      disp('Export successful');
      successful = true;
    catch me
      mydispexception(me);
      disp('Export failed');
    end
  case 'bmp'
    try
      imwrite(im,fname,'bmp');
      successful = true;
      disp('Export successful');
    catch me
      mydispexception(me)
      disp('Export failed');
    end
  case {'tif','tiff'}
    try
      imwrite(im,fname,'tif');
      successful = true;
      disp('Export successful');
    catch me
      mydispexception(me);
      disp('Export failed');
    end
    
  case 'pacs'
    makeimagedicom(im,fname,DATA.ViewPanels(DATA.CurrentPanel));
    try
      successful = pacsaccess('savetopacs_helper',{fname},getpreferencespath);
    catch me
      successful = 0;
      mydispexception(me);
      myfailed(me.message);
    end
    try
      status = 1;
      delete(fname)
    catch me
      mydispexception(me)
      status = 0;
    end
    
    if successful
      stri = dprintf('Save to PACS success');
    else
      stri = dprintf('Save to PACS failed');
    end
    
    if status == 1
      stri = [stri,dprintf(', succesfully removed report from computer.')];
    else
      stri = [stri,dprintf(', failed to remove report from computer.')];
    end
    mymsgbox(stri)
end

%---------------------------------
function img = getimage
%---------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;
img = [];
selectedbtn = gui.handles.selectcropregioncuibuttongroup.SelectedObject;
switch selectedbtn.Tag
  case 'imagepanelradiobutton'
    try
      img = gui.handles.screenshotimgh.CData;
    catch me  
      mydispexception(me);
    end
  case 'userselectionradiobutton'
    [sizeX, sizeY,~] = size(gui.handles.screenshotimgh.CData);
    gui.startx  = min(sizeY,max(1,round(gui.startx)));
    gui.starty = min(sizeX,max(1,round(gui.starty)));
    gui.endx = min(sizeY,max(1,round(gui.endx)));
    gui.endy = min(sizeX,max(1,round(gui.endy)));
    try
      img = gui.handles.screenshotimgh.CData(gui.starty:gui.endy,gui.startx:gui.endx,:);
    catch me  
      mydispexception(me);
    end
end

%---------------------------------
function selectgui_Callback
%---------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;

%find all opened GUIs
indices = ~structfun(@isempty, DATA.GUI);
guis = fieldnames(DATA.GUI);
openedguis = guis(indices);
toremove = strcmp(openedguis,'CropScreenshot'); %remove screeshot GUI
openedguis(toremove) = [];

%prompt user to choose a GUI
n = 1;
f(n).Field = 'selectedgui';
f(n).Label = dprintf('Select a GUI for the screenshot:');
f(n).Default = openedguis;
[outs,ok] = myinputstruct(f,dprintf('Select GUI'),15);
if ok
  ind = outs.selectedgui;
  g = openedguis(ind);
else
  % cancel was clicked
  g = {'Segment'}; %main GUI is default
end

%assign selected gui
gui.SelectedGUI = g{1};
setimage;
resetcropboundaries;
updatecropbox;

%---------------------------------
function setimage
%---------------------------------
global DATA
gui = DATA.GUI.CropScreenshot;
selectedbtn = gui.handles.selectcropregioncuibuttongroup.SelectedObject;
try
  switch selectedbtn.Tag
    case 'imagepanelradiobutton'
      h = DATA.Handles.boxaxes;
      imgframe = mygetframe(h);
      boxvis = 'off';    
    case 'userselectionradiobutton'
      % get whole figure of segment as a frame 
      selectedGUI = gui.SelectedGUI;
      sgui = DATA.GUI.(selectedGUI);
      imgframe = mygetframe(sgui.fig);
      imgframe.colormap = colormap(sgui.fig); 
      boxvis = 'on';
  end
  img = frame2im(imgframe);
  gui.handles.screenshotimgh.CData = img;
  set([gui.handles.cropbox,gui.handles.resizeu,gui.handles.resizel],'Visible',boxvis);
  resetframe;
catch me
  mydispexception(me);
end