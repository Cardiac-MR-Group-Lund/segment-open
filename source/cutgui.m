function [varargout] = cutgui(varargin)
%Gui for square cuts of x mm size.

if nargin==0
  varargin = {'init'};
end

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%-------------
function init
%-------------
%initialize the cut GUI

global DATA SET NO
DATA.GUI.CutGui = mygui('cutgui.fig');
gui = DATA.GUI.CutGui;
gui.no = NO;

%this is the box handle which is shown in the image.
gui.handles.box=[];
gui.boxplaced = 0;

%Start by showing the middle image in the preview axes.
gui.slice = SET(gui.no).CurrentSlice; %1+SET(gui.no).ZSize-max(1,SET(gui.no).CurrentSlice);
%gui.slice = round(SET(gui.no).ZSize/2);
cmap = gray(256);
im = SET(gui.no).IM(:,:,1,gui.slice);
c = SET(gui.no).IntensityMapping.Contrast;
b = SET(gui.no).IntensityMapping.Brightness;
rim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,1),c,b);%segment('remap',im,cmap(:,1),c,b);
gim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,2),c,b);%segment('remap',im,cmap(:,2),c,b);
bim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,3),c,b);%segment('remap',im,cmap(:,3),c,b);
im = cat(3,rim, gim, bim);
gui.handles.imagehandle = image(im,'parent',gui.handles.previewaxes);
set(gui.handles.previewaxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [])
%axis(gui.handles.previewaxes,'image','off');

%Set the slider so that it matches the middle slice
if SET(gui.no).ZSize==1
  set(gui.handles.cutslider,'visible','off')
else
  set(gui.handles.cutslider,'Max',SET(gui.no).ZSize)
  set(gui.handles.cutslider,'Value',1+SET(gui.no).ZSize-gui.slice)
  set(gui.handles.cutslider,'Min',1)  
  set(gui.handles.cutslider,'Sliderstep',[1/(SET(gui.no).ZSize-1),0.1])
end

%Set the 150 radio button as indented by default
setsize(150)

set(gui.handles.imagehandle,'buttondownfcn','cutgui(''placebox'')')
set(gui.fig,'windowbuttonmotionfcn','cutgui(''motion'')')


%-------------------
function updateimage
%-------------------
%update image display

global DATA SET
gui = DATA.GUI.CutGui;
cmap = gray(256);
%always look in first timeframe for now
im = SET(gui.no).IM(:,:,1,gui.slice);
c = SET(gui.no).IntensityMapping.Contrast;
b = SET(gui.no).IntensityMapping.Brightness;
rim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,1),c,b);%segment('remap',im,cmap(:,1),c,b);
gim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,2),c,b);%segment('remap',im,cmap(:,2),c,b);
bim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,3),c,b);%segment('remap',im,cmap(:,3),c,b);
im = cat(3,rim, gim, bim);
set(gui.handles.imagehandle,'cdata',im,'parent',gui.handles.previewaxes);


%-------------------
function setsize(x)
%------------------
%set sixe of cop box

global DATA SET
gui = DATA.GUI.CutGui;

if nargin == 0
  x = get(gui.handles.sizeedit,'string');
  x = str2double(x);
  if isnan(x) %str2double returns nan if any text is included in the string. 
    return
  end
end

switch x
  case 10
    set(gui.handles.tenmmradiobutton,'Value',1)
    set(gui.handles.onehundredmmradiobutton,'Value',0)
    set(gui.handles.onefiftymmradiobutton,'Value',0)
    set(gui.handles.twofiftymmradiobutton,'Value',0)
    set(gui.handles.sizeedit,'string','10');

  case 100
    set(gui.handles.tenmmradiobutton,'Value',0)
    set(gui.handles.onehundredmmradiobutton,'Value',1)
    set(gui.handles.onefiftymmradiobutton,'Value',0)
    set(gui.handles.twofiftymmradiobutton,'Value',0)
    set(gui.handles.sizeedit,'string','100');

  case 150
    set(gui.handles.tenmmradiobutton,'Value',0)
    set(gui.handles.onehundredmmradiobutton,'Value',0)
    set(gui.handles.onefiftymmradiobutton,'Value',1)
    set(gui.handles.twofiftymmradiobutton,'Value',0)
    set(gui.handles.sizeedit,'string','150');

  case 250
    set(gui.handles.tenmmradiobutton,'Value',0)
    set(gui.handles.onehundredmmradiobutton,'Value',0)
    set(gui.handles.onefiftymmradiobutton,'Value',0)
    set(gui.handles.twofiftymmradiobutton,'Value',1)
    set(gui.handles.sizeedit,'string','250');

  otherwise
    set(gui.handles.tenmmradiobutton,'Value',0)
    set(gui.handles.onehundredmmradiobutton,'Value',0)
    set(gui.handles.onefiftymmradiobutton,'Value',0)
    set(gui.handles.twofiftymmradiobutton,'Value',0)
    
end

%store the side length so it is accesible for draw functions outside this
%one
gui.x = x;

%check if boxplaced then redraw placed box
if gui.boxplaced
  %get the pixel radi in the image
  xrad = gui.x/2/SET(gui.no).ResolutionX;
  yrad = gui.x/2/SET(gui.no).ResolutionY;
  set(gui.handles.box,'ydata',[gui.xc-xrad,gui.xc+xrad,gui.xc+xrad,gui.xc-xrad,gui.xc-xrad],...
    'xdata',[gui.yc-yrad,gui.yc-yrad,gui.yc+yrad,gui.yc+yrad,gui.yc-yrad],'visible','on');
end
%-----------------------
function slider_callback %#ok<DEFNU>
%-----------------------
%callback from slice slider to change image slice in display

global DATA SET
gui = DATA.GUI.CutGui;
gui.slice = 1+SET(gui.no).ZSize-max(1,round(get(gui.handles.cutslider,'Value')));

updateimage

%-----------------------
function close_callback
%----------------------
%close the crop GUI

global DATA
close(DATA.GUI.CutGui);

%---------------
function motion
%---------------
%move the crop box
global DATA
gui = DATA.GUI.CutGui;

handleAddress=hittest(gui.fig);
if isequal(handleAddress,gui.handles.imagehandle)
  drawbox;
elseif ~isempty(gui.handles.box)
  try 
  delete(gui.handles.box)
  catch
  end
  gui.box = [];
end

%-------------------------
function [xc,yc] = drawbox
%-------------------------
%draw the crop box in image display
global DATA SET
gui = DATA.GUI.CutGui;
box = gui.handles.box;
ima = gui.handles.previewaxes;

%get the pixel radi in the image
xrad = gui.x/2/SET(gui.no).ResolutionX;
yrad = gui.x/2/SET(gui.no).ResolutionY;

%Draw box around cursor
if isempty(box) ||~ishandle(box) 
  hold(ima,'on')
  [yc,xc] = mygetcurrentpoint(ima);
  gui.handles.box = plot(ima,...
    [yc-yrad,yc-yrad,yc+yrad,yc+yrad,yc-yrad],...
    [xc-xrad,xc+xrad,xc+xrad,xc-xrad,xc-xrad],'y');
  hold(ima,'off')
else
  [yc,xc] = mygetcurrentpoint(ima);
  set(box,'ydata',[xc-xrad,xc+xrad,xc+xrad,xc-xrad,xc-xrad],...
    'xdata',[yc-yrad,yc-yrad,yc+yrad,yc+yrad,yc-yrad],'visible','on');
end

%-----------------
function placebox
%----------------
%first box is placed then if apply is clicked cut is made

global DATA
gui = DATA.GUI.CutGui;

gui.boxplaced = 1;
ima = gui.handles.previewaxes;
[xc,yc] = drawbox;
gui.xc =round(xc); 
gui.yc =round(yc);

set(gui.handles.imagehandle,'buttondownfcn','cutgui(''unlockbox'')')
set(gui.fig,'windowbuttonmotionfcn','')

%------------------
function unlockbox
%------------------
%unlock the crop box

global DATA
gui = DATA.GUI.CutGui;
motion
gui.boxplaced = 0;
set(gui.handles.imagehandle,'buttondownfcn','cutgui(''placebox'')')
set(gui.fig,'windowbuttonmotionfcn','cutgui(''motion'')')

%--------------------
function cut_callback %#ok<DEFNU>
%--------------------
%perform the crop 

global DATA SET
gui = DATA.GUI.CutGui;

if ~gui.boxplaced
  mywarning('Box must be placed before applying crop. Click in the preview where you want to crop to place the box.')
  return
end

xc = gui.xc;
yc = gui.yc;
xrad = round(gui.x/2/SET(gui.no).ResolutionX);
yrad = round(gui.x/2/SET(gui.no).ResolutionY);

%make the cut
xl = max(1,xc-xrad);
yl = max(1,yc-yrad);
xu = min(SET(gui.no).XSize,xc+xrad);
yu = min(SET(gui.no).YSize,yc+yrad);

%This switch makes the cut correct in the internal coordinate system.
xind = xl:xu;
yind = yl:yu;

%update interface
nos = tools('crophelper',gui.no,xind,yind);
tools('cropupdate',nos);

%update values in SET struct

%close the image
close_callback

