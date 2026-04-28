function varargout = croplvsaxstack(varargin)
%CROPLVSAXSTACK GUI to crop short-axis stack before performing 
%LV segmentation. Opens if crop is necessary for the current stack.
%#ok<*GVMIS> 
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%----------------------------------------------------
function [go,xall,yall] = init(no,cropped,xall,yall,isonlyboundingbox) 
%----------------------------------------------------
%Initate GUI for cropping LV stack.
%nos - nos to include in GUI
%cropped - boolean vector of whether stack is to be cropped or not
%xall - cell of x coordinates of suggested cropbox for each stack
%yall - cell of y coordinates of suggested cropbox for each stack
global DATA SET 

if nargin < 2
  cropped = true;
end
if nargin < 4 || isempty(xall) || isempty(yall)
  xall = round(SET(no).XSize/4):round(SET(no).XSize*3/4);
  yall = round(SET(no).YSize/4):round(SET(no).YSize*3/4);
end
if nargin < 5
  isonlyboundingbox = false;
end

gui = mygui('croplvsaxstack.fig');
DATA.GUI.CropLVSAXStack = gui;
gui.state = false;

set(gui.handles.saxaxes,'Visible','off');
gui.yall = yall;
gui.xall = xall;
gui.xallinit = xall;
gui.yallinit = yall;
gui.no = no;
gui.cropped = cropped;
if isonlyboundingbox
  % update figure name and button
  set(gui.fig,'Name',dprintf('Place bounding box'));
  set(gui.handles.defaultautocroppushbutton,'String',dprintf('Reset'));
  % add information text
  stri = dprintf('Make sure the LV is placed inside the yellow bounding box.');
  set(gui.handles.textinfo,'String',stri);
end

if cropped
  im = calcfunctions('remapuint8', ...
    SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),no);
  imagesc(im,'Parent',gui.handles.saxaxes);
  colormap(gui.handles.saxaxes,'gray');
  hold(gui.handles.saxaxes,'on');
  set(gui.handles.saxaxes,'Visible','on');
  axis(gui.handles.saxaxes,'image','off');
  y = yall;
  x = xall;
  gui.frameh = plot(gui.handles.saxaxes,[y(1) y(1) y(end) y(end) y(1)] , ...
    [x(1) x(end) x(end) x(1) x(1)],'y');
  gui.resizeh = plot(gui.handles.saxaxes,y(end),x(end),'yo');
  gui.translate = plot(gui.handles.saxaxes,y(1),x(1),'yo');
  set(gui.resizeh,'ButtonDownFcn',@(hObject,eventdata,handles)croplvsaxstack('resize'));
  set(gui.translate,'ButtonDownFcn',@(hObject,eventdata,handles)croplvsaxstack('move'));
end


if DATA.Testing
  gui.state = true;
else
  uiwait(gui.fig);
end
go = gui.state;
xall = gui.xall;
yall = gui.yall;
closegui;


%----------
function ok_Callback 
%----------
global DATA
gui = DATA.GUI.CropLVSAXStack;
gui.state = true;
uiresume(gui.fig);


%--------------
function cancel_Callback 
%--------------
global DATA

try
  gui = DATA.GUI.CropLVSAXStack;
  gui.state = false;
  uiresume(gui.fig);
catch
end


%----------------
function move
%----------------
global DATA
gui = DATA.GUI.CropLVSAXStack;
set(gui.fig,'WindowButtonMotionFcn',@(hObject,eventdata,handles)croplvsaxstack('motion','move'));
set(gui.fig,'WindowButtonUpFcn',@(hObject,eventdata,handles)croplvsaxstack('buttonup'));


%------------------
function resize
%------------------
global DATA
gui = DATA.GUI.CropLVSAXStack;
set(gui.fig,'WindowButtonMotionFcn',@(hObject,eventdata,handles)croplvsaxstack('motion','resize'));
set(gui.fig,'WindowButtonUpFcn',@(hObject,eventdata,handles)croplvsaxstack('buttonup'));


%----------------------
function motion(arg) 
%----------------------
global DATA
gui = DATA.GUI.CropLVSAXStack;
[y,x] = mygetcurrentpoint(gui.handles.saxaxes);
ax = axis(gui.handles.saxaxes);
ax = ax + 0.5*[1 -1 1 -1];

switch arg
    %NOT moving but resizing from the left upper corner; 
    %needs renaming/refactoring in the future
    case 'move' 
    xmax = gui.xall(end);
    ymax = gui.yall(end);
    y = round(max(1,min(ymax-2,y)));
    x = round(max(1,min(xmax-2,x)));
    set(gui.frameh,'XData',[y y ymax ymax y], ...
      'YData',[x xmax xmax x x]);
    set(gui.resizeh,'XData',ymax,'YData',xmax);
    set(gui.translate,'XData',y,'YData',x);

  case 'resize'
    xmin = gui.xall(1);
    ymin = gui.yall(1);   
    y = round(max(ymin+2,min(ax(2),y)));
    x = round(max(xmin+2,min(ax(4),x)));
    set(gui.frameh,'XData',[ymin ymin y y ymin], ...
      'YData',[xmin x x xmin xmin]);
    set(gui.resizeh,'XData',y,'YData',x);
    set(gui.translate,'XData',ymin,'YData',xmin);
end


%--------------------
function buttonup 
%--------------------
global DATA
gui = DATA.GUI.CropLVSAXStack;
x = get(gui.frameh,'YData');
y = get(gui.frameh,'XData');
gui.xall = x(1):x(3);
gui.yall = y(1):y(3);
set(gui.fig,'WindowButtonMotionFcn',[]);
set(gui.fig,'WindowButtonUpFcn',[]);


%--------------------
function default_Callback 
%--------------------
global DATA
gui = DATA.GUI.CropLVSAXStack;
gui.xall = gui.xallinit;
gui.yall = gui.yallinit;
if gui.cropped
  y = gui.yall;
  x = gui.xall;
  set(gui.frameh,'XData',[y(1) y(1) y(end) y(end) y(1)], ...
    'YData',[x(1) x(end) x(end) x(1) x(1)]);
  set(gui.resizeh,'XData',y(end),'YData',x(end));
  set(gui.translate,'XData',y(1),'YData',x(1));
end


%----------------
function closegui
%----------------
%Close GUI
global DATA

try
  DATA.GUI.CropLVSAXStack = close(DATA.GUI.CropLVSAXStack);
catch   %#ok<CTCH>
  DATA.GUI.CropLVSAXStack = [];
  delete(gcf);
end
