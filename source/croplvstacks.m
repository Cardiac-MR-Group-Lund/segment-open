function varargout = croplvstacks(varargin)
%CROPLVSTACKS GUI to crop short-axis and long-axis stacks before performing 
%LV segmentation. Opens if crop is necessary for the current stacks.
macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%----------------------------------------------------
function [go,xall,yall] = init(nos,cropped,xall,yall) %#ok<DEFNU>
%----------------------------------------------------
%Initate GUI for cropping LV stacks.
%nos - nos to include in GUI
%cropped - boolean vector of whether stack is to be cropped or not
%xall - cell of x coordinates of suggested cropbox for each stack
%yall - cell of y coordinates of suggested cropbox for each stack
global DATA SET

numnos = numel(nos);
if nargin < 4
  xall = cell(1,numnos);
  yall = xall;
  for i = 1:numnos
    xall{i} = round(SET(nos(i)).XSize/3):round(SET(nos(i)).XSize*2/3);
    yall{i} = round(SET(nos(i)).YSize/3):round(SET(nos(i)).YSize*2/3);
  end
  if nargin < 2
    cropped = true(1,numnos);
  end
end

gui = mygui('croplvstacks.fig');
DATA.GUI.CropLVStacks = gui;
gui.state = false;

gui.axh = [gui.handles.saxaxes ...
  gui.handles.ch2axes ...
  gui.handles.ch3axes ...
  gui.handles.ch4axes];
set(gui.axh,'Visible','off');
gui.yall = yall;
gui.xall = xall;
gui.xallinit = xall;
gui.yallinit = yall;
gui.nos = nos;
gui.cropped = cropped;

for i = 1:numnos
  no = nos(i);
  if cropped(i)
    im = calcfunctions('remapuint8', ...
      SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice),no);
    imagesc(im,'Parent',gui.axh(i));
    colormap(gui.axh(i),'gray');
    hold(gui.axh(i),'on');
    set(gui.axh(i),'Visible','on');
    axis(gui.axh(i),'image','off');
    y = yall{i};
    x = xall{i};
    gui.frameh(i) = plot(gui.axh(i),[y(1) y(1) y(end) y(end) y(1)] , ...
      [x(1) x(end) x(end) x(1) x(1)],'y');
    gui.resizeh(i) = plot(gui.axh(i),y(end),x(end),'yo');
    gui.translate(i) = plot(gui.axh(i),y(1),x(1),'yo');
    set(gui.resizeh(i),'ButtonDownFcn',@(hObject,eventdata,handles)croplvstacks('resize',i));
    set(gui.translate(i),'ButtonDownFcn',@(hObject,eventdata,handles)croplvstacks('move',i));
  end
end

uiwait(gui.fig);
go = gui.state;
xall = gui.xall;
yall = gui.yall;
closegui;

%----------
function ok_Callback %#ok<DEFNU>
%----------
global DATA
gui = DATA.GUI.CropLVStacks;
gui.state = true;
uiresume(gui.fig);

%--------------
function cancel_Callback %#ok<DEFNU>
%--------------
global DATA

try
  gui = DATA.GUI.CropLVStacks;
  gui.state = false;
  uiresume(gui.fig);
catch
end

%----------------
function move(ix) %#ok<DEFNU>
%----------------
global DATA
gui = DATA.GUI.CropLVStacks;
set(gui.fig,'WindowButtonMotionFcn',@(hObject,eventdata,handles)croplvstacks('motion','move',ix));
set(gui.fig,'WindowButtonUpFcn',@(hObject,eventdata,handles)croplvstacks('buttonup',ix));

%------------------
function resize(ix) %#ok<DEFNU>
%------------------
global DATA
gui = DATA.GUI.CropLVStacks;
set(gui.fig,'WindowButtonMotionFcn',@(hObject,eventdata,handles)croplvstacks('motion','resize',ix));
set(gui.fig,'WindowButtonUpFcn',@(hObject,eventdata,handles)croplvstacks('buttonup',ix));

%----------------------
function motion(arg,ix) %#ok<DEFNU>
%----------------------
global DATA
gui = DATA.GUI.CropLVStacks;
[y,x] = mygetcurrentpoint(gui.axh(ix));
ax = axis(gui.axh(ix));
ax = ax + 0.5*[1 -1 1 -1];
% needed for the old proper move
% xmin = gui.xall{ix}(1);
% ymin = gui.yall{ix}(1);
% xdiff = gui.xall{ix}(end)-xmin;
% ydiff = gui.yall{ix}(end)-ymin;
switch arg
    %NOT moving but resizing from the left upper corner; 
    %needs renaming/refactoring in the future
    case 'move' 
    xmax = gui.xall{ix}(end);
    ymax = gui.yall{ix}(end);
    y = round(max(1,min(ymax-2,y)));
    x = round(max(1,min(xmax-2,x)));
    set(gui.frameh(ix),'XData',[y y ymax ymax y], ...
      'YData',[x xmax xmax x x]);
    set(gui.resizeh(ix),'XData',ymax,'YData',xmax);
    set(gui.translate(ix),'XData',y,'YData',x);
        % old case with proper move
%   case 'move' 
%     y = round(min(max(y,ax(3)),ax(4)-ydiff));
%     x = round(min(max(x,ax(1)),ax(2)-xdiff));
%     set(gui.frameh(ix),'XData',[y y y+ydiff y+ydiff y], ...
%       'YData',[x x+xdiff x+xdiff x x]);
%     set(gui.resizeh(ix),'XData',y+ydiff,'YData',x+xdiff);
%     set(gui.translate(ix),'XData',y,'YData',x);
  case 'resize'
    xmin = gui.xall{ix}(1);
    ymin = gui.yall{ix}(1);   
    y = round(max(ymin+2,min(ax(2),y)));
    x = round(max(xmin+2,min(ax(4),x)));
    set(gui.frameh(ix),'XData',[ymin ymin y y ymin], ...
      'YData',[xmin x x xmin xmin]);
    set(gui.resizeh(ix),'XData',y,'YData',x);
    set(gui.translate(ix),'XData',ymin,'YData',xmin);
end
%--------------------
function buttonup(ix) %#ok<DEFNU>
%--------------------
global DATA
gui = DATA.GUI.CropLVStacks;
x = get(gui.frameh(ix),'YData');
y = get(gui.frameh(ix),'XData');
gui.xall{ix} = x(1):x(3);
gui.yall{ix} = y(1):y(3);
set(gui.fig,'WindowButtonMotionFcn',[]);
set(gui.fig,'WindowButtonUpFcn',[]);


%--------------------
function default_Callback %#ok<DEFNU>
%--------------------
global DATA
gui = DATA.GUI.CropLVStacks;
gui.xall = gui.xallinit;
gui.yall = gui.yallinit;
for i = 1:numel(gui.nos)
  if gui.cropped(i)
    y = gui.yall{i};
    x = gui.xall{i};
    set(gui.frameh(i),'XData',[y(1) y(1) y(end) y(end) y(1)], ...
      'YData',[x(1) x(end) x(end) x(1) x(1)]);
    set(gui.resizeh(i),'XData',y(end),'YData',x(end));
    set(gui.translate(i),'XData',y(1),'YData',x(1));
  end
end


%----------------
function closegui
%----------------
%Close GUI
global DATA

try
  DATA.GUI.CropLVStacks = close(DATA.GUI.CropLVStacks);
catch   %#ok<CTCH>
  DATA.GUI.CropLVStacks = [];
  delete(gcf);
end
