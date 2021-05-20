function [varargout] = makeisotropic(varargin) 
%Make isotropic GUI

%Einar Heiberg

if nargin==0
  varargin = {'init_Callback'};
end

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%--------------------------
function init_Callback(arg) %#ok<DEFNU>
%--------------------------
%Init function

global DATA SET NO

gui = mygui('makeisotropic.fig');

DATA.GUI.MakeIsotropic = gui;

%Store operation
gui.no = NO;
gui.res = [SET(NO).ResolutionX SET(NO).ResolutionY (SET(NO).SliceThickness+SET(NO).SliceGap)];
gui.newres = median(gui.res);
gui.oldsize = [SET(NO).XSize SET(NO).YSize SET(NO).ZSize];
gui.init3dp = false; %if true init 3dp after makeisotropic operation
    
if nargin>0
  if isequal(arg,'init3DP')
    gui.init3dp = true;
  end
end

update;

%--------------
function update 
%--------------
global DATA

gui = DATA.GUI.MakeIsotropic;

gui.newsize = round(gui.oldsize.*gui.res/gui.newres); 
set(gui.handles.resolutiontext,'String',sprintf('%0.3g*%0.3g*%0.3g',gui.res(1),gui.res(2),gui.res(3)));
set(gui.handles.resolutionedit,'String',sprintf('%0.3g',gui.newres));
set(gui.handles.sizetext,'String',sprintf('%d*%d*%d',gui.oldsize(1),gui.oldsize(2),gui.oldsize(3)));
set(gui.handles.numvoxelstext,'String',sprintf('%d (100 %%)',prod(gui.oldsize)));
set(gui.handles.newsizetext,'String',sprintf('%d*%d*%d',gui.newsize(1),gui.newsize(2),gui.newsize(3)));
set(gui.handles.newnumvoxelstext,'String',sprintf('%d (%d %%)',prod(gui.newsize),round(100*prod(gui.newsize)/prod(gui.oldsize))));
gui.out = false;

%-------------------------------
function resolutionedit_Callback %#ok<DEFNU>
%-------------------------------
%edit box

global DATA

gui = DATA.GUI.MakeIsotropic;

s = get(gui.handles.resolutionedit,'String');
v = str2double(s);
if ~isnan(v)
  gui.newres = v;
  update;
else
  update;
end

%------------------------------
function makeisotropic_Callback %#ok<DEFNU>
%------------------------------
%Computes to isotropic, uses routines in tools

global DATA

gui = DATA.GUI.MakeIsotropic;

numel = prod(gui.newsize);
if numel>(512*512*768)
  if yesno('Very large volume size. Abort?')
    return;
  end
end

myworkon;

if abs(gui.res(3)/gui.newres-1)>1e-4
  %resample slices
  f = gui.res(3)/gui.newres;
  disp('Upsamling slices.');
  tools('upsampleslices_Callback',f);
end

if abs(gui.res(1)/gui.newres-1)>1e-4
  %resample inplane
  f = gui.res(1)/gui.newres;
  disp('Upsamling in plane.');
  tools('upsampleimage_Callback',f,gui.no);  
end

if gui.init3dp
  DATA.Handles.toggleiconholder.indent('ribbon3dp',1)
end

myworkoff;

close_Callback;

%----------------------
function close_Callback
%----------------------
%Close the gui

global DATA

try
  gui = DATA.GUI.MakeIsotropic;
  delete(gui.fig);
catch
end
DATA.GUI.MakeIsotropic = [];

