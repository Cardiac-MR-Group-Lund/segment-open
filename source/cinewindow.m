function cinewindow(varargin)
%Creates a cine window that loops until it is closed.

%Einar Heiberg

if nargin==0
  init;
else
  macro_helper(varargin{:});
  feval(varargin{:});
end;

%------------
function init
%------------
%Initialize a plotting window.

global DATA  

if isa(DATA.CineTimer,'timer')
  return;
else
    
  DATA.CineTimer = timer(...
    'Name','updatetimer',...
    'Period',0.05,...
    'StartDelay',0.05,...
    'executionmode','fixeddelay',...
    'TimerFcn',@update);

  update; %Call to initiate.

  %Start the timer
  start(DATA.CineTimer);

end;

%------------------------
function update(varargin)
%------------------------
global DATA SET NO

persistent s

if length(varargin)==1
  %Called with a kill. When the function is called from the timer the
  %length is 2. When called on setup the length is zero.
 
  try
    set(DATA.Handles.cinetoolicon,'state','off');
    stop(DATA.CineTimer);
  catch %#ok<CTCH>
  end;
  DATA.CineTimer = [];
  
  try
    delete(DATA.Handles.cineaxes);
  catch %#ok<CTCH>
  end;

  s = [];
  return;
end;

if isempty(s) || isempty(varargin)
  %First call => init
  s.no = NO;
  s.slice = SET(s.no).CurrentSlice;
  s.tf = 1;
  s.count = 1;
  
  %Create new handle
  if isfield(DATA.Handles,'cineaxes')
    try
      delete(DATA.Handles.cineaxes);
    catch
    end;
  end;
  
  DATA.Handles.cineaxes = axes('position',...
    [0.04 0.56 0.3 0.3],... %[0.05 0.68 0.3 0.3],...
    'parent',DATA.imagefig);

  %Create image
  s.im = calcfunctions('remapuint8', SET(s.no).IM,s.no,...
    calcfunctions('returnmapping',s.no),...
    SET(s.no).IntensityMapping.Contrast,...
    SET(s.no).IntensityMapping.Brightness);
  s.imagehandle = image(s.im(:,:,s.tf,s.slice),'parent',DATA.Handles.cineaxes);
  
  %Draw line
  x = SET(s.no).YSize;
  y = SET(s.no).XSize;  
  hold(DATA.Handles.cineaxes,'on');
  plot([1 x x 1 1],[1 1 y y 1],'y-');
  hold(DATA.Handles.cineaxes,'on');
  
  axis(DATA.Handles.cineaxes,'image','off');
  
  %When gray is current colormap, this causes all images to go dark
  if isempty(SET(s.no).Colormap)
    cinemap = gray(DATA.GUISettings.ColorMapSize);
  else
    cinemap = SET(s.no).Colormap;
  end
  colormap(DATA.Handles.cineaxes,cinemap);
  
  %Buttondown
  set(s.imagehandle,'ButtondownFcn','cinewindow(''update'',''kill'')');  
  return;
end;

%next update
s.count = s.count+1;
  
if SET(s.no).TSize>1
  %Get next frame
  f = 1.2;
  s.tf = 1+round((SET(s.no).TSize-1)*rem(now*3600*24*f,1));
  s.slice = SET(s.no).CurrentSlice;
else
  %Get next slice
  f = 0.25;
  s.slice = 1+round((SET(s.no).ZSize-1)*rem(now*3600*24*f,1));  
end;

if s.count>40000;
  update('kill');
  return;
end;

%Check if is handle
if ishandle(s.imagehandle)
  set(s.imagehandle,'cdata',s.im(:,:,s.tf,s.slice));
else
  update('kill');
  return;
end;
