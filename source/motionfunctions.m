function varargout = motionfunctions(varargin)
% Functions for cursor motion
% Klas
%Invoke subfunction

%#ok<*GVMIS>
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%---------------------------
function select_motion(panel)
%---------------------------
%Function that highlights the selected panel if montage also highlight the frame
global DATA SET

no = DATA.ViewPanels(panel);
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

slice = viewfunctions('clickedslice',panel,x,y);

%clicked in black area of montage
if isempty(slice)
  return
end

if slice == SET(no).CurrentSlice
  SET(no).StartSlice = slice;
  SET(no).EndSlice = slice;
elseif slice > SET(no).CurrentSlice
  SET(no).EndSlice = slice;
elseif slice < SET(no).CurrentSlice
  SET(no).StartSlice = slice;
end

drawfunctions('drawselectedslice',panel)

%----------------------------------------------
function [contrastout,brightnessout] = contrast_motion(panel,slice,xstart,ystart,xsize,ysize)
%----------------------------------------------
%motion function for contrast image
global DATA SET

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

%we need to update the delta contrast and delta brightness based on dx and
%dy
deltacontrast = (x-xstart)/xsize;
deltabrightness = -(y-ystart)/ysize;

no = DATA.ViewPanels(panel);

if ~isempty(SET(no).IntensityScaling) && ~isempty(SET(no).IntensityMapping)
  [window,level] = calcfunctions('con2win',...
    SET(no).IntensityMapping.Contrast+deltacontrast,...
    SET(no).IntensityMapping.Brightness+deltabrightness,no);
  [contrast,brightness] = calcfunctions('win2con',window,level,no);
else
  contrast = SET(no).IntensityMapping.Contrast+deltacontrast;
  brightness = SET(no).IntensityMapping.Brightness+deltabrightness;
end

if any(strcmp(DATA.ViewPanelsType{panel},{'montage','montagerow','montagesegmented'}))
  if strcmp(DATA.ViewPanelsType{panel},'montagesegmented')
    segmentedonly=true;
  else
    segmentedonly=false;
  end
  oneextraslice = true;
  slicestoinclude = [];
  usezoomstate = true;
  im = calcfunctions('calcmontageviewim', ...
    no,fliplr(DATA.ViewPanelsMatrix{panel}),segmentedonly,calcfunctions('returnmapping',no),contrast,brightness, ...
    [],SET(no).CurrentTimeFrame,oneextraslice,slicestoinclude,usezoomstate);
  set(DATA.Handles.imagehandles(panel),'cdata',squeeze(im));
else
  %update the image with the new contrast and brightness settings
  im = calcfunctions('remapuint8',SET(no).IM(:,:,SET(no).CurrentTimeFrame,slice),no,...
    calcfunctions('returnmapping',no),contrast,brightness);
  scale = viewfunctions('getscale',panel);
  im = imresize(im,scale);
  if ~isempty(SET(no).Colormap)|| ndims(im)==3
    set(DATA.Handles.imagehandles(panel),'cdata',im);
  else
    set(DATA.Handles.imagehandles(panel),'cdata',cat(3,im,im,im));
  end


  %update colorbar
  % if DATA.GUISettings.ShowColorbar
  %   h = colorbar(DATA.Handles.imageaxes(panel),'EastOutside');
  %     if isempty(SET(no).Colormap)
  %         colormap(h,'gray');
  %     else
  %         colormap(h,SET(no).Colormap);
  %     end
  %     h.Color = DATA.GUISettings.ForegroundColor;
  %     [win level]=calcfunctions('con2win',SET(no).IntensityMapping.Contrast, SET(no).IntensityMapping.Brightness)
  %     set(DATA.Handles.imageaxes(panel),'CLim',[(level-win/2) (level+win/2)])
  % end
  %

end

%update colorbar
if DATA.GUISettings.ShowColorbar
  createfunctions('addcolorbar',panel,contrast, brightness)
  %   h = colorbar(DATA.Handles.imageaxes(panel),'East');
  %     [win, level]=calcfunctions('con2win',contrast, brightness);
  %     set(DATA.Handles.imageaxes(panel),'CLim',[(level-win/2) (level+win/2)])
  %     if isempty(SET(no).Colormap)
  %         colormap(h,'gray');
  %     else
  %         colormap(h,SET(no).Colormap);
  %     end
  %     h.Color = DATA.GUISettings.ForegroundColor;

end
%
% set(DATA.Handles.imagehandles(panel),'cdata',cat(3,im,im,im));

if nargout == 2
  brightnessout = brightness;
  contrastout = contrast;
end

%----------------------------------------------
function [Xout, Yout] = crop_motion(panel,xstart,ystart)
%----------------------------------------------
%Motion function for cropping.

global DATA

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

%This always makes a box
X = [xstart,xstart,x,x,xstart];
Y = [ystart,y,y,ystart,ystart];

DATA.CursorX = X;
DATA.CursorY = Y;

DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;

if nargout == 2
  Xout = X;
  Yout = Y;
end

%----------------------------------------------
function measure_motion(panel,pointind)
%----------------------------------------------
%Motion function for measures. Has panel and point index of the currently dragged
%measurement as input.

global DATA SET NO

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
DATA.CursorX(pointind) = x;
DATA.CursorY(pointind) = y;
DATA.CursorZ(pointind) = SET(NO).CurrentSlice;

DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;


%----------------------------------------------
function measuretext_motion(txthandle,panel,measureind)
%----------------------------------------------
% Motion function for measures' text
global DATA SET

[xpos,ypos] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

set(txthandle,'Position', [xpos ypos 0])
textcounter = get(txthandle,'UserData');
scale = viewfunctions('getscale',panel);
no = DATA.ViewPanels(panel);
if contains(DATA.ViewPanelsType{panel},'montage')
  [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},SET(no).Measure(measureind).Z(end));
  imdim = zoomfunctions.getxysize(no,panel);
  xt = (xl-1)*imdim.XSize - imdim.XStart +1;
  yt = (yl-1)*imdim.YSize - imdim.YStart +1;
else
  xt = 0;
  yt = 0;
end
offsety = xpos/scale -yt - SET(no).Measure(measureind).Y(end);
offsetx = ypos/scale -xt - SET(no).Measure(measureind).X(end);
SET(no).Measure(measureind).Offset = [offsety offsetx];
set(DATA.Handles.measurementtextline(panel,textcounter),'XData',[scale*(SET(no).Measure(measureind).Y(end)+yt), xpos],'YData',[scale*(SET(no).Measure(measureind).X(end)+xt), ypos]);


%----------------------------------------------
function putroi_motion(panel)
%----------------------------------------------
%This function expands a balloon using otsu and polar transform
global DATA

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

DATA.CursorX = DATA.CursorX-mean(DATA.CursorX)+x;
DATA.CursorY = DATA.CursorY-mean(DATA.CursorY)+y;

DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;


%----------------------------------------------
function interp_motion(panel,type,sliceind,pointind)
%----------------------------------------------
%Motion function input is panel type sliceind and pointind. The types handled here
%are {EndoInterp,EpiInterp,RVendoInterp,RVepiInterp}. The indexes are
%determined in the buttondown so that the correct point is shifted.

global DATA SET

[no,slice,timeframe,x,y] = helperfunctions('getinterpparameters',panel);

%Check object ind in case of General Pen
currentobjectind = 1;
if strcmp(type,'GeneralPenInterp')
  currentobjectind = DATA.GeneralPenSettings.getcurrentobject;
end

helperfunctions('assignsetfield',no,type,'X',x,currentobjectind,timeframe,sliceind,pointind);
helperfunctions('assignsetfield',no,type,'Y',y,currentobjectind,timeframe,sliceind,pointind);

%write the results back to the contour field in the SET struct if we've
%placed more than two points
threshold = 2;
helperfunctions('storeinterptocontour',no,panel,type,timeframe,slice,currentobjectind,threshold);

drawfunctions('drawinterp',panel,type)

if strcmp(type(1:end-6),'LA') && ismember(timeframe,[SET(no).EST, SET(no).EDT])
  drawfunctions('drawatrialparameters',type(1:end-6),SET(no).ImageViewPlane);
end

%---------------------------------
function glarotatehandle_Motion(panel)
%---------------------------------
%Motion fcn for using GLA view rotation handle
global SET DATA

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
x = x/scale;
y = y/scale;
glaangle = mod(atan2(SET(no).ResolutionX*(y-SET(no).HLA.slice), ...
  SET(no).ResolutionY*(x-SET(no).VLA.slice))+pi/2,pi)-pi/2;

%We draw a line with markers the endpoint markers can not be seen as they
%are outside the field of view
DATA.Handles.cursor.XData = scale*(x+1000/SET(no).ResolutionY*cos(glaangle)*[1 0 -1]);
DATA.Handles.cursor.YData = scale*(y+1000/SET(no).ResolutionX*sin(glaangle)*[1 0 -1]);
%------------------------------
function [Xout,Yout,incr] = scale_motion(panel,startrad)
%-----------------------
global DATA

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system

X = DATA.CursorX;
Y = DATA.CursorY;

xc = mean(X);
yc = mean(Y);
a = norm([x,y]-[xc,yc]);

incr = a/startrad;

X = xc+(X-xc)*incr;
Y = yc+(Y-yc)*incr;

DATA.Handles.cursor.YData = X;
DATA.Handles.cursor.XData = Y;

if nargout>0
  Xout=X;
  Yout=Y;
end

%------------------------------
function [Xout,Yout] = translate_motion(panel,xstart,ystart)
%-----------------------
global DATA

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
X = DATA.CursorX-xstart+x;
Y = DATA.CursorY-ystart+y;

DATA.Handles.cursor.YData = X;
DATA.Handles.cursor.XData = Y;

if nargout>0
  Xout = X;
  Yout = Y;
end

%----------------------------------------------
function point_motion(panel)
%----------------------------------------------
%New motion function that handles points.

global DATA

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

%plot clicked position
DATA.Handles.cursor.XData = x;
DATA.Handles.cursor.YData = y;

%----------------------------------------------
function pan_motion(panel,xwstart,ywstart)
%----------------------------------------------
%New motion function for panning. Uses cursorX and cursor Y to store the
%initial positions of XLim and YLim

global DATA SET NO 

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

xw = x-DATA.Handles.imageaxes(panel).XLim(1);
yw = y-DATA.Handles.imageaxes(panel).YLim(1);

dx = xw - xwstart;
dy = yw - ywstart;

%move
xc = DATA.CursorX - dx;
yc = DATA.CursorY - dy;

%this is to halt execution assert that the above is set before running
panelslinked = find(ismember(DATA.ViewPanels,SET(DATA.ViewPanels(panel)).Linked));
paneltype = DATA.ViewPanelsType{DATA.CurrentPanel};

if ismember(paneltype,{'trans3DP','sag3DP','cor3DP'})
  panelslinked = [DATA.CurrentPanel segment3dp.viewfunctions('findspeedim')];
end
if isequal(paneltype,'speedim')
  switch SET(NO).LevelSet.Pen.Color
    case 'r'
      type = 'trans3DP';
    case 'g'
      type = 'sag3DP';
    case 'b'
      type = 'cor3DP';
  end
  panelslinked = [find(strcmp(DATA.ViewPanelsType,type)) DATA.CurrentPanel];
end

for actpanel = panelslinked
  set(DATA.Handles.imageaxes(actpanel),'XLim' ,xc,'YLim',yc)
end

drawfunctions('drawselectedframe',panel);

%----------------------------------------------
function pen_motion(panel)
%----------------------------------------------
%New motion function that handles all the different pen tools.

global DATA
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

DATA.CursorX = cat(1,DATA.CursorX,x);
DATA.CursorY = cat(1,DATA.CursorY,y);

DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;

