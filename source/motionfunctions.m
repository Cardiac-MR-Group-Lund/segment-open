 function varargout = motionfunctions(varargin)
% Functions for cursor motion 
% Klas 
%Invoke subfunction

macro_helper(varargin{:}); %future macro recording use
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%---------------------------
function select_motion(panel) %#ok<DEFNU>
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
function [contrastout,brightnessout] = contrast_motion(panel,slice,xstart,ystart,xsize,ysize) %#ok<DEFNU>
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
    im = calcfunctions('calcmontageviewim', ...
        no,fliplr(DATA.ViewPanelsMatrix{panel}),segmentedonly,calcfunctions('returnmapping',no),contrast,brightness,[],SET(no).CurrentTimeFrame);
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
function measure_motion(panel,pointind) %#ok<DEFNU>
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
function interp_motion(panel,type,sliceind,pointind) %#ok<DEFNU>
%----------------------------------------------
%Motion function input is panel type sliceind and pointind. The types handled here
%are {EndoInterp,EpiInterp,RVendoInterp,RVepiInterp}. The indexes are
%determined in the buttondown so that the correct point is shifted.

global DATA SET
no = DATA.ViewPanels(panel);

%get closest interpolation point
[yclick,xclick] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

scale = viewfunctions('getscale',panel);
slice = viewfunctions('clickedslice',panel,yclick,xclick);
slices = viewfunctions('slicesinpanel',panel);

%normalize clicked position to contour domain
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
x = xclick/scale - (xl-1)*SET(no).XSize;
y = yclick/scale - (yl-1)*SET(no).YSize;

SET(no).([type,'X']){SET(no).CurrentTimeFrame,sliceind}(pointind) = x;
SET(no).([type,'Y']){SET(no).CurrentTimeFrame,sliceind}(pointind) = y;

if SET(no).([type(1:end-6),'InterpOngoing'])
  numinterppoints = length(SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice});
  numpoints = floor(DATA.NumPoints/15)*(numinterppoints-1);
  if numpoints > DATA.NumPoints
    numpoints = DATA.NumPoints;
  end
else
  numpoints = DATA.NumPoints;
end

%here we need to interpolate and store curve into contours
%removes duplicate points and resamples the contour
[x,y] = calcfunctions('resamplecurve',SET(no).([type,'X']){SET(no).CurrentTimeFrame,sliceind},...
    SET(no).([type,'Y']){SET(no).CurrentTimeFrame,sliceind},numpoints-1);

%write the results back to the SET struct
if length(x)>2
  switch type
    case 'Roi'
        %would be nice to introduce RoiInterpX...
    otherwise
      if SET(no).([type(1:end-6),'InterpOngoing'])
        expectedlength = length(SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice));
        if numpoints < expectedlength
          % fill up with NaNs        
          x = cat(2,x,nan(1,expectedlength-numpoints));
          y = cat(2,y,nan(1,expectedlength-numpoints));
        end
      end
      SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,sliceind)= [x,x(1)];
      SET(no).([type(1:end-6),'Y'])(:,SET(no).CurrentTimeFrame,sliceind)= [y,y(1)];
      drawfunctions('drawcontours',panel,type(1:end-6))
  end
end
drawfunctions('drawinterp',panel,type)

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
function [Xout,Yout,incr] = scale_motion(panel,startrad) %#ok<DEFNU>
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
function [Xout,Yout] = translate_motion(panel,xstart,ystart) %#ok<DEFNU>
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

%--------------------------------
function viewportpoints_motion(ind) %#ok<DEFNU>
%--------------------------------
%Motion function for viewport and points

global DATA SET NO

coord = getclickedcoord(DATA.LevelSet.ViewPort);

if ~isnan(coord(1))
  x = coord(1)+DATA.LevelSet.Box.xmin-1; %This adds if the the volume is cropped
  y = coord(2)+DATA.LevelSet.Box.ymin-1;
  z = coord(3)+DATA.LevelSet.Box.zmin-1;

  SET(NO).Point.X(ind) = x;
  SET(NO).Point.Y(ind) = y;
  SET(NO).Point.Z(ind) = z;
  
end

segment3dp.tools('addpointstoviewport');

%See if we should update the line
if ~isempty(SET(NO).Line3D.Points)
  if ismember(ind,SET(NO).Line3D.Points{1})
    segment3dp.linetools('createline','',SET(NO).Line3D.Points{1}); %Later possibility to have more lines
  end
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
function pan_motion(panel,xwstart,ywstart) %#ok<DEFNU>
%----------------------------------------------
%New motion function for panning. Uses cursorX and cursor Y to store the
%initial positions of XLim and YLim

global DATA SET

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
  panelslinked = [DATA.CurrentPanel viewfunctions('findspeedim')];
end
for actpanel = panelslinked 
  set(DATA.Handles.imageaxes(actpanel),'XLim' ,xc,'YLim',yc) 
end

drawfunctions('drawselectedframe',panel);

%----------------------------------------------
function pen_motion(panel) %#ok<DEFNU>
%----------------------------------------------
%New motion function that handles all the different pen tools.

global DATA
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

DATA.CursorX = cat(1,DATA.CursorX,x);
DATA.CursorY = cat(1,DATA.CursorY,y);

DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY; 

%--------------------------
function select3dp_motion  %#ok<DEFNU>
%--------------------------
%Motion function for select
global DATA

segment3dp.tools('storeclickedposition');

%segment3dp.tools('update3DP')
for p = find(DATA.ViewPanels)  
  if ~isequal(DATA.ViewPanelsType{p},'viewport') %no need to draw viewport
    drawfunctions('drawpanel',p);
  end
end

for loop = 1:length(DATA.ViewPanels)
  if ~isequal(DATA.ViewPanelsType{loop},'speedim')
    drawfunctions('drawtext',loop);
  end
end

%--------------------------
function fastmarching_motion  %#ok<DEFNU>
%--------------------------
%Motion function for graphical updates using fastmarching tool
global DATA SET NO

if DATA.Run
  return
end

[r,g,b] = segment3dp.tools('getclickedposition3DP');

SET(NO).LevelSet.View.RSlice = r;
SET(NO).LevelSet.View.GSlice = g;
SET(NO).LevelSet.View.BSlice = b;

%Convert to Segment indices
[x,y,z] = segment3dp.tools('rgb2xyz',r,g,b);

%Extrac arrival time at mouse click
arrivaltime = DATA.LevelSet.ArrivalTime(x,y,z);

if isnan(arrivaltime)
  arrivaltime = max(DATA.LevelSet.ArrivalTime(:));
end

imh_r = DATA.Handles.imagehandles(strcmp(DATA.ViewPanelsType,'trans3DP'));
imh_g = DATA.Handles.imagehandles(strcmp(DATA.ViewPanelsType,'sag3DP'));
imh_b = DATA.Handles.imagehandles(strcmp(DATA.ViewPanelsType,'cor3DP'));
%Red
im = segment3dp.tools('getimage','r');
bw = segment3dp.tools('getimagehelper',SET(NO).LevelSet.BW,'r')>uint8(127)  |  (segment3dp.tools('getimagehelper',DATA.LevelSet.ArrivalTime,'r')<arrivaltime);
set(imh_r,'cdata',segment3dp.tools('levelsetremapandoverlay',im,bw,uint8(0*bw)));

%Green
im = segment3dp.tools('getimage','g');
bw = segment3dp.tools('getimagehelper',SET(NO).LevelSet.BW,'g')>uint8(127)  |  (segment3dp.tools('getimagehelper',DATA.LevelSet.ArrivalTime,'g')<arrivaltime);
set(imh_g,'cdata',segment3dp.tools('levelsetremapandoverlay',im,bw,uint8(0*bw)));

%Blue
im = segment3dp.tools('getimage','b');
bw  = segment3dp.tools('getimagehelper',SET(NO).LevelSet.BW,'b')>uint8(127)  |  (segment3dp.tools('getimagehelper',DATA.LevelSet.ArrivalTime,'b')<arrivaltime);
set(imh_b,'cdata',segment3dp.tools('levelsetremapandoverlay',im,bw,uint8(0*bw)));

segment3dp.tools('updatergbsliders');

%----------------------
function pen3dp_motion(tool,update)  %#ok<DEFNU>
%----------------------
%Manual draw

global SET NO DATA
persistent counter

if isempty(counter)
  counter = 0;
end

%second input asserts update of all panels
if nargin<2
    update=0;
end

[r,g,b] = segment3dp.tools('getclickedposition3DP');

[x,y,z] = segment3dp.tools('rgb2xyz',r,g,b);

%Old code to ensure it
%Check if valid position & move
%x = min(max(round(x),SET(NO).LevelSet.Pen.XSize),SET(NO).XSize-SET(NO).LevelSet.Pen.XSize);
%y = min(max(round(y),SET(NO).LevelSet.Pen.YSize),SET(NO).YSize-SET(NO).LevelSet.Pen.YSize);
%z = min(max(round(z),SET(NO).LevelSet.Pen.ZSize),SET(NO).ZSize-SET(NO).LevelSet.Pen.ZSize);

x = min(max(round(x),1),SET(NO).XSize);
y = min(max(round(y),1),SET(NO).YSize);
z = min(max(round(z),1),SET(NO).ZSize);

%Convert back to RGB for storage
[r,g,b] = segment3dp.tools('xyz2rgb',x,y,z);
SET(NO).LevelSet.View.RSlice = r;
SET(NO).LevelSet.View.GSlice = g;
SET(NO).LevelSet.View.BSlice = b;

%Find index
xind = SET(NO).LevelSet.Pen.X+x;
yind = SET(NO).LevelSet.Pen.Y+y;
zind = SET(NO).LevelSet.Pen.Z+z;

%Find inside positions
logind = ...
  (xind>0) & (xind<=SET(NO).XSize) & ...
  (yind>0) & (yind<=SET(NO).YSize) & ...
  (zind>0) & (zind<=SET(NO).ZSize);
  
%Remove outside
xind = xind(logind);
yind = yind(logind);
zind = zind(logind);

try
  motionindex = sub2ind([SET(NO).XSize SET(NO).YSize SET(NO).ZSize],xind,yind,zind);
catch
  disp('outside image');
  set(DATA.fig,'windowbuttonmotionfcn',@DATA.toggleplaceholdermotion)
  return;
end

if isequal(tool,'localthreshold')
  %Normal draw/erase
  SET(NO).LevelSet.BW(motionindex) = segment3dp.tools('calcthreshold3DP',SET(NO).LevelSet.SpeedIM(motionindex));
else
  switch DATA.CurrentTool
    case 'rubber'
      DATA.LevelSet.Man(motionindex) = int8(-1);
      SET(NO).LevelSet.BW(motionindex) = min(SET(NO).LevelSet.BW(motionindex),SET(NO).LevelSet.Pen.Value(logind));      
    case 'draw'
      DATA.LevelSet.Man(motionindex) = int8(1);
      SET(NO).LevelSet.BW(motionindex) = max(SET(NO).LevelSet.BW(motionindex),SET(NO).LevelSet.Pen.Value(logind));
  end
end

DATA.toggleplaceholdermotion

%For speed sake only update every 4th
counter = mod(counter+1,8);
if counter==1 || update
  DATA.LevelSet.bwupdated = true;
  segment3dp.tools('update3DP')  
  %segment3dp.tools('bwupdated'); %this is way to slow!
else
  drawfunctions('drawimages',DATA.CurrentPanel)
end

segment3dp.tools('updatergbsliders');


