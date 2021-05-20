function [varargout] = LVrotation(varargin)
%function to define the LV rotation


macro_helper(varargin{:});
if nargin == 0 || isempty(varargin{1})
  varargin{1} = 'init';
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard


%------------
function init %#ok<DEFNU>
%------------
%open the gui to define LV rotation
global DATA NO SET

gui = mygui('LVrotation.fig');
DATA.GUI.LVrotation = gui;
gui.no = NO;
gui.slice = round(0.5*(SET(gui.no).StartSlice+SET(gui.no).EndSlice)); 

%set rotation and slice
set(gui.handles.rotationslider,'value',SET(gui.no).SectorRotation);
set(gui.handles.slicetext,'String',dprintf('Slice %d',gui.slice));
if SET(NO).ZSize > 1
  set(gui.handles.sliceslider,'Min',1,'Max',...
    SET(NO).ZSize,'Value',SET(NO).ZSize-gui.slice+1,'SliderStep',...
    [1 3]/(SET(NO).ZSize));
else
  set(gui.handles.sliceslider,'Visible','off');
  set(gui.handles.slicetext,'Visible','off');
end
segment('recursekeypressfcn',gui.fig,@(hObject,eventdata)LVrotation('keypressed',eventdata))
%plot image
plotimage;

if ~DATA.Testing
  %block any further analysis until this interface is closed
  uiwait(gui.fig);
else
  ok_Callback;
end

%--------------------------
function requestfocus 
%--------------------------
global DATA
gui = DATA.GUI.LVrotation;
warning off
jFig = get(gui.fig,'JavaFrame'); 
jFig.requestFocus; 
warning on;

%-------------------------
function keypressed(evt) %#ok<DEFNU>
%-------------------------
%%Keypress function
global DATA
gui = DATA.GUI.LVrotation;
% take away focus from sliders
requestfocus;
switch evt.Key
  case {'downarrow','uparrow'}
    %move slice
    oldvalue = mygetslider(gui.handles.sliceslider);
    h = gui.handles.sliceslider;
    stepvalue = (h.Max - h.Min)*h.SliderStep(1);
    if contains(evt.Key,'up')
      newvalue = oldvalue + stepvalue;
    else
      newvalue = oldvalue - stepvalue;
    end
    if newvalue > h.Min && newvalue < h.Max
      gui.handles.sliceslider.Value = newvalue;
    elseif newvalue > h.Max
      gui.handles.sliceslider.Value = h.Max;
    elseif newvalue < h.Min
      gui.handles.sliceslider.Value = h.Min;
    end
    sliceslider_Callback
  
    
  case {'rightarrow','leftarrow'}
    oldvalue = mygetslider(gui.handles.rotationslider);
    h = gui.handles.rotationslider;
    stepvalue = (h.Max - h.Min)*h.SliderStep(1);
    if contains(evt.Key,'right')
      newvalue = oldvalue + stepvalue;
    else
      newvalue = oldvalue - stepvalue;
    end
    if newvalue > h.Min && newvalue < h.Max
      gui.handles.rotationslider.Value = newvalue;
    elseif newvalue > h.Max
      gui.handles.rotationslider.Value = h.Max;
    elseif newvalue < h.Min
      gui.handles.rotationslider.Value = h.Min;
    end
    rotationslider_Callback
  
  otherwise
    return

end

%-----------------
function plotimage
%-----------------
%plot the image

global DATA SET
gui = DATA.GUI.LVrotation;

tf = SET(gui.no).CurrentTimeFrame;
SET(gui.no).SectorRotation = mygetvalue(gui.handles.rotationslider);

%Plot image
temp = SET(gui.no).IM(:,:,tf,gui.slice);
temp = min(max(temp,0),1);
imsize = size(temp);

%force true color
image(calcfunctions('remapuint8',temp,gui.no,calcfunctions('returnmapping',gui.no,true)),...
  'parent',gui.handles.imageaxes);
axis(gui.handles.imageaxes,'image','off');

%Upsample model
if not(DATA.Pref.RadialProfiles == DATA.NumPoints)
  [endox,endoy] = calcfunctions('resamplemodel',SET(gui.no).EndoX(:,tf,gui.slice),SET(gui.no).EndoY(:,tf,gui.slice),DATA.Pref.RadialProfiles);
  if ~isempty(SET(gui.no).EpiX)
    [epix,epiy] = calcfunctions('resamplemodel',SET(gui.no).EpiX(:,tf,gui.slice), SET(gui.no).EpiY(:,tf,gui.slice),DATA.Pref.RadialProfiles);
  end
else
  endox = SET(gui.no).EndoX(:,tf,gui.slice);
  endoy = SET(gui.no).EndoY(:,tf,gui.slice);
  if ~isempty(SET(gui.no).EpiX)
    epix = SET(gui.no).EpiX(:,tf,gui.slice);
    epiy = SET(gui.no).EpiY(:,tf,gui.slice);
  end
end
if isempty(SET(gui.no).EpiX)
  epix = NaN;
  epiy = NaN;
end

%Plot contours
hold(gui.handles.imageaxes,'on');
gui.handles.endocontour = plot(gui.handles.imageaxes,endoy,endox,'r-');
set(gui.handles.endocontour,'linewidth',DATA.Pref.LineWidth);

if ~isempty(SET(gui.no).EpiX)
  gui.handles.epicontour = plot(gui.handles.imageaxes,epiy,epix,'g-');
  set(gui.handles.epicontour,'linewidth',DATA.Pref.LineWidth);
end


%Plot sectors
hold(gui.handles.imageaxes,'on');
if isnan(epix(1))
  %--- Only endo exist => draw only endo
  %get positions.
  [gui.meanx,gui.meany,gui.sectors] = calcfunctions('findmeaninsectorslice','endo',DATA.Pref.RadialProfiles,...
    tf,gui.slice,1,gui.no);
  
  xpos = endoy(gui.sectors(1));
  ypos = endox(gui.sectors(1));
  %Draw an extra long line
  h = plot(gui.handles.imageaxes,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,1.5*xpos-0.5*gui.meany))], ...
    [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,1.5*ypos-0.5*gui.meany))],'y-');
  set(h,'linewidth',DATA.Pref.LineWidth);

else
  %--- Epi exists => draw both
  %get positions.
  [gui.meanx,gui.meany,gui.sectors] = calcfunctions('findmeaninsectorslice','epi',DATA.Pref.RadialProfiles,...
    tf,gui.slice,1,gui.no);  
  xpos = epiy(gui.sectors(1));  
  ypos = epix(gui.sectors(1));
  %Draw an extra long line
  gui.handles.line = plot(gui.handles.imageaxes,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,1.5*xpos-0.5*gui.meany))], ...
    [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,1.5*ypos-0.5*gui.meanx))],'y-');
  set(gui.handles.line,'linewidth',DATA.Pref.LineWidth);
end
hold off



%-------------------------------
function rotationslider_Callback %#ok<DEFNU>
%-------------------------------
%Callback for rotation slider. Update slice image.
plotimage;



%----------------------------
function sliceslider_Callback 
%----------------------------
%Callback for slider to toggle slice
global DATA SET
gui = DATA.GUI.LVrotation;
requestfocus;

gui.slice = SET(gui.no).ZSize-round(mygetvalue(gui.handles.sliceslider))+1;
set(gui.handles.sliceslider,'Value',round(mygetvalue(gui.handles.sliceslider)));
set(gui.handles.slicetext,'String',dprintf('Slice %d',gui.slice));
plotimage;


%---------------------------------------
function rotationfromannotation_Callback %#ok<DEFNU>
%---------------------------------------
%Finds rotation by looking at RV insertion points.

global DATA SET
gui = DATA.GUI.LVrotation;
requestfocus;

v = get(gui.handles.rotationfromannotationcheckbox,'value');
if v
  pos = rotationfromannotationhelper(gui.no);
  if isempty(pos)
    set(gui.handles.rotationfromannotationcheckbox,'value',0);
  end
end

%update slider
set(gui.handles.rotationslider,'value',SET(gui.no).SectorRotation);
%call to graphically update rotation & bullseye
plotimage;

%--------------------------------------
function pos = rotationfromannotationhelper(no)
%--------------------------------------
%Find suitable sector rotation based on RV insertion points

global SET

%Find slices with RV insertion points
slices = false(1,SET(no).ZSize);
for loop = 1:length(SET(no).Point.Z)
  if isequal(SET(no).Point.Label{loop},'RV insertion') || isequal(SET(no).Point.Label{loop},'P1')|| isequal(SET(no).Point.Label{loop},'P2')
    slices(SET(no).Point.Z(loop)) = true;
  end
end

%Find slices
pos = find(slices);
if isempty(pos)
  mywarning(dprintf('No RV points found. Define two RV annotation points.')); 
  return
end
pos = reportbullseye('sectorrotationhelper',no);


%-------------------
function ok_Callback
%-------------------
%start the bullseye analsysia and close the LV rotation interface

global DATA SET
gui = DATA.GUI.LVrotation;

%if user click ok for 0 in LV rotation, define it as 0.1 degrees in order
%to remember that the user have checked the LV rotation
if SET(gui.no).SectorRotation == 0
  SET(gui.no).SectorRotation = 0.1;
end

close_Callback;
reportbullseye('startbullseye');

%------------------------
function cancel_Callback %#ok<DEFNU>
%-----------------------
%cancel the analysis and close the LV rotation interface
close_Callback;


%----------------------
function close_Callback 
%----------------------
%Properly close the GUI.

global DATA

try
  DATA.GUI.LVrotation = close(DATA.GUI.LVrotation);  %close the gui
catch %#ok<CTCH>
  delete(gcbf)
end

DATA.GUI.LVrotation= [];