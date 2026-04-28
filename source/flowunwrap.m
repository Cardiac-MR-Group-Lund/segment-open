function [varargout] = flowunwrap(varargin)
%---------------------------------
% GUI for doing manual and automatic phase unwrapping in 2D and 3D
% flow datasets.

% Started by Johannes Toger 090710

%#ok<*GVMIS> 

if nargin == 0
  varargin = {'init'};
end

[varargout{1:nargout}] = feval(varargin{:});

%-------------------
function init(force,no)
%-------------------
%Initialize the Flow Unwrap GUI
%

global SET DATA NO 

if nargin==0
  force = false;%force is only used by maketest to be able to test without user interaction
end
if nargin < 2
  no = NO;
end

% Check for [] and 
if no == 0
  myfailed('No data loaded. Please load flow data.')
  return
end

availableflows = checkflowdata(no);

if isempty(availableflows)
  myfailed('No flows available in current stack. Please load and select flow data.')
  return
end

myworkon;

% Launch GUI
DATA.GUI.FlowUnwrap = mygui('flowunwrap.fig');
gui = DATA.GUI.FlowUnwrap;

gui.currentslice = SET(no).CurrentSlice;
gui.totalslices = SET(no).ZSize;
gui.currenttimeframe = SET(no).CurrentTimeFrame;
gui.totaltimeframes = SET(no).TSize;
gui.TIncr = SET(no).TIncr;
gui.currentflowno = 1;
gui.playlock = false; % Lock so that only one play is done at one time
gui.zoomlevel = 0;
gui.zoomfactor = 1.2; % Length scaling factor for zoom

% Compute max zoom factor so that there is at least 20 pixels in each
% direction
minpixels = 20;
zoomfactorx = floor(-log(minpixels/SET(no).XSize)/log(gui.zoomfactor));
zoomfactory = floor(-log(minpixels/SET(no).YSize)/log(gui.zoomfactor));

gui.maxzoomlevel = min(zoomfactorx, zoomfactory);

% Hide some deprecated things
%initializehidejumpsizeselect;

% Make a working copy of the flow data.
gui.workingcopy = struct;
gui.workingcopy.VENC = SET(no).VENC;
gui.workingcopy.SETnumbers = availableflows;
gui.workingcopy.flows = cell(size(availableflows));

for flownumber = 1:length(availableflows)
  gui.workingcopy.flows{flownumber} = SET(availableflows(flownumber)).IM;
end

% Rescale to -1 = -VENC, +1 = +VENC
for flownumber = 1:length(gui.workingcopy.flows)
  gui.workingcopy.flows{flownumber} = ...
                            2 * (gui.workingcopy.flows{flownumber} - 1/2);
end

% If this dataset has been auto-unwrapped with the old algorithm, the VENC
% will have changed. This new algorithm saves the original venc in
% SET().Flow.OriginalVENC, but the old one didn't. Therefore, we have to
% ask for it if it's not there, and then correct for it. Bah.
%
if not(force)
  try
    originalVENC = SET(no).Flow.OriginalVENC;
  catch exception %#ok<NASGU>
    vencquestion = dprintf('What was the original VENC in this dataset?');
    dlgtitle = dprintf('Please set the original VENC.');
    default = num2str(SET(no).VENC);
    userinput = myinputdlg({vencquestion}, dlgtitle, 1, {default});
    if ~isempty(userinput)
      originalVENC = str2double(userinput{1});
    else
      closewindow; 
      return
    end
  end
else
  originalVENC=SET(no).VENC;
end

gui.workingcopy.OriginalVENC = originalVENC;

% Restore to original range: -1 = -originalVENC, +1 = +originalVENC
vencfactor = gui.workingcopy.VENC / gui.workingcopy.OriginalVENC;
for flownumber = 1:length(gui.workingcopy.flows)
  gui.workingcopy.flows{flownumber} = gui.workingcopy.flows{flownumber}*vencfactor;
end

% The working copy is stored in units of gui.workingcopy.VENC, which is
% equal to gui.workingcopy.OriginalVENC*3.
gui.workingcopy.VENC = 3*gui.workingcopy.OriginalVENC;
for flownumber = 1:length(gui.workingcopy.flows)
  gui.workingcopy.flows{flownumber} = gui.workingcopy.flows{flownumber} / 3;
end

% Initialize current point
midx = round(SET(no).XSize / 2);
midy = round(SET(no).YSize / 2);
gui.currentimagepoint = [midx midy];

% Set the expected jump size (in OriginalVENCs)
% Don't change this in code, since there is no updater to reflect this to
% the GUI (the radio button group 'Jump size').
% The only reason for allowing this to be changed is that we had some
% Hadamard 3Dflow datasets with 1-venc jumps (unknown reason).
%
% Don't change this. Really. /JT.
gui.expectedjumpsize = 2;

% Initialize gui
initializegui;

% Set keypressed callback. slider4 is an invisible control used only for
% this. See function refocus for details.
for uielement = [gui.fig gui.handles.slider4]
  set(uielement, 'KeyPressFcn', ...
      @(src, event)flowunwrap('keypress_Callback', src, event));
end

% Set velocity image button down callback. This must be done here, since 
% the image object doesn't exist until after initialization
% (initializetemporalplot)
set(gui.handles.velocityimage, ...
    'ButtonDownFcn', 'flowunwrap(''velocityimage_ButtonDown'')');

% Button group for jump size callback. 
% For some reason you can't set this in GUIDE.
set(gui.handles.jumpsizebuttongroup, 'SelectionChangeFcn', ...
    @(src,event)flowunwrap('jumpsizebuttongroup_Callback', src, event));

% Update all gui elements (first kick)
updatefullgui;
refocus;
myworkoff;
set(gui.fig,'pointer','arrow');



%---------------------
function initializegui
%---------------------
% Performs startup tasks, like setting slider limits. A full update is
% also done after this initialization, so there is no need to duplicate
% code in both.
%

initializehideunnecessary;
initializesliders;
initializevelocityimage;
initializetemporalplot;
initializeflowdirectionlistbox;

%------------------------------------
function initializehidejumpsizeselect
%------------------------------------
% Hides the jump size selector.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

set(gui.handles.jumpsizebuttongroup, 'Visible', 'off');


%---------------------------------
function initializehideunnecessary
%---------------------------------
% Hide unnecessary UI elements:
% 
% * If there is only one slice, hide slice selection.
%
% * If zoom is 0, hide pan sliders.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

if gui.totalslices == 1
  set(gui.handles.sliceslider, 'Visible', 'off');
  set(gui.handles.slicetext, 'Visible', 'off');
end

if gui.zoomlevel == 0
  set(gui.handles.panimagehorizontalslider, 'Visible', 'off');
  set(gui.handles.panimageverticalslider, 'Visible', 'off');
end


%-------------------------
function initializesliders
%-------------------------
% Initializes the sliders.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

timeframesliderstep = [1 5] / gui.totaltimeframes;

set(gui.handles.timeframeslider, 'Min', 1);
set(gui.handles.timeframeslider, 'Max', gui.totaltimeframes);
set(gui.handles.timeframeslider, 'SliderStep', timeframesliderstep);

% No slice slider needed if there is only one slice.
if gui.totalslices > 1 
  slicesliderstep = [1 5] / gui.totalslices;

  set(gui.handles.sliceslider, 'Min', 1);
  set(gui.handles.sliceslider, 'Max', gui.totalslices);
  set(gui.handles.sliceslider, 'SliderStep', slicesliderstep);
end

%-------------------------------
function initializevelocityimage
%-------------------------------
% Initializes the velocity image
%

global DATA
gui = DATA.GUI.FlowUnwrap;

gui.handles.velocityimage = image(getcurrentvelocityslice, ...
                  'parent', gui.handles.velocityimageaxes, ...
                  'CDataMapping','scaled');

hold(gui.handles.velocityimageaxes, 'on');

% Set color scale so that white = +2*OriginalVENC
% and black = -2*OriginalVENC
% Remember: working copy has VENC = 3*originalVENC.
vencratio = gui.workingcopy.OriginalVENC/gui.workingcopy.VENC;
colorlimits =  [-1 1]*vencratio;

set(gui.handles.velocityimageaxes, 'CLim', colorlimits);
set(gui.fig, 'Colormap', gray);                

set(gui.handles.velocityimageaxes, 'XTick', []);
set(gui.handles.velocityimageaxes, 'YTick', []);
set(gui.handles.velocityimageaxes, 'DataAspectRatio', [1 1 1]);

% Create current pixel overlay
[X, Y] = getoverlayplotcoords;

overlayhandle = plot(X, Y, 'r-', 'parent', gui.handles.velocityimageaxes);
gui.handles.currentpointoverlayplot = overlayhandle;

%------------------------------
function initializetemporalplot
%------------------------------
% Initializes the temporal plot
%

global DATA
gui = DATA.GUI.FlowUnwrap;

% Main plot
pixel_over_time = getcurrentpixelovertime;
timevector = (0:(gui.totaltimeframes-1))*gui.TIncr*1e3; % plot in ms

plots = plot(timevector, pixel_over_time, 'k+', ...
             timevector, pixel_over_time, 'b-', ...
             'parent', gui.handles.temporalplotaxes);

set(gui.handles.temporalplotaxes, 'YLimMode', 'manual','YColor',DATA.GUISettings.ForegroundColor);
set(gui.handles.temporalplotaxes, 'YLim', [-2 2]);

tlimits = [0 (gui.totaltimeframes-1)*gui.TIncr*1e3];
set(gui.handles.temporalplotaxes, 'XLimMode', 'manual','XColor',DATA.GUISettings.ForegroundColor);
set(gui.handles.temporalplotaxes, 'XLim', tlimits);

gui.handles.temporalplots = plots;

% Labels
timestring = makeunitstring(dprintf('Time'),'ms');
set(get(gui.handles.temporalplotaxes, 'XLabel'), 'String', timestring);
set(get(gui.handles.temporalplotaxes, 'Ylabel'), 'String', dprintf('Velocity (original VENC)'));

% Current timeframe overlay (vertical line)
hold(gui.handles.temporalplotaxes, 'on');

[X, Y] = gettimeframemarkercoords;
timemarkerhandle = plot(X, Y, 'r--', ...
                        'parent', gui.handles.temporalplotaxes);
                      
gui.handles.timeframemarkerplot = timemarkerhandle;

%---------------------------------------
function initializeflowdirectionlistbox
%--------------------------------------
% Initializes the listbox where flow direction can be selected.
%

global DATA SET
gui = DATA.GUI.FlowUnwrap;

directiontexts = cell(size(gui.workingcopy.SETnumbers));
directiontexttemplate = '%d: %s, %s';

for textnumber = 1:length(directiontexts)
  thistext = sprintf(directiontexttemplate, textnumber, ...
                     SET(gui.workingcopy.SETnumbers(textnumber)).ImageType, ...
                     SET(gui.workingcopy.SETnumbers(textnumber)).ImageViewPlane);
                   
  directiontexts{textnumber} = thistext;
end

set(gui.handles.flowdirectionlistbox, 'String', directiontexts);

%---------------------
function updatefullgui
%---------------------
% Updates all gui elements to reflect current state.
%

updatevelocityimage;
updatetemporalplotaxes;
updatesliders;
updatetextelements;

%---------------------------
function updatevelocityimage
%---------------------------
% Updates the axes showing the velocity image
%

global DATA;
gui = DATA.GUI.FlowUnwrap;

set(gui.handles.velocityimage, 'CData', getcurrentvelocityslice);

%-------------------------------------
function updatevelocityimagezoomandpan
%-------------------------------------
% Updates zoom and pan in the velocity image
%

global DATA
gui = DATA.GUI.FlowUnwrap;

% Zoom and pan
XData = get(gui.handles.velocityimage, 'XData');
YData = get(gui.handles.velocityimage, 'YData');
sizex = XData(2);
sizey = YData(2);

xpan = mygetslider(gui.handles.panimagehorizontalslider);
ypan = mygetslider(gui.handles.panimageverticalslider);

if gui.zoomlevel == 0
  XLim = XData + 0.5;
  YLim = YData + 0.5;
else
  numpixelsx = floor(sizex / power(gui.zoomfactor, gui.zoomlevel));
  numpixelsy = floor(sizey / power(gui.zoomfactor, gui.zoomlevel));

  topleftx = (sizex - numpixelsx)*xpan;
  toplefty = (sizey - numpixelsy)*(1 - ypan);
  
  XLim = round([topleftx, topleftx + numpixelsx]) + 0.5;
  YLim = round([toplefty, toplefty + numpixelsy]) + 0.5;
end

set(gui.handles.velocityimageaxes, 'XLim', XLim);
set(gui.handles.velocityimageaxes, 'YLim', YLim);

%------------------------------
function updatetemporalplotaxes
%------------------------------
% Updates the temporal plot of the current pixel.
%

global DATA;
gui = DATA.GUI.FlowUnwrap;

pixel_over_time = getcurrentpixelovertime;

for plt = gui.handles.temporalplots
  set(plt, 'YData', pixel_over_time);
end

[X, ~] = gettimeframemarkercoords; 
set(gui.handles.timeframemarkerplot, 'XData', X);

%-----------------------------
function updatepanimagesliders
%-----------------------------
% Updates the pan image sliders - invisible if zoomlevel is 0
%

global DATA;
gui = DATA.GUI.FlowUnwrap;

for panslider = [gui.handles.panimageverticalslider ...
                 gui.handles.panimagehorizontalslider]
  if gui.zoomlevel == 0
    set(panslider, 'Visible', 'off')
  else
    set(panslider, 'Visible', 'on');
  end
end

%---------------------
function updatesliders
%---------------------
% Updates the time and slice sliders to reflect current positions.
%

global DATA;
gui = DATA.GUI.FlowUnwrap;

set(gui.handles.sliceslider, 'Value', gui.currentslice);
set(gui.handles.timeframeslider, 'Value', gui.currenttimeframe);

updatetextelements;

%--------------------------
function updatetextelements
%--------------------------
% Updates all text elements in the GUI.
%

global DATA;
gui = DATA.GUI.FlowUnwrap;


slicetext = sprintf('%s: %d/%d', ...
                    dprintf('Slice'),...
                    gui.currentslice,gui.totalslices);
                        
timeframetext = sprintf('%s: %d/%d (%1.0f ms)', ...
                dprintf('Time frame'),...
                gui.currenttimeframe, gui.totaltimeframes, ...
                getcurrenttime);

set(gui.handles.timeframetext, 'String', timeframetext);
set(gui.handles.slicetext, 'String', slicetext);

%---------------------------------
function updatecurrentpointoverlay
%---------------------------------
% Updates the current point overlay (red square in velocity image)
%

global DATA
gui = DATA.GUI.FlowUnwrap;

[X, Y] = getoverlayplotcoords;

set(gui.handles.currentpointoverlayplot, 'XData', X);
set(gui.handles.currentpointoverlayplot, 'YData', Y);


%--------------------------------
function timeframeslider_Callback
%---------------------------------
% Timeframe slider callback
%

global DATA
gui = DATA.GUI.FlowUnwrap;

newtimeframe = mygetslider(gui.handles.timeframeslider);
newtimeframe = round(newtimeframe);
gui.currenttimeframe = newtimeframe;

updatesliders;
updatevelocityimage;
updatetemporalplotaxes;
refocus;

%-----------------------------
function sliceslider_Callback
%-----------------------------
% 
% Slice slider callback

global DATA
gui = DATA.GUI.FlowUnwrap;

newslice = mygetslider(gui.handles.sliceslider);
newslice = round(newslice);

gui.currentslice = newslice;

updatesliders;
updatevelocityimage;
updatetemporalplotaxes;
refocus;

%-------------------------------------
function flowdirectionlistbox_Callback 
%-------------------------------------
% 
% Callback for the flow direction listbox.
%
global DATA
gui = DATA.GUI.FlowUnwrap;

currentflowno = mygetlistbox(gui.handles.flowdirectionlistbox);

gui.currentflowno = currentflowno;

updatevelocityimage;
updatetemporalplotaxes;
refocus;

%--------------------------------
function velocityimage_ButtonDown 
%--------------------------------
% Callback for clicks in velocityimage.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

[y,x]= mygetcurrentpoint(gui.handles.velocityimageaxes); %EH: Note switched conventions.

x = round(x);
y = round(y);

sizes = size(gui.workingcopy.flows{1});
xsize = sizes(1);
ysize = sizes(2);

% Bounds check
if (x < 0.5) || (y < 0.5) || (x > xsize + 0.49) || (y > ysize + 0.49)
  return
end

gui.currentimagepoint = [x y];

updatetemporalplotaxes;
updatecurrentpointoverlay;
refocus;

%-------------------------------------
function keypress_Callback(src, event) 
%------------------------------------
% Handles keypress events
%
%
% Up      30  Slice -
% Down    31  Slice +
% Left    28  Timestep +
% Right   29  Timestep -
%
% wW      119 Pixel selector up
% aA      97  Pixel selector left  
% sS      115 Pixel selector down
% dD      100 Pixel selector right
%
% space   32  Play current slice
%
% shift       Wrap up
% control     Wrap down
%

global DATA
gui = DATA.GUI.FlowUnwrap;

% Ignore keypresses from other sources
if ~(src == gui.handles.slider4)
  return
end

switch event.Key
  case 'uparrow' % Up, slice -
    newslice = gui.currentslice - 1;
    
    if newslice < 1
      newslice = 1;
    end

    gui.currentslice = newslice;
    updatetemporalplotaxes;
    updatevelocityimage;
    updatesliders;
    
  case 'downarrow' % Down, slice +
    newslice = gui.currentslice + 1;
    
    if newslice > gui.totalslices
      newslice = gui.totalslices;
    end
    
    gui.currentslice = newslice;
    updatetemporalplotaxes;
    updatevelocityimage;
    updatesliders;
    
  case 'leftarrow' % Left, timestep -
    newtimeframe = gui.currenttimeframe - 1;
    
    if newtimeframe < 1
      newtimeframe = gui.totaltimeframes;
    end
    
    gui.currenttimeframe = newtimeframe;
    updatetemporalplotaxes;
    updatevelocityimage;
    updatesliders;
    
  case 'rightarrow' % Right, timestep +
    newtimeframe = gui.currenttimeframe + 1;
    
    if newtimeframe > gui.totaltimeframes
      newtimeframe = 1;
    end
    
    gui.currenttimeframe = newtimeframe;
    updatetemporalplotaxes;
    updatevelocityimage;
    updatesliders;
  
  case {'w', 'a', 's', 'd'} % WASD, pixel selection
    switch event.Key
      case 'w'
        movecurrentpixel('up');
      case 'a'
        movecurrentpixel('left');
      case 's'
        movecurrentpixel('down');
      case 'd'
        movecurrentpixel('right');
    end
    
    updatecurrentpointoverlay;
    updatetemporalplotaxes;
    
  case 'space' % space, play current slice
    playcurrentslice;
    
  case 'shift' % shift, wrap up
    unwrapcurrentpixel('up');
    updatetemporalplotaxes;
    updatevelocityimage;
    
  case 'control' % control, wrap down
    unwrapcurrentpixel('down');
    updatetemporalplotaxes;
    updatevelocityimage;
    
end
%------------------------------------------------
function jumpsizebuttongroup_Callback(src, event) 
%------------------------------------------------
% Jump size selector callback.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

if ~(src == gui.handles.jumpsizebuttongroup) % Only handle callbacks from the right figure.
  return
end

switch get(event.NewValue, 'Tag')
  case 'venc1radiobutton'
    gui.expectedjumpsize = 1;
  case 'venc2radiobutton'
    gui.expectedjumpsize = 2;
end

refocus;

%-------------------------------------
function autounwrapregionpushbutton_Callback 
%-------------------------------------
% Perform automatic phase unwrapping regionally

%Einar Heiberg

global DATA
gui = DATA.GUI.FlowUnwrap;

for flownumber = 1:length(gui.workingcopy.flows)
  if(ndims(gui.workingcopy.flows{flownumber}) < 3)
    disp('This flow is not time-resolved, can''t do auto-unwrap.');
    return;
  end
  
  currentflow = gui.workingcopy.flows{flownumber};
  
  jumpsize = gui.expectedjumpsize * gui.workingcopy.OriginalVENC / gui.workingcopy.VENC;
  jumpfraction = 0.5;
  unwrapinputcheck = false; 
  waitbarmessagetemplate = dprintf('Please wait, automatic phase unwrapping %d/%d.');
    
  switch ndims(currentflow)
    case 3 % 2D-flow
      myworkon;      
      roisize = 15; %"radius" of cube
      point = gui.currentimagepoint;
      regionx = (max(1,point(1)-roisize)):(min(size(currentflow,1),point(1)+roisize)); %Take roisize around centre pixel
      regiony = (max(1,point(2)-roisize)):(min(size(currentflow,2),point(2)+roisize)); %Take roisize around centre pixel
      
      npixels = numel(regionx)*numel(regiony);
      waitbarmessage = sprintf(waitbarmessagetemplate, ...
                               flownumber, length(gui.workingcopy.flows));
                             
      h = mywaitbarstart(npixels, waitbarmessage,[]);
      
      for xn = 1:numel(regionx)
        for yn = 1:numel(regiony)
          currentflow(regionx(xn), regiony(yn), :) = ...
                unwraptimeseries(currentflow(regionx(xn), regiony(yn), :), jumpsize, ...
                                 jumpfraction, unwrapinputcheck);

          h = mywaitbarupdate(h);
        end
      end
      mywaitbarclose(h);
      myworkoff;
      
    case 4 % 3D-flow
      myfailed('Not supported for 4D flow.');
      return;      
  end
  
  gui.workingcopy.flows{flownumber} = currentflow;  
end

updatevelocityimage;
updatetemporalplotaxes;
refocus;

%-------------------------------------
function autounwrappushbutton_Callback 
%-------------------------------------
% Perform automatic phase unwrapping on all flows.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

for flownumber = 1:length(gui.workingcopy.flows)
%for flownumber = 1:1
  if(ndims(gui.workingcopy.flows{flownumber}) < 3)
    disp('This flow is not time-resolved, can''t do auto-unwrap.');
    return;
  end
  
  currentflow = gui.workingcopy.flows{flownumber}; 
  
  jumpsize = gui.expectedjumpsize * gui.workingcopy.OriginalVENC / gui.workingcopy.VENC;
  jumpfraction = 0.5;
  unwrapinputcheck = false; 
  waitbarmessagetemplate = dprintf('Please wait, automatic phase unwrapping %d/%d.');
    
  switch ndims(currentflow)
    case 3 % 2D-flow
      myworkon;
      npixels = size(currentflow,1)*size(currentflow,2);
      waitbarmessage = sprintf(waitbarmessagetemplate, ...
                               flownumber, length(gui.workingcopy.flows));
                             
      h = mywaitbarstart(npixels, waitbarmessage,[]);
      
      for xn = 1:size(currentflow, 1)
        for yn = 1:size(currentflow, 2)
          currentflow(xn, yn, :) = ...
                unwraptimeseries(currentflow(xn, yn, :), jumpsize, ...
                                 jumpfraction, unwrapinputcheck);

          h = mywaitbarupdate(h);
        end
      end
      mywaitbarclose(h);
      myworkoff;
      
    case 4 % 3D-flow
      myworkon;
      npixels = size(currentflow, 1)*size(currentflow,2);
      waitbarmessage = sprintf(waitbarmessagetemplate, ...
                               flownumber, length(gui.workingcopy.flows));
      
      h = mywaitbarstart(npixels, waitbarmessage,[]);
      
      for xn = 1:size(currentflow, 1)
        for yn = 1:size(currentflow, 2)
          for sln = 1:size(currentflow, 4)
            currentflow(xn, yn, :, sln) = ...
                unwraptimeseries(currentflow(xn, yn, :, sln), jumpsize, ...
                                 jumpfraction, unwrapinputcheck);
            
            
          end
          h = mywaitbarupdate(h);
        end
      end
      mywaitbarclose(h);
      myworkoff;
  end
  
  gui.workingcopy.flows{flownumber} = currentflow;  
end

updatevelocityimage;
updatetemporalplotaxes;
refocus;

%_--------------------------------------
function applyandexitpushbutton_Callback 
%---------------------------------------
% Apply and Exit push button.
% 
% Saves the phase images to where they came from, overwriting.
% and using the VENC, double relative the original.
%

global DATA SET
gui = DATA.GUI.FlowUnwrap;

% Save data to the SETs it came from.

for flownumber = 1:length(gui.workingcopy.flows)
  setnumber = gui.workingcopy.SETnumbers(flownumber);
  
  currentflow = gui.workingcopy.flows{flownumber};
  currentflow_segment = currentflow/2 + 1/2; % Convert to segment units
  
  SET(setnumber).IM = currentflow_segment;
  SET(setnumber).Flow.OriginalVENC = gui.workingcopy.OriginalVENC;
  SET(setnumber).VENC = gui.workingcopy.VENC;
end

magnitudenumber = SET(setnumber).Flow.MagnitudeNo;
SET(magnitudenumber).VENC = gui.workingcopy.VENC;
SET(magnitudenumber).Flow.OriginalVENC = gui.workingcopy.OriginalVENC;

closewindow;
%drawfunctions('drawimageno');
viewfunctions('setview'); %drawfunctions('drawall')

%---------------------------------
function cancelpushbutton_Callback 
%---------------------------------
% Cancel pushbutton: Exit window, save nothing.
%

closewindow;

%-------------------------------------
function zoompushbutton_Callback(inout)
%--------------------------------------
% Zoom pushbuttons. Change zoom level, update GUI.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

switch inout
  case 'in'
    gui.zoomlevel = gui.zoomlevel + 1;
  case 'out'
    gui.zoomlevel = gui.zoomlevel - 1;
end

% Enforce max/min zoom levels
if gui.zoomlevel < 0
  gui.zoomlevel = 0;
end

if gui.zoomlevel > gui.maxzoomlevel % arbitrary
  gui.zoomlevel = gui.maxzoomlevel;
end

updatepanimagesliders;
updatevelocityimagezoomandpan;
refocus;

%-------------------------------
function panimageslider_Callback 
%--------------------------------
% Callback for the pan sliders.
%

updatevelocityimagezoomandpan;
refocus;


%------------------------
function playcurrentslice
%-------------------------
% Plays the current slice
%

global DATA
gui = DATA.GUI.FlowUnwrap;

if ~gui.playlock
  gui.playlock = true;
  oldtimeframe = gui.currenttimeframe;

  for framenumber = 1:gui.totaltimeframes
    gui.currenttimeframe = framenumber;
    updatetemporalplotaxes;
    updatevelocityimage;
    drawnow;
  end

  gui.currenttimeframe = oldtimeframe;
  updatetemporalplotaxes;
  updatevelocityimage;

  gui.playlock = false;
end

%----------------------------------
function unwrapcurrentpixel(updown)
%----------------------------------
% Unwrap current pixel in specified direction
%

global DATA
gui = DATA.GUI.FlowUnwrap;

switch updown
  case 'up'
    sgn = +1;
  case 'down'
    sgn = -1;
  otherwise
    myfailed('unwrapcurrentpixel must be called with either ''up'' or ''down'' as its only argument.')
    return
end

expectedjumpsize_orgvenc = gui.expectedjumpsize; % in OriginalVENC
expectedjumpsize = expectedjumpsize_orgvenc * ...
      gui.workingcopy.OriginalVENC / gui.workingcopy.VENC;

currentflow = gui.workingcopy.flows{gui.currentflowno};

x = gui.currentimagepoint(1);
y = gui.currentimagepoint(2);
t = gui.currenttimeframe;
sl = gui.currentslice;

currentflow(x,y,t,sl) = currentflow(x,y,t,sl) + sgn*expectedjumpsize;

gui.workingcopy.flows{gui.currentflowno} = currentflow;

%----------------------------------
function movecurrentpixel(direction)
%--------------------------------
% Moves the current pixel in requested direction, but stops at boundaries.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

x = gui.currentimagepoint(1);
y = gui.currentimagepoint(2);

switch direction
  case 'up'
    x = x - 1;
  case 'down'
    x = x + 1;
  case 'left'
    y = y - 1;
  case 'right'
    y = y + 1;
end

% Bounds check
sizes = size(gui.workingcopy.flows{1});
xsize = sizes(1);
ysize = sizes(2);

x = max(1, x);
x = min(xsize, x);

y = max(1, y);
y = min(ysize, y);

gui.currentimagepoint = [x,y];

%----------------------------
function tms = getcurrenttime
%----------------------------
% Returns the current time in ms.
%

global DATA;
gui = DATA.GUI.FlowUnwrap;

tms = 1000 * gui.TIncr * (gui.currenttimeframe-1);

%-------------------
function closewindow
%-------------------
% Everything that is common to the 'Apply and Exit' and 'Cancel' buttons.
%

global DATA

myworkoff(DATA.fig);

try
  DATA.GUI.FlowUnwrap = close(DATA.GUI.FlowUnwrap);
catch %#ok<CTCH>
  close(gcf);
end

%------------------------------------------
function flows = checkflowdata(NO_to_check)
%-----------------------------------------
% Check what flow data, if any, exists and return an array of stack
% numbers.
% 

global SET

flows=[];
try
  flowstruct = SET(NO_to_check).Flow;
catch ex %#ok<NASGU> 
  % No flow struct. Return empty list.
 % flows = [];
  return;
end

if isempty(flowstruct)
  myfailed('No Flow data found.');
  return;
end
flows = [flowstruct.PhaseNo flowstruct.PhaseX flowstruct.PhaseY];

%---------------------------------------
function slice = getcurrentvelocityslice
%---------------------------------------
% Returns the current velocity slice
%

global DATA;
gui = DATA.GUI.FlowUnwrap;

currentflow = gui.workingcopy.flows{gui.currentflowno};

slice = currentflow(:, :, gui.currenttimeframe, gui.currentslice);
slice = squeeze(slice);

%-------------------------------------------------
function pixel_over_time = getcurrentpixelovertime
%------------------------------------------------
% Returns the currently selected pixel over time, in units of the
% N.B: ORIGINAL VENC.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

x = gui.currentimagepoint(1);
y = gui.currentimagepoint(2);

currentflow = gui.workingcopy.flows{gui.currentflowno};
currentsliceno = gui.currentslice;

pixel_over_time = currentflow(x, y, :, currentsliceno);
pixel_over_time = squeeze(pixel_over_time);

pixel_over_time = pixel_over_time * gui.workingcopy.VENC / gui.workingcopy.OriginalVENC;

%-------------------------------------
function [X, Y] = getoverlayplotcoords
%-------------------------------------
% Computes the coordinate vectors for the overlay plot
%
% Corners:  (x-0.5, y+0.5)
%           (x+0.5, y+0.5)
%           (x+0.5, y-0.5)
%           (x-0.5, y-0.5)
%           (1st again)

global DATA
gui = DATA.GUI.FlowUnwrap;
x = gui.currentimagepoint(2);
y = gui.currentimagepoint(1);

X = [x-0.5, x+0.5, x+0.5, x-0.5, x-0.5];
Y = [y+0.5, y+0.5, y-0.5, y-0.5, y+0.5];

%-----------------------------------------
function [X, Y] = gettimeframemarkercoords
%-----------------------------------------
% Computes the timeframe marker (vertical line in temporal plot) 
% coordinates.
%

x = getcurrenttime;

X = [x x];
Y = [-100 100];

%---------------
function refocus
%----------------
% Gives back focus to slider4 - which has keypressedfcn set to the callback
% handling all keyboard shortcuts. Of course we would like to give focus to
% the figure itself and let the figure handle the keypresses, but matlab won't 
% let us do that. Yes, this is a hack to get around a matlab limitation.
%
% slider4 is placed at (-100, -100), so it shouldn't cause any
% trouble.
%

global DATA
gui = DATA.GUI.FlowUnwrap;

% If this is done to gui.fig, it creates a pushbutton in the lower left
% corner of the figure instead of giving focus to the figure (?!).
uicontrol(gui.handles.slider4);

%-----------------------------------------------------------------
function unwrapped = unwraptimeseries(wrapped, jumpsize, varargin)
%-----------------------------------------------------------------
% function unwrapped = unwraptimeseries(wrapped, jumpsize, jumpfraction,
%                                       debug)
%
% Unwraps a timeseries of samples.
%
% unwrapped:  The unwrapped phase
%
% wrapped:      The wrapped phase
% jumpsize:     Size of the expected jumps in the data
% jumpfraction: Fraction of jumpsize to consider a jump. Default is 0.5.
% debug:        Perform sanity checks on input. Default is true, but can be
%               set to false for speed.
%
% Example usage, when SET(NO) is a phase stack:
%
% SET(NO).IM(x,y,:,slice) = ...
%       unwraptimeseries(SET(NO).IM(x,y,:,slice), 1, 0.5)
% 


if isempty(varargin)
  jumpfraction = 0.5;
end

if length(varargin) >= 1
  if isempty(varargin{1})
    jumpfraction = 0.5;
  else
    jumpfraction = varargin{1};
  end
end

if length(varargin) == 2
  dochecks = varargin{2};
else
  dochecks = true;
end

wrapped = wrapped(:);

if dochecks
  assert(length(varargin) <= 2, 'Too many input arguments.')
  assert(islogical(dochecks), 'debug must be either true or false')
  
  fracmsg = 'Invalid jump fraction, must have 0 < jumpfraction < 1.';
  assert(isnumeric(jumpfraction), fracmsg);
  assert(length(jumpfraction) == 1, fracmsg);
  assert(jumpfraction > 0, fracmsg);
  assert(jumpfraction < 1, fracmsg);

  jumpsizemsg = 'Invalid jump size, must be > 0';
  assert(isnumeric(jumpsize), jumpsizemsg);
  assert(length(jumpsize) == 1, jumpsizemsg);
  assert(jumpsize > 0, jumpsizemsg);

  sizemsg = 'wrapped phase must be a single timeseries.';
  assert(ndims(wrapped) == 2, sizemsg);
  assert(sum(size(wrapped) > 1) == 1, sizemsg);
end

% OK, do the unwrap.
limit = jumpsize * jumpfraction;
diffwrapped = conv2(wrapped, [1;-1], 'valid');
ups = find(diffwrapped > (limit));
downs = find(diffwrapped < (-limit));
unwrapped = wrapped;

if (length(ups)==1) && (length(downs)==1) %Exactly one up and one down
  if ups>downs
    unwrapped((downs+1):(ups)) = (unwrapped((downs+1):(ups)) + jumpsize);
  end

  if ups<downs
    unwrapped((ups+1):(downs)) = (unwrapped((ups+1):(downs)) - jumpsize);
  end
end



