function varargout = plugin_calibrate(fcn,varargin)
% Plugin to calibrate image resolutions.
%
% Main plugin function. Acts as switchyard, also
% reports dependencys and register menu items.

%Written by Jonatan Wulcan. Slightly commented by Einar Heiberg

global DATA

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end;

switch fcn
  case 'getname'
    varargout = cell(1,1);

    %Register callbacks
    uimenu(varargin{1},'Label','Calibrate resolution','Callback','plugin_calibrate(''opengui_Callback'')');
    set(varargin{1},'Callback','');

    %Return title
    varargout{1} = 'Calibrate';

    % Init DATA.GUI.plugin_calibrate
    DATA.GUI.PluginCalibrate = [];
  case 'getdependencies'
    %Here: List all depending files. This is required if your plugin should
    %be possible to compile to the stand-alone version of Segment.
    varargout = cell(1,4);

    %M-files
    varargout{1} = {};

    %Fig-files
    varargout{2} = {'plugin_calibrate.fig'};

    %Mat-files
    varargout{3} = {};

    %Mex-files
    varargout{4} = {};

  otherwise
    macro_helper(fcn,varargin{:});
    [varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard
end;

%------------------------
function opengui_Callback %#ok<DEFNU>
%------------------------
% Callback that is executed when user press the button in the menu.
% Opens the calibrate gui.

global DATA SET NO

% Check if there is data loaded
if not(DATA.DataLoaded)
  myfailed('No data loaded.');
  return;
end;

% Check if gui is already open
if not(isempty(DATA.GUI.PluginCalibrate))
  return
end

% Create gui
DATA.GUI.PluginCalibrate = mygui('plugin_calibrate.fig');
gui = DATA.GUI.PluginCalibrate;
axis(gui.handles.imageaxes,'off');
axis(gui.handles.imageaxes,'equal');

% Create polygon
xl = SET(NO).XSize*0.3;
xh = SET(NO).XSize*0.7;
yl = SET(NO).YSize*0.3;
yh = SET(NO).YSize*0.7;
DATA.GUI.PluginCalibratePoly = [yl xl; yl xh; yh xh; yh xl];
DATA.GUI.PluginCalibrateActive = 1;

draw();

%--------------------------
function guidelete_Callback %#ok<DEFNU>
%--------------------------
% Executes when the calibrate gui is closed. Cleans up some
% variables.

global DATA

DATA.GUI.PluginCalibratePoly = [];
DATA.GUI.PluginCalibrate = [];

%--------------------------
function guicancel_Callback %#ok<DEFNU>
%--------------------------
% Executes when the user press cancel in the gui.

global DATA

close(DATA.GUI.PluginCalibrate);

%----------------------
function guiok_Callback %#ok<DEFNU>
%----------------------
% Executes when user press ok in the gui.
% Calculates the correct resolution and
% sets the SET struct accordingly. Then
% redraws the images.

global DATA SET NO

correct_size = str2double(mygetedit(DATA.GUI.PluginCalibrate.handles.sizeedit));
if isnan(correct_size)
  mymsgbox('Size contains illegal characters', 'Error', DATA.GUI.PluginCalibrate);
  return;
end

pos = DATA.GUI.PluginCalibratePoly;
area = polyarea(pos(:, 1), pos(:, 2));
SET(NO).ResolutionX = sqrt(correct_size/area);
SET(NO).ResolutionY = sqrt(correct_size/area);
close(DATA.GUI.PluginCalibrate);

%Update graphically in Segment to get grid around image panels correct.
drawfunctions('drawimageno');

%------------
function draw
%------------
% Redraws the images axes in the gui. Called when the polygon has
% been moved.

global DATA SET NO

img = single([]);
img(:, :, 1) = SET(NO).IM(:, :, 1, 1);
img(:, :, 2) = SET(NO).IM(:, :, 1, 1);
img(:, :, 3) = SET(NO).IM(:, :, 1, 1);

image(img, ...
  'parent', DATA.GUI.PluginCalibrate.handles.imageaxes, ...
  'ButtonDownFcn', @moveactive_ButtonDown);
axis(DATA.GUI.PluginCalibrate.handles.imageaxes,'off');
axis(DATA.GUI.PluginCalibrate.handles.imageaxes,'equal');
hold on;
plot([DATA.GUI.PluginCalibratePoly(:, 1); ...
  DATA.GUI.PluginCalibratePoly(1, 1)], ...
  [DATA.GUI.PluginCalibratePoly(:, 2); ...
  DATA.GUI.PluginCalibratePoly(1, 2)], ...
  'r-', 'LineWidth', 3);

for n=1:4
  if n == DATA.GUI.PluginCalibrateActive
    color = 'k';
  else
    color = 'r';
  end
  plot(DATA.GUI.PluginCalibratePoly(n, 1), ...
    DATA.GUI.PluginCalibratePoly(n, 2), ...
    [color 'o'], 'MarkerFaceColor', color, 'ButtonDownFcn', ...
    {@selectpoint_Buttondown, n});
end
hold off;

%-----------------------------------------------
function selectpoint_Buttondown(h, event, point) %#ok<INUSL>
%-----------------------------------------------
% Called when user press a vertex in the polygon.
% Sets the Active point accordingly.

global DATA

seltype = get(DATA.GUI.PluginCalibrate.fig, 'selectiontype');

if strcmp(seltype, 'normal')
  DATA.GUI.PluginCalibrateActive = point;
  draw();
end

%---------------------------------------
function moveactive_ButtonDown(h, event) %#ok<INUSD>
%---------------------------------------
% Called when user clicks the image axes.
% Moves the active point to the new location.

global DATA

seltype = get(DATA.GUI.PluginCalibrate.fig,'selectiontype');
[cp(1),cp(2)] = mygetcurrentpoint(DATA.GUI.PluginCalibrate.handles.imageaxes);
if strcmp(seltype, 'normal')
  DATA.GUI.PluginCalibratePoly(DATA.GUI.PluginCalibrateActive, :) = cp;
  draw();
end
