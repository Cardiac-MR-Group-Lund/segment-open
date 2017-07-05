function weightingslider(arg)
%Adjust slider for MO threshold

%Einar Heiberg

global DATA

if nargin == 0
  arg = 'init';
end;

switch arg
  case 'autopushbutton_Callback'
    autopushbutton_Callback;
  case 'init'
    %init
    DATA.GUI.WeightingSlider = mygui('weightingslider.fig','blocking'); %Create object and store in global the variable DATA    
    gui = DATA.GUI.WeightingSlider;
     recursekeypressfcn(gui.fig);
     
  case 'slider'
    try
      disp('Weightingslider moved.');
      gui = DATA.GUI.WeightingSlider;

      %Update edit box
      v = get(gui.handles.weightingslider,'value'); %get slider value
      stri = sprintf('%0.2g',v); %make string
      set(gui.handles.weightingedit,'string',stri); %update edit box      
      %Update scar calculation and display
      update;      
    catch me
      mydispexception(me);
    end;
  case 'edit'
    try
      disp('Weightingeditbox.');
      gui = DATA.GUI.WeightingSlider;
      
      stri = get(gui.handles.weightingedit,'String'); %get string
      v = str2double(stri); %parse the string to a number
      
      if ~isnan(v)        
        oldv = v;
        v = min(max(0.001,v),1); %Ensure in range
        
        if ~isequal(oldv,v)
          stri = sprintf('%0.2g',v); %make string
          set(gui.handles.weightingedit,'string',stri); %update edit box if changed to be inside range
        end;
        
        set(gui.handles.weightingslider,'value',v);
        %Update scar calculation and display
        update;
      end;
    catch me
      mydispexception(me);
    end;
  case 'close'
    disp('trying to close weightingslider');
    try
      DATA.GUI.WeightingSlider = close(DATA.GUI.WeightingSlider);
    catch
      delete(gcbf);
      DATA.GUI.WeightingSlider = [];
    end
  case 'up'    
    up;
  case 'down'
    down;
end;

%---------------
function update
%---------------
%Updates infarct info

global SET NO DATA

try
  gui = DATA.GUI.WeightingSlider; %get  handle
  v = get(gui.handles.weightingslider,'value'); %get value

  %Perform update if scar
  if ~isempty(SET(NO).Scar)                       
      
    SET(NO).Scar.MR.WeightingScale = v; %Store in SET struct
    
    viability('viabilitycalc');
    drawfunctions('drawimagepanel',DATA.CurrentPanel);
    
    figure(gui.fig);
  end;
catch me
  mydispexception(me);
end

%-----------------------------
function recursekeypressfcn(h)
%-----------------------------
global DATA

%Start with the current handle.
try
  set(h,'keypressfcn',@DATA.keypressed);
catch %#ok<CTCH>
end;

children = get(h,'children');
for loop=1:length(children)
  recursekeypressfcn(children(loop));
end;

%----------
function up
%----------
changehelper(0.05);

%----------
function down
%----------
changehelper(-0.05);

%-------------------------
function changehelper(d)
%-------------------------
global DATA SET NO

try
  gui = DATA.GUI.WeightingSlider; %get  handle
  v = get(gui.handles.weightingslider,'value'); %get value
catch
  %Likely the slider is not active => return;
  return;
end;

if isempty(SET(NO).Scar)
  return;
end;

v = v + d;
v = min(max(0.001,v),1); %Ensure in range

%Update graphically
stri = sprintf('%0.2g',v); %make string
set(gui.handles.weightingedit,'string',stri); %update edit box if changed to be inside range
set(gui.handles.weightingslider,'value',v);

%Update is setstruct
SET(NO).Scar.MR.WeightingScale = v; %Store in SET struct

%update graphically
viability('viabilitycalc');
drawfunctions('drawimagepanel',DATA.CurrentPanel);

%Make figure visible if not
figure(gui.fig);

%-----------------------------------
function autopushbutton_Callback
%-----------------------------------
global DATA SET NO

gui = DATA.GUI.WeightingSlider; %get  handle

if isempty(SET(NO).Scar)
  return;
end;

ifsize = inputdlg('Infarctsize','Input value');

if isempty(ifsize)
  myfailed('Aborted.');
  return;
end;

ifsize = str2double(ifsize{1});
if isnan(ifsize)
  myfailed('Invalid number.');
  return;
end;

maxw = 1;
minw = 0.001;
w = 0.5*(maxw+minw);

h = waitbar(0,'Please wait.');
for loop = 1:12
  SET(NO).Scar.MR.WeightingScale = w;
  viability('viabilitycalc');
  p =   SET(NO).Scar.Percentage;
  
  if p>ifsize
    maxw = w;
    w = 0.5*(maxw+minw);
  end;
  
  if p<ifsize
    minw = w;
    w = 0.5*(maxw+minw);
  end;

  waitbar(loop/12,h);
end;
close(h);

SET(NO).Scar.MR.WeightingScale = w;

%Update graphically
stri = sprintf('%0.2g',w); %make string
set(gui.handles.weightingedit,'string',stri); %update edit box if changed to be inside range
set(gui.handles.weightingslider,'value',w);

%update graphically
viability('viabilitycalc');
drawfunctions('drawimagepanel',DATA.CurrentPanel);
