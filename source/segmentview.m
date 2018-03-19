function [varargout] = segmentview(varargin)
%GUI to select and store views

%TODO
%- rename
%- double click to view
%- store hotkeys 
%- hotkeys in keypressed
%- DONE extra GUI to select hotkey

%Written by Einar and Nisse. Viewing class by Nisse.

if nargin==0
  varargin = {'init'};
end;

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init
%------------
%Initialize GUI
global DATA

%This turns off the icon in Segment
%Change this line when icon in Segment...
%set(DATA.Handles.reportlongaxisicon,'state','off');

%Check if no data loaded
if ~DATA.DataLoaded
  myfailed('You need to have data loaded to run this GUI.');
  return;
end;

%Real initialization
DATA.GUI.SegmentView = mygui('segmentview.fig'); %Create object and store in global the variable DATA
gui = DATA.GUI.SegmentView; %make the shorter variable name gui point to the object also

gui.viewsfile = fullfile(getpreferencespath,'.segment_views.mat');
gui.windowsfile = fullfile(getpreferencespath,'.segment_conbri.mat');
gui.hotkey = '';

set(gui.handles.storeuipanel,'visible','off');
set(gui.handles.viewlistbox,'visible','on');

%Call update to fill listbox
update;
updatecontrast;

%--------------
function update
%--------------
%Fill the listbox
global DATA

gui = DATA.GUI.SegmentView;

%Check for existence & load
if exist(gui.viewsfile,'file')
  load(gui.viewsfile,'viewsarray');
else
  viewsarray = [];
end

if isempty(viewsarray)
  views = {};
  set([gui.handles.renamepushbutton ...
    gui.handles.deletepushbutton ...
    gui.handles.viewpushbutton],'enable','off');
else
  views = cell(1,numel(viewsarray));
  for i = 1:numel(viewsarray)
    hk = viewsarray(i).hotkey;
    if hk
      hk = sprintf('(%s)',hk);
    end
    views{i} = sprintf('%s %s',viewsarray(i).name,hk);
  end
  set([gui.handles.renamepushbutton ...
    gui.handles.deletepushbutton ...
    gui.handles.viewpushbutton],'enable','on');
end;

%Update handles
set(gui.handles.viewlistbox,'String',views);

%Store the viewsarray
gui.viewsarray = viewsarray;
DATA.GUI.SegmentView = gui;

%----------------------
function updatecontrast
%----------------------
%Fill the listbox
global DATA

gui = DATA.GUI.SegmentView;

%Check for existence & load
if exist(gui.windowsfile,'file')
  load(gui.windowsfile,'windowsarray');
else
  windowsarray = [];
end

if isempty(windowsarray)
  windows = {};
  set([gui.handles.renamecontrastpushbutton ...
    gui.handles.deletecontrastpushbutton ...
    gui.handles.applypushbutton],'enable','off');
else
  windows = cell(1,numel(windowsarray));
  for i = 1:numel(windowsarray)
    hk = windowsarray(i).hotkey;
    if hk
      hk = sprintf('(%s)',hk);
    end
    windows{i} = sprintf('%s %s',windowsarray(i).name,hk);
  end
  set([gui.handles.renamecontrastpushbutton ...
    gui.handles.deletecontrastpushbutton ...
    gui.handles.applypushbutton],'enable','on');
end;

%Update handles
set(gui.handles.windowlistbox,'String',windows);

%Store the windowsarray
gui.windowsarray = windowsarray;
DATA.GUI.SegmentView = gui;

%---------------------
function view_Callback %#ok<DEFNU>
%---------------------
global DATA

gui = DATA.GUI.SegmentView;

%Extract selected
v = get(gui.handles.viewlistbox,'value');

%Call load
loadviewfrompreset(gui.viewsarray(v));

%-----------------------------
function viewcontrast_Callback
%-----------------------------
global DATA

gui = DATA.GUI.SegmentView;

%Extract selected
v = get(gui.handles.windowlistbox,'value');

%Call load
loadcontrastfrompreset(gui.windowsarray(v));

%----------------------------------
function loadcontrastfrompreset(cv)
%----------------------------------
global SET NO
if isfield(cv.conbri,'Contrast')
  SET(NO).IntensityMapping.Contrast = cv.conbri.Contrast;
  SET(NO).IntensityMapping.Brightness = cv.conbri.Brightness;
else
  [contrast, brightness] = calcfunctions('win2con',...
    cv.conbri.window,cv.conbri.level);
  SET(NO).IntensityMapping.Contrast = contrast;
  SET(NO).IntensityMapping.Brightness = brightness;
end
drawfunctions('drawcontrastimage',NO);

%-----------------------
function delete_Callback %#ok<DEFNU>
%-----------------------
global DATA

gui = DATA.GUI.SegmentView;

%Extract selected
v = get(gui.handles.viewlistbox,'value');
if isempty(v)
  return
end

%Call delete
deleteviewfrompreset(v);

%Update current element
v = get(gui.handles.viewlistbox,'value');
if ~isequal(v,1)  
  set(gui.handles.viewlistbox,'value',v-1);
end;

%Call update to reflect changes
update;

%-------------------------------
function deletecontrast_Callback %#ok<DEFNU>
%-------------------------------
global DATA

gui = DATA.GUI.SegmentView;

%Extract selected
v = get(gui.handles.windowlistbox,'value');
if isempty(v)
  return
end

%Call delete
deletewindowfrompreset(v);

%Update current element
v = get(gui.handles.windowlistbox,'value');
if ~isequal(v,1)  
  set(gui.handles.windowlistbox,'value',v-1);
end;

%Call update to reflect changes
updatecontrast;

%----------------------
function close_Callback  %#ok<DEFNU>
%----------------------
global DATA

try
  DATA.GUI.SegmentView = close(DATA.GUI.SegmentView);
catch %#ok<CTCH>
  DATA.GUI.SegmentView = [];
  close(gcbf);
end

%-----------------------
function rename_Callback %#ok<DEFNU>
%-----------------------
global DATA

gui = DATA.GUI.SegmentView;

set(gui.handles.storeuipanel,'visible','on');
h = [...
  gui.handles.viewlistbox ...
  gui.handles.viewtext ...
  ];
set(h,'visible','off');

v = get(gui.handles.viewlistbox,'value');
set(gui.handles.nameedit,'String',gui.viewsarray(v).name);
set([gui.handles.shiftcheckbox ...
  gui.handles.ctrlcheckbox ...
  gui.handles.altcheckbox],'value',false);
keystri = gui.viewsarray(v).hotkey;
if ~isempty(keystri)
  set(gui.handles.shiftcheckbox,'value',~isempty(regexpi(keystri,'shift')));
  set(gui.handles.ctrlcheckbox,'value',~isempty(regexpi(keystri,'ctrl')));
  set(gui.handles.altcheckbox,'value',~isempty(regexpi(keystri,'alt')));
  f_ix = regexpi(keystri,'f');
  set(gui.handles.hotkeylistbox,'value',str2double(keystri(f_ix(end)+1:end)));
end
set(gui.handles.storepushbutton,'Callback','segmentview(''gorename_Callback'')');
saveupdate_Callback;

%-------------------------------
function renamecontrast_Callback %#ok<DEFNU>
%-------------------------------
global DATA

gui = DATA.GUI.SegmentView;

set(gui.handles.storecontrastuipanel,'visible','on');
h = [...
  gui.handles.windowlistbox ...
  gui.handles.windowtext ...
  ];
set(h,'visible','off');

v = get(gui.handles.windowlistbox,'value');
set(gui.handles.contrastnameedit,'String',gui.windowsarray(v).name);
set(gui.handles.windowedit,'String',sprintf('%0.4f',gui.windowsarray(v).conbri.window));
set(gui.handles.leveledit,'String',sprintf('%0.4f',gui.windowsarray(v).conbri.level));
set([gui.handles.shiftcontrastcheckbox ...
  gui.handles.ctrlcontrastcheckbox ...
  gui.handles.altcontrastcheckbox],'value',false);
keystri = gui.windowsarray(v).hotkey;
if ~isempty(keystri)
  set(gui.handles.shiftcontrastcheckbox,'value',~isempty(regexpi(keystri,'shift')));
  set(gui.handles.ctrlcontrastcheckbox,'value',~isempty(regexpi(keystri,'ctrl')));
  set(gui.handles.altcontrastcheckbox,'value',~isempty(regexpi(keystri,'alt')));
  f_ix = regexpi(keystri,'f');
  set(gui.handles.contrasthotkeylistbox,'value',str2double(keystri(f_ix(end)+1:end)));
end
set(gui.handles.storecontrastpushbutton,'Callback','segmentview(''gorenamecontrast_Callback'')');
savecontrastupdate_Callback;

%---------------------
function save_Callback %#ok<DEFNU>
%---------------------
%Brings upp new GUI until abort is clicked
global DATA

gui = DATA.GUI.SegmentView;

set(gui.handles.storeuipanel,'visible','on');
h = [...
  gui.handles.viewlistbox ...
  gui.handles.viewtext];
set(h,'visible','off');

set(gui.handles.nameedit,'String','Name of View');
set([gui.handles.shiftcheckbox ...
  gui.handles.ctrlcheckbox ...
  gui.handles.altcheckbox],'Value',false);
set(gui.handles.hotkeylistbox,'Value',1);

set(gui.handles.storepushbutton,'Callback','segmentview(''gosave_Callback'')');
saveupdate_Callback;

%-----------------------------
function savecontrast_Callback
%-----------------------------
%Brings upp new GUI until abort is clicked
global DATA

gui = DATA.GUI.SegmentView;

set(gui.handles.storecontrastuipanel,'visible','on');
h = [...
  gui.handles.windowlistbox ...
  gui.handles.windowtext];
set(h,'visible','off');

set(gui.handles.contrastnameedit,'String','Name of Window');
set([gui.handles.shiftcontrastcheckbox ...
  gui.handles.ctrlcontrastcheckbox ...
  gui.handles.altcontrastcheckbox],'Value',false);
set(gui.handles.contrasthotkeylistbox,'Value',1)

%Set window/level values
[window, level] = calcfunctions('con2win');
set(gui.handles.windowedit,'String',sprintf('%0.4f',window));
set(gui.handles.leveledit,'String',sprintf('%0.4f',level));

set(gui.handles.storecontrastpushbutton,'Callback','segmentview(''gosavecontrast_Callback'')');
savecontrastupdate_Callback;

%----------------------
function abort_Callback 
%----------------------
%Brings upp new GUI until abort is clicked
global DATA

gui = DATA.GUI.SegmentView;

set(gui.handles.storeuipanel,'visible','off');
set(gui.handles.viewlistbox,'visible','on');
set(gui.handles.viewtext,'visible','on');

%------------------------------
function abortcontrast_Callback
%------------------------------
%Brings upp new GUI until abort is clicked
global DATA

gui = DATA.GUI.SegmentView;

set(gui.handles.storecontrastuipanel,'visible','off');
set(gui.handles.windowlistbox,'visible','on');
set(gui.handles.windowtext,'visible','on');

%---------------------------
function saveupdate_Callback %#ok<DEFNU>
%---------------------------
global DATA

gui = DATA.GUI.SegmentView;

%Extract name
name = get(gui.handles.nameedit,'string');
isshift = get(gui.handles.shiftcheckbox,'value');
isctrl = get(gui.handles.ctrlcheckbox,'value');
isalt = get(gui.handles.altcheckbox,'value');
v = get(gui.handles.hotkeylistbox,'value');
hotkeystring = get(gui.handles.hotkeylistbox,'String');
hotkeystring = hotkeystring{v};

hotkey = '';

%Please note these need to be in correct order
if isshift
  hotkey = [hotkey 'Shift-'];
end;
if isctrl
  hotkey = [hotkey 'Ctrl-'];
end;
if isalt
  hotkey = [hotkey 'Alt-'];
end;
hotkey = [hotkey hotkeystring];

set(gui.handles.hotkeytext,'string',hotkey);

%-----------------------------------
function savecontrastupdate_Callback 
%-----------------------------------
global DATA

gui = DATA.GUI.SegmentView;

%Extract name
name = get(gui.handles.contrastnameedit,'string');
isshift = get(gui.handles.shiftcontrastcheckbox,'value');
isctrl = get(gui.handles.ctrlcontrastcheckbox,'value');
isalt = get(gui.handles.altcontrastcheckbox,'value');
v = get(gui.handles.contrasthotkeylistbox,'value');
hotkeystring = get(gui.handles.contrasthotkeylistbox,'String');
hotkeystring = hotkeystring{v};

hotkey = '';

%Please note these need to be in correct order
if isshift
  hotkey = [hotkey 'Shift-'];
end;
if isctrl
  hotkey = [hotkey 'Ctrl-'];
end;
if isalt
  hotkey = [hotkey 'Alt-'];
end;
hotkey = [hotkey hotkeystring];

set(gui.handles.contrasthotkeytext,'string',hotkey);

%-----------------------
function gosave_Callback %#ok<DEFNU>
%-----------------------
global DATA

gui = DATA.GUI.SegmentView;

hotkey = get(gui.handles.hotkeytext,'string');
name = get(gui.handles.nameedit,'string');

%Call save
saveviewtopreset(name,hotkey);

%Call update to reflect changes
update;

abort_Callback; %This resets the view

%-------------------------------
function gosavecontrast_Callback %#ok<DEFNU>
%-------------------------------
global DATA

gui = DATA.GUI.SegmentView;

hotkey = get(gui.handles.contrasthotkeytext,'string');
name = get(gui.handles.contrastnameedit,'string');

%Call save
savecontrasttopreset(name,hotkey);

%Call update to reflect changes
updatecontrast;

abortcontrast_Callback; %This resets the view

%-------------------------
function gorename_Callback
%--------------------------
global DATA

gui = DATA.GUI.SegmentView;

hotkey = get(gui.handles.hotkeytext,'string');
name = get(gui.handles.nameedit,'string');

%Rename
v = get(gui.handles.viewlistbox,'value');
[ok,viewsarray] = checkexistence(name,hotkey,'view',v);

if ~ok
  return
end
% if ismember(hotkey,{viewsarray([1:v-1 v+1:end]).hotkey})
%   myfailed('Hotkey is already used. Please select another key combination.')
%   return
% end
viewsarray(v).name = name;
viewsarray(v).hotkey = hotkey;
save(gui.viewsfile,'viewsarray');

%Call update to reflect changes
update;

abort_Callback; %This resets the view

%-------------------------
function gorenamecontrast_Callback
%--------------------------
global DATA

gui = DATA.GUI.SegmentView;

hotkey = get(gui.handles.contrasthotkeytext,'string');
name = get(gui.handles.contrastnameedit,'string');
window = str2double(get(gui.handles.windowedit,'String'));
level = str2double(get(gui.handles.leveledit,'String'));
  
%Rename
v = get(gui.handles.windowlistbox,'value');
[ok,windowsarray] = checkexistence(name,hotkey,'window',v);

if ~ok
  return
end
% if ismember(hotkey,{viewsarray([1:v-1 v+1:end]).hotkey})
%   myfailed('Hotkey is already used. Please select another key combination.')
%   return
% end
windowsarray(v).name = name;
windowsarray(v).hotkey = hotkey;
windowsarray(v).conbri.window = window;
windowsarray(v).conbri.level = level;
save(gui.windowsfile,'windowsarray');

%Call update to reflect changes
updatecontrast;

abortcontrast_Callback; %This resets the view

%-----------------------------------------
function savecontrasttopreset(name,hotkey)
%-----------------------------------------
global DATA

gui = DATA.GUI.SegmentView;

v = get(gui.handles.windowlistbox,'value');
if isempty(v)
  set(gui.handles.windowlistbox,'value',1);
end

%Filename
[ok,windowsarray,saveix] = checkexistence(name,hotkey,'window');

if ~ok
  return
end

try
  window = str2double(get(gui.handles.windowedit,'String'));
  level = str2double(get(gui.handles.leveledit,'String'));
  savestruct = struct('name',name,'hotkey',hotkey,...
      'conbri',struct('window',window,'level',level));
  if isempty(saveix)
    windowsarray = [windowsarray savestruct]; %#ok<NASGU>
  else
    windowsarray(saveix) = savestruct; %#ok<NASGU>
  end
catch me
  mydispexception(me)
  myfailed('Could not save contrast/brightness window')
  return
end
save(gui.windowsfile,'windowsarray');

%-------------------------------------
function saveviewtopreset(name,hotkey)
%-------------------------------------
global DATA

gui = DATA.GUI.SegmentView;

v = get(gui.handles.viewlistbox,'value');
if isempty(v)
  set(gui.handles.viewlistbox,'value',1);
end
[ok,viewsarray,saveix] = checkexistence(name,hotkey,'view');

if ~ok
  return
end

try
  if isempty(saveix)
    viewsarray = [viewsarray presetview(name,hotkey)];
  else
    viewsarray(saveix) = presetview(name,hotkey); %presetview is a class
  end
catch me
  myfailed('Could not save view')
  return
end
save(gui.viewsfile,'viewsarray');

%---------------------------------------------------------------
function [ok,array,saveix] = checkexistence(name,hotkey,mode,ix)
%---------------------------------------------------------------
global DATA

if nargin < 4
  ix = [];
end

gui = DATA.GUI.SegmentView;

%Filename
ok = true;
saveix = [];
emptystruct = struct('hotkey',{});

if exist(gui.viewsfile,'file')
  load(gui.viewsfile,'viewsarray');
else
  viewsarray = emptystruct;
end

if exist(gui.windowsfile,'file')
  load(gui.windowsfile,'windowsarray');
else
  windowsarray = struct('name',{},'hotkey',{},'conbri',{});
end

switch mode
  case 'view'
    array = viewsarray;
    hotkeylist = {windowsarray.hotkey};
  case 'window'
    array = windowsarray;
    hotkeylist = {viewsarray.hotkey};
end

if ~isempty(array)
  trueix = setdiff(1:numel(array),ix);
  hotkeylist = [hotkeylist {array(trueix).hotkey}];
  if ismember(name,{array(trueix).name})
    if nargout < 3 || ~yesno('View with this name already exists. Do you want to overwrite?')
      myfailed('Name already exists, saving aborted');
      ok = false;
      return
    end
    saveix = find(strcmp(name,{array.name}),1);
  end
end

%Check for hotkey
if ~isempty(hotkey) && ismember(hotkey,hotkeylist)
  if ~yesno('Hotkey is already used. Do you want to overwrite?')
    ok = false;
    return
  end
  hkix = find(strcmp(hotkey,{array.hotkey}),1);
  if ~isempty(hkix)
    array(hkix).hotkey = '';
  else
    switch mode
      case 'view'
        hkix = find(strcmp(hotkey,{windowsarray.hotkey}),1);
        windowsarray(hkix).hotkey = '';
        save(gui.windowsfile,'windowsarray');
        updatecontrast;
      case 'window'
        hkix = find(strcmp(hotkey,{viewsarray.hotkey}),1);
        viewsarray(hkix).hotkey = '';
        save(gui.viewsfile,'viewsarray');
        update;
    end
  end
end

%------------------------------
function deleteviewfrompreset(v)
%-------------------------------
%Delete a saved view

global DATA

gui = DATA.GUI.SegmentView;

thisname = gui.viewsarray(v).name;
if yesno(dprintf('Delete view ''%s''?',thisname))
  n = numel(gui.viewsarray);
  gui.viewsarray = gui.viewsarray([1:v-1 v+1:n]);
  viewsarray = gui.viewsarray; %#ok<NASGU>
  save(gui.viewsfile,'viewsarray'); %Store to file
else
  myfailed('Operation cancelled.');
end

%---------------------------------
function deletewindowfrompreset(v)
%---------------------------------
%Delete a saved view

global DATA

gui = DATA.GUI.SegmentView;

thisname = gui.windowsarray(v).name;
if yesno(dprintf('Delete contrast/brightness window ''%s''?',thisname))
  n = numel(gui.windowsarray);
  gui.windowsarray = gui.windowsarray([1:v-1 v+1:n]);
  windowsarray = gui.windowsarray; %#ok<NASGU>
  save(gui.windowsfile,'windowsarray'); %Store to file
else
  myfailed('Operation cancelled.');
end

%------------------------------
function loadviewfrompreset(pv)
%------------------------------
%Code to bring it up!
nos = pv.match;
shape = pv.matrix;
drawfunctions('drawimageview',nos,shape,pv.panelstype);
% DATA.ViewPanelsType = pv.panelstype;
% DATA.ViewPanels = nos;
% DATA.ViewIM = cell(1,numel(nos));
% drawfunctions('drawall',shape(1),shape(2));

%-----------------------
function keypressed(key) %#ok<DEFNU>
%-----------------------
viewsfile = fullfile(getpreferencespath,'.segment_views.mat');
windowsfile = fullfile(getpreferencespath,'.segment_conbri.mat');

arrayloaded = false;
if exist(viewsfile,'file')
  load(viewsfile,'viewsarray');
  arrayloaded = true;
end
if exist(windowsfile,'file')
  load(windowsfile,'windowsarray');
  arrayloaded = true;
end
if ~arrayloaded
  return
end

if exist('viewsarray','var')
  viewkeys = {viewsarray.hotkey};
  keymatch = find(strcmpi(key,viewkeys),1);
  if keymatch
    loadviewfrompreset(viewsarray(keymatch));
    return
  end
end;

if exist('windowsarray','var')
  windowkeys = {windowsarray.hotkey};
  keymatch = find(strcmpi(key,windowkeys),1);
  if keymatch
    loadcontrastfrompreset(windowsarray(keymatch));
    return
  end
end;
