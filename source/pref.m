function varargout = pref(varargin)
%GUI for preferences
% Panel P1, P2 etc are for different tabs. 
%To make a new tab, do a new panel in the fig file and place your
%preferences in the panel. Add a new tab to the tabgroup and add the panel
%to the tab in initgui.

%Fanny Månefjord, Medviso, 2024

%#ok<*GVMIS>

if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-----------------------------------------------
function initgui
%-----------------------------------------------
%Init GUI
global DATA

%--- Create figure
DATA.GUI.Preferences = mygui('preffig.fig');
gui = DATA.GUI.Preferences;
guih = gui.handles;
set(gui.fig, 'Name', dprintf('Preferences'));
%tab group
gui.handles.tgroup = uitabgroup('Parent', gui.fig,'TabLocation', 'left', 'SelectionChangedFcn',@tabgroup_Callback);
% long name of initpanel to prevent tabs to change size when getting bold
gui.handles.tabinit = uitab('Parent', gui.handles.tgroup, 'Title', dprintf('                               '), 'Tag', 'tabinit'); 
gui.handles.tab1 = uitab('Parent', gui.handles.tgroup, 'Title', dprintf('AI AutoMate'));
% gui.handles.tab2 = uitab('Parent', gui.handles.tgroup, 'Title', dprintf('Patient database'));
% gui.handles.tab3 = uitab('Parent', gui.handles.tgroup, 'Title', dprintf('Services'));
%Add new tabs here
% gui.handles.tab4 = uitab('Parent', gui.handles.tgroup, 'Title', dprintf('My new tab'));
tabgroup_Callback(gui.handles.tgroup);

%Place panels into each tab
set(guih.P1,'Parent',gui.handles.tab1)
% set(guih.P2,'Parent',gui.handles.tab2)
% set(guih.P3,'Parent',gui.handles.tab3)

set(guih.initpanel,'Parent',gui.handles.tabinit)

%Reposition each panel to same location as panel 1
% set(guih.P2,'Position',get(guih.P1,'Position'));
% set(guih.P3,'Position',get(guih.P1,'Position'));
% set(guih.panel2,'Position',get(guih.panel1,'Position'));
% set(guih.panel3,'Position',get(guih.panel1,'Position'));
set(guih.title1, 'String', dprintf('AI AutoMate - default values'));
set(guih.savebutton, 'String', dprintf('Save'), 'Callback', 'pref(''save_Callback'')');
set(guih.closebutton, 'String', dprintf('Close'), 'Callback', 'pref(''close_Callback'')');
set(guih.resetbutton, 'String', dprintf('Reset'), 'Callback', 'pref(''reset_Callback'')');
handles =  fieldnames(gui.handles);
for n = 1:numel(handles)
  objname = handles{n};
  set(gui.handles.(objname),'Units','normalized');
end
set(gui.fig,'Units','pixels');
% Set position for panels and tabgroup
set(guih.P1, 'Position', [0 0 1 1]);
set(guih.initpanel, 'Position', [0 0 1 1]);
set(gui.handles.tgroup, 'Position', [0 0 1 1]);

initpanel1;

%-------------------------------------
function tabgroup_Callback(tabgrouphandle,eventdata)
%-------------------------------------
%Callback for tab group
global DATA

gui = DATA.GUI.Preferences;
handles = gui.handles;

% Get the selected tab
if nargin > 1
    selectedTab = eventdata.NewValue;
else %initialisation
  selectedTab = gui.handles.tab1;
  tabgrouphandle.SelectedTab = selectedTab;
end
if contains(selectedTab.Tag, 'tabinit')
  %special case for tab with empty name
  %Set normal text to all tabs
  alltabs = tabgrouphandle.Children;
  for loop = 1:numel(alltabs)
    htmltxt1 = '<html><strong>';
    htmltxt2 = '</strong>';
    % Remove the htmltxt from Title. It is needed to be done in two steps
    renewtitle = regexprep(alltabs(loop).Title, htmltxt1, '');
    cleanedtitle = regexprep(renewtitle, htmltxt2, '');
    set(alltabs(loop),'Title',cleanedtitle);
  end
else
  strainmitt.strainguifunctions('updateselectedtab',tabgrouphandle);
end
% Move buttons based on selected tab
    if selectedTab == handles.tab1
        set(handles.savebutton, 'Parent', handles.P1);
        set(handles.closebutton, 'Parent', handles.P1);
        set(handles.resetbutton, 'Parent', handles.P1);
        set(handles.title1, 'Parent', handles.P1, 'string', dprintf('AI AutoMate - default values'));
%     elseif selectedTab == handles.tab2
%         set(handles.savebutton, 'Parent', handles.P2);
%         set(handles.closebutton, 'Parent', handles.P2);
%         set(handles.resetbutton, 'Parent', handles.P2);
%         set(handles.title1, 'Parent', handles.P2, 'string', dprintf('Patient database'));
%     elseif selectedTab == handles.tab3
%         set(handles.savebutton, 'Parent', handles.P3);
%         set(handles.closebutton, 'Parent', handles.P3);
%         set(handles.resetbutton, 'Parent', handles.P3);
%         set(handles.title1,'Parent', handles.P3, 'string', dprintf('Services'));
    end
    
    % Ensure buttons are visible
    set(handles.savebutton, 'Visible', 'on');
    set(handles.closebutton, 'Visible', 'on');
    set(handles.resetbutton, 'Visible', 'on');

%-------------------------------------
function initpanel1
%-------------------------------------
%init panel 1, AI AutoMate preferences
global DATA 
gui = DATA.GUI.Preferences;
guih = gui.handles;

%Analysis panel
set(guih.P1_analysis, 'Title', dprintf('Analysis'));
set(guih.P1_checkboxLV, 'String', dprintf('LV'), 'Parent', guih.P1_analysis, 'UserData', 'DATA.Pref.AutomateDefault.LV');
set(guih.P1_checkboxRV, 'String', dprintf('RV'), 'UserData', 'DATA.Pref.AutomateDefault.RV');
set(guih.P1_checkboxLGELV, 'String', dprintf('LV in LGE'), 'UserData', 'DATA.Pref.AutomateDefault.LGELV');
set(guih.P1_checkboxSAX, 'String', dprintf('Strain SAX'), 'UserData', 'DATA.Pref.AutomateDefault.StrainSax');
set(guih.P1_checkboxLAX, 'String', dprintf('Strain LAX'), 'UserData', 'DATA.Pref.AutomateDefault.StrainLax');
panelpos = get(guih.P1_analysis,'position');

%Load panel
panelpos_load = get(guih.P1_load, 'Position');
panelpos_load(3) = panelpos(3); %set same width
set(guih.P1_load, 'Title', dprintf('Loading'), 'Position', panelpos_load);
set(guih.P1_checkboxscout, 'String', dprintf('Exclude scout images'), 'UserData', 'DATA.Pref.AutomateDefault.SkipScout');
set(guih.P1_checkboxsec, 'String', dprintf('Exclude secondary captures'), 'UserData', 'DATA.Pref.AutomateDefault.SkipSecCap');

%Path panel
panelpos_path = get(guih.P1_path, 'Position');
panelpos_path(3) = panelpos(3); %set same width
set(guih.P1_path, 'Title', dprintf('Create Segment files from DICOM'), 'Position', panelpos_path);
inputpathstr = sprintf('%s:',dprintf('Input Path'));
outputpathstr = sprintf('%s:',dprintf('Output Path'));
browsestr = dprintf('Browse');
set(guih.P1_stext1, 'String', inputpathstr);
set(guih.P1_stext2, 'String', outputpathstr);
set(guih.P1_text1, 'String', DATA.Pref.AutomateDefault.PathInput, 'UserData', 'DATA.Pref.AutomateDefault.PathInput');
set(guih.P1_text2, 'String', DATA.Pref.AutomateDefault.PathOutput, 'UserData', 'DATA.Pref.AutomateDefault.PathOutput');
set(guih.P1_pushbuttonset1, 'String', browsestr, 'Callback', 'pref(''browse_Callback'', ''p1_set1'')');
set(guih.P1_pushbuttonset2, 'String', browsestr, 'Callback', 'pref(''browse_Callback'', ''p1_set2'')');

panelpos_dbpath = get(guih.P1_databasepath, 'Position');
panelpos_dbpath(3) = panelpos(3); %set same width
set(guih.P1_databasepath, 'Title', dprintf('Import DICOM with AI AutoMate to database'), 'Position', panelpos_dbpath);
set(guih.P1_stext3, 'String', inputpathstr);
set(guih.P1_text3, 'String', DATA.Pref.AutomateDefault.PathDatabase, 'UserData', 'DATA.Pref.AutomateDefault.PathDatabase');
set(guih.P1_pushbuttonset3, 'String', browsestr, 'Callback', 'pref(''browse_Callback'', ''p1_set3'')');

init = true;
looppanels(init);

%-----------------------------------------------
function initfromdatapref(settingcontrols)
%-----------------------------------------------
%Init values in panel from DATA.Pref
global DATA

for setting = settingcontrols(:)'
  userdatastr = get(setting, 'UserData');
  if ~isempty(userdatastr)
    fields = strsplit(userdatastr, '.');
    savedvalue = getfield(DATA, fields{2:end});
    if strcmp(get(setting, 'Style'), 'text')
      set(setting, 'String', savedvalue);
    else
      set(setting, 'Value', savedvalue);
    end
  end
end


%-----------------------------------------------
function close_Callback
%-----------------------------------------------
%Close the GUI
global DATA

try
  close(DATA.GUI.Preferences);
catch
  close(gcbf);
end

%-----------------------------------------------
function save_Callback
%-----------------------------------------------
%Get all current settings (all tabs)
%Save all current settings to DATA.Pref. 
%Save DATA.Pref to mat file
global DATA
init = false;
looppanels(init);
%save to default_preferences.mat
silent = true;
segpref('save_Callback', silent);

%-----------------------------------------------
function reset_Callback
%-----------------------------------------------
%Reset values to DATA.Pref

init = true;
looppanels(init);

%-----------------------------------------------
function looppanels(init)
%-----------------------------------------------
%Loop the panels (tabs) in the GUI. If init, init from DATA.Pref, otherwise update DATA.Pref

global DATA 
gui = DATA.GUI.Preferences;
guih = gui.handles;

panelnames = {'panel1', 'panel2', 'panel3'}; %Add more names when new tabs are added

for ii = 1:length(panelnames)
  currentpanelname = panelnames{ii}; %different tabs

  if isfield(guih, currentpanelname) % Check if the field exists
    currentpanelhandles = guih.(currentpanelname);
    panelchildren = currentpanelhandles.Children;
    for childind = 1 : length(panelchildren)
      subpanel = panelchildren(childind);
      children = subpanel.Children;
      if init
        initfromdatapref(children)
      else
        updatedatapref(children);
      end
    end
  end
end


%-----------------------------------------------
function updatedatapref(settingcontrols)
%-----------------------------------------------
%Update DATA.Pref to choices in GUI. 

global DATA

for setting = settingcontrols(:)'
  userdatastr = get(setting, 'UserData');
  if ~isempty(userdatastr)
    fields = strsplit(userdatastr, '.');
    if strcmp(get(setting, 'Style'), 'text')
      value = get(setting, 'String'); %for text boxes, get the string
    else
      value = get(setting, 'Value');
    end
    DATA = setfield(DATA, fields{2:end}, value);
  end
end


%-----------------------------------------------
function getdatapref(settingcontrols)
%-----------------------------------------------
%Update choices in GUI from DATA.Pref

global DATA

for setting = settingcontrols(:)'
  userdatastr = get(setting, 'UserData');
  if ~isempty(userdatastr)
    fields = strsplit(userdatastr, '.');
    value = getfield(DATA, fields{2:end}, value);
    set(setting, 'Value', value)
  end
end

%-----------------------
function browse_Callback(button)
%-----------------------
global DATA

gui = DATA.GUI.Preferences;
guih = gui.handles;

switch button
  case 'p1_set1'
    titlestr = dprintf('Select default input path');
    temp = myuigetdir(DATA.SegmentFolder,titlestr);
    if ~(isempty(temp)||isequal(temp,0))
      set(guih.P1_text1, 'String', temp);
    end
  case 'p1_set2'
    titlestr = dprintf('Select default output path');
    temp = myuigetdir(DATA.SegmentFolder,titlestr);
    if ~(isempty(temp)||isequal(temp,0))
      set(guih.P1_text2, 'String', temp);
    end
  case 'p1_set3'
    titlestr = dprintf('Select default input path to the patient database');
    temp = myuigetdir(DATA.SegmentFolder,titlestr);
    if ~(isempty(temp)||isequal(temp,0))
      set(guih.P1_text3, 'String', temp);
    end
end
 
