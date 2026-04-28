function varargout = generalpenfunctions(varargin)
%This contains all callbacks related to general pen objects.
%
% Justine Le Douaron, Medviso, 2023

%#ok<*GVMIS>

if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%----------------------------------------------
function setlabel_Callback
%----------------------------------------------
global DATA SET NO

labelsdictionnary = DATA.GeneralPenSettings.getstandardlabels;
labels = labelsdictionnary.values;
labels = cat(2,'Define a custom label...',labels); %add an entry for user-defined label

n = 1;
field = 'GeneralPenLabels';
inputstruct(n).Field = field;
inputstruct(n).Label = dprintf('Select a label');
inputstruct(n).Default = labels;
[outs,ok] = myinputstruct(inputstruct,dprintf('Set object''s label'),10);
if ok
  userchoice = outs.(field);
else
  %user cancelled
  return
end

if userchoice == 1
  defaultlabel = 'Object';
  label = myinputdlg({'Enter a label'},'',1,{defaultlabel});
  if isempty(label)
    label = defaultlabel;
  else
    label = label{1};
  end
end

if ((userchoice > 1) || (userchoice < 1)) && userchoice > 0
  label = labels{userchoice};
end

objind = DATA.GeneralPenSettings.getcurrentobject;
SET(NO).GeneralPenObjects(objind).setlabel(label);
SET(NO).GeneralPenObjects(objind).updateicon;
indentobjecticon(objind);

%Graphical update
drawfunctions('drawcontoursgeneralpen',DATA.CurrentPanel);
addicontotoolbar(SET(NO).GeneralPenObjects(objind).Icon);

%----------------------------------------------
function setcolor_Callback
%----------------------------------------------
global DATA SET NO

colorssdictionary = DATA.GeneralPenSettings.getgeneralpencolors;
colorkeys = colorssdictionary.keys;

n = 1;
field = 'GeneralPenColors';
inputstruct(n).Field = field;
inputstruct(n).Label = dprintf('Select a color');
inputstruct(n).Default = colorkeys;
[outs,ok] = myinputstruct(inputstruct,dprintf('Set General Pen color'),10);
if ~ok %user cancelled
  return
end

userchoice = outs.(field);
colorvalue = colorssdictionary(colorkeys{userchoice});
objind = DATA.GeneralPenSettings.getcurrentobject;
SET(NO).GeneralPenObjects(objind).setcolor(colorvalue);
SET(NO).GeneralPenObjects(objind).updateicon;
indentobjecticon(objind);

%Graphical update
drawfunctions('drawcontoursgeneralpen',DATA.CurrentPanel);
addicontotoolbar(SET(NO).GeneralPenObjects(objind).Icon);

%---------------------------------------------------------
function deleteobject_Callback
%---------------------------------------------------------
%Delete selected General Pen object
global DATA SET NO

if isempty(DATA.GeneralPenSettings.getcurrentobject)
  myfailed('No object selected.',DATA.GUI.Segment);
  return
else
  objecttodeleteind = DATA.GeneralPenSettings.getcurrentobject;
end

%Remove object from SET and from objects' toolbar
removeiconfromobjectstoolbar(SET(NO).GeneralPenObjects(objecttodeleteind).Icon);
SET(NO).GeneralPenObjects(objecttodeleteind) = [];

%Update the currently selected object
if DATA.GeneralPenSettings.getnumobjects > 0
  newobjectind = DATA.GeneralPenSettings.getnumobjects;
  indentobjecticon(newobjectind);
else
  newobjectind = [];
end
DATA.GeneralPenSettings.setcurrentobject(newobjectind);

%Graphical update
set(DATA.Handles.generalpentext(DATA.CurrentPanel,objecttodeleteind),'Position',[nan nan]);
drawfunctions('clearcontoursgeneralpen',DATA.CurrentPanel,objecttodeleteind);

%---------------------------------------------------------
function deleteallobjects_Callback
%---------------------------------------------------------
%Deleted all General Pen objects
global DATA SET NO

%Remove objects from SET and from objects' toolbar
SET(NO).GeneralPenObjects(:) = [];
DATA.GeneralPenSettings.setcurrentobject([]);
clearobjectstoolbar;

%Graphical update
set(DATA.Handles.generalpentext(:,DATA.CurrentPanel),'Position',[nan nan]);
drawfunctions('drawno',NO);

%---------------------------------------------------------
function varargout = createnewobject
%---------------------------------------------------------
%Create a new General Pen object.
global DATA

%Instantiate a new object
newobjectind = DATA.GeneralPenSettings.createnewobject;

%Indent the icon that is associated with the new object
indentobjecticon(newobjectind);

if nargout == 1
  varargout{1} = newobjectind;
end

%---------------------------------------------------------
function createnewobject_Callback
%---------------------------------------------------------
%Create a new General Pen object and perform graphical update of the panel.
global DATA

maxnumobjects = DATA.GeneralPenSettings.getmaxnumobjects;
if DATA.GeneralPenSettings.getnumobjects >= maxnumobjects
  str1 = dprintf('Cannot create a new object.');
  str2 = dprintf('Maximum number of objects allowed:');
  str = sprintf('%s %s %d.',str1,str2,maxnumobjects);
  myfailed(str);
  return
end

%Instantiate a new object and add the new object to the toolbar
createnewobject;

%Graphical update
drawfunctions('unselectcontoursgeneralpen',DATA.CurrentPanel);

%---------------------------------------------------------
function indentobjecticon(objectind)
%---------------------------------------------------------
%Indent the icon associated wiht the current object
global DATA SET NO
executeicon = false;
iconname = SET(NO).GeneralPenObjects(objectind).Icon.name;
indent(DATA.Handles.configiconholder,iconname,executeicon);

%---------------------------------------------------------
function setcurrentobject_Callback
%---------------------------------------------------------
%Update the current object. Called when clicking an object's icon the 
%objects' toolbar.
global DATA SET NO

%Get the clicked icon
clickedicon = DATA.Handles.configiconholder.clickedicon.name;

%Find the corresponding object
for objectind = 1:length(SET(NO).GeneralPenObjects)
  if strcmp(clickedicon,SET(NO).GeneralPenObjects(objectind).Icon.name)
    break
  end
end

%Select the corresponding object
DATA.GeneralPenSettings.setcurrentobject(objectind);

%Graphical update
drawfunctions('drawcontoursgeneralpen',DATA.CurrentPanel);

%---------------------------------------------------------
function addicontotoolbar(icon)
%---------------------------------------------------------
%Update the toolbar with an object's icon
global DATA

generalpeniconcell = DATA.Icons.generalpeniconcell;

numberoficons = numel(generalpeniconcell);
iconind = numberoficons + 1;
for iconloop = 1:numberoficons
  name = generalpeniconcell{iconloop}.name;
  if isequal(deblank(name),icon.name)
    iconind = iconloop;
    break;
  end
end

generalpeniconcell{1,iconind} = icon;
DATA.Icons.generalpeniconcell = generalpeniconcell;
DATA.Handles.configiconholder.add(DATA.Icons.generalpeniconcell);

%---------------------------------------------------------
function clearobjectstoolbar
%---------------------------------------------------------
%Clear the objects' toolbar
global DATA

generalpeniconcell = DATA.Icons.generalpeniconcell;

if isempty(generalpeniconcell)
  return %returning here saves some time if generalpeniconcell is not used, such in Segment 3DPrint
end

iconstodelete = [];
numberoficons = numel(generalpeniconcell);
for iconloop = 1:numberoficons
  name = generalpeniconcell{iconloop}.name;
  if contains(deblank(name),"generalpenobject")
    iconstodelete = [iconstodelete iconloop]; %#ok<AGROW> 
  end
end

generalpeniconcell(iconstodelete) = [];
DATA.Icons.generalpeniconcell = generalpeniconcell;
DATA.Handles.configiconholder.add(DATA.Icons.generalpeniconcell);

%---------------------------------------------------------
function removeiconfromobjectstoolbar(icon)
%---------------------------------------------------------
%Remove an icon from the toolbar
global DATA

generalpeniconcell = DATA.Icons.generalpeniconcell;
numberoficons = numel(generalpeniconcell);
for iconloop = 1:numberoficons
  name = generalpeniconcell{iconloop}.name;
  if isequal(deblank(name),icon.name)
    iconind = iconloop;
    break;
  end
end

generalpeniconcell(iconind) = [];
DATA.Icons.generalpeniconcell = generalpeniconcell;
DATA.Handles.configiconholder.add(DATA.Icons.generalpeniconcell);

%---------------------------------------------------------
function initialiseobjectstoolbar
%---------------------------------------------------------
%Populate the objects' toolbar
global DATA SET NO

%Verify that the current stack has any general pen object
if isempty(SET(NO).GeneralPenObjects)
  return
end

%Get all the objects' icons
generalpeniconcell = DATA.Icons.generalpeniconcell;
numberofpermanenticons = numel(generalpeniconcell);
for objectind = 1:length(SET(NO).GeneralPenObjects)
  icon = SET(NO).GeneralPenObjects(objectind).Icon;
  generalpeniconcell{1,numberofpermanenticons+objectind} = icon;
end
indentobjecticon(objectind);

%Update the toolbar
DATA.Icons.generalpeniconcell = generalpeniconcell;
DATA.Handles.configiconholder.add(DATA.Icons.generalpeniconcell);

%---------------------------------------------------------
function updateobjectstoolbar
%---------------------------------------------------------
%Update the objeects' toolbar, used when switching image stack
global DATA

if strcmp(DATA.CurrentTheme,'generalpen')
  clearobjectstoolbar;
  initialiseobjectstoolbar;
end