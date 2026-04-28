function varargout = generatefilenames(varargin)
%---------------------------------------------
% Generate new names GUI
%Generating file names with a prefix and suffix. Prefix can be chosen,
%random or non-exsistent. Suffix is made by numbers and and start number
%can be chosen, or the suffix can be random or non-existent.
%
%Fanny Månefjord, Medviso, 2023

%#ok<*GVMIS>

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%--------------------------------------------------
function init
%--------------------------------------------------
%Initialize the generate file names GUI 
global DATA
DATA.GUI.NameGenerator = mygui('generatefilenames.fig');
gui = DATA.GUI.NameGenerator;

%Set numbers in pop up menu to 3-9 and let 5 be default value
numbers = string(3:9);
set(gui.handles.prefix_popup,'String',numbers);
set(gui.handles.prefix_popup,'Value',3);
set(gui.handles.suffix_popup,'String',numbers);
set(gui.handles.suffix_popup,'Value',3);

%Random will be default for both prefix and suffix
set(gui.handles.randomradiobutton,'Value',1);
set(gui.handles.randomsuffix_radiobutton,'Value',1);

%Set keypressfunctions for edit boxes
set(gui.handles.prefixedit,'KeyPressFcn','generatefilenames(''writeprefix_Callback'')');
set(gui.handles.suffix_edit,'KeyPressFcn','generatefilenames(''writenumber_Callback'')');

%Set callbacks for pressing the popup menus
set(gui.handles.prefix_popup,'Callback','generatefilenames(''ppopup_Callback'')');
set(gui.handles.suffix_popup,'Callback','generatefilenames(''spopup_Callback'')');

%Set callbacks for pushbuttons
set(gui.handles.pushbuttoncancel,'Callback','generatefilenames(''close_Callback'')');
set(gui.handles.pushbutton_generate,'Callback','generatefilenames(''generate_Callback'')');

%--------------------------------------------------
function ok_Callback
%--------------------------------------------------
global DATA
gui = DATA.GUI.NameGenerator;

%--------------------------------------------------
function close_Callback %called from GUI
%--------------------------------------------------
% close figure
global DATA

try
  DATA.GUI.NameGenerator = close(DATA.GUI.NameGenerator);
catch %#ok<CTCH>
  DATA.GUI.NameGenerator = [];
  delete(gcbf);
end

%-----------------------------------------------
function generate_Callback
%-----------------------------------------------
%Create file names with prefix and suffix

global WL DATA

gui = DATA.GUI.NameGenerator;
reset(RandStream.getGlobalStream,sum(100*clock)); %reset random seed

%find which files to anonymize from worklist
anonymize=zeros(size(WL,1),1);
[anonymize(:)]=[WL{:,1}];
indextoanonymize = find(anonymize == 1);
numfiles = length(indextoanonymize);
if numfiles == 0
  myfailed('No valid data to pseudononymize');
  return;
end

newnames = cell(numfiles ,2);

%Choices for prefix and suffix
randomprefix = get(gui.handles.randomradiobutton,'value');
chooseradiotton = get(gui.handles.choose_radiobutton,'value');
randomsuffix = get(gui.handles.randomsuffix_radiobutton,'value');
choosesuffix = get(gui.handles.choosesuffix_radiobutton,'value');
noprefix = get(gui.handles.noprefix_radiobutton,'value');
nosuffix = get(gui.handles.nosuffix_radiobutton,'value');

%It should not be possible to not have a prefix nor suffix
if noprefix && nosuffix
  myfailed('Please choose at least prefix or suffix');
  return;
end

%---Prefix

if (randomprefix) 
  letters = (get(gui.handles.prefix_popup,'string'));
  index = (get(gui.handles.prefix_popup,'value'));
  nbrletters = str2num(letters{index}); %how long the random prefix should be
  for j = 1:numfiles
    prefix = char(randi([65,90], 1, nbrletters));
    newnames{j,1} = prefix;
  end

elseif (chooseradiotton) %user has chosen a prefix
  prefix = gui.handles.prefixedit.String;
  if isempty(prefix)
    myfailed('Please choose a prefix.');
    return;
  end
  for j = 1:numfiles
  newnames{j,1} = prefix;
  end
end

%---Suffix

if randomsuffix  

  %Number of digits in the popup menu
  digits = (get(gui.handles.suffix_popup,'string'));  
  index = (get(gui.handles.suffix_popup,'value'));
  nbrdigits = str2num(digits{index});
  

  for j=1:numfiles
    str = num2str(randi([0,9], 1, nbrdigits));
    str = str(~isspace(str));
     newnames{j,2}=str;
  end

elseif choosesuffix     %starting point of suffix is chosen 
  start = gui.handles.suffix_edit.String;
  len = length(start);
  if len == 0
    myfailed('Please choose a starting number.');
    return;
  end
  start = str2num(start);
  if isempty(start) || ~isreal(start)
    myfailed('Choose a number');
    return;
  end
  
  ind=start;
  for j = 1:numfiles
    suf = pad(num2str(ind),len,'left','0'); %pad with zeros, e.g. 0001
    newnames{j,2} = suf;
    ind = ind+1;
  end
end

%concatenate prefix and suffix and add to worklist
 for j=1:numfiles
   name = strcat(newnames{j,1},newnames{j,2});
   WL{indextoanonymize(j),2}=name;
 end

anonymization('updateoutput');
anonymization('updatelist');
close_Callback;

%-----------------------------------------------
function writenumber_Callback
%-----------------------------------------------
%Mark radiobutton when writing in the edit field
global DATA

gui = DATA.GUI.NameGenerator;
set(gui.handles.choosesuffix_radiobutton,'Value',1);

%-----------------------------------------------
function writeprefix_Callback
%-----------------------------------------------
%Mark radiobutton when writing in the edit field
global DATA

gui = DATA.GUI.NameGenerator;
set(gui.handles.choose_radiobutton,'Value',1);

%-----------------------------------------------
function ppopup_Callback
%-----------------------------------------------
%Mark radiobutton when clicking in the popup menu
global DATA

gui = DATA.GUI.NameGenerator;
set(gui.handles.randomradiobutton,'Value',1);

%-----------------------------------------------
function spopup_Callback
%-----------------------------------------------
%Mark radiobutton when clicking in the popup menu
global DATA

gui = DATA.GUI.NameGenerator;
set(gui.handles.randomsuffix_radiobutton,'Value',1);

