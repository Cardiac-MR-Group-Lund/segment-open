function myinfostruct(s,tit,width,fighandle,illustrationfile,illustrationsize)
%Create input dialog Medviso style. Edited version of myinputstruct
%
%2024
%
%String will shown in a scrollable box.
%Char will be shown as static text.
%
%MYINFOSTRUCT(S, [TITLE],width,[FIGHANDLE],[ILLUSTRATIONFILE],[ILLUSTRATIONSIZE])
%- S = Struct, see below for example.
%- TIT = Title string.
%- FIGHANDLE = optional handle to figure, when supplied used to position input
%- ILLUSTRATIONFILE = optional filename that will be showed as an illustration
%- ILLUSTRATIONSIZE = optional on how big the illustration is, scale 0..1
%
%Note that the function does not call translation and you need to do this
%manually before calling the function.
%
%Example:
%s = []; %reset
%analysedstr = dprintf('This file was automatically loaded and segmented. Please verify the segmentations.');
%s(n).Text = analysedstr;
%n = n + 1;
%s(n).Text = warningstring;
%height = min(amountwarnings, 16);
%s(n).Height = height; %Height of string box

%myinfostruct(s,'Info struct example');

%#ok<*GVMIS> 

global DATA 

if nargin<1
  error('Invalid input.');
end

%Text that should be shown
if ~isfield(s,'Text')
  for loop = 1:length(s)
    s(loop).Text = '';
  end
end

%Max
if ~isfield(s,'Max')
  for loop = 1:length(s)
    s(loop).Max = []; %for now, fixed later
  end
end

%Min
if ~isfield(s,'Min')
  for loop = 1:length(s)
    s(loop).Min = []; %for now, fixed later
  end
end

%Height
if ~isfield(s,'Height')
  s(1).Height = []; %This sets all others to empty as well. 
end

if nargin<2
  tit = '';
end

if nargin<3
  width = 10; %characters
end

if nargin < 5
  fighandle = [];
end

if nargin < 6
  illustrationfile = '';
end

if nargin < 7
  illustrationsize = 1;
end

%Call creation function
[fig, s] = createit(s,dprintf('%s',tit),width,fighandle,illustrationfile,illustrationsize); %There is a uiwait so it returns when user clicked ok

if ~isempty(DATA)
  testmode = DATA.Testing;
else
  testmode = false;
end
if ~ishandle(fig)
  figstatus = 'Cancelled';
else
  figstatus = fig.UserData.Status;
end

if (~ishandle(fig) || isequal(figstatus,'Cancelled')) && ~testmode
  ok = false;
  %Copy default and setup
else
  s = fig.UserData.S; %Ensure get the update one 
end

%Close the figure
if ishandle(fig)
  delete(fig);
end

%--------------------------
function im = fixcolors(im)
%--------------------------
%Recomputes from black and white to correct background details
%Im is an rgb image both in and out.

global DATA

%recompute, assume black and white
bc = DATA.GUISettings.BackgroundColor;
fg = DATA.GUISettings.ForegroundColor;

%Recompute
im = single(im(:,:,1))/255;

%Fix colors
imr = fg(1)*(1-im)+bc(1)*im;
img = fg(2)*(1-im)+bc(2)*im;
imb = fg(3)*(1-im)+bc(3)*im;
im = cat(3,imr,img,imb);

%----------------------------------------------------------------------------------------------------
function [fig, s] = createit(s,tit,minchars,fighandle,illustrationfile,illustrationsize)
%----------------------------------------------------------------------------------------------------
%Creates the figure and populates it

global DATA

%get color settings
if  ~isempty(DATA)
  bc = DATA.GUISettings.BackgroundColor;
  fc = DATA.GUISettings.ForegroundColor;
else
  % default values
  bc = [0.9400 0.9400 0.9400];
  fc = [0 0 0];
end

UserData = [];
UserData.Status = 'Cancelled';
UserData.S = s;
UserData.Outs = [];
UserData.KeyPressed = false;

fig = figure('Name',tit,...
  'Color',bc,...
  'UserData',UserData,...
  'MenuBar','none',...
  'NumberTitle','off'); %set userData to cancelled, if not then change later
  %'WindowStyle','modal',...% ensures that figure reacts to uiwait/uiresume

%Set up Segment icon
setupicon(fig);

%Count number of elements
n = length(s);

%margins
topmarg = 20;
bottommarg = 20;
ymarg = 15; %distance between lines

leftmarg = 10;
rightmarg = 10;
xmarg = 10; %distance between columns

%constants
fontsize = 13;
fontf = 8; %12; %width of character in pixels
rowheight = 24;

%Loop over all to compute width and height
boxwidth = zeros(1,n+1);
height = zeros(1,n+1);
type = cell(1,n+1);
defaultv = zeros(1,n+1);

for loop = 1:n
  switch class(s(loop).Text)
    case 'char'
      %Normal
      boxwidth(loop) = max(length(s(loop).Text),minchars)*fontf;
      height(loop) = rowheight;
      type{loop} = 'char';
    case 'string'
      boxwidth(loop) = max(length(s(loop).Text),minchars)*fontf;
      if isempty(s(loop).Height)
        height(loop) = 5*rowheight;
      else
        height(loop) = s(loop).Height*rowheight;
      end
      type{loop} = 'string';
  end
end

%Add button
boxwidth(n+1) = 0;
height(n+1) = 20;

maxboxwidth = max(boxwidth(1:n))+xmarg;
cumheight = [0 cumsum(height+ymarg)];

%Required width and height
reqw = leftmarg+xmarg+maxboxwidth+rightmarg;
reqh = topmarg+sum(height)+n*ymarg+bottommarg;

%Load the illustration
if ~isempty(illustrationfile)
  try
    illustration = imread(illustrationfile);
  catch
    logdisp(sprintf('Could not read illustration %s',illustrationfile));
    illustration = [];
  end

  if ~isempty(illustration)
    illustration = fixcolors(illustration);
  end
else
  illustration = [];
end

%Add illustration size
oldreqw = reqw;
if ~isempty(illustration)

  rows = size(illustration,1);
  cols = size(illustration,2);

  addwidth = (cols/rows)*reqh*illustrationsize; %illustrationsize is 0..1
  reqw = reqw + addwidth; 
  reqw = round(reqw);
end

% set possible height and width, so always the whole figure is visible
[px, py, reqw, reqh] = helperfunctions('getfigposition',reqw,reqh,fighandle);

%set size
fig.Position = [px py reqw reqh]; 

%---  Add axes
if ~isempty(illustration)
  ax = axes(fig);

  ax.Units = 'pixels';
  axesheight = (reqh-20)*illustrationsize;
  ax.Position = [oldreqw+20 reqh-axesheight-10 addwidth-40 axesheight];
  
  image(ax,illustration);
  axis(ax,'off','image');
end

%--- Add fields
numfields = length(s);
u = zeros(1,numfields);
for loop = 1:numfields
  
  %Compute y position
  posy = reqh-(topmarg+cumheight(loop))-bottommarg-height(loop)+rowheight;
  posx = leftmarg+xmarg; 
  
  shouldfocus = false; %EiH: Default (I think)

  %Fix stri
  switch type{loop}     
    case {'char'}
      stri = s(loop).Text;
      u(loop) = uicontrol(...
        'parent',fig,...
        'Style','text',...
        ...'UserData',s(loop).Field,...
        ...'Tag',s(loop).Field,...
        'TooltipString','',...
        'Position',[posx posy maxboxwidth height(loop)],...
        'HorizontalAlignment','left',...
        'String',stri,...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize,...
        'Enable','on');
     case 'string'
      stri = s(loop).Text;
      u(loop) = uicontrol(...
        'parent',fig,...
        'Style','edit',...
        ...'UserData',s(loop).Field,...
        ...'Tag',s(loop).Field,...
        'TooltipString','',...
        'Position',[posx posy maxboxwidth height(loop)],...
        'HorizontalAlignment','left',...
        'Max',3,'Min',1,...
        'String',stri,...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize,...
        'Enable','inactive');
%     case 'cell'
%       if (numfields == 1) && isequal(s(loop).Max,s(loop).Min) % if max == min, then it is single choice
%         keyfcn = @keyhelper;
%         shouldfocus = true;
%       else
%         keyfcn = '';
%         shouldfocus = false;
%       end
%       if ~isempty(s(loop).Height) && s(loop).Height == 1
%         % maked dropdown
%         uicontrolstyle = 'popupmenu';
%       else
%         uicontrolstyle = 'listbox';
%       end
%       fig.KeyPressFcn = keyfcn;
%       stri = s(loop).Text;
%       u(loop) = uicontrol(...
%         'parent',fig,...
%         'Style',uicontrolstyle,...
%         'UserData',s(loop).Field,...
%         'Tag',s(loop).Field,...
%         'TooltipString','',...
%         'Position',[posx posy maxboxwidth height(loop)],...
%         'HorizontalAlignment','left',...
%         'String',stri,...
%         'Value',s(loop).Value,...
%         'Max',s(loop).Max,...
%         'Min',s(loop).Min,...        
%         'ForegroundColor',fc,...
%         'BackgroundColor',bc,...
%         'Fontsize',fontsize,...
%         'Enable','on',...
%         'KeyPressFcn',keyfcn);
  end

end  %for loop
  
posx = leftmarg;
posy = reqh-(topmarg+cumheight(loop+1))-bottommarg;
  
%Add ok button
%[reqw-50-rightm posy-7 50 height(n+1)]
uicontrol(fig,...
  'style','pushbutton',...
  'Position',[reqw-rightmarg-50 posy 40 height(loop+1)],...
  'Callback', @(u,v)ok_Callback(u,v),...
  'String',dprintf('OK'),...
  'ForegroundColor',fc,...
  'BackgroundColor',bc,...
  'Fontsize',fontsize);

if ~isempty(DATA) &&  DATA.Testing
  v = popfrombuffer('KeyStroke');
  if isempty(v)
    msgstr = 'KeyStroke buffer is empty.';
    myfailed(msgstr);
    return;
  end
  if strcmpi(v,'OK')    
    ind = popfrombuffer('ListboxValue');
    if ~isempty(ind)
      % find listbox
      lb = findobj(fig.Children,'style','listbox');
      if ~isempty(lb)
        lb.Value = ind;
      end
    end
    
    if ~isempty(DATA.Buffer.EditString)
      % find editboxes
      editboxes = findobj(fig.Children,'style','edit');
      if ~isempty(editboxes)
        for loop = 1:numel(editboxes)
          editboxes(loop).String = popfrombuffer('EditString');
        end
      end
    end
    s(1).ok = true;
  else
    s(1).ok = false;
  end
else
  if shouldfocus
    % shift focus to the listbox
%     hlistbox = findobj(fig.Children,'Style','listbox');
%     uicontrol(hlistbox);
  end

  %activate uiwait to wait for user selection until user clicks cancel or ok.
  uiwait(fig);
end

%------------------------
function ok_Callback(h,~)
%------------------------
%Callback function, just resume

uiresume(h.Parent);
UserData = h.Parent.UserData;
UserData.Status = 'Ok';
h.Parent.UserData = UserData;

%------------------------------
function stri = mydouble2str(d)
%------------------------------
%Convert numeric to string

stri = strtrim(sprintf('%11.4g',d));