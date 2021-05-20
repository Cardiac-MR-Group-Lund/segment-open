function [outs,ok] = myinputstruct(s,tit,width)
%Automatically create input dialogs with Medviso style
%This function replaces the old inputstruct.
%
%MYINPUTSTRUCT(S,TITLE,WIDTH)
%- S = Struct, see below for example.
%- TIT = Title string.
%- WIDTH = Width of editbox (optional, default = 10)
%
%Note that the function does not call translation and you need to do this
%manually before calling the function.
%
%Example:
%s = []; %reset
%s(1).Field = 'Name';
%s(1).Label = dprintf('Name'); %Use dprintf only if you want the input to be translated.
%s(1).Default = 'Box';
%s(2).Field = 'H';
%s(2).Label = dprintf('Height [mm]');
%s(2).Default = 10;
%s(2).Limits = [0,10];
%s(3).Field = 'W';
%s(3).Label = dprintf('Width [mm]');
%s(3).Default = 20;
%s(3).Limits = [0,30];
%s(4).Field = 'Visible';
%s(4).Label = dprintf('Visible');
%s(4).Default = true;
%[outs,ok] = myinputstruct(s,'Input box parameters',10); %10 is 10 characters minimum
%
%If user press cancel or close then ok = false;
%
%Optional fields
%- Tooltip, sets Tooltip string

%REMOVED FUNCTIONALITY
%- Limits (default [-inf inf])
%- ValueDisplayFormat (default '%11.4g')
%- LowerLimitInclusive (default 'on')
%- UpperLimitInclusive (default 'on')
%
%Hint to input decimals use any of the int/uint classes

%Einar Heiberg

%left to do:
%- put limits back.
%- listbox support
%- DONE ok button
%- DONE correct colors
%- DONE add translation

outs = [];

if nargin<1
  error('Invalid input.');
end

%Error checks
if ~isfield(s,'Field')
  disp('Required fieldname Field is missing.');
  error('Invalid input.');
end

if ~isfield(s,'Label')
  disp('Required fieldname Label is missing.');
  error('Invalid input.');
end

if ~isfield(s,'Default')
  for loop = 1:length(s)
    s(loop).Default = '';
  end
end

if ~isfield(s,'Tooltip')
  for loop = 1:length(s)
    s(loop).Tooltip = '';
  end
end

% if ~isfield(s,'ValueDisplayFormat')
%   for loop = 1:length(s)
%     s(loop).ValueDisplayFormat = '%11.4g';
%   end
% end
% 
% if ~isfield(s,'Limits')
%   for loop = 1:length(s)
%     s(loop).Limits = [-inf inf];
%   end
% end
% if ~isfield(s,'LowerLimitInclusive')
%   for loop = 1:length(s)
%     s(loop).LowerLimitInclusive = 'on';
%   end
% end
% 
% if ~isfield(s,'UpperLimitInclusive')
%   for loop = 1:length(s)
%     s(loop).UpperLimitInclusive = 'on';
%   end
% end
% 
% %Check validity
% for loop = 1:length(s)
%   if isempty(s(loop).Limits)
%     s(loop).Limits = [-inf inf];
%   end
% end
% 
% for loop = 1:length(s)
%   if isempty(s(loop).LowerLimitInclusive)
%     s(loop).ToolTip = 'on';
%   end
% end
% 
% for loop = 1:length(s)
%   if isempty(s(loop).UpperLimitInclusive)
%     s(loop).ToolTip = 'on';
%   end
% end
% 
% for loop = 1:length(s)
%   if isempty(s(loop).ValueDisplayFormat)
%     s(loop).ValueDisplayFormat = '%11.4g';
%   end
% end

for loop = 1:length(s)
  if isempty(s(loop).Tooltip)
    s(loop).ToolTip = '';
  end
end

for loop = 1:length(s)
  if isempty(s(loop).Label)
    s(loop).Label = s(loop).Field;
  end
end

if nargin<2
  tit = '';
end

if nargin<3
  width = 10; %characters
end

%Call creation function
fig = createit(s,dprintf('%s',tit),width); %There is a uiwait so it returns when user clicked ok

%Get data
if ~ishandle(fig) || isequal(fig.UserData,'Cancelled')
  ok = false;
  %Copy default and setup
  for loop = 1:length(s)
    outs.(s(loop).Field) = s(loop).Default;
  end

else
  [outs,ok] = parse(fig,s);
end

%Close the figure
if ishandle(fig)
  delete(fig);
end

%--------------------------------
function [outs,ok] = parse(fig,s)
%--------------------------------
%Gets the edited results

outs = [];
ok = true;

%Copy default and setup
for loop = 1:length(s)
  outs.(s(loop).Field) = s(loop).Default;
end

c = fig.Children;

errors = false;
message = '';
%Loop over elements to get them
for loop = 1:length(c)
  ud = c(loop).UserData;
  if ~isempty(ud)
    
    %Loop over to find. A bit uggly but who cares, fast enough
    id = 1;
    for sloop = 1:length(s)
      if isequal(s(sloop).Field,ud)
        id = sloop;
      end
    end
    
    %Assign it
    if isnumeric(s(id).Default)
      outs.(ud) = str2double(c(loop).String);
      if isnan(outs.(ud))
        errors = true;
        message = dprintf('Invalid response for %s',s(id).Label);       
      end
    end
    
    if islogical(s(id).Default)
      outs.(ud) = logical(c(loop).Value);
    end
    
    if ischar(s(id).Default)
      outs.(ud) = c(loop).String;
    end
    
  end
end

if errors
  ok = false;
  myfailed(message);
end

%--------------------------------------
function fig = createit(s,tit,minchars)
%--------------------------------------
%Creates the figure and populates it

global DATA

%get color settings
bc = DATA.GUISettings.BackgroundColor;
fc = DATA.GUISettings.ForegroundColor;

fig = figure('Name',tit,...
  'Color',bc,...
  'UserData','Cancelled',...
  'MenuBar','none',...
  'NumberTitle','off'); %set userData to cancelled, if not then change later

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
fontsize = 14;
fontf = 12; %width of character in pixels
rowheight = 25;

%Loop over all to compute width and height
labelwidth = zeros(1,n+1);
boxwidth = zeros(1,n+1);
height = zeros(1,n+1);
type = cell(1,n+1);

for loop = 1:n
  switch class(s(loop).Default)
    case {'single','double','int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
      %Numeric
      labelwidth(loop) = length(s(loop).Label)*fontf;
      boxwidth(loop) = max(length(strtrim(sprintf('%11.4g',s(loop).Default))),minchars)*fontf;
      height(loop) = rowheight;
      type{loop} = 'numeric';
    case 'logical'
      labelwidth(loop) = (2+length(s(loop).Label))*fontf;
      boxwidth(loop) = 0;
      height(loop) = rowheight;     
      type{loop} = 'logical';
    case 'char'
      labelwidth(loop) = length(s(loop).Label)*fontf;
      boxwidth(loop) = max(length(s(loop).Default),minchars)*fontf;
      height(loop) = rowheight;
      type{loop} = 'char';
  end
end

%Add button
labelwidth(n+1) = 300;
boxwidth(n+1) = 0;
height(n+1) = 20;

maxlabelwidth = max(labelwidth(1:n));
maxboxwidth = max(boxwidth(1:n))+xmarg;
cumheight = [0 cumsum(height+ymarg)];

%Required width and height
reqw = leftmarg+maxlabelwidth+xmarg+maxboxwidth+rightmarg;
reqh = topmarg+sum(height)+n*ymarg+bottommarg;

%Get and set size
p = fig.Position;
fig.Position = [p(1) p(2) reqw reqh]; %later center

u = zeros(1,length(s));
for loop = 1:length(s)
  
  %Compute y position
  posy = reqh-(topmarg+cumheight(loop))-bottommarg;
  posx = leftmarg+maxlabelwidth+xmarg;
  
  %Fix stri
  switch type{loop}
    case 'numeric'
      %Numeric
      stri = sprintf('%11.4g',s(loop).Default);
      stri = stri(stri~=' ');
      u(loop) = uicontrol(...
        'parent',fig,...
        'Style','edit',...
        'UserData',s(loop).Field,...
        'TooltipString',sprintf('%s',s(loop).Tooltip),...
        'Position',[posx posy maxboxwidth height(loop)],...
        'HorizontalAlignment','center',...
        'String',stri,...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize);     
    case 'logical'
      stri = s(loop).Label;
      u(loop) = uicontrol(...
        'parent',fig,...
        'Style','checkbox',...
        'UserData',s(loop).Field,...
        'TooltipString',sprintf('%s',s(loop).Tooltip),...
        'Position',[leftmarg posy maxlabelwidth height(loop)],...
        'String',stri,...
        'Value',s(loop).Default,...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize);
    case 'char'
      stri = s(loop).Default;
      u(loop) = uicontrol(...
        'parent',fig,...
        'Style','edit',...
        'UserData',s(loop).Field,...
        'TooltipString',sprintf('%s',s(loop).Tooltip),...
        'Position',[posx posy maxboxwidth height(loop)],...
        'HorizontalAlignment','left',...
        'String',stri,...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize);      
  end
  
  %Add the label name
  if ~isequal(type{loop},'logical')
    uicontrol(fig,...
      'Style','text',... %note do not add userdata here
      'TooltipString',sprintf('%s',s(loop).Tooltip),...
      'Position',[leftmarg posy maxlabelwidth+xmarg height(loop)],...
      'Stri',sprintf('%s',s(loop).Label),...
      'HorizontalAlignment','left',...
      'ForegroundColor',fc,...
      'BackgroundColor',bc,...
      'Fontsize',fontsize);
  end
  
end

posx = leftmarg;
posy = reqh-(topmarg+cumheight(loop+1))-bottommarg;

%Add cancel button
uicontrol(fig,...
  'Style','pushbutton',...
  'Position',[posx posy 100 height(loop+1)],...
  'Callback', @(u,v)cancel_Callback(u,v),...
  'String',dprintf('Cancel'),...
  'ForegroundColor',fc,...
  'BackgroundColor',bc,...
  'Fontsize',fontsize);
  
%Add ok button
%[reqw-50-rightm posy-7 50 height(n+1)]
uicontrol(fig,...
  'style','pushbutton',...
  'Position',[reqw-rightmarg-50 posy 50 height(loop+1)],...
  'Callback', @(u,v)ok_Callback(u,v),...
  'Stri',dprintf('Ok'),...
  'ForegroundColor',fc,...
  'BackgroundColor',bc,...
  'Fontsize',fontsize);

%activate uiwait to wait for user selection until user clicks cancel or ok.
uiwait(fig);

%----------------------------
function cancel_Callback(h,~)
%----------------------------
%If cancel then ...

uiresume(h.Parent);
h.Parent.UserData = 'Cancelled';

%------------------------
function ok_Callback(h,~)
%------------------------
%Callback function, just resume

uiresume(h.Parent);
h.Parent.UserData = '';
