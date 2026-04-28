function [outs,ok] = myinputstruct(s,tit,width,editcallback,fighandle,illustrationfile,illustrationsize,fontsize)
%Automatically create input dialogs with Medviso style
%This function replaces the old inputstruct.
%
%MYINPUTSTRUCT(S, [TITLE], [WIDTH], [EDITCALLBACK],[FIGHANDLE],[ILLUSTRATIONFILE],[ILLUSTRATIONSIZE],[FONTSIZE])
%- S = Struct, see below for example.
%- TIT = Title string.
%- WIDTH = Width of editbox (optional, default = 10)
%- EDITCALLBACK = optional callback that is called when users edits a field.
%- FIGHANDLE = optional handle to figure, when supplied used to position input
%- ILLUSTRATIONFILE = optional filename that will be showed as an illustration (normally the file is something like ['+3segment3dp' filesep 'illustration_something.png'])
%- ILLUSTRATIONSIZE = optional on how big the illustration is, scale 0..1
%- FONTSIZE = optional fontsize scaling, default is 1
%
%Note that the function does not call translation and you need to do this
%manually before calling the function.
%
%Example:
%
%s = []; %reset
%s(1).Field = 'Name';
%s(1).Label = dprintf('Name'); %Use dprintf only if you want the input to be translated.
%s(1).Default = 'Box';
%s(2).Field = 'H';
%s(2).Label = dprintf('Height [mm]');
%s(2).Default = 10;
%s(3).Field = 'W';
%s(3).Label = dprintf('Width [mm]');
%s(3).Default = 20;
%s(4).Field = 'Action';
%s(4).Label = dprintf('Action');
%s(4).Default = {'Fill','Erase','Invert','Ignore'}; %Listbox, returned as value.
%s(4).Value = 2; %Set 'Erase' as default
%s(4).Height = 2; %Set the height as 1, i.e only show two options at the time.
%s(5).Field = 'Visible';
%s(5).Label = dprintf('Visible');
%s(5).Default = true;
%s(6).Field = 'Comment';
%s(6).Label = dprintf('Comment');
%s(6).Default = "Use string to get multi-line input";
%s(6).Height = 2;
%s(7).Field = 'Dir';
%s(7).Label = dprintf('Directory');
%s(7).Default = @dir; %Use a function handle to create a pushbutton. If value is [] or zero then pushing does not end the input.
%s(7).Max = 1; %if 1 then quit input after exection of the callback, else then wait for user to press ok or cancel before exiting. The output of the callback (if any) is stored in outs as for the other fields.
%s(7).Tooltip = 'Displays the files in the current directory and ends the input';
%
%%It is also possible to manually force the type of the input to be of a
%%specific type by adding the field Type. Examples are 'pushbutton', or
%%'text' that creates a non-editable line of text.
%
%[outs,ok] = myinputstruct(s,'Input box parameters',10); %10 is 10 characters minimum
%
%If user press cancel or close then ok = false;
%
%--- Special fieldnames ---
%- PathName if this then a browse button is shown
%
%--- Optional fields ---
%- Tooltip, sets Tooltip string
%- Enable, allow to enable / disable certain options, used in conjunction with editcallback option.
%- Visible allow to hide objects, often used to avoid persistent variables inputupdate callbacks. Should be true/false.
%- Value, set updated value, selects different options.
%- Max and Min, if Max > Min for cell input, then this allows for multiple choice in the listbox
%- Height, sets the height of the option. Especially useful for creating dropdown lists. This has no effect on single row inputs such as scalars.
%
%Example on editcallback function usage:
%---------------------------------------
%How to call it:
%
%[outs,ok] = myinputstruct(s,'Input box parameters',10,@editsettings_Callback); %10 is 10 characters minimum
%
%%Sample function definition:
%
% function [s,ok] = editsettings_Callback(outs,s,clickedon)
% %Called from myinputstruct when user have made a change. Return the 
% %modified setting struct s. outs is the parsed output that can be 
% %helpful, and clicked on is the name of the Field for what the user
% %clicked.
%
% ok = false; %when this is true, then it is as the user pressed ok and myinputstruct exits
%
% if isequal(clickedon,'H')
%   n = myinputstructfindfield('H',s)
%   s(n).Value = outs.H; %Set the width equal to the H
% end
%
% nvisible = myinputstructfindfield('Visible',s)
% naction = myinputstructfindfield('Action',s)
% if s(nvisible).Value %Can also write this as outs.Visible
%   s(naction).Enable = 'on'; %If Visible checked then enable Action
% else
%   s(nactioin).Enable = 'off'; %If Visible no checked then disable field Action
% end

%REMOVED FUNCTIONALITY (FOR NOW)
%- Limits (default [-inf inf])
%- ValueDisplayFormat (default '%11.4g')
%- LowerLimitInclusive (default 'on')
%- UpperLimitInclusive (default 'on')
%
%Hint to input decimals use any of the int/uint classes

%Einar Heiberg

%Plans for the future
%- Automatic determination of how many many rows there are and then possible split into multiple screens where a user can go back and forth.
%- DONE Possibility to have graphical screen to help users. For instance to show graphically what parameters mean such as height, width, depths etc...
%- DONE Nicer dropdown behaviour

%#ok<*GVMIS> 

global DATA 

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

%Default
if ~isfield(s,'Default')
  for loop = 1:length(s)
    s(loop).Default = '';
  end
end

%Tooltip
if ~isfield(s,'Tooltip')
  for loop = 1:length(s)
    s(loop).Tooltip = '';
  end
end

%Enable
if ~isfield(s,'Enable')
  for loop = 1:length(s)
    s(loop).Enable = 'on';
  end
end

%Visible
if ~isfield(s,'Visible')
  for loop = 1:length(s)
    s(loop).Visible = true;
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

%Type
if ~isfield(s,'Type')
  for loop = 1:length(s)
    s(loop).Type = ''; %if empty then guess based on data type from Default
  end
end

%Fix good values for enable
for loop = 1:length(s)
  if isequal(s(loop).Enable,false)
    s(loop).Enable = 'off';
  end
  if isequal(s(loop).Enable,true) || isempty(s(loop).Enable)
    s(loop).Enable = 'on';
  end
end

%Fix good values for visible
for loop = 1:length(s)
  if isempty(s(loop).Visible)
    s(loop).Visible = true;
  end
end

%Value
if ~isfield(s,'Value')
  s(1).Value = [];
end
for loop = 1:length(s)
  if isempty(s(loop).Value)
    switch class(s(loop).Default)
      case 'function_handle'
        s(loop).Value = 0;
      case 'cell'
        s(loop).Value = 1;
      otherwise
        s(loop).Value = s(loop).Default; %Copy the default to Value
    end
  end
end

%Height
if ~isfield(s,'Height')
  s(1).Height = []; %This sets all others to empty as well. 
end

%Ensure min, max have good values
for loop = 1:length(s)

  if isempty(s(loop).Max)
    if isnumeric(s(loop).Value)
      s(loop).Max = inf;
    end
    if iscell(s(loop).Default)
      s(loop).Max = 0;
    end
  end

   if isempty(s(loop).Min)
    if isnumeric(s(loop).Value)
      s(loop).Min = -inf;
    end
    if iscell(s(loop).Default)
      s(loop).Min = 0;
    end
  end

end

%Label
for loop = 1:length(s)
  if isempty(s(loop).Label)
    s(loop).Label = s(loop).Field;
  end
end

if nargin<2
  tit = '';
else
  % check if exact the same myinputstruct is open, then close it
  fh = findobj('Name',tit,'Type','figure');
  if ~isempty(fh)
    delete(fh);
  end
end

if nargin<3
  width = 10; %characters
end

if nargin<4 || isempty(editcallback)
  editcallback = @(a,b,c)validate(a,b,c);
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

if nargin < 8 
  fontsize = 13;
end

%Call creation function
[fig, s] = createit(s,dprintf('%s',tit),width,editcallback,fighandle,illustrationfile,illustrationsize,fontsize); %There is a uiwait so it returns when user clicked ok

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
  for loop = 1:length(s)
    outs.(s(loop).Field) = s(loop).Default;
  end
else
  s = fig.UserData.S; %Ensure get the update one
  [outs,ok,~] = parse(fig,s); %This is where the input is parsed and converted to outs
end

%Close the figure
if ishandle(fig)
  delete(fig);
end

%----------------------------------
function [outs,ok,s] = parse(fig,s)
%----------------------------------
%Gets the edited results

outs = [];
try
  ok = s(1).ok; %Not sure about this line
catch
  ok = true; %default
end

%Copy default and setup
for loop = 1:length(s)
  outs.(s(loop).Field) = s(loop).Default;
end

c = fig.Children;

%Loop over elements to get them
for loop = 1:length(c)
  ud = c(loop).UserData;

  if ~isempty(ud) %if not isempty user data 

    if ~ishghandle(ud)

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
        s(id).Value = outs.(ud);
      end

      if islogical(s(id).Default)
        outs.(ud) = logical(c(loop).Value);
        s(id).Value = outs.(ud);
      end

      if ischar(s(id).Default) || isstring(s(id).Default)
        outs.(ud) = c(loop).String;
        s(id).Value = c(loop).String;
      end

      if iscell(s(id).Default)
        outs.(ud) = c(loop).Value;
        s(id).Value = outs.(ud);
      end      

    end

  end %~isempty(ud)
end

%------------------------------
function s = validate(outs,s,~) %last argument is clickedon, unused
%------------------------------
%Validate the input and returns output

%Loop over in the struct
errors = false;
message = '';
for loop = 1:length(s)
  fname = s(loop).Field;
  if isnumeric(outs.(fname)) % && ~isnan(outs.(fname))
    v = s(loop).Value;
    if isempty(v)
      v = NaN;
    end
    if (length(v) <= 1) && isnan(v) && ~isequal(s(loop).Type,'text')
      errors = true;
      message = dprintf('Invalid response for %s',s(loop).Label);
      s(loop).Value = s(loop).Default;
    end
    if  ~iscell(s(loop).Default) % exclude cell from checking, if max-min>0 then multiple choice, otherwise single choice
      if v<s(loop).Min
        errors = true;
        message = sprintf('%s. %s %0.5g',dprintf('Too low value'),dprintf('Minimum value is'),s(loop).Min);
        s(loop).Value = s(loop).Min;
      end
      if v>s(loop).Max
        errors = true;
        message = sprintf('%s. %s %0.5g',dprintf('Too large value'),dprintf('Maximum value is'),s(loop).Max);
        s(loop).Value = s(loop).Max;
      end
    end
  end
end

if errors
  myfailed(message);
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

%Filter
%f = [1 2 1];
%f = f/sum(f);
%im = conv2(conv2(im,f,'same'),f(:),'same');

%Fix colors
imr = fg(1)*(1-im)+bc(1)*im;
img = fg(2)*(1-im)+bc(2)*im;
imb = fg(3)*(1-im)+bc(3)*im;
im = cat(3,imr,img,imb);

%-------------------------------------------------------------------------------------------------------------
function [fig, s] = createit(s,tit,minchars,editcallback,fighandle,illustrationfile,illustrationsize,fontsize)
%-------------------------------------------------------------------------------------------------------------
%Creates the figure and populates it

global DATA

%get color settings
if  ~isempty(DATA)
  bc = DATA.GUISettings.BackgroundColor;
  fc = DATA.GUISettings.ForegroundColor;
else
  % default values
%   bc = [0.2118 0.2353 0.2824];
%   fc = [0.9400 0.9400 0.9400];
  bc = [0.9400 0.9400 0.9400];
  fc = [0 0 0];
end

UserData = [];
UserData.Status = 'Cancelled';
UserData.Callback = editcallback;
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
fontf = fontsize-1; %width of character in pixels
rowheight = 24;

%Loop over all to compute width and height
boxwidth = zeros(1,n+1);
height = zeros(1,n+1);
datatype = cell(1,n+1);
defaultv = zeros(1,n+1);

for loop = 1:n
  switch class(s(loop).Default)
    case {'single','double','int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
      %Numeric
      boxwidth(loop) = max(length(mydouble2str(s(loop).Value)),minchars)*fontf;
      height(loop) = rowheight;
      datatype{loop} = 'numeric';
    case 'logical'
      boxwidth(loop) = 0;
      height(loop) = rowheight;
      datatype{loop} = 'logical';
    case 'char'
      if contains(s(loop).Field,'PathName')
        %special case with pathname
        boxwidth(loop) = max(max(length(s(loop).Value),20),minchars)*fontf;
        height(loop) = rowheight; %discuss if height should affect this
        datatype{loop} = 'pathname';
      else
        %Normal
        boxwidth(loop) = max(length(s(loop).Value),minchars)*fontf;
        height(loop) = rowheight;
        datatype{loop} = 'char';
      end
    case 'string'
      boxwidth(loop) = max(length(s(loop).Value),minchars)*fontf;
      if isempty(s(loop).Height)
        height(loop) = 5*rowheight;
      else
        height(loop) = s(loop).Height*rowheight;
      end

      datatype{loop} = 'string';
    case 'cell'
      maxlen = 0;
      defaultv(loop) = 1; %Default value for the selection
      if iscell(s(loop).Value)
        s(loop).Value = 1;
      end
      numentries = length(s(loop).Default);
      for subloop = 1:numentries
        if length(s(loop).Default{subloop})>maxlen
          maxlen = length(s(loop).Default{subloop});
        end
        temp = s(loop).Default{subloop};
        if ~isempty(temp)
          if any(temp=='*')
            defaultv(loop) = subloop; %Set default value for the selection
            s(loop).Value = subloop;
          end
        end        
      
      end
      boxwidth(loop) = max(maxlen,minchars)*fontf;
      maxentries = 15;
      if numentries > maxentries
        numentries = maxentries;
      end
      if ~isempty(s(loop).Height)
        height(loop) = rowheight*s(loop).Height; %Use user supplied height
      else
        height(loop) = rowheight*numentries; %User have not supplied height
      end
      datatype{loop} = 'cell';
    case 'function_handle'
      boxwidth(loop) = 0;
      height(loop) = rowheight;
      datatype{loop} = 'function_handle';
  end
  
  %If not visible, then no height
  if ~s(loop).Visible
    height(loop) = 0;
  end

end

%Add button
boxwidth(n+1) = 0;
height(n+1) = 20;

maxlabelwidth = getmaxlabelwidth(fig,fontsize,s);
maxboxwidth = max(boxwidth(1:n))+xmarg;
ymargvector = repmat(ymarg,1,length(height));
ymargvector(height==0) = 0;
cumheight = [0 cumsum(height+ymargvector)];
numvisible = sum((height>0));

%Required width and height
reqw = leftmarg+maxlabelwidth+xmarg+maxboxwidth+rightmarg;
reqh = topmarg+sum(height)+numvisible*ymarg+bottommarg;

%Load the illustration
if ~isempty(illustrationfile)

  %The illustration file is normally something like ['+3segment3dp' filesep 'illustration_something.png']

  %Fix folder depending on location
  if isdeployed
    foldername = DATA.getsoftwarenamenospace;
    if length(foldername) > 12
      foldername = foldername(1:12);
    end
    illustrationfilename = [ctfroot filesep foldername filesep illustrationfile];
  else
    illustrationfilename = [DATA.SegmentFolder filesep illustrationfile];
  end

  try
    illustration = imread(illustrationfilename);
  catch
    logdisp(sprintf('Could not read illustration %s',illustrationfilename));
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
  %ax.Units = 'Normalized';
  %ax.Position = [0.5 0.1 0.55 0.8];

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
  posx = leftmarg+maxlabelwidth+xmarg;
  
  shouldfocus = false; %EiH: Default (I think)

  %--- Find what type of object we want
  if ~isempty(s(loop).Type)
    datatype{loop} = lower(s(loop).Type);
  end

  %Fix stri
  if s(loop).Visible
    switch datatype{loop}
      case 'numeric'
        %Numeric
        stri = mydouble2str(s(loop).Value);
        stri = stri(stri~=' ');

        u(loop) = uicontrol(...
          'parent',fig,...
          'Style','edit',...
          'UserData',s(loop).Field,...
          'Tag',s(loop).Field,...
          'TooltipString',sprintf('%s',s(loop).Tooltip),...
          'Position',[posx posy maxboxwidth height(loop)],...
          'HorizontalAlignment','left',...
          'String',stri,...
          'ForegroundColor',fc,...
          'BackgroundColor',bc,...
          'Fontsize',fontsize,...
          'Enable',s(loop).Enable,...
          'Callback',@editcallbackhelper);
      case 'logical'
        u(loop) = uicontrol(...
          'parent',fig,...
          'Style','checkbox',...
          'UserData',s(loop).Field,...
          'Tag',s(loop).Field,...
          'TooltipString',sprintf('%s',s(loop).Tooltip),...
          'Position',[posx posy maxboxwidth height(loop)],...
          'String','',...
          'Value',s(loop).Value,...
          'ForegroundColor',fc,...
          'BackgroundColor',bc,...
          'Fontsize',fontsize,...
          'Enable',s(loop).Enable,...
          'Callback',@editcallbackhelper);

      case {'char','pathname'}
        stri = s(loop).Value;
        u(loop) = uicontrol(...
          'parent',fig,...
          'Style','edit',...
          'UserData',s(loop).Field,...
          'Tag',s(loop).Field,...
          'TooltipString',sprintf('%s',s(loop).Tooltip),...
          'Position',[posx posy maxboxwidth height(loop)],...
          'HorizontalAlignment','left',...
          'String',stri,...
          'ForegroundColor',fc,...
          'BackgroundColor',bc,...
          'Fontsize',fontsize,...
          'Enable',s(loop).Enable,...
          'Callback',@editcallbackhelper);
      case 'string'
        stri = s(loop).Value;
        u(loop) = uicontrol(...
          'parent',fig,...
          'Style','edit',...
          'UserData',s(loop).Field,...
          'Tag',s(loop).Field,...
          'TooltipString',sprintf('%s',s(loop).Tooltip),...
          'Position',[posx posy maxboxwidth height(loop)],...
          'HorizontalAlignment','left',...
          'Max',3,'Min',1,...
          'String',stri,...
          'ForegroundColor',fc,...
          'BackgroundColor',bc,...
          'Fontsize',fontsize,...
          'Enable',s(loop).Enable,...
          'Callback',@editcallbackhelper);
      case 'cell'
        if (numfields == 1) && isequal(s(loop).Max,s(loop).Min) % if max == min, then it is single choice
          callback = @listbox_Callback;
          keyfcn = @keyhelper;
          shouldfocus = true;
        else
          callback = @editcallbackhelper;
          keyfcn = '';
          shouldfocus = false;
        end
        if ~isempty(s(loop).Height) && s(loop).Height == 1
          % maked dropdown
          uicontrolstyle = 'popupmenu';
        else
          uicontrolstyle = 'listbox';
        end
        fig.KeyPressFcn = keyfcn;
        stri = s(loop).Default;
        u(loop) = uicontrol(...
          'parent',fig,...
          'Style',uicontrolstyle,...
          'UserData',s(loop).Field,...
          'Tag',s(loop).Field,...
          'TooltipString',sprintf('%s',s(loop).Tooltip),...
          'Position',[posx posy maxboxwidth height(loop)],...
          'HorizontalAlignment','left',...
          'String',stri,...
          'Value',s(loop).Value,...
          'Max',s(loop).Max,...
          'Min',s(loop).Min,...
          'ForegroundColor',fc,...
          'BackgroundColor',bc,...
          'Fontsize',fontsize,...
          'Enable',s(loop).Enable,...
          'KeyPressFcn',keyfcn,...
          'Callback',callback);
      case 'function_handle'
        stri = s(loop).Label;
        u(loop) = uicontrol(...
          'parent',fig,...
          'Style','pushbutton',...
          'UserData',s(loop).Field,...
          'Tag',s(loop).Field,...
          'TooltipString',sprintf('%s',s(loop).Tooltip),...
          'Position',[leftmarg posy maxlabelwidth height(loop)],...
          'String',stri,...
          'Value',s(loop).Value,...
          'Max',s(loop).Max,...
          'ForegroundColor',fc,...
          'BackgroundColor',bc,...
          'Fontsize',fontsize,...
          'Enable',s(loop).Enable,...
          'Callback', @(u,v)pushbutton_Callback(u,v,s(loop).Default));
      case 'pushbutton'
        stri = s(loop).Label;

        u(loop) = uicontrol(...
          'parent',fig,...
          'Style','pushbutton',...
          'UserData',s(loop).Field,...
          'Tag',s(loop).Field,...
          'TooltipString',sprintf('%s',s(loop).Tooltip),...
          'Position',[leftmarg posy maxlabelwidth height(loop)],...
          'String',stri,...
          'ForegroundColor',fc,...
          'BackgroundColor',bc,...
          'Fontsize',fontsize,...
          'Enable',s(loop).Enable,...
          'Callback', @editcallbackhelper);
    end
  end
  
  %Add the label name
  switch datatype{loop}
    case 'pathname'
      uicontrol(fig,...
        'Style','pushbutton',... %note do not add userdata here
        'TooltipString',sprintf('%s',s(loop).Tooltip),...
        'Position',[leftmarg posy maxlabelwidth+0*xmarg height(loop)],...
        'Stri',sprintf('%s',s(loop).Label),...
        'HorizontalAlignment','left',...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize,...
        'Enable',s(loop).Enable,...
        'Callback',@(h,a)pathnamecallbackhelper(h,a,loop),...
        'UserData',u(loop)); %Link to it its field
    case 'text'
      uicontrol(fig,...
        'Style','text',... %note do not add userdata here
        'TooltipString',sprintf('%s',s(loop).Tooltip),...
        'Position',[leftmarg posy maxlabelwidth+0*xmarg height(loop)],...
        'Stri',sprintf('%s',s(loop).Label),...
        'HorizontalAlignment','left',...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize,...
        'Enable',s(loop).Enable,...
        'UserData',s(loop).Field); %Link to it its field        
    case {'function_handle','pushbutton'}
      %Do nothing
    otherwise
      uicontrol(fig,...
        'Style','text',... %note do not add userdata here
        'TooltipString',sprintf('%s',s(loop).Tooltip),...
        'Position',[leftmarg posy maxlabelwidth+xmarg height(loop)],...
        'Stri',sprintf('%s',s(loop).Label),...
        'HorizontalAlignment','left',...
        'ForegroundColor',fc,...
        'BackgroundColor',bc,...
        'Fontsize',fontsize,...
        'Enable',s(loop).Enable,...
        'UserData',u(loop)); %Link to it its field
  end

end  %for loop
  
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

    if ~isempty(DATA.Buffer.ListboxValue)
      % find listbox
      lb = findobj(fig.Children,'style','listbox');
      for loop = 1:numel(lb)
        if ~isempty(lb(loop))
          lbv = popfrombuffer('ListboxValue');
          if ~isempty(lbv)
            lb(loop).Value = lbv;
          end
        end
      end
    end
    
    if ~isempty(DATA.Buffer.EditString)
      % find editboxes
      editboxes = findobj(fig.Children,'style','edit');
      if ~isempty(editboxes)
        for loop = 1:numel(editboxes)
          edstri = popfrombuffer('EditString');
          if ~isempty(edstri)
            editboxes(loop).String = edstri;
          end
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
    hlistbox = findobj(fig.Children,'Style','listbox');
    if ~isempty(hlistbox)
      uicontrol(hlistbox);
    end
  end

  %activate uiwait to wait for user selection until user clicks cancel or ok.
  uiwait(fig);
end

%----------------------------
function cancel_Callback(h,~)
%----------------------------
%If cancel then ...

uiresume(h.Parent);
UserData = h.Parent.UserData ;
UserData.Status = 'Cancelled';
h.Parent.UserData = UserData;

%------------------------
function ok_Callback(h,~)
%------------------------
%Callback function, just resume

uiresume(h.Parent);
UserData = h.Parent.UserData;
UserData.Status = 'Ok';
h.Parent.UserData = UserData;

%----------------------------------
function pushbutton_Callback(h,a,f)
%----------------------------------
%Pushbutton function, the user have pressed on some of the user defined pushbutton

%Evaluate the callback f
n = nargout(f);  %check output arguments that f gives
if n > 0 || n < 0 %varargout produces negative number
  output = feval(f);
else
  feval(f);
  output = [];
end

%Get the value to see if we should exit or not
value = h.Max;

if isequal(value,1)
  %Then we should exist after
  uiresume(h.Parent);

  UserData = h.Parent.UserData;
  UserData.Status = 'Ok';
  h.Parent.UserData = UserData;
end

%--- Store output, if no output then store [] as then output is set to [] above
UserData = h.Parent.UserData;
s = UserData.S;
tag = h.Tag;
for loop = 1:length(s)
  if isequal(tag,s(loop).Field)
    s(loop).Default = output; %Store the output
  end
end

%Store back
UserData.S = s;
h.Parent.UserData = UserData;
editcallbackhelper(h,a);

%-------------------------------
function editcallbackhelper(h,a)
%-------------------------------
%User have pressed edit.

%Figure out handle to parent
fig = h.Parent;

%Get userdata
UserData = fig.UserData;

%Check if there is a callback, if not, no need to do anything.
if isempty(UserData.Callback)
  return
end

%Get the s
[outs,~,s] = parse(fig,UserData.S);
UserData.Outs = outs;
UserData.S = s;

%Get what clicked on
clickedon = get(a.Source,'Tag');

%Store to figure
fig.UserData = UserData;

%Update interface, to reflect users edits
update(fig);

%Call callback to allow coder make changes
ok = false;
try %This is an ugly method to handle extra optional output argument
  [UserData.S,ok] = feval(UserData.Callback,outs,UserData.S,clickedon);
catch
  [UserData.S] = feval(UserData.Callback,outs,UserData.S,clickedon);
end

%Store to figure
fig.UserData = UserData;

if ok
  ok_Callback(h); %This ends it
end

%Update interface with new changes
update(fig);

%-------------------------------------
function pathnamecallbackhelper(h,~,n)
%-------------------------------------
%User have pressed pathname (i.e browse).

%Figure out handle to parent
fig = h.Parent;

%Get userdata
UserData = fig.UserData;

if isempty(UserData)
  return
end

%Get s
s = UserData.S;
pathname = s(n).Default;

%Call function to select
selpath = myuigetdir(pathname,s(n).Label);

%Check if user cancelled
if isequal(selpath,0)
  selpath = pathname;
end

%Store to UserData struct
s(n).Default = selpath;
s(n).Value = selpath;
UserData.S = s;

%Store to figure
fig.UserData = UserData;

%Update interface with new changes
update(fig);

%-------------------------------
function s = getvisiblestring(v)
%-------------------------------
%Convert logical to on / off

if v
  s = 'on';
else
  s = 'off';
end

%-------------------
function update(fig)
%-------------------

%Get userdata
UserData = fig.UserData;

%Get s struct
s = UserData.S;

%Get Children
c = fig.Children;

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
    
    if ~ishghandle(ud)
      %Normal
      set(c(loop),'Enable',s(id).Enable);

      if isnumeric(s(id).Default)
        stri = mydouble2str(s(id).Value);
        set(c(loop),...
          'String',stri,...
          'Visible',getvisiblestring(s(id).Visible));
          %...'Enable',s(id).Enable ...
      end
    
      if islogical(s(id).Default)
        set(c(loop),...
          'Value',s(id).Value...
          ...'Enable',s(id).Enable ...
          );
      end
    
      if ischar(s(id).Default) || isstring(s(id).Default)
        set(c(loop),...
          'String',s(id).Value,...
          'Visible',getvisiblestring(s(id).Visible));
          %...'Enable',s(id).Enable ...
      end
    
      if iscell(s(id).Default)
        set(c(loop),...
          'String',s(id).Default,...
          'Value',s(id).Value...
          ...'Enable',s(id).Enable ...
          );
      end

      if isa(s(id).Default, 'function_handle')
        set(c(loop),...
          'String',s(id).Label...
          ...'Enable',s(id).Enable ...
          );
      end

    end % ~ishghandle(ud)

  end % ~isempty(ud)

end %loop over handles

%Update the text object(s), text objects have handle to corresponding edit
%box as userdata

for loop = 1:length(s)
  ud = c(loop).UserData;
  if ~isempty(ud) 
    if ishghandle(ud) 
      if (ud>0)
        enablestatus = get(ud,'Enable'); %get if the corresponding edit box is enabled
        set(c(loop),'Enable',enablestatus)
      end
    end
  end
end

%------------------------------
function stri = mydouble2str(d)
%------------------------------
%Convert numeric to string

stri = strtrim(sprintf('%11.4g',d));

%----------------------------------------
function listbox_Callback(h,~)
%----------------------------------------

if strcmp(h.Type,'figure')
  hlistbox = findobj(h.Children,'Style','listbox');
  hfig = h;
  if isempty(hlistbox)
    return
  end
else
  hlistbox = h;
  hfig = h.Parent;
end

keypressed =  hfig.UserData.KeyPressed;
if keypressed
  % hfig.UserData.KeyPressed
  hfig.UserData.KeyPressed = false;
else
  ok_Callback(hlistbox);
end


%----------------------------------------
function keyhelper(h,event)
%----------------------------------------

if strcmp(h.Type,'figure')
  hlistbox = findobj(h.Children,'Style','listbox');
  hfig = h;
  if isempty(hlistbox)
    return
  end
else
  hlistbox = h;
  hfig = h.Parent;
end

if strcmp(event.EventName,'KeyPress')
  key = event.Key;
  hfig.UserData.KeyPressed = true;
else
  return
end

%--- If pressed return, then take currently selected.
if isequal(key,'return')
  ok_Callback(hlistbox);
  return;
end

if isempty(key)
  return;
end

%--- Check for special keys
switch key
  case 'escape'
    cancel_Callback(hlistbox)
    return;
end

%----------------------------------------
function maxwidth = getmaxlabelwidth(fig,fontsize,s)
%----------------------------------------
% Function to get extent width for the longest label

% Set invisible uicontrol to determine extent
extentcontrol = uicontrol(...
  'parent',fig,...
  'Style','pushbutton',...
  'Visible','off',...
  'Fontsize',fontsize...
  );

% Extract 'Label' field values into a cell array
labels = {s.Label};

% Replace string for field 'PathName'
indpath = matches({s.Field},'PathName');
labels(indpath) = {dprintf('Browse')};

% Find strings containing newline
splitind = contains(labels, newline);
if any(splitind)
  % Split strings on newline
  for n = find(splitind)
    labels{n} = strsplit(labels{n}, newline);
  end
  % Merge all into one cell array
  labels = cat(2, labels{:});
end

% Find the string with the maximum length
celllength = cellfun(@length,labels);
[~,ind] = max(celllength);

% Set 'extentcontrol' string and get its width extent
extentcontrol.String = ['aa',labels{ind}]; % add 'aa' to make string a bit longer
ext = extentcontrol.Extent;
maxwidth = ext(3);
