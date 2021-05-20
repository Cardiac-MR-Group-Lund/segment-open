function [varargout] = mymenumultiplechoice(header,varargin)
%m = MYMENU(header,items,defaultstring)
%m = MYMENU(header,items,defaultstring,maxselections)
%m = MYMENU(header,item1,item2,item3,...)
%
%Same as menu, but with the following additions:
%- more elegant! uses keyboard to select items.
%- modal display
%- can use default string
%- last optional argument is figure handle 

%Einar Heiberg
global DATA

persistent handles

if nargout>0
  varargout = cell(1,1);
end

if isa(DATA, 'maingui') && DATA.Testing
  testing = DATA.Testing;
else
  testing = false;
end


if testing
  %If testing then also DATA.Buffer should exist
  v = popfrombuffer('Mymenu');
  if isempty(v)
    myfailed('Menu buffer is empty.');
    return;
  elseif isnan(v) 
    %Test of warning and nbr of options (one less than cell size)
    pushtobuffer('Warnings',sprintf('%d %s',numel(varargin{1})-1,header));
    v = 2;
  end
  varargout{1}=v;
  return;
end

%Check if special commands
switch header
%   //case 'ok'
%     print('ok pressed')
  case {'cancel'}
    %--- User has pressed selection or has pressed OK button to confirm
    %selection
    closehelper(handles);
    handles.m = [];
    if nargout>0
      varargout{1} = handles.m;
    end
    return;
  case '_close'
    closehelper(handles);
    if nargout>0
      varargout{1} = handles.m; %return;
    end
    return
  case {'_selection'}
    %--- User has pressed selection or has pressed OK button to confirm
    %selection
    handles = selectionhelper(handles);
    if nargout>0
      varargout{1} = handles.m; %return;
    end
    return;
  case '_key'
    handles = keyhelper(handles,varargin{1});
    if nargout>0
      varargout{1} = handles.m; %return;
    end
    return;
end %End of case clause

%--- Init if got to here
if nargin<2
  myfailed('Too few input arguments.');
  return;
end

if length(varargin)>1
  %Last argument is handle to figure.
  n = length(varargin);
  if ishandle(varargin{n})
    fighandle = varargin{n};
    varargin = varargin(1:(n-1)); %Remove last...
  else
    fighandle =[];
  end
else
  if isa(varargin{1},'cell')
    n = length(varargin{1});
    if ishandle(varargin{1}{n})
      fighandle = varargin{1}{n};
      varargin = varargin{1}(1:(n-1)); %Remove last...
    else
      fighandle = [];
    end
  else
  fighandle = [];
  end
end

%Extract paramters
if iscell(varargin{1})
  options = varargin{1};
  if length(varargin)==2
    default = varargin{2};
  else
    default = '';
  end
else
  options = cell(1,length(varargin));
  [options{:}] = deal(varargin{:});
  default = '';
end

% %translation
% header = translation.dictionary(header);
% for k=1:numel(options)
%   options{k} = translation.dictionary(options{k});
% end

%Call initialization code.
handles = init(header,fighandle,options,default);
try
  uiwait(handles.fig); %wait until terminates
catch %#ok<CTCH>
end

if nargout>0
  varargout{1} = handles.m; %return;
  showout(handles.m);
end
flushlog;

%------------------------------------------
function handles = selectionhelper(handles)
%------------------------------------------
%Get selection
handles.m = mygetlistbox(handles.menulistbox);
% closehelper(handles);

%----------------------------------------
function handles = keyhelper(handles,key)
%----------------------------------------

%--- If pressed return, then take currently selected.
if isequal(key,'return')
  handles = selectionhelper(handles);
  closehelper(handles);
  return;
end

if isempty(key) || contains(key,'shift') || contains(key,'ctrl')
  return;
end

%--- Check for special keys
switch key
  case 'uparrow'
    v = mygetlistbox(handles.menulistbox);
    v = max(v-1,1);
    set(handles.menulistbox,'value',v);
    return;
  case 'downarrow'
    v = mygetlistbox(handles.menulistbox);
    v = min(v+1,length(handles.options));
    set(handles.menulistbox,'value',v);
    return;
  case 'escape'
    closehelper(handles);
    handles.m = [];
    return;
end

%--- Find menu items beginning with the same letter.
ind = [];
for loop=(mygetlistbox(handles.menulistbox)+1):length(handles.options)
  temp = handles.options{loop};
  if (~isempty(temp))
    if isequal(lower(temp(1)),key(1))
      ind = [ind loop];
    end
  end
end

if isempty(ind)
  %Try to loop from begining
  for loop=1:length(handles.options)
    temp = handles.options{loop};
    if (~isempty(temp))
      if isequal(lower(temp(1)),key(1))
        ind = [ind loop];
      end
    end
  end
end

if isempty(ind)
  return;
end

set(handles.menulistbox,'value',ind(1));

%----------------------------
function closehelper(handles)
%----------------------------
set(handles.fig,'closerequestfcn','');
try
  delete(handles.fig);
catch %#ok<CTCH>
end
    
%-------------------------------
function keypressed(fignum,evnt) %#ok<INUSL>
%-------------------------------
key = getkey(evnt);
mymenumultiplechoice('_key',key);

%----------------------------------------------
function handles = init(header,fighandle,options,default)
%----------------------------------------------

fig = openfig('mymenu.fig','reuse','invisible');%menu not visble until it is horisontally aligned 
setupicon(fig);
setinterfacecolor(fig);
if ~isempty(fighandle)
  myadjust(fig,fighandle);
end

%Extract handles
handles = guihandles(fig);
handles.fig = fig;
handles.header = header;
handles.options = options;
handles.default = default;
handles.m = 0; %default...

%Set handles.
set(handles.fig,'name',translation.dictionary(header));
% set(handles.fig,'name','Select');
hintstr = dprintf('Use CTRL + click to choose multiple options');
set(handles.titletext,'String',hintstr);
set(handles.menulistbox,'String',options,'Max',10,...
  'Callback','mymenumultiplechoice(''_selection'','''')');
set(handles.menulistbox,'keypressfcn',@keypressed);
set(handles.fig,'keypressfcn',@keypressed,'CloseRequestFcn','mymenumultiplechoice(''_close'')');
set(handles.okpushbutton,'Callback','mymenumultiplechoice(''_close'')');
set(handles.cancelpushbutton,'Callback','mymenumultiplechoice(''cancel'')');

translation.translatealllabels(fig);

%Adjust size
p = get(handles.fig,'position');
extrapixels = 0;
if length(handles.options)>5
  extrapixels = 50;
end
if length(handles.options)>10
  extrapixels = 150;
end
if length(handles.options)>15
  extrapixels = 250;
end
if length(handles.options)>20
  extrapixels = 400;
end
p(4) = p(4)+extrapixels;
p(2) = p(2)-extrapixels;
set(handles.fig,'position',p);

pmenu=get(handles.menulistbox,'position');
ptitle=get(handles.titletext,'position');

pmenu(4)=p(4)-ptitle(4)-pmenu(2)-20;
pmenu(2)=pmenu(2);%+extrapixels;
set(handles.menulistbox,'position',pmenu);

ptitle(2)=p(4)-ptitle(4)-10;
set(handles.titletext,'position',ptitle);

try
  if ~isempty(fighandle)
    myadjust(handles.fig,fighandle); %Horisontal alignment.
  end
catch %#ok<CTCH>
end
set(fig,'visible','on');%menu not visble until it is horisontally aligned 

%Find default and mark it
if ~isempty(default)
  defaultpos = [];
  for loop=1:length(options)
    if isequal(options{loop},default)
      defaultpos = loop;
    end
  end
  
  if ~isempty(defaultpos)
    set(handles.menulistbox,'value',defaultpos);
  end
end

%------------------
function showout(d)
%------------------
global DATA

if isa(DATA, 'maingui') && DATA.RecordMacro
  recordmacro = DATA.RecordMacro;
else
  recordmacro = false;
end

if recordmacro
  macro_helper('put',sprintf('pushtobuffer(''Mymenu'',%d); %%select menu item',d));
  macro_helper('switchorder'); %We need to store data in buffer before the callback
end
