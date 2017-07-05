function result = yesno(s,arg,fighandle)
%YESNO(STRI)
%  Asks the user an yes/no question. STRI is the question
%  string.
%
%YESNO(STRI,[],FIGHANDLE)
%  Same as above, but fighandle indicates alignment.
%
%See also MYFAILED, MYWARNING, MYWAITBARSTART, MYMSGBOX.

%Einar Heiberg

global DATA

persistent state 

s = translation.dictionary(s);

state = false;

if (nargin==1)||(nargin==3)
  
  keystroke = popfrombuffer('KeyStroke');
  if ~isempty(keystroke)
    switch keystroke
      case 'yes' 
        result = true;
        mydisp(dprintf('yesno: %s. Yes',s));
        if ~isempty(DATA) && DATA.RecordMacro
          macro_helper('put','pushtobuffer(''KeyStroke'',''yes''); %yes in yesno');
          macro_helper('switchorder'); %We need to store data in buffer before the callback
        end;
        return;
      case 'no'
        result = false;
        mydisp(dprintf('yesno: %s. No',s));
        if ~isempty(DATA) && DATA.RecordMacro
          macro_helper('put','pushtobuffer(''KeyStroke'',''no''); %no in yesno');
          macro_helper('switchorder'); %We need to store data in buffer before the callback
        end;
        return;
      otherwise
        mydisp(dprintf('Yesno got:%s',keystroke));
        error('Expected either ''yes'' or ''no'', got %s.',keystroke);
    end;
  else
    if nargin==3
      handles = initgui(s,fighandle);
    else
      handles = initgui(s);
    end;
  end;
  uiwait(handles.fig);
else
  if isa(DATA, 'maingui') && DATA.RecordMacro
    recordmacro = DATA.RecordMacro;
  else
    recordmacro = false;
  end;
  switch arg
    case 'yes'
      state = true;
      if recordmacro
        macro_helper('put','pushtobuffer(''KeyStroke'',''yes''); %yes in yesno');
        macro_helper('switchorder'); %We need to store data in buffer before the callback
      end;
      mydisp('yes');
    case 'no'
      state = false;
      if recordmacro
        macro_helper('put','pushtobuffer(''KeyStroke'',''no''); %no in yesno');
        macro_helper('switchorder'); %We need to store data in buffer before the callback
      end;
      mydisp('no');      
    case 'yeskeypressed'
      yeskeypressed(s,arg);
    case 'nokeypressed'
      nokeypressed(s,arg);
    case 'return'
      state = false;
  end;
  
  uiresume(s);
  close(s);
end;

result = state;
flushlog;

%--------------------------------------
function handles = initgui(s,fighandle)
%--------------------------------------
%Init the GUI

fig = openfig('yesno.fig','reuse');
translation.translatealllabels(fig);
set(fig,'visible','off');%dlgbox not visble until it is horisontally aligned 
set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

%Horisontal alignment
if nargin>1
  myadjust(fig,fighandle); 
else
  myadjust(fig);
end;

%Extract handles
handles = guihandles(fig);
handles.fig = fig;

%Add text
stri = textwrap({s},50); %was 70 changed to 50 for Mac
numlinestopad = round((5-size(stri,1))/2);
if numlinestopad>0
  stri = [repmat({''},numlinestopad,1);s];
else
end;

set(handles.questiontext,'String',stri);

mydisp(s);

try
  myalign(handles.fig,fighandle); %Align horisontally
catch %#ok<CTCH>
end;
set(fig,'visible','on');%dlgbox not visble until it is horisontally aligned 

%Display question image
load('dialogicons.mat','questIconData','questIconMap');
% questIconMap(256,:)=get(handles.fig,'color');
image(questIconData,'parent',handles.imageaxes);
colormap(handles.imageaxes,questIconMap);
axis(handles.imageaxes,'off');

%Add keypress callback
set(handles.fig,'keypressfcn',@keypressed);
set(handles.yespushbutton','keypressfcn',@yeskeypressed);
set(handles.nopushbutton','keypressfcn',@nokeypressed);

%-------------------------------
function keypressed(fignum,evnt)
%-------------------------------
key = getkey(evnt);

switch key
  case 'return'
    yesno(fignum,'return');
  case {'y','Y'}
    yesno(fignum,'yes');
  case {'n','N'}
    yesno(fignum,'no');
end;

%----------------------------------
function yeskeypressed(fignum,evnt) %#ok<INUSL>
%----------------------------------
key = getkey(evnt);
switch key
  case {'return','y','Y'}
    yesno(gcbf,'yes');
  otherwise
    yesno(gcbf,'no');
end;
      
%----------------------------------
function nokeypressed(fignum,evnt) %#ok<INUSL>
%----------------------------------
key = getkey(evnt);

switch key
  case {'y','Y'}
    yesno(gcbf,'yes');
  otherwise
    yesno(gcbf,'no');
end;