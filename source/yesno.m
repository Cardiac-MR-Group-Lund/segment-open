function result = yesno(s,arg,fighandle)
%YESNO(STRI)
%  Asks the user an yes/no question. STRI is the question
%  string.
%
%YESNO(STRI,[],FIGHANDLE)
%  Same as above, but fighandle indicates alignment.
%
%  If DATA.Silent then message is just displayed in the console
%  and exectution continues as if user had pressed ok.
%
%See also MYFAILED, MYWARNING, MYWAITBARSTART, MYMSGBOX, MYINPUTSTRUCT.

%Einar Heiberg

global DATA

persistent state

%If DATA.Silent then we can assume that user answers yes, which is the
%default.
try
  if DATA.Silent
    disp(sprintf('Yesno: %s : Yes',s)); %#ok<DSPS>
    result = true;
    return;
  end
catch
  %do nothing as error here depends on DATA not initialized
end

s = translation.dictionary(s);

state = false;

if (nargin==1)||(nargin==3)
  
  keystroke = popfrombuffer('KeyStroke');
  if ~isempty(keystroke)
    switch keystroke
      case 'yes'
        result = true;
        disp(sprintf('yesno: %s. Yes',s)); %#ok<DSPS>
        if ~isempty(DATA) && DATA.RecordMacro
          macro_helper('put','pushtobuffer(''KeyStroke'',''yes''); %yes in yesno');
          macro_helper('switchorder'); %We need to store data in buffer before the callback
        end
        return;
      case 'no'
        result = false;
        disp(sprintf('yesno: %s. No',s)); %#ok<DSPS>
        if ~isempty(DATA) && DATA.RecordMacro
          macro_helper('put','pushtobuffer(''KeyStroke'',''no''); %no in yesno');
          macro_helper('switchorder'); %We need to store data in buffer before the callback
        end
        return;
      otherwise
        disp(sprintf('Yesno got:%s',keystroke)); %#ok<DSPS>
        error('Expected either ''yes'' or ''no'', got %s.',keystroke);
    end
  else
    if nargin==3
      handles = initgui(s,fighandle);
    else
      handles = initgui(s);
    end
  end
  uiwait(handles.fig);
else
  if isa(DATA, 'maingui') && DATA.RecordMacro
    recordmacro = DATA.RecordMacro;
  else
    recordmacro = false;
  end
  switch arg
    case 'yes'
      state = true;
      if recordmacro
        macro_helper('put','pushtobuffer(''KeyStroke'',''yes''); %yes in yesno');
        macro_helper('switchorder'); %We need to store data in buffer before the callback
      end
      mydisp('yes');
    case 'no'
      state = false;
      if recordmacro
        macro_helper('put','pushtobuffer(''KeyStroke'',''no''); %no in yesno');
        macro_helper('switchorder'); %We need to store data in buffer before the callback
      end
      mydisp('no');
    case 'yeskeypressed'
      yeskeypressed(s,arg);
    case 'nokeypressed'
      nokeypressed(s,arg);
    case 'return'
      state = false;
  end
  
  uiresume(s);
  close(s);
end

result = state;
flushlog;

%--------------------------------------
function handles = initgui(s,fighandle)
%--------------------------------------
%Init the GUI
global DATA

fig = openfig('yesno.fig','invisible','reuse');
try
  set(fig,'Color',DATA.GUISettings.BackgroundColor);
catch
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
end

setupicon(fig); %set up Segment icon

%Horisontal alignment
if nargin>1
  myadjust(fig,fighandle);
else
  myadjust(fig);
end

%Extract handles
handles = guihandles(fig);
handles.fig = fig;

%Add text
stri = textwrap({s},50); %was 70 changed to 50 for Mac
numlinestopad = round((5-size(stri,1))/2);
if numlinestopad>0
  stri = [repmat({''},numlinestopad,1);s];
else
end

try
  set(handles.questiontext,'String',stri,'BackgroundColor',DATA.GUISettings.BackgroundColor,...
    'ForegroundColor',DATA.GUISettings.ForegroundColor,'FontSize',16);
catch
  set(handles.questiontext,'String',stri);
end

mydisp(s);

try
  myalign(handles.fig,fighandle); %Align horisontally
catch %#ok<CTCH>
end
set(fig,'visible','on');%dlgbox not visble until it is horisontally aligned

%Display question image
try
  switch DATA.GUISettings.BackgroundColor(1)
    case 0.2118
      load('dialogiconsblue.mat','questIconDataBlue','questIconMap');
      image(questIconDataBlue,'parent',handles.imageaxes);
    case 0
      load('dialogiconsblack.mat','questIconDataBlack','questIconMap');
      image(questIconDataBlack,'parent',handles.imageaxes);
    otherwise
      load('dialogicons.mat','questIconData','questIconMap');
      % questIconMap(256,:)=get(handles.fig,'color');
      image(questIconData,'parent',handles.imageaxes);
  end
catch
  load('dialogicons.mat','questIconData','questIconMap');
  % questIconMap(256,:)=get(handles.fig,'color');
  image(questIconData,'parent',handles.imageaxes);
end
colormap(handles.imageaxes,questIconMap);
axis(handles.imageaxes,'off');
set(handles.yespushbutton,'String',dprintf('Yes'));
set(handles.nopushbutton,'String',dprintf('No'));
handles.fig.Name = dprintf('Are you sure?');

%Add keypress callback
set(handles.fig,'keypressfcn',@keypressed);
try
  set(handles.yespushbutton','keypressfcn',@yeskeypressed,'BackgroundColor',...
    DATA.GUISettings.BackgroundColor,'ForegroundColor',DATA.GUISettings.ForegroundColor);
catch
  set(handles.yespushbutton','keypressfcn',@yeskeypressed,'FontSize',12);
end
try set(handles.nopushbutton','keypressfcn',@nokeypressed,'BackgroundColor',...
    DATA.GUISettings.BackgroundColor,'ForegroundColor',DATA.GUISettings.ForegroundColor);
catch
  set(handles.nopushbutton','keypressfcn',@nokeypressed,'FontSize',12);
end

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
end

%----------------------------------
function yeskeypressed(fignum,evnt) %#ok<INUSL>
%----------------------------------
key = getkey(evnt);
switch key
  case {'return','y','Y'}
    yesno(gcbf,'yes');
  otherwise
    yesno(gcbf,'no');
end

%----------------------------------
function nokeypressed(fignum,evnt) %#ok<INUSL>
%----------------------------------
key = getkey(evnt);

switch key
  case {'y','Y'}
    yesno(gcbf,'yes');
  otherwise
    yesno(gcbf,'no');
end
