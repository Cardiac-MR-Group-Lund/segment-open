function myfailed(stri,fighandle)
global DATA %#ok<*GVMIS> 
%MYFAILED(STRI,FIGHANDLE)
%  Displays an error message STRI. FIGHANDLE is 
%  an optional alignment indicating alignment. 
%  Difference to errordlg is the alignment and 
%  possibility to hit return to okey the message.
%  Alignment string may be a fig handle or an MYGUI object.
%
%  If DATA.Silent then message is just displayed in the console
%  and exectution continues as if user had pressed ok.
%
%See also MYWARNING, MYMSGBOX, MYADJUST, MYWAITBARSTART.

%Einar Heiberg

if nargin < 1
  myfailed('Expected one input argument.');
  return;
end

%If silent then just "write to log-file" and return

% make sure that the string we display is in English to be able to read and
% understand our log file
try
  logstr = translation.dictionary(stri,'English',DATA.Pref.Language);
catch
  logstr = stri;
end

%For autoloader, add warning to warninglist instead of showing it
try 
  if ~isempty(DATA) && DATA.Autoloader
    fprintf('Error: %s\n',logstr); 
    warningfunctions('addwarning', logstr);
    return;
  end
catch
  %Do nothing as if error then depends on DATA not initialised
end

try
  if ~isempty(DATA) && DATA.Silent
    fprintf('Myfailed: %s\n',logstr); 
    return;
  end
catch
  %Do nothing as if error then depends on DATA not initialised
end

logdisp(['Error: ' logstr])
%fprintf('Error: %s\n',logstr);

stri = translation.dictionary(stri);
keystroke = popfrombuffer('KeyStroke');

if isequal(keystroke,'ok')
  pushtobuffer('Warnings',stri);
  return
end

h = errordlg(stri,dprintf('Operation failed.'),'invisible');
fontsize = 10;
try
  if ~isempty(DATA)
    htext = findobj(h, 'Type', 'Text');  %find text control in dialog
    set(htext,'Color',DATA.GUISettings.ForegroundColor,'FontSize',fontsize);
    if length(stri) > 25
      % get current extent of the text message
      ext = htext.Extent;
      pos = h.Position;
      % adjust width with enough space for the message plus 10 pixels margin
      pos(3) = ext(1)+ext(3);
      h.Position = pos;
      set(htext,'Units','normalized')
    end
  end
catch me
  mydispexception(me)
end
setupicon(h);
  
try
  set(h,'Color',DATA.GUISettings.BackgroundColor);
  kids = h.Children;
  for i = 1:length(kids)
    try set(kids(i),'BackgroundColor',DATA.GUISettings.BackgroundColor);catch, end
    try set(kids(i),'Color',DATA.GUISettings.BackgroundColor);catch, end
    try set(kids(i),'ForegroundColor',DATA.GUISettings.ForegroundColor);catch, end
    try
      kidstyle = kids(i).Style;
      if strcmp(kidstyle,'pushbutton')
        set(kids(i),'Units','normalized','FontSize',fontsize)
      end
    catch
    end
  end
catch
end

if length(stri) > 25
  ps = get(h,'Position');
  set(h,'Position',[ps(1) ps(2) ps(3)*1.2 ps(4)*1.2])
end
try
  if nargin>1
    myadjust(h,fighandle);
  else
    myadjust(h);
  end
catch %#ok<CTCH>
end

try
  if ~isempty(DATA) && DATA.Testing
    DATA.GUI.Dialog = h;
    pause(0.1)
    oktimer = timer('ExecutionMode','singleShot',...
      'TimerFcn', 'clickonbutton', ...
      'StartDelay',0.5,'Name','okbuttonclick');
    start(oktimer)
    pushtobuffer('Warnings',stri);
  end
catch me
  mydispexception(me)
end

set(h,'windowstyle','modal','visible','on');
flushlog;
uiwait(h);
try
  macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from myfailed');
  macro_helper('switchorder'); %We need to store data in buffer before the callback
catch %#ok<CTCH>
end