function myfailed(stri,fighandle)
global DATA
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

if nargin<1
  myfailed('Expected one input argument.');
  return;
end

%If silent then just "write to log-file" and return
try
  if DATA.Silent
    disp(sprintf('Error: %s',stri)); %#ok<DSPS>
    return;
  end
catch
  %Do nothing as if error then depends on DATA not initialised
end

stri = translation.dictionary(stri);

mydisp(dprintf('Error:%s\n',stri));

keystroke = popfrombuffer('KeyStroke');

if isequal(keystroke,'ok')
  pushtobuffer('Warnings',stri);
  return
end

h = errordlg(stri,dprintf('Operation failed.'),'invisible');
try
  if ~isempty(DATA)
    htext = findobj(h, 'Type', 'Text');  %find text control in dialog
    set(htext,'Color',DATA.GUISettings.ForegroundColor,'FontSize',9,'Units','normalized');
  end
catch me
  mydispexception(me)
end
setupicon(h);
  
try
  set(h,'Color',DATA.GUISettings.BackgroundColor);
  kids = h.Children;
  for i=1:length(kids)
    try set(kids(i),'BackgroundColor',DATA.GUISettings.BackgroundColor);catch, end
    try set(kids(i),'Color',DATA.GUISettings.BackgroundColor);catch, end
    try set(kids(i),'ForegroundColor',DATA.GUISettings.ForegroundColor);catch, end
    try
      kidstyle = kids(i).Style;
      if strcmp(kidstyle,'pushbutton')
        set(kids(i),'Units','normalized','FontSize',9)
      end
    catch
    end
  end
catch
end
if length(stri) > 40
  ps = get(h,'Position');
  set(h,'Position',[ps(1) ps(2) ps(3)*1.2 ps(4)])
end
try
  if nargin>1
    myadjust(h,fighandle);
  else
    myadjust(h);
  end
catch %#ok<CTCH>
end

set(h,'windowstyle','modal');
set(h,'visible','on');
flushlog;
uiwait(h);
try
  macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from myfailed');
  macro_helper('switchorder'); %We need to store data in buffer before the callback
catch %#ok<CTCH>
end