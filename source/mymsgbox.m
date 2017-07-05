function mymsgbox(stri,title,fighandle)
%MYMSGBOX Helper function, works as MSGBOX
%
%MYMSGBOX(STRI,TITLE,[ALIGNMENT])
%  STRI String to display
%  TITLE Title of messagebox
%  FIGHANDLE Optional figure handle indicating alignment
%
%See also MYWARNING, MYMSGBOX, MYADJUST, MYWAITBARSTART.

%Einar Heiberg

global DATA
if nargin<3
  fighandle=[];
end
if nargin<2
  title = '';
  fighandle=[];
end;
  
stri = translation.dictionary(stri);
try 
  if DATA.Pref.DoNotAsk
    %If do not ask then also do not display
    mydisp(dprintf('Message: %s\n',stri));
    return;
  end;
catch %#ok<CTCH>
end;

try
  mydisp(dprintf('Message: %s\n',stri));
catch %#ok<CTCH>
  disp(sprintf('Message: %s\n',stri));
end;

keystroke = popfrombuffer('KeyStroke');

if isempty(keystroke)
  h = msgbox(stri,title);
  myadjust(h,fighandle);
  set(h,'windowstyle','modal');
  uiwait(h);
  try
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mymsgbox');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
  catch %#ok<CTCH>
    %Do nothing if this fails. Most likely DATA variable is not
    %initialized.
  end;
else
  if isequal(lower(keystroke),'ok')
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mymsgbox');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
    return;
  else
    error('Expected ''ok'' as keystroke');
  end;
end;
flushlog;
