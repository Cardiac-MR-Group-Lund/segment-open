function myfailed(stri,fighandle)
%MYFAILED(STRI,FIGHANDLE)
%  Displays an error message STRI. FIGHANDLE is 
%  an optional alignment indicating alignment. 
%  Difference to errordlg is the alignment and 
%  possibility to hit return to okey the message.
%  Alignment string may be a fig handle or an MYGUI object.
%
%See also MYWARNING, MYMSGBOX, MYADJUST, MYWAITBARSTART.

%Einar Heiberg

if nargin<1
  myfailed('Expected one input argument.');
  return;
end;

stri = translation.dictionary(stri);

mydisp(dprintf('Error:%s\n',stri));

keystroke = popfrombuffer('KeyStroke');

if isequal(keystroke,'ok')
  pushtobuffer('Warnings',stri);
  return
end

h = errordlg(stri,translation.dictionary('Operation failed.'));

try
  if nargin>1
    myadjust(h,fighandle);
  else
    myadjust(h);
  end;
catch %#ok<CTCH>
end;

set(h,'windowstyle','modal');
flushlog;
uiwait(h);
try
  macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from myfailed');
  macro_helper('switchorder'); %We need to store data in buffer before the callback
catch %#ok<CTCH>
end;