function mywarning(stri,fighandle)
%MYWARNING(STRI,HANDLE) Displays a warning message
%  User needs to click OK to discard message.
%
%See also  MYWARNINGNOBLOCK

%%%% If changing this file, please also change mywarning no block, since
%%%% they should be identical except for the uiwait(h); line.
global DATA

%Einar Heiberg

stri = translation.dictionary(stri);

try
  if DATA.Pref.DoNotAsk
    %Just print warning message in window if DoNotAsk mode.
    mydisp(dprintf('Warning:%s\n',stri));
    return;
  end;
catch %#ok<CTCH>
end;

keystroke = popfrombuffer('KeyStroke');

if isempty(keystroke)
  mydisp(dprintf('Warning:%s\n',stri));
  h = warndlg(stri,'Warning:');

  %Adjust horisontally
  if nargin>1
    myadjust(h,fighandle);
  else
    myadjust(h);
  end;
  flushlog;
  uiwait(h); %This is the only difference to mywarningnoblock.
  try
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mywarning');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
  catch %#ok<CTCH>
    %Ignore if failed, likely due to DATA is cleared or not present
  end;
else
  mydisp(dprintf('Warning:%s\n',stri));  
  if isequal(lower(keystroke),'ok')
    pushtobuffer('Warnings',stri);
    macro_helper('put','pushtobuffer(''KeyStroke'',''ok''); %ok from mywarning');
    macro_helper('switchorder'); %We need to store data in buffer before the callback
    return;
  else
    error(sprintf('Expected ''ok'' as keystroke, got %s',keystroke));
  end;
end;
flushlog;
