function mywarningnoblock(stri,fighandle)
%MYWARNING(STRI,HANDLE) Displays a warning message
%  User needs to click OK to discard message.
%
%  If DATA.Silent then message is just displayed in the console
%  and exectution continues as if user had pressed ok.
%
%See also  MYWARNING

global DATA

%Einar Heiberg

try
  if DATA.Pref.DoNotAsk || DATA.Silent
    %Just print warning message in window if DoNotAsk mode.
    mydisp(dprintf('Warning:%s\n',stri));
    return;
  end
catch %#ok<CTCH>
end

mydisp(dprintf('Warning:%s\n',stri));
h = warndlg(stri,'Warning:');

%Adjust horisontally
if nargin>1
  myadjust(h,fighandle);
else
  myadjust(h);
end

%uiwait(h); %This is the only difference to mywarning. If changing
%mywarning, please also update this file.
