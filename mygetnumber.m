function [z,ok] = mygetnumber(promptstri,titlestri,default,minval,maxval)
%Input dialog function that get a scalar value. 
%Syntax
%  MYGETNUMBER(STRI,DEFAULT,MIN,MAX)
%
%Where:
% - PROMPTSTRI prompt for dialog.
% - TITLESTRI title of dialog box.
% - DEFAULT is default value.
% - <MIN>, optional minimal value.
% - <MAX>, optional maximal value.

%Einar Heiberg

global DATA

if nargin<3
  error('Expected at least three input arguments (prompt,title,default).');
end;

if nargin==4
  error('When minimal value is specified, then also maximal value needs to be specified.');
end;

if DATA.Testing
  testing = DATA.Testing;
else
  testing = false;
end;

if DATA.RecordMacro
  recordmacro = DATA.RecordMacro;
else
  recordmacro = false;
end;

if testing
  %Take from buffer
  z = popfrombuffer('Number');
  if ~isempty(z)
    ok = true;
  else
    error('No number in buffer.');
  end;
else
  %Ask user
  z = [];
  s = inputdlg({promptstri},titlestri,1,{sprintf('%0.9g',default)});
  if isempty(s)
    ok = false;
    return;
  else
    [z,ok] = str2num(s{1}); %#ok<ST2NM>
    if not(ok)
      return;
    end;
    if length(z)>1 || isempty(z)
      ok = false;
      myfailed('Need to be a scalar value.');
      return;
    end;
    if nargin>3
      if z<minval
        myfailed(dprintf('Too small value. [%0.5g %0.5g]',minval,maxval));
        ok = false;
        return;
      end;
      if z>maxval
        myfailed(dprintf('Too large value. [%0.5g %0.5g]',minval,maxval));
        ok = false;
        return;
      end;
    end;
  end;
end;

if recordmacro
  macro_helper('put',sprintf('pushtobuffer(''Number'',[ %0.9g DATA.Buffer.Number]); %%add to buffer',z));
  macro_helper('switchorder'); %We need to store data in buffer before the callback
end;
