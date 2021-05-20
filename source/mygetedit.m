function stri = mygetedit(h)
%Get string from an edit box. Use this instead of get(h,'String') to make
%macro evaluation possible.

%Einar Heiberg
global DATA

if ~DATA.Testing
  stri = get(h,'string');
else
  if isempty(DATA.Buffer.EditString)
    stri = get(h,'string');
  else
    stri = popfrombuffer('EditString');
    set(h,'String',stri); %Store value
  end  
end

if DATA.RecordMacro
  tagstri = get(h,'tag');
  macro_helper('put',sprintf('pushtobuffer(''EditString'',''%s''); %% %s',stri,tagstri));
  macro_helper('switchorder'); %We need to position the slider data before the callback
end
