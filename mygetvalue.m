function v = mygetvalue(h)
%Get value from a slider handle. Used instead of normal get to also be able
%to use macros in Segment.

%Einar Heiberg
global DATA

if ~DATA.Testing
  v = get(h,'value');
else
  if isempty(DATA.Buffer.SliderValue)
    v = get(h,'value');
  else
    v = popfrombuffer('SliderValue');
    set(h,'value',v); %Store value
  end;  
end;

if DATA.RecordMacro
  stri = get(h,'tag');
  namestri = inputname(1);
  macro_helper('put',sprintf('pushtobuffer(''SliderValue'',%0.5g); %% %s %s',v,stri,namestri));
  macro_helper('switchorder'); %We need to position the slider data before the callback
end;