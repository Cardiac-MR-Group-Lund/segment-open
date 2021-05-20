function v = mygetlistbox(h)
%Get value from a listbox handle. Used instead of normal get to also be able
%to use macros in Segment.

%Einar Heiberg
global DATA

try
if ~DATA.Testing
  v = get(h,'value');
else
  if isempty(DATA.Buffer.ListboxValue)
    v = get(h,'value');
  else
    v = popfrombuffer('ListboxValue');
    set(h,'value',v); %Store the value
  end
end

if DATA.RecordMacro
  stri = get(h,'tag');
  macro_helper('put',sprintf('pushtobuffer(''ListboxValue'',%0.5g); %% %s',v,stri));
  macro_helper('switchorder'); %We need to position the slider data before the callback
end
catch  %#ok<CTCH>
  %This catch is to catch when DATA struct is not initialized
  v = get(h,'value');
end