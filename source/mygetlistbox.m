function v = mygetlistbox(h)
%Get value from a listbox handle. Used instead of normal get(handles,'Value') 
%to facilitate writing testing scripts in the softwares.

%Einar Heiberg

global DATA %#ok<*GVMIS> 

try
if ~DATA.Testing
  if ~isempty(DATA.Buffer.ListboxValue) %used for batch export
    v = popfrombuffer('ListboxValue');
    set(h,'value',v); %Store the value
  else
    v = get(h,'value');
  end
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