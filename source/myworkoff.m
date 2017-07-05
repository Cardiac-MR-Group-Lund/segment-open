function myworkoff
%MYWORKOFF Graphically show that calculation is
%  finished by restoring pointer.

%Einar Heiberg

global DATA
%resume last pointer

if ~isempty(DATA)
  set(gcf,'pointer',DATA.LastPointer);
else
  return;
end;

if ~isempty(DATA)
  if isequal(DATA.LastPointer,'custom')
    set(DATA.imagefig,'pointer',DATA.LastPointer,...
      'pointershapecdata',DATA.LastPointerShapeCData);
  end;
end;
flushlog;
      