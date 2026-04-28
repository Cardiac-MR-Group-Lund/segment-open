function myworkoff(fignr)
%MYWORKOFF Graphically show that calculation is
%  finished by restoring pointer.

%Einar Heiberg

global DATA %#ok<*GVMIS> 
%resume last pointer

if nargin==0
  fignr = gcf;  
end

if ~isempty(DATA)
  set(fignr,'pointer',DATA.LastPointer);
else
  close(fignr)
  return;
end

if ~isempty(DATA)
  if isequal(DATA.LastPointer,'custom')
    set(fignr,'pointer',DATA.LastPointer,...
      'pointershapecdata',DATA.LastPointerShapeCData);
  end
end
flushlog;
      