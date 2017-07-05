function myworkon
%MYWORKON Graphically indicate that Segment is 
%  busy by showing watch pointer.

%Einar Heiberg
global DATA

if ~isempty(DATA)
  temp = get(DATA.imagefig,'pointer');
else
  return;
end;

if isequal(temp,'watch')
  return;
else
  if ~isempty(DATA)
    DATA.LastPointer = temp;
    DATA.LastPointerShapeCData = get(DATA.imagefig,'pointershapecdata');
  end;
  set(gcf,'pointer','watch');
end;
drawnow;