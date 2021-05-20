function myworkon(fignr)
%MYWORKON Graphically indicate that Segment is 
%  MYWORKON(<FIG>);
%  busy by showing watch pointer.

%Einar Heiberg
global DATA

if ~isempty(DATA)
  temp = get(DATA.imagefig,'pointer');
else
  return;
end

if nargin==0
  fignr = gcf;  
end

if isequal(temp,'watch')
  return;
else
  if ~isempty(DATA)
    DATA.LastPointer = temp;
    DATA.LastPointerShapeCData = get(DATA.imagefig,'pointershapecdata');
  end
  set(fignr,'pointer','watch');
end
drawnow;