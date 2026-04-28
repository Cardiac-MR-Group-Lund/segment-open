function myworkon(fignr)
%MYWORKON Graphically indicate that Segment is 
%  MYWORKON(<FIG>);
%  busy by showing watch pointer.

%Einar Heiberg

global DATA %#ok<*GVMIS> 

if nargin==0
  fignr = gcf;  
end

if ~isempty(DATA)
  temp = get(fignr,'pointer'); %Was DATA.imagefig
else
  return;
end

if isequal(temp,'watch')
  return;
else
  if ~isempty(DATA)
    DATA.LastPointer = temp;
    DATA.LastPointerShapeCData = get(fignr,'pointershapecdata');
  end
  set(fignr,'pointer','watch');
end
drawnow limitrate;