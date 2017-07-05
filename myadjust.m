function myadjust(objecthandle,fighandle)
%MYADJUST objecthandle and fighandle/mygui object.
%Adjust so that a message box gets to the correct screen.
%fighandle is a figure handle

%Einar Heiberg
global DATA

try
  if nargin==1
    fighandle=DATA.GUI.Segment.fig;
  elseif isempty(fighandle)
      fighandle=DATA.GUI.Segment.fig;
  elseif isa(fighandle,'mygui')
      fighandle = fighandle.fig;
  end
catch
  fighandle=[];
end

if not(isempty(fighandle))

  figunits=get(fighandle,'units');
  objectunits=get(objecthandle,'units');
  set(fighandle,'units','pixels');
  set(objecthandle,'units','pixels');

  %get position for gui and the handle h
  figpos=get(fighandle,'position');%figposition [x(g) y(g) width(g) height(g)]
  objectpos=get(objecthandle,'position');%objectposition [x(h) y(h) width(h) height(h)]

  %calculate new position for the handle h
  temp=figpos(1:2)+(figpos(3:4)-objectpos(3:4))/2;%[x y](g) + ([width height](g) -[width height](h))/2
  objectpos(1:2)=temp;

  %set position
  set(objecthandle,'position',objectpos);

  set(fighandle,'units',figunits);
  set(objecthandle,'units',objectunits);
end

