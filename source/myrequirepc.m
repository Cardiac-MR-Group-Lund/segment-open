function z = myrequirepc
%Z = MYREQUIREPC Returns true if PC platform. 
%  If not PC platform, then an error messsage is
%  displayed. Used to block certain functions 
%  that is not implemented for other platforms 
%  than PC.

%Einar Heiberg
global DATA

if ispc
  z = true;
  return;
else
  z = false;
  myfailed('This functionality is not yet implemented for Linux/Mac OS.',DATA.GUI.Segment);
  return;
end;