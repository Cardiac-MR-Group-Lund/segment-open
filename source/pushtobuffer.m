function pushtobuffer(f,o)
%PUSHTOBUFFER(F). Where F is the field to push to. 
%
%See also CLEARBUFFER, POPFROMBUFFER.

%Einar Heiberg

global DATA

if nargin<2
  error('Expected two input arguments.');
end;

buffer = DATA.Buffer.(f);
isacell = iscell(buffer);
if isacell
  buffer = [{o} buffer];
else
  buffer = [o buffer];
end;

%Store it.
DATA.Buffer.(f) = buffer;
