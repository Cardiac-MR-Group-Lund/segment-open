function [name,s,id,cmd,args] = decodemessage(stri)
%DECODEMESSAGE Decodes a client/server string
%
%[name,s,id,cmd,arg] = DECODEMESSAGE(stri);
%
%See also ENCODEMESSAGE, MYCLIENTSERVER.

%Einar Heiberg

%Note this function should not be changed unless specifically checked with
%CVQ. Please ask Einar before modifying this file.

%Syntax
%|name|s|id|cmd|arg|

%Default output
name = '';
s = -1;
id = 0;
cmd = '';
args = '';

if isempty(stri)
  args = 'Empty';
  return;
end;

log = stri=='|';
numbars = sum(log);
if numbars<6
  args = 'Inconsitent message.';
  return;
end;

pos = find(log);
if pos(1)~=1
  args = 'Should start with |';
  return;
end;

%if ~isequal(pos(end),length(stri))
%  args = 'Should end with |';
%  return;
%end;
  
%Extract from string
name = stri((pos(1)+1):(pos(2)-1));
s = stri((pos(2)+1):(pos(3)-1));
id = stri((pos(3)+1):(pos(4)-1));
cmd = stri((pos(4)+1):(pos(5)-1));
args = stri((pos(5)+1):(pos(6)-1));

%Convert to double for s and id
s = str2double(s);
if isnan(s)
  args = 'Invalid status';
  return;
end;
id = str2double(id);
if isnan(id)
  args = 'Invalid id';
  return;
end;

%Make it nicer (could be Empty string 1-by-0)
if isempty(cmd)
  cmd = '';
end;
if isempty(args)
  args = '';
end;
