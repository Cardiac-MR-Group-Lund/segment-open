function stri = encodemessage(name,s,id,cmd,args)
%ENCODEMESSAGE Encodes a string from a status, id, cmd, and arg.
%
%STRI = ENCODEMESSAGE(NAME,S,ID,CMD,ARG)
%
%NAME - Name of client/sever pair.
%S    - Status (integer, 0 is ok).
%ID   - ID of message (integer). Currently unused, unless for CVQ.
%CMD  - Command to send.
%ARGS - Optional args to send.
%
%See also DECODEMESSAGE, MYCLIENTSERVER.

%NOTE: This function should not be modified unless it is checked with CVQ.
%Ask Einar for furhter details.

%Einar Heiberg

if nargin<3
  error('Requires at least three input arguments.');
end;
if nargin<4
  cmd = '';
end;
if nargin<5
  args = '';
end;

stri = sprintf('|%s|%d|%d|%s|%s|\n',name,s,id,cmd,args);

%%---------------------------------
%function stri = clearnewline(stri)
%---------------------------------
%Removes all newline chars from the string
%stri(stri==sprintf('\n')) = ' ';