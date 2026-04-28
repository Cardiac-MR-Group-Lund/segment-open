%-------------------------------
function writelog(message,varargin) % skiptranslation
%-------------------------------
%Store a message to the log file.
%Syntax:
%  log('open',pathname)
%  log(message)
%  log('close',pathname)

persistent filename

if nargin > 1
  filename = varargin{1};
end
% 'a'     open or create file for writing; append data to end of file
logfid = fopen(filename,'at');
  
%Display and write message
disp(message);
try
  fprintf(logfid,'%s\n',message);
catch me
  dispexception(me)
end
fclose(logfid);
