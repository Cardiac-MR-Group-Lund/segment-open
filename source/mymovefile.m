function [varargout] = mymovefile(source,destination,mode)
%[STATUS,MESSAGE,MESSAGEID] = MOVEFILE(SOURCE,DESTINATION,MODE) 
%moves the file in source to dest.
%
%Essentially a wrapper for MOVEFILE, kept for naming consistency.
%
%See also MYCOPYFILE, MOVEFILE.

%Einar Heiberg
if nargin < 3
  [status,message] = movefile(source,destination);
else
  [status,message] = movefile(source,destination,mode);
end

varargout = cell(1,nargout);
if nargout>0
  varargout{1} = status;
end
if nargout>1
  varargout{2} = message;
end
