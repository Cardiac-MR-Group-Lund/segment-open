function [varargout] = mydel(source)
%MYDEL Delete file from disk
%
%[OK] = MYDEL(filename)
%
%Should be platform independent.
%
%See also MYCOPY, RECYCLE.

%Einar Heiberg

ok = true;
message = '';

if ispc && not(isempty(findstr(source,'*'))) %using wildcards
 [ok,message]= dos(sprintf('del "%s"',source));
 if isempty(message),
   ok=true;
 end
elseif isunix && not(isempty(findstr(source,'*')))
 [ok,message]= system(sprintf('rm "%s"',source));
else
  if exist(source,'file')
    delete(source);
  else
    ok = false;
    message = 'File not found.';
  end
end

varargout = cell(1,nargout);
if nargout>0
  varargout{1} = ok;
end
if nargout>1
  varargout{2} = message;
end

%Alternative code:
%if ispc
%  dos(sprintf('del "%s"',source));
%else
%  system(sprintf('rm "%s"',source));
%end