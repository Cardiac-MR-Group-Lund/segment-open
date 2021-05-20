function [varargout] = mycopyfile(source,dest)
%MYCOPYFILE Copies a file from source to dest
%
%MYCOPYFILE(SOURCE,DEST)
%
%Should be platform independent and faster than making a system
%call. If really large files are used then do a system call instead.
%Essentially a wrapper to COPYFILE, kept for naming consistency.
%
%See also MYDEL, MYMKDIR.

[sucess,message] = copyfile(source,dest);

varargout = cell(1,nargout);
if nargout>0 
  varargout{1} = sucess;
end
if nargout>1
  varargout{2} = message;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Alternative code Einar Heiberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Einar Heiberg
% sucess = false;
% message = '';
% 
% if nargin<3
%   siz = inf;
% end;
% 
% if nargin<2
%   error('Expected two input arguments.');
%   return;
% end;
% 
% fidi = fopen(source);
% if fidi==-1
%   message = sprintf('Could not open %s for input.',source);
%   return;
% end;
% 
% fido = fopen(dest,'w','n');
% if fido==-1
%   message = sprintf('Could not open %s for output.',source);
%   return;
% end;
% 
% [temp,count] = fread(fidi,siz,'uint8=>uint8');
% if isinf(siz)
%   siz = count;
% end;
% if count~=siz
%   message = sprintf('Could not read whole file:%s.',source);
%   fclose(fidi);
%   fclose(fido);
%   return;
% end;
% count = fwrite(fido,temp);
% if count~=siz
%   message = sprintf('Could not write whole file:%s.',dest);
%   fclose(fidi);
%   fclose(fido);
%   return;
% end;  
% 
% fclose(fidi);
% fclose(fido);
% sucess=true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alternative code Written by Johan Ugander %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %---------------------------
% function mycopy(source,dest)
% %---------------------------
% if ispc
%   dos(sprintf('copy "%s" "%s"',source,dest));
% else
%   system(sprintf('cp %s %s',source,dest));
% end
% disp(sprintf('copy "%s" "%s"',source,dest));