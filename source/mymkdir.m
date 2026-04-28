function [varargout] = mymkdir(d)
%MYMKDIR(NEWDIR) 
%[SUCESS] = MYMKDIR(NEWDIR)
%
%Make directory, difference to MKDIR is that 
%is does not issue an error message
%if the folder already exists.
%
%See also MYCOPY, MKDIR

%Einar Heiberg

sucess = true;
if ~exist(d,'dir')
  sucess = mkdir(d);
end

if ~sucess   %try to elevate in case no admin rights
  commandstri = sprintf('cmd.exe /c mkdir "%s"',d);
   [status,result] = mysystemadmin(commandstri);
  
    if isequal(status,0)
      %Ok check that it is really there
      if exist(d,'dir')
        sucess = true;
      else
        sucess = false;
      end
    else
      sucess = false;
    end
end

varargout = cell(1,nargout);
if nargout>0
  varargout{1} = sucess;
end
if nargout>1
  error('No more than one output argument.');
end

