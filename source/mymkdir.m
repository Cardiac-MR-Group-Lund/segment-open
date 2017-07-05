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
end;

varargout = cell(1,nargout);
if nargout>0
  varargout{1} = sucess;
end;
if nargout>1
  error('No more than one output argument.');
end;

