function s = myexecutableext
%MYEXECUTABLEEXT Returns proper file extension for exectutables.
%
%See also, MYREQUIREPC, MYMKDIR, MYMOVEFILE, MYCOPYFILE, MYDEL, MYDIR.

%Einar Heiberg

if ispc
  s = '.exe';
else
  temp = mexext;
  s = temp(4:end);
  s = ['_' s];
  if isequal(s, '_maci64')
    s = '_maci';
  end
end