function name = getcredsfile(inpath)
%Returns the name of the credentials file. If not existing, then return
%empty. Function looks in the default directory. If multiple files are
%existing, the return them in a cell array.nn 
%
%See also READCREDENTIALS

%Einar Heiberg

name = '';

if nargin == 0
  f = dir('*.transfercreds');
else
  f = dir([inpath filesep '*.transfercreds']);
end;

if isempty(f)
  return; 
end;

if length(f)>1
  name = cell(1,length(f));
  for loop = 1:length(f)
    name{loop} = f(loop).name;
  end;
else
  name = f(1).name;
end;