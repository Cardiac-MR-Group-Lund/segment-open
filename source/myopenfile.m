function myopenfile(filename)
% open file in an appropriate program

if ~exist(filename, 'file')
  myfailed(dprintf('Could not find %s',filename));
  return;
end

if ispc
  winopen(filename);
end