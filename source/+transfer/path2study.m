function r = path2study(path)
%Returns a cell array with the paths to studies
  tt = dir(path);
  r = {};
  for i=1:numel(tt)
    if isequal(tt(i).name, '.')
      continue
    end
    if isequal(tt(i).name, '..')
      continue
    end
    if isequal(tt(i).name, '.svn')
      continue
    end
    if not(tt(i).isdir)
      continue
    end
    r{end+1} = transfer.path2files([path filesep tt(i).name filesep]);
  end
end