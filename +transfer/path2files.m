function files = path2files(path)
%Return the files in a path

  dd = dir(path);
  files = {};
  
  for i=1:numel(dd) 
    if(dd(i).isdir)
      continue
    end
    if isequal(dd(i).name, 'dicom.cache')
      continue
    end
    files(end+1) = {[path dd(i).name]};
  end
end