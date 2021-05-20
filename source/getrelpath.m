function relpath = getrelpath(basepath,path)
%Return relative path
%basepath is the basepath, and path is the pathname that should be
%shortened to a relative path

%Einar Heiberg

pos = length(basepath);
extractpath = path(1:pos);
if ~isequal(path(pos+1),filesep)   
  %Extra save attempt if last character is fileseparator.
  if isequal(path(pos),filesep)
    pos=pos-1;
  else
    error(sprintf('Expected file separator at this position. Basepath:%s Path:%s',basepath,path)); %#ok<SPERR>
  end  
end
relpath = path(pos+2:end);
if ~isequal(extractpath,basepath)
  error(sprintf('Basepath and the begining of relpath should be the same. Should not occur. Basepath:%s Path:%s',basepath,path)); %#ok<SPERR>
end
