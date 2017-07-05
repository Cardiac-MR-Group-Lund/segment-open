function pathname = getpreferencespath
%Get path to preferences folders
global DATA
if isa(DATA,'maingui')
  pathname = DATA.getpreferencespath;
  return
end
foldername = 'Segment';

if ispc
  %Is Windows, check for better
  temp = getenv('APPDATA');
  if not(isempty(temp))
    pathname = temp;
  else
    temp = getenv('USERPROFILE');
    if not(isempty(temp))
      pathname = temp;
    else
      temp = getenv('HOMEPATH');
      if not(isempty(temp))
        pathname = temp;
      end;
    end;
  end;
  
  %Check if subdirectory exists
  if not(exist([temp filesep foldername],'dir'))
    temppath = pwd;
    cd(pathname);
    disp(sprintf('Creating new folder %s%s%s',temp,filesep,foldername));
    suc = mkdir(foldername);
    cd(temppath);
    if not(suc)
      myfailed(sprintf('Could not create %s%s%s',temp,filesep,foldername));
    end;
  end;
  
  %Add Segment to the path
  pathname= [pathname filesep foldername];
  
else
  pathname = pwd;
end; %is pc