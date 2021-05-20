function pathname = getpreferencespath
%Get path to preferences folders
% global DATA
persistent localprefpath
% if isa(DATA,'maingui')
%   pathname = DATA.getpreferencespath;
%   return
% end
foldername = 'Segment';
if isempty(localprefpath)
  if ispc
    % check if Siemens version
%     s = commandlinehelper('getparameters');
%     if isfield(s,'LogFolder') && (~isempty(s.LogFolder))
%       % Siemens version
%       pathname = s.LogFolder;
%     else
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
          end
        end
      end
%     end
    %Check if subdirectory exists
    dirpath = [pathname filesep foldername]; 
    if not(exist(dirpath,'dir'))
      try
        fprintf('Creating new folder %s\n',dirpath);
        suc = mkdir(dirpath);
        if not(suc)
          myfailed(dprintf('Could not create %s',dirpath));
        end  
      catch
        myfailed(dprintf('Could not create %s',dirpath));
      end
  %     temppath = DATA.SegmentFolder;
  %     cd(pathname);
  %     disp(sprintf('Creating new folder %s%s%s',temp,filesep,foldername));
  %     suc = mkdir(foldername);
  %     cd(temppath);
  %     if not(suc)
  %       myfailed(sprintf('Could not create %s%s%s',temp,filesep,foldername));
  %     end
    end

    %Add Segment to the path
    localprefpath = dirpath;

  else
    localprefpath = pwd;
  end %is pc
end
pathname = localprefpath;