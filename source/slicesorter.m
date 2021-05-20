function slicesorter
%Sorts a folder of mixed files so that they are sorted into subfolders of
%according to slicelocation. 

%Einar Heiberg

global DATA

%If Segment is running
try
  pathname = DATA.Pref.datapath;
catch %#ok<CTCH>
  pathname = DATA.SegmentFolder;
end

%If loader is running
try
  if isa(DATA.GUI.OpenFile,'mygui')
    f = DATA.Preview.FileList;
    v = get(DATA.Preview.Handles.pathlistbox,'value');
    f = f(v);
    if (length(f) == 1) && f(1).isdir
      pathname = [DATA.Preview.PathName filesep f(1).name];
    end
  end
catch %#ok<CTCH>
end

[pathname,ok] = myuigetdir(pathname,'Select folder to sort');
if ~ok
%   myfailed('Aborted');
  return;
end

outpathname = [pathname '-sorted'];
ok = mymkdir(outpathname);
if ~ok
  myfailed(dprintf('Could not create output path %s',outpathname));
  return;  
end

f = dir(pathname);

h = mywaitbarstart(length(f),'Please wait, processing files.');
for loop = 1:length(f)
  if ~f(loop).isdir
    %Not a folder => filename
    filename = [pathname filesep f(loop).name];
    if guessifdicom(filename,1) %1=heuristic mode
      
      info = fastdicominfo(filename);
      if ischar(info) || isequal(info.TransferSyntaxUID,'jpeg') || isequal(info.TransferSyntaxUID,'1.2.840.10008.1.2.4.90')
        dicomconvert(filename); %try to convert on the fly
        info = fastdicominfo(filename);
        if ischar(info)
          error(info);
        end
      end
      %Use this instead of sliceLocation which sometimes is buggy
      normal = cross(info.ImageOrientation(1:3),info.ImageOrientation(4:6));
      sl = sum(normal.*info.ImagePosition);
      
      %z = sum(self.getnormal.*self.getposition); %Use this instead to look for the tag. It is buggy under certain conditions, this is more stable.

      if ~isempty(info.SeriesDescription)
        slicestring = [strtrim(info.SeriesDescription) '_'];    
      else
        slicestring = '';
      end
      if sl<0
        slicestring = [slicestring '-']; %#ok<AGROW>
      else
        slicestring = [slicestring '']; %#ok<AGROW>
      end     
      slicestring = [slicestring sprintf('%03d',floor(abs(sl)))]; %#ok<AGROW>
      slicestring = [slicestring '.']; %#ok<AGROW>
      slicestring = [slicestring sprintf('%03d',round(1000*rem(abs(sl),1)))]; %#ok<AGROW>          
      
      newfolder = [outpathname filesep slicestring];
      if ~exist(newfolder,'dir')
        mkdir(newfolder);
      end
      [status,msg] = copyfile(filename,[newfolder filesep f(loop).name]);
      if ~isequal(status,1)
        disp(sprintf('Failed moving file %s. Error message %s',f(loop).name,msg)); 
      end
    end
  end
  h = mywaitbarupdate(h);
end
mywaitbarclose(h);
