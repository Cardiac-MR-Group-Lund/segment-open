function gantrytilttool(files2load)
%Tool to fix gantry tool stacks problems
%Takes a cell of files and creates a new folder with fixed files. If not
%call with any input argument, then uses information from the file loader
%(openfile).

%Einar Heiberg

global DATA %#ok<*GVMIS> 
 
if nargin==0
  %Check if there is gantrydetectortilt
  try
    files2load = openfile('getfiles2load');
  catch
    files2load = DATA.Preview.Files2Load;
    if isempty(files2load)
      logdisp('Gantry tilt tool: Could not find files to load. Openfile might not be open.');
      return;
    end
  end
end

%Check gantry tilt
logdisp('Gantry tilt tool');

%Safety check to ensure that all files are in the same folder
[pathname,~,~] = fileparts(files2load{1});
dowarn = false;
for loop = 2:length(files2load)
  [thispathname,~,~] = fileparts(files2load{loop});
  if ~isequal(pathname,thispathname)
    dowarn = true;
  end
end
if dowarn
  if yesno('Multiple series in gantry tilt, is not supported. Abort?')
    return
  end
end

%Check a middle file about info
file2check = max(1,round(length(files2load)/2));
dinfo = dicominfo(files2load{file2check});

try
  gantrydetectortilt = dinfo.GantryDetectorTilt;
catch
  gantrydetectortilt = 0;
end
if isequal(gantrydetectortilt,0)
  myfailed('Gantrydetector is already zero, nothing to do.');
  return;
end

%Reserve memory
filenames = cell(1,length(files2load));
slicelocations = zeros(1,length(files2load));
imageorientations = filenames;
imagepositions = filenames;
imcell = filenames;
dinfocell = filenames;
ignored = false(size(filenames));

%Load information
h = mywaitbar(dprintf('Please wait.'),gcf);
for loop = 1:length(files2load)
  try
    filenumber = openfile('getfilenumber',files2load{loop});
    filenames{loop} = sprintf('%d',filenumber);
    if length(filenames{loop})>10
      filenames{loop} = filenames{loop}((end-4):end);
    end
    dinfo = dicominfo(files2load{loop});
    dinfocell{loop} = dinfo;
    imcell{loop} = dicomread(dinfo);

    %Extract position
    try
      imageorientations{loop} = dinfo.ImageOrientationPatient;
    catch
      imageorientations{loop} = dinfo.ImageOrientation;
    end

    try
      imagepositions{loop} = dinfo.ImagePositionPatient;
    catch
      imagepositions{loop} = dinfo.ImagePosition;
    end
  catch
    %Something failed
    ignored(loop) = true;
    logdisp(sprintf('Ignored %s',files2load{loop}));
  end

  if isequal(rem(loop,10),0)
    set(h,loop/length(files2load));
  end

end
close(h); %waitbar

%Remove ignored (if any)
files2load = files2load(~ignored);
slicelocations = slicelocations(~ignored);
imagepositions = imagepositions(~ignored);
imageorientations = imageorientations(~ignored);
imcell = imcell(~ignored);
dinfocell = dinfocell(~ignored);

%Compute slicelocations pass 1
for loop = 1:length(files2load)
  zdir = cross(imageorientations{loop}(1:3),imageorientations{loop}(4:6));
  slicelocations(loop) = sum(imagepositions{loop}.*zdir);
end

%Find true-zdir
[~,ind] = sort(slicelocations);
newzdir = imagepositions{ind(end)}-imagepositions{ind(1)};
newzdir = newzdir./(sqrt(sum(newzdir.*newzdir)));

%Compute slicelocations pass 2, this time with zdir based on block
for loop = 1:length(files2load)
  slicelocations(loop) = sum(imagepositions{loop}.*newzdir);
end

[~,ind] = sort(slicelocations);
minslicelocation = min(slicelocations);

%Find gantrydetectortilt angle. Note we should take this as according to
%DICOM specs the code in DICOM files may only be approximate....
gantrydetectortiltfromdicom = gantrydetectortilt;
gantrydetectortilt = acos(sum(newzdir.*zdir))/pi*180;

s = []; %reset
s(1).Field = 'GantryDetectorTilt';
s(1).Label = ['GantryDetectorTilt [' sprintf('%0.5g',gantrydetectortiltfromdicom) ' ' dprintf('in DICOM file') ']']; 
s(1).Default = gantrydetectortilt;
[outs,ok] = myinputstruct(s,'GantryDetectorTilt',10); %10 is 10 characters minimum

if ~ok
  return
end
gantrydetectortilt = outs.GantryDetectorTilt;

%find new folder name
filename = files2load{1};
pos = find(filename==filesep);
if length(pos)<2
  myfailed('Failed.');
  return
end
basepathname = filename(1:(pos(end-1)-1));
foldername = filename((pos(end-1)+1):(pos(end)-1));
newfoldername = [foldername '-gantrymodified-' datestr(now,'HHMMSS')];

%Create directory
mymkdir([basepathname filesep newfoldername]);

%find start positon 
startposition = imagepositions{ind(1)};
imageorientation = imageorientations{ind(1)};

%pre compute tangens
tang = tan(gantrydetectortilt/180*pi); %tangens gantrytilt

%Comput how much the box should be translated
fixtranslation =  -tang*(max(slicelocations)-min(slicelocations))/2;

%Set up struct to move images
T = affine2d(eye(3));

imageorientation = imageorientation(:); %Ensure row vector

%Loop over files
h = mywaitbar(dprintf('Please wait.'),gcf); 
createmode = 'copy';

for loop = 1:length(files2load)
  
  newfilename = sprintf('out%06d.dcm',loop-1);
  outfilename = [basepathname filesep newfoldername filesep newfilename];

  %Compute total distance this slice has from the first slice
  dist = slicelocations(loop)-minslicelocation;
  newposition = startposition+dist*zdir; %obs zdir as slices are orthogonal to zdir
  newposition = newposition-fixtranslation*(imageorientation(4:6)); %Also add the fix translation we move
  newslicelocation = sum(newposition.*zdir);
  
  %Retrieve dinfo
  dinfo = dinfocell{loop};
  im = imcell{loop};
  
  %Assign new position
  if isfield(dinfo,'ImagePosition')
    dinfo.ImagePosition = newposition(:);
    dinfo.ImageOrientation = imageorientation; %ensure all are the same, take first
  else
    dinfo.ImagePositionPatient = newposition(:);
    dinfo.ImageOrientationPatient = imageorientation;
  end
  if isfield(dinfo,'SliceLocation')
    dinfo.SliceLocation = newslicelocation;
  end  
  
  %Compute distance to move image
  moveimagemm = dist*tang+fixtranslation;
  moveimage = moveimagemm/dinfo.PixelSpacing(2);
  
  fillvalues = min(im(:));
  
  %Set translation
  R = imref2d(size(im),[1 size(im,2)],[1 size(im,1)]+moveimage);
  movedim = imwarp(im,T,'OutputView',R,'FillValues',fillvalues); %translate image
          
  %Set zero gantry tilt
  dinfo.GantryDetectorTilt = 0;

  %write
  didfailed = false;
  try
    dicomwrite(movedim,outfilename,dinfo,'WritePrivate',true,'CreateMode',createmode); %First try with createmode
  catch
    if isequal(createmode,'copy') %See if we should try create instead
      createmode = 'create';
      try
        dicomwrite(movedim,outfilename,dinfo,'WritePrivate',true,'CreateMode',createmode);
      catch
        didfailed = true;
      end
    else
      didfailed = true;
    end
  end

  if didfailed
     myfailed('Failed to write file.')
     close(h); %waitbar
     return
  end

  set(h,loop/length(files2load)); %waitbar

end

close(h); %waitbar

%--- Try to delete thumbs file
try
  thumbsfile = [DATA.Preview.PathName filesep 'thumbs.cache'];
  if exist(thumbsfile,'file')
    delete(thumbsfile);
    logdisp('Thumbs file deleted');
  end
catch me
  mydispexception(me)
  myfailed('Could not delete thumbnail file');
  logdisp('Could not delete thumbnail file');
end

%Refresh for filereader
try
  openfile('refresh_Callback');
catch
end

%Refresh for the patient database
try
  patientdatabase('studylistbox_Callback')
catch
end