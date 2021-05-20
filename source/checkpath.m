function checkpath(pathname)
%Ensures that the pathname for Segment is on matlab's path. This is utilized in
%initialize to ensure that the Segment root directory in on the search
%path.

%Einar Heiberg

global DATA

if isdeployed()
 
  %Check if in correct path
  if iscorrectpath(pathname)
    return;
  end
  
  %--- Ok. Not correct look in program files
  programfilesdir = getenv('programfiles');
  
  programnamenospace = DATA.ProgramName;
  programnamenospace = programnamenospace(programnamenospace~=' ');
  
  softwarepath = [programfilesdir filesep 'Medviso' filesep programnamenospace filesep 'application'];
  if iscorrectpath(softwarepath)
    %Found it change path to this folder
    cd(softwarepath);
    return
  end
  
  myfailed('Software started from another folder than where it was installed. Please adjust shortcut path.');
  
else
  
  if nargin<1
    [pathname] = fileparts(which('segment_main'));
  end
  
  %get current path
  p = path;
  
  if ~contains(p,pathname)
    %Add to the path;
    addpath(pathname);
    if exist([pathname filesep 'tensorarray'],'dir')
      addpath([pathname filesep 'tensorarray']); %Required by Vortex methods
      addpath([pathname filesep 'tensorarray' filesep 'VortexMethods']); %Required by Vortex methods
    end
  end
  
end

%-----------------------------------------
function correct = iscorrectpath(pathname)
%-----------------------------------------
%Check if we are in the correct folder
global DATA

correct = false;

%--- Check if correct folder by finding ProgramName.exe
namenospace = DATA.ProgramName;
namenospace = namenospace(namenospace~=' ');
wantfile = [namenospace '.exe'];
f = dir(pathname);
foundit = false;
for loop = 1:length(f)
  if ~f(loop).isdir
    if strcmp(f(loop).name,wantfile)
      foundit = true;
    end
  end
end

if foundit
  %We seem to be in the correct folder
  correct = true;
end
