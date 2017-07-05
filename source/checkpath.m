function checkpath(pathname)
%Ensures that the pathname for Segment is on matlab's path. This is utilized in
%initialize to ensure that the Segment root directory in on the search
%path.

%Einar Heiberg
if nargin<1
  [pathname] = fileparts(which('segment_main'));
end;

%get current path
p = path;

if isempty(strfind(p,pathname))
  %Add to the path;
  addpath(pathname);
  if exist([pathname filesep 'tensorarray'],'dir')
    addpath([pathname filesep 'tensorarray']); %Required by Vortex methods
    addpath([pathname filesep 'tensorarray' filesep 'VortexMethods']); %Required by Vortex methods
  end
end;
