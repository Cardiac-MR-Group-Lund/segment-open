function f = mydir(pathname,checkfilter)
%MYDIR Extracts directory information.
%. & .. are always first on the list
% thumbs.cache,folders.cache,dicomdir excluded.
% folders comes first

%Einar Heiberg
global DATA

filter = false; %Assume without filter

if nargin<2
  checkfilter = true;
end;

if checkfilter
  try    
    %get(DATA.Preview.Handles); This line disabled filter action
    filter = get(DATA.Preview.Handles.filtercheckbox,'value');
  catch %#ok<CTCH>
    filter = false;
  end;
end;

%DATA might not be correctly set => put inside catch clause
try
  if DATA.Preview.LoadAll
    filter = false;
  end
catch %#ok<CTCH>
  filter = false;
end;
  
if ~filter
  f = dir(pathname);
else
  %--- Use filter(s)
  filtersstri = mygetedit(DATA.Preview.Handles.filteredit);

  %Extract filters
  numfilters = sum(filtersstri==';')+1;
  filters = cell(1,numfilters+1);
  filters{1} = '*.'; %This takes directories
  loop = 2;
  while sum(filtersstri==';')>0
    pos = find(filtersstri==';');
    pos = pos(1);
    filters{loop} = filtersstri(1:(pos-1));
    filtersstri = filtersstri((pos+1):end);
    loop = loop+1;
  end;
  filters{loop} = filtersstri;

  %Set up template
  f = [];

  for loop=1:length(filters)
    f = cat(1,f,dir([pathname filesep filters{loop}]));
  end;

end;

% Windows 10(?) bugfix: Sometimes isdir is not set for . and ..
for floop = 1:length(f)
  fname = f(floop).name;
  if isequal(fname, '.') || isequal(fname, '..')
    f(floop).isdir = true;
  end
end

%Sort so that folders comes first
[~,ind] = sort(-cat(1,f(:).isdir));
f = cat(1,f(ind));

%Remove forbidden files
ind = true(size(f));
for loop=1:length(f)
  if isequal(f(loop).name,'thumbs.cache')
    ind(loop) = false;
  end;
  if isequal(f(loop).name,'folders.cache')
    ind(loop) = false;
  end;
  if isequal(f(loop).name,'dicom.cache')
    ind(loop) = false;
  end;
  if isequal(lower(f(loop).name),'dicomdir')
    ind(loop) = false;
  end;
	if isequal(f(loop).name,'.svn'),
		ind(loop) = false;
	end
end;

% Exclude hidden files on unix, i.e. starting with . but isn't . or ..
if DATA.Pref.HideFilesUnix && isunix
  for loop = 1:length(f)
    if ~isequal(f(loop).name, '.') && ~isequal(f(loop).name, '..')
      thisname = f(loop).name;
      if isequal(thisname(1), '.')
        ind(loop) = false;
      end
    end
  end
end

% Filter forbidden files and hidden files on unix.
f = cat(1,f(ind));