function showdicominfo(filename)
%SHOWDICOMINFO shows DICOM tags in a window.
%
%See also DICOMINFO

%Einar Heiberg, 2003-12-22.
global DATA

if nargin==0  
  filename = DATA.Preview.PreviewFile;
end

%Check if not a filename
pos = strfind(filename,'_Callback');
if ~isempty(pos)
  eval(filename); %it is not a filename, it is a function
  return;
end

%Check if copy_Callback 
if isequal(filename,'copy_Callback')
  copy_Callback;
  return;
end

%Temporary switch
stable = true;

if ~exist(filename,'file')
  myfailed(dprintf('Could not find file %s.',filename));
  return;
end

%debugging
%system(sprintf('dcmdump%s "%s" > temp.txt',platformextension,filename));
set(gcf,'pointer','watch');
%info = dicominfo(filename)

%Add this line for "real" DICOM debugging
%system(sprintf('dcmdump.exe "%s"',filename));

try
  if stable
    dinfo = dicominfo(filename);
  else
    dinfo = fastdicominfo(filename);
    if ischar(dinfo) || isequal(dinfo.TransferSyntaxUID,'jpeg') || isequal(dinfo.TransferSyntaxUID,'1.2.840.10008.1.2.4.90')
      dicomconvert(filename); %try to convert on the fly
      dinfo = fastdicominfo(filename);
      if ischar(dinfo)
        error(dinfo)
      end
    end
  end
catch me
  myfailed('Major problems parsing file.');
  mydispexception(me);
  dinfo.Error = 'Could not parse the file using ordninary parser. DICOM?';
end

fields = fieldnames(dinfo);
[~,ind] = sort(fields);
fields = fields(ind);
res = cell(size(fields));

%Find longest tagname
maxlen = 0;
for loop=1:length(fields)
  if length(fields{loop})>maxlen
    maxlen = length(fields{loop});
  end
end

priv = false(length(fields),1);
for loop=1:length(fields)
  if findstr(fields{loop},'Private')
    priv(loop) = true;
  end
  name = [fields{loop} repmat(' ',1,maxlen-length(fields{loop}))];
  f = getfield(dinfo,fields{loop});
  res{loop} = writefield(f,name);
end

res = res([find(not(priv));find(priv)]);
set(gcf,'pointer','arrow');

gui = mygui('dinfo.fig');

DATA.GUI.ShowDICOMInfo = gui;

%Store it
gui.dinfo = dinfo;

% Use system color scheme for figure:
%handles = guihandles(fig);
set(gui.handles.dinfolistbox,'String',res,'fontname',get(0,'FixedWidthFontName'));

%--------------------------------
function res = writefield(f,name)
%--------------------------------
%Helper function to write a field

switch class(f)
  case 'double'
    if numel(f)<=10
      res = sprintf('%s:%s',name,mynum2str(f));
    else
      res = sprintf('%s:[double size %s]',name,num2str(size(f)));
    end    
  case 'single'
    if numel(f)<=10
      res = sprintf('%s:%s',name,mynum2str(f));
    else
      res = sprintf('%s:[single size %s]',name,num2str(size(f)));
    end        
  case 'char'
    res = sprintf('%s:%s',name,f);
  case {'int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
    if numel(f)<5
      res = sprintf('%s:%s',name,sprintf('%d ',f));
    else
      res = sprintf('%s:[%d] ... %s',name,numel(f),sprintf('%d ',f(1:5)));
    end
  case 'struct'
    if length(fieldnames(f))>2
      res = sprintf('%s:[Struct]',name);
    else
      tempstri = '';
      subnames = fieldnames(f);
      for subloop=1:length(subnames)
        if isa(getfield(f,subnames{subloop}),'char')
          tempstri = [tempstri sprintf('%s:%s',subnames{subloop},getfield(f,subnames{subloop})) ' ']; %#ok<AGROW>
        else
          tempstri = [tempstri sprintf('%s:%s',subnames{subloop},class(getfield(f,subnames{subloop}))) ' ']; %#ok<AGROW>
        end
      end
      res = sprintf('%s=>%s',name,tempstri);
    end
  otherwise
    res = sprintf('%s:[%s]',name,class(f));
end
  
%-------------------------
function stri = mynum2str(vec)
%-------------------------
stri = num2str(vec(:)');

%-------------------------
function compare_Callback %#ok<DEFNU>
%-------------------------
global DATA

gui = DATA.GUI.ShowDICOMInfo;

fnames = fieldnames(gui.dinfo);
for loop = 1:length(fnames)
  if isequal(gui.dinfo.(fnames{loop}),gui.olddinfo.(fnames{loop}))
  else
    disp(sprintf('New: %s',writefield(gui.dinfo.(fnames{loop}),fnames{loop}))); %#ok<DSPS>
    disp(sprintf('Old: %s',writefield(gui.olddinfo.(fnames{loop}),fnames{loop}))); %#ok<DSPS>
  end
end

%----------------------
function store_Callback %#ok<DEFNU>
%----------------------
global DATA

gui = DATA.GUI.ShowDICOMInfo;

gui.olddinfo = gui.dinfo;

mymsgbox('Store for future comparison.');

%---------------------
function copy_Callback
%---------------------
global DATA

gui = DATA.GUI.ShowDICOMInfo;

v = get(gui.handles.dinfolistbox,'Value');
cstri = get(gui.handles.dinfolistbox,'String');
stri = cstri{v};

%Cut off after first colon
pos = find(stri==':',1);
if ~isempty(pos)
  stri = stri(pos+1:end);
end

segment('cell2clipboard',{stri});

%----------------------
function close_Callback
%----------------------
%Close showdicominfo GUI
global DATA

try
  DATA.GUI.ShowDICOMInfo = close(DATA.GUI.ShowDICOMInfo);
catch   %#ok<CTCH>
  DATA.GUI.ShowDICOMInfo=[];
  delete(gcbf);
end
