function showdicominfo(varargin)
%SHOWDICOMINFO shows DICOM tags in a window.
%
%See also DICOMINFO

%#ok<*GVMIS>
%Einar Heiberg, 2003-12-22. 
%Rewritten by Fanny Manefjord 2023-01-25, ticket #2584

global DATA

if nargin==0  
  init(DATA.Preview.PreviewFile);

elseif nargin==1 && ~contains(varargin{1},'Callback','IgnoreCase',true)
  init(varargin{1});
  return
else
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
end

%--------------------------------
function init(filename)
%--------------------------------
%Initialize GUI, input argument is filename of the dicom file.
global DATA

gui = mygui('dinfo.fig');
DATA.GUI.ShowDICOMInfo = gui;

worklist = getdicominfo(filename);
[pathname,file] = fileparts(filename);
gui.startpathname = pathname;
table = gui.handles.dicomtable;
set(table,'Data',worklist);
%set column width as half of the available size 
table.Units = 'pixels';
width = table.Position(3)-31; % subtracting by 31 to account for row header
ncol = size(table.Data, 2);
table.ColumnWidth = num2cell(ones(1,ncol)*width/ncol);
table.Units = 'normalized';
table.ColumnName = {'DICOM tag' file};

%--------------------------------
function worklist = getdicominfo(filename)
%--------------------------------
%Write dicom info to a cell array

worklist = {''};
try
  dinfo = dicominfo(filename);

catch me
  myfailed('Major problems parsing file.');
  mydispexception(me);
  worklist{1,1} = filename;
  worklist{1,2} = 'Could not parse the file using ordninary parser. Is it a DICOM file?';
  return
end

fields = fieldnames(dinfo);
[~,ind] = sort(fields);
fields = fields(ind);

for loop=1:length(fields)
  name = fields{loop}; %take out the dicom tag name

  if loop == 1
    worklist{end,1} = name;
  else
    worklist{end+1,1} = name; %#ok<AGROW> 
  end
  value = dinfo.(name);
  worklist = worklist_helper(value,worklist,0);
end

%Put the private tags at the end
privind = contains(worklist(:,1),'private', IgnoreCase=true);
ind = [find(not(privind));find(privind)];
worklist = worklist(ind,1:2);

%--------------------------------
function worklist = worklist_helper(value, worklist, index)
%--------------------------------
%Helper function to make a worklist of dicom tags

if isstruct(value)
  index = index + 1;
  names = fieldnames(value);
  for structloop = 1 : numel(names)
    arrow = repmat('> ',1,index);
    name = names{structloop};
    str = sprintf(' %s %s',arrow, name);
    worklist{end+1,1} = str; %#ok<AGROW> 
    v = value.(name);
    worklist = worklist_helper(v,worklist,index);
  end
else
  worklist{end,2} = writefield(value);
end

%--------------------------------
function res = writefield(f)
%--------------------------------
%Helper function to write a field

switch class(f)
  case 'double'
    if numel(f)<=10
      res = sprintf('%s',mynum2str(f));
    else
      res = sprintf('[double size %s]',num2str(size(f)));
    end    
  case 'single'
    if numel(f)<=10
      res = sprintf('%s',mynum2str(f));
    else
      res = sprintf('[single size %s]',num2str(size(f)));
    end        
  case 'char'
    if (contains(f,char(10))) %take away newline
      f = strrep(f,char(10),' ');
    end
    if (contains(f,char(13))) 
      f = strrep(f,char(13),' ');
    end
    res = sprintf('%s',f);
  case {'int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
    if numel(f)<5
      res = sprintf('%s',sprintf('%d ',f));
    else
      res = sprintf('[%d] ... %s',numel(f),sprintf('%d ',f(1:5)));
    end
  case 'struct'
    if length(fieldnames(f))>2
      res = sprintf('[Struct]');
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
      res = sprintf('>%s',tempstri);
    end
  otherwise
    res = sprintf('[%s]',class(f));
end
  
%-------------------------
function stri = mynum2str(vec)
%-------------------------
stri = num2str(vec(:)');

%-------------------------
function compare_Callback
%-------------------------
%Opens a browse window to choose a new file and puts the new dicom values next
%to the dicom values of the old file

global DATA

myworkon
gui = DATA.GUI.ShowDICOMInfo;

%select file
[filename{1},filepath{1}] = myuigetfile('*.*','Select dicom-file.',gui.startpathname);
if isequal(filepath{1},0)||isequal(filename,0)
  myfailed('Function aborted. No file selected.',DATA.GUI.Segment);
  myworkoff
  return;
end

file = strcat(char(filepath{1}), char(filename{1}));

firstwl = getlist; %dicom tags of the old file
secondwl = getdicominfo(file); %dicom tags of the new file

alltags = firstwl; %list of all tags

%go through all dicom tags in the new file and add them to the end of the list if they're not alreday there
for tagnumber = 1 : length(secondwl)
  tagname = secondwl{tagnumber,1};
  tagvalue = secondwl{tagnumber,2};
  istagpresent = matches(alltags(:,1),tagname);
  added = 0;
  if ~any(istagpresent) %dicom tag not in alltags, put it at the end
    alltags{end+1,1} = tagname; %#ok<AGROW> 
    alltags{end,3} = tagvalue; %put value in 3rd column
    added = 1;
  else %dicom tag already in alltags, put value in 3rd column
    ind = find(istagpresent);

    %tag is present more than 1 time, probably a nested dicom tag
    if (numel(ind)>1)
      
      % how many parent tags
      numparent = length(strfind(tagname, '>'));
      parentfound = zeros(1,numparent);
      index = tagnumber - 1; %go up to find parents 
      currentarrows = numparent;
      parent = strings(numparent, 1);
      for loop = 1:numparent
        while(~parentfound(loop))
        parent(loop,1) = secondwl{index,1};
        arrows = length(strfind(parent(loop), '>'));
        if (arrows == currentarrows - 1)        
          parentfound(loop) = 1;
          currentarrows = currentarrows - 1;
          index = index - 1;
        else
          index = index - 1;
        end
        end
      end
      
      parent = flip(parent);
      parent(end + 1) = tagname;
      p = matches(alltags(:,1),parent(1)); %start with the first parent
      parentindex = find(p);
      childindex = parentindex + 1; %go down to find next parent/child
      childfound = zeros(numparent,1);
      child = 1;

      if (childindex>length(alltags))
        childfound(end) = 1;
      end
      
      for loop = 1 : length(childfound)
        while (~childfound(loop) && child)
          current = alltags{childindex,1};
          if strcmp(current,parent(loop+1))
            childfound(loop) = 1;
            if (sum(childfound) == length(childfound))
              alltags{childindex,3} = tagvalue;
              added = 1;
            end
           
          else %not the correct child
            childindex = childindex + 1;
            if ~contains(current,'>') %we have moved on and it's not a child anymore
              child = 0;
            end
          end
        end
      end

       if (childindex>length(alltags))
        alltags{childindex,1} = tagname;
        alltags{childindex,3} = tagvalue;
        added = 1;
        continue;
      end

      %parent exists in alltags, but not the child
      %add child in alltags on next index and shift everything down
      if ~added
        newline = {tagname, [], tagvalue};
        alltags = [alltags(1:parentindex,:);newline;alltags(parentindex + 1:end,:)];
      end

    else
    %not present more than 1 time, just add the value in the right place
      alltags{ind,3} = char(secondwl{tagnumber,2});
    end
  end
end

%---Put new worklist (2 dicom files) in figure
table = gui.handles.dicomtable;
set(table,'Data',alltags);
table.Units = 'pixels';
width = table.Position(3)-31; % subtracting by 31 to account for row header
ncol = size(table.Data, 2);
table.ColumnWidth = num2cell(ones(1,ncol)*width/ncol);
table.Units = 'normalized';
columnnames = table.ColumnName;
table.ColumnName = {columnnames{1}; columnnames{2}; filename{1}};

myworkoff

%---------------------
function copy_Callback
%---------------------
%Copies tags and values to the clipboard

global DATA
gui = DATA.GUI.ShowDICOMInfo;
table = gui.handles.dicomtable;
columnnames = table.ColumnName;
worklist = getlist;

str = [];
newline = sprintf('\n');
len = length(worklist);

%add headers in first row
for column = 1 : length(columnnames)
  firstrow = sprintf('%s\t', columnnames{column});
  str = [str firstrow];
end

str = [str newline];

%add all values, tabs between columns and newline before next row
for i = 1:len
  row = sprintf('%s\t', worklist{i,:});
  row(end) = newline;
  str = [str row];  %#ok<AGROW>
end

segment('cell2clipboard',{str});

%---------------------
function worklist = getlist
%---------------------
global DATA

gui = DATA.GUI.ShowDICOMInfo;
table = gui.handles.dicomtable;
worklist = get(table,'Data');%'fontname',get(0,'FixedWidthFontName'));

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
