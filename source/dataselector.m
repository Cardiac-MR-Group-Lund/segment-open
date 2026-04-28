function varargout = dataselector(varargin)
%Functions to show data as thumbnails and able to select/unselect, datatype can
%be screenshots or thumbnails of SET struct
%#ok<*GVMIS>

if nargin==0
  varargin = {'init','thumbnail'};
end
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%----------------------------------
function varargout = init(datatype)
%----------------------------------
% open figure and adjust according to the data type

global DATA

%Init
DATA.GUI.SeriesSelector = mygui('dataselector.fig');
gui = DATA.GUI.SeriesSelector;
switch datatype
  case 'screenshot'
    [datastructure, filedir] = getscreenshots;
    cm = uicontextmenu(gui.fig);
    gui.handles.deleteuimenu = uimenu(cm,'Text',dprintf('Delete'),'Callback',@deletefile);
    gui.fig.ContextMenu = cm;
    instr01 = dprintf('Left click to select/unselect screenshots to include into report.');
    instr02 = dprintf('Right click to delete screenshots from disc.');
    instructstr = sprintf('%s\n\n%s',instr01,instr02);
    buttonstr = dprintf('Include');
    figname = dprintf('Select screenshot');
  case 'thumbnail'
    % this feature is not implemented yet
    filedir = [];
    instructstr = '';
    buttonstr = dprintf('Close stacks');
    figname = dprintf('Close Multiple Image Stacks');
  otherwise
    return
end
if isempty(datastructure)
  varargout{1} = filedir;
  close_Callback;
else
  % set instruction
  set(gui.fig, 'Name',figname)
  set(gui.handles.instructiontext,'String',instructstr);
  set(gui.handles.loadpushbutton,'String',buttonstr, ...
    'BackgroundColor',DATA.GUISettings.HighlightedButtonColor, ...
    'ForegroundColor',DATA.GUISettings.HighlightedButtonText);

  gui.datastructure = datastructure;
  gui.datatype = datatype;
  gui.filedir = filedir;

  drawpreview;
  uiwait(gui.fig)
  if ~isempty(gui.filedir)
    outdata = gui.datastructure;
    ind = [outdata.selected];
    varargout{1} = gui.filedir(ind);    
  else
    varargout{1} = gui.filedir;
  end
  close_Callback;
end



%----------------------
function drawpreview
%----------------------
% initialize all handles and draw th original preview
gui = getgui;
datastructure = gui.datastructure;

numimages = length(datastructure);
imageaxesratio = 2;
minnumimages = 8;
minNbrOfGridPoints = ceil(sqrt(ceil(minnumimages/imageaxesratio)));
nbrOfGridPoints = [1 imageaxesratio]*max(ceil(sqrt(ceil(numimages/imageaxesratio))),minNbrOfGridPoints);

thumbnailsize = length(datastructure(1).preview);
x = (0:nbrOfGridPoints(1)-1)*thumbnailsize+1;
y = (0:nbrOfGridPoints(2)-1)*thumbnailsize+1;
[gridPointX, gridPointY] = meshgrid(x,y);
gui.handles.thumbnailIm = zeros(thumbnailsize*nbrOfGridPoints);
gridPointIm = zeros(thumbnailsize*nbrOfGridPoints);
for index = 1:numimages
  x = gridPointX(index):gridPointX(index)+thumbnailsize-1;
  y = gridPointY(index):gridPointY(index)+thumbnailsize-1;
  gridPointIm(x,y) = index;

  %--- ensure that the preview image is the right size

  %If larger in some dimension, downsample
  if (size(datastructure(index).preview,1)>thumbnailsize) || (size(datastructure(index).preview,2)>thumbnailsize)
    factor = max(size(datastructure(index).preview)/thumbnailsize);
    datastructure(index).preview = imresize(datastructure(index).preview,floor(size(datastructure(index).preview)/factor));
  end

  %Zeropad if required
  temppreview = zeros(thumbnailsize,thumbnailsize,3);
  temppreview(1:size(datastructure(index).preview,1),1:size(datastructure(index).preview,2),:) = datastructure(index).preview;

  datastructure(index).preview = temppreview;

  %store the preview image
  try
    gui.handles.thumbnailIm(x,y,1:3) = rescale(datastructure(index).preview);
  catch
    disp('rfsd');
  end

end
gui.gridPointIm = gridPointIm;

gui.handles.image = imagesc(gui.handles.thumbnailIm,'parent',gui.handles.imageaxes);
set(gui.handles.imageaxes,'xtick',[],'ytick',[]);
colormap(gui.handles.imageaxes,gray);
set(gui.handles.image,'buttondownfcn','dataselector(''buttondown'')');
set(gui.fig,'windowbuttonmotionfcn','dataselector(''buttonmotion'')');

cols = imageaxesratio*nbrOfGridPoints(1);

%Draw load selection box around previews
gui.handles.loadselectbox = zeros(1,numimages);
hold(gui.handles.imageaxes,'on');
for ind = 1:numimages
  x1 = 4+mod((ind-1),cols)*thumbnailsize;
  x2 = x1+thumbnailsize-5;
  y1 = 4+floor((ind-1)/cols)*thumbnailsize;
  y2 = y1+thumbnailsize-5;
  gui.handles.loadselectbox(ind) = plot(gui.handles.imageaxes,...
    [x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'y-');
  if datastructure(ind).selected
    visiblestate = 'on';
  else
    visiblestate = 'off';
  end
  set(gui.handles.loadselectbox(ind),'visible',visiblestate,'linewidth',2);
end

% %Draw group selection box around previews
% gui.handles.groupselectcolor = [zeros(1,numimages);linspace(1,0,numimages);linspace(0,1,numimages)]';
% hsvmap = hsv(256);
% hsvmap = hsvmap(60:230,:); %limits to avoid yellows and reds
% gui.handles.groupselectcolor = hsvmap(round(linspace(1,size(hsvmap,1),numimages)),:);
% gui.handles.groupselectbox = zeros(1,numimages);
% gui.handles.groupselectbox = zeros(1,numimages);
% hold(gui.handles.imageaxes,'on');
for ind=1:numimages
  x1 = 6+mod((ind-1),cols)*thumbnailsize;
%   x2 = x1+thumbnailsize-9;
  y1 = 6+floor((ind-1)/cols)*thumbnailsize;
%   y2 = y1+thumbnailsize-9;
%   gui.handles.groupselectbox(ind) = plot(gui.handles.imageaxes,...
%     [x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'b--');
%   set(gui.handles.groupselectbox(ind),'visible','off','linewidth',2);
  gui.handles.groupnumber(ind) = text(x1,y1,sprintf('#%d',ind), ...
      'Color','w','VerticalAlignment','top');
end

%--------------------
function buttonmotion
%--------------------



%------------------
function buttondown
%-------------------

gui = getgui;

type = get(gui.fig,'SelectionType');

%get clicked coordinates
[y,x] = mygetcurrentpoint(gui.handles.imageaxes); %Obs note switched conventions.

x = round(x);
y = round(y);
imgnumber = gui.gridPointIm(x,y);
if imgnumber > 0
  % make adjustements only when a proper image was selected
  switch type
    case 'normal'
      select(imgnumber)
    case 'alt'
      [figx,figy] = mygetcurrentpoint(gui.fig);
      str = sprintf('%s #%d',dprintf('Delete'),imgnumber);
      gui.imgnumber = imgnumber;
      set(gui.handles.deleteuimenu,'Text',str);
      set(gui.fig.ContextMenu,'Visible','on','Position' ,[figx,figy]);
  end
  update;
end

%--------------
function update
%--------------
% update frames around images based on their visibility

gui = getgui;
datastructure = gui.datastructure;

%update so that selectionboxes are visible
for n = 1:length(datastructure)
  if datastructure(n).selected
    visiblestate ='on';
  else
    visiblestate = 'off';
  end
  set(gui.handles.loadselectbox(n),'visible',visiblestate);
  set(gui.handles.groupnumber(n),'String',sprintf('#%d',n));
  %   if datastructure(n).groupselected
  %     set(gui.handles.groupselectbox(n),'visible','on',...
  %       'color',[1 0 0]);
  %   else
  %     set(gui.handles.groupselectbox(n),'visible','off');
  %   end
  %   if datastructure(n).groupedto ~= 0
  %     groups = unique([datastructure.groupedto]);
  %     groups = groups(groups > 0);
  %     %find current group as fraction of all used groups, to index color
  %     groupfrac = (find(groups == datastructure(n).groupedto,1)-1) / ...
  %       (numel(groups)-1) * (numel(datastructure)-1) + 1;
  %     if isnan(groupfrac)
  %       groupfrac = 1;
  %     end
  %     set(gui.handles.groupselectbox(n),'visible','on',...
  %       'color',gui.handles.groupselectcolor(round(groupfrac),:));
  %     allmembers = find([datastructure.groupedto] == datastructure(n).groupedto);
  %     memberindex = find(allmembers == n);
  %     set(gui.handles.groupnumber(n),'String',num2str(memberindex));
  %   else
  %     set(gui.handles.groupnumber(n),'String','');
  %   end
end

%--------------------
function select(nums)
%--------------------
gui = getgui;
datastructure = gui.datastructure;

for n = nums
  endstate = ~datastructure(n).selected;
  datastructure(n).selected = endstate;
  %   for j = 1:numel(datastructure)
  %     if j == num
  % %         || datastructure(j).groupedto ~= 0 && ...
  % %         datastructure(j).groupedto == datastructure(num).groupedto
  %       datastructure(j).selected = endstate;
  %     end
  %   end
end
gui.datastructure = datastructure;

%----------------
function buttonup
%----------------
%Buttonup function when dragging to select stacks.
gui = getgui;

set(gui.fig,...
  'WindowButtonMotionFcn','seriesselector(''buttonmotion'')')

% switch gui.dragtype
%   case 'normal'
%     field = 'selected';
%   case 'extend'
%     field = 'groupselected';
%   otherwise
%     return
% end

%----------------------
function ok_Callback
%----------------------
gui = getgui;
uiresume(gui.fig)

%----------------------
function selectall
%----------------------
selectstate = true;
selecthelper(selectstate);
update;

%----------------------
function unselectall
%----------------------

selectstate = false;
selecthelper(selectstate);
update;

function selecthelper(selectstate)
gui = getgui;
datastructure = gui.datastructure;

for n = 1:numel(datastructure)
  datastructure(n).selected = selectstate;
end
gui.datastructure = datastructure;


%----------------------
function close_Callback
%----------------------
global DATA

try
  close(DATA.GUI.SeriesSelector);
catch
  close(gcbf);
end
DATA.GUI.SeriesSelector = [];
%----------------------
function gui = getgui
%----------------------
global DATA

gui = DATA.GUI.SeriesSelector;

%----------------------------------------------------
function [datastructure, filestruct] = getscreenshots
%----------------------------------------------------
%Get list of screenshots

gui = getgui;

[pathname, filestruct] = getfilenames;
gui.pathname = pathname;

if isempty(filestruct)
  datastructure = [];
  msgstr = sprintf('%s\n%s',dprintf('No valid files found in directory.'),pathname);
  myfailed(msgstr)
  return
end
datastructure = getscreenshotpreviews(filestruct);

%--------------------------------------------
function [pathname,filestruct] = getfilenames
%--------------------------------------------
% get screenshot file names

global DATA SET

name = SET(1).Report.Name;
pathname = fullfile(DATA.Pref.Pacs.ReportsheetPath, ...
  name,'Screenshots');
filestruct = dir(fullfile(pathname,'*.png'));
if ~isempty(filestruct)
  % sort files by date oldest first
  [~,ind] = sort([filestruct.datenum]);
  filestruct = filestruct(ind);
end

%----------------------
function datastructure = getscreenshotpreviews(filestruct)
%----------------------
% get preview and further information for screenshots from input file
% structure
px = 256; % size of
numfiles = numel(filestruct);
tmpstruct.filename = '';
tmpstruct.preview = [];
tmpstruct.previewsize = [];
tmpstruct.selected = true;
datastructure = repmat(tmpstruct,numfiles,1);
for n = 1:numel(filestruct)
  fname = fullfile(filestruct(n).folder,filestruct(n).name);
  im = imread(fname);
  szimage = size(im);
  szheight = szimage(1);
  szwidth = szimage(2);
  szdiff = abs((szwidth-szheight))/2;
  if szwidth >= szheight
    padpre = [floor(szdiff) 0];
    padpost = [ceil(szdiff) 0];
  else
    padpre = [0 floor(szdiff)];
    padpost = [0 ceil(szdiff)];
  end
  im = padarray(im,padpre,'pre');
  im = padarray(im,padpost,'post');

  im = imresize(im, [px px]);
  datastructure(n).filename = fname;
  datastructure(n).preview = im;
  datastructure(n).previewsize = size(im);
  datastructure(n).selected = true;
end

%--------------------------
function deleteall_Callback
%--------------------------
%Delete all screenshots

gui = getgui;
filedir = gui.pathname;
if yesno('Are you sure?')
  delete([filedir filesep '*.png']);
  logdisp('All screenshots deleted');
  ok_Callback();
end

%----------------------
function deletefile(~,~)
%----------------------
% delete specified file from disc

gui = getgui;
imgnumber = gui.imgnumber;
filedir = gui.filedir;
fname = fullfile(filedir(imgnumber).folder,filedir(imgnumber).name);
askstr = sprintf('%s #%d\n%s', dprintf('Delete'),imgnumber,fname);
if yesno(askstr)
  delete(fname);
  % update
  numimages = numel(gui.datastructure);
  ind = true(1,numimages);
  ind(imgnumber) = false;
  gui.datastructure = gui.datastructure(ind);
  gui.filedir = filedir(ind);
  if ~isempty (gui.filedir)
    drawpreview;
    gui.imgnumber = [];
  else    
    uiresume(gui.fig)
  end
end
