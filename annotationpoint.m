function varargout = annotationpoint(varargin)
% ANNOTATIONPOINT
% Functions to place and modify annotation points.

% Moved out from segment_main.m by Nisse Lundahl

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%---------------------------------
function pointmove_Callback(dx,dy) %#ok<DEFNU>
%---------------------------------
%Helper function to move points.
global SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

pos = pointfind(true); %Call it in silent mode

if isnan(pos)
  return;
end;

tools('enableundo',no);

SET(no).Point.X(pos) = SET(no).Point.X(pos)+dx;
SET(no).Point.Y(pos) = SET(no).Point.Y(pos)+dy;

drawfunctions('drawsliceno');

%-----------------------------
function pointforward_Callback %#ok<DEFNU>
%-----------------------------
%Track point forward in time.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

pos = pointfind;

if isnan(SET(no).Point.T(pos))
  myfailed('Not a time-resolved point.',DATA.GUI.Segment);
  return;
end;

x = SET(no).Point.X(pos);
y = SET(no).Point.Y(pos);
z = SET(no).Point.Z(pos);
label = SET(no).Point.Label{pos};

%Find corresponding point
t = SET(no).CurrentTimeFrame;
t = t+1;
if t>SET(no).TSize
  t = 1;
end;

%Possible matches
ind = find((round(SET(no).Point.Z)==z)&(SET(no).Point.T==t));
score = 1e10; %Smaller the better
foundind = NaN;
for loop=1:length(ind)
  if isequal(SET(no).Point.Label{ind(loop)},label)
    temp = sqrt(...
      (SET(no).Point.X(ind(loop))-x).^2+...
      (SET(no).Point.Y(ind(loop))-y).^2);
    if temp<score
      score = temp;
      foundind = ind(loop);
    end;
  end;
end;

tools('enableundo',no);

if isnan(foundind)
  %Found no corresponding point, create
  newind = length(SET(no).Point.X)+1;
  SET(no).Point.X(newind) = x;
  SET(no).Point.Y(newind) = y;  
  SET(no).Point.Z(newind) = z;
  SET(no).Point.T(newind) = t;
  SET(no).Point.Label{newind} = label;  
else
  %Found point, modify.
  SET(no).Point.X(foundind) = x;
  SET(no).Point.Y(foundind) = y;
end;
  
SET(no).CurrentTimeFrame = t;

drawfunctions('drawimageno');


%----------------------
function pointclearall
%----------------------
%Clear all points.
global SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;
tools('enableundo',no);

SET(no).Point = [];
SET(no).Point.X = [];
SET(no).Point.Y = [];
SET(no).Point.T = [];
SET(no).Point.Z = [];
SET(no).Point.Label = {};

%-------------------------------
function pointexportall_Callback %#ok<DEFNU>
%-------------------------------
%Export all point data
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

if isempty(SET(no).Point)
  pointclearall;
  myfailed('No points marked',DATA.GUI.Segment);
  return;
end;

if length(SET(no).Point.X)<1
  myfailed('No points marked',DATA.GUI.Segment);
  return;
end;

outdata = cell(10,10);
outdata{1,1} = 'Filename';
outdata{1,2} = SET(no).FileName;

%Non time resolved points
if sum(isnan(SET(no).Point.T))>0
  outdata{3,1} = 'Non time resolved points';
  outdata{4,1} = 'Label';
  outdata{4,2} = 'X';
  outdata{4,3} = 'Y';
  outdata{4,4} = 'Z';
  pos = find(isnan(SET(no).Point.T));
  offset = 2+length(pos)+1;
  for loop=1:length(pos)
    outdata{4+loop,1} = SET(no).Point.Label{pos(loop)};
    outdata{4+loop,2} = SET(no).Point.X(pos(loop))*SET(no).ResolutionX;
    outdata{4+loop,3} = SET(no).Point.Y(pos(loop))*SET(no).ResolutionY;
    outdata{4+loop,4} = SET(no).Point.Z(pos(loop))*(SET(no).SliceThickness+SET(no).SliceGap);    
  end;
else
  offset = 0;
end;

%--- Time resolved points
if sum(not(isnan(SET(no).Point.T)))>0
  %Basic titles 
  outdata{3+offset,1} = 'Time resolved points';
  outdata{4+offset,1} = 'Time[ms]';
  for tloop=1:SET(no).TSize
    outdata{4+offset+tloop,1} = (tloop-1)*SET(no).TIncr*1000;
  end;
    
  %extract data
  ind = not(isnan(SET(no).Point.T));
  labels = SET(no).Point.Label(ind);
  
  %Find unique names
  labelstodo = unique(labels);
  
  for loop=1:length(labelstodo)
    outdata{4+offset,2+(loop-1)*4} = 'Label';
    outdata{4+offset,3+(loop-1)*4} = 'X';
    outdata{4+offset,4+(loop-1)*4} = 'Y';
    outdata{4+offset,5+(loop-1)*4} = 'Z';    
  end;
  
  for lloop=1:length(labelstodo)
    for loop=1:length(SET(no).Point.X)
      if not(isnan(SET(no).Point.T(loop)))&&isequal(labelstodo{lloop},SET(no).Point.Label{loop})
        if not(isempty(outdata{4+offset+SET(no).Point.T(loop),2+(lloop-1)*4}))
          mywarning(dprintf(...
            'Multiple points with same label, ambigous result. Timeframe %d, slice %d ',...
            SET(no).Point.T(loop),SET(no).Point.Z(loop)),DATA.GUI.Segment);
        end;
        outdata{4+offset+SET(no).Point.T(loop),2+(lloop-1)*4} = SET(no).Point.Label{loop};
        outdata{4+offset+SET(no).Point.T(loop),3+(lloop-1)*4} = SET(no).Point.X(loop)*SET(no).ResolutionX;
        outdata{4+offset+SET(no).Point.T(loop),4+(lloop-1)*4} = SET(no).Point.Y(loop)*SET(no).ResolutionY;
        outdata{4+offset+SET(no).Point.T(loop),5+(lloop-1)*4} = SET(no).Point.Z(loop)*(SET(no).SliceThickness+SET(no).SliceGap);
      end;
    end;
  end;
end;

segment('cell2clipboard',outdata);

%-------------------------------
function pointdeletethis_Callback %#ok<DEFNU>
%-------------------------------
%Delete point

global SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;
tools('enableundo',no);

pos = pointfind;
if isnan(pos)
  myfailed('No point found.');
  return;
end
ind = true(1,length(SET(no).Point.X));
ind(pos) = false;

SET(no).Point.X = SET(no).Point.X(ind);
SET(no).Point.Y = SET(no).Point.Y(ind);
SET(no).Point.T = SET(no).Point.T(ind);
SET(no).Point.Z = SET(no).Point.Z(ind);
SET(no).Point.Label = SET(no).Point.Label(ind);
drawfunctions('drawimageno');

%------------------------------------
function pointcleartemplate_Callback %#ok<DEFNU>
%------------------------------------
%Clear points using naming template.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

if isempty(SET(no).Point)||isempty(SET(no).Point.X)
  myfailed('No points to delete.',DATA.GUI.Segment)
  return;
end;

s = inputdlg({'Delete points labeled as:'},'Template',1,{''});
if isempty(s)
  myfailed('Invalid template.',DATA.GUI.Segment);
  return;
else
  s = s{1};
end;

tools('enableundo',no);

ind = true(1,length(SET(no).Point.Z));
for loop=1:length(SET(no).Point.Z)
  if isequal(SET(no).Point.Label{loop},s)
    ind(loop) = false;
  end;
end;

if ~yesno(dprintf('Deleting %d points are you sure?',sum(not(ind))),[],DATA.GUI.Segment);
  return;
end;

SET(no).Point.X = SET(no).Point.X(ind);
SET(no).Point.Y = SET(no).Point.Y(ind);
SET(no).Point.T = SET(no).Point.T(ind);
SET(no).Point.Z = SET(no).Point.Z(ind);
SET(no).Point.Label = SET(no).Point.Label(ind);
drawfunctions('drawimageno');

%------------------------------------------
function pointmaketimeresolvedthis_Callback %#ok<DEFNU>
%------------------------------------------
%Converts a none time resolved point to a time resolved point.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;
tools('enableundo',no);

pos = pointfind;

if isnan(pos)
  myfailed('No point found');
 return;
end

ind = true(1,length(SET(no).Point.X));
ind(pos) = false;

%Backup 
x = SET(no).Point.X(pos);
y = SET(no).Point.Y(pos);
z = SET(no).Point.Z(pos);
label = SET(no).Point.Label(pos);

%Remove old
SET(no).Point.X = SET(no).Point.X(ind);
SET(no).Point.Y = SET(no).Point.Y(ind);
SET(no).Point.T = SET(no).Point.T(ind);
SET(no).Point.Z = SET(no).Point.Z(ind);
SET(no).Point.Label = SET(no).Point.Label(ind);

%Add new
x = repmat(x,[1 SET(no).TSize]);
y = repmat(y,[1 SET(no).TSize]);
z = repmat(z,[1 SET(no).TSize]);
t = 1:SET(no).TSize;
label = repmat(label,[1 SET(no).TSize]);
SET(no).Point.X = cat(2,SET(no).Point.X,x);
SET(no).Point.Y = cat(2,SET(no).Point.Y,y);
SET(no).Point.Z = cat(2,SET(no).Point.Z,z);
SET(no).Point.T = cat(2,SET(no).Point.T,t);
SET(no).Point.Label = cat(2,SET(no).Point.Label,label);

mymsgbox('Point now timeresolved.','Done!',DATA.GUI.Segment);

%Update image
drawfunctions('drawimageno');

%------------------------------------
function pointrenametemplate_Callback %#ok<DEFNU>
%------------------------------------
%Rename points according to a renaming template.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

if isempty(SET(no).Point)||isempty(SET(no).Point.X)
  pointclearall;
  myfailed('No points exist.',DATA.GUI.Segment);
  return;
end;

s = inputdlg({'Rename points labeled as:'},'Template',1,{''});
if isempty(s)
  myfailed('Invalid template.',DATA.GUI.Segment);
  return;
else
  s = s{1};
end;

stri = inputdlg({'New name:'},'Newname',1,{''});
if isempty(stri)
  myfailed('Invalid new name.',DATA.GUI.Segment);
  return;
else
  stri = stri{1};
end;

tools('enableundo',no);

for loop=1:length(SET(no).Point.Z)
  if isequal(SET(no).Point.Label{loop},s)
    SET(no).Point.Label{loop} = stri;
  end;
end;

drawfunctions('drawimageno');

%--------------------------------
function pointrenamethis_Callback %#ok<DEFNU>
%--------------------------------
%Rename point 
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;
tools('enableundo',no);

ind = pointfind;

if isnan(ind)
  myfailed('No point found.');
  return;
end

s = inputdlg({'Enter name'},'Name',1,{sprintf('%s',SET(no).Point.Label{ind})});
if isempty(s)
  myfailed('Invalid name.',DATA.GUI.Segment);
  return;
end;
SET(no).Point.Label(ind) = s;

drawfunctions('drawimageno');

%------------------------------
function pointclearall_Callback %#ok<DEFNU>
%------------------------------
%Clear all points.
global SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;
tools('enableundo',no);

SET(no).Point.X = [];
SET(no).Point.Y = [];
SET(no).Point.T = [];
SET(no).Point.Z = [];
SET(no).Point.Label = {};
drawfunctions('drawimageno');

%----------------------
function point_Buttonup %#ok<DEFNU>
%----------------------
%Button up function for points.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;
tools('enableundo',no);

%Restore so no motion is called
set(DATA.imagefig,'WindowButtonMotionFcn','');

%Restore main buttonup function
set(DATA.imagefig,'WindowButtonUpFcn',...
  sprintf('%s(''buttonup_Callback'')','segment'));

SET(no).Point.X(DATA.MeasureN) = DATA.MeasureX;
SET(no).Point.Y(DATA.MeasureN) = DATA.MeasureY;

drawfunctions('drawimageno');



%--------------------
function point_Motion %#ok<DEFNU>
%--------------------
%Motion function for points.
global DATA SET NO

%If straintagging initiated adjust LVupdated
if ~isempty(SET(NO).StrainTagging) && isfield(SET(NO).StrainTagging, 'LVupdated')
  SET(NO).StrainTagging.LVupdated = 1;
end

[x,y,slice] = segment('getclickedcoords');

% first line is for montage
% second line is for one-view.
if (slice==SET(NO).CurrentSlice) && ...
   (x>0.5) && (y>0.5) && (x<SET(NO).YSize+0.5) && (y<SET(NO).XSize+0.5)
  DATA.MeasureX = y;
  DATA.MeasureY = x;

  set(DATA.Handles.pointp{DATA.CurrentPanel}(DATA.MeasureN),...
    'xdata',DATA.MeasureOffsetX+DATA.MeasureY,...
    'ydata',DATA.MeasureOffsetY+DATA.MeasureX);
  set(DATA.Handles.pointo{DATA.CurrentPanel}(DATA.MeasureN),...
    'xdata',DATA.MeasureOffsetX+DATA.MeasureY,...
    'ydata',DATA.MeasureOffsetY+DATA.MeasureX);
  set(DATA.Handles.pointtext{DATA.CurrentPanel}(DATA.MeasureN),'position',[...
    DATA.MeasureOffsetX+DATA.MeasureY+2 ...
    DATA.MeasureOffsetY+DATA.MeasureX ...
    0]);
end
%-------------------------------
function ind = pointfind(silent) %#ok<INUSD>
%-------------------------------
%Find neareast point.
global DATA SET NO

%Get clicked coordinate
[x,y,slice] = segment('getclickedcoords');

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

%Find correct point
ind = NaN;
mindist = 1e10;
for loop=1:length(SET(no).Point.X)
  if slice==round(SET(no).Point.Z(loop))
    if (isnan(SET(no).Point.T(loop))) || (SET(no).Point.T(loop)==SET(no).CurrentTimeFrame)
      dist = min(sqrt(...
        (SET(no).Point.X(loop)-y).^2+...
        (SET(no).Point.Y(loop)-x).^2));
      if dist<mindist
        ind = loop;
        mindist = dist;
      end;
    end;
  end;
end;

if nargin==0
  if isnan(ind)
    disp('Could not find measurement point, should not occur.');
    return;
  end;
end;

%--------------------------------------
function pointat_Buttondown(panel,type) %#ok<DEFNU>
%--------------------------------------
%Buttondown function when clicked at a point.
global DATA SET NO

if nargin==1
  type = get(DATA.imagefig,'SelectionType');
end;

segment('switchtopanel',panel);

[~,~,slice] = segment('getclickedcoords'); 
if (slice>SET(NO).ZSize)
  return;
end
segment('switchtoslice',slice);

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end;

ind = pointfind;

%Prepare DATA.Measure info
DATA.MeasureN = ind;
DATA.MeasureX = SET(no).Point.X(ind);
DATA.MeasureY = SET(no).Point.Y(ind);
switch DATA.ViewPanelsType{panel}
  case {'one','mmodespatial','ortho'}
    DATA.MeasureOffsetX = 0;
    DATA.MeasureOffsetY = 0;
  case {'montage','montagerow','montagefit','sax3'}
    DATA.MeasureOffsetY = SET(no).XSize*floor((slice-1)/DATA.ViewPanelsMatrix{panel}(2));
    DATA.MeasureOffsetX = SET(no).YSize*mod(slice-1,DATA.ViewPanelsMatrix{panel}(2));
  otherwise
    return;
end;

switch type
  case 'alt'
    DATA.contextmenu;
  otherwise
    %Call point motion
    set(DATA.imagefig,'WindowButtonMotionFcn',sprintf('%s(''point_Motion'');',mfilename));
    set(DATA.imagefig,'WindowButtonUpFcn',...
      sprintf('%s(''point_Buttonup'')',mfilename));
end;


%----------------------------
function filterpoints_Callback %#ok<DEFNU>
%----------------------------
%Filter point in time using a Kalman filter

%Einar Heiberg

global SET NO

no = NO;
 
%Find points to filter
z = SET(no).Point.Z;
labels = SET(no).Point.Label;
labels = union(labels,{}); %Remove duplicates.

c = mymenu('Select point to filter',labels{:});
if isequal(c,0)
  myfailed('Aborted.');
  return;
end;

name = labels{c}; %This is the name to be found.

%---Extract coordinates.

%Find points
ind = false(1,length(SET(no).Point.X));
for loop = 1:length(SET(no).Point.X)
  if isequal(SET(no).Point.Label{loop},name)
    ind(loop) = true;
  end;
end;

%Extract
x = SET(no).Point.X;
y = SET(no).Point.Y;
x = x(ind);
y = y(ind);

%Check
if ~isequal(length(x),SET(no).TSize)
  myfailed('Number of points is not equal to number of timeframes. Duplicated point names?');
  return;  
end;

%Ask for standard deviation
s.Noise = 0.3;
[s,ok] = inputstruct(s,'Enter Noise estimate');
if ~ok
  myfailed('Aborted or illegal value.');
  return;
end;

%Extract it
objectnoise = s.Noise;

%Apply filter
xnew = kalman([x x x],objectnoise);
ynew = kalman([y y y],objectnoise);

%Plot to ask user
figure(99);
t = [SET(no).TimeVector SET(no).TimeVector+SET(no).TimeVector(end) SET(no).TimeVector+2*SET(no).TimeVector(end)];
clf;

%X coordinate
subplot(2,1,1);
plot(t,[x x x],'b.');
hold on;
plot(t,xnew,'b-');
hold off;
title('X-coordinate.');
xlabel('Time [s]');
ylabel('Position (pixel)');

%X coordinate
subplot(2,1,2);
plot(t,[y y y],'b.');
hold on;
plot(t,ynew,'b-');
hold off;
title('Y-coordinate.');
xlabel('Time [s]');
ylabel('Position (pixel)');

set(99,'numbertitle','off','name','Position over time');

%Crop it back from triplicate
xnew = xnew((length(x)+1):2*length(x));
ynew = ynew((length(y)+1):2*length(y));

if yesno('Do you want to apply filter?');
  %Prepare to store by backup...
  tools('enableundo',no);

  %Store the result
  SET(no).Point.X(ind) = xnew;
  SET(no).Point.Y(ind) = ynew;
  
  mymsgbox('Filter applied.');  
end;

close(99);


%----------------------------------------------
function importpoint_Callback(no) %#ok<DEFNU>
%----------------------------------------------
%Import point from another image stack.
%Imports to current image stack NO from no or if called with no input
%arguments user is asked.

global DATA SET NO

tools('enableundo');

if length(SET)<2
  myfailed('Only one image stack in memory, import from file instead (under File menu).',DATA.GUI.Segment);
  return;
end;

if nargin==0
  %Find what imagestack
  menuitems = cell(1,numel(SET)-1);
  impstacks = setdiff(1:numel(SET),NO);
  for n = 1:numel(impstacks)
    nn = impstacks(n);
    menuitems{n} = sprintf('%d. %s',nn,[SET(nn).ImageType ' / ' SET(nn).ImageViewPlane]);
  end
  s = mymenu('Select which stack to import from',menuitems);

  if s == 0 %operation cancelled
    return;
  else
    no = impstacks(s);
  end;
end;

if no==NO
  myfailed('Cannot import from same image stack.',DATA.GUI.Segment);
  return;
end;
if (no>length(SET))||(no<1)
  myfailed('Invalid image stack selected.',DATA.GUI.Segment);
  return;
end;

%Check what to do - later optionally take from input arguments.
if ~isempty(SET(no).Point)    
    importpointhelper(NO,no);
else
    return;
end


segment('updatemodeldisplay');
segment('updatevolume');
segment_main('viewrefreshall_Callback');
drawfunctions('drawallslices');

%----------------------------------------------------------------------------------
function importpointhelper(tono,fromno)
%----------------------------------------------------------------------------------
%Helper function to segmentimportsegmention Callback. 
%- tono is destination of segmentation.
%- fromno is source.
%
%The function is capable of handeling slice offsets and different
%pixelssizes as well as situations when number of timeframes differ. When
%destination is not timeresolved and source is timeresolved then user is
%asked from what timeframe to take the segmentation from. 
%
%Work horse in importing. This function would benefit from anti-cut 
%and paste treatment.

global SET

[sourceslices,destslices,sourcetime,desttime,zdirsource,zdirdest] = ...
  segmentation('findmatchingslices',tono,fromno,0,0,0,0);

%Loop over the number of slices in destination images
if ~isempty(SET(tono).Point)
    currentind = length(SET(tono).Point.X);
else
    currentind = 0;
end
for zloop = 1:length(destslices)
  
  %Match slices, see above.
  sourceslice = sourceslices(zloop);
  destslice = destslices(zloop);
  
  %Match positions
  cornerpossource = SET(fromno).ImagePosition-zdirsource*(zloop-1)*(SET(fromno).SliceThickness+SET(fromno).SliceGap);
  cornerposdest = SET(tono).ImagePosition-zdirdest*(destslice-1)*(SET(tono).SliceThickness+SET(tono).SliceGap);
  
  xofs = (cornerpossource-cornerposdest)*(SET(tono).ImageOrientation(4:6)'); %Project on one axis in mm
  yofs = (cornerpossource-cornerposdest)*(SET(tono).ImageOrientation(1:3)'); %Project on the other axis in mm

  xofs = xofs / SET(tono).ResolutionX; %Now in pixels in destination coordinates
  yofs = yofs / SET(tono).ResolutionY; %Now in pixels in destination coordinates  
  
  %factor between
  fx = SET(fromno).ResolutionX/SET(tono).ResolutionX;
  fy = SET(fromno).ResolutionY/SET(tono).ResolutionY;  
  
  ind = find(SET(fromno).Point.Z == sourceslice);
  if ~isempty(ind)
      %Equal time resolution
      if isequal(desttime,sourcetime)
          SET(tono).Point.X(currentind+1:currentind+length(ind)) = min(SET(tono).XSize,max(1,SET(fromno).Point.X(ind)*fx+xofs));
          SET(tono).Point.Y(currentind+1:currentind+length(ind)) = min(SET(tono).YSize,max(1,SET(fromno).Point.Y(ind)*fy+yofs));
          SET(tono).Point.T(currentind+1:currentind+length(ind)) = SET(fromno).Point.T(ind);
          SET(tono).Point.Z(currentind+1:currentind+length(ind)) = ones(1,length(ind))*destslice;
          for indloop = 1:length(ind)
              SET(tono).Point.Label{currentind+indloop} = SET(fromno).Point.Label{ind(indloop)};
          end
          currentind =  currentind+length(ind);
      end
  end
    
end;
