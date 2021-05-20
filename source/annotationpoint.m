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

%--------------------------------
function pointclearall_helper(no)
%--------------------------------
%Helperfunction to clear all points

global SET

SET(no).Point = [];
SET(no).Point.X = [];
SET(no).Point.Y = [];
SET(no).Point.T = [];
SET(no).Point.Z = [];
SET(no).Point.Label = {};

%----------------------
function pointclearall
%----------------------
%Clear all points.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end
tools('enableundo',no);

%Update in 2D
for p = 1:length(DATA.ViewPanels)
  viewfunctions('updatedrawlist',p)
  drawfunctions('drawpanel',p)
end

try
  if segment3dp.isviewportalive
    DATA.LevelSet.ViewPort.setpoints([],[],[])
  end
catch
end

pointclearall_helper(no);

%-------------------------------
function pointexportall_Callback %#ok<DEFNU>
%-------------------------------
%Export all point data
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

if isempty(SET(no).Point)
  pointclearall;
  myfailed('No points marked',DATA.GUI.Segment);
  return;
end

if length(SET(no).Point.X)<1
  myfailed('No points marked',DATA.GUI.Segment);
  return;
end

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
  end
else
  offset = 0;
end

%--- Time resolved points
if sum(not(isnan(SET(no).Point.T)))>0
  %Basic titles 
  outdata{3+offset,1} = 'Time resolved points';
  outdata{4+offset,1} = 'Time[ms]';
  for tloop=1:SET(no).TSize
    outdata{4+offset+tloop,1} = (tloop-1)*SET(no).TIncr*1000;
  end
    
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
  end
  
  for lloop=1:length(labelstodo)
    for loop=1:length(SET(no).Point.X)
      if not(isnan(SET(no).Point.T(loop)))&&isequal(labelstodo{lloop},SET(no).Point.Label{loop})
        if not(isempty(outdata{4+offset+SET(no).Point.T(loop),2+(lloop-1)*4}))
          mywarning(dprintf('Multiple points with same label, ambigous result. Timeframe %d, slice %d ',...
            SET(no).Point.T(loop),SET(no).Point.Z(loop)),DATA.GUI.Segment);
        end
        outdata{4+offset+SET(no).Point.T(loop),2+(lloop-1)*4} = SET(no).Point.Label{loop};
        outdata{4+offset+SET(no).Point.T(loop),3+(lloop-1)*4} = SET(no).Point.X(loop)*SET(no).ResolutionX;
        outdata{4+offset+SET(no).Point.T(loop),4+(lloop-1)*4} = SET(no).Point.Y(loop)*SET(no).ResolutionY;
        outdata{4+offset+SET(no).Point.T(loop),5+(lloop-1)*4} = SET(no).Point.Z(loop)*(SET(no).SliceThickness+SET(no).SliceGap);
      end
    end
  end
end

segment('cell2clipboard',outdata);

%------------------------------------
function pointcleartemplate_Callback %#ok<DEFNU>
%------------------------------------
%Clear points using naming template.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

if isempty(SET(no).Point)||isempty(SET(no).Point.X)
  myfailed('No points to delete.',DATA.GUI.Segment)
  return;
end

s = myinputdlg({'Delete points labeled as:'},'Template',1,{''});
if isempty(s)
  myfailed('Invalid template.',DATA.GUI.Segment);
  return;
else
  s = s{1};
end

tools('enableundo',no);

ind = true(1,length(SET(no).Point.Z));
for loop=1:length(SET(no).Point.Z)
  if isequal(SET(no).Point.Label{loop},s)
    ind(loop) = false;
  end
end

if ~yesno(dprintf('Deleting %d points are you sure?',sum(not(ind))),[],DATA.GUI.Segment)
  return;
end

SET(no).Point.X = SET(no).Point.X(ind);
SET(no).Point.Y = SET(no).Point.Y(ind);
SET(no).Point.T = SET(no).Point.T(ind);
SET(no).Point.Z = SET(no).Point.Z(ind);
SET(no).Point.Label = SET(no).Point.Label(ind);
viewfunctions('setview');  %drawfunctions('drawimageno');

%------------------------------------
function pointrenametemplate_Callback %#ok<DEFNU>
%------------------------------------
%Rename points according to a renaming template.
global DATA SET NO

%Use to point to mag data set
no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

if isempty(SET(no).Point)||isempty(SET(no).Point.X)
  pointclearall;
  myfailed('No points exist.',DATA.GUI.Segment);
  return;
end

s = myinputdlg({'Rename points labeled as:'},'Template',1,{''});
if isempty(s)
  myfailed('Invalid template.',DATA.GUI.Segment);
  return;
else
  s = s{1};
end

stri = myinputdlg({'New name:'},'Newname',1,{''});
if isempty(stri)
  myfailed('Invalid new name.',DATA.GUI.Segment);
  return;
else
  stri = stri{1};
end

tools('enableundo',no);

for loop=1:length(SET(no).Point.Z)
  if isequal(SET(no).Point.Label{loop},s)
    SET(no).Point.Label{loop} = stri;
  end
end

viewfunctions('setview');  %drawfunctions('drawimageno');

%----------------------------
function filterpoints_Callback %#ok<DEFNU>
%----------------------------
%Filter point in time using a Kalman filter

%Einar Heiberg

global SET NO

no = NO;
 
%Find points to filter
labels = SET(no).Point.Label;
labels = union(labels,{}); %Remove duplicates.

c = mymenu('Select point to filter',labels{:});
if isequal(c,0)
%   myfailed('Aborted.');
  return;
end

name = labels{c}; %This is the name to be found.

%---Extract coordinates.

%Find points
ind = false(1,length(SET(no).Point.X));
for loop = 1:length(SET(no).Point.X)
  if isequal(SET(no).Point.Label{loop},name)
    ind(loop) = true;
  end
end

%Extract
x = SET(no).Point.X;
y = SET(no).Point.Y;
x = x(ind);
y = y(ind);

%Check
if ~isequal(length(x),SET(no).TSize)
  myfailed('Number of points is not equal to number of timeframes. Duplicated point names?');
  return;  
end

%Ask for standard deviation
s.Noise = 0.3;
[s,ok] = inputstruct(s,'Enter Noise estimate');
if ~ok
  myfailed('Aborted or illegal value.');
  return;
end

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
title(dprintf('X-coordinate.'));
xlabel(dprintf('Time [s]'));
ylabel(dprintf('Position (pixel)'));

%X coordinate
subplot(2,1,2);
plot(t,[y y y],'b.');
hold on;
plot(t,ynew,'b-');
hold off;
title(dprintf('Y-coordinate.'));
xlabel(dprintf('Time [s]'));
ylabel(dprintf('Position (pixel)'));

set(99,'numbertitle','off','name',dprintf('Position over time'));

%Crop it back from triplicate
xnew = xnew((length(x)+1):2*length(x));
ynew = ynew((length(y)+1):2*length(y));

if yesno('Do you want to apply filter?')
  %Prepare to store by backup...
  tools('enableundo',no);

  %Store the result
  SET(no).Point.X(ind) = xnew;
  SET(no).Point.Y(ind) = ynew;
  
  mymsgbox('Filter applied.');  
end

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
end

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
  end
end

if no==NO
  myfailed('Cannot import from same image stack.',DATA.GUI.Segment);
  return;
end
if (no>length(SET))||(no<1)
  myfailed('Invalid image stack selected.',DATA.GUI.Segment);
  return;
end

%Check what to do - later optionally take from input arguments.
if ~isempty(SET(no).Point)    
    importpointhelper(NO,no);
else
    return;
end

segment('updatevolume');
viewfunctions('setview');

%----------------------------------------------------------------------------------
function importpointhelper(tono,fromno)
%----------------------------------------------------------------------------------
%Helper function to segmentimportsegmention Callback. 
%- tono is destination of segmentation.
%- fromno is source.

global SET

x = SET(fromno).Point.X;
y = SET(fromno).Point.Y;
z = SET(fromno).Point.Z;
  
%Convert coordinates
pos = calcfunctions('xyz2rlapfh',fromno,x,y,z);
pos = calcfunctions('rlapfh2xyz',tono,pos(:,1),pos(:,2),pos(:,3));

%Copy the points
for loop = 1:size(pos,2)  
  x = pos(1,loop);
  y = pos(2,loop);
  z = pos(3,loop);
  if ...
      (x>0) && (x<=SET(tono).XSize) && ...
      (y>0) && (y<=SET(tono).YSize) && ...
      (z>0) && (z<=SET(tono).ZSize)
    SET(tono).Point.X(end+1) = x;
    SET(tono).Point.Y(end+1) = y;
    SET(tono).Point.Z(end+1) = z;
    SET(tono).Point.T(end+1) = SET(fromno).Point.T(loop);
    SET(tono).Point.Label{end+1} = SET(fromno).Point.Label{loop};    
  end
 
end

