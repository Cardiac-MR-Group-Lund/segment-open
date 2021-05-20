function [varargout] = atrialscar(varargin)
%Module for atrial scar quantification.

%Einar Heiberg

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%Todo:
%- resample DATA.NumPoints
%- data summary (before visualization)
%- auto scar
%
%- DONE add to menu
%- DONE visualization
%- DONE add drawslicepanelone
%- DONE finish DATA.Handles with atrial scar...
%- DONE hide/show scar 

%-------------------------
function clearall_Callback %#ok<DEFNU>
%-------------------------
%Callback to clear all atrial scar data.

global SET NO DATA

SET(NO).AtrialScar = [];
drawfunctions('drawpanel',DATA.CurrentPanel);
segment_main('updatevolume');

%--------------------
function initscar(no)
%--------------------
global SET

SET(no).AtrialScar = [];
SET(no).AtrialScar.Manual = int8(zeros(size(SET(no).RVEndoX)));
SET(no).AtrialScar.Auto = false(size(SET(no).RVEndoX));
SET(no).AtrialScar.Result = false(size(SET(no).RVEndoX));
SET(no).AtrialScar.IM = nan(size(SET(no).RVEndoX));
SET(no).AtrialScar.Percentage = 0;
SET(no).AtrialScar.TotArea = 0;

%------------------
function update(no)
%------------------
%Updates atrial volume and graphical update

global DATA SET

%Update scar struct
SET(no).AtrialScar.Result = SET(no).AtrialScar.Auto; %Start with auto
SET(no).AtrialScar.Result = SET(no).AtrialScar.Result | (SET(no).AtrialScar.Manual>0); %Add manual adds (scar)
SET(no).AtrialScar.Result = SET(no).AtrialScar.Result & (SET(no).AtrialScar.Manual>=0); %Remove rubber regions

%Extract data
x = squeeze(SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,:));
y = squeeze(SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,:));

%Update percentage
slicethickness = (SET(no).SliceThickness+SET(no).SliceGap)/10; %in cm
totarea = 0;
totscararea = 0;
for slice = 1:SET(no).ZSize
  if ~isnan(x(1,slice))
    xr = x(:,slice)*SET(no).ResolutionX;
    yr = y(:,slice)*SET(no).ResolutionY;
    rr = SET(no).AtrialScar.Result(:,slice);
    
    %Calc length
    dist = sqrt((xr-circshift(xr,[1 0])).^2+(yr-circshift(yr,[1 0])).^2)/10; %10 to convert to cm
    totarea = totarea+sum(dist(:))*slicethickness;
    totscararea = totscararea+sum(dist(:))*slicethickness*(sum(rr))/length(rr);
  end
end

SET(no).AtrialScar.TotArea = totarea;
SET(no).AtrialScar.Percentage = (totscararea/totarea)*100;

%drawfunctions('drawimagepanel',DATA.CurrentPanel);
drawfunctions('drawpanel',DATA.CurrentPanel);
segment_main('updatevolume');

%-----------------------------------------
function manualdraw_Buttonup(no,type,xr,yr) %#ok<DEFNU>
%-----------------------------------------
%Callback called from buttonupfunctions.m

global SET

if isempty(SET(no).AtrialScar)
  initscar(no); %Initalize data structure
end

%Create a mask
mask = segment('createmask',[SET(no).XSize SET(no).YSize],xr,yr);

%Interpolate RV contour
x = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
y = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
z = interp2(double(mask),y,x);

ind = find(z>0.5); %logical index

%Assign
switch lower(type)
  case 'scar'
    SET(no).AtrialScar.Manual(ind,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = int8(1);    
  case 'rubberpen'
    SET(no).AtrialScar.Manual(ind,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = int8(-1);
end

update(no);

%---------------------------
function visualize_Callback %#ok<DEFNU>
%---------------------------
%Plots a 3D model of the atrial with colored surface with scar

global SET NO

if isempty(SET(NO).AtrialScar)
  myfailed('No Atrial Scar data found.');
  return;
end

%Extract data
tf = SET(NO).CurrentTimeFrame;
x = SET(NO).RVEndoX;
y = SET(NO).RVEndoY;

if isempty(x)
  myfailed('Need delineation using RV Endo tool');
  return
end

%Extract only valid slices
slices = find(~isnan(squeeze(x(1,tf,:))));

%Check if slices are contigous
contcheck = sum(diff(slices)~=1);
if contcheck>0
  myfailed('Slices are not contigous. Missing segmentation in some slice.');
  return;
end

if length(slices)<2 
  myfailed('Expected at least two slices for atrial scar visualization.');
  return;
end

%Extract tf,slices
x = squeeze(x(:,tf,slices));
y = squeeze(y(:,tf,slices));
result = squeeze(SET(NO).AtrialScar.Result(:,tf,slices));

%Create mask for surface
mask = zeros(SET(NO).XSize,SET(NO).YSize,length(slices));
for loop = 1:length(slices)
  tempmask = segment('createmask',[SET(NO).XSize SET(NO).YSize],x(:,loop),y(:,loop));
  mask(:,:,loop) = tempmask;
end

%--- Create smoothing filter
n = 7;
xi = linspace(-2,2,n);
f = exp(-xi.^2);
f = f./sum(f(:));

%figure(83);
%plot(xi,f);

%--- Smooth mask
mask = econv3(mask,f);
mask = econv3(mask,f');
mask = econv3(mask,reshape(f,[1 1 n]));

%--- Calculate isosurface
fv = isosurface(mask,0.5);

%--- Create scar mask
mask = zeros(SET(NO).XSize,SET(NO).YSize,length(slices));
for loop = 1:length(slices)
  
  %Extract contour
  xr = x(:,loop);
  yr = y(:,loop);
  mx = mean(xr);
  my = mean(yr);
  %make the line thicker to compensate for smoothed myocardium
  for k = linspace(0.8,1.2,20)
    xr2 = xr-mx;
    yr2 = yr-my;
    xr2 = mx+xr2*k;
    yr2 = my+yr2*k;
    
    %Take only infarct points
    xr2 = xr2(result(:,loop)>0.5);
    yr2 = yr2(result(:,loop)>0.5);
    
    ind = sub2ind([SET(NO).XSize SET(NO).YSize length(slices)],round(xr2),round(yr2),repmat(loop,size(xr2)));
    mask(ind) = 1;
  end
end

%--- Smooth scar mask
% originalmask = mask;
% mask = econv3(mask,f);
% mask = econv3(mask,f');
% mask = econv3(mask,reshape(f,[1 1 n]));
    
%Interpolate to find scar region
cdata = interp3(mask,...
  fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3),'linear');

%Rescale
maxv = max(cdata(:));
% cdata(cdata<(maxv/2)) = 0;
if maxv==0
  maxv = 1;
end
cdata = 0.2+0.8*cdata/maxv;

%--- Scale vertices before visualization
fv.vertices(:,1) = fv.vertices(:,1)*SET(NO).ResolutionX;
fv.vertices(:,2) = fv.vertices(:,2)*SET(NO).ResolutionY;
fv.vertices(:,3) = fv.vertices(:,3)*(SET(NO).SliceThickness+SET(NO).SliceGap);

%--- Display it
fig = figure(12);
clf;
h = patch(fv);
set(h,'cdata',cdata,'facecolor','interp','facealpha',0.7,'edgealpha',0);
%set(h,'facecolor',[1 0 0],'edgecolor',[1 1 0],'facealpha',0.7);
%set(h,'facealpha',0.7);
%set(h,'edgealpha',0);
axis off image;
colormap(hot);
cameratoolbar(fig);
set(gca,'clim',[0 1]);

