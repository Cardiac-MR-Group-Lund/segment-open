function varargout = buttonupfunctions(varargin)
% Functions for buttonups

% Broken out by Klas

%Invoke subfunction
%#ok<*GVMIS>
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-----------------------------------
function glarotatehandle_Buttonup(panel)
%-----------------------------------
%Update view according to new rotation
global SET DATA

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

x = x/scale;
y = y/scale;

glaangle = mod(atan2(SET(no).ResolutionX*(y-SET(no).HLA.slice), ...
  SET(no).ResolutionY*(x-SET(no).VLA.slice))+pi/2,pi)-pi/2;
SET(no).GLA.angle = glaangle;

%also do gla center update this updates the x0 y0 coordinates. Which is
%needed for calculating the plane intersection.
viewfunctions('setglacenter',no);

set(DATA.fig,'WindowButtonMotionFcn',[], ...
  'WindowButtonUpFcn',[]);

DATA.Handles.orthoanglehandle.XData = DATA.Handles.cursor.XData(2);
DATA.Handles.orthoanglehandle.YData = DATA.Handles.cursor.YData(2);

%clear cursor
DATA.Handles.cursor.XData = nan;
DATA.Handles.cursor.YData = nan;
DATA.Handles.cursor.Marker = 'none';

%update the ortho panel and the gla panel intersection
drawfunctions('drawplaneintersections');

%and gla panel
glapanel = find(strcmp(DATA.ViewPanelsType,'gla'));

%update gla segmentation intersections
for t = 1:SET(no).TSize
  calcfunctions('segmentationintersection_helper',glapanel,t) %this updates stored segmentation intersections for remaining timeframes no graphical update needed.
end

DATA.ViewIM{glapanel} = []; %forces viewim update of the glapanel
drawfunctions('drawpanel',glapanel);

%------------------------------
function crop_buttonup(panel,xstart,ystart,slice)
%-----------------------
%This function adds the interpolated contour to all timeframes if all
%timeframes mode is selected it also removes motion and buttonupfunction
%settings in DATA.fig

global DATA SET
no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

[Xout, Yout] = motionfunctions('crop_motion',panel,xstart,ystart);

%get translation from montage case
slices = viewfunctions('slicesinpanel',panel);
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice));
imdim = zoomfunctions.getxysize(no,panel);
xt = (yl-1)*imdim.YSize - imdim.YStart +1;
yt = (xl-1)*imdim.XSize - imdim.XStart +1;

%this removes any translation from montage and scaling so the format is the
%same as we store coordinates in.
xstart = xstart/scale-xt;
ystart = ystart/scale-yt;
xlast = Xout(3)/scale-xt;
ylast = Yout(3)/scale-yt;

%Get min max points in box and provide to crophelper
yind = max(round([1,min([xlast,xstart])])):min(round([SET(no).YSize,max([xlast,xstart])]));
xind = max(round([1,min([ylast,ystart])])):min(round([SET(no).XSize,max([ylast,ystart])]));

%remove buttonup and motionfunctions
buttonup_Callback

qststr = dprintf('Apply crop?');
answ = strainmitt.strainmitt('existstraininnos',no);
if answ
  qststr = dprintf('Strain calculation exist in the stack %d and will be deleted after cropping. Apply crop?',no);
end

if ~yesno(qststr,[],DATA.GUI.Segment)
  %reset cursor
  DATA.Handles.cursor.XData = nan;
  DATA.Handles.cursor.YData = nan;
  return;
end

%reset cursor
DATA.Handles.cursor.XData = nan;
DATA.Handles.cursor.YData = nan;

tools('crophelper',no,xind,yind); %after this we are done

tools('disableundo')
DATA.ViewIM{DATA.CurrentPanel}=[];
drawfunctions('drawno',no)
drawfunctions('drawselectedframe',panel)

%also update thumbnail!
drawfunctions('drawthumbnails',1,0)

if ~isempty(SET(no).StrainMitt)
  strainmitt.strainmitt('clearstrain_helper',no);
  %   type = 'all';
  %   updateredo(SET(no).StrainMitt,type);
end

%------------------------------
function interp_buttonup(panel,type,slice)
%-----------------------
%This function adds the interpolated contour to all timeframes if all
%timeframes mode is selected it also removes motion and buttonupfunction
%settings in DATA.fig
global DATA SET

[no,~,ctf,~,~] = helperfunctions('getinterpparameters',panel);

%Check object ind in case of General Pen
currentobjectind = 1;
if strcmp(type,'GeneralPenInterp')
  currentobjectind = DATA.GeneralPenSettings.getcurrentobject;
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

%This applies the contours and interpolation points to all timeframes
redostrainmitt = false;
if not(DATA.ThisFrameOnly) %all frames mode
  redostrainmitt = true;
  %Get contour
  typeX = helperfunctions('parsesetfield',SET(no),type(1:end-6),'X',currentobjectind); %get contour
  typeY = helperfunctions('parsesetfield',SET(no),type(1:end-6),'Y',currentobjectind); %get contour
  typeInterpX = helperfunctions('parsesetfield',SET(no),type,'X',currentobjectind); %get interp contour
  typeInterpY = helperfunctions('parsesetfield',SET(no),type,'Y',currentobjectind); %get interp contour
  for tf = 1:SET(no).TSize
    %Assign contour at current time frame to all timeframes
    helperfunctions('assignsetfield',no,type(1:end-6),'X',typeX(:,ctf,slice),currentobjectind,tf,slice);
    helperfunctions('assignsetfield',no,type(1:end-6),'Y',typeY(:,ctf,slice),currentobjectind,tf,slice);
    helperfunctions('assignsetfield',no,type,'X',typeInterpX{ctf,slice},currentobjectind,tf,slice);
    helperfunctions('assignsetfield',no,type,'Y',typeInterpY{ctf,slice},currentobjectind,tf,slice);
  end
end

if (redostrainmitt || isequal(ctf,SET(no).EDT)) && ~isempty(SET(no).StrainMitt)
  updateredo(SET(no).StrainMitt,type);
end

buttonup_Callback;
isinterpongoing = helperfunctions('isinterpongoing',SET(no),type,currentobjectind);
if isinterpongoing
  drawfunctions('drawinterp',panel,type);
else
  drawfunctions('drawno',no);
end

calcfunctions('updatemarandscar',no);

%update result panel
switch type
  case {'EndoInterp','EpiInterp'}
    if ~SET(no).([type(1:end-6),'InterpOngoing'])
      if isempty(DATA.LVNO)
        DATA.LVNO = no;
      end
      segment('updatevolume')
    end
  case {'RVEndoInterp','RVEpiInterp'}
    if ~SET(no).([type(1:end-6),'InterpOngoing'])
      if isempty(DATA.RVNO)
        DATA.RVNO = no;
      end
      segment('updatevolume')
    end
end

if (strcmp(type(1:end-6),'LA') || strcmp(type(1:end-6),'RA')) && ismember(ctf,[SET(no).EST, SET(no).EDT]) 
  segment('updatemeasurement');
  if isequal(ctf,SET(no).EST)
    SET(no).(type(1:end-6)).isManualContourinES = true;
  end
end

%------------------------------
function scale_buttonup(panel,type,objectind,slice,startrad)
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);
if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  magno = SET(no).Flow.MagnitudeNo;
else
  magno = no;
end
scale = viewfunctions('getscale',panel);

xc = mean(DATA.CursorX);
yc = mean(DATA.CursorY);


%Update cursor coordinates then right back into SET struct
[DATA.CursorX,DATA.CursorY, scalevalue] = motionfunctions('scale_motion',panel,startrad);
switch type
  case {'Endo','Epi'}
    stepsize = (scalevalue-1)/0.02;
    lv('segmentexpandcontract_Callback',stepsize,lower(type),slice);

    %If straintagging initiated adjust LVupdated
    if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
      SET(no).StrainTagging.LVupdated = 1;
    end
    if ~isempty(SET(no).StrainMitt)
      updateredo(SET(no).StrainMitt,type);
    end
    calcfunctions('updatemarandscar',no);
  case {'RVEndo','RVEpi'}
    stepsize = (scalevalue-1)/0.02;
    rv('expandcontract_Callback',stepsize,lower(type),slice);

  case 'Roi'
    if DATA.ThisFrameOnly  %not(findindented(DATA.Handles.hideiconholder,'allframesmode'))
      slices = viewfunctions('slicesinpanel',panel);
      [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
      imdim = zoomfunctions.getxysize(no,panel);

      xt = (xl-1)*imdim.XSize - imdim.XStart +1;
      yt = (yl-1)*imdim.YSize - imdim.YStart +1;

      %single frame mode
      DATA.DoThisFrameOnly = true; %flag for calc functions
      xx = DATA.CursorX/scale - xt;
      yy = DATA.CursorY/scale - yt;
      [xx,yy] = calcfunctions('correctoutboundcontour',xx,yy,no);

      SET(magno).Roi(objectind).X(:,SET(no).CurrentTimeFrame) = xx;
      SET(magno).Roi(objectind).Y(:,SET(no).CurrentTimeFrame) = yy;
    else
      tf = 1:SET(magno).TSize;
      mx = repmat(mean(SET(magno).Roi(SET(magno).RoiCurrent).X(:,tf)),size(SET(magno).Roi(SET(magno).RoiCurrent),1),1);
      my = repmat(mean(SET(magno).Roi(SET(magno).RoiCurrent).Y(:,tf)),size(SET(magno).Roi(SET(magno).RoiCurrent),1),1);
      xx = bsxfun(@plus,mx,scalevalue*bsxfun(@minus,SET(magno).Roi(SET(magno).RoiCurrent).X(:,tf),mx));
      yy = bsxfun(@plus,my,scalevalue*bsxfun(@minus,SET(magno).Roi(SET(magno).RoiCurrent).Y(:,tf),my));
      [xx,yy] = calcfunctions('correctoutboundcontour',xx,yy,magno);

      SET(magno).Roi(objectind).X(:,tf) = xx;
      SET(magno).Roi(objectind).Y(:,tf) = yy;
    end
    %KG: update Area and Mean of the ROI
    [~,SET(magno).Roi(objectind).Area] = ...
      calcfunctions('calcroiarea',magno,objectind);
    [m,sd]=calcfunctions('calcroiintensity',magno,objectind);
    SET(magno).Roi(objectind).Mean = m;
    SET(magno).Roi(objectind).StD = sd;
    %Update Flow ROI
    segment('updateflow');
    DATA.DoThisFrameOnly = false; %flag for calc functions
end

drawfunctions('drawno',no)

DATA.Handles.cursor.YData = nan;
DATA.Handles.cursor.XData = nan;
DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

%----------------------------------------------
function contrast_buttonup(panel,slice,xstart,ystart,xsize,ysize)
%----------------------------------------------
global DATA SET
no = DATA.ViewPanels(panel);
[contrast,brightness] = motionfunctions('contrast_motion',panel,slice,xstart,ystart,xsize,ysize);
SET(no).IntensityMapping.Contrast = contrast;
SET(no).IntensityMapping.Brightness = brightness;
DATA.ViewIM{panel} = [];
drawfunctions('drawimages',panel);

DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

%------------------------------
function translate_buttonup(panel,type,objectind,xstart,ystart,slice)
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);
if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  magno = SET(no).Flow.MagnitudeNo;
else
  magno = no;
end
scale = viewfunctions('getscale',panel);

slices = viewfunctions('slicesinpanel',panel);
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
imdim = zoomfunctions.getxysize(no,panel);

xt = ((xl-1)*imdim.XSize- imdim.XStart +1);
yt = ((yl-1)*imdim.YSize- imdim.YStart +1);

%Update cursor coordinates then right back into SET struct
[DATA.CursorX,DATA.CursorY] = motionfunctions('translate_motion',panel,xstart,ystart);

x = DATA.CursorX/scale - xt;
y = DATA.CursorY/scale - yt;

%check if any part of the contour is out of bounds then place it along the
%edge of the image.
x = min(x,SET(no).XSize-2);
x = max(x,2);

y = min(y,SET(no).YSize-2);
y = max(y,2);

switch type
  case {'Endo','Epi','RVEndo','RVEpi'}
    %objectind is not used for contours
    if not(DATA.ThisFrameOnly)  %findindented(DATA.Handles.hideiconholder,'allframesmode')
      %all frames mode
      diffx = SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice) - DATA.CursorX/scale + xt;
      diffy = SET(no).([type,'Y'])(:,SET(no).CurrentTimeFrame,slice) - DATA.CursorY/scale + yt;

      newx = SET(no).([type,'X'])(:,:,slice) - diffx;
      newx = min(newx,SET(no).XSize-2);
      newx = max(newx,2);
      newy = SET(no).([type,'Y'])(:,:,slice) - diffy;

      newy = min(newy,SET(no).YSize-2);
      newy = max(newy,2);
      SET(no).([type,'X'])(:,:,slice) = newx;%        repmat(x,[1,SET(no).TSize]);
      SET(no).([type,'Y'])(:,:,slice) = newy;%        repmat(y,[1,SET(no).TSize]);
    else
      SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice) = x;
      SET(no).([type,'Y'])(:,SET(no).CurrentTimeFrame,slice) = y;
    end

    %If straintagging initiated adjust LVupdated
    if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
      SET(no).StrainTagging.LVupdated = 1;
    end
    if ~isempty(SET(no).StrainMitt)
      updateredo(SET(no).StrainMitt,type);
    end
    calcfunctions('updatemarandscar',no);
    %update volumes in case of cropped by border
    segment('updatevolume')

  case 'Roi'
    if DATA.ThisFrameOnly  %~findindented(DATA.Handles.hideiconholder,'allframesmode')
      %single frame mode
      DATA.DoThisFrameOnly = true; %flag for calc functions
      SET(magno).Roi(objectind).X(:,SET(no).CurrentTimeFrame) = x;
      SET(magno).Roi(objectind).Y(:,SET(no).CurrentTimeFrame) = y;
    else
      diffx = SET(magno).Roi(objectind).X(:,SET(no).CurrentTimeFrame) - DATA.CursorX/scale + xt;
      diffy = SET(magno).Roi(objectind).Y(:,SET(no).CurrentTimeFrame) - DATA.CursorY/scale + yt;
      newx = SET(magno).Roi(objectind).X - diffx;
      newx = min(newx,SET(no).XSize-2);
      newx = max(newx,2);
      newy = SET(magno).Roi(objectind).Y - diffy;

      newy = min(newy,SET(no).YSize-2);
      newy = max(newy,2);
      SET(magno).Roi(objectind).X = newx;%repmat(x,[1,SET(no).TSize]);
      SET(magno).Roi(objectind).Y = newy;%repmat(y,[1,SET(no).TSize]);
    end

    %update mean area std in case of cropped by border
    [~,SET(magno).Roi(objectind).Area] = ...
      calcfunctions('calcroiarea',magno,objectind);
    [m,sd]=calcfunctions('calcroiintensity',magno,objectind);
    SET(magno).Roi(objectind).Mean = m;
    SET(magno).Roi(objectind).StD = sd;

  case 'Measure'
    measurex = x;
    measurey = y;
    switch DATA.ViewPanelsType{panel}
      case 'hla'
        SET(no).Measure(objectind).X = ones(size(measurey))*SET(no).HLA.slice;
        SET(no).Measure(objectind).Y = measurey;
        SET(no).Measure(objectind).Z = measurex;
      case 'vla'
        SET(no).Measure(objectind).X = measurey;
        SET(no).Measure(objectind).Y = ones(size(measurey))*SET(no).VLA.slice;
        SET(no).Measure(objectind).Z = measurex;
      case 'gla'
        [SET(no).Measure(objectind).X,SET(no).Measure(measureN).Y,...
          SET(no).Measure(objectind).Z] = calcfunctions('gla2sax',measurex,measurey,no);
      otherwise
        SET(no).Measure(objectind).X = measurex;
        SET(no).Measure(objectind).Y = measurey;
        SET(no).Measure(objectind).Z = ones(size(measurey))*SET(no).CurrentSlice;
    end
    %Calc length in case of cropped by border
    L = sum(sqrt(...
      (SET(no).ResolutionX*diff(measurex)).^2+...
      (SET(no).ResolutionY*diff(measurey)).^2+...
      ((SET(no).SliceThickness+SET(no).SliceGap)*diff(SET(no).Measure(objectind).Z)).^2));

    SET(no).Measure(objectind).Length = L;

  case 'Point'
    SET(no).Point.X(objectind) = x;
    SET(no).Point.Y(objectind) = y;
  case 'Center'
    SET(no).CenterX = x;
    SET(no).CenterY = y;
end

%draw changes
drawfunctions('drawno',no)

%update result panel
switch type
  case {'Endo','Epi','RVEndo','RVEpi'}
    segment('updatevolume');
  case 'Roi'
    segment('updateflow');
    DATA.DoThisFrameOnly = false; %flag for calc functions
  otherwise
    segment('updatemeasurement');
end

DATA.Handles.cursor.YData = nan;
DATA.Handles.cursor.XData = nan;

DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';


%------------------------------
function translateall_buttonup(panel,xstart,ystart,slice)
%-----------------------
global DATA SET

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

%translate the contours by xstart-x transferred to SET contour domain.
tx = (xstart-x)/scale;
ty = (ystart-y)/scale;



types = {'Endo','Epi','RVEndo','RVEpi'};

for i = 1:length(types)
  type = types{i};
  if ~isempty(SET(no).([type,'X'])) && ~all(isnan(SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice)))
    contourx = SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice)-tx;
    contoury = SET(no).([type,'Y'])(:,SET(no).CurrentTimeFrame,slice)-ty;

    %check if any part of the contour is out of bounds then place it along the
    %edge of the image.
    contourx = min(contourx,SET(no).XSize-2);
    contourx = max(contourx,2);
    contoury = min(contoury,SET(no).YSize-2);
    contoury = max(contoury,2);

    SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice) = contourx;
    SET(no).([type,'Y'])(:,SET(no).CurrentTimeFrame,slice) = contoury;
  end
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

if ~isempty(SET(no).StrainMitt)
  updateredo(SET(no).StrainMitt,type);
end

calcfunctions('updatemarandscar',no);
drawfunctions('drawno',no)

DATA.CursorX = [];
DATA.CursorY = [];
DATA.Handles.cursor.YData = nan;
DATA.Handles.cursor.XData = nan;

DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

%----------------------------------------------
function point_buttonup(panel,pointind,slice)
%----------------------------------------------
%Point button up function. Has the current panel and index of the
%dragged point we are updating as input.

global DATA SET

if nargin == 2
  slice = [];
end

no = DATA.ViewPanels(panel);
%Restore/Kill all action functions
DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

%check if placed in 3dp view then we need to store it in a different way
if any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','speedim','cor3DP'}))

  [xsz,ysz,~,~] = segment3dp.tools('viewsizeandres',SET(no).LevelSet.Pen.Color,no);

  px = min(x,xsz-2);
  px = max(px,2);
  py = min(y,ysz-2);
  py = max(py,2);

  [r,g,b] = segment3dp.tools('xy2rgb',SET(no).LevelSet.Pen.Color,px,py);
  [x,y,z] = segment3dp.tools('rgb2xyz',r,g,b);

  %Add the point
  segment3dp.tools('addpoint_helper',x,y,z,no)

  for p = 1:length(DATA.ViewPanels)
    if ~isequal(DATA.ViewPanelsType{p},'viewport')
      viewfunctions('updatedrawlist',p)
      drawfunctions('drawpanel',p)
    end
  end

else
  %Not 3DP
  slices = viewfunctions('slicesinpanel',panel);
  scale = viewfunctions('getscale',panel);
  %get montage coordinates of clicked position
  [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice));
  imdim = zoomfunctions.getxysize(no,panel);
  xt = (xl-1)*imdim.XSize - imdim.XStart +1;
  yt = (yl-1)*imdim.YSize - imdim.YStart +1;

  %this removes any translation from montage and scaling so the format is the
  %same as we store coordinates in.
  px = x/scale -xt;
  py = y/scale -yt;

  %check if any part of the contour is out of bounds then place it along the
  %edge of the image.
  px = min(px,SET(no).XSize-2);
  px = max(px,2);

  py = min(py,SET(no).YSize-2);
  py = max(py,2);

  SET(no).Point.X(pointind) = px;
  SET(no).Point.Y(pointind) = py;
  SET(no).Point.Z(pointind) = slice;
  if length(SET(no).Point.Label)<pointind
    SET(no).Point.Label{pointind} = sprintf('P%d',pointind);
    if strcmpi(DATA.CurrentTheme,'strain') && strcmpi(SET(no).ImageViewPlane, '4CH')
      pointname = 'TV plane';
      if (sum(contains(SET(no).Point.Label, pointname))<2)
        SET(no).Point.Label{pointind} = pointname;
      end
    end
    if strcmpi(DATA.CurrentTheme,'strain') && ...
        (strcmpi(SET(no).ImageViewPlane, '2CH') || strcmpi(SET(no).ImageViewPlane, '3CH'))...
        && (sum(contains(SET(no).Point.Label, 'MV plane'))<2)
      SET(no).Point.Label{pointind} = 'MV plane';
    end
  end

  if not(DATA.ThisFrameOnly)
    %all frames mode
    SET(no).Point.T(pointind) = nan;
  else
    SET(no).Point.T(pointind) = SET(no).CurrentTimeFrame;
  end

  drawfunctions('drawno',no);

  %Check if strain need to be recalculated
  if strcmp(SET(no).Point.Label{pointind},'MV plane') && ~isempty(SET(no).StrainMitt)
    updateredo(SET(no).StrainMitt,'LV');
    updateredo(SET(no).StrainMitt,'LA');
  end
  if strcmp(SET(no).Point.Label{pointind},'TV plane') && ~isempty(SET(no).StrainMitt)
    updateredo(SET(no).StrainMitt,'RV');
    updateredo(SET(no).StrainMitt,'RA');
  end

end

DATA.Handles.cursor.XData = nan;
DATA.Handles.cursor.YData = nan;

%----------------------------------------------
function measuretext_buttonup(panel,pointind)
%----------------------------------------------
% Button-up function for measures' text
global DATA
drawfunctions('drawmeasures',panel)
%Restore/Kill all action functions
DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

%----------------------------------------------
function measure_buttonup(panel,measureN,slice)
%----------------------------------------------
%Measure button up function. Has the current panel index of the currently
%dragged point and finally the index of the measure we are updating as
%input.

global DATA SET

no = DATA.ViewPanels(panel);
%Restore/Kill all action functions
DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

%get the slice and normalize the coordinates
scale = viewfunctions('getscale',panel);

slices = viewfunctions('slicesinpanel',panel);

%get montage coordinates of clicked position
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
if isempty(yl)
  yl = 1;
  xl = 1;
end
imdim = zoomfunctions.getxysize(no,panel);
xt = (xl-1)*imdim.XSize - imdim.XStart +1;
yt = (yl-1)*imdim.YSize - imdim.YStart +1;
measurey = DATA.CursorX/scale-yt;
measurex = DATA.CursorY/scale-xt;

%check if any part of the contour is out of bounds then place it along the
%edge of the image.
measurex = min(measurex,SET(no).XSize-2);
measurex = max(measurex,2);

measurey = min(measurey,SET(no).YSize-2);
measurey = max(measurey,2);

switch DATA.ViewPanelsType{panel}
  case 'hla'
    SET(no).Measure(measureN).X = ones(size(measurey))*SET(no).HLA.slice;
    SET(no).Measure(measureN).Y = measurey;
    SET(no).Measure(measureN).Z = measurex;
  case 'vla'
    SET(no).Measure(measureN).X = measurey;
    SET(no).Measure(measureN).Y = ones(size(measurey))*SET(no).VLA.slice;
    SET(no).Measure(measureN).Z = measurex;
  case 'gla'
    [SET(no).Measure(measureN).X,SET(no).Measure(measureN).Y,...
      SET(no).Measure(measureN).Z] = calcfunctions('gla2sax',measurex,measurey,no);
  otherwise
    SET(no).Measure(measureN).X = measurex;
    SET(no).Measure(measureN).Y = measurey;
    SET(no).Measure(measureN).Z = DATA.CursorZ;%ones(size(measurey))*SET(no).CurrentSlice;
end

%Calc length
L = sum(sqrt(...
  (SET(no).ResolutionX*diff(measurex)).^2+...
  (SET(no).ResolutionY*diff(measurey)).^2+...
  ((SET(no).SliceThickness+SET(no).SliceGap)*diff(SET(no).Measure(measureN).Z)).^2));

SET(no).Measure(measureN).Length = L;
SET(no).Measure(measureN).T = SET(no).CurrentTimeFrame;

if ~isfield(SET(no).Measure(measureN),'Name') || ~isfield(SET(no).Measure(measureN),'LongName')
  SET(no).Measure(measureN).Name = '';
  SET(no).Measure(measureN).LongName = '';
end

if ~isfield(SET(no).Measure(measureN),'Offset')
  SET(no).Measure(measureN).Offset = [10 0]; % default offset for text
end

DATA.Handles.cursor.XData = nan;
DATA.Handles.cursor.YData = nan;
DATA.CursorX = [];
DATA.CursorY = [];
DATA.CursorZ = [];

%Remove markers from cursor object
DATA.Handles.cursor.Marker = 'none';

if any(strcmp(DATA.ViewPanelsType{panel},{'orth','hla','gla','vla'}))
  for p = 1:length(DATA.ViewPanels)
    drawfunctions('drawmeasures',p)
  end
else
  drawfunctions('drawmeasures',panel)
end

buttondownfunctions('updatebuttondowns','Measure');
viewfunctions('updatedrawlist',panel);
drawfunctions('drawpanel',panel);
if DATA.issegmentcmr
  % show renaming tool directly
  if isempty(SET(no).Measure(measureN).Name) || isempty(SET(no).Measure(measureN).LongName)
    callbackfunctions('measurerenamethis_Callback',measureN)
  end
end
segment('updatemeasurement');


%------------------------------------------------------------------------
function select_buttonup
%------------------------------------------------------------------------
%Function that highlights the selected panel if montage also highlight the frame
global DATA

DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';
drawfunctions('drawplaneintersections')


%----------------------------------------------
function balloon_buttonup(panel,type,slice,tf)
%----------------------------------------------
%Generic buttonup that only clears motion and buttonup function in fig.
global DATA SET

no = DATA.ViewPanels(panel);

if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  no = SET(no).Flow.MagnitudeNo;
end

y = DATA.CursorX;
x = DATA.CursorY;
if isrow(x)
  x = x';
end
if isrow(y)
  y = y';
end

if isempty(x)
  DATA.Handles.cursor.XData = nan;
  DATA.Handles.cursor.YData = nan;
  DATA.CursorX = [];
  DATA.CursorY = [];

  DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
  DATA.fig.WindowButtonUpFcn = '';
  return
end

%Determine which slice we are drawing in and interpolation mode to translate the cursor contour .
scale = viewfunctions('getscale',panel);
slices = viewfunctions('slicesinpanel',panel);
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice));
imdim = zoomfunctions.getxysize(no,panel);
xt = (xl-1)*imdim.XSize*scale - imdim.XStart +1;
yt = (yl-1)*imdim.YSize*scale - imdim.YStart +1;
%this removes any translation from montage and scaling so the format is the
%same as we store coordinates in.
x = (x-xt)/scale;
y = (y-yt)/scale;

%Adjust curve direction.
[x,y]=mypoly2cw(x,y);

switch type
  case 'Roi'
    %we do not merge rois
    x_comp = [];
    y_comp = [];
  otherwise
    if ~isempty(SET(no).([type,'X']))
      x_comp = SET(no).([type,'X'])(:,tf,slice);
      y_comp = SET(no).([type,'Y'])(:,tf,slice);
    else
      x_comp = [];
      y_comp = [];
    end
end


%this will trigger if exists segmentation in slice and timeframe
if ~isempty(x_comp) && ~all(isnan(x_comp))
  %determine if polygon is inside other polygon. If all is inside we keep the
  % new contour.
  mnew = poly2mask(x,y,SET(no).XSize,SET(no).YSize);
  mold = poly2mask(x_comp,y_comp,SET(no).XSize,SET(no).YSize);

  %this means that there is overlap between the old balloon and the new and an addition in the new blob that we need to
  %incorporate
  if any(mnew(:).*mold(:)>0) && any(mnew(:)-mold(:)>0)
    m = mold+mnew>0;
    [mi,mj] = ind2sub(size(m),find(m,1));
    X = bwtraceboundary(m,[mi,mj],'W');

    xm = zeros(length(X),1);
    ym = zeros(length(X),1);

    %idea is to find the closest point of all points which are not
    %inside a polygon
    [in1,~] = inpolygon(x,y,x_comp,y_comp);
    [in2,~] = inpolygon(x_comp,y_comp,x,y);
    x_cat = cat(1,x(~in1),x_comp(~in2));
    y_cat = cat(1,y(~in1),y_comp(~in2));

    %find the closest points on the combined polygon of the contour and
    %the drawn. Then let resamplecurve work removing any duplicates.
    for i = 1:length(X)
      [~,ind] = min((X(i,2)-x_cat).^2+(X(i,1)-y_cat).^2);
      xm(i) = x_cat(ind);
      ym(i) = y_cat(ind);
    end
    x = xm;
    y = ym;
  end
end

%This might fail for fast clicks where xy contains duplicates or is empty.
%if so clear the cursor object and motion buttonupfunctions and return before writing anything back to SET.
try
  % check for duplicate points
  nonduplicates = unique([x,y],'rows','stable');
  % remove duplicate points
  x = nonduplicates(:,1);
  y = nonduplicates(:,2);
  %removes duplicate points and resamples the contour
  [x,y] = calcfunctions('resamplecurve',x,y,DATA.NumPoints-1);

  %Close the contour
  x = [x,x(1)];
  y = [y,y(1)];

catch
  DATA.Handles.cursor.XData = nan;
  DATA.Handles.cursor.YData = nan;
  DATA.CursorX = [];
  DATA.CursorY = [];
  DATA.Handles.cursor.LineStyle = '-';

  DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
  DATA.fig.WindowButtonUpFcn = '';
  return
end

switch type
  case 'Roi'
    %For the balloon case the we always create a new Roi
    emptyroiinds = find(cellfun(@isempty,{SET(no).Roi.X}));
    if ~isempty(emptyroiinds)
      roiind = emptyroiinds(1);
    else
      roiind =length(SET(no).Roi)+1;
    end

    if not(DATA.ThisFrameOnly)  %findindented(DATA.Handles.hideiconholder,'allframesmode')
      %all frames mode
      SET(no).Roi(roiind).X = repmat(x',[1,SET(no).TSize]);
      SET(no).Roi(roiind).Y = repmat(y',[1,SET(no).TSize]);
    else
      DATA.DoThisFrameOnly = true; %flag for calc functions
      tmpx = nan(length(x),SET(no).TSize);
      tmpy = nan(length(x),SET(no).TSize);
      tmpx(:,tf) = x;
      tmpy(:,tf) = y;
      SET(no).Roi(roiind).X = tmpx;
      SET(no).Roi(roiind).Y = tmpy;
    end

    SET(no).Roi(roiind).T = 1:SET(no).TSize;
    SET(no).Roi(roiind).Z = slice;

    %Set roicurrent
    SET(no).RoiCurrent = roiind;

    %Increase number of rois
    SET(no).RoiN = SET(no).RoiN + 1;

    %Update roi name

    SET(no).Roi(roiind).Sign = 1;
    SET(no).Roi(roiind).Name = roi('roisetnewname',no);
    newlinecolor = roi('roisetcolorbasedonposition',no);
    SET(no).Roi(roiind).LineSpec = [newlinecolor,'-'];

    %Calculate area and intensity of ROI
    [~, roiarea] = calcfunctions('calcroiarea',no,roiind);
    threshold = 0.01*(SET(no).ResolutionX *SET(no).ResolutionY);% area of 1 pixel in cm^2
    ind = ~isnan(roiarea);
    isvalidroi = all(roiarea(ind) >= threshold);
    SET(no).Roi(roiind).Area = roiarea;
      
    [m,sd]=calcfunctions('calcroiintensity',no,roiind);
    SET(no).Roi(roiind).Mean = m;
    SET(no).Roi(roiind).StD = sd;

    %Update result panel
    if isvalidroi
      segment('updateflow');
    end
    DATA.DoThisFrameOnly = false; %flag for calc functions
    if ~isvalidroi
      roi('roidelete_Callback');
    end

  otherwise

    if strcmp(type, 'RVEndo')
      [x,y] = calcfunctions('calcpointsoutsideLV',x,y,no);
    end

    if isempty(SET(no).([type,'X']))
      SET(no).([type,'Y']) = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
      SET(no).([type,'X']) = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
    end

    if not(DATA.ThisFrameOnly)  %findindented(DATA.Handles.hideiconholder,'allframesmode')
      %all frames mode
      SET(no).([type,'Y'])(:,:,slice)=repmat(y',[1,SET(no).TSize]);
      SET(no).([type,'X'])(:,:,slice)=repmat(x',[1,SET(no).TSize]);
    else
      SET(no).([type,'Y'])(:,tf,slice)=y;
      SET(no).([type,'X'])(:,tf,slice)=x;
    end
    %If straintagging initiated adjust LVupdated
    if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
      SET(no).StrainTagging.LVupdated = 1;
    end

    if ~isempty(SET(no).StrainMitt)
      updateredo(SET(no).StrainMitt,type);
    end

    switch type
      case {'Endo','Epi'}
        if isempty(DATA.LVNO)
          DATA.LVNO = no;
        end
      case {'RVEndo','RVEpi'}
        if isempty(DATA.RVNO)
          DATA.RVNO = no;
        end
    end

    calcfunctions('updatemarandscar',no);
    segment('updatevolume')
end

% %generate image from normals of the endo contour
% x = SET(no).([type,'Y'])(:,tf,slice);
% y = SET(no).([type,'X'])(:,tf,slice);
%
% %use 100 points
% [x,y] = calcfunctions('resamplecurve',x,y,150);%round(DATA.NumPoints/2));
% x = x';
% y = y';
%
% dy = conv([x(end-1);x],[1 0 -1]/2,'valid');%diff(SET(no).([type,'Y'])(:,tf,slice),);
% dx = conv([y(end-1);y],[1 0 -1]/2,'valid');%diff(SET(no).([type,'X'])(:,tf,slice));
%
% %get normals
% normals = -[0 -1 ; 1 0]*[dx,dy]'./sqrt(dx.*dx+dy.*dy)';
% X = [y(1:end-1),x(1:end-1)]';
% %extract lines from along normals in image
% A = linspace(0,40,200);
% lines =X + normals.*reshape(repmat(A,2,1),2,1,200);
%
% imA = zeros(size(lines,2),size(lines,3));
%
% %image we retrieve values from
% imq = histeq(SET(no).IM(:,:,tf,slice),100);
%
% for i =1:size(lines,2)
%   xinds = round(squeeze(lines(1,i,:)));
%   yinds = round(squeeze(lines(2,i,:)));
%   inds = xinds>0 & xinds<=SET(no).YSize & yinds>0 & yinds<=SET(no).XSize;
%   xinds = xinds(inds);
%   yinds = yinds(inds);
%   indq = sub2ind(size(imq),xinds,yinds);
%   imA(i,inds) = imq(indq);
% end
%
% %we can bin the bloodpool content and fit a gaussian to this which we draw
% %the content of imA from.
% padim = padarray(imA,[2,2],'replicate');
% A = conv2(padim,ones(5)/25,'valid');
%
% %alternative idea
% %[~,medmed] = max(abs(diff(mean(A))));
% %medianres = ones(1,size(lines,2))*medmed;
%
% sample = A(:,1:6);
% mu = mean2(sample);
% sigma = std(sample(:));
% y = normpdf(A(:),mu,sigma);
% draw_imA = reshape(y,size(imA));
% level = graythresh(draw_imA);
%
% %multithresh
% %levels = multithresh(draw_imA,2);
% %a = single(draw_imA>levels(1) & draw_imA<levels(2));
%
% a = single(draw_imA>level);
% a(a==0) = -inf;
% t = cumsum(a,2);
%
% %take the local median or less around a row intervall
% medianres = max(t,[],2);
%
% medianres(isinf(medianres)) = nan;
% %medianres = medfilt1(maxvals,5);
% medmed = round(nanmedian(medianres));
%
% %just apply median for all
% medianres(:) = medmed;
%
% if isempty(SET(no).EpiX)
%   SET(no).EpiY = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
%   SET(no).EpiX = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
% end
%
% %apply this median index as the number of steps along the normal of the endo contour
% tmpy = zeros(1,size(lines,2));
% tmpx = zeros(1,size(lines,2));
%
% for i = 1:size(lines,2)
%   tmpy(i)=lines(2,i,medianres(i));
%   tmpx(i)=lines(1,i,medianres(i));
% end
%
% %apply this median index as the number of steps along the normal of the endo contour
% [x,y] = calcfunctions('resamplecurve',tmpx',tmpy',DATA.NumPoints-1);
%
% %Close the contour
% x = [x,x(1)];
% y = [y,y(1)];
%
% SET(no).EpiY(:,tf,slice)=y;
% SET(no).EpiX(:,tf,slice)=x;
%
% %do refine then smooth twice
% %epismoothsegmentation_helper(no,tf,slice)
% %lvpeter('segmentrefineepi_Callback')
% %epismoothsegmentation_helper(no,tf,slice)
%

% figure; plot(SET(no).EpiY(:,tf,slice),SET(no).EpiX(:,tf,slice),'g')
% hold on
% plot(SET(no).EndoY(:,tf,slice),SET(no).EndoX(:,tf,slice),'r')
% find
% figure; imagesc(draw_imA<level);



% figure;
% hold on
% for i =1:size(lines,2)
%   plot(squeeze(lines(1,i,:)),squeeze(lines(2,i,:)))
% end

% figure;
% hold on
% plot(SET(no).([type,'X'])(1:end-1,tf,slice)+normals(1,:)'*5,SET(no).([type,'Y'])(1:end-1,tf,slice)+normals(2,:)'*5,'*')
% plot(SET(no).([type,'X'])(:,tf,slice),SET(no).([type,'Y'])(:,tf,slice))

drawfunctions('drawno',no)

DATA.Handles.cursor.XData = nan;
DATA.Handles.cursor.YData = nan;
DATA.CursorX = [];
DATA.CursorY = [];
DATA.Handles.cursor.LineStyle = '-';

DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';


%----------------------------------------------
function buttonup_Callback
%----------------------------------------------
%Generic buttonup that only clears motion and buttonup function in fig.
global DATA

%Restore/Kill all action functions
DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

%----------------------------------------------
function pan_buttonup(panel)
%----------------------------------------------
global DATA SET
panelslinked = find(ismember(DATA.ViewPanels,SET(DATA.ViewPanels(panel)).Linked));
for actpanel = panelslinked
  no = DATA.ViewPanels(actpanel);

  %Update text position and show text again
  viewfunctions('updatetextposition',actpanel)
  scale = viewfunctions('getscale',actpanel);

  switch DATA.ViewPanelsType{actpanel}
    case 'trans3DP'
      SET(no).LevelSet.View.RZoomState(1:2) = DATA.Handles.imageaxes(actpanel).XLim/scale;
      SET(no).LevelSet.View.RZoomState(3:4) = DATA.Handles.imageaxes(actpanel).YLim/scale;
    case 'sag3DP'
      SET(no).LevelSet.View.GZoomState(1:2) = DATA.Handles.imageaxes(actpanel).XLim/scale;
      SET(no).LevelSet.View.GZoomState(3:4) = DATA.Handles.imageaxes(actpanel).YLim/scale;
    case 'cor3DP'
      SET(no).LevelSet.View.BZoomState(1:2) = DATA.Handles.imageaxes(actpanel).XLim/scale;
      SET(no).LevelSet.View.BZoomState(3:4) = DATA.Handles.imageaxes(actpanel).YLim/scale;
    otherwise
      %update the normalzoomstate field
      SET(no).NormalZoomState(1:2) = DATA.Handles.imageaxes(panel).XLim/scale;
      SET(no).NormalZoomState(3:4) = DATA.Handles.imageaxes(panel).YLim/scale;
  end

  %Draw all existing texts by calling drawpanel
  drawfunctions('drawpanel',actpanel)
end

buttonup_Callback

%----------------------------------------------
function pen_buttonup(panel,type,slice,tf,doall)
%----------------------------------------------
%New button up function input is type and panel. The types handled here
%are {Endo,Epi,RVendo,RVepi,Roi}.
global DATA SET

no = DATA.ViewPanels(panel);

if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  no = SET(no).Flow.MagnitudeNo;
end

if nargin<4
  tf = SET(no).CurrentTimeFrame;
end

if nargin<5
  doall = 0;
end

if nargin < 3
  slice = SET(no).CurrentSlice;
end

y = DATA.CursorX;
x = DATA.CursorY;

if isempty(x) || ismember(DATA.ViewPanelsType(panel),{'hla','vla','gla'})
  DATA.Handles.cursor.XData = nan;
  DATA.Handles.cursor.YData = nan;
  DATA.CursorX = [];
  DATA.CursorY = [];

  DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
  DATA.fig.WindowButtonUpFcn = '';
  return
end

%Adjust curve direction.
[x,y] = mypoly2cw(x,y);

%Determine which slice we are drawing in and interpolation mode to translate the cursor contour .
scale = viewfunctions('getscale',panel);
slices = viewfunctions('slicesinpanel',panel);
imdim = zoomfunctions.getxysize(no,panel);
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice));

xt = scale*((xl-1)*imdim.XSize- imdim.XStart +1);
yt = scale*((yl-1)*imdim.YSize- imdim.YStart +1);

%this removes any translation from montage and scaling so the format is the
%same as we store coordinates in.
x = (x-xt)/scale;
y = (y-yt)/scale;

%check if any part of the contour is out of bounds then place it along the
%edge of the image.
x = min(x,SET(no).XSize-2);
x = max(x,2);
y = min(y,SET(no).YSize-2);
y = max(y,2);

%two possibilities, either the segmentation overwrites the entire contour
%field in the SET struct or it is appended to it. this is determined in the
%function "appendcontour".
[x,y,newroi,roiind] = appendcontour(no,type,x,y,slice,tf,panel);

% check for duplicate points
nonduplicates = unique([x,y],'rows','stable');

if size(nonduplicates,1) < 2
  % if the non-duplicate points contains just 1 point or less
  % reset this drawing attempt
  DATA.Handles.cursor.XData = nan;
  DATA.Handles.cursor.YData = nan;
  DATA.CursorX = [];
  DATA.CursorY = [];

  DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
  DATA.fig.WindowButtonUpFcn = '';
  return
end
% remove duplicate points
x = nonduplicates(:,1);
y = nonduplicates(:,2);
% resamples the contour
numpoints = tools('getnumpointsforno',no);
opencontour = false;
if any(strcmp(type,{'LA','RA'}))
  opencontour = true;
end
if ~opencontour
  numpointstoresample = numpoints - 1;
else
  numpointstoresample = numpoints;
end
[x,y] = calcfunctions('resamplecurve',x,y,numpointstoresample,opencontour);
if ~opencontour
  x = [x,x(1)];
  y = [y,y(1)];
end

isonlyannotation = 1;

%write the results back to the SET struct
switch type
  case 'Roi'
    if not(DATA.ThisFrameOnly) || doall  % findindented(DATA.Handles.hideiconholder,'allframesmode')
      %all frames mode
      SET(no).Roi(roiind).X = repmat(x',[1,SET(no).TSize]);
      SET(no).Roi(roiind).Y = repmat(y',[1,SET(no).TSize]);
    elseif newroi
      DATA.DoThisFrameOnly = true; %flag for calc functions
      tmpx = nan(length(x),SET(no).TSize);
      tmpy = nan(length(x),SET(no).TSize);
      tmpx(:,tf) = x;
      tmpy(:,tf) = y;
      SET(no).Roi(roiind).X = tmpx;
      SET(no).Roi(roiind).Y = tmpy;
    else
      DATA.DoThisFrameOnly = true; %flag for calc functions
      SET(no).Roi(roiind).X(:,tf) = x;
      SET(no).Roi(roiind).Y(:,tf) = y;
    end

    SET(no).Roi(roiind).T = 1:SET(no).TSize;
    SET(no).Roi(roiind).Z = slice;

    %Set roicurrent
    SET(no).RoiCurrent = roiind;

    %Update Roi name
    if newroi
      SET(no).RoiN = SET(no).RoiN + 1;
      SET(no).Roi(roiind).Sign = 1;
      newname = roi('roisetnewname',no);
      SET(no).Roi(roiind).Name = newname;
      if strcmp(newname,'Static tissue')  || strcmp(newname,'Non-static tissue')
        clr = roi('roisetcolorbasedonname',no);
      else
        clr = roi('roisetcolorbasedonposition',no);
      end
      SET(no).Roi(roiind).LineSpec = [clr,'-'];
    end

    %Calculate area and intensity of ROI
    [~,SET(no).Roi(roiind).Area] = calcfunctions('calcroiarea',no,roiind);
    [m,sd]=calcfunctions('calcroiintensity',no,roiind);
    SET(no).Roi(roiind).Mean = m;
    SET(no).Roi(roiind).StD = sd;

  case {'Scar','MO'}
    %Check if it is atrial scar drawing, or normal LV scar

    %If atria scar then call special functions and exit otherwise continue
    if segment('doatrialscar',no)
      atrialscar('manualdraw_Buttonup',no,type,y,x);
      return;
    end

    if isempty(SET(no).Scar)
      viability('viabilityreset_Callback');
      SET(no).Scar.Mode = 'manual';
    end

    tempmask = segment('createmask',[SET(no).XSize SET(no).YSize],y,x) & SET(no).Scar.MyocardMask(:,:,slice);
    temp = SET(no).Scar.Manual(:,:,slice);

    if isequal(type,'Scar')
      temp(tempmask) = int8(1); %Mark manual scar as 1
    else
      temp(tempmask) = int8(2); %Mark manually no reflow as 2
    end

    SET(no).Scar.Manual(:,:,slice) = temp;

    if SET(no).Scar.UpdateDirectly
      force = true;  %force not to ask since this causes problem in manual drawing
      viability('viabilitycalc',force);
    end
    %--------------------------------------------------------------------------------------------------------------------------

  case {'CalciumPen'}
    %add code here
    isonlyannotation=0;
    if isempty(SET(no).CT.CalciumPenMask)
      viability('viabilityreset_Callback');
      SET(no).CT.CalciumPenMask = 'manual';
    end
    tempmask = uint8(segment('createmask',[SET(no).XSize SET(no).YSize],y,x) & SET(no).CT.CalciumMask(:,:,slice));

    sizesetstruct=size(SET);

    if sizesetstruct(1,2)==1
      m=1;
    else
      m=2;
    end
    temp=[];

    for n =1:sizesetstruct(1,m)
      temp=[temp isempty(SET(n).EndoX)==0];
    end

    if sum(temp)==1
      seg= find(temp);
    elseif sum(temp)==2
      findindex1=find(temp);
      seg = mymenu('Choose segmentation stack',['Image stack no ' num2str(SET(findindex1(1)).Linked)], ['Image stack no ' num2str(SET(findindex1(2)).Linked)]);
    end

    if SET(no).XSize==SET(seg).XSize
      SET(no).CT.CalciumPenMask(:,:,slice) = tempmask;
    else
      SET(no).CT.CalciumPenMask(:,:,slice) = imresize(tempmask, [SET(seg).XSize SET(seg).YSize]);
    end

    %   %  if SET(no).Scar.UpdateDirectly
    %       force = true;  %force not to ask since this causes problem in manual drawing
    %       viability('viabilitycalc',force);
    % %   %  end

    CalciumMask=SET(no).CT.CalciumMask(:,:,slice);

    fives=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==5);
    sixs=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==6);
    sevens=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==7);
    eights=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==8);

    if ~isempty(fives)
      CalciumMask(fives)=3;
    end
    if ~isempty(sixs)
      CalciumMask(sixs)=3;
    end
    if ~isempty(sevens)
      CalciumMask(sevens)=3;
    end
    if ~isempty(eights)
      CalciumMask(eights)=3;
    end

    SET(no).CT.CalciumMask(:,:,slice)=CalciumMask;

    DATA.ViewIM=[];

  case {'CalciumPenBlue'}
    %add code here
    isonlyannotation=0;
    if isempty(SET(no).CT.CalciumPenMask)
      viability('viabilityreset_Callback');
      SET(no).CT.CalciumPenMask = 'manual';
    end

    tempmask = uint8(segment('createmask',[SET(no).XSize SET(no).YSize],y,x) & SET(no).CT.CalciumMask(:,:,slice));

    sizesetstruct=size(SET);

    if sizesetstruct(1,2)==1
      m=1;
    else
      m=2;
    end

    temp=[];

    for n=1:sizesetstruct(1,m)
      temp=[temp isempty(SET(n).EndoX)==0];
    end

    if sum(temp)==1
      seg= find(temp);
    elseif sum(temp)==2
      findindex1=find(temp);
      seg = mymenu('Choose segmentation stack',['Image stack no ' num2str(SET(findindex1(1)).Linked)], ['Image stack no ' num2str(SET(findindex1(2)).Linked)]);
    end

    if SET(no).XSize==SET(seg).XSize
      SET(no).CT.CalciumPenMask(:,:,slice) = tempmask;
    else
      SET(no).CT.CalciumPenMask(:,:,slice) = imresize(tempmask, [SET(seg).XSize SET(seg).YSize]);
    end

    %   %  if SET(no).Scar.UpdateDirectly
    %       force = true;  %force not to ask since this causes problem in manual drawing
    %       viability('viabilitycalc',force);
    % %   %  end

    CalciumMask=SET(no).CT.CalciumMask(:,:,slice);

    threes=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==3);
    fives=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==5);
    sevens=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==7);
    eights=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==8);

    if ~isempty(threes)
      CalciumMask(threes)=6;
    end
    if ~isempty(fives)
      CalciumMask(fives)=6;
    end
    if ~isempty(sevens)
      CalciumMask(sevens)=6;
    end
    if ~isempty(eights)
      CalciumMask(eights)=6;
    end

    SET(no).CT.CalciumMask(:,:,slice)=CalciumMask;

    DATA.ViewIM=[];

  case {'CalciumPenRed'}
    %add code here
    isonlyannotation=0;
    if isempty(SET(no).CT.CalciumPenMask)
      viability('viabilityreset_Callback');
      SET(no).CT.CalciumPenMask = 'manual';
    end

    tempmask = uint8(segment('createmask',[SET(no).XSize SET(no).YSize],y,x) & SET(no).CT.CalciumMask(:,:,slice));

    sizesetstruct=size(SET);

    if sizesetstruct(1,2)==1
      m=1;
    else
      m=2;
    end

    temp=[];

    for n=1:sizesetstruct(1,m)
      temp=[temp isempty(SET(n).EndoX)==0];
    end

    if sum(temp)==1
      seg= find(temp);
    elseif sum(temp)==2
      findindex1=find(temp);
      seg = mymenu('Choose segmentation stack',['Image stack no ' num2str(SET(findindex1(1)).Linked)], ['Image stack no ' num2str(SET(findindex1(2)).Linked)]);
    end

    if SET(no).XSize==SET(seg).XSize
      SET(no).CT.CalciumPenMask(:,:,slice) = tempmask;
    else
      SET(no).CT.CalciumPenMask(:,:,slice) = imresize(tempmask, [SET(seg).XSize SET(seg).YSize]);
    end

    %   %  if SET(no).Scar.UpdateDirectly
    %       force = true;  %force not to ask since this causes problem in manual drawing
    %       viability('viabilitycalc',force);
    % %   %  end

    CalciumMask=SET(no).CT.CalciumMask(:,:,slice);

    threes=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==3);
    fives=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==5);
    sixs=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==6);
    eights=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==8);


    if ~isempty(threes)
      CalciumMask(threes)=7;
    end
    if ~isempty(fives)
      CalciumMask(fives)=7;
    end
    if ~isempty(sixs)
      CalciumMask(sixs)=7;
    end
    if ~isempty(eights)
      CalciumMask(eights)=7;
    end

    SET(no).CT.CalciumMask(:,:,slice)=CalciumMask;

    DATA.ViewIM=[];

  case {'CalciumPenLilac'}
    %add code here
    isonlyannotation=0;
    if isempty(SET(no).CT.CalciumPenMask)
      viability('viabilityreset_Callback');
      SET(no).CT.CalciumPenMask = 'manual';
    end

    tempmask = uint8(segment('createmask',[SET(no).XSize SET(no).YSize],y,x) & SET(no).CT.CalciumMask(:,:,slice));

    sizesetstruct=size(SET);

    if sizesetstruct(1,2)==1
      m=1;
    else
      m=2;
    end

    temp=[];

    for n=1:sizesetstruct(1,m)
      temp=[temp isempty(SET(n).EndoX)==0];
    end

    if sum(temp)==1
      seg= find(temp);
    elseif sum(temp)==2
      findindex1=find(temp);
      seg = mymenu('Choose segmentation stack',['Image stack no ' num2str(SET(findindex1(1)).Linked)], ['Image stack no ' num2str(SET(findindex1(2)).Linked)]);
    end

    if SET(no).XSize==SET(seg).XSize
      SET(no).CT.CalciumPenMask(:,:,slice) = tempmask;
    else
      SET(no).CT.CalciumPenMask(:,:,slice) = imresize(tempmask, [SET(seg).XSize SET(seg).YSize]);
    end

    %   %  if SET(no).Scar.UpdateDirectly
    %       force = true;  %force not to ask since this causes problem in manual drawing
    %       viability('viabilitycalc',force);
    % %   %  end

    CalciumMask=SET(no).CT.CalciumMask(:,:,slice);

    threes=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==3);
    fives=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==5);
    sixs=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==6);
    sevens=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==7);

    if ~isempty(threes)
      CalciumMask(threes)=8;
    end
    if ~isempty(fives)
      CalciumMask(fives)=8;
    end
    if ~isempty(sixs)
      CalciumMask(sixs)=8;
    end
    if ~isempty(sevens)
      CalciumMask(sevens)=8;
    end

    SET(no).CT.CalciumMask(:,:,slice)=CalciumMask;

    DATA.ViewIM=[];

  case {'CalciumPenRemove'}
    %add code here
    isonlyannotation=0;
    if isempty(SET(no).CT.CalciumPenMask)
      viability('viabilityreset_Callback');
      SET(no).CT.CalciumPenMask = 'manual';
    end

    tempmask = uint8(segment('createmask',[SET(no).XSize SET(no).YSize],y,x) & SET(no).CT.CalciumMask(:,:,slice));

    sizesetstruct=size(SET);

    if sizesetstruct(1,2)==1
      m=1;
    else
      m=2;
    end

    temp=[];

    for n=1:sizesetstruct(1,m)
      temp=[temp isempty(SET(n).EndoX)==0];
    end

    if sum(temp)==1
      seg= find(temp);
    elseif sum(temp)==2
      findindex1=find(temp);
      seg = mymenu('Choose segmentation stack',['Image stack no ' num2str(SET(findindex1(1)).Linked)], ['Image stack no ' num2str(SET(findindex1(2)).Linked)]);
    end

    if SET(no).XSize==SET(seg).XSize
      SET(no).CT.CalciumPenMask(:,:,slice) = tempmask;
    else
      SET(no).CT.CalciumPenMask(:,:,slice) = imresize(tempmask, [SET(seg).XSize SET(seg).YSize]);
    end

    %   %  if SET(no).Scar.UpdateDirectly
    %       force = true;  %force not to ask since this causes problem in manual drawing
    %       viability('viabilitycalc',force);
    % %   %  end

    threes=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==3);
    sixs=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==6);
    sevens=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==7);
    eights=find (single(tempmask).*single(SET(no).CT.CalciumMask(:,:,slice))==8);

    CalciumMask=SET(no).CT.CalciumMask(:,:,slice);
    if ~isempty(threes)
      CalciumMask(threes)=5;
    end
    if ~isempty(sixs)
      CalciumMask(sixs)=5;
    end
    if ~isempty(sevens)
      CalciumMask(sevens)=5;
    end
    if ~isempty(eights)
      CalciumMask(eights)=5;
    end

    SET(no).CT.CalciumMask(:,:,slice)=CalciumMask;

    DATA.ViewIM=[];

    %---------------------------------------------------------------------------------------------------------------------------------
  case 'ScarRubber'
    %Check if atrial scar then call special code, otherwise continue
    if segment('doatrialscar',no)
      atrialscar('manualdraw_Buttonup',no,type,y,x);
      return;
    end

    if isempty(SET(no).Scar)
      myfailed('You need to enter viability view mode before drawing infarct regions.',DATA.GUI.Segment);
      return;
    end

    %Create mask
    tempmask = segment('createmask',[SET(no).XSize SET(no).YSize],y,x)&SET(no).Scar.MyocardMask(:,:,slice);

    %Update scar
    temp = SET(no).Scar.Manual(:,:,slice);
    temp(tempmask) = int8(-1);
    SET(no).Scar.Manual(:,:,slice) = temp;

    %Update noreflow
    temp = SET(no).Scar.NoReflow(:,:,slice);
    temp(tempmask) = 0;
    SET(no).Scar.NoReflow(:,:,slice) = temp;

    if SET(no).Scar.UpdateDirectly
      viability('viabilitycalc');
    end

  case 'MaR'
    mar('createmyocardmask');
    SET(no).MaR.Mode = 'manual';
    tempmask = segment('createmask',[SET(no).XSize SET(no).YSize],y,x)&SET(no).MaR.MyocardMask(:,:,tf,slice);
    temp = SET(no).MaR.Manual(:,:,tf,slice);
    temp(tempmask) = int8(1); %Mark manual scar as 1
    SET(no).MaR.Manual(:,:,tf,slice) = temp;
    mar('update');

  case 'MaRRubber'
    mar('createmyocardmask');
    tempmask = segment('createmask',[SET(no).XSize SET(no).YSize],y,x)&SET(no).MaR.MyocardMask(:,:,tf,slice);
    temp = SET(no).MaR.Manual(:,:,tf,slice);
    temp(tempmask) = int8(-1);
    SET(no).MaR.Manual(:,:,tf,slice) = temp;
    mar('update');
  case 'MORubber'
    %Create the mask
    tempmask = segment('createmask',[SET(no).XSize SET(no).YSize],y,x) & SET(no).Scar.MyocardMask(:,:,slice);  
    temp = SET(no).Scar.NoReflow(:,:,slice);
    temp(tempmask) = false; %Remove the MO
    SET(no).Scar.NoReflow(:,:,slice) = temp;

    %Also remove from Manual
    temp = SET(no).Scar.Manual(:,:,slice);
    mask = (temp==2) & tempmask;
    temp(mask) = int8(1);
    SET(no).Scar.Manual(:,:,slice) = min(SET(no).Scar.Manual(:,:,slice),temp);

    if SET(no).Scar.UpdateDirectly
      viability('viabilitycalc');
    end

  case 'GeneralPen'
    currentobjectind = DATA.GeneralPenSettings.getcurrentobject;
    if isempty(currentobjectind) %object list is empty, create a new object
      currentobjectind = generalpen.generalpenfunctions('createnewobject');
    end    
    numpoints = SET(no).GeneralPenObjects(currentobjectind).getnumpoints;
    if isempty(SET(no).GeneralPenObjects(currentobjectind).X)
      SET(no).GeneralPenObjects(currentobjectind).Y = nan(numpoints,SET(no).TSize,SET(no).ZSize);
      SET(no).GeneralPenObjects(currentobjectind).X = nan(numpoints,SET(no).TSize,SET(no).ZSize);
    end
    [x,y] = tools('distributepointsalongcontour',x,y,numpoints);

    if ~DATA.ThisFrameOnly || doall
      %all frames mode
      if isrow(x)
        x = x';
      end
      if isrow(y)
        y = y';
      end
      SET(no).GeneralPenObjects(currentobjectind).Y(:,:,slice) = repmat(y,[1,SET(no).TSize]);
      SET(no).GeneralPenObjects(currentobjectind).X(:,:,slice) = repmat(x,[1,SET(no).TSize]);
    else
      SET(no).GeneralPenObjects(currentobjectind).Y(:,tf,slice) = y;
      SET(no).GeneralPenObjects(currentobjectind).X(:,tf,slice) = x;
    end

  case {'LA','RA'}
    if isempty(SET(no).(type)) %object is empty, create a new object
      generalpen.atriumpenfunctions('createnewobject',lower(type));
    end    
    numpoints = SET(no).(type).getnumpoints;
    if isempty(SET(no).(type).X)
      SET(no).(type).Y = nan(numpoints,SET(no).TSize,SET(no).ZSize);
      SET(no).(type).X = nan(numpoints,SET(no).TSize,SET(no).ZSize);
    end
    [x,y] = tools('distributepointsalongcontour',x,y,numpoints,opencontour);

    if ~DATA.ThisFrameOnly || doall
      %all frames mode
      if isrow(x)
        x = x';
      end
      if isrow(y)
        y = y';
      end
      SET(no).(type).Y(:,:,slice) = repmat(y,[1,SET(no).TSize]);
      SET(no).(type).X(:,:,slice) = repmat(x,[1,SET(no).TSize]);
    else
      SET(no).(type).Y(:,tf,slice) = y;
      SET(no).(type).X(:,tf,slice) = x;
    end
    %Update Strain flag
    if isequal(tf,SET(no).EDT) && ~isempty(SET(no).StrainMitt)
      updateredo(SET(no).StrainMitt,type);
    end
    if isequal(tf,SET(no).EST)
      SET(no).(type).isManualContourinES = true;
    end

  otherwise
    if isempty(SET(no).([type,'X']))
      datanumpoints = tools('getnumpointsforno',no);
      SET(no).([type,'Y']) = nan(datanumpoints,SET(no).TSize,SET(no).ZSize);
      SET(no).([type,'X']) = nan(datanumpoints,SET(no).TSize,SET(no).ZSize);
    end
    [x,y] = tools('distributepointsalongcontour',x,y,numpoints);

    if not(DATA.ThisFrameOnly) || doall  %findindented(DATA.Handles.hideiconholder,'allframesmode')
      %all frames mode
      if isrow(x)
        x = x';
      end
      if isrow(y)
        y = y';
      end
      SET(no).([type,'Y'])(:,:,slice) = repmat(y,[1,SET(no).TSize]);
      SET(no).([type,'X'])(:,:,slice) = repmat(x,[1,SET(no).TSize]);
    else
      SET(no).([type,'Y'])(:,tf,slice) = y;
      SET(no).([type,'X'])(:,tf,slice) = x;
    end

    %If straintagging initiated adjust LVupdated
    if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
      SET(no).StrainTagging.LVupdated = 1;
    end
    if isequal(tf,SET(no).EDT) && ~isempty(SET(no).StrainMitt)
      updateredo(SET(no).StrainMitt,type);
    end

    % redo scar and MaR caclcualtions since LV was changed
    calcfunctions('updatemarandscar',no);
%     calcfunctions('calclvvolume',no);
end

%draw changes
drawfunctions('drawno',no,isonlyannotation)

%update result panel
switch type
  case {'Endo','Epi'}
    if isempty(DATA.LVNO)
      DATA.LVNO = no;
    end
    segment('updatevolume');
    segment('updatemeasurement')
  case {'RVEndo','RVEpi'}
    if isempty(DATA.RVNO)
      DATA.RVNO = no;
    end
    segment('updatevolume');
  case 'Roi'
    segment('updateflow');
    DATA.DoThisFrameOnly = false; %flag for calc functions
  otherwise
    segment('updatemeasurement');
end

DATA.Handles.cursor.XData = nan;
DATA.Handles.cursor.YData = nan;
DATA.CursorX = [];
DATA.CursorY = [];
DATA.Handles.cursor.LineStyle = '-';

DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
DATA.fig.WindowButtonUpFcn = '';

%----------------------------------------------
function [x,y,newroi,roiind] = appendcontour(no,type,x,y,slice,tf,panel)
%----------------------------------------------
%two possibilities, either the segmentation overwrites the entire contour
%field in the SET struct or it is appended to it. this is determined by
%checking if it is closer to the startpoint of the drawing or the closest
%point on the contour. As the orientation of the contour is adjusted so
%both have a counterclockwise orientation we only need to handle two cases
%if we are going to append.
%These cases are if the start index of the prior is removed in the append
%or not. We set up both these cases and evaluate the two potential closed
%curves and then choose the one which constitutes the largest area.
global DATA SET

newroi = true;
roiind = length(SET(no).Roi)+1;
isopencontour = false;
numpointsalongcontour = 0;
switch type
  case {'Scar','MaR','MO','ScarRubber','MaRRubber','MORubber','CalciumPen', 'CalciumPenRemove', 'CalciumPenBlue', 'CalciumPenRed', 'CalciumPenLilac'}
    %This will trigger new contour case
    x_close = [];
    y_close = [];
    startdist = inf;
    enddist = inf;
    drawdist = 0;

  case 'Roi'
    if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
      [startdist,roiind1,ind1] = findfunctions('closestroi',panel,x(1),y(1),[],[],no);
      [enddist,roiind2,ind2] = findfunctions('closestroi',panel,x(end),y(end),[],[],no);
    else
      [startdist,roiind1,ind1] = findfunctions('closestroi',panel,x(1),y(1));
      [enddist,roiind2,ind2] = findfunctions('closestroi',panel,x(end),y(end));
    end
    drawdist = sqrt((x(1)-x(end)).^2+(y(1)-y(end)).^2);

    [~,closest] = min([startdist,enddist,drawdist]);
    switch closest
      case 1 %the starting point has a roi that is very close
        x_close = SET(no).Roi(roiind1).X(:,tf);
        y_close = SET(no).Roi(roiind1).Y(:,tf);
        newroi = false;
        roiind = roiind1;
      case 2 %the end point has a roi that is very close
        x_close = SET(no).Roi(roiind2).X(:,tf);
        y_close = SET(no).Roi(roiind2).Y(:,tf);
        newroi = false;
        roiind = roiind2;
      case 3 %the starting point is closest to the endpoint of the roi.
        x_close = [];
        y_close = [];
        newroi = true;
        emptyroiinds = find(cellfun(@isempty,{SET(no).Roi.X}));
        if ~isempty(emptyroiinds)
          roiind = emptyroiinds(1);
        else
          roiind = length(SET(no).Roi)+1;
        end
    end

  otherwise
    %Endo, Epi, RVEndo, RVEpi, LA, RA, GeneralPen
    newcontour = true;
    if strcmp(type,'GeneralPen')
      %GeneralPen
      objind = DATA.GeneralPenSettings.getcurrentobject;
      if ~isempty(objind) && ...
          ~isempty(SET(no).GeneralPenObjects(objind).X) && ...
          ~all(isnan(SET(no).GeneralPenObjects(objind).X(:,tf,slice)))
        newcontour = false;
        x_close = SET(no).GeneralPenObjects(objind).X(:,tf,slice);
        y_close = SET(no).GeneralPenObjects(objind).Y(:,tf,slice);
      end
    elseif any(matches({'LA','RA'},type,IgnoreCase=true))
      %LA, RA
      if ~isempty(SET(no).(type)) && ...
          ~isempty(SET(no).(type).X) && ...
          ~all(isnan(SET(no).(type).X(:,tf,slice)))
        newcontour = false;
        isopencontour = true;
        numpointsalongcontour = SET(no).(type).NumPoints;
        x_close = SET(no).(type).X(:,tf,slice);
        y_close = SET(no).(type).Y(:,tf,slice);
      end
    else
      %Endo, Epi, RVEndo, RVEpi
      if ~isempty(SET(no).([type,'X'])) && ...
          ~all(isnan(SET(no).([type,'X'])(:,tf,slice)))
        newcontour = false;
        x_close = SET(no).([type,'X'])(:,tf,slice);
        y_close = SET(no).([type,'Y'])(:,tf,slice);
      end
    end

    if ~newcontour
      %distance to and index of closestcontour point from drawing startpoint
      [startdist,ind1] = min((y_close - y(1)).^2 + (x_close-x(1)).^2);

      %distance to and index of closestcontour point from drawing endpoint
      [enddist,ind2] = min((y_close - y(end)).^2 + (x_close-x(end)).^2);

      %distance between drawing start and endpoint
      drawdist = (y(1)-y(end)).^2 + (x(1)-x(end)).^2;
    else
      %This will trigger new contour case
      x_close = [];
      y_close = [];
      startdist = inf;
      enddist = inf;
      drawdist = 0;
    end
end

if isrow(x_close)
  x_close = x_close';
  y_close = y_close';
end

%this is the new contour decision
if drawdist<startdist && drawdist<enddist
  y1 = [flipud(y);y(end)];

  x1 = [flipud(x);x(end)];

  y2 = [y;y(1)];

  x2 = [x;x(1)];

  %The following handles how concatenation should be done
  %depending on the starting index of the contour that should be
  %merged with the new contour
elseif ind1 < ind2 || (ind1 == ind2 && ind1 == numpointsalongcontour) % second condition to ensure appending at the end of open contour
  if isopencontour && ind2 == numpointsalongcontour
    % this case appends to the end of the open contour
    y1 = cat(1,y_close(ind1:ind2),flipud(y));

    x1 = cat(1,x_close(ind1:ind2),flipud(x));

    y2 = cat(1,y_close(1:ind1),y);

    x2 = cat(1,x_close(1:ind1),x);
  else
    y1 = cat(1,y(1),y_close(ind1:ind2),flipud(y));

    x1 = cat(1,x(1),x_close(ind1:ind2),flipud(x));

    y2 = cat(1,y_close(1:ind1),y,y_close(ind2:end));

    x2 = cat(1,x_close(1:ind1),x,x_close(ind2:end));
  end
  
else
  if isopencontour && ind2 == 1
    % this case append to the beginnning of the open contour
    y1 = cat(1,y_close(ind2:ind1),y);

    x1 = cat(1,x_close(ind2:ind1),x);

    y2 = cat(1,flipud(y),y_close(ind1:end));

    x2 = cat(1,flipud(x),x_close(ind1:end));
  else
    y1 = cat(1,y(end),y_close(ind2:ind1),y);

    x1 = cat(1,x(end),x_close(ind2:ind1),x);

    y2 = cat(1,y_close(1:ind2),flipud(y),y_close(ind1:end));

    x2 = cat(1,x_close(1:ind2),flipud(x),x_close(ind1:end));
  end
 
end

%Comparison of area choose maximum area contour
a1 = calcfunctions('stablepolyarea',x1,y1);
a2 = calcfunctions('stablepolyarea',x2,y2);
[~,i] = max([a1,a2]);

if i==1
  x = x1;
  y = y1;
else
  x = x2;
  y = y2;
end