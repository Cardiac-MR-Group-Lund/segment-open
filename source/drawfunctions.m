function [varargout] = drawfunctions(varargin)
% Functions for drawing in panels
% Klas

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-----------------------
function drawthumbnails(calculatepreview,sliderupdated) %#ok<DEFNU>
%-----------------------
%Draw all thumbnails. Calculatepreview is a boolean
%indicating if thumbnails needs to be redrawn.

global DATA SET
persistent setlength

if DATA.Silent || (~DATA.DataLoaded) || isempty(DATA.Handles.datasetaxes)
  return;
end

if isempty(DATA.VisibleThumbnails)
  DATA.VisibleThumbnails=1:min(DATA.Pref.NumberVisibleThumbnails,length(SET));
  setlength=length(SET);
end

thumbsize=DATA.GUISettings.ThumbnailSize;
if nargin<2 ||(nargin==2 && not(sliderupdated))
  segment('updateslider');
end

if nargin<1
  calculatepreview=0;
end

if setlength~=length(SET)
  calculatepreview=1;
end
setlength=length(SET);

if calculatepreview
  delete(DATA.Handles.datasetaxes.Children)
  calcfunctions('calcdatasetpreview');
end

DATA.Handles.datasetpreviewimage = image(DATA.DATASETPREVIEW,'parent',DATA.Handles.datasetaxes);
set(DATA.Handles.datasetpreviewimage,'ButtondownFcn',...
  sprintf('%s(''thumbnail_Buttondown'')','segment'));
%colormap(DATA.Colormap); %no longer needed, now truecolor /JU
axis(DATA.Handles.datasetaxes,'image','off','ij');
% axis(DATA.Handles.datasetaxes,)
ylim=[(DATA.VisibleThumbnails(1)-1) DATA.VisibleThumbnails(end)]*thumbsize+1;
set(DATA.Handles.datasetaxes,'ylim',ylim ,'dataaspectratio', [1 1 1]);
hold(DATA.Handles.datasetaxes,'on');

%draw frames
drawthumbnailframes;
%print number of thumbnail in the upper left corner of the image
DATA.printthumbnailnumber(thumbsize);

hold(DATA.Handles.datasetaxes,'off');

%---------------------------
function drawthumbnailframes
%---------------------------
%Draw frames around thumbnail images
global DATA SET NO

if isempty(DATA.Handles.datasetaxes)
  return
end

thumbsize=DATA.GUISettings.ThumbnailSize;
try
  delete(DATA.Handles.datasetpreviewline);
  delete(DATA.Handles.datasetflowline);
catch %#ok<CTCH>
end

hold(DATA.Handles.datasetaxes,'on');
%Make frames used for linked images
%DATA.Handles.datasetselectline = [];
DATA.Handles.datasetflowline = [];
visibleandlinked = unique([SET(DATA.VisibleThumbnails).Linked]);
for loop =  visibleandlinked%1:length(SET)
  ypos = (loop-1)*thumbsize+1;
  DATA.Handles.datasetflowline =  ...
    [DATA.Handles.datasetflowline ...
    plot(DATA.Handles.datasetaxes,...
    [1    1                thumbsize        thumbsize 1   ],...
    [ypos ypos+thumbsize-1 ypos+thumbsize-1 ypos      ypos],...
    'color',DATA.GUISettings.ThumbFlowLineColor)...
    ];
end
set(DATA.Handles.datasetflowline,'visible','off');

%Show frames around linked images
if length(SET(NO).Linked) > 1
  linkies = SET(NO).Linked(SET(NO).Linked ~= NO);
  ind = ismember(visibleandlinked,linkies);
  set(DATA.Handles.datasetflowline(ind),'visible','on');
end

%draw frame around current image
ypos = (NO-1)*thumbsize+1;
DATA.Handles.datasetpreviewline =  plot(DATA.Handles.datasetaxes,...
  [1    1                thumbsize        thumbsize 1   ],...
  [ypos ypos+thumbsize-1 ypos+thumbsize-1 ypos      ypos],...
  'color',DATA.GUISettings.ThumbLineColor);

%---------------------
function killhandles %#ok<DEFNU>
%-----------------------
global DATA

if ~isfield(DATA.Handles,'imageaxes')
  return
end

graphics = [];
for i = 1:length(DATA.Handles.imageaxes)
  graphics = [graphics;DATA.Handles.imageaxes(i).Children];
end

delete(graphics)
delete(DATA.Handles.boxaxes)
delete(DATA.Handles.imageaxes)

%---------------------
function setxynan(panels)%#ok<DEFNU>
%-----------------------
%this function hides all overlay and axes by nan setting all objects.
global DATA

%all imageaxes needs to be checked.
if nargin == 0
  panels = 1:length(DATA.Handles.imageaxes);
end

%If a panel contains any image information then it needs to be cleared as
%it is likely that there should be a new stack in that panel after the
%configuration of the panels.
graphics = [];
for i = panels
  if ~isempty(DATA.Handles.imagehandles(i).CData)
    graphics = [graphics;DATA.Handles.imageaxes(i).Children]; %#ok<AGROW>
  end
end

%nan set all objects in panel image is always last so fish by setting it as
%empty
for i = 1:length(graphics)
  if strcmp(graphics(i).Type,'contour')
    %This is the contour case
    %tmp = [1,0;0,0];
    %set(graphics(i),'ZData',tmp,'XData',tmp,'YData',tmp);
    %graphics(i).Visible = 'off'; %This is needed in 3DPrint due to Matlab bug
    
    graphics(i).ZData = [1,0;0,0];
    graphics(i).Visible = 'off'; %This is needed in 3DPrint due to Matlab bug
        
  elseif strcmp(graphics(i).Type,'image')
    graphics(i).CData = [];
  elseif isprop(graphics(i),'XData')
    graphics(i).XData = nan;
    graphics(i).YData = nan;
  elseif isprop(graphics(i),'Position')
    graphics(i).Position = [nan nan];
  end
end
        
%---------------------
function drawtext(panel) %#ok<DEFNU>
%-----------------------
%drawtext

global DATA SET

if isempty(panel)
  return
end

no = DATA.ViewPanels(panel);
slice = [];
switch DATA.ViewPanelsType{panel}
  case 'trans3DP'
    slice = SET(no).LevelSet.View.RSlice;
  case 'sag3DP'
    slice = SET(no).LevelSet.View.GSlice;
  case 'cor3DP'
    slice = SET(no).LevelSet.View.BSlice;
  case {'one','montage','montagerow'}
    slice = SET(no).CurrentSlice;
end

%we do not need to do the full update here it is sufficient if we update
%the current image and corresponding lines in the others

stri = '';
if ~isempty(SET(DATA.ViewPanels(panel)).DICOMImageType)
  stri = sprintf('%s\n',SET(DATA.ViewPanels(panel)).DICOMImageType);
end

if ~isempty(SET(DATA.ViewPanels(panel)).SeriesDescription)
  stri = [stri sprintf('%s\n',SET(DATA.ViewPanels(panel)).SeriesDescription)];
end

if ~isnan(SET(DATA.ViewPanels(panel)).TimeVector(SET(DATA.ViewPanels(panel)).CurrentTimeFrame))
  stri = [stri dprintf('Time: %.0f ms\n',1000*SET(DATA.ViewPanels(panel)).TimeVector(SET(DATA.ViewPanels(panel)).CurrentTimeFrame))];
end

if (SET(DATA.ViewPanels(panel)).ZSize>1) && ~isempty(slice)
  stri = [stri dprintf('Slice: %d',slice)];
end

DATA.Handles.text(panel).String = stri;

%---------------------
function drawmeasures(panel) %#ok<DEFNU>
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);

scale = viewfunctions('getscale',panel);
slices = viewfunctions('slicesinpanel',panel);
if not(isempty(SET(no).Measure))
  [measure,slice] = viewfunctions('getmeasurecoords',panel);
  x = nan;
  y = nan;
  xo = nan; %used to set other points not in this slice (measurementoutsideplane)
  yo = nan;
  
  %This makes measurement text rendering compact and index will depend on
  %number of measurements in the current timeframe and slice
  if DATA.Pref.BackgroundColor
    bgcolor = 'k';
  else
    bgcolor = 'none';
  end
  textcounter = 1;
  for loop=1:length(SET(no).Measure)
    ziv = round(measure(loop).Z);
    ziv = min(ziv):max(ziv);
    if any(ismember(slices,ziv)) && ...
        (SET(no).Measure(loop).T==SET(no).CurrentTimeFrame)||(isnan(SET(no).Measure(loop).T))
      
      %sometimes the column format is not used
      if ~iscolumn(measure(loop).Y)
        measure(loop).Y = measure(loop).Y';
        measure(loop).X = measure(loop).X';
      end
      
      if ~all(ismember(ziv , slices))
        DATA.Handles.measurement(panel).LineStyle = '--';
      else
        DATA.Handles.measurement(panel).LineStyle = '-';        
      end
      
      %if the current view is montage we need to translate the
      %contour to the correct slice if the entire
      if any(strcmp(DATA.ViewPanelsType{panel},{'montage','montagesegmented','montagerow'}))
        for i = 1:length(ziv)
          [x1,y1] = ind2sub(DATA.ViewPanelsMatrix{panel},find(ziv(i)==slices,1));
          x = [x;nan;scale*(measure(loop).Y+(x1-1)*SET(no).YSize)];
          y = [y;nan;scale*(measure(loop).X+(y1-1)*SET(no).XSize)];
        end
      else
        x = [x;nan;scale*(measure(loop).Y)];
        y = [y;nan;scale*(measure(loop).X)];
        curslice = SET(DATA.ViewPanels(panel)).CurrentSlice; %get current slice
        logind = (round(measure(loop).Z) ~= curslice); %logical index to not in this slice
        if sum(logind)>0
          xo = [xo;nan;scale*(measure(loop).Y(logind))];
          yo = [yo;nan;scale*(measure(loop).X(logind))];          
        end
        %This is used to place montage text correctly
        x1 = 1;
        y1 = 1;
      end
      
      if ~DATA.Run
        [ymax,ix] = max(measure(loop).Y);
        
        DATA.Handles.measurementtext(panel,textcounter).Position = scale*[...
          ymax+1+(x1-1)*SET(no).YSize...
          measure(loop).X(ix)+(y1-1)*SET(no).XSize];
        DATA.Handles.measurementtext(panel,textcounter).String = sprintf('%s\n%0.1f [mm]',SET(no).Measure(loop).Name,SET(no).Measure(loop).Length);
        DATA.Handles.measurementtext(panel,textcounter).BackgroundColor = bgcolor;
        textcounter = textcounter + 1;
      end
    end
  end
  
  if ~DATA.Run
    set(DATA.Handles.measurementtext(panel,textcounter:end),...
      'Position',[nan nan]);
  end
  DATA.Handles.measurement(panel).XData = x;
  DATA.Handles.measurement(panel).YData = y;
  DATA.Handles.measurementoutsideplane(panel).XData = xo;
  DATA.Handles.measurementoutsideplane(panel).YData = yo;
  
end


%---------------------
function draworthoanglehandle(panel) %#ok<DEFNU>
%-----------------------
%This function adds the orthoangle handle to the ortho panel
global DATA SET

scale = viewfunctions('getscale',panel);
no = DATA.ViewPanels(panel);

if isnan(DATA.Handles.orthoanglehandle.XData)
  %default is quarter from the image border
  [x,y] = calcfunctions('calcplaneintersections',no,no,'orth','gla');
  
  %this is the intersection point of all planeintersections
  xc = SET(no).HLA.slice;
  yc = SET(no).VLA.slice;
  
  %then we find in which direction on the gla line there is the most
  %space and place the anglehandle point here.
  val1 = norm([xc,yc]-[x(1),y(1)]);
  val2 = norm([xc,yc]-[x(2),y(2)]);
  
  if val1<val2
    x = xc + (xc-x(1))/2;
    y = yc + (yc-y(1))/2;
  else
    x = xc + (xc-x(end))/2;
    y = yc + (yc-y(end))/2;
  end
  
  DATA.Handles.orthoanglehandle(panel).Parent = DATA.Handles.imageaxes(panel);
  DATA.Handles.orthoanglehandle(panel).YData = scale*x;%[nan y(i([1 end]))*[0.3;0.7]];
  DATA.Handles.orthoanglehandle(panel).XData = scale*y;%[nan x(i([1 end]))*[0.3;0.7]];
end

%---------------------
function drawpanel(panel)
%-----------------------
%does what it's told from drawlist
global DATA

if (~DATA.Silent) && (~isempty(panel))
  for i = 1:length(DATA.drawlist{panel})
    fcn = DATA.drawlist{panel}{i};
    if isa(fcn,'function_handle')
      fcn();
    else
      eval(fcn);
    end
  end
end

%---------------------
function drawroi(panel,colortypes) %#ok<DEFNU>
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);

if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  originalno = no;% no actually shown in panel
  no = SET(no).Flow.MagnitudeNo;
  % resetting that current time frame is always the one that is shown in
  % the panel
  SET(no).CurrentTimeFrame = SET(originalno).CurrentTimeFrame;
end

scale = viewfunctions('getscale',panel);
slicestoinclude = viewfunctions('slicesinpanel',panel);

%Define colortypes. These are then used to concatenate into the correct
%index in the x,y variables which finally are dealt into the different
%color roi graphics objects.'
textreset = false;
if nargin == 1%colortypes
  colortypes = 'cgbrwkym';
  set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
  textreset = true;
end

x = cell(1,length(colortypes));
y = cell(1,length(colortypes));

%only work with rois that exist in the included slices
numroisinslice = nnz(cellfun(@(x,y) any(~isempty(x)) && any(x==slicestoinclude),{SET(no).Roi.Z}));
roistodo = find(cellfun(@(x,y,z) numroisinslice ~= 0 && any(~isempty(x)) && any(x==slicestoinclude) && ismember(y(1),colortypes) && any(~isempty(z)),{SET(no).Roi.Z},{SET(no).Roi.LineSpec},{SET(no).Roi.Area}));
if ~textreset
   if numroisinslice <= 10 
    if length(DATA.Handles.roitext(panel,:))> numroisinslice
      %reset all roitexts that are left after using ctrl-z
      set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
    else
      set(DATA.Handles.roitext(panel,1:numroisinslice),'Position',[nan nan]);
    end 
  else
    set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
  end
end
linewidth = DATA.Pref.LineWidth;
if isempty(linewidth)
  linewidth = 1;
end
if DATA.Pref.BackgroundColor
  bgcolor = 'k';
else
  bgcolor = 'none';
end
%KG: 
for loop = 1:length(roistodo)
  %KG: 
  [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slicestoinclude==SET(no).Roi(roistodo(loop)).Z));
  yl = (yl-1)*SET(no).XSize;
  xl = (xl-1)*SET(no).YSize;
  if roistodo(loop) == SET(no).RoiCurrent
    DATA.Handles.roicurrent(panel).XData = scale*(xl+SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
    DATA.Handles.roicurrent(panel).YData = scale*(yl+SET(no).Roi(roistodo(loop)).X(:,SET(no).CurrentTimeFrame));
    DATA.Handles.roicurrent(panel).Color = SET(no).Roi(roistodo(loop)).LineSpec(1);
    DATA.Handles.roicurrent(panel).LineWidth = linewidth + 1;
  else 
    cind = regexp(colortypes,SET(no).Roi(roistodo(loop)).LineSpec(1));
    x{cind} = cat(1,x{cind},nan,xl+SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
    y{cind} = cat(1,y{cind},nan,yl+SET(no).Roi(roistodo(loop)).X(:,SET(no).CurrentTimeFrame));
    DATA.Handles.([colortypes(cind),'roi'])(panel).LineWidth = linewidth;
  end
end

if ~DATA.Run
  %We only display text in one and orth viewmode
  if any(strcmp(DATA.ViewPanelsType{panel},{'one','orth'}))
    textloop = 1;
    for loop = 1:length(roistodo)
      %if ~all(isnan(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame)))
      %get sign for text      
      %KG: 
      if SET(no).Roi(roistodo(loop)).Sign > 0
        roisign = '';
      else
        roisign = ' (-)';
      end
% % %       cind = regexp(colortypes,SET(no).Roi((loop)).LineSpec(1));
      
      %ROI position
      
      % we also need to check if the roitext is within the panel limits
      % otherwise we dont plot the text
      %[ymin,ix] = min(SET(no).Roi((loop)).Y(:,SET(no).CurrentTimeFrame));
      if numroisinslice <= 10 
      %This checks if it is a flow roi we are dealing with (do we have the output generated from calcflow in all timeframes). If so
      %do not plot the Std and Mean intensities.      
      if ~isempty(SET(no).Flow)
        [ymin,ix] = min(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
        ypos = ymin-1;
        xpos = SET(no).Roi(roistodo(loop)).X(ix,SET(no).CurrentTimeFrame);
        halign = 'right';
        if isempty(SET(no).Roi(roistodo(loop)).Area)
          labelstr = sprintf('%s%s',SET(no).Roi(roistodo(loop)).Name,roisign);
        else
          labelstr = sprintf('%s%s\n%3.1f [cm^2]',...
                     SET(no).Roi(roistodo(loop)).Name,roisign, ...
                     SET(no).Roi(roistodo(loop)).Area(SET(no).CurrentTimeFrame));
        end
      else
        ypos = mean(SET(no).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame));
        xpos = mean(SET(no).Roi(roistodo(loop)).X(:,SET(no).CurrentTimeFrame));
        halign = 'left'; 
        labelstr = sprintf('%s%s\n%3.1f [cm^2]\n%3.1f +/- %3.1f',...
                  SET(no).Roi(roistodo(loop)).Name,roisign, ...
                  SET(no).Roi(roistodo(loop)).Area(SET(no).CurrentTimeFrame), ...
                  SET(no).Roi(roistodo(loop)).Mean(SET(no).CurrentTimeFrame), ...
                  SET(no).Roi(roistodo(loop)).StD(SET(no).CurrentTimeFrame));
      end
      DATA.Handles.roitext(panel,(textloop)).Position = scale*[ypos xpos];
      DATA.Handles.roitext(panel,(textloop)).HorizontalAlignment = halign;
      DATA.Handles.roitext(panel,(textloop)).String = labelstr;
      DATA.Handles.roitext(panel,(textloop)).BackgroundColor = bgcolor;
      
      textloop = textloop+1;
      end
    end
  else
    %set the different colortypes
    set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
  end
end

%set the different colortypes
for i = 1:length(colortypes)
  DATA.Handles.([colortypes(i),'roi'])(panel).XData = scale*x{i};
  DATA.Handles.([colortypes(i),'roi'])(panel).YData = scale*y{i};
end

%------------------------------
function drawplaneintersections %#ok<DEFNU>
%------------------------------
global DATA SET

panels = find(DATA.ViewPanels);
for p1 = panels %this is the nos panel which we are going to plot in
  if ~any(strcmp(DATA.ViewPanelsType{p1},{'montage','montagesegmented','montagerow'})) %dont plot plane intersections within montage views
    scale = viewfunctions('getscale',p1);    
    for p2 = panels
      if p1~= p2
        [x,y] = calcfunctions('calcplaneintersections',...
          DATA.ViewPanels(p1),DATA.ViewPanels(p2),DATA.ViewPanelsType{p1},...
          DATA.ViewPanelsType{p2});%,'one','one',DATA.slices(p1),DATA.slices(p2));

        DATA.Handles.planeintersection(p1,p2).YData = scale*x;
        DATA.Handles.planeintersection(p1,p2).XData = scale*y;
        if p2 == DATA.CurrentPanel
          DATA.Handles.planeintersection(p1,p2).Color = 'y';
        else
          DATA.Handles.planeintersection(p1,p2).Color = 'w';
        end
      end
    end
  else
    set(DATA.Handles.planeintersection(p1,:),'YData', nan, 'XData', nan);
  end
end

%----------------------
function drawspeedimage %#ok<DEFNU>
%----------------------
%Ensure speedimage is updated.
global DATA

drawimages(find(strcmp(DATA.ViewPanelsType,'speedim')));

%-----------------------
function im = speedimage %#ok<DEFNU>
%-----------------------
%generates the 3dp speed image in rgb format ready for display in segment
global DATA SET NO

%--- generate 2D 'speed' image

if isempty(SET(NO).LevelSet.SpeedIM)
  segment3dp.tools('updatespeed');
end

tempimage = getimage(SET(NO).LevelSet.Pen.Color);

%---Remap the image

%bug in fastremap and int16, then need to ensure it is in the range
if isa(tempimage,'int16')
  tempremapped = fastremap(min(max(tempimage,int16(SET(NO).minValue)),int16(SET(NO).maxValue)),SET(NO).LevelSet.Speed.IntensityMap,int16(SET(NO).minValue),int16(SET(NO).maxValue)); 
else
  tempremapped = fastremap(tempimage,SET(NO).LevelSet.Speed.IntensityMap);
  
end

tempremapped = single(tempremapped)/4000+0.5; %scale it 0...1
tempremapped = uint8(tempremapped*255)+1; %scale it to 1..256

%convert colormap to uint8
cmap = DATA.LevelSet.colormap;
cmap = uint8(255*cmap);
temprgb = cat(2,...
  cmap(tempremapped,1),...
  cmap(tempremapped,2),...
  cmap(tempremapped,3));

im = reshape(temprgb,[size(tempimage,1), size(tempimage,2), 3]);

%---------------------
function drawimages(panel) 
%-----------------------
global DATA SET

% if isempty(DATA.ViewIM{panel}) && any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','cor3DP','speedim'}))
%     im = createfunctions('createviewim',panel);
%     DATA.Handles.imagehandles(panel).CData = im;
%     return
if isempty(panel)
  return
end

if panel>length(DATA.ViewIM) || isempty(DATA.ViewIM{panel})
  createfunctions('createviewim',panel);
end

%if any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','cor3DP','speedim'}))
switch DATA.ViewPanelsType{panel}
  case 'trans3DP'
    im = drawfunctions('getoverlayimage','r');
    set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]);    
    return
  case 'sag3DP'
    im = drawfunctions('getoverlayimage','g');
    set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]); 
    return
  case 'cor3DP'
    im = drawfunctions('getoverlayimage','b');    
    set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]);
    return
  case 'speedim'
    im = drawfunctions('speedimage');
    set(DATA.Handles.imagehandles(panel),'CData',im,'XData',[1 size(im,2)],'YData',[1 size(im,1)]); 
    return
  case 'viewport'
    return
  case 'montagesegmented'
    im = createfunctions('createviewim',panel);
end

if ndims(DATA.ViewIM{panel})==5
  im=squeeze(DATA.ViewIM{panel}(:,:,SET(DATA.ViewPanels(panel)).CurrentTimeFrame,:,:));
  DATA.Handles.imagehandles(panel).CData = im;
elseif ~isempty(SET(DATA.ViewPanels(panel)).Colormap)
  DATA.Handles.imagehandles(panel).CData = ind2rgb(im,SET(DATA.ViewPanels(panel)).Colormap);
else
  im = DATA.ViewIM{panel}(:,:,SET(DATA.ViewPanels(panel)).CurrentTimeFrame);  
  DATA.Handles.imagehandles(panel).CData = cat(3,im,im,im);
end

%Fix for interpolated. More work is required to verify this.
updatexydata('imagehandles',panel,size(im));
  
%--------------------------------
function drawselectedframe(panel) %#ok<DEFNU>
%--------------------------------
%Function that highlights the selected panel

global DATA

%turn all frames off
for p = 1:length(DATA.ViewPanels)
  DATA.Handles.selectedframe(p).Visible = 'off';
end

xl = DATA.Handles.imageaxes(panel).XLim;
yl = DATA.Handles.imageaxes(panel).YLim;
DATA.Handles.selectedframe(panel).XData =[xl,fliplr(xl),xl(1)];
DATA.Handles.selectedframe(panel).YData =[yl(1),yl(1),yl(2),yl(2),yl(1)];
DATA.Handles.selectedframe(panel).Visible = 'on';

%------------------------------------------------------------------------
function drawselectedslice(panel)
%------------------------------------------------------------------------
%Function that highlights the selected slice
global DATA SET

no = DATA.ViewPanels(panel);

%if panels are in mode 'one'/'orth' or
%if user unselected all slices (start and end slices are empty),
%then no frame is shown
if any(strcmp(DATA.ViewPanelsType{panel},{'one','orth'}))||((isempty(SET(no).StartSlice))&&isempty(SET(no).EndSlice))
  DATA.Handles.selectedslice(panel).XData = nan;
  DATA.Handles.selectedslice(panel).YData = nan;
  return
end

slicestoinclude = viewfunctions('slicesinpanel',panel);

%create boxes
x = [];
y = [];

ss = find(SET(no).StartSlice==slicestoinclude);
se = find(SET(no).EndSlice==slicestoinclude);

for zloop=ss:se
  [x1,y1] = ind2sub(DATA.ViewPanelsMatrix{panel},zloop);
  x1 = (x1-1)*SET(no).YSize + 1.5;
  y1 = (y1-1)*SET(no).XSize + 1.5;
  x2 = -1.5+x1+SET(no).YSize;
  y2 = -1.5+y1+SET(no).XSize;
  x = [x nan x1 x2 x2 x1 x1];%nan separation creates separate box
  y = [y nan y1 y1 y2 y2 y1];
end

DATA.Handles.selectedslice(panel).XData = x;
DATA.Handles.selectedslice(panel).YData = y;

%fix so that linked images slices are plotted
nos = SET(no).Linked;
%nos(nos==no)=[]; %remove the current panel. EH: why?
panelstodo=find(ismember(DATA.ViewPanels,nos));
for loop=1:length(panelstodo)
  panel = panelstodo(loop);
  switch DATA.ViewPanelsType{panel}
    case {'montage','montagerow','montagefit','montagesegmented'}
      DATA.Handles.selectedslice(panel).XData = x;
      DATA.Handles.selectedslice(panel).YData = y;
      DATA.Handles.selectedslice(panel).LineStyle = '--';
      DATA.Handles.selectedslice(panel).Visible = 'on';
  end
end

%----------------------------
function drawviability(panel) %#ok<DEFNU>
%----------------------------
%Function to draw viability contour on screen used from drawimageslice,
%drawimagemontage
global DATA SET

hidecell={'hidescar'};
stateandicon = viewfunctions('iconson',hidecell);
state=[stateandicon{:,1}];
if state(1)
  return
end

type = DATA.ViewPanelsType{panel};
scale = viewfunctions('getscale',panel);

no = DATA.ViewPanels(panel);
%Get mask
switch type
  case {'montage','montagerow','montagefit','sax3'}
    result = viewfunctions('reshape2layout',SET(no).Scar.Result,panel);
  case {'one','orth'}
    result = SET(no).Scar.Result(:,:,SET(no).CurrentSlice);
  otherwise
    myfailed('Unknown viewtype in drawviabilityhelper.');
end

if sum(result(:))>0
  
  %Normal scar outline
  scarim = imresize(double(result),scale,'bilinear');
  DATA.Handles.scarcontour(panel).ZData = scarim;
  DATA.Handles.scarcontour(panel).Visible = 'on';
  %Fix for interpolated
  updatexydata('scarcontour',panel,size(scarim))
  %Weighted scar outline
  if SET(no).Scar.UseWeighting
    
    weighting = viability('viabilityweight',no); %This is a vector with weight and position as corresponding pixels in the whole volume.
    
    if ~isempty(weighting)
      %Get weighting
      temp = zeros(size(SET(no).Scar.Result));
      
      temp(SET(no).Scar.Result) = weighting;
      temp(SET(no).Scar.NoReflow) = 1; %Make sure MO is always included weighted graphically.
      
      %Sort the weighting
      sortedweighting = sort(weighting);
      weighted = sum(weighting(:)); %Calculate weighted percentage
      total = sum(SET(no).Scar.Result(:));
      f = weighted/total;
      
      %Find threshold so that it would include f% of the area
      f = 1-f; %Since we want the pixelse that are larger than the threshold
      thres = sortedweighting(min(max(round(length(sortedweighting)*f),1),length(sortedweighting)));
      
      %Convert to right size
      switch type
        case {'montage','montagerow','montagefit','sax3'}
          weighting = viewfunctions('reshape2layout',temp,panel);
        case {'one','orth'}
          weighting = temp(:,:,SET(no).CurrentSlice);
        otherwise
          myfailed('Unknown viewtype in drawviabilityhelper.');
      end
      if max(weighting(:)) < thres
        %no pink region to show
        DATA.Handles.weightedscarcontour(panel).ZData = [];
        DATA.Handles.weightedscarcontour(panel).XData = [];
        DATA.Handles.weightedscarcontour(panel).YData = [];
        DATA.Handles.weightedscarcontour(panel).Visible = 'off';
      else
        weightedscarim = imresize(weighting,scale,'bilinear');
        DATA.Handles.weightedscarcontour(panel).ZData = weightedscarim;
        DATA.Handles.weightedscarcontour(panel).Visible = 'on';
        DATA.Handles.weightedscarcontour(panel).LevelList = [thres thres];
        %Fix for interpolated
        updatexydata('weightedscarcontour',panel,size(weightedscarim));
      end      
    else
      DATA.Handles.weightedscarcontour(panel).ZData = [];
      DATA.Handles.weightedscarcontour(panel).XData = [];
      DATA.Handles.weightedscarcontour(panel).YData = [];
      DATA.Handles.weightedscarcontour(panel).Visible = 'off';
    end %nonempty weighting
  else
    DATA.Handles.weightedscarcontour(panel).ZData = [];
    DATA.Handles.weightedscarcontour(panel).XData = [];
    DATA.Handles.weightedscarcontour(panel).YData = [];
    DATA.Handles.weightedscarcontour(panel).Visible = 'off';
  end %use weighting
end

%Update taken mocontour and moextentcontour graphically
if sum(SET(no).Scar.NoReflow(:))>0
  
  if ~isfield(SET(no).Scar,'MOThreshold')
    SET(no).Scar.MOThreshold = 1.5;
  end
  
  %Extract and convert to right size
  %Take image and mask with no reflow.
  doit = true;
  
  switch type%DATA.ViewPanelsType{panel}
    case {'montage','montagerow','montagefit','sax3'}
      contourim = SET(no).Scar.IM;
      contourim(~SET(no).Scar.NoReflow) = 1; %Mask with no reflow. 1 is max.
      contourim = viewfunctions('reshape2layout',contourim,panel,NaN); %NaN is outsideelement.
    case {'one','orth'}
      if existfunctions('anyall',SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice))
        contourim = SET(no).Scar.IM(:,:,SET(no).CurrentSlice);
        contourim(~SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice)) = 1; %1 is max
      else
        doit = false;
      end
    otherwise
      myfailed('Unknown viewtype in drawviabilityhelper.');
  end
  
  if (not(isequal(length(SET(no).Scar.mthreshold),SET(no).ZSize)||isequal(length(SET(no).Scar.mthreshold),1)))&& doit
    %Length is not correct.
    viability('viabilitycalc'); %this will updat the size correctly.
  end
  
  if doit
    %If we need to do, then doit.
    
    %--- Find suitable threshold
    
    if isequal(SET(no).Scar.mthreshold,0)
      %Manual mode
      mthreshold = 1;
      logslice = SET(no).Scar.Result(:,:,SET(no).CurrentSlice);
      logslice = logslice & (~SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice));
      if sum(logslice(:))>0
        im = SET(no).Scar.IM(:,:,SET(no).CurrentSlice);
        mthreshold = min(im(logslice(:)));
        if mthreshold<0.01
          mthreshold=0.01;
        end
      end
    else
      %Some automated mode
      mthreshold = SET(no).Scar.mthreshold;
    end
    
    %Use remote intensity to calculate good threshold
    if length(mthreshold)<2
      thres = mthreshold;
    else
      thres = mthreshold(SET(no).CurrentSlice);
    end
    thres = thres*SET(no).Scar.MOThreshold; %OBS if you change here, you need ALSO to change in viabilitycalcvolume.
    mocoim = imresize(contourim,scale,'bilinear');
    DATA.Handles.mocontour(panel).ZData = mocoim;
    DATA.Handles.mocontour(panel).Visible = 'on';
    DATA.Handles.mocontour(panel).LevelList = [thres thres];
    DATA.Handles.moextentcontour(panel).ZData = imresize(contourim,scale,'bilinear');
    DATA.Handles.moextentcontour(panel).Visible = 'on';
    %Fix for interpolated
    szim = size(mocoim);
    updatexydata('mocontour',panel,szim);
    updatexydata('moextentcontour',panel,szim);    
  end
end

%----------------------
function drawmar(panel) %#ok<DEFNU>
%----------------------
%Function to draw MaR contour on screeen used from drawimageslice,
%drawimagemontage

global DATA SET
no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

%Get mask
switch DATA.ViewPanelsType{panel}
  case {'montage','montagerow','montagefit','sax3'}
    result = viewfunctions('reshape2layout',squeeze(SET(no).MaR.Result(:,:,SET(no).CurrentTimeFrame,:)),panel);
  case {'one','orth'}
    result = SET(no).MaR.Result(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
  otherwise
    myfailed('Unknown viewtype in mar-drawhelper.');
end
%Normal MaR outline
%Contour rendering not triggered for constant Zdata
if all(result(:) == 0)
  result(1)=1;
end
marimg = imresize(double(result),scale,'bilinear');
DATA.Handles.marcontour(panel).ZData = marimg;
updatexydata('marcontour',panel,size(marimg));
stateandicon = viewfunctions('iconson','hidemar');
if ~stateandicon{1}
  DATA.Handles.marcontour(panel).Visible = 'on';
end
DATA.Handles.marcontour(panel).LevelList = [0.5 0.5];

%------------------------------
function drawcentercross(panel) %#ok<DEFNU>
%------------------------------
global SET DATA
no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);
slices = viewfunctions('slicesinpanel',panel);

[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},slices-slices(1)+1);
DATA.Handles.centercross(panel).XData = scale*(SET(no).CenterY +(yl-1)*SET(no).YSize);
DATA.Handles.centercross(panel).YData = scale*(SET(no).CenterX +(xl-1)*SET(no).XSize);

%------------------------------
function drawno(no, isonlyannotation) %#ok<DEFNU>
%------------------------------
%This function is intended to be used after calls that alter a no and may
%effect multiple panels and linked images. This will be used in all
%callbacks altering contours and rois for now. If more speed is needed then we will
%need to write specific calls to draw functions at each callback.

global DATA SET NO

if DATA.Silent
  return
end

if nargin < 1
  no = NO;
end

if nargin < 2
  isonlyannotation = false;
end


%get all panels which dont have no in it or are linked these are the ones
%that possible have intersections with the updated no and therefore needs
%updating.
panelstodo = find(~ismember(DATA.ViewPanels,SET(no).Linked) & DATA.ViewPanels>0);
panelslinked = find(ismember(DATA.ViewPanels,SET(no).Linked));

%Then we udpate intersections in intersecting panels
for p = panelstodo
  drawintersections(p)
end

%update altered panels fully
for p = panelslinked
  if strcmp(DATA.ViewPanelsType{p},'one') || strcmp(DATA.ViewPanelsType{p},'orth')
    DATA.ViewIM{p} = [];
  end
  viewfunctions('updatedrawlist',p,isonlyannotation)
  drawpanel(p)
  if isonlyannotation
    % add "drawimage" again, since it was skipped
    viewfunctions('addtodrawlist',p,sprintf('drawfunctions(''drawimages'',%d)',p));
  end
  if isempty(DATA.ViewIM{p})
    createfunctions('createviewim',p);
  end
end

%if montage then we should update the selected slice.
if any(strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},{'montage','montagesegmented','montagerow'}))
  drawselectedslice(DATA.CurrentPanel)
end

%------------------------------
function drawpoint3D(panel) %#ok<DEFNU>
%------------------------------
global DATA NO
%this draws and updates the panel image intersections.
drawfunctions('drawplaneintersections')

if nargin == 1
  % add points to all panels
  scale = viewfunctions('getscale',panel);
  switch DATA.ViewPanelsType{panel}
    case 'hla'
      [y,z,x] = segment('getclickedcoords');
      DATA.Handles.point3D(panel).XData(end+1) = scale*y;
      DATA.Handles.point3D(panel).YData(end+1) = scale*z;
    case 'vla'
      [x,z,y] = segment('getclickedcoords');
      DATA.Handles.point3D(panel).XData(end+1) = scale*x;
      DATA.Handles.point3D(panel).YData(end+1) = scale*z;
    case 'gla'
      [ytemp,xtemp] = segment('getclickedcoords');
      [x,y,z] = calcfunctions('gla2sax',xtemp,ytemp,NO);
      [xo,yo] = calcfunctions('sax2gla',x,y,z,NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*xo;
      DATA.Handles.point3D(panel).YData(end+1) = scale*yo;
    case {'one','orth'}
      [y,x,z] = segment('getclickedcoords');
      DATA.Handles.point3D(panel).XData(end+1) = y;
      DATA.Handles.point3D(panel).YData(end+1) = x;
      x = x/scale;
      y = y/scale;
   case {'montage'}
      [y,x,z] = segment('getclickedcoords');
      [xofs,yofs] = calcfunctions('calcoffset',z,[],NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*(y+yofs);
      DATA.Handles.point3D(panel).YData(end+1) = scale*(x+xofs);
    case {'montagerow'}
      [y,x,z] = segment('getclickedcoords');
      [xofs,yofs] = calcfunctions('calcoffset',z,[],NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*(y+xofs);
      DATA.Handles.point3D(panel).YData(end+1) = scale*(x+yofs);
    otherwise
      [y,x,z] = segment('getclickedcoords');
      [xofs,yofs] = calcfunctions('calcoffset',z,[],NO);
      DATA.Handles.point3D(panel).XData(end+1) = scale*(y+yofs);
      DATA.Handles.point3D(panel).YData(end+1) = scale*(x+xofs);
  end 
  pos = calcfunctions('xyz2rlapfh',NO,x,y,z);
  for panel = find(DATA.ViewPanels ~= NO & DATA.ViewPanels > 0)
    no = DATA.ViewPanels(panel);
    scale = viewfunctions('getscale',panel);
    xyz = calcfunctions('rlapfh2xyz',no,pos(:,1),pos(:,2),pos(:,3));
    switch DATA.ViewPanelsType{panel}
      case {'one','orth'}
        newx = scale*xyz(2,:);
        newy = scale*xyz(1,:);
      case {'montage'}
        [xofs,yofs] = calcfunctions('calcoffset',xyz(3,:),[],no, panel);
        newx = scale*(xyz(2,:))+yofs;
        newy = scale*(xyz(1,:))+xofs;
      otherwise
        [xofs,yofs] = calcfunctions('calcoffset',xyz(3,:),[],no, panel);
        newx = scale*(xyz(2,:)+yofs);
        newy = scale*(xyz(1,:)+xofs);
    end
    
    DATA.Handles.point3D(panel).XData(end+1) = newx;%scale*xyz(2,:);
    DATA.Handles.point3D(panel).YData(end+1) = newy;%scale*xyz(1,:);
  end
else
  % delete all points be setting all values to nana
  for panel = find(DATA.ViewPanels > 0)
    DATA.Handles.point3D(panel).XData(:) = nan;
    DATA.Handles.point3D(panel).YData(:) = nan; 
  end
  
end

%------------------------------
function drawpoint(panel) %#ok<DEFNU>
%------------------------------
global DATA SET

no = DATA.ViewPanels(panel);

%clear all text.
if ~DATA.Run
  set(DATA.Handles.pointtext(panel,:),'Position', [nan nan]);
  textoffset = DATA.Handles.point(panel).MarkerSize/2;
end

if not(isempty(SET(no).Point)) 
  
  if DATA.Pref.BackgroundColor
    bgcolor = 'k';
  else
    bgcolor = 'none';
  end
  %check if placed in 3dp view then we need to draw x and y in a different
  %way
  
  if any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','speedim','cor3DP'}))
    
    if DATA.Handles.configiconholder.findindented('showpoint')
      [r,g,b] = segment3dp.tools('xyz2rgb',SET(no).Point.X,SET(no).Point.Y,SET(no).Point.Z);
      switch DATA.ViewPanelsType{panel}
        case 'trans3DP'
          [~,~,pointstodo] = findfunctions('closestpoint3dp',panel,1,1,SET(no).LevelSet.View.RSlice);
          x = g(pointstodo);
          y = b(pointstodo);
        case 'sag3DP'
          [~,~,pointstodo] = findfunctions('closestpoint3dp',panel,1,1,SET(no).LevelSet.View.GSlice);
          x = b(pointstodo);
          y = r(pointstodo);
        case 'cor3DP'
          [~,~,pointstodo] = findfunctions('closestpoint3dp',panel,1,1,SET(no).LevelSet.View.BSlice);
          x = g(pointstodo);
          y = r(pointstodo);
        otherwise
          pointstodo = [];
          x = nan;
          y = nan;
      end
    else
      pointstodo = [];
      x = [];
      y = [];      
    end
    
  else
    %Normal
    slices = viewfunctions('slicesinpanel',panel);
    scale = viewfunctions('getscale',panel);
    
    %this is the montageindex of the
    pointstodo = find(ismember(SET(no).Point.Z,slices)&...
      (ismember(SET(no).Point.T,SET(no).CurrentTimeFrame)|isnan(SET(no).Point.T))); %nan means constant over time
    
    [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},SET(no).Point.Z(pointstodo)-slices(1)+1);
    
    y = scale*(SET(no).Point.X(pointstodo)+(xl-1)*SET(no).XSize);
    x = scale*(SET(no).Point.Y(pointstodo)+(yl-1)*SET(no).YSize);
  end
  
  DATA.Handles.point(panel).XData = x;
  DATA.Handles.point(panel).YData = y;
  if ~DATA.Run
    for i = 1:length(x)
      DATA.Handles.pointtext(panel,i).Position = [x(i)+textoffset,y(i)];
      DATA.Handles.pointtext(panel,i).String = SET(no).Point.Label(pointstodo(i));
      DATA.Handles.pointtext(panel,i).BackgroundColor = bgcolor;
    end
  end
end

%---------------------
function drawinterp(panel,type) %#ok<DEFNU>
%-----------------------
%Draws interpolation points in panel.
global DATA SET


no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);
if nargin == 1
  type = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
else
  type = {type};
end

slicestoinclude = viewfunctions('slicesinpanel',panel);

%when playing we dont want to see interpolation points
if DATA.Run == 1
  slicestoinclude = [];
end

for i = 1:length(type)
  x = nan;
  y = nan;
  if SET(no).([type{i}(1:end-6),'InterpOngoing']) && ~strcmp(DATA.Handles.([lower(type{i}(1:end-6)),'contour'])(panel).LineStyle, 'none')
    updateinterpolationsettings(panel,type{i})
  end
  for j = 1:length(slicestoinclude)
    if ~isempty(SET(no).([type{i},'X'])) && ~isempty(SET(no).([type{i},'X']){SET(no).CurrentTimeFrame,slicestoinclude(j)})
      [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},j);
      yl = (yl-1)*SET(no).XSize;
      xl = (xl-1)*SET(no).YSize;
      x = cat(1,x,nan,xl+SET(no).([type{i},'Y']){SET(no).CurrentTimeFrame,slicestoinclude(j)});
      y = cat(1,y,nan,yl+SET(no).([type{i},'X']){SET(no).CurrentTimeFrame,slicestoinclude(j)});
    end
  end
  DATA.Handles.(lower(type{i}))(panel).XData = scale*x;
  DATA.Handles.(lower(type{i}))(panel).YData = scale*y;
end

%---------------------
function drawcontours(panel,type) %#ok<DEFNU>
%-----------------------
%This function draws contours in the designated panel of the input type. It
%depends on that the naming of the coordinate fields of the new contour is homogenous to the
%classical EndoX EndoY representation.
global DATA SET

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

%To reduce code amount we use dynamic calling
if nargin == 1
  type = {'Endo','Epi','RVEndo','RVEpi'};
else
  type = {type};
end

slicestoinclude = viewfunctions('slicesinpanel',panel);

for i = 1:length(type)
  x = [];
  y = [];
  if ~isempty(SET(DATA.ViewPanels(panel)).([type{i},'X']))
    for j = 1:length(slicestoinclude)
      [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},j);
      yl = (yl-1)*SET(no).XSize;
      xl = (xl-1)*SET(no).YSize;
      x = cat(1,x,nan,xl+SET(no).([type{i},'Y'])(:,SET(no).CurrentTimeFrame,slicestoinclude(j)));
      y = cat(1,y,nan,yl+SET(no).([type{i},'X'])(:,SET(no).CurrentTimeFrame,slicestoinclude(j)));
    end
    
  end
  %check if interpolationis ongoing in current NO/stack and adjust
  %contours appearance
  if  ~SET(no).([type{i},'InterpOngoing'])
    tools('resetinterpolationcontours',panel,type{i})
  else
    updateinterpolationsettings(panel,type{i})
  end
  DATA.Handles.([lower(type{i}),'contour'])(panel).XData = scale*x;
  DATA.Handles.([lower(type{i}),'contour'])(panel).YData = scale*y;
end

%---------------------
function drawintersections(panel)
%-----------------------
global DATA SET
persistent markersize
if nargin < 1
  panels = 1:length(DATA.ViewPanels);
else
  panels = panel;
end
  
for panel = panels
  no = DATA.ViewPanels(panel);
  if not(no==0)
    scale = viewfunctions('getscale',panel);
    t = SET(no).CurrentTimeFrame;
    [endointersectionx,endointersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'endo',t,DATA.ViewPanelsType{panel});
    %    calcfunctions('calcsegmentationintersections',DATA.ViewPanels(panel),'endo',t,DATA.ViewPanelsType{panel});
    [epiintersectionx,epiintersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'epi',t,DATA.ViewPanelsType{panel});
    %    calcfunctions('calcsegmentationintersections',DATA.ViewPanels(panel),'epi',t,DATA.ViewPanelsType{panel});
    
    %
    [rvendointersectionx,rvendointersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'rvendo',t,DATA.ViewPanelsType{panel});
    [rvepiintersectionx,rvepiintersectiony] = ...
      calcfunctions('calcsegmentationintersections',panel,'rvepi',t,DATA.ViewPanelsType{panel});
    if isempty(markersize) ||(markersize ~= DATA.Pref.MarkerSize)
      markersize = DATA.Pref.MarkerSize;
   
      set([DATA.Handles.endocontourintersection(panel),...
           DATA.Handles.epicontourintersection(panel),...
           DATA.Handles.rvendocontourintersection(panel),...
           DATA.Handles.rvepicontourintersection(panel)],'markersize',markersize);
    end
    
    DATA.Handles.endocontourintersection(panel).XData = scale*endointersectiony;
    DATA.Handles.endocontourintersection(panel).YData = scale*endointersectionx;
    DATA.Handles.epicontourintersection(panel).XData = scale*epiintersectiony;
    DATA.Handles.epicontourintersection(panel).YData = scale*epiintersectionx;
    DATA.Handles.rvendocontourintersection(panel).XData = scale*rvendointersectiony;
    DATA.Handles.rvendocontourintersection(panel).YData = scale*rvendointersectionx;
    DATA.Handles.rvepicontourintersection(panel).XData = scale*rvepiintersectiony;
    DATA.Handles.rvepicontourintersection(panel).YData = scale*rvepiintersectionx;
  end
end



% DATA.Handles.endocontourintersection(panel).XData = scale*DATA.endointersectiony{panel}{SET(DATA.ViewPanels(panel)).CurrentTimeFrame};
% DATA.Handles.endocontourintersection(panel).YData = scale*DATA.endointersectionx{panel}{SET(DATA.ViewPanels(panel)).CurrentTimeFrame};
% DATA.Handles.epicontourintersection(panel).XData = scale*DATA.epiintersectiony{panel}{SET(DATA.ViewPanels(panel)).CurrentTimeFrame};
% DATA.Handles.epicontourintersection(panel).YData = scale*DATA.epiintersectionx{panel}{SET(DATA.ViewPanels(panel)).CurrentTimeFrame};

%------------------------------
function showviabilityedits(panel) %#ok<DEFNU>
%------------------------------
%Show viability edits on screen as a temporary overlay.
global DATA SET

no = DATA.ViewPanels(panel);
hidecell={'hidescar'};
stateandicon = viewfunctions('iconson',hidecell);
state=[stateandicon{:,1}];

if isempty(SET(no).Scar)|| state(1)
  return;
  %viability('viabilityreset_Callback');
end

if isequal(get(DATA.Handles.viabilityshowinfarctaswhitemenu,'checked'),'on')
  showaswhite = true;
else
  showaswhite = false;
end

if ~isempty(SET(no).Scar.GreyZone.map) && ...
    isequal(get(DATA.Handles.viabilityshowgrayzonemenu,'checked'),'on')
  showgreyzone = true;
else
  showgreyzone = false;
end

% %If neither then we can just safely exit.
% if not(showaswhite || showgreyzone)
%     return;
% end

tempnos=no;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

panelstodo = find(DATA.ViewPanels==no);

for panel=panelstodo
  %Ok lets draw it
  switch DATA.ViewPanelsType{panel}
    case {'one','mmodespatial','orth'}
      
      temp = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      if showgreyzone
        temp = imresize(temp,size(SET(no).Scar.GreyZone.map(:,:,1)),'bilinear');
      end
      colmap = SET(no).Colormap;
      if isempty(colmap)
        colmap = colormap('gray');
      end
      tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp,no,colmap));
      
      sz=size(tempuint8);
      scarimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
      tmp = SET(no).Scar.Result(:,:,SET(no).CurrentSlice);
      tmp = logical(tmp(:));
      
      %Show infarct
      if DATA.Pref.LineWidth>0
        if showaswhite
          scarimrgb(tmp,:) = uint8(255);
          %Greyzone
        elseif showgreyzone
          greytmp = (SET(no).Scar.GreyZone.map(:,:,SET(no).CurrentSlice) == 1);
          coretmp = (SET(no).Scar.GreyZone.map(:,:,SET(no).CurrentSlice) == 2);
          scarimrgb(greytmp(:),:) = repmat(uint8([127 127 0]),sum(greytmp(:)),1);
          scarimrgb(coretmp(:),:) = repmat(uint8([127 0 0]),sum(coretmp(:)),1);
        end
      end
      
      if not(isfield(SET(no).Scar,'NoReflow'))
        SET(no).Scar.NoReflow = SET(no).Scar.Auto;
      end
      
      
      tmp=SET(no).Scar.Manual(:,:,SET(no).CurrentSlice);
      tmp=tmp(:);
      scarimrgb((tmp==int8( 1)),2) = uint8(255);
      scarimrgb((tmp==int8(-1)),3) = uint8(255);
      tmp=SET(no).Scar.NoReflow(:,:,SET(no).CurrentSlice);
      tmp=tmp(:);
      scarimrgb(logical(tmp),1) = uint8(255);
      
      scarimrgb=reshape(scarimrgb,[sz(1:2) 3]); %EH: added 3
      
    case {'montage','montagerow','montagefit','sax3'}
      % Convert to 2D and layout
      scale=1;
      colmap = SET(no).Colormap;
      if isempty(colmap)
        colmap = colormap('gray');
      end
      tempuint8 = calcfunctions('remapuint8',...
        viewfunctions('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:)),panel),...
        no,colmap);
      
      tempuint8 = min(uint8(255),uint8(1)+tempuint8);
      sz = size(tempuint8);
      scarimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
      
      tmp = viewfunctions('reshape2layout',SET(no).Scar.Result,panel);%used to be NO but this may cause error
      tmp = logical(tmp(:));
      
      %Show infarct
      if DATA.Pref.LineWidth>0
        if showaswhite
          scarimrgb(tmp,:) = uint8(255);
          %Greyzone
        elseif showgreyzone
          temp = imresize(SET(no).Scar.GreyZone.map,[SET(no).XSize SET(no).YSize],'bilinear');
          greytmp = viewfunctions('reshape2layout',temp,panel);
          gztmp = greytmp(:) == 1;
          coretmp = greytmp(:) == 2;
          
          scarimrgb(coretmp,:) = repmat(uint8([127 0 0]),sum(coretmp),1);
          scarimrgb(gztmp(:),:) = repmat(uint8([127 127 0]),sum(gztmp(:)),1);
        end
      end
      
      if not(isfield(SET(no).Scar,'NoReflow'))
        SET(no).Scar.NoReflow = repmat(uint8(1),size(SET(no).Scar.Auto));
      end
      
      
      tmp=viewfunctions('reshape2layout',SET(no).Scar.Manual,panel);
      tmp=tmp(:);
      scarimrgb((tmp==int8( 1)),2) = uint8(255);
      scarimrgb((tmp==int8(-1)),3) = uint8(255);
      tmp=viewfunctions('reshape2layout',SET(no).Scar.NoReflow,panel);
      tmp=tmp(:);
      scarimrgb(logical(tmp),1) = uint8(255);
      
      scarimrgb=reshape(scarimrgb,sz);
  end
  if not(DATA.Silent)
    scale = viewfunctions('getscale',panel);
    imxsz = size(scarimrgb,1)*scale;
    imysz = size(scarimrgb,2)*scale;
    im = imresize(scarimrgb,[imxsz imysz],'bilinear');
    set(DATA.Handles.imagehandles(panel),'CData',im);
  end
end

%--------------------------------
function draw3dpline(panel,type) %#ok<DEFNU>
%--------------------------------
%Draws stack intersections in 3dp.

global DATA SET
no = DATA.ViewPanels(panel);

if DATA.Handles.configiconholder.findindented('showcross')
  %Find max size
  %[rmax,gmax,bmax] = segment3dp.tools('xyz2rgb',SET(no).XSize,SET(no).YSize,SET(no).ZSize);
  %just use 8000 instead never used image that large
  
  switch type
    case 'br'
      x = [0, 8000];
      y = [SET(no).LevelSet.View.RSlice SET(no).LevelSet.View.RSlice];
    case 'gr'
      x = [0, 8000];
      y = [SET(no).LevelSet.View.RSlice SET(no).LevelSet.View.RSlice];
    case 'gb'
      x = [SET(no).LevelSet.View.BSlice SET(no).LevelSet.View.BSlice];
      y = [0, 8000];
    case 'rb'
      x = [0, 8000];
      y = [SET(no).LevelSet.View.BSlice SET(no).LevelSet.View.BSlice];
    case 'rg'
      x = [SET(no).LevelSet.View.GSlice SET(no).LevelSet.View.GSlice];
      y = [0, 8000];
    case 'bg'
      x = [SET(no).LevelSet.View.GSlice SET(no).LevelSet.View.GSlice];
      y = [0, 8000];
  end
else
  x = NaN;
  y = NaN;
end

DATA.Handles.([type,'line'])(panel).XData = x;
DATA.Handles.([type,'line'])(panel).YData = y;

%--------------------------------
function displaypoints(view,ima) %#ok<DEFNU>
%--------------------------------
%Draws annotation points in 3dp imageaxes ima according to view.

global SET DATA NO

if ~DATA.Handles.configiconholder.findindented('showpoint')
  %Hide
  delete(DATA.Handles.([view,'text3dp']))%This removes old points
else
  
  %add points
  [r,g,b] = segment3dp.tools('xyz2rgb',SET(NO).Point.X,SET(NO).Point.Y, SET(NO).Point.Z);
  switch view
    case 'r'
      pointstoplot = round(r)==SET(NO).LevelSet.View.RSlice;
      x = g;
      y = b;
    case 'g'
      pointstoplot = round(g)==SET(NO).LevelSet.View.GSlice;
      x = b;
      y = r;
    case 'b'
      pointstoplot = round(b)==SET(NO).LevelSet.View.BSlice;
      x = g;
      y = r;
  end
  
  delete(DATA.Handles.([view,'text3dp']))%This removes old points
  
  DATA.Handles.([view,'text3dp']) = [];
  
  if any(pointstoplot)
    for loop = find(pointstoplot)
      DATA.Handles.([view,'text3dp'])=[DATA.Handles.([view,'text3dp']),text('parent',ima,...
        'position',[x(loop) y(loop)],'HorizontalAlignment','center',...
        'string','+','color',[1,1,1]), text('parent',ima,...
        'position',[x(loop)+5 y(loop)],'HorizontalAlignment','left',...
        'string',SET(NO).Point.Label{loop},'color',[1,1,1])];
      set(DATA.Handles.([view,'text3dp']),'ButtonDownFcn','segment3dp.tools(''delete3dppoints'')');
    end
  end
end

%-----------------------------------
function draw3dpoutline(panel,type) %#ok<DEFNU>
%------------------------------------
%Shows contour of 3d segmentation if outline button is indented and motion function isnt running.
global DATA SET NO

% if speed im we want to use the current selected imagestack for contour.
%if strcmp(type,'s')
%  im = double(getimagehelper(SET(NO).LevelSet.BW,SET(NO).LevelSet.Pen.Color));
%else
%  im = double(getimagehelper(SET(NO).LevelSet.BW,type));
%end

switch DATA.ViewPanelsType{panel}
  case 'trans3DP'
    im = double(getimagehelper(SET(NO).LevelSet.BW,'r'));
  case 'sag3DP'
    im = double(getimagehelper(SET(NO).LevelSet.BW,'g'));
  case 'cor3DP'
    im = double(getimagehelper(SET(NO).LevelSet.BW,'b'));
  case 'speedim'
    im = double(getimagehelper(SET(NO).LevelSet.BW,SET(NO).LevelSet.Pen.Color));
end

%no segmentation
if sum(im(:)>127)==0 %all(im==0)
  DATA.Handles.threedpcontour(panel).ZData = [127 0 ; 0 0];
  return
end

DATA.Handles.threedpcontour(panel).ZData = im;
DATA.Handles.threedpcontour(panel).Visible = 'on'; 

%-----------------------------------------
function im = getimagehelper(vol,view)
%----------------------------------------
%Similar as getimage, but do assme no TSize as stored locally in BW.

global DATA SET NO

switch lower(DATA.LevelSet.imageorientation)
  case 'transversal'
    switch view
      case 'r'
        im = squeeze(vol(:,:,SET(NO).LevelSet.View.RSlice));
      case 'g'
        im = squeeze(vol(:,SET(NO).LevelSet.View.GSlice,:))';
      case 'b'
        im = squeeze(vol(SET(NO).LevelSet.View.BSlice,:,:))';
    end
  case 'sagittal'
    switch view
      case 'r'
        im = squeeze(vol(SET(NO).LevelSet.View.RSlice,:,:));
      case 'g'
        im = squeeze(vol(:,:,SET(NO).LevelSet.View.GSlice));
      case 'b'
        im = squeeze(vol(:,SET(NO).LevelSet.View.BSlice,:));
    end
  case 'coronal'
    switch view
      case 'r'
        im = squeeze(vol(SET(NO).LevelSet.View.RSlice,:,:))';
      case 'g'
        im = squeeze(vol(:,SET(NO).LevelSet.View.GSlice,:));
      case 'b'
        im = squeeze(vol(:,:,SET(NO).LevelSet.View.BSlice));
    end
end

%------------------------------
function im = getimage(view)
%-----------------------------
%Get image for view given imageorientation, r = transversal, g = sagittal,
%b=coronal

global DATA SET NO

switch lower(DATA.LevelSet.imageorientation)
  case 'transversal'
    switch view
      case 'r'
        im = squeeze(SET(NO).IM(:,:,1,SET(NO).LevelSet.View.RSlice));
      case 'g'
        im = squeeze(SET(NO).IM(:,SET(NO).LevelSet.View.GSlice,1,:))';
      case 'b'
        im = squeeze(SET(NO).IM(SET(NO).LevelSet.View.BSlice,:,1,:))';
    end
  case 'sagittal'
    switch view
      case 'r'
        im = squeeze(SET(NO).IM(SET(NO).LevelSet.View.RSlice,:,1,:));
      case 'g'
        im = squeeze(SET(NO).IM(:,:,1,SET(NO).LevelSet.View.GSlice));
      case 'b'
        im = squeeze(SET(NO).IM(:,SET(NO).LevelSet.View.BSlice,1,:));
    end
  case 'coronal'
    switch view
      case 'r'
        im = squeeze(SET(NO).IM(SET(NO).LevelSet.View.RSlice,:,1,:))';
      case 'g'
        im = squeeze(SET(NO).IM(:,SET(NO).LevelSet.View.GSlice,1,:));
      case 'b'
        im = squeeze(SET(NO).IM(:,:,1,SET(NO).LevelSet.View.BSlice));
    end
end

%----------------------------------
function im = getoverlayimage(view) %#ok<DEFNU>
%----------------------------------
%Get image including overlay, view is either of 'r','g','b'.

global DATA SET NO

im = getimage(view);
bw = getimagehelper(SET(NO).LevelSet.BW,view);
if ~isfield(DATA.LevelSet,'Man') || isempty(DATA.LevelSet.Man)
  DATA.LevelSet.Man = repmat(int8(0),size(SET(NO).LevelSet.BW));
end
man = getimagehelper(DATA.LevelSet.Man,view);

im = segment3dp.tools('levelsetremapandoverlay',im,bw,man);

%--------------------------------
function drawintensitymapping %#ok<DEFNU>
%--------------------------------
%Draws intensity mapping curve and add gray colorbar

global DATA SET NO

c = [0 0 0];

offset = SET(NO).LevelSet.IntensityOffset;
slope = SET(NO).LevelSet.IntensitySlope;
[minvalue,maxvalue] = segment3dp.tools('getminmax');
windowmin = SET(NO).LevelSet.WindowCenter-SET(NO).LevelSet.WindowWidth/2;
windowmax = SET(NO).LevelSet.WindowCenter+SET(NO).LevelSet.WindowWidth/2;

%Plot mapping--
DATA.Handles.intensityline = plot(DATA.Handles.intensityaxes,...
  [minvalue-abs(minvalue) -offset/slope (0.5-offset)/slope (1-offset)/slope maxvalue+abs(maxvalue)],...
  [0 0 0.5 1 1],...
  'b-');

set(DATA.Handles.intensityaxes,'ytick',[0 1],...
  'yticklabel',{'0%' ,'100%'});

if isequal(SET(NO).Modality,'CT') %&& ~isempty(regexpi(class(SET(NO).IM),'int'))
  mintext = sprintf('%d HU',round(windowmin));
  maxtext = sprintf('%d HU',round(windowmax));
  set(DATA.Handles.contrastlistbox,'visible','on');
else
  mintext = sprintf('%d',windowmin);
  maxtext = sprintf('%d',windowmax);
  set(DATA.Handles.contrastlistbox,'visible','off');
end

set(DATA.Handles.intensityaxes,...
  'xtick',[windowmin windowmax],...
  'xticklabel',{mintext maxtext});

%Add markers to move the lines
hold(DATA.Handles.intensityaxes,'on');

%x coordinates for the three points are:
x1 = max(windowmin,(0-offset)/slope);
x3 = min(windowmax,(1-offset)/slope);
x2 = (x1+x3)/2;

%use linear algebra to find y coordinates, y = kx+m
y1 = x1*slope+offset;
y2 = x2*slope+offset;
y3 = x3*slope+offset;

DATA.Handles.intensitypoint1 = plot(DATA.Handles.intensityaxes,windowmin,y1,'bo');
DATA.Handles.intensitypoint2 = plot(DATA.Handles.intensityaxes,(windowmin+windowmax)/2,y2,'bo');
DATA.Handles.intensitypoint3 = plot(DATA.Handles.intensityaxes,windowmax,y3,'bo');

set(DATA.Handles.intensitypoint1,'markersize',8,'buttondownfcn','segment3dp.tools(''intensitypoint1buttondown'')');
set(DATA.Handles.intensitypoint2,'markersize',8,'buttondownfcn','segment3dp.tools(''intensitypoint2buttondown'')');
set(DATA.Handles.intensitypoint3,'markersize',8,'buttondownfcn','segment3dp.tools(''intensitypoint3buttondown'')');

hold(DATA.Handles.intensityaxes,'off');

extra = (windowmax-windowmin)*0.1;
set(DATA.Handles.intensityaxes,'ylim',[-0.1 1.1],'xlim',[windowmin-extra windowmax+extra],...
  'xcolor',c,'ycolor',c,'zcolor',c);

segment3dp.tools('updateintensityline');

%--------------------------------------------------------
function [x1,y1,x2,y2,x3,y3] = getpointcoordinates(no,x,y)
%--------------------------------------------------------
%Get points on speedmapping

global SET

center = SET(no).LevelSet.Speed.Center;
width = SET(no).LevelSet.Speed.Width;
type = SET(no).LevelSet.Speed.MappingMode;

%Get x-coordinate for points
switch type
  case {'positiveslope','negativeslope'}
    x1 = center-width/2;
    x2 = center;
    x3 = center+width/2;
  case 'gaussian'
    [~,ind] = find(y>=0,1,'first');
    x1 = x(ind);
    x2 = center;
    [~,ind] = find(y>=0,1,'last');
    x3 = x(ind);
end

%Ensure visible on screen
if x1<(center-width/2)
  x1 = center-width/2;
end

if x3>(center+width/2)
  x3 = center+width/2;
end

%Find y-coordinates
y1 = interp1(x,y,x1);
y2 = interp1(x,y,x2);
y3 = interp1(x,y,x3);

%-------------------------
function updatemapping %#ok<DEFNU>
%-------------------------
%update the mapping graphically

global DATA SET NO

%y = k(x-offset)

[x,y] = segment3dp.tools('getmappingline',NO);

SET(NO).LevelSet.Speed.IntensityMap = int16(2000*y);

[x1,y1,x2,y2,x3,y3] = getpointcoordinates(NO,x,y);

set(DATA.Handles.speedimline,'xdata',x,'ydata',y);
set(DATA.Handles.speedpoint1,'xdata',x1,'ydata',y1);
set(DATA.Handles.speedpoint2,'xdata',x2,'ydata',y2);
set(DATA.Handles.speedpoint3,'xdata',x3,'ydata',y3);

%----------------------
function drawmapping %#ok<DEFNU>
%----------------------
%Prepare the mapping area for drawing

global DATA SET NO

fgc = [0 0 0];
bgc = [0.94 0.94 0.94];

windowmin = min(SET(NO).LevelSet.WindowCenter-SET(NO).LevelSet.WindowWidth/2,SET(NO).LevelSet.Speed.Center-SET(NO).LevelSet.Speed.Width/2);
windowmax = max(SET(NO).LevelSet.WindowCenter+SET(NO).LevelSet.WindowWidth/2,SET(NO).LevelSet.Speed.Center+SET(NO).LevelSet.Speed.Width/2);

%Plot mapping
[x,y] = segment3dp.tools('getmappingline',NO);
DATA.Handles.speedimline = plot(DATA.Handles.mappingaxes,x,y,'b-');

set(DATA.Handles.mappingaxes,'ytick',[-1 0 1],...
  'yticklabel',{'0%' '50%','100%'});

if isequal(SET(NO).Modality,'CT') %&& ~isempty(regexpi(class(SET(NO).IM),'int'))
  mintext = sprintf('%d HU',round(windowmin));
  maxtext = sprintf('%d HU',round(windowmax));
else
  mintext = sprintf('%d',windowmin);
  maxtext = sprintf('%d',windowmax);
end

set(DATA.Handles.mappingaxes,'xtick',[windowmin windowmax],...
  'xticklabel',{mintext maxtext});

%Plot center line
extra = (windowmax-windowmin)*0.1;
hold(DATA.Handles.mappingaxes,'on');
plot(DATA.Handles.mappingaxes,[windowmin-extra windowmax+extra],[0 0],'k:');

%Get coordinates of the points
[x1,y1,x2,y2,x3,y3] = getpointcoordinates(NO,x,y);

%Add markers to move the lines
DATA.Handles.speedpoint1 = plot(DATA.Handles.mappingaxes,x1,y1,'bo');
DATA.Handles.speedpoint2 = plot(DATA.Handles.mappingaxes,x2,y2,'bo');
DATA.Handles.speedpoint3 = plot(DATA.Handles.mappingaxes,x3,y3,'bo');
set(DATA.Handles.speedpoint1,'markersize',8,'buttondownfcn','segment3dp.tools(''speedpoint1buttondown'')');
set(DATA.Handles.speedpoint2,'markersize',8,'buttondownfcn','segment3dp.tools(''speedpoint2buttondown'')');
set(DATA.Handles.speedpoint3,'markersize',8,'buttondownfcn','segment3dp.tools(''speedpoint3buttondown'')');

hold(DATA.Handles.mappingaxes,'off');

%Adjust color and limits
set(DATA.Handles.mappingaxes,'xlim',[windowmin-extra windowmax+extra],'ylim',[-1.05 1.05],...
  'xcolor',fgc,...
  'ycolor',fgc,...
  'zcolor',fgc);

hv = [...
  DATA.Handles.decreaseintensitypushbutton ...
  DATA.Handles.increaseintensitypushbutton ...
  DATA.Handles.decreaseoffsetpushbutton ...
  DATA.Handles.centeredit ...
  DATA.Handles.increaseoffsetpushbutton ...
  DATA.Handles.slopetext ...
  DATA.Handles.widthedit ...
  DATA.Handles.positivesloperadiobutton ...
  DATA.Handles.negativesloperadiobutton ...
  DATA.Handles.gaussianradiobutton ...
  DATA.Handles.pickpushbutton ...
  DATA.Handles.predefinedlistbox ...
  DATA.Handles.windowtext ...
  DATA.Handles.centertext ...
  ];
set(hv,'Backgroundcolor',bgc,'Foregroundcolor',fgc);

%----------------------------------------
function drawintersectionpoints(no,panel)
%----------------------------------------
%Initiate handles and draw intersections with other contours
global DATA

if isempty(DATA.Handles.hideothercontouricon) || isequal(get(DATA.Handles.hideothercontouricon,'state'),'off')
  [endointersectline,maxintersect] = segment('getendointersection',no);
  viewtype = DATA.ViewPanelsType{panel};
  [endox,endoy] = calcfunctions('calcsegmentationintersections',no,'endo',viewtype);
  [epix,epiy] = calcfunctions('calcsegmentationintersections',no,'epi',viewtype);
  [rvendox,rvendoy] = calcfunctions('calcsegmentationintersections',no,'rvendo',viewtype);
  [rvepix,rvepiy] = calcfunctions('calcsegmentationintersections',no,'rvepi',viewtype);
  
  DATA.Handles.endointersectionline{panel} = zeros(1,maxintersect);
  for i=1:maxintersect
    if(i<=length(endointersectline))
      n = endointersectline(i).NPoints;
      DATA.Handles.endointersectionline{panel}(i) = plot(DATA.Handles.imageaxes(panel),...
        endointersectline(i).Y(1:n),endointersectline(i).X(1:n),'r');
    else
      % create empty handle for intersections in other time frames
      DATA.Handles.endointersectionline{panel}(i) = plot(DATA.Handles.imageaxes(panel),...
        [0 1],[0 1],'r');
      set(DATA.Handles.endointersectionline{panel}(i),'visible','off');
    end
  end
  DATA.Handles.endointersectionpoints{panel} = plot(DATA.Handles.imageaxes(panel),...
    endoy,endox,'r.');
  DATA.Handles.epiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    epiy,epix,'g.');
  DATA.Handles.rvendointersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    rvendoy,rvendox,'m.');
  DATA.Handles.rvepiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    rvepiy,rvepix,'c.');
else
  for i = 1:length(DATA.Handles.endointersectionline{panel})
    DATA.Handles.endointersectionline{panel}(i) = plot(DATA.Handles.imageaxes(panel),...
      NaN,NaN,'r');
  end
  DATA.Handles.endointersectionpoints{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'r.');
  DATA.Handles.epiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'g.');
  DATA.Handles.rvendointersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'m.');
  DATA.Handles.rvepiintersection{panel} = plot(DATA.Handles.imageaxes(panel),...
    NaN,NaN,'c.');
end

if ~DATA.Pref.BlackWhite
  set(DATA.Handles.endointersectionpoints{panel},'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.epiintersection{panel},'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvendointersection{panel},'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvepiintersection{panel},'markersize',DATA.Pref.MarkerSize);
else
  set(DATA.Handles.endointersectionline{panel},'color',[1 1 1]);
  set(DATA.Handles.endointersectionpoints{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
  set(DATA.Handles.epiintersection{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
  set(DATA.Handles.rvendointersection{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
  set(DATA.Handles.rvepiintersection{panel},'markersize',DATA.Pref.MarkerSize,'color',[1 1 1]);
end

%---------------------------------------
function updateintersectionpoints(panel)
%---------------------------------------
%Draw intersection with segmentation in other image stacks

global DATA
no = DATA.ViewPanels(panel(1));
viewtype = DATA.ViewPanelsType{panel};
if ~isequal(get(DATA.Handles.hideothercontouricon,'state'),'on')
  endointersect = segment('getendointersection',no);
  [endox,endoy] = calcfunctions('calcsegmentationintersections',no,'endo',viewtype);
  [epix,epiy] = calcfunctions('calcsegmentationintersections',no,'epi',viewtype);
  [rvendox,rvendoy] = calcfunctions('calcsegmentationintersections',no,'rvendo',viewtype);
  [rvepix,rvepiy] = calcfunctions('calcsegmentationintersections',no,'rvepi',viewtype);
  set(DATA.Handles.endointersectionline{panel},'visible','off');
  for i=1:length(endointersect)
    n = endointersect(i).NPoints;
    set(cellref(DATA.Handles.endointersectionline(panel),i),'xdata',endointersect(i).Y(1:n),...
      'ydata',endointersect(i).X(1:n),'visible','on');
  end
  if ~isempty([endox epix rvendox rvepix])
    set([DATA.Handles.endointersectionpoints{panel}], ...
      'XData',endoy,'YData',endox);
    set([DATA.Handles.epiintersection{panel}], ...
      'XData',epiy,'YData',epix);
    set([DATA.Handles.rvendointersection{panel}], ...
      'XData',rvendoy,'YData',rvendox);
    set([DATA.Handles.rvepiintersection{panel}], ...
      'XData',rvepiy,'YData',rvepix);
  end
end
%---------------------------------------
function updateinterpolationsettings(panel,type)
%---------------------------------------
global DATA
DATA.Handles.([lower(type),'contour'])(panel).LineStyle = 'none';
DATA.Handles.([lower(type),'contour'])(panel).Marker = '.';
DATA.Handles.([lower(type),'contour'])(panel).MarkerSize = 4; 

%---------------------------------------
function updatexydata(handlename,panel,szim)
%---------------------------------------
% this function updates XData and YData depending on the scale and image
% size
global DATA

if DATA.Pref.ViewInterpolated && strcmp(DATA.ViewPanelsType{panel},'one') 
    %Assume scale = 2
  DATA.Handles.(handlename)(panel).XData = 1.5:1:szim(2)+0.5; %XData is corner coordinates for image
  DATA.Handles.(handlename)(panel).YData = 1.5:1:szim(1)+0.5;
else
  DATA.Handles.(handlename)(panel).XData = 1:1:szim(2); %XData is corner coordinates for image
  DATA.Handles.(handlename)(panel).YData = 1:1:szim(1);
end
