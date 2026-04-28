function buttondownfunctions(varargin)
% Functions for buttondowns

% Broken out by Klas

%Remove all buttondownfunctions calls, unecessary overhead

%#ok<*GVMIS> 
%Invoke subfunction
feval(varargin{:}); % FEVAL switchyard

%-------------------------------------
function buttondown(panel,currenttool)
%-------------------------------------
global DATA SET

persistent ismaxnumrois

if DATA.issegment3dp
  segment3dp.buttondownfunctions('buttondown',panel,currenttool)
  return
end

scale = viewfunctions('getscale',panel);
no = DATA.ViewPanels(panel);
clicktype = get(DATA.fig,'SelectionType');

%store the current panel after graphical update and return. The
%procedure is that first the panel is selected then you can start selecting
%slices and rois.
if ~isequal(DATA.CurrentPanel,panel)
  viewfunctions('switchpanel',panel)
  DATA.CurrentPanel = panel;

  return
end

%if montage select slice
if any(strcmp(DATA.ViewPanelsType{panel},{'montagesegmented','montagerow','montage'}))
  [x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
  slice = viewfunctions('clickedslice',panel,x,y);
  if ~isempty(slice)
    if ~isequal(slice,SET(no).CurrentSlice) && ismember(currenttool,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp','GeneralPenInterp'})
      SET(no).StartSlice = slice;
      SET(no).EndSlice = slice;
      SET(no).CurrentSlice = slice;
      drawfunctions('drawno',no);
      return;
    else
      SET(no).StartSlice = slice;
      SET(no).EndSlice = slice;
      SET(no).CurrentSlice = slice;
      drawfunctions('drawno',no);
    end
  else
    return
  end
end

%Only measures available for these panels
if any(strcmp(DATA.ViewPanelsType{panel},{'hla','gla','vla'}))
  if ~strcmp(currenttool,'Measure')
    if ~strcmp(currenttool,'select')
      return
    end
  end
end

%This enables undo.
tools('enableundo',no);

switch clicktype
  case 'normal'
    %these are used to store information about rightclick. Type tells us
    %which context menu to use objectind tells us which index we right
    %clicked next to.
    DATA.LastObjectType = [];
    DATA.LastObject = [];
  case 'open'
    %this is the double click
    % check the current tool, since double click option is only valid for
    % select tool
    if ~strcmp(currenttool,'select') 
      return
    end
    if isempty(DATA.LastView)
      % save the last view
      DATA.LastView.ViewPanels = DATA.ViewPanels;
      DATA.LastView.ViewPanelsType = DATA.ViewPanelsType;
      DATA.LastView.ViewMatrix = DATA.ViewMatrix;
      DATA.LastView.CurrentPanel = DATA.CurrentPanel;
      DATA.LastView.NormalZoomState = {SET.NormalZoomState};
      if numel(DATA.ViewPanels) == 1
        % check whether to switch to one panel view
        switch DATA.ViewPanelsType{1}
          case 'one'
            newview = 'montage';
          case {'montage', 'montagerow'}
            newview = 'one';
          otherwise
            return
        end
        viewfunctions('setviewtype',newview);
      else
        % go to one view
        DATA.ViewMatrix = [1,1];
        viewfunctions('setview',DATA.ViewMatrix(1),DATA.ViewMatrix(2),[]);
      end
      return
    else
      %reset to previous view
      try
        DATA.ViewPanels = DATA.LastView.ViewPanels;
        DATA.ViewPanelsType = DATA.LastView.ViewPanelsType;
        DATA.ViewMatrix = DATA.LastView.ViewMatrix;
        DATA.CurrentPanel = DATA.LastView.CurrentPanel;
        panelstodo = find(DATA.LastView.ViewPanels);
        for actpanel = panelstodo
          no = DATA.ViewPanels(actpanel);
          SET(no).NormalZoomState = DATA.LastView.NormalZoomState{no};
        end
        viewfunctions('setview',DATA.ViewMatrix(1),DATA.ViewMatrix(2),DATA.ViewPanels);
      catch me
        mydispexception(me);
      end
      DATA.LastView = [];
    end

  case 'alt'
    [lvrvongoing,laraongoing,generalpenongoing] = helperfunctions('isanyinterpongoing',no);
    if lvrvongoing
      switch currenttool
        case {'EndoInterp','EpiInterp'}
          if isempty(DATA.LVNO)
            DATA.LVNO = no;
          end
        case {'RVEndoInterp','RVEpiInterp'}
          if isempty(DATA.RVNO)
            DATA.RVNO = no;
          end
      end
      tools('connectinterpolation',no,currenttool);
      return
    elseif laraongoing
      if isempty(SET(no).(currenttool(1:end-6))) %object is empty, create a new object
        generalpen.atriumpenfunctions('createnewobject',lower(currenttool(1:end-6)));
      end
      tools('connectlarainterpolation',no,currenttool);
      return
    elseif generalpenongoing
      tools('connectgeneralpeninterpolation',no);
      return
    end
    [y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
    slice = viewfunctions('clickedslice',panel,y,x);
    slices = viewfunctions('slicesinpanel',panel);

    %normalize clicked position to contour domain
    [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
    imdim = zoomfunctions.getxysize(no,panel);
    x = x/scale - ((xl-1)*imdim.XSize - imdim.XStart +1);
    y = y/scale - ((yl-1)*imdim.YSize - imdim.YStart +1);
    stateandicon = viewfunctions('iconson','point');
    pointstate = stateandicon{1};
    modifier = get(gcbf, 'Currentmodifier');
    if pointstate && not(isempty(modifier)) && strcmp(modifier,'control')
      % if place annotation point is on then place point where the right
      % click was done
      [~,pointind] = findfunctions('closestpoint',panel,x,y);

      %This way if you click close to a measure you drag that point of the
      %measure.
      DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''point_buttonup'',%d,%d,%d)',panel,pointind,slice);
      return
    else
      DATA.LastClickType = clicktype;% helper in findfunctions(closestmeasure)
      [type,objectind] = findfunctions('closestobject',panel,x,y);
      %place context menu. which one is triggered depends on the
      %closestobject
      DATA.LastObjectType = type;
      DATA.LastObject = [objectind,slice,SET(no).CurrentTimeFrame];
      viewfunctions('placecontextmenu',panel,type,objectind);
      DATA.LastClickType = [];
      return
    end

  case 'extend'
    if SET(no).EndoInterpOngoing || SET(no).EpiInterpOngoing ||...
        SET(no).RVEndoInterpOngoing || SET(no).RVEpiInterpOngoing
      tools('connectinterpolation',no,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})
      return
    end
    if ~strcmp(DATA.CurrentTool,'Measure')% && ~strcmp(DATA.CurrentTool,'spliteraser')
      if ~any(strcmp(DATA.ViewPanelsType{panel},{'montagesegmented','montagerow','montage'}))
        %Do panning of current image
        pan_buttondown(panel);
        return
      end
    end
end

%check if number of points, rois or measurements in current slice is larger than 20 then give
%message that this is not supported.

%Rois
maxnumroitext = viewfunctions('getnumberofmaxroitext',panel);
roitools = {'roi', 'roiput', 'roiballoon'};
if ismember(lower(currenttool),roitools)
  %check number of rois
  if strcmp(DATA.ViewPanelsType{panel},'one')
    % number of ROIs in the current slice
    numrois = sum(cellfun(@(x) any(~isempty(x)) && x==SET(no).CurrentSlice,{SET(no).Roi.Z}));
  else
    % number of overall ROIs
    numrois = length(SET(no).Roi);
  end
  if numrois >= maxnumroitext
    if isempty (ismaxnumrois)
      ismaxnumrois = true;
      msg1 = dprintf('Cannot have more than %d ROI labels in one slice view.',12);
      msg2 = dprintf('Cannot have more than %d ROI labels in montage view.',35);
      mymsgbox([msg1,newline,msg2],'Maximum description reached', DATA.GUI.Segment);
%       if (isempty(regexpi('roiput',currenttool, 'once')))
        return
%       end
    end
    set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
  end
end

%Points
if ~isempty(regexpi('point',currenttool)) && strcmp(DATA.ViewPanelsType{panel},'one')
  %check number of points in slice current timeframe if in single frame
  %mode.
  %check maximum number of points in slice all timeframes.
  pinslices= SET(no).Point.Z==SET(no).CurrentSlice;

  if not(DATA.ThisFrameOnly)    %findindented(DATA.Handles.hideiconholder,'allframesmode')
    %all frames mode
    uniquetf = unique(SET(no).Point.T(pinslices));
    numpoints = 0;
    for t = uniquetf(~isnan(uniquetf))
      tmp = sum(ismember(SET(no).Point.T(pinslices),[t,nan]));
      if numpoints<tmp
        numpoints = tmp;
      end
    end
  else
    numpoints = sum(SET(no).Point.T(pinslices)==SET(no).CurrentTimeFrame);
  end

  if numpoints>=length(DATA.Handles.pointtext)
    myfailed('Cannot have more than twenty points in a slice and timeframe')
    return
  end
end

%Measurements
if  ~isempty(SET(no).Measure) && ~isempty([SET(no).Measure.X])
  [measure,slice] = viewfunctions('getmeasurecoords',panel);
  nummeasures = 0;
  if not(DATA.ThisFrameOnly)  %findindented(DATA.Handles.hideiconholder,'allframesmode')
    %all frames mode
    for tloop = 1:SET(no).TSize
      nummeasuresintf = 0;
      for loop=1:length(SET(no).Measure)
        ziv = round(measure(loop).Z);
        ziv = min(ziv):max(ziv);

        if ismember(slice,ziv) && ...
            (SET(no).Measure(loop).T==tloop)||(isnan(SET(no).Measure(loop).T))
          nummeasuresintf = nummeasuresintf + 1;
        end
      end

      if nummeasuresintf>nummeasures
        nummeasures = nummeasuresintf;
      end
    end
  else
    for loop=1:length(SET(no).Measure)
      ziv = round(measure(loop).Z);
      ziv = min(ziv):max(ziv);
      if ismember(slice,ziv) && ...
          (SET(no).Measure(loop).T==SET(no).CurrentTimeFrame)||(isnan(SET(no).Measure(loop).T))
        nummeasures = nummeasures + 1;
      end
    end
  end

  if nummeasures>=length(DATA.Handles.measurementtext)
    myfailed('Cannot have more than twenty measurements in a slice and timeframe')
    return
  end
end

%Here we toggle to different buttondowns.
switch currenttool
  case 'scale'
    buttondownfunctions('scale_buttondown',panel)
  case 'move'
    buttondownfunctions('translate_buttondown',panel)
  case 'select'
    buttondownfunctions('select_buttondown',panel)
  case {'Endo','Epi','RVEndo','RVEpi','LA','RA','Roi','Scar','MO','MaR','ScarRubber','MORubber','MaRRubber','GeneralPen', 'CalciumPen', 'CalciumPenRemove', 'CalciumPenBlue', 'CalciumPenRed', 'CalciumPenLilac' }
    buttondownfunctions('pen_buttondown',panel,currenttool)
  case {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp','GeneralPenInterp','LAInterp','RAInterp'}
    buttondownfunctions('interp_buttondown',panel,currenttool)
  case 'Measure'
    buttondownfunctions('measure_buttondown',panel)
  case 'Contrast'
    buttondownfunctions('contrast_buttondown',panel)
  case 'RoiPut'
    buttondownfunctions('putroi_buttondown',panel,'Roi')
  case 'RoiBalloon'
    buttondownfunctions('balloon_buttondown',panel,'Roi')
  case 'EndoBalloon'
    buttondownfunctions('balloon_buttondown',panel,'Endo')
  case 'RVEndoBalloon'
    buttondownfunctions('balloon_buttondown',panel,'RVEndo')
  case 'EpiBalloon'
    buttondownfunctions('scale_buttondown',panel)
    %     buttondownfunctions('balloon_buttondown',panel,'Epi')
  case 'crop'
    %can only crop in one and montage viewmodes
    if ~any(strcmp(DATA.ViewPanelsType{panel},{'one','montage','montagesegmented','montagerow'}))
      return
    end
    buttondownfunctions('crop_buttondown',panel)
  case 'Point'
    buttondownfunctions('point_buttondown',panel)
    %case 'cutvessel'
    %Not implemented
    %segment3dp.tools('enableundo');
    %segment3dp.buttondownfunctions('normalclicksurface',panel); %cutvessel_buttondown
  case 'moveall'
    buttondownfunctions('translateall_buttondown',panel)
  case 'click3d'
    drawfunctions('drawpoint3D',panel)
  case 'addcalcium'
    addcalcium_buttondown;
  case 'addaortacalcium'
    addaortacalcium_buttondown;
  case 'addmitralcalcium'
    addmitralcalcium_buttondown;
  case 'addcoronarycalcium'
    addcoronarycalcium_buttondown;
end

%---------------------------------------------
function updatebuttondowns(currenttool,panels) 
%---------------------------------------------
%The buttondown is always used within a panel. A buttondown should
%therefore always have the panel as the first argument.

global DATA

if nargin == 0 || isempty(currenttool)
  currenttool = DATA.CurrentTool;
end

%if still is empty then it hasnt been set yet then assume select.
if isempty(currenttool)
  currenttool = 'select';
end

% make sure to undent contrast icon
if ~strcmpi(currenttool,'contast') && ~DATA.issegment3dp
  runicon = 0;
  undent(DATA.Handles.permanenticonholder,'contrastbrightness',runicon);
end

if nargin<2
  panels = 1:length(DATA.ViewPanels);
end

%connect interpolation
tools('connectallinterpolation');

%Set the correct cursor type
switch currenttool
  case {'scale','move','moveall','Contrast'}
    load('pointers.mat','pointer');
    set(DATA.imagefig,...
      'pointer','custom',...
      'pointershapecdata',1+pointer.(lower(currenttool)),...
      'pointershapehotspot',[7 7]);
    
  case {'select'}
    set(DATA.imagefig,'pointer','arrow');
    
  case {'RoiPut'}
    load('pointers.mat','pointer');
    set(DATA.imagefig,...
      'pointer','custom',...
      'pointershapecdata',1+pointer.point,...
      'pointershapehotspot',[7 7]);
    
  case {'Endo','Epi','RVEndo','RVEpi','Roi','Scar','MO','MaR','crop',...
      'ScarRubber','MORubber','MaRRubber','LA','RA','LAInterp','RAInterp',...
      'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp','GeneralPenInterp',...
      'Measure','RoiBalloon','EndoBalloon', 'EpiBalloon','GeneralPen',...
      'CalciumPen','CalciumPenBlue','CalciumPenRed','CalciumPenLilac','CalciumPenRemove',...
      'addcalcium','addaortacalcium','addmitralcalcium','addcoronarycalcium'}
    if isa(DATA.GUISettings.DrawPointer,'string')||isa(DATA.GUISettings.DrawPointer,'char')
      set(DATA.imagefig,'pointer',DATA.GUISettings.DrawPointer);
    else
      set(DATA.imagefig,'pointer','custom',...
        'pointershapecdata',DATA.GUISettings.DrawPointer.cdata,...
        'pointershapehotspot',DATA.GUISettings.DrawPointer.hotspot);
    end
end

viewfunctions('updatetoolhidestate',currenttool); 
%Set the buttondowns on each panel
for p = panels
  axeschildren = DATA.Handles.imageaxes(p).Children;
  axeschildren = findobj(axeschildren,'-not','Tag','measurementtext'); % this excludes measurement text
  set(axeschildren,'ButtonDownFcn',sprintf('buttondownfunctions(''buttondown'',%d, ''%s'')',p,currenttool))
  set(DATA.Handles.imageaxes(p),'ButtonDownFcn',sprintf('buttondownfunctions(''buttondown'',%d,''%s'')',p,currenttool))
end

%Reset previous motion and buttonup functions just for sure
%set(DATA.fig,'WindowButtonMotionFcn','');
set(DATA.fig,'WindowButtonMotionFcn','segment(''toggleplaceholdermotion'')');
set(DATA.fig,'WindowButtonUpFcn','');
      
%Set the CurrentTool field
DATA.CurrentTool = currenttool;

%--------------------------------------
function measuretextbuttondown(panel,h)
%--------------------------------------
global DATA SET

if ~isequal(DATA.CurrentPanel,panel)
  viewfunctions('switchpanel',panel)
  DATA.CurrentPanel = panel;
  if ~DATA.issegment3dp
    DATA.updatevolumeaxes;
  end
  return
end
clicktype = get(DATA.fig,'SelectionType');
no = DATA.ViewPanels(panel);
type = 'Measure';

str = h.String;
objectind = h.UserData;
if isempty(objectind) ||length(objectind) > 1
  measlength = compose("%.1f", [SET(no).Measure.Length]);
  measstr = strsplit(str{2}); % separates mm from the value
  tmpind = find(contains(measlength,measstr{1}));
  if isempty(objectind)
    objectind = tmpind;
  else
    objectind = intersect(tmpind,objectind);
  end
end
if isempty(objectind)
  objectind = callbackfunctions('measureaskhelper');
  if objectind == 0
    return
  end
end
switch clicktype
  case 'alt'
    % right-click, then show measurement context menu
    if isempty(DATA.LastObject)
      DATA.LastObject = objectind;
    else
      DATA.LastObject(1) = objectind;
    end
    viewfunctions('placecontextmenu',panel,type,objectind);
  case 'normal'
    % left click, allow to move text
    set(DATA.fig,'WindowButtonMotionFcn',@(obj,event)motionfunctions('measuretext_motion',h,panel,objectind));
    set(DATA.fig,'WindowButtonUpFcn',sprintf('buttonupfunctions(''measuretext_buttonup'',%d,%d)',panel,objectind));
  otherwise
end

%-----------------------------------
function contrast_buttondown(panel) 
%-----------------------------------
%Activated by contrast tool, different from resetlight_Callback (above).

global DATA

%first we get the clicked position
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,x,y);

%get x and y size of panel to normalize distance with
xsize = get(DATA.Handles.imageaxes(panel),'xlim');
xsize = xsize(end)-xsize(1);
ysize = get(DATA.Handles.imageaxes(panel),'ylim');
ysize = ysize(end)-ysize(1);

%we set the motion function and the buttonupfunction with input start
%position
set(DATA.fig,'WindowButtonMotionFcn',sprintf('motionfunctions(''contrast_motion'',%d,%d,%f,%f,%f,%f)',panel,slice,x,y,xsize,ysize));
set(DATA.fig,'WindowButtonUpFcn',sprintf('buttonupfunctions(''contrast_buttonup'',%d,%d,%f,%f,%f,%f)',panel,slice,x,y,xsize,ysize));

%-----------------------------------------------------------------------------------
function [interpx,interpy] = interpdrawhelper(interpx,interpy,contourx,contoury,x,y,slice)
%-----------------------------------------------------------------------------------
%Helper fcn to interpdrawbuttondown
%Side effects is update of DATA.Pin.
global DATA SET NO

% new=true;
tf = SET(NO).CurrentTimeFrame;
if nargin < 7
  slice = SET(NO).CurrentSlice;
end

%New interpolation point. User has not clicked close to old point.
%Prepare to store point
if isempty(interpx)||...
    size(interpx,1)<SET(NO).TSize||...
    size(interpx,2)<SET(NO).ZSize
  interpx = cell(SET(NO).TSize,SET(NO).ZSize);
  interpy = cell(SET(NO).TSize,SET(NO).ZSize);
end

newcontourx = [];
newcontoury = [];
origx=x;
origy=y;
if (~isempty(contourx)) && (~isnan(contourx(1,tf,slice)))
  %Segmentation exist.
  if isempty(interpx{tf,slice})
    %Contour exist but no previous points>=add many points.
    %Interpolate points
    newcontourx = contourx(:,tf,slice);
    newcontoury = contoury(:,tf,slice);
    n = size(contourx,1);
    nstep = n/DATA.Pref.NumInterpPoints; %Later take from preferences.
    newcontourx = newcontourx(round(1:nstep:n));
    newcontoury = newcontoury(round(1:nstep:n));
    newcontourx = newcontourx(:);
    newcontoury = newcontoury(:);
    y = []; %Do not add the point
    x = []; %Do not add the point
  else
    %Contour and previous points exist => find closest gap
    pinx = [interpx{tf,slice}; y];
    piny = [interpy{tf,slice}; x];
    contx = contourx(:,tf,slice);
    conty = contoury(:,tf,slice);
    
    pinxrep=repmat(pinx',[length(contx) 1]);
    contxrep=repmat(contx,[1 length(pinx)]);
    pinyrep=repmat(piny',[length(conty) 1]);
    contyrep=repmat(conty,[1 length(piny)]);
    pindist2cont = (pinxrep-contxrep).^2+(pinyrep-contyrep).^2;
    [~,mindistindex] =min(pindist2cont);
    [~,sortindex] =sort(mindistindex);
    pinx=pinx(sortindex);
    piny=piny(sortindex);
    
    interpx{tf,slice}=pinx;
    interpy{tf,slice}=piny;
    
    y = []; %Do not add the point
    x = []; %Do not add the point
  end
end

%Add points
interpx{tf,slice} = [...
  interpx{tf,slice} ; ...
  newcontourx;y];
interpy{tf,slice} = [...
  interpy{tf,slice} ; ...
  newcontoury;x];

%----------------------------------------------
function crop_buttondown(panel)
%----------------------------------------------
%New buttondown function for cropping. You cannot crop in montage viewmodes only in one.
global DATA

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,x,y);
%5 coords make a box
X = [x,x,x,x,x];
Y = [y,y,y,y,y];

DATA.CursorX = X;
DATA.CursorY = Y;

DATA.Handles.cursor.Color = 'y';
DATA.Handles.cursor.LineStyle = '--';
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);
DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''crop_motion'',%d,%f,%f)',panel,x,y);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''crop_buttonup'',%d, %f, %f, %d)',panel,x,y,slice);

%----------------------------------------------
function interp_buttondown(panel,type) 
%----------------------------------------------
%New buttondown function input is panel and type.
arguments
  panel
  type {mustBeMember(type,{'GeneralPenInterp','LAInterp','RAInterp', ...
    'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})}
end

global DATA SET

[no,slice,timeframe,x,y] = helperfunctions('getinterpparameters',panel);
closethecurve = true; %used later when resampling contour
opencontour = false;

%determine number of points for the contour
currentobjectind = 1;
switch type
  case 'GeneralPenInterp'
    currentobjectind = DATA.GeneralPenSettings.getcurrentobject;
    if isempty(currentobjectind) %object list is empty, create a new object
      currentobjectind = generalpen.generalpenfunctions('createnewobject');
    end
    datanumpoints = SET(no).GeneralPenObjects(currentobjectind).getnumpoints;
  case {'LAInterp','RAInterp'}
    if isempty(SET(no).(type(1:end-6))) %object is empty, create a new object
      generalpen.atriumpenfunctions('createnewobject',lower(type(1:end-6)));
    end
    datanumpoints = SET(no).(type(1:end-6)).getnumpoints;
    opencontour = true;
    closethecurve = false;
  otherwise
    datanumpoints = tools('getnumpointsforno',no);
end

%then we need to create a nan matrix with size [numpoints,TSize,ZSize]
tempnanmatrix = nan(datanumpoints,SET(no).TSize,SET(no).ZSize);
typeX = helperfunctions('parsesetfield',SET(no),type(1:end-6),'X',currentobjectind); %get contour
if isempty(typeX)
  helperfunctions('assignsetfield',no,type(1:end-6),'X',tempnanmatrix,currentobjectind);
  helperfunctions('assignsetfield',no,type(1:end-6),'Y',tempnanmatrix,currentobjectind);
end
%then we need to create a cell with size [TSize,ZSize]
tempcell = cell(SET(no).TSize,SET(no).ZSize);
typeInterpX = helperfunctions('parsesetfield',SET(no),type,'X',currentobjectind); %get interp contour
if isempty(typeInterpX)
  helperfunctions('assignsetfield',no,type,'X',tempcell,currentobjectind);
  helperfunctions('assignsetfield',no,type,'Y',tempcell,currentobjectind);
end

%Double check so that there is a cell with the correct size for the used
%type
tmpX = cell(SET(no).TSize,SET(no).ZSize);
tmpY = cell(SET(no).TSize,SET(no).ZSize);
typeInterpX = helperfunctions('parsesetfield',SET(no),type,'X',currentobjectind); %get interp contour
if ~all([size(typeInterpX,1) == SET(no).TSize, size(typeInterpX,2) == SET(no).ZSize])
  sz = size(typeInterpX);
  tmpX(1:sz(1),1:sz(2)) = typeInterpX(1:sz(1),1:sz(2));
  tmpY(1:sz(1),1:sz(2)) = typeInterpX(1:sz(1),1:sz(2));
  helperfunctions('assignsetfield',no,type,'X',tmpX,currentobjectind);
  helperfunctions('assignsetfield',no,type,'Y',tmpY,currentobjectind);
end

%Old files may have nan vectors in them remove them if so
typeInterpX = helperfunctions('parsesetfield',SET(no),type,'X',currentobjectind); %get interp contour
if all(isnan(typeInterpX{timeframe,slice}))
  helperfunctions('assignsetfield',no,type,'X',[],currentobjectind,timeframe,slice);
  helperfunctions('assignsetfield',no,type,'Y',[],currentobjectind,timeframe,slice);
end
switch type
  case 'GeneralPenInterp'
    [mindist,ind] = findfunctions('closestgeneralpeninterp',no,x,y,slice,timeframe);
  case {'LAInterp','RAInterp'}
    [mindist,ind] = findfunctions('closestlarainterp',panel,type,x,y,slice,timeframe);
  otherwise
    [mindist,ind] = findfunctions('closestinterp',panel,type,x,y,slice,timeframe);
end
typeX = helperfunctions('parsesetfield',SET(no),type(1:end-6),'X',currentobjectind); %get contour
typeY = helperfunctions('parsesetfield',SET(no),type(1:end-6),'Y',currentobjectind); %get contour
typeInterpX = helperfunctions('parsesetfield',SET(no),type,'X',currentobjectind); %get interp contour
typeInterpY = helperfunctions('parsesetfield',SET(no),type,'Y',currentobjectind); %get interp contour

if ~isempty(mindist) && mindist<DATA.Pref.ContourAdjustDistance
  %adjust old point
  minind = ind;
elseif isempty(typeInterpX{timeframe,slice}) ...
  && ~isempty(typeX) ...
  && ~all(isnan(typeX(:,timeframe,SET(no).CurrentSlice)))
  %Contour exist but no previous points>=add many points.
  %Interpolate points
  newcontourx = typeX(:,timeframe,slice);
  newcontoury = typeY(:,timeframe,slice);
  n = size(typeX,1);
  if strcmpi(type,'endointerp') && any(contains({'4CH','2CH','3CH'},SET(no).ImageViewPlane,"IgnoreCase",true))
    numinterpoints = max(54,DATA.Pref.NumInterpPoints);
  else
    % Take from preferences.
    numinterpoints = DATA.Pref.NumInterpPoints;
  end
  nstep = n/numinterpoints;
  valuetoassignX = newcontourx(round(1:nstep:n));
  valuetoassignY = newcontoury(round(1:nstep:n));
  if opencontour
    % make sure to add the last point from the contour to open contour
    valuetoassignX = [valuetoassignX; newcontourx(end)];
    valuetoassignY = [valuetoassignY; newcontoury(end)];
  end
  helperfunctions('assignsetfield',no,type,'X',valuetoassignX,currentobjectind,timeframe,slice);
  helperfunctions('assignsetfield',no,type,'Y',valuetoassignY,currentobjectind,timeframe,slice);

  %do graphical update and return
  viewfunctions('updatedrawlist',panel);
  drawfunctions('drawinterp',panel,type)
  return
else
  %closest point on the contour then use this to find the two closest
  %interpolation points in different directions along the contour this
  %should remove the case where you have close turns on the contour i.e the
  %horseshoe
  isinterpongoing = helperfunctions('isinterpongoing',SET(no),type,currentobjectind);
  [val,~] = min(sqrt((x-typeX(:,timeframe,slice)).^2 + (y-typeY(:,timeframe,slice)).^2));
  if ~isnan(val) && ~isinterpongoing%val< DATA.Pref.ContourAdjustDistance
    %First we find the closest point on the contour then we consider
    %the distance to each point from the closest points. Since contour
    %is equidistantly sampled we can use indexes as our metric. To
    %ensure that there is a point close when clicking next to the
    %contour we upsample the contour to a thousand points with linear interpolation.
    datax = typeX(:,timeframe,slice);
    datay = typeY(:,timeframe,slice);
    
    valueind = ~isnan(datax);
    datax = datax(valueind);
    datay = datay(valueind);
    closethecurve = false; %curve is already closed
    [contourx,contoury] = calcfunctions('resamplecurve',datax,datay,1000,opencontour,closethecurve);
    [~,ind] = min(sqrt((x-contourx).^2+(y-contoury).^2));
    minind2contour = ind; %for LA/RA to determine if point should be placed adjacent to first/last point without taking into account the second shortest distance
    p2pcontour = zeros(1,length(contourx));
    p2pcontour(1:ind) = ind - 1:-1:0;
    p2pcontour(ind+1:end) = 1:length(p2pcontour)-ind;
    p2pcontour(end-(floor(length(p2pcontour)/2)-ind)+1:end) = floor(length(p2pcontour)/2):-1:ind + 1;
    %then we map each interpolation point to a distance on p2pcontour
    %by finding there equivalent on the contour
    interpdistalongcontourmap = zeros(1,length(typeInterpY{timeframe,slice})); 
    values = interpdistalongcontourmap;
    for i = 1:length(typeInterpY{timeframe,slice})
      [val,ind] = min((typeInterpX{timeframe,slice}(i) - contourx).^2 ...
        + (typeInterpY{timeframe,slice}(i) - contoury).^2);
      interpdistalongcontourmap(i) = ind;
      values(i) = val;
    end
    %find the two minimal distance points
    [~,inds] = sort(p2pcontour(interpdistalongcontourmap));
    if ~opencontour
      %JB-2020-08-06 attempt to solve #2107
      % use the closest point ind(1) -> determine index of the
      indstart = inds(1);
      pstart = [typeInterpX{timeframe,slice}(indstart),...
        typeInterpY{timeframe,slice}(indstart)];
      if indstart == 1
        indbefore = length(inds);
      else
        indbefore = indstart-1;
      end
      pbefore = [typeInterpX{timeframe,slice}(indbefore),...
        typeInterpY{timeframe,slice}(indbefore)];
      if indstart == length(inds)
        indafter = 1;
      else
        indafter = indstart+1;
      end
      pafter = [typeInterpX{timeframe,slice}(indafter),...
        typeInterpY{timeframe,slice}(indafter)];
      pnew = [x,y];
      % calculate angle between closest point-newpoint-left neighbor
      %     angleft = atan2d((det([pnew-pstart;pleft-pnew])),dot(pnew-pstart,pleft-pnew))
      P0 = [pnew 0];
      P1 = [pbefore 0];
      P2 = [pstart 0];
      n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
      n2 = (P1 - P0) / norm(P1 - P0);
      angleft = atan2d(norm(cross(n1, n2)), dot(n1, n2));

      P1 = [pafter 0];
      n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
      n2 = (P1 - P0) / norm(P1 - P0);
      angright = atan2d(norm(cross(n1, n2)), dot(n1, n2));

      anglelimit = 60;
      if (indbefore == inds(2) && angleft < anglelimit)
        % switch the order of the points
        inds(2) = indafter;
      elseif(indafter == inds(2) && angright < anglelimit)
        inds(2) = indbefore;
      end
      % end of attempt to solve #2107
    end
    %this should capture the case between the last and first endpoint
    if abs(inds(1)-inds(2)) > 1 && ~opencontour
      indl = 0;
      valuetoassignX = [x;typeInterpX{timeframe,slice}];
      valuetoassignY = [y;typeInterpY{timeframe,slice}];
    elseif inds(1) == 1 && opencontour && minind2contour == 1
      %this is the case for LA/RA when clicking close to first point
      indl = 0;
      valuetoassignX = [x;typeInterpX{timeframe,slice}];
      valuetoassignY = [y;typeInterpY{timeframe,slice}];
    elseif inds(1) == length(inds) && opencontour && minind2contour == 1000
      %this is the case for LA/RA when clicking close to last point
      indl = 0;
      valuetoassignX = [typeInterpX{timeframe,slice};x];
      valuetoassignY = [typeInterpY{timeframe,slice};y];
    else
      indl = min(inds(1:2));
      valuetoassignX = [typeInterpX{timeframe,slice}(1:indl); x; typeInterpX{timeframe,slice}(indl+1:end)];
      valuetoassignY = [typeInterpY{timeframe,slice}(1:indl); y; typeInterpY{timeframe,slice}(indl+1:end)];
    end
    helperfunctions('assignsetfield',no,type,'X',valuetoassignX,currentobjectind,timeframe,slice);
    helperfunctions('assignsetfield',no,type,'Y',valuetoassignY,currentobjectind,timeframe,slice);
    minind = indl+1;
  else
    minind = length(typeInterpX{timeframe,slice}) + 1;
    helperfunctions('assignsetfield',no,type,'X',x,currentobjectind,timeframe,slice,minind);
    helperfunctions('assignsetfield',no,type,'Y',y,currentobjectind,timeframe,slice,minind);
  end
end

%write the results back to the contour field in the SET struct if we've
%placed more than one point
threshold = 1;
helperfunctions('storeinterptocontour',no,panel,type,timeframe,slice,currentobjectind,threshold);

drawfunctions('drawinterp',panel,type);

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''interp_motion'',%d,''%s'',%d,%d)',panel,type,slice,minind);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''interp_buttonup'',%d,''%s'',%d)',panel,type,slice);

%-------------------------------------
function orthoview_buttondown(panel,x,y)
%-------------------------------------
%Rotate GLA view using handle on intersection line
global DATA SET

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);

switch DATA.ViewPanelsType{panel}
  case 'orth'
    %get the angle handle position
    xhan = DATA.Handles.orthoanglehandle.XData;
    yhan = DATA.Handles.orthoanglehandle.YData;
    
    %If sufficiently close to rotation handle trigger rotation else update center position
    if norm([xhan yhan]-[x,y])<3
      set(DATA.Handles.cursor,'markersize',20,'marker', '.', 'color', 'w','parent',DATA.Handles.imageaxes(panel))
      set(DATA.fig,'WindowButtonMotionFcn',sprintf('motionfunctions(''glarotatehandle_Motion'',%d)',panel));
      set(DATA.fig,'WindowButtonUpFcn',sprintf('buttonupfunctions(''glarotatehandle_Buttonup'',%d)',panel));
      return
    end
    
    SET(no).HLA.slice = round(y/scale);
    SET(no).VLA.slice = round(x/scale);
    viewfunctions('setglacenter',no);
  case 'hla'
    SET(no).VLA.slice = round(x/scale);
    SET(no).CurrentSlice = round(y/scale);
    viewfunctions('setglacenter',no);
  case 'gla'
    [vslice,hslice] = calcfunctions('gla2sax',y/scale,x/scale,no);
    SET(no).HLA.slice = round(vslice);
    SET(no).VLA.slice = round(hslice);
    SET(no).CurrentSlice = round(y/scale);
    viewfunctions('setglacenter',no);
  case 'vla'
    SET(no).HLA.slice = round(x/scale);
    SET(no).CurrentSlice = round(y/scale);
    viewfunctions('setglacenter',no);
end

%updates the plane intersections
drawfunctions('drawplaneintersections');

%update the angle handle for orthoview
DATA.Handles.orthoanglehandle.XData = nan;
DATA.Handles.orthoanglehandle.YData = nan; % this forces the orthoangle handle to refresh after translation
drawfunctions('draworthoanglehandle',find(strcmp('orth',DATA.ViewPanelsType)))

%empty the DATA.viewim to force creation of new images in all panels but
%the current as it is unnecessary.
for p = 1:length(DATA.ViewPanels)
  if p~=panel
    calcfunctions('segmentationintersection_helper',p,SET(no).CurrentTimeFrame) %this updates segmentation intersections for the selected panel and current timeframe so it is necessary that intersections are updated for all files at buttonup at buttonup.
    DATA.ViewIM{p} = [];
    drawfunctions('drawpanel',p)
  end
end

%force graphics update
drawnow

%after doing the current timeframe update we do it for all timeframes
for p = 1:length(DATA.ViewPanels)
  if p~=panel
    for t = 1:SET(no).TSize
      calcfunctions('segmentationintersection_helper',p,t) %this updates stored segmentation intersections for remaining timeframes no graphical update needed.
    end
  end
end

%------------------------------
function scale_buttondown(panel) 
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);
if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  magno = SET(no).Flow.MagnitudeNo;
else
  magno = no;
end
scale = viewfunctions('getscale',panel);

%this makes the clicked roi current roi
select_buttondown(panel)
[ystart,xstart] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',panel,ystart,xstart);
slices = viewfunctions('slicesinpanel',panel);
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
imdim = zoomfunctions.getxysize(no,panel);

xt = ((xl-1)*imdim.XSize- imdim.XStart +1);
yt = ((yl-1)*imdim.YSize- imdim.YStart +1);

x = xstart/scale - xt;
y = ystart/scale - yt;

[type,objectind] = findfunctions('closestobject',panel,x,y,slice);

%below follows the scalables
switch type
  case {'Endo','Epi','RVEndo','RVEpi','EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'}
    %objectind is not used for contours
    if contains(type,'Interp')
      tools('connectinterpolation',no,{type});
      %delete Interp at the end
      type = type(1:end-6); 
    else
      tools('connectinterpolation',no,{[type,'Interp']});
    end
    DATA.CursorX = xt +...
      scale*SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice);
    DATA.CursorY = yt + ...
      scale*SET(no).([type,'Y'])(:,SET(no).CurrentTimeFrame,slice);
    
  case 'Roi'
    DATA.CursorX = xt + ...
      scale*SET(magno).Roi(objectind).X(:,SET(no).CurrentTimeFrame);
    DATA.CursorY = yt + ...
      scale*SET(magno).Roi(objectind).Y(:,SET(no).CurrentTimeFrame);
    
  otherwise
    return
end
xc = mean(DATA.CursorX);
yc = mean(DATA.CursorY);
startrad = norm([xstart,ystart]-[xc,yc]);

%Set marker color
DATA.Handles.cursor.Color = 'w';
DATA.Handles.cursor.Marker = 'none';

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.Handles.cursor.YData = DATA.CursorX;
DATA.Handles.cursor.XData = DATA.CursorY;

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''scale_motion'',%d,%f)',panel,startrad);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''scale_buttonup'',%d,''%s'',%d,%d,%f)',panel,type,objectind,slice,startrad);

%---------------------
function translate_buttondown(panel) 
%-----------------------
global DATA SET
persistent paninmontagemsg; %Variable for the panning message in the montage view.

scale = viewfunctions('getscale',panel);
no = DATA.ViewPanels(panel);

%fix so possible to translate roi in phase image stack
if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  magno = SET(no).Flow.MagnitudeNo;
else
  magno = no;
end

[ystart,xstart] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',panel,ystart,xstart);
slices = viewfunctions('slicesinpanel',panel);
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
imdim = zoomfunctions.getxysize(no,panel);

xt = ((xl-1)*imdim.XSize- imdim.XStart +1);
yt = ((yl-1)*imdim.YSize- imdim.YStart +1);

x = xstart/scale - xt;
y = ystart/scale - yt;

[type,objectind] = findfunctions('closestobject',panel,x,y,slice);

%below follows the translateables
switch type
  case {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'} %does the same as the later but need to handle dynamic calling differently
    tools('connectinterpolation',no,{type});
    %objectind is not used for contours
    DATA.CursorX = xt +...
      scale*SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice);
    DATA.CursorY = yt + ...
      scale*SET(no).([type(1:end-6),'Y'])(:,SET(no).CurrentTimeFrame,slice);
    
    %When translating we remove the interpolation points
    if ~isempty(SET(no).([type,'X']))
      SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice} = [];
      SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice} = [];
      drawfunctions('drawinterp',panel)
    end
    
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = 'none';
    
  case {'Endo','Epi','RVEndo','RVEpi'}
    tools('connectinterpolation',no,{[type,'Interp']});
    %objectind is not used for contours
    DATA.CursorX = xt +...
      scale*SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice);
    DATA.CursorY = yt + ...
      scale*SET(no).([type,'Y'])(:,SET(no).CurrentTimeFrame,slice);
    
    %When translating we remove the interpolation points
    if ~isempty(SET(no).([type,'InterpX']))
      SET(no).([type,'InterpX']){SET(no).CurrentTimeFrame,slice} = [];
      SET(no).([type,'InterpY']){SET(no).CurrentTimeFrame,slice} = [];
      drawfunctions('drawinterp',panel)
    end
    
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = 'none';
    
  case 'Roi'
    DATA.CursorX = xt + ...
      scale*SET(magno).Roi(objectind).X(:,SET(no).CurrentTimeFrame);
    DATA.CursorY = yt + ...
      scale*SET(magno).Roi(objectind).Y(:,SET(no).CurrentTimeFrame);
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = 'none';
  case 'Measure'
    %measures
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = '+';
    [measure,~] = viewfunctions('getmeasurecoords',panel);
    DATA.CursorX = xt + scale * measure(objectind).X;
    DATA.CursorY = yt + scale * measure(objectind).Y;
  case 'Point'
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = '+';
    DATA.CursorX = xt + scale * SET(no).Point.X(objectind);
    DATA.CursorY = yt + scale * SET(no).Point.Y(objectind);
  case 'Center'
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = '+';
    DATA.CursorX = xt + scale * SET(no).CenterX;
    DATA.CursorY = yt + scale * SET(no).CenterY;
    %When translating we remove the center point
    if ~isempty(SET(no).([type,'X']))
      % set to a very high number so that the original center cross is not
      % seen during motion
      SET(no).([type,'X']) = 10^6;
      SET(no).([type,'Y']) = 10^6;
      drawfunctions('drawcentercross',panel)
    end
  otherwise
    %cannot pan in montage viewmodes
    if any(strcmp(DATA.ViewPanelsType{panel},{'montage','montagesegmented','montagerow'}))
      if isempty(paninmontagemsg)
        mywarning('Image panning is not available in the montage view yet')
      end
      paninmontagemsg = true;
      return
    end
    pan_buttondown(panel);
    return;
end
%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.Handles.cursor.YData = DATA.CursorX;DATA.Handles.cursor.XData = DATA.CursorY;

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''translate_motion'',%d,%f,%f)',panel,xstart,ystart);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''translate_buttonup'',%d,''%s'',%d,%f,%f,%d)',panel,type,objectind,xstart,ystart,slice);


%----------------------------------------------
function balloon_buttondown(panel,type) 
%----------------------------------------------
%New buttondown function input is type. The types handled here
%are {Endo,Epi,RVendo,RVepi,Roi}. All temporary drawing is made using the
%cursor object.

global DATA SET

no = DATA.ViewPanels(panel);

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,y,x);
slices = viewfunctions('slicesinpanel',panel);
scale = viewfunctions('getscale',panel);
%Set the color of the cursor object which draws the temporary line according to the selected type
switch type
  case 'Endo'
    DATA.Handles.cursor.Color = 'r';
    DATA.Handles.cursor.LineStyle = '--';
    DATA.Handles.cursor.Marker = 'none';
    useconvhull = 1;
    thetasz = 200;
    tools('connectinterpolation',no,{'EndoInterp'})
  case 'Epi'
    DATA.Handles.cursor.Color = 'g';
    DATA.Handles.cursor.LineStyle = '--';
    DATA.Handles.cursor.Marker = 'none';
    useconvhull = 1;
    tools('connectinterpolation',no,{'EpiInterp'})
  case 'RVEndo'
    %check if Epi segmentation exists
    if isempty(SET(no).EpiX)|| any(isnan(SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      myfailed('LV Epi segmentation must exist in order to use this tool');
      DATA.Handles.configiconholder.undent('balloonrvendo',0);
      DATA.Handles.configiconholder.indent('rvendopen',1);
      return
    end
    tools('connectinterpolation',no,{'RVEndoInterp'})
    DATA.Handles.cursor.Color = 'm';
    DATA.Handles.cursor.LineStyle = '--';
    DATA.Handles.cursor.Marker = 'none';
    useconvhull = 1;
    thetasz = 200;
  case 'RVEpi'
    DATA.Handles.cursor.Color = 'c';
    DATA.Handles.cursor.LineStyle = '--';
    DATA.Handles.cursor.Marker = 'none';
    useconvhull = 1;
    tools('connectinterpolation',no,{'RVEpiInterp'})
  case 'Roi'
    if ~DATA.ThisFrameOnly && strcmpi(type,'roi')
      clr = 'm';
    else
      clr = 'b';
    end    
    DATA.Handles.cursor.Color = clr;
    DATA.Handles.cursor.LineStyle = '--';
    DATA.Handles.cursor.Marker = 'none';
    useconvhull = 1;
    thetasz = 200;
end

%clear interpolation points
%When drawing we remove the interpolation points
if ~any(strcmp(type,{'Roi'})) && ~isempty(SET(no).([type,'InterpX']))
  SET(no).([type,'InterpX']){SET(no).CurrentTimeFrame,slice} = [];
  SET(no).([type,'InterpY']){SET(no).CurrentTimeFrame,slice} = [];
  drawfunctions('drawinterp',panel)
end

%normalize clicked position to contour domain
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
imdim = zoomfunctions.getxysize(no,panel);
x = x/scale - ((xl-1)*imdim.XSize - imdim.XStart +1);
y = y/scale - ((yl-1)*imdim.YSize - imdim.YStart +1);

%create polar image that can be used in motion and
r = round(min([abs(x-SET(no).XSize),x-1,abs(y-SET(no).YSize),y-1]));
im = SET(no).IM(:,:,SET(no).CurrentTimeFrame,slice);

%Get polar matrix
DATA.ImP = calcfunctions('imtopolar',im(round(x)-r:round(x)+r,round(y)-r:round(y)+r),0,1,r,thetasz);

%Draw initial circle with radius 5 using Cursor before the algorithm has found anything
d = 1;
theta = linspace(0,2*pi,DATA.NumPoints);
DATA.CursorX = (d*cos(theta)+y+(yl-1)*imdim.YSize)*scale;
DATA.CursorY = (d*sin(theta)+x+(xl-1)*imdim.XSize)*scale;
DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.fig.WindowButtonMotionFcn = sprintf('lvsegmentation(''balloon_motion'',%d,%d,%f,%f,%d)',panel,slice,x,y,useconvhull);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''balloon_buttonup'',%d,''%s'',%d,%d)',panel,type,slice,SET(no).CurrentTimeFrame);


%----------------------------------------------
function putroi_buttondown(panel,type) 
%----------------------------------------------
%New buttondown function input is type. The types handled here
%are {Endo,Epi,RVendo,RVepi,Roi}. All temporary drawing is made using the
%cursor object.

global DATA SET
no = DATA.ViewPanels(panel);
[ystart,xstart] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,ystart,xstart);
scale = viewfunctions('getscale',panel);

%Set the color of the cursor object which draws the temporary line according to the selected type
if ~DATA.ThisFrameOnly && strcmpi(type,'roi')
  clr = 'm';
else
  clr = 'b';
end
DATA.Handles.cursor.Color = clr;
DATA.Handles.cursor.LineStyle = '--';

%Draw initial circle with radius 5 using Cursor before the algorithm has found anything
d = 10;
theta = linspace(0,2*pi,DATA.NumPoints);
DATA.CursorX = (d*cos(theta)'/SET(no).ResolutionY)*scale +ystart;
DATA.CursorY = (d*sin(theta)'/SET(no).ResolutionX)*scale +xstart;
DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''putroi_motion'',%d)',panel);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pen_buttonup'',%d,''%s'',%d,%d)',panel,type,slice,SET(no).CurrentTimeFrame);

%----------------------------------------------
function pen_buttondown(panel,type) 
%----------------------------------------------
%New buttondown function input is type. The types handled here
%are {Endo,Epi,RVendo,RVepi}. All temporary drawing is made using the
%cursor object.

global DATA SET
no = DATA.ViewPanels(panel);

switch type
  case 'MaRRubber'
    if isempty(SET(no).MaR)
      return
    end
    case 'MORubber'
    if isempty(SET(no).Scar)
      return
    end
    case 'ScarRubber'
    if isempty(SET(no).Scar)
      return
    end
end
%connect interpolation
tools('connectinterpolation',no,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'});
tools('connectlarainterpolation',no,'LAInterp');
tools('connectlarainterpolation',no,'RAInterp');

viewfunctions('updatetoolhidestate',type);

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,y,x);
objectind = [];
%Set the color of the cursor object which draws the temporary line according to the selected type
switch type
  case 'GeneralPen'
    objectind = DATA.GeneralPenSettings.getcurrentobject;
    if ~isempty(objectind)
      color = SET(no).GeneralPenObjects(objectind).getcolor();
    else
      color = DATA.GeneralPenSettings.getdefaultcolor;
    end
    DATA.Handles.cursor.Color = color;
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'Endo'
    DATA.Handles.cursor.Color = 'r';
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'Epi'
    DATA.Handles.cursor.Color = 'g';
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'RVEndo'
    DATA.Handles.cursor.Color = 'm';
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'RVEpi'
    DATA.Handles.cursor.Color = 'c';
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case {'LA','RA'}
    color = DATA.GeneralPenSettings.getatriumcolor(lower(type));
    DATA.Handles.cursor.Color = color;
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'Roi'
    if DATA.ThisFrameOnly
      clr = 'b';
    else
      clr = 'm';
    end
    DATA.Handles.cursor.Color = clr;
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'Scar'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'MO'
    DATA.Handles.cursor.Color = 'r';
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'MaR'
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.LineStyle = '-';
    DATA.Handles.cursor.Marker = 'none';
  case 'ScarRubber'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.Marker = 'none';
  case 'MORubber'
    DATA.Handles.cursor.Color = 'r';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.LineWidth = 2;
    DATA.Handles.cursor.Marker = 'none';
  case 'MaRRubber'
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.Marker = 'none';
  case 'CalciumPen'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.LineWidth = 2;
    DATA.Handles.cursor.Marker = 'none';    
  case 'CalciumPenBlue'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.LineWidth = 2;
    DATA.Handles.cursor.Marker = 'none';   
  case 'CalciumPenRed'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.LineWidth = 2;
    DATA.Handles.cursor.Marker = 'none';   
  case 'CalciumPenLilac'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.LineWidth = 2;
    DATA.Handles.cursor.Marker = 'none';   
  case 'CalciumPenRemove'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.LineWidth = 2;
    DATA.Handles.cursor.Marker = 'none';    
end

%Clear the cursorX and cursorY field
DATA.CursorX = [];
DATA.CursorY = [];

%If scar then the stack must not be timeresolved
if any(strcmp(type,{'Scar','MO','ScarRubber'}))&& SET(no).TSize>1
  myfailed('Tool not enabled for time-resolved data. Perhaps delete time frames?')
  indent(DATA.Handles.configiconholder,'select',1);
  return
end

%If scar then endocardium should exist
if any(strcmp(type,{'Scar','MO', 'MaR'}))  && isempty(SET(no).EndoX)
  myfailed('Endocardium segmentation must exist.')  
  return
end


%clear interpolation points
%When drawing we remove the interpolation points
closeinterpolationcontour(no,panel,type,slice,objectind);

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

%Update Callbacks
DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''pen_motion'',%d)',panel);
if strcmp(type,'MaR')
  DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pen_buttonup'',%d,''%s'')',panel,type);
else
  DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pen_buttonup'',%d,''%s'',%d,%d)',panel,type,slice,SET(no).CurrentTimeFrame);%last input is doall timeframes
end

%----------------------------------------------
function closeinterpolationcontour(no,panel,type,slice,objectind)
%----------------------------------------------
%Helper function to close interpolation contour when calling pen_buttondown
global SET
if any(strcmp(type,{'Scar','MaR','MO','MORubber','ScarRubber','MaRRubber','Roi',...
    'CalciumPen','CalciumPenBlue','CalciumPenRed','CalciumPenLilac','CalciumPenRemove'}))
  return
end

timeframe = SET(no).CurrentTimeFrame;
switch type
  case 'GeneralPen'
    if (~isempty(SET(no).GeneralPenObjects(objectind)) && ...
        ~isempty(SET(no).GeneralPenObjects(objectind).InterpX))
      SET(no).GeneralPenObjects(objectind).InterpX{timeframe,slice} = [];
      SET(no).GeneralPenObjects(objectind).InterpY{timeframe,slice} = [];
    end
  case {'LA','RA'}
    if (~isempty(SET(no).(type)) && ~isempty(SET(no).(type).InterpX))
      SET(no).(type).InterpX{timeframe,slice} = [];
      SET(no).(type).InterpY{timeframe,slice} = [];
    end
  otherwise
    if ~isempty(SET(no).([type,'InterpX']))
      SET(no).([type,'InterpX']){timeframe,slice} = [];
      SET(no).([type,'InterpY']){timeframe,slice} = [];
    end
end
drawfunctions('drawinterp',panel);

%----------------------------------------------
function measureadd_buttondown(panel,measureind,slice) 
%----------------------------------------------
%New buttondown function input is type. The types handled here
%are {Endo,Epi,RVendo,RVepi}. All temporary drawing is made using the
%cursor object.

global DATA 

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,x,y);
DATA.CursorX = [DATA.CursorX,x];
DATA.CursorY = [DATA.CursorY,y];
DATA.CursorZ = [DATA.CursorZ,slice];

if strcmp(DATA.fig.SelectionType, 'normal')
  buttonupfunctions('measure_buttonup',panel,measureind,slice);
  return
end

pointind = length(DATA.CursorY);
DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;
DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''measure_motion'',%d,%d)',panel,pointind);

%----------------------------------------------
function translateall_buttondown(panel) 
%----------------------------------------------
%button down for translating all existing contours

global DATA SET

%retrieve clicked point no in panel, clicked slice, slices in panel and scale
[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
no = DATA.ViewPanels(panel);
slice = viewfunctions('clickedslice',panel,y,x);
scale = viewfunctions('getscale',panel);
slices = viewfunctions('slicesinpanel',panel);

%get parameters for transferring into panel coordinates
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
imdim = zoomfunctions.getxysize(no,panel);
%set cursor color
DATA.Handles.cursor.Color = 'w';
DATA.Handles.cursor.LineStyle = '--';
DATA.Handles.cursor.Marker = 'none';
%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

%Empty out any scrap in CursorX and CursorY container.
DATA.CursorX = [];
DATA.CursorY = [];

%Retrieve all the existing contours transfer there coordinates to the panel
%coordinates and store into DATA.CursorX/Y
types = {'Endo','Epi','RVEndo','RVEpi'};

for i = 1:length(types)
  if ~isempty(SET(no).([types{i},'X'])) && ~all(isnan(SET(no).([types{i},'X'])(:,SET(no).CurrentTimeFrame,slice)))
    xt = (xl-1)*imdim.XSize - imdim.XStart +1;
    yt = (yl-1)*imdim.YSize - imdim.YStart +1;
    DATA.CursorX =[DATA.CursorX;nan; (SET(no).([types{i},'X'])(:,SET(no).CurrentTimeFrame,slice) + xt)*scale];
    DATA.CursorY =[DATA.CursorY;nan; (SET(no).([types{i},'Y'])(:,SET(no).CurrentTimeFrame,slice) + yt)*scale];
  end
end
types = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
for type = 1:length(types)
%When translating we remove the interpolation points
    if ~isempty(SET(no).([types{type},'X']))
      SET(no).([types{type},'X']){SET(no).CurrentTimeFrame,slice} = [];
      SET(no).([types{type},'Y']){SET(no).CurrentTimeFrame,slice} = [];
      drawfunctions('drawinterp',panel)
    end
end

DATA.Handles.cursor.XData = DATA.CursorY;
DATA.Handles.cursor.YData = DATA.CursorX;

DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''translateall_buttonup'',%d,%f,%f,%d)',panel,x,y,slice);
DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''translate_motion'',%d,%f,%f)',panel,x,y);

%-----------------------------
function pan_buttondown(panel)
%-----------------------------
%button down for panning

global DATA

%get starting position relative to window this is so that we can obtain the
%movement undependent of the xlim ylim
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

xw = x-DATA.Handles.imageaxes(panel).XLim(1);
yw = y-DATA.Handles.imageaxes(panel).YLim(1);

%Place the original boundaries in Cursor information fields
DATA.CursorX = DATA.Handles.imageaxes(panel).XLim;%-mean(DATA.Handles.imageaxes(panel).XLim);
DATA.CursorY = DATA.Handles.imageaxes(panel).YLim;%-mean(DATA.Handles.imageaxes(panel).YLim);

%hide all texts
DATA.Handles.text(panel).Position = [nan,nan];
% % % for c = 'cygmkwrb'
% % %   set(DATA.Handles.([c,'roitext'])(panel,:),'Position',[nan,nan])
% % % end
if ~DATA.issegment3dp
  set(DATA.Handles.roitext(panel,:),'Position',[nan,nan])
  set(DATA.Handles.measurementtext(panel,:),'Position',[nan,nan])
  set(DATA.Handles.pointtext(panel,:),'Position',[nan,nan])
end

%set motion and buttonup function
DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''pan_motion'',%d,%f,%f)',panel,xw,yw);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pan_buttonup'',%d)',panel);

%----------------------------------------------
function point_buttondown(panel) %,measureind)
%----------------------------------------------
%button down for points

global DATA SET

%Get global graphics handles
handles = DATA.Handles;
no = DATA.ViewPanels(panel);

%Get current point and slice
[y,x] = mygetcurrentpoint(handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,y,x);

%Set the color of the cursor object which draws the temporary line according to the selected type
set(handles.cursor,'Color','w','Marker','+','MarkerSize',8);

%set parent of the cursor
set(handles.cursor,'Parent',handles.imageaxes(panel));

%For 3dp the slices is empty
if ~any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','speedim','cor3DP'}))
    
  slices = viewfunctions('slicesinpanel',panel);
  scale = viewfunctions('getscale',panel);
  %get montage coordinates of clicked position
  [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
  imdim = zoomfunctions.getxysize(no,panel);

  xt = (xl-1)*imdim.XSize - imdim.XStart +1;
  yt = (yl-1)*imdim.YSize - imdim.YStart +1;

  x = x/scale -xt;
  y = y/scale -yt;
  
  %Find closest point in current timeframe and slice
  try
    [pointdist,pointind] = findfunctions('closestpoint',...
      panel,x, y,slice);
  catch
    pointdist = Inf;
    pointind = [];
  end
else
  pointdist = [];
  pointind = [];
  %if strcmp(DATA.ViewPanelsType{panel},'speedim')
  %  pointdist = [];
  %  pointind = [];
  %else
  %  %Find closest 3dp point
  %  [pointdist,pointind] = findfunctions('closestpoint3dp',panel,x,y,slice);
  %end

  %Update cursor position
  view = SET(no).LevelSet.Pen.Color;
  [r,g,b] = segment3dp.tools('xy2rgb',view,x,y); %
  [x1,y1,z1] = segment3dp.tools('rgb2xyz',r,g,b);
  segment3dp.graphics('setcursor',x1,y1,z1);
  %plot clicked position
  set(handles.cursor,'XData',y,'YData',x);
end

%If it after the above still is empty create new measure
if isempty(pointind) || pointdist>DATA.Pref.ContourAdjustDistance
  pointind = length(SET(no).Point.X)+1;
end

%This way if you click close to a measure you drag that point of the
%measure.
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''point_buttonup'',%d,%d,%d)',panel,pointind,slice);
DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''point_motion'',%d)',panel);

%-------------------------------
function cutvessel_buttondown(~) 
%-------------------------------
%Buttondown function for cutting vessels

myfailed('3D Vessel Cut not yet implemented in 2D view.');

%----------------------------------------------
function centercross_buttondown(panel) %,measureind)
%----------------------------------------------
%New buttondown function input is type. All temporary drawing is made using the
%cursor object.

global DATA

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));

DATA.Handles.cursor.Color = 'w';
DATA.Handles.cursor.Marker = '+';
DATA.Handles.cursor.MarkerSize = 8;

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.CursorX = [x,x];
DATA.CursorY = [y,y];

%----------------------------------------------
function measure_buttondown(panel) %,measureind)
%----------------------------------------------
%New buttondown function input is type. The types handled here
%are {Endo,Epi,RVendo,RVepi}. All temporary drawing is made using the
%cursor object.

global DATA SET

no = DATA.ViewPanels(panel);

[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,x,y);
slices = viewfunctions('slicesinpanel',panel);
scale = viewfunctions('getscale',panel);

%get montage coordinates of clicked position
[xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));

%Set the color of the cursor object which draws the temporary line according to the selected type
DATA.Handles.cursor.Color = 'w';
DATA.Handles.cursor.Marker = 'o';
DATA.Handles.cursor.MarkerSize = 5;

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

%get measure number if new measurenN = nummeasures+1;
[measureind,pointind,measureX,measureY, measureZ] = findfunctions('closestmeasure',panel,x/scale - (xl-1)*SET(no).XSize...
  ,y/scale - (yl-1)*SET(no).YSize);

%If it after the above still is empty create new measure
if isempty(pointind)
  measureind = length(SET(no).Measure)+1;
  DATA.CursorX = [x,x];
  DATA.CursorY = [y,y];
  DATA.CursorZ = [slice,slice];%clicked slice
  pointind = 2;
else %this case handles grabbing of measure point aswell as extend clicking since pointind = number of coordinates + 1
  measureX(pointind) = x/scale;
  measureY(pointind) = y/scale;
  DATA.CursorY = measureX*scale; %the usual switch for correct rendering
  DATA.CursorX = measureY*scale;
  DATA.CursorZ = measureZ;
end

%This way if you click close to a measure you drag that point of the
%measure.
if strcmp(DATA.fig.SelectionType, 'normal')
  DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''measure_buttonup'',%d,%d,%d)',panel,measureind,slice);
else
  set(DATA.Handles.imageaxes(panel).Children,'ButtonDownFcn',sprintf('buttondownfunctions(''measureadd_buttondown'',%d,%d,%d)',panel,measureind,slice))
  set(DATA.Handles.imageaxes(panel),'ButtonDownFcn',sprintf('buttondownfunctions(''measureadd_buttondown'',%d,%d,%d)',panel,measureind,slice))
end

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''measure_motion'',%d,%d)',panel,pointind);

%--------------------------------
function select_buttondown(panel)
%--------------------------------
%Function that highlights the selected panel if montage also highlight the frame
global DATA SET
% %store the current panel after graphical update and return. The
% %procedure is that first the panel is selected then you can start selecting
% %slices and rois.

%no = DATA.ViewPanels(panel);

%orthoview, montage and montage segmented procedure
switch DATA.ViewPanelsType{panel}
  case 'one'
    %select and update roi
    viewfunctions('selectroi',panel);
    viewfunctions('selectgeneralpen',panel);
  case {'montage','montagesegmented','montagerow'}
    %select and update roi
    viewfunctions('selectroi',panel);
    [x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
    slice = viewfunctions('clickedslice',panel,x,y);
    
    DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''select_motion'',%d)',panel);
    DATA.fig.WindowButtonUpFcn = 'buttonupfunctions(''select_buttonup'')';
    
  case {'orth','hla','gla','vla'}
    %select and update roi
    viewfunctions('selectroi',panel);
    [x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
    orthoview_buttondown(panel,x,y);
    
  case {'trans3DP','speedim','sag3DP','cor3DP'}
    clicktype = get(DATA.fig,'SelectionType');

    switch clicktype
      case 'normal'
        %Normal click
        segment3dp.tools('storeclickedposition')
        segment3dp.graphics('update2D')
        for loop = 1:length(DATA.ViewPanels)
          if ~isequal(DATA.ViewPanelsType{loop},'speedim')
            drawfunctions('drawtext',loop);
          end
        end
      case 'alt'
        %Right click => none
    end
end

if strcmp(DATA.ProgramName, 'Segment') % only executed for segment research version
  scale = viewfunctions('getscale',panel);
  no = DATA.ViewPanels(panel);
  [ystart,xstart] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
  slice = viewfunctions('clickedslice',panel,ystart,xstart);
  slices = viewfunctions('slicesinpanel',panel);
  [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
  imdim = zoomfunctions.getxysize(no,panel);

  xt = (xl-1)*imdim.XSize - imdim.XStart +1;
  yt = (yl-1)*imdim.YSize - imdim.YStart +1;
  x = xstart/scale -xt;
  y = ystart/scale -yt;

  [type,objectind] = findfunctions('closestobject',panel,x,y,slice);
  if strcmp(type, 'Center')
    state = viewfunctions('iconson','hideplus');
    if not(state{1}) %not indented
      %Set marker color
      DATA.Handles.cursor.Color = 'w';
      DATA.Handles.cursor.Marker = '+';
      DATA.CursorX = xt + scale * SET(no).CenterX;
      DATA.CursorY = yt + scale * SET(no).CenterY;
      %When translating we remove the center point
      if ~isempty(SET(no).([type,'X']))
        % set to a very high number so that the original center cross is not
        % seen during motion
        SET(no).([type,'X']) = 10^6;
        SET(no).([type,'Y']) = 10^6;
        drawfunctions('drawcentercross',panel)
      end
      %set parent of the cursor
      DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);
      DATA.Handles.cursor.YData = DATA.CursorX;
      DATA.Handles.cursor.XData = DATA.CursorY;

      DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''translate_motion'',%d,%f,%f)',panel,xstart,ystart);
      DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''translate_buttonup'',%d,''%s'',%d,%f,%f,%d)',panel,type,objectind,xstart,ystart,slice);
     end
   end
end

%-----------------------------
function addcalcium_buttondown
%-----------------------------
%buttondown when user clicks at object

global DATA SET 

no = DATA.ViewPanels(DATA.CurrentPanel);

if isempty(SET(no).CT)
  return;
end

if ~isfield(SET(no).CT,'CalciumMask')
  return;
end

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(DATA.CurrentPanel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',DATA.CurrentPanel,y,x);

x = max(min(round(x),SET(no).XSize),1);
y = max(min(round(y),SET(no).YSize),1);
slice = round(slice);

%
[sx, sy, st, sslice]=size(SET(no).CT.CalciumMask);
calculatedpixelid=[];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y,slice)];

calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y-1,slice)];

binim=SET(no).CT.CalciumMask>0;
binim = squeeze(binim);
bw = bwconncomp(binim);

for i=1:length(bw.PixelIdxList)
  for j=1:length(calculatedpixelid)
  if ismember(calculatedpixelid(j),bw.PixelIdxList{1,i})%==1
    rightlist=bw.PixelIdxList{1,i};
    for k=1:length(rightlist)
    SET(no).CT.CalciumMask(rightlist(k))=uint8(3);
    end   
  end
  end
end
%
%SET(no).CT.CalciumMask(x,y,slice) = uint8(3);

DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel);
%Maybe later the bwconn should be precalculated and stored
%im = calcfunctions('calctruedata',SET(no).IM,no); %g?r om till Hounsfieldunits
%outputmask = uint8((im>=130));
%imbin = (im>=130); %allt ?ver 130 Hu r?knas som kalk
%imbin = squeeze(imbin);
%bw = bwconncomp(imbin);

%-----------------------------
function addmitralcalcium_buttondown
%-----------------------------
%buttondown when user clicks at object

global DATA SET 

no = DATA.ViewPanels(DATA.CurrentPanel);

if isempty(SET(no).CT)
  return;
end

if ~isfield(SET(no).CT,'CalciumMask')
  return;
end

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(DATA.CurrentPanel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',DATA.CurrentPanel,y,x);

x = max(min(round(x),SET(no).XSize),1);
y = max(min(round(y),SET(no).YSize),1);
slice = round(slice);

%
[sx, sy, st, sslice]=size(SET(no).CT.CalciumMask);
calculatedpixelid=[];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y,slice)];

calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y-1,slice)];


binim=SET(no).CT.CalciumMask>0;
binim = squeeze(binim);
bw = bwconncomp(binim);

for i=1:length(bw.PixelIdxList)
 for j=1:length(calculatedpixelid)
  if ismember(calculatedpixelid(j),bw.PixelIdxList{1,i})%==1
    rightlist=bw.PixelIdxList{1,i};
    for k=1:length(rightlist)
    SET(no).CT.CalciumMask(rightlist(k))=uint8(7);
    end   
  end
 end
end
%

DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel);

%-----------------------------
function addaortacalcium_buttondown
%-----------------------------
%buttondown when user clicks at object

global DATA SET 

no = DATA.ViewPanels(DATA.CurrentPanel);

if isempty(SET(no).CT)
  return;
end

if ~isfield(SET(no).CT,'CalciumMask')
  return;
end

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(DATA.CurrentPanel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',DATA.CurrentPanel,y,x);

x = max(min(round(x),SET(no).XSize),1);
y = max(min(round(y),SET(no).YSize),1);
slice = round(slice);

%
[sx, sy, st, sslice]=size(SET(no).CT.CalciumMask);
calculatedpixelid=[];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y,slice)];

calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y-1,slice)];

binim=SET(no).CT.CalciumMask>0;
binim = squeeze(binim);
bw = bwconncomp(binim);

for i=1:length(bw.PixelIdxList)
  for j=1:length(calculatedpixelid)
  if ismember(calculatedpixelid(j),bw.PixelIdxList{1,i})%==1
    rightlist=bw.PixelIdxList{1,i};
    for k=1:length(rightlist)
    SET(no).CT.CalciumMask(rightlist(k))=uint8(6);
    end   
  end
  end
end
%

DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel);

%-----------------------------
function addcoronarycalcium_buttondown
%-----------------------------
%buttondown when user clicks at object

global DATA SET 

no = DATA.ViewPanels(DATA.CurrentPanel);

if isempty(SET(no).CT)
  return;
end

if ~isfield(SET(no).CT,'CalciumMask')
  return;
end

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(DATA.CurrentPanel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',DATA.CurrentPanel,y,x);

x = max(min(round(x),SET(no).XSize),1);
y = max(min(round(y),SET(no).YSize),1);
slice = round(slice);

%
[sx, sy, st, sslice]=size(SET(no).CT.CalciumMask);
calculatedpixelid=[];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y,slice)];

calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x+1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x-1,y-1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y+1,slice)];
calculatedpixelid= [calculatedpixelid sub2ind(size(squeeze(SET(no).CT.CalciumMask)),x,y-1,slice)];
binim=SET(no).CT.CalciumMask>0;
binim = squeeze(binim);
bw = bwconncomp(binim);

for i=1:length(bw.PixelIdxList)
  for j=1:length(calculatedpixelid)
  if ismember(calculatedpixelid(j),bw.PixelIdxList{1,i})%==1
    rightlist=bw.PixelIdxList{1,i};
    for k=1:length(rightlist)
    SET(no).CT.CalciumMask(rightlist(k))=uint8(8);
    end   
  end
  end
end
%

DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel);

