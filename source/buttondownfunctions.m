function buttondownfunctions(varargin)
% Functions for buttondowns

% Broken out by Klas

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
feval(varargin{:}); % FEVAL switchyard

%-------------------------------------
function buttondown(panel,currenttool) %#ok<DEFNU>
%-------------------------------------
global DATA SET
persistent ismaxnumrois

scale = viewfunctions('getscale',panel);
no = DATA.ViewPanels(panel);
clicktype = get(DATA.fig,'SelectionType');

%store the current panel after graphical update and return. The
%procedure is that first the panel is selected then you can start selecting
%slices and rois.
if DATA.CurrentPanel ~= panel
  viewfunctions('switchpanel',panel)
  DATA.updatevolumeaxes;
  return
end

%if montage select slice
if any(strcmp(DATA.ViewPanelsType{panel},{'montagesegmented','montagerow','montage'}))
  [x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
  slice = viewfunctions('clickedslice',panel,x,y);
  if ~isempty(slice)
    if slice ~= SET(no).CurrentSlice && ismember(currenttool,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp','GeneralPenInterp'})
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
if any(ismember(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','cor3DP','speedim'}))
  %segment3dp.tools('enableundo')
else
  tools('enableundo',no);
end

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
    if ~strcmp(currenttool,'select')|| strcmp(DATA.CurrentTheme,'3dp') %|| isdeployed
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
    if ~strcmp(DATA.CurrentTheme,'3dp')
      if SET(no).EndoInterpOngoing || SET(no).EpiInterpOngoing ||...
          SET(no).RVEndoInterpOngoing || SET(no).RVEpiInterpOngoing     
        tools('connectinterpolation',no);
        return        
      end
      [y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
      slice = viewfunctions('clickedslice',panel,y,x);
      slices = viewfunctions('slicesinpanel',panel);
      
      %normalize clicked position to contour domain
      [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
      x = x/scale - (xl-1)*SET(no).XSize;
      y = y/scale - (yl-1)*SET(no).YSize;
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
    else
      %3DP and right click
      p = get(DATA.fig,'CurrentPoint');
      [r,g,b] = segment3dp.tools('getclickedposition3DP');
      %SET(no).LevelSet.View.RSlice = r; %Store
      %SET(no).LevelSet.View.GSlice = g;
      %SET(no).LevelSet.View.BSlice = b;
      
      [x,y,z] = segment3dp.tools('rgb2xyz',r,g,b);
      %[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel))%this needs to transformed to the stored coordinate system         
      [type,objectind] = findfunctions('closestobject',panel,x,y,z);
      switch type
        case 'Image'
          set(DATA.Handles.view2dcontextmenu,'Position',p,'Visible','on');
        case 'Point'
          DATA.LastObject = [objectind z 1]; %1=timeframe
          set(DATA.Handles.point3dpcontextmenu,'Position',p,'Visible','on');
      end
      return
    end
  case 'extend'
    if ~strcmp(DATA.CurrentTheme,'3dp')
      if SET(no).EndoInterpOngoing || SET(no).EpiInterpOngoing ||...
          SET(no).RVEndoInterpOngoing || SET(no).RVEpiInterpOngoing        
        tools('connectinterpolation',no)
        return        
      end
    end
    if ~strcmp(DATA.CurrentTool,'Measure')% && ~strcmp(DATA.CurrentTool,'spliteraser')
      %Do panning of current image
      pan_buttondown(panel);
      return
    end
end

%check if number of points, rois or measurements in current slice is larger than 20 then give
%message that this is not supported.

%Rois
if (not(isempty(regexpi('roi',currenttool))) || not(isempty(regexpi('roiput',currenttool)))...
    || not(isempty(regexpi('roiballoon',currenttool)))) && strcmp(DATA.ViewPanelsType{panel},'one')
  %check number of rois in current slice
  numrois = sum(cellfun(@(x) any(~isempty(x)) && x==SET(no).CurrentSlice,{SET(no).Roi.Z}));
  if numrois >=10 %isfield(DATA.Handles,'roitext') && numrois>=length(DATA.Handles.roitext)
    if isempty (ismaxnumrois)
      ismaxnumrois = true;
      mymsgbox('Cannot have more than 10 ROIs desriptions.','Maximum description reached', DATA.GUI.Segment);
      if (isempty(regexpi('roiput',currenttool)))
        return
      end
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
if  ~isempty(SET(no).Measure)
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
  case {'Endo','Epi','RVEndo','RVEpi','Roi','Scar','MO','MaR','ScarRubber','MORubber','MaRRubber','GeneralPen', 'CalciumPen', 'CalciumPenRemove', 'CalciumPenBlue', 'CalciumPenRed', 'CalciumPenLilac' }
    buttondownfunctions('pen_buttondown',panel,currenttool)
  case {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp','GeneralPenInterp'}
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
  case 'cutvessel'
    segment3dp.tools('enableundo');
    segment3dp.buttondownfunctions('cutvessel_buttondown',panel)
  case 'moveall'
    buttondownfunctions('translateall_buttondown',panel)
  case 'fastmarching'
    segment3dp.tools('enableundo');
    segment3dp.buttondownfunctions('fastmarching_buttondown')
  case 'wand'
    segment3dp.buttondownfunctions('wand_buttondown')
  case 'localwand'
    segment3dp.tools('enableundo');
    segment3dp.tools('storeclickedposition');
    segment3dp.tools('wandexpand_Callback');
  case 'localsmooth'
    segment3dp.tools('enableundo');
    segment3dp.buttondownfunctions('localsmooth_buttondown')
  case 'localclose'
    segment3dp.tools('enableundo');
    segment3dp.buttondownfunctions('localclose_buttondown')
  case 'localpatch'
    segment3dp.tools('enableundo');
    segment3dp.buttondownfunctions('localpatch_buttondown')
  case 'localseparate'
    segment3dp.tools('enableundo');
    segment3dp.buttondownfunctions('localseparate_buttondown')
  case 'localfill'
    segment3dp.tools('enableundo');
    segment3dp.buttondownfunctions('localfill_buttondown') 
  case 'keep'
    segment3dp.tools('keepclick_Callback')
  case {'localthreshold','draw','rubber'}
    segment3dp.tools('enableundo');
    if ~strcmp(clicktype, 'alt')
      segment3dp.buttondownfunctions('pen3dp_buttondown',currenttool)
    else
      if segment3dp.tools('is2d')
        segment3dp.tools('split2d_Callback');
      else
        segment3dp.tools('split3d_Callback');
      end
    end
  case 'crop3dp'
    crop3dp_buttondown(panel)
  case 'box3dp'
    segment3dp.tools('enableundo');
    box3dp_buttondown(panel)
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

%----------------------------------
function updatebuttondowns(currenttool,panels) %#ok<DEFNU>
%----------------------------------
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

if nargin<2
  panels = 1:length(DATA.ViewPanels);
end

%Set the correct cursor type
switch currenttool
  case {'scale','move','Contrast'}
    load('pointers.mat','pointer');
    set(DATA.imagefig,...
      'pointer','custom',...
      'pointershapecdata',1+pointer.(lower(currenttool)),...
      'pointershapehotspot',[7 7]);
    
  case 'select'
    set(DATA.imagefig,'pointer','arrow');
  case {'RoiPut','Point'}
    load('pointers.mat','pointer');
    set(DATA.imagefig,...
      'pointer','custom',...
      'pointershapecdata',1+pointer.point,...
      'pointershapehotspot',[7 7]);
    
  case {'Endo','Epi','RVEndo','RVEpi','Roi','Scar',...
      'MO','MaR','ScarRubber','MORubber','MaRRubber',...
      'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp','GeneralPenInterp',...
      'Measure','RoiBalloon','EndoBalloon', 'EpiBalloon','crop','GeneralPen','addcalcium', 'CalciumPen' , 'CalciumPenBlue', 'CalciumPenRed', 'CalciumPenLilac', 'CalciumPenRemove', 'addaortacalcium', 'addmitralcalcium', 'addcoronarycalcium'}
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
  set(DATA.Handles.imageaxes(p).Children,'ButtonDownFcn',sprintf('buttondownfunctions(''buttondown'',%d, ''%s'')',p,currenttool))
  set(DATA.Handles.imageaxes(p),'ButtonDownFcn',sprintf('buttondownfunctions(''buttondown'',%d,''%s'')',p,currenttool))
end

%Reset previous motion and buttonup functions just for sure
%set(DATA.fig,'WindowButtonMotionFcn','');
set(DATA.fig,'WindowButtonMotionFcn','segment(''toggleplaceholdermotion'')');
set(DATA.fig,'WindowButtonUpFcn','');
      
%Set the CurrentTool field
DATA.CurrentTool = currenttool;

%-----------------------------------
function contrast_buttondown(panel) %#ok<DEFNU>
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
function box3dp_buttondown(panel)
%----------------------------------------------
%New buttondown function for cropping. You cannot crop in montage viewmodes only in one.
global DATA
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
segment3dp.tools('boxmotion',x,y);
set(DATA.fig,'WindowButtonMotionFcn','segment3dp.tools(''boxmotion'')');
set(DATA.fig,'WindowButtonUpFcn','segment3dp.tools(''levelsetboxbuttonup'')');

%----------------------------------------------
function crop3dp_buttondown(panel)
%----------------------------------------------
%New buttondown function for cropping. You cannot crop in montage viewmodes only in one.
global DATA
[x,y] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
segment3dp.tools('boxmotion',x,y);
set(DATA.fig,'WindowButtonMotionFcn','segment3dp.tools(''boxmotion'')');
set(DATA.fig,'WindowButtonUpFcn','segment3dp.tools(''cropboxbuttonup'')');

%----------------------------------------------
function interp_buttondown(panel,type) %#ok<DEFNU>
%----------------------------------------------
%New buttondown function input is panel and type. The types handled here
%are {EndoInterp,EpiInterp,RVendoInterp,RVepiInterp}.

global DATA SET

no = DATA.ViewPanels(panel);
scale = viewfunctions('getscale',panel);
%get closest interpolation point
[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
slice = viewfunctions('clickedslice',panel,y,x);
slices = viewfunctions('slicesinpanel',panel);

%normalize clicked position to contour domain
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices == slice,1));
x = x/scale - (xl-1)*SET(no).XSize;
y = y/scale - (yl-1)*SET(no).YSize;

%then we need to create a nan matrix with size [numpoints,TSize,ZSize]
if isempty(SET(no).([type(1:end-6),'X']))
  SET(no).([type(1:end-6),'X']) = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
  SET(no).([type(1:end-6),'Y']) = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
end

%then we need to create a cell with size [TSize,ZSize]
if isempty(SET(no).([type,'X']))
  SET(no).([type,'X']) = cell(SET(no).TSize,SET(no).ZSize);
  SET(no).([type,'Y']) = cell(SET(no).TSize,SET(no).ZSize);
end

%Double check so that there is a cell with the correct size for the used
%type
if ~all([size(SET(no).([type,'X']),1)==SET(no).TSize,size(SET(no).([type,'X']),2)==SET(no).ZSize])
  tmpX = cell(SET(no).TSize,SET(no).ZSize);
  tmpY = cell(SET(no).TSize,SET(no).ZSize);
  sz = size(SET(no).([type,'X']));
  tmpX(1:sz(1),1:sz(2)) = SET(no).([type,'X'])(1:sz(1),1:sz(2));
  tmpY(1:sz(1),1:sz(2)) = SET(no).([type,'Y'])(1:sz(1),1:sz(2));
  SET(no).([type,'X']) = tmpX;
  SET(no).([type,'Y']) = tmpY;
end


%Old files may have nan vectors in them remove them if so
if all(isnan(SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}))
  SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice} = [];
  SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice} = [];
end

[mindist,ind] = findfunctions('closestinterp',panel,type,x,y,slice);
if ~isempty(mindist) && mindist<DATA.Pref.ContourAdjustDistance
  %adjust old point
  minind = ind;
elseif isempty(SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}) && ~isempty(SET(no).([type(1:end-6),'X'])) && ~all(isnan(SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
  %Contour exist but no previous points>=add many points.
  %Interpolate points
  newcontourx = SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice);
  newcontoury = SET(no).([type(1:end-6),'Y'])(:,SET(no).CurrentTimeFrame,slice);
  n = size(SET(no).([type(1:end-6),'X']),1);
  nstep = n/DATA.Pref.NumInterpPoints; %Later take from preferences.
  SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}= newcontourx(round(1:nstep:n));
  SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice} = newcontoury(round(1:nstep:n));
  
  %do graphical update and return
  viewfunctions('updatedrawlist',panel);
  drawfunctions('drawinterp',panel)
  return
else
  %closest point on the contour then use this to find the two closest
  %interpolation points in different directions along the contour this
  %should remove the case where you have close turns on the contour i.e the
  %horseshoe
  
  [val,~] = min(sqrt((x-SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice)).^2+(y-SET(no).([type(1:end-6),'Y'])(:,SET(no).CurrentTimeFrame,slice)).^2));
  if ~isnan(val) && ~SET(no).([type(1:end-6),'InterpOngoing'])%val< DATA.Pref.ContourAdjustDistance
    %First we find the closest point on the contour then we consider
    %the distance to each point from the closest points. Since contour
    %is equidistantly sampled we can use indexes as our metric. To
    %ensure that there is a point close when clicking next to the
    %contour we upsample the contour to a thousand points with linear interpolation.
    datax = SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice);
    datay = SET(no).([type(1:end-6),'Y'])(:,SET(no).CurrentTimeFrame,slice);
    
    valueind = ~isnan(datax);
    datax = datax(valueind);
    datay = datay(valueind);
    [contourx,contoury] = calcfunctions('resamplecurve',datax,datay,1000);
    [~,ind] = min(sqrt((x-contourx).^2+(y-contoury).^2));
    p2pcontour = zeros(1,length(contourx));
    p2pcontour(1:ind) = ind - 1:-1:0;
    p2pcontour(ind+1:end) = 1:length(p2pcontour)-ind;
    p2pcontour(end-(floor(length(p2pcontour)/2)-ind)+1:end) = floor(length(p2pcontour)/2):-1:ind + 1;
    %then we map each interpolation point to a distance on p2pcontour
    %by finding there equivalent on the contour
    interpdistalongcontourmap = zeros(1,length(SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice})); 
    values = interpdistalongcontourmap;
    for i = 1:length(SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice})
      [val,ind] = min((SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}(i)...
        -contourx).^2+...
        (SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}(i)...
        -contoury).^2);
      interpdistalongcontourmap(i) = ind;
      values(i) = val;
    end
    %find the two minimal distance points
    [~,inds] = sort(p2pcontour(interpdistalongcontourmap));
    %JB-2020-08-06 attempt to solve #2107    
    % use the closest point ind(1) -> determine index of the
    indstart = inds(1);
    pstart = [SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}(indstart),...
              SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}(indstart)];
    if indstart == 1
      indbefore = length(inds);      
    else
      indbefore = indstart-1;
    end
    pbefore = [SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}(indbefore),...
             SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}(indbefore)];
    if indstart == length(inds)
      indafter = 1;
    else
      indafter = indstart+1;
    end
    pafter = [SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}(indafter),...
             SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}(indafter)];
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
    %this should capture the case between the last and first endpoint
    if abs(inds(1)-inds(2))>1
      indl = 0;
      SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice} = [x;SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}];
      SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice} = [y;SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}];
    else
      indl = min(inds(1:2));
      SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice} = [SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}(1:indl);x;SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}(indl+1:end)];
      SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice} = [SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}(1:indl);y;SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}(indl+1:end)];
    end
    minind = indl+1;
  else
    minind = length(SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice})+1;
    SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice}(minind,1)=x;
    SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice}(minind,1)=y;
  end
end
if SET(no).([type(1:end-6),'InterpOngoing'])
  numinterppoints = length(SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice});
  numpoints = floor(DATA.NumPoints/15)*(numinterppoints-1);
  if numpoints > DATA.NumPoints
    numpoints = DATA.NumPoints;
  end
else
  numpoints = DATA.NumPoints;
end
%here we need to interpolate and store curve into contours.
%The below function removes duplicate points and resamples the contour.
[x,y] = calcfunctions('resamplecurve',SET(no).([type,'X']){SET(no).CurrentTimeFrame,slice},...
  SET(no).([type,'Y']){SET(no).CurrentTimeFrame,slice},numpoints-1);

%write the results back to the contour field in the SET struct if we've
%placed more than one point
if length(x)>1
  switch type
    case 'Roi'
      %would be nice to introduce RoiInterpX...
    otherwise
      if SET(no).([type(1:end-6),'InterpOngoing'])
        expectedlength = length(SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice));
        if numpoints < expectedlength
          % fill up with NaNs        
          x = cat(2,x,nan(1,expectedlength-numpoints));
          y = cat(2,y,nan(1,expectedlength-numpoints));
        end
      end
      SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice)= [x,x(1)];
      SET(no).([type(1:end-6),'Y'])(:,SET(no).CurrentTimeFrame,slice)= [y,y(1)];
      drawfunctions('drawcontours',panel)
  end
elseif length(x) == 1
  SET(no).([type(1:end-6),'InterpOngoing']) = true;
  drawfunctions('updateinterpolationsettings',panel,type);
end

drawfunctions('drawinterp',panel);

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
function scale_buttondown(panel) %#ok<DEFNU>
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

xl = floor(xstart/scale/SET(no).XSize)*SET(no).XSize;
yl = floor(ystart/scale/SET(no).YSize)*SET(no).YSize;
x = xstart/scale - xl;
y = ystart/scale - yl;

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
    DATA.CursorX = xl +...
      scale*SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice);
    DATA.CursorY = yl + ...
      scale*SET(no).([type,'Y'])(:,SET(no).CurrentTimeFrame,slice);
    
  case 'Roi'
    DATA.CursorX = xl + ...
      scale*SET(magno).Roi(objectind).X(:,SET(no).CurrentTimeFrame);
    DATA.CursorY = yl + ...
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
function translate_buttondown(panel) %#ok<DEFNU>
%-----------------------
global DATA SET

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
xl = floor(xstart/scale/SET(no).XSize)*SET(no).XSize;
yl = floor(ystart/scale/SET(no).YSize)*SET(no).YSize;
x = xstart/scale - xl;
y = ystart/scale - yl;

[type,objectind] = findfunctions('closestobject',panel,x,y,slice);

%below follows the translateables
switch type
  case {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'} %does the same as the later but need to handle dynamic calling differently
    tools('connectinterpolation',no,{type});
    %objectind is not used for contours
    DATA.CursorX = xl +...
      scale*SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,slice);
    DATA.CursorY = yl + ...
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
    DATA.CursorX = xl +...
      scale*SET(no).([type,'X'])(:,SET(no).CurrentTimeFrame,slice);
    DATA.CursorY = yl + ...
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
    DATA.CursorX = xl + ...
      scale*SET(magno).Roi(objectind).X(:,SET(no).CurrentTimeFrame);
    DATA.CursorY = yl + ...
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
    DATA.CursorX = xl + scale * measure(objectind).X;
    DATA.CursorY = yl + scale * measure(objectind).Y;
  case 'Point'
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = '+';
    DATA.CursorX = xl + scale * SET(no).Point.X(objectind);
    DATA.CursorY = yl + scale * SET(no).Point.Y(objectind);
  case 'Center'
    %Set marker color
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.Marker = '+';
    DATA.CursorX = xl + scale * SET(no).CenterX;
    DATA.CursorY = yl + scale * SET(no).CenterY;
    %When translating we remove the center point
    if ~isempty(SET(no).([type,'X']))
      % set to a very high number so that the original center cross is not
      % seen during motion
      SET(no).([type,'X']) = 10^6;
      SET(no).([type,'Y']) = 10^6;
      drawfunctions('drawcentercross',panel)
    end
  otherwise
    pan_buttondown(panel);
    return;
end
%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.Handles.cursor.YData = DATA.CursorX;DATA.Handles.cursor.XData = DATA.CursorY;

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''translate_motion'',%d,%f,%f)',panel,xstart,ystart);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''translate_buttonup'',%d,''%s'',%d,%f,%f,%d)',panel,type,objectind,xstart,ystart,slice);


%----------------------------------------------
function balloon_buttondown(panel,type) %#ok<DEFNU>
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
x = x/scale - (xl-1)*SET(no).XSize;
y = y/scale - (yl-1)*SET(no).YSize;

%create polar image that can be used in motion and
r = round(min([abs(x-SET(no).XSize),x-1,abs(y-SET(no).YSize),y-1]));
im = SET(no).IM(:,:,SET(no).CurrentTimeFrame,slice);

%Get polar matrix
DATA.ImP = calcfunctions('imtopolar',im(round(x)-r:round(x)+r,round(y)-r:round(y)+r),0,1,r,thetasz);

%Draw initial circle with radius 5 using Cursor before the algorithm has found anything
d = 1;
theta = linspace(0,2*pi,DATA.NumPoints);
DATA.CursorX = (d*cos(theta)+y+(yl-1)*SET(no).YSize)*scale;
DATA.CursorY = (d*sin(theta)+x+(xl-1)*SET(no).XSize)*scale;
DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.fig.WindowButtonMotionFcn = sprintf('lvsegmentation(''balloon_motion'',%d,%d,%f,%f,%d)',panel,slice,x,y,useconvhull);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''balloon_buttonup'',%d,''%s'',%d,%d)',panel,type,slice,SET(no).CurrentTimeFrame);


%----------------------------------------------
function putroi_buttondown(panel,type) %#ok<DEFNU>
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
if ~DATA.ThisFrameOnly && strcmpi(type,'roi')
  clr = 'm';
else
  clr = 'b';
end
DATA.Handles.cursor.Color = clr;
DATA.Handles.cursor.LineStyle = '--';

%normalize clicked position to contour domain
[yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
x = x/scale - (xl-1)*SET(no).XSize;
y = y/scale - (yl-1)*SET(no).YSize;

%Draw initial circle with radius 5 using Cursor before the algorithm has found anything
d = 10;
theta = linspace(0,2*pi,DATA.NumPoints);
DATA.CursorX = (d*cos(theta)'/SET(no).ResolutionY+y + (yl-1)*SET(no).YSize)*scale;
DATA.CursorY = (d*sin(theta)'/SET(no).ResolutionX+x + (xl-1)*SET(no).XSize)*scale;
DATA.Handles.cursor.XData = DATA.CursorX;
DATA.Handles.cursor.YData = DATA.CursorY;

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''putroi_motion'',%d)',panel);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pen_buttonup'',%d,''%s'',%d,%d)',panel,type,slice,SET(no).CurrentTimeFrame);


%----------------------------------------------
function pen_buttondown(panel,type) %#ok<DEFNU>
%----------------------------------------------
%New buttondown function input is type. The types handled here
%are {Endo,Epi,RVendo,RVepi}. All temporary drawing is made using the
%cursor object.

global DATA SET
no = DATA.ViewPanels(panel);
switch type
  case {'Endo','Epi','RVEndo','RVEpi'}
    tools('connectinterpolation',no,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})
end
viewfunctions('updatetoolhidestate',type);

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,y,x);
%Set the color of the cursor object which draws the temporary line according to the selected type
switch type
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
  case 'MaRRubber'
    DATA.Handles.cursor.Color = 'w';
    DATA.Handles.cursor.LineStyle = ':';
    DATA.Handles.cursor.Marker = 'none';
  case 'GeneralPen'
    DATA.Handles.cursor.Color = 'y';
    DATA.Handles.cursor.LineStyle = '-';
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
  return;
end

%If scar then endocardium should exist
if any(strcmp(type,'Scar')) && isempty(SET(no).EndoX)
  myfailed('Endocardium segmentation must exist to draw scar.')  
  return;
end

%clear interpolation points
%When drawing we remove the interpolation points
if ~any(strcmp(type,{'Scar','MaR','MO','ScarRubber','MaRRubber', 'Roi','CalciumPen','CalciumPenBlue', 'CalciumPenRed', 'CalciumPenLilac', 'CalciumPenRemove'})) && ~isempty(SET(no).([type,'InterpX']))
  SET(no).([type,'InterpX']){SET(no).CurrentTimeFrame,slice} = [];
  SET(no).([type,'InterpY']){SET(no).CurrentTimeFrame,slice} = [];
  drawfunctions('drawinterp',panel)
end

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''pen_motion'',%d)',panel);
if strcmp(type,'MaR')
  DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pen_buttonup'',%d,''%s'')',panel,type);
else
  DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pen_buttonup'',%d,''%s'',%d,%d)',panel,type,slice,SET(no).CurrentTimeFrame);%last input is doall timeframes
end
%----------------------------------------------
function measureadd_buttondown(panel,measureind,slice) %#ok<DEFNU>
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
function translateall_buttondown(panel) %#ok<DEFNU>
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
    DATA.CursorX =[DATA.CursorX;nan; (SET(no).([types{i},'X'])(:,SET(no).CurrentTimeFrame,slice) + (xl-1)*SET(no).XSize)*scale];
    DATA.CursorY =[DATA.CursorY;nan; (SET(no).([types{i},'Y'])(:,SET(no).CurrentTimeFrame,slice) + (yl-1)*SET(no).YSize)*scale];
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

%----------------------------------------------
function pan_buttondown(panel)
%----------------------------------------------
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
set(DATA.Handles.roitext(panel,:),'Position',[nan,nan])
set(DATA.Handles.measurementtext(panel,:),'Position',[nan,nan])
set(DATA.Handles.pointtext(panel,:),'Position',[nan,nan])

%set motion and buttonup function
DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''pan_motion'',%d,%f,%f)',panel,xw,yw);
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''pan_buttonup'',%d)',panel);

%----------------------------------------------
function point_buttondown(panel)%#ok<DEFNU> %,measureind)
%----------------------------------------------
%button down for points

global DATA SET

no = DATA.ViewPanels(panel);
[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,y,x);

%Set the color of the cursor object which draws the temporary line according to the selected type
DATA.Handles.cursor.Color = 'w';
DATA.Handles.cursor.Marker = '+';
DATA.Handles.cursor.MarkerSize = 8;

%set parent of the cursor
DATA.Handles.cursor.Parent = DATA.Handles.imageaxes(panel);

%For 3dp the slices is empty
if ~any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','speedim','cor3DP'}))
    
  slices = viewfunctions('slicesinpanel',panel);
  scale = viewfunctions('getscale',panel);
  %get montage coordinates of clicked position
  [yl,xl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slices==slice,1));
  
  %Find closest point in current timeframe and slice
  try
    [pointdist,pointind] = findfunctions('closestpoint',panel,x/scale - (xl-1)*SET(no).XSize,...
                                         y/scale - (yl-1)*SET(no).YSize,slice);
  catch
    pointdist = Inf;
    pointind = [];
  end
else
  %Find closest 3dp point
  [pointdist,pointind] = findfunctions('closestpoint3dp',panel,x,y,slice);
end
%If it after the above still is empty create new measure
if isempty(pointind) || pointdist>DATA.Pref.ContourAdjustDistance
  pointind = length(SET(no).Point.X)+1;
end

%plot clicked position
DATA.Handles.cursor.XData = y;
DATA.Handles.cursor.YData = x;

%This way if you click close to a measure you drag that point of the
%measure.
DATA.fig.WindowButtonUpFcn = sprintf('buttonupfunctions(''point_buttonup'',%d,%d,%d)',panel,pointind,slice);
DATA.fig.WindowButtonMotionFcn = sprintf('motionfunctions(''point_motion'',%d)',panel);

%-------------------------------
function cutvessel_buttondown(~) %#ok<DEFNU>
%-------------------------------
%Buttondown function for cutting vessels

myfailed('3D Vessel Cut not yet implemented in 2D view.');

%----------------------------------------------
function centercross_buttondown(panel)%#ok<DEFNU> %,measureind)
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
function measure_buttondown(panel)%#ok<DEFNU> %,measureind)
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
        segment3dp.tools('update3DP')
        for loop = 1:length(DATA.ViewPanels)
          if ~isequal(DATA.ViewPanelsType{loop},'speedim')
            drawfunctions('drawtext',loop);
          end
        end
        %set(DATA.fig,'WindowButtonMotionFcn','motionfunctions(''select3dp_motion'')')
        %set(DATA.fig,'WindowButtonUpFcn','buttonupfunctions(''buttonup_Callback'')');
      case 'alt'
        %Right click => was cut
        %Now done by context menu
        %if segment3dp.tools('is2d')
        %  segment3dp.tools('split2d_Callback');
        %else
        %  segment3dp.tools('split3d_Callback');
        %end
    end
end

if strcmp(DATA.ProgramName, 'Segment') % only executed for segment research version
  scale = viewfunctions('getscale',panel);
  no = DATA.ViewPanels(panel);
  [ystart,xstart] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));%this needs to transformed to the stored coordinate system
  slice = viewfunctions('clickedslice',panel,ystart,xstart);
  xl = floor(xstart/scale/SET(no).XSize)*SET(no).XSize;
  yl = floor(ystart/scale/SET(no).YSize)*SET(no).YSize;
  x = xstart/scale - xl;
  y = ystart/scale - yl;

  [type,objectind] = findfunctions('closestobject',panel,x,y,slice);
  if strcmp(type, 'Center')
    state = viewfunctions('iconson','hideplus');
    if not(state{1}) %not indented
      %Set marker color
      DATA.Handles.cursor.Color = 'w';
      DATA.Handles.cursor.Marker = '+';
      DATA.CursorX = xl + scale * SET(no).CenterX;
      DATA.CursorY = yl + scale * SET(no).CenterY;
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

