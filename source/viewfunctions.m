function [varargout] = viewfunctions(varargin)
% Functions for querying about the view such as slices in panel and
% zoomstate

%#ok<*GVMIS>

% Klas
%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-------------------------------------
function slices = slicesinpanel(panel)
%-------------------------------------
%find the number of slices in panel
global DATA SET 

if panel > length(DATA.ViewPanels)
  %Should not happen
  slices = 0;
  return
end

no = DATA.ViewPanels(panel);

if isequal(no,0)
  slices = 0;
  return
end

if isempty(DATA.ViewPanelsType{panel})
  DATA.ViewPanelsType{panel} = 'one';
end
switch DATA.ViewPanelsType{panel}
  case 'montagesegmented'
    slices = segment_main('getmontagesegmentedslices',no);
  case {'montage','montagerow'}
    slices = 1:SET(no).ZSize;
  case {'one','orth'}
    slices = SET(no).CurrentSlice;
  case 'gla'
    slices = 0;
  case 'vla'
    slices = SET(no).VLA.slice;
  case 'hla'
    slices = SET(no).HLA.slice;
  otherwise
    slices = [];
end

%-------------------------------
function scale = getscale(panel)
%-------------------------------
%Returns if the panel should be interpolated or not

global DATA

if ~DATA.Pref.ViewInterpolated || contains(DATA.ViewPanelsType{panel},'montage')
  scale = 1;
else
  scale = 2;
end

%---------------------------------------
function updatezoomandaspectratio(panel)
%---------------------------------------
%Update zoom state

global DATA SET

%Get global graphics handles
handles = DATA.Handles;

no = DATA.ViewPanels(panel);

scale = getscale(panel);
if length(DATA.ViewIM)< panel
  return
else
  %get the xsz and ysz of the images. Overwritten if 3dp view see case below
  xsz = size(DATA.ViewIM{panel},1)/scale;
  ysz = size(DATA.ViewIM{panel},2)/scale;

  switch DATA.ViewPanelsType{panel}
    case {'montage','montagesegmented','montagerow'}
      zoomstate = [0.5;ysz-0.5;0.5;xsz-0.5];
      xres = SET(no).ResolutionX;
      yres = SET(no).ResolutionY;
    case {'one','orth'}
      zoomstate = SET(no).NormalZoomState;
      xres = SET(no).ResolutionX;
      yres = SET(no).ResolutionY;
    case 'hla'
      xres = SET(no).SliceThickness + SET(no).SliceGap;
      yres = SET(no).ResolutionY;
      zoomstate = [0.5;ysz-0.5;0.5;xsz-0.5];
    case 'vla'
      xres = SET(no).SliceThickness + SET(no).SliceGap;
      yres = SET(no).ResolutionX;
      zoomstate = [0.5;ysz-0.5;0.5;xsz-0.5];
    case 'gla'
      %General longaxis. glaangle = 0 <=> HLA; glaangle = pi/2 <=> VLA
      glaangle = SET(no).GLA.angle;
      xres = SET(no).SliceThickness + SET(no).SliceGap;
      yres = SET(no).ResolutionY*cos(glaangle)+SET(no).ResolutionX*abs(sin(glaangle));
      zoomstate = updateglazoomstate(no,ysz);
  end
end

if isempty(zoomstate)
  zoomstate = getnewzoomstate(panel,no);
else
  % this is the previous zoom state
  if ~isempty(SET(no).NormalZoomState) %&& strcmp(DATA.ViewPanelsType{panel},'one')
    zoomstateold = SET(no).NormalZoomState;
    xsizeold = zoomstateold(2)-zoomstateold(1);
    ysizeold = zoomstateold(4)-zoomstateold(3);
    zoomfold = xsizeold/ysizeold;

    %recalculate zoom state without updating SET struct
    zoomstatenew = getnewzoomstate(panel,no,false);
    xsizenew = zoomstatenew(2)-zoomstatenew(1);
    ysizenew = zoomstatenew(4)-zoomstatenew(3);
    zoomfnew = xsizenew/ysizenew;
    if ~isapproxequal(zoomfnew,zoomfold) && ~isnan(zoomfnew)
      % if the zoom ratio has changed, update to the newest
      SET(no).NormalZoomState = zoomstatenew;
      if strcmp(DATA.ViewPanelsType{panel},'one')
        zoomstate = zoomstatenew;
      end
    end
  end
end

plotboxaspectratio = [ysz/yres xsz/xres 1];

minres = min(xres,yres);
dataaspectratio = [xres/minres  yres/minres 1];

if ~DATA.Silent
  set(handles.imageaxes(panel),...  
  'XLim',scale*[zoomstate(1) zoomstate(2)],...
  'YLim',scale*[zoomstate(3) zoomstate(4)],...
  'DataAspectRatio',dataaspectratio,...
  'PlotBoxAspectRatio',plotboxaspectratio...
  );
end

%---------------------
function updatedrawlist(panel,isonlyannotation)
%-----------------------
%this function should be used at the end of a addition of annotation for
%a panel.
%Linked stacks are taken care of within annotationinpanel
do = annotationinpanel(panel);
if nargin < 2
  isonlyannotation = false;
end
updatedrawlistfcn_helper(panel,do,isonlyannotation)

%--------------------------------
function addtodrawlist(panel,fcn)
%--------------------------------
%Adds fcn to drawlist for panel

global DATA

%%only add if not already there, ismember did not work as mix of strings
%%and function handles.
count = 1;
found = false;
drawlist = DATA.drawlist{panel};
while (count<=length(drawlist)) && ~found
  if isequal(drawlist{count},fcn)
    found = true;
  end
  count = count+1;
end

%add it
if not(found)
  DATA.drawlist{panel}{end+1} = fcn;
end

%--------------------------------------------------------
function updatedrawlistfcn_helper(panel,do,isonlyannotation)
%--------------------------------------------------------
%helper function for update the drawlist function

global DATA

%Drawlist is a cell of cells
DATA.drawlist{panel} = {};

if nargin < 3
  isonlyannotation = false;
end

%image are always drawn if not viewport viewpaneltype
if ~strcmp(DATA.ViewPanelsType{panel},'viewport')
  if ~isonlyannotation
    addtodrawlist(panel,@()drawfunctions('drawimages',panel));
  end
else
  if DATA.issegment3dp
    addtodrawlist(panel,@()segment3dp.graphics('helprender'));
  else
    %Consider if we need something here
    addtodrawlist(panel,@()angioview('update'));
    
  end
end

%Centerpoint always drawn if segment
if strcmp(DATA.ProgramName, 'Segment')
  addtodrawlist(panel,@()drawfunctions('drawcentercross',panel));
end

if do.point3D
  addtodrawlist(panel,@()drawfunctions('drawpoint3D',panel));
end

switch DATA.ViewPanelsType{panel}
  case {'one','orth'}
    
    %If orth view then we need to draw the orthanglehandle
    if strcmp(DATA.ViewPanelsType{panel},'orth')
      addtodrawlist(panel,@()drawfunctions('draworthoanglehandle',panel));
    end
    
    %Contours
    types = {'Endo','Epi','RVEndo','RVEpi'};
    for i = 1:length(types)
      if do.([lower(types{i}),'contour'])
        addtodrawlist(panel,@()drawfunctions('drawcontours',panel,types{i}));
      end
    end
    if do.generalpencontour
      addtodrawlist(panel,@()drawfunctions('drawcontoursgeneralpen',panel));
    end
    types = {'LA','RA'};
    for loop = 1:length(types)
      if do.([lower(types{loop}),'contour'])
        addtodrawlist(panel,@()drawfunctions('drawcontourslara',panel,types{loop}));
      end
    end    

    %Interpolation
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen','LA','RA'};
    for i = 1:length(types)
      if do.([lower(types{i}),'interp'])
        addtodrawlist(panel,@()drawfunctions('drawinterp',panel,[types{i},'Interp']));
      end
    end

    %Others
    if do.measures
      addtodrawlist(panel,@()drawfunctions('drawmeasures',panel));
    end
    
    if do.text
      addtodrawlist(panel,@()drawfunctions('drawtext',panel));
    end
    
    if do.intersections
      addtodrawlist(panel,@()drawfunctions('drawintersections',panel));
    end
    
    drawrois = false;
    for c = 'cgymkwrb'
      if do.([c,'roi'])
        drawrois = true;
        break
      end
    end
    if drawrois
      addtodrawlist(panel,@()drawfunctions('drawroi',panel));
    end
    
    if do.point
      addtodrawlist(panel,@()drawfunctions('drawpoint',panel));
    end
    
    if do.scar
      addtodrawlist(panel,@()drawfunctions('drawviability',panel));
    else
      do.manualinteraction = 0;
      %DATA.ViewIM{panel} = [];
    end
    
    if do.manualinteraction
      addtodrawlist(panel,@()drawfunctions('showviabilityedits',panel));
    else
      %DATA.ViewIM{panel} = [];
    end
    
    if do.mar
      addtodrawlist(panel,@()drawfunctions('drawmar',panel));
    end
    
  case {'montage','montagerow','montagesegmented'}
    %Contours
    types = {'Endo','Epi','RVEndo','RVEpi'};
    for i = 1:length(types)
      if do.([lower(types{i}),'contour'])
        addtodrawlist(panel,@()drawfunctions('drawcontours',panel,types{i}));
      end
    end
    if do.generalpencontour
      addtodrawlist(panel,@()drawfunctions('drawcontoursgeneralpen',panel));
    end
    types = {'LA','RA'};
    for loop = 1:length(types)
      if do.([lower(types{loop}),'contour'])
        addtodrawlist(panel,@()drawfunctions('drawcontourslara',panel,types{loop}));
      end
    end

    %Interpolation
    types = {'Endo','Epi','RVEndo','RVEpi','LA','RA'};
    for i = 1:length(types)
      if do.([lower(types{i}),'interp'])
        addtodrawlist(panel,@()drawfunctions('drawinterp',panel,[types{i},'Interp']));
      end
    end
%     if do.generalpeninterp
%       addtodrawlist(panel,@()drawfunctions('drawinterpgeneralpen',panel));
%     end
    
    %Others
    if do.measures
      addtodrawlist(panel,@()drawfunctions('drawmeasures',panel));
    end
    
    if do.text
      addtodrawlist(panel,@()drawfunctions('drawtext',panel));
    end
    
    doroitext = false;
    for c = 'cgymkwrb'
      if do.([c,'roi'])
        addtodrawlist(panel,@()drawfunctions('drawroi',panel,c));
        doroitext = true;
      end
    end
    if doroitext
      addtodrawlist(panel,@()drawfunctions('drawroitext',panel));
    end
    
    if do.point
      addtodrawlist(panel,@()drawfunctions('drawpoint',panel));
    end

    % Setting ViewIM to empty, triggers recalculation
    % of ViewIM for each time frame. TRY NOT to use
    
    if do.scar
      addtodrawlist(panel,@()drawfunctions('drawviability',panel));
    else
      do.manualinteraction = 0;
      %DATA.ViewIM{panel} = [];
    end
    
    if do.manualinteraction
      addtodrawlist(panel,@()drawfunctions('showviabilityedits',panel));
    else
      %DATA.ViewIM{panel} = [];
    end
    
    if do.mar
      addtodrawlist(panel,@()drawfunctions('drawmar',panel));
    end
  case {'hla','gla','vla'}
    if do.text
      addtodrawlist(panel,@()drawfunctions('drawtext',panel));
    end
    
    if do.measures
      addtodrawlist(panel,@()drawfunctions('drawmeasures',panel));
    end
    
    if do.intersections
      addtodrawlist(panel,@()drawfunctions('drawintersections',panel));
    end
  case {'trans3DP','sag3DP','cor3DP','speedim'}
    
    types = {'gb','bg','gr','br','rb','rg'};
    for i = 1:length(types)
      if do.([types{i} 'line'])
        addtodrawlist(panel,@()drawfunctions('draw3dpline',panel,types{i}));
      end
    end
    
    if do.threedpcontour
      addtodrawlist(panel,@()drawfunctions('draw3dpoutline',panel,'threedp')); %remove 'threedp'?
      addtodrawlist(panel,@()drawfunctions('draw3dpsurfaceoutline',panel));
    end
    
    if do.point
      addtodrawlist(panel,@()drawfunctions('drawpoint',panel));
    end
    
end
%--------------------------------------------------------
function updatedrawlist_helper(panel,do,isonlyannotation)
%--------------------------------------------------------
%helper function for update the drawlist

global DATA

%Drawlist is a cell of cells
DATA.drawlist{panel} = {};

if nargin < 3
  isonlyannotation = false;
end

%image are always drawn if not viewport viewpaneltype
if ~strcmp(DATA.ViewPanelsType{panel},'viewport')
  if ~isonlyannotation
    addtodrawlist(panel,sprintf('drawfunctions(''drawimages'',%d)',panel));
  end
else
  addtodrawlist(panel,sprintf('segment3dp.graphics(''helprender'');'));
end

%Centerpoint always drawn if segment
if strcmp(DATA.ProgramName, 'Segment')
  addtodrawlist(panel,sprintf('drawfunctions(''drawcentercross'',%d)',panel));
end

if do.point3D
  addtodrawlist(panel,sprintf('drawfunctions(''drawpoint3D'',%d)',panel));
end

switch DATA.ViewPanelsType{panel}
  case {'one','orth'}
    
    %If orth view then we need to draw the orthanglehandle
    if strcmp(DATA.ViewPanelsType{panel},'orth')
      addtodrawlist(panel,sprintf('drawfunctions(''draworthoanglehandle'',%d)',panel));
    end
    
    %Contours
    types = {'Endo','Epi','RVEndo','RVEpi'};
    for i = 1:length(types)
      if do.([lower(types{i}),'contour'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawcontours'',%d,''%s'')',panel,types{i}));
      end
    end
    if do.generalpencontour
      addtodrawlist(panel,sprintf('drawfunctions(''drawcontoursgeneralpen'',%d)',panel));
    end
    types = {'LA','RA'};
    for loop = 1:length(types)
      if do.([lower(types{loop}),'contour'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawcontourslara'',%d,''%s'')',panel,types{loop}));
      end
    end  
    
    %Interpolation
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'interp'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawinterp'',%d,''%s'')',panel,[types{i},'Interp']));
      end
    end
%     if do.generalpeninterp
%       addtodrawlist(panel,sprintf('drawfunctions(''drawinterpgeneralpen'',%d)',panel));
%     end

    %Others
    if do.measures
      addtodrawlist(panel,sprintf('drawfunctions(''drawmeasures'',%d)',panel));
    end
    
    if do.text
      addtodrawlist(panel,sprintf('drawfunctions(''drawtext'',%d)',panel));
    end
    
    if do.intersections
      addtodrawlist(panel,sprintf('drawfunctions(''drawintersections'',%d)',panel));
    end
    
%     for c = 'cgymkwrb'
%       if do.([c,'roi'])
%         DATA.drawlist{panel}{end+1} = sprintf('drawfunctions(''drawroi'',%d,''%s'')',panel,c);
%       end
%     end
    drawrois = false;
    for c = 'cgymkwrb'
      if do.([c,'roi'])
        drawrois = true;
        break
      end
    end
    if drawrois
      addtodrawlist(panel,sprintf('drawfunctions(''drawroi'',%d)',panel));
    end
    
    if do.point
      addtodrawlist(panel,sprintf('drawfunctions(''drawpoint'',%d)',panel));
    end
    
    if do.scar
      addtodrawlist(panel,sprintf('drawfunctions(''drawviability'',%d)',panel));
    else
      do.manualinteraction = 0;
      %DATA.ViewIM{panel} = [];
    end
    
    if do.manualinteraction
      addtodrawlist(panel,sprintf('drawfunctions(''showviabilityedits'',%d)',panel));
    else
      %DATA.ViewIM{panel} = [];
    end
    
    if do.mar
      addtodrawlist(panel,sprintf('drawfunctions(''drawmar'',%d)',panel));
    end
    
  case {'montage','montagerow','montagesegmented'}
    
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'contour'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawcontours'',%d,''%s'')',panel,types{i}));
      end
    end
    
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'interp'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawinterp'',%d,''%s'')',panel,[types{i},'Interp']));
      end
    end
    
    if do.measures
      addtodrawlist(panel,sprintf('drawfunctions(''drawmeasures'',%d)',panel));
    end
    
    if do.text
      addtodrawlist(panel,sprintf('drawfunctions(''drawtext'',%d)',panel));
    end
    doroitext= false;
    for c = 'cgymkwrb'
      if do.([c,'roi'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawroi'',%d,''%s'')',panel,c));
        doroitext= true;
      end
    end
    if doroitext
      addtodrawlist(panel,sprintf('drawfunctions(''drawroitext'',%d)',panel));
    end
    
    if do.point
      addtodrawlist(panel,sprintf('drawfunctions(''drawpoint'',%d)',panel));
    end
    
    if do.scar
      addtodrawlist(panel,sprintf('drawfunctions(''drawviability'',%d)',panel));
    else
      do.manualinteraction = 0;
      DATA.ViewIM{panel} = [];
    end
    
    if do.manualinteraction
      addtodrawlist(panel,sprintf('drawfunctions(''showviabilityedits'',%d)',panel));
    else
      DATA.ViewIM{panel} = [];
    end
    
    if do.mar
      addtodrawlist(panel,sprintf('drawfunctions(''drawmar'',%d)',panel));
    end
  case {'hla','gla','vla'}
    if do.text
      addtodrawlist(panel,sprintf('drawfunctions(''drawtext'',%d)',panel));
    end
    
    if do.measures
      addtodrawlist(panel,sprintf('drawfunctions(''drawmeasures'',%d)',panel));
    end
    
    if do.intersections
      addtodrawlist(panel,sprintf('drawfunctions(''drawintersections'',%d)',panel));
    end
  case {'trans3DP','sag3DP','cor3DP','speedim'}
    
    types = {'gb','bg','gr','br','rb','rg'};
    for i = 1:length(types)
      if do.([types{i} 'line'])
        addtodrawlist(panel,sprintf('drawfunctions(''draw3dpline'',%d,''%s'')',panel,types{i}));
      end
    end
    
    if do.threedpcontour
      addtodrawlist(panel,sprintf('drawfunctions(''draw3dpoutline'',%d,''threedp'')',panel));
      addtodrawlist(panel,sprintf('drawfunctions(''draw3dpsurfaceoutline'',%d)',panel));      
    end
    
    if do.point
      addtodrawlist(panel,sprintf('drawfunctions(''drawpoint'',%d)',panel));
    end
    
end

%---------------------------------------------------------
function do = annotationinpanel(panel) 
%-----------------------------------------------------------
%this function sets graphics objects to nan in panel if there is no
%occurence in that slice of the annotation pertaining to that object.
%Returns struct with booleans indicating occurence of annotation.
global DATA SET

%Get global graphics handles
handles = DATA.Handles;

%this is done on panel level so if no segmentatio found in panel then it is
%not displayed in that panel;
slices = slicesinpanel(panel);
no = DATA.ViewPanels(panel);

%Assume all exists then turn off later in code. These fields are all
%lowercase and contours have the field name which matches the SET field
%name but lower case and contour after.
do.endocontour = 1;
do.epicontour = 1;
do.rvendocontour = 1;
do.rvepicontour = 1;
do.generalpencontour = 1;
do.lacontour = 1;
do.racontour = 1;
do.endointerp = 1;
do.epiinterp = 1;
do.rvendointerp = 1;
do.rvepiinterp = 1;
do.generalpeninterp = 1;
do.lainterp = 1;
do.rainterp = 1;
do.measures = 1;
do.point = 1;
do.point3D = 1;
do.croi = 1;
do.yroi = 1;
do.mroi = 1;
do.kroi = 1;
do.groi = 1;
do.rroi = 1;
do.broi = 1;
do.wroi = 1;
do.text = 1;
do.intersections = 1;
do.scar = 1;
do.mar = 1;
do.manualinteraction = 1;

%3dp objects
do.rgline = 1;
do.bgline = 1;
do.rbline = 1;
do.gbline = 1;
do.brline = 1;
do.grline = 1;
do.threedpcontour = 1;

%3D point is not needed by default
do.point3D = 0;

switch DATA.ViewPanelsType{panel}
  case {'hla','gla','vla'}
    %Somethings are never shown in hla gla vla
    do.endocontour = 0;
    do.epicontour = 0;
    do.rvendocontour = 0;
    do.rvepicontour = 0;
    do.generalpencontour = 0;
    do.lacontour = 0;
    do.racontour = 0;
    do.endointerp = 0;
    do.epiinterp = 0;
    do.rvendointerp = 0;
    do.rvepiinterp = 0;
    do.generalpeninterp = 0;
    do.lainterp = 0;
    do.rainterp = 0;
    do.measures = 1;
    do.point = 0;
    do.croi = 0;
    do.yroi = 0;
    do.mroi = 0;
    do.kroi = 0;
    do.groi = 0;
    do.rroi = 0;
    do.broi = 0;
    do.wroi = 0;
    do.text = 1;
    do.intersections = 1;
    do.scar = 0;
    do.mar = 0;
    do.manualinteraction = 0;
    
    %3dp objects
    do.rgline = 0;
    do.bgline = 0;
    do.rbline = 0;
    do.gbline = 0;
    do.brline = 0;
    do.grline = 0;
    do.threedpcontour = 0;
    %do.rcontour = 0;
    %do.bcontour = 0;
    %do.gcontour= 0;
    %do.scontour= 0;
    
    %measures-----------------------------------------------
    if ~isempty(SET(no).Measure)
      [measure,slice] = viewfunctions('getmeasurecoords',panel);
      measureinslice = zeros(1,length(SET(no).Measure));
      for loop=1:length(SET(no).Measure)
        ziv = round(measure(loop).Z);
        ziv = min(ziv):max(ziv);
        if ismember(slice,ziv)
          measureinslice(loop) = 1;
          break;
        end
      end
    else
      measureinslice =0;
    end
    
    if ~any(measureinslice)
      do.measures = 0;
      set([handles.measurement(panel) handles.measurementoutsideplane(panel)],...
        'XData', nan, 'YData', nan);
      set(handles.measurementtext(panel,:),'position',[nan nan]);
    end
    
  case {'trans3DP','sag3DP','cor3DP','speedim'}
    %Somethings are never shown in 3DP views
    do.endocontour = 0;
    do.epicontour = 0;
    do.rvendocontour = 0;
    do.rvepicontour = 0;
    do.generalpencontour = 0;
    do.lacontour = 0;
    do.racontour = 0;
    do.endointerp = 0;
    do.epiinterp = 0;
    do.rvendointerp = 0;
    do.rvepiinterp = 0;
    do.generalpeninterp = 0;
    do.lainterp = 0;
    do.rainterp = 0;
    do.measures = 0;
    do.point = 1;
    do.croi = 0;
    do.yroi = 0;
    do.mroi = 0;
    do.kroi = 0;
    do.groi = 0;
    do.rroi = 0;
    do.broi = 0;
    do.wroi = 0;
    do.text = 1;
    do.intersections = 0;
    do.scar = 0;
    do.mar = 0;
    do.manualinteraction = 0;
    
    paneltype = DATA.ViewPanelsType{panel};
    if isequal(paneltype,'speedim')
      switch SET(no).LevelSet.Pen.Color
        case 'r'
          paneltype = 'trans3DP';
        case 'g'
          paneltype = 'sag3DP';
        case 'b'
          paneltype = 'cor3DP';
      end
    end
    
    switch paneltype
      case 'trans3DP'
        do.rgline = 1;
        do.rbline = 1;
        do.grline = 0;
        do.gbline = 0;
        do.brline = 0;
        do.bgline = 0;        
      case 'sag3DP'
        do.rgline = 0;
        do.rbline = 0;
        do.grline = 1;
        do.gbline = 1;
        do.brline = 0;
        do.bgline = 0;        
      case 'cor3DP'
        do.rgline = 0;
        do.rbline = 0;
        do.grline = 0;
        do.gbline = 0;
        do.brline = 1;
        do.bgline = 1;        
    end
    
    do.threedpcontour = true; %outlineon*panels(currentind);
    
    types = {'rg','gr','rb','br','bg','gb'};
    for i =  1:length(types)
      if ~do.([types{i},'line'])
        set(handles.([types{i},'line'])(panel),'XData',nan,'YData',nan);
      end
    end
    
    %if ~do.([currenttype,'contour'])
    if ~do.threedpcontour
      handles.threedpcontour(panel).ZData = [127 0 ; 0 0];
    end
    
  case 'viewport'
    %the viewport panel has no annotations whatsoever.
    do.endocontour = 0;
    do.epicontour = 0;
    do.rvendocontour = 0;
    do.rvepicontour = 0;
    do.generalpencontour = 0;
    do.lacontour = 0;
    do.racontour = 0;
    do.endointerp = 0;
    do.epiinterp = 0;
    do.rvendointerp = 0;
    do.rvepiinterp = 0;
    do.generalpeninterp = 0;
    do.lainterp = 0;
    do.rainterp = 0;
    do.measures = 0;
    do.point = 0;
    do.croi = 0;
    do.yroi = 0;
    do.mroi = 0;
    do.kroi = 0;
    do.groi = 0;
    do.rroi = 0;
    do.broi = 0;
    do.wroi = 0;
    do.text = 0;
    do.intersections = 0;
    do.scar = 0;
    do.mar = 0;
    do.manualinteraction = 0;
    
    %3dp objects
    do.rgline = 0;
    do.bgline = 0;
    do.rbline = 0;
    do.gbline = 0;
    do.brline = 0;
    do.grline = 0;
    do.threedpcontour = 0;
    
  otherwise
    %3dp objects
    do.rgline = 0;
    do.bgline = 0;
    do.rbline = 0;
    do.gbline = 0;
    do.brline = 0;
    do.grline = 0;
    do.threedpcontour = 0;
    
    types = {'Endo','Epi','RVEndo','RVEpi'};
    %contours-------------------------------------------------
    for i = 1:length(types)
      ind = findfunctions('findframeswithsegmentation',types{i},no,slices);
      
      if ~any(ind)
        do.([lower(types{i}),'contour']) = 0;
        set(handles.([lower(types{i}),'contour'])(panel),'XData',nan,'YData',nan);
      end
    end
    
    types = {'Endo','Epi','RVEndo','RVEpi','LA','RA','GeneralPen'};
    %interpolationpoints--------------------------------------
    for i = 1:length(types)
      ind = findfunctions('findframeswithinterpolationpoints',types{i},no,slices);
      
      if ~any(ind)
        do.([lower(types{i}),'interp']) = 0;
        set(handles.([lower(types{i}),'interp'])(panel),'XData',nan,'YData',nan);
      end
    end

    %la, ra contours
    types = {'LA','RA'};
    for loop = 1:length(types)
      if isempty(SET(no).(types{loop}))
        set(handles.([lower(types{loop}),'contour'])(panel),'XData',nan,'YData',nan);
        set([handles.laavplane(panel), handles.laaxis(panel)],'Visible','off');
        set([handles.raavplane(panel), handles.laaxis(panel)],'Visible','off');
        do.([lower(types{loop}),'contour']) = 0;
      end
    end

    %general pen objects
    if isempty(SET(no).GeneralPenObjects)
      set(handles.generalpencontour(panel,:),'XData',nan,'YData',nan);
      set(handles.generalpentext(panel,:),'Position',[nan nan]);
      do.generalpencontour = 0;
    end    
    
    %rois----------------------------------------------------
    linkednos = SET(no).Linked;
    colors = [];
    for i = 1:length(linkednos)
      roisinslice = cellfun(@(x) any(~isempty(x)) && any(x==slices),{SET(linkednos(i)).Roi.Z});
      ls = [SET(linkednos(i)).Roi(roisinslice).LineSpec];
      ls = ls(1:2:end);
      colors = [colors,unique(ls)];

      %The roicurrent is set to nan if it isnt in the viewed slices

      if ~isempty(SET(linkednos(i)).RoiCurrent) && ~any(SET(linkednos(i)).RoiCurrent == find(roisinslice))...
          || isempty(SET(linkednos(i)).RoiCurrent) && isempty(find(roisinslice, 1))
        set(handles.roicurrent(panel),'XData', nan,'YData',nan);
        set(handles.roitext(panel,:),'Position',[nan nan]);
      end
    end
    
    colorlist = 'cgykmwrb';
    notincludedcolors = colorlist(~ismember(colorlist,colors));
    
    for c = notincludedcolors
      set(handles.([c,'roi'])(panel),'XData',nan,'YData',nan);
%       set(DATA.Handles.roitext(panel,:),'Position', [nan,nan]);
      do.([c,'roi']) = 0;
    end
  
    
    %points---------------------------------------------------
    if ~any(ismember(SET(no).Point.Z,slices))
      do.point = 0;
      set(handles.point(panel),'XData',nan,'YData',nan);
      set(handles.pointtext(panel,:),'position',[nan nan]);
    end
    
    %Scar------------------------------------------------------
    %For the contour objects the equivalent of nan rendering is placing a zero matix with a one in a corner.
%     visibleslices = viewfunctions('slicesinpanel', panel)
    if isempty(SET(no).Scar)
      tmp = [1,0;0,0];
      set([handles.scarcontour(panel),handles.weightedscarcontour(panel),...
           handles.moextentcontour(panel),handles.mocontour(panel)],...
           'ZData',tmp,'YData',tmp,'XData',tmp);
      do.scar = 0;
      do.manualinteraction = 0;
    elseif ~any(SET(no).Scar.Result(:,:,slices),'all')
      tmp = [1,0;0,0];
      set([handles.scarcontour(panel),handles.weightedscarcontour(panel),...
           handles.moextentcontour(panel),handles.mocontour(panel)],...
           'ZData',tmp,'YData',tmp,'XData',tmp);
      do.scar = 0;
      do.manualinteraction = 0;
    elseif ~any(SET(no).Scar.NoReflow(:,:,slices),'all')
      tmp = [1,0;0,0];
      set([handles.moextentcontour(panel),handles.mocontour(panel)],...
        'ZData',tmp,'YData',tmp,'XData',tmp);      
    end
    
    %Since manual interaction is incorporated into the panel
    %imagehandle we need to check the button state here
    stateandicon = iconson({'hidescarmanual','click3D'});
    if stateandicon{1} %state for hidescarmanual
      do.manualinteraction = 0;
    end

    if ~stateandicon{2} %state for click3D
      do.point3D = 0;
    end

    
    %Mar--------------------------------------------------------
    if isempty(SET(no).MaR)
      %no one can see it anyway
      tmp = [1,0;0,0];
      set(handles.marcontour(panel),'ZData',tmp,'YData',tmp,'XData',tmp);
      do.mar = 0;
      %JB 2019-11-06 commented out so users can draw mar even in time-resolved
      %data
%     elseif ~any(SET(no).MaR.Result(:,:,visibleslices),'all')
%       %no one can see it anyway
%       tmp = [1,0;0,0];
%       DATA.Handles.marcontour(panel).ZData = tmp;
%       do.mar = 0;
    end
    
    %measures---------------------------------------------------------------------
    if ~isempty(SET(no).Measure)
      [measure,slice] = viewfunctions('getmeasurecoords',panel);
      measureinslice = zeros(1,length(SET(no).Measure));
      for loop=1:length(SET(no).Measure)
        ziv = round(measure(loop).Z);
        ziv = min(ziv):max(ziv);
        if ismember(slice,ziv)
          measureinslice(loop) = 1;
          break;
        end
      end
    else
      measureinslice =0;
    end
    
    if ~any(measureinslice)
      do.measures = 0;
      set([handles.measurement(panel) handles.measurementoutsideplane(panel) handles.measurementtextline(panel,:)],...
        'XData',nan,'YData',nan);
      set(handles.measurementtext(panel,:),'position',[nan nan]);
    end
end

%------------------------------------
function saxpanels = getpanelswithsax
%------------------------------------
% Get all panels that contain SAX images
global DATA SET
lvnopanel = [];
rvnopanel = [];
saxpanels = [];
if ~isempty(DATA.LVNO)
  lvnopanel = find(DATA.ViewPanels ==DATA.LVNO);
end
if ~isempty(DATA.RVNO)
  rvnopanel = find(DATA.ViewPanels ==DATA.RVNO);
end
ind = nonzeros(DATA.ViewPanels);
viewplanes = {SET(ind).ImageViewPlane};
saxpanels = find(contains(viewplanes,'short-axis','IgnoreCase',true));
saxpanels = [lvnopanel,rvnopanel,saxpanels];

%---------------------
function updatevisibility(hideallexceptsax)
%-----------------------
%update the visibility of overlays in the images
global DATA
if nargin == 0
  hideallexceptsax = shouldhideallexceptsax;
end
%Get global variable
h = DATA.Handles;
handles = [h.marcontour h.scarcontour h.weightedscarcontour ...
  h.moextentcontour h.mocontour h.endocontourintersection ...
  h.epicontourintersection h.rvendocontourintersection ...
  h.rvepicontourintersection h.endointerp h.epiinterp h.rvendointerp ...
  h.rvepiinterp h.endocontour h.epicontour h.rvendocontour ...
  h.rvepicontour h.lacontour h.racontour, h.lainterp, h.rainterp, ...
  h.laavplane h.laaxis ...
  h.roicurrent h.croi h.rroi h.mroi h.kroi h.groi h.wroi ...
  h.yroi h.broi h.measurement h.point h.point3D h.text h.centercross];
  
set(handles,'Visible', 'on');
set(h.roitext(:),'Visible', 'on');
set(h.planeintersection,'Visible', 'on');
set([h.measurementtext h.measurementtextline h.pointtext],'Visible', 'on');
set([h.generalpencontour h.generalpeninterp],'Visible', 'on');
set(h.generalpentext(:),'Visible','on');

%The do structures contains information of which updates to use.
tmp = cell(1,length(DATA.ViewPanels));
for i  = 1:length(tmp)
  tmp{i} = 1;
end

do = struct('endocontour',tmp,'epicontour',tmp,'endointerp',tmp,'epiinterp',tmp,...
  'rvendocontour',tmp,'rvepicontour',tmp,'rvendointerp',tmp,'rvepiinterp',tmp,...
  'lacontour',tmp,'racontour',tmp,'lainterp',tmp,'rainterp',tmp,...
  'generalpencontour',tmp,'generalpeninterp',tmp,...
  'measures',tmp,'point',tmp,'point3D',tmp,'croi',tmp,'mroi',tmp,'groi',tmp,'rroi',tmp,...
  'broi',tmp,'wroi',tmp,'yroi',tmp,'kroi',tmp,'text',tmp,...
  'intersections',tmp,'scar',tmp,'mar',tmp,'manualinteraction',tmp,...
  'rgline',tmp,'bgline',tmp,'rbline',tmp,'gbline',tmp,'brline',tmp,...
  'grline',tmp,'threedpcontour',tmp); 

%Then we perform a check if there exists objects in any of the timeframes
%for contours interpolation points measures rois viability 
%this is done on panel level so if no segmentatio found in panel then it is
%not displayed in that panel;
for p = find(DATA.ViewPanels)
  do(p) = annotationinpanel(p);
end

%Then hide icons are incorporated. The hide state is set for all panels so
%if a state is hidden then that is applied for all panels.
hidecell={'hideplus','hidemar','hidescar',...
  'hideintersections','hideothercontour',...
  'hideinterp','hidelv','hiderv','hideroi',...
  'hidemeasure','hidepoint','hidetext','hidescarmanual',...
  'hidegeneralpenall','hidela','hidera','hidesegmentation'};

stateandicon = iconson(hidecell);
%If icon not found, then put its state to the same state as "Hide/show all
%overlays". The hide/show functions are activated only if the do structure agrees.
%This is useful in Segment 3DP for ex., where hide/show icons are not
%displayed in the main GUI.
state = iconson('hideall'); %get "Hide/show all overlays" state
for loop = 1:length(stateandicon)
  if ~isa(stateandicon{loop,2},'myicon') %then is nan
    stateandicon{loop,1} = state{1,1};
  end
end
state = [stateandicon{:,1}];

if hideallexceptsax
  saxpanels = getpanelswithsax;
  if ~isempty(saxpanels)
    nonsaxpanels = setdiff(1:numel(DATA.ViewPanels),saxpanels);
    if ~isempty(nonsaxpanels)
      state(5) = true; % contours' intersections
      state(7:8) = true; %LV/RV contours
      state(15:16) = true; %LA/RA contours
    end
  end
end

% First state corresponds to centercross
if state(1)
  set(h.centercross,'Visible', 'off');
end

%Second state corresponds to the marcontour.
if state(2)
  set(h.marcontour,'Visible', 'off');
  %remake images without manual interaction
  DATA.ViewIM = cell(size(DATA.ViewIM));
  [do.mar] = deal(0);
end

%Third state corresponds to the scarcontour.
if state(3)
  handles = [h.scarcontour h.weightedscarcontour h.moextentcontour h.mocontour];
  set(handles,'Visible', 'off');
  %remake images without manual interaction visually turned off at bottom
  %of function
  DATA.ViewIM = cell(size(DATA.ViewIM));
  [do.scar] = deal(0);
  [do.manualinteraction] = deal(0);
end

%Fourth state corresponds to the plane intersections.
if state(4)
  set(h.planeintersection,'Visible','off');
end

%Fifth state corresponds to the contour intersections.
if state(5)
  if ~hideallexceptsax
    handles = [h.endocontourintersection h.epicontourintersection ...
      h.rvendocontourintersection h.rvepicontourintersection];
    set(handles,'Visible', 'off');
  elseif hideallexceptsax && ~isempty(saxpanels)
    % hide contour intersections in SAX
    handles = [h.endocontourintersection(saxpanels) h.epicontourintersection(saxpanels) ...
      h.rvendocontourintersection(saxpanels) h.rvepicontourintersection(saxpanels)];
    set(handles,'Visible', 'off');
  end
  [do.intersections] = deal(0);
end

%Sixth state corresponds to the interpolationpoints.
if state(6)
  iph = mycat(2,h.endointerp, ...
    h.epiinterp, ...
    h.rvendointerp, ...
    h.rvepiinterp, ...
    h.lainterp, ...
    h.rainterp);
  set(iph,'Visible', 'off');
  set(h.generalpeninterp,'Visible', 'off');
  [do.interp] = deal(0);
end

%Seventh state corresponds to LV contours
if state(7)
  if hideallexceptsax && ~isempty(nonsaxpanels)
    lvh = mycat(2, h.endocontour(nonsaxpanels), ...
      h.epicontour(nonsaxpanels),...
      h.endointerp(nonsaxpanels), ...
      h.epiinterp(nonsaxpanels) ...
      );
  else
    lvh = mycat(2, h.endocontour, ...
      h.epicontour,...
      h.endointerp, ...
      h.epiinterp, ...
      h.endocontourintersection, ...
      h.epicontourintersection);
  end
  set(lvh,'Visible', 'off');
  [do.intersections] = deal(0);
end

%Eigth state corresponds to RV contours
if state(8)
  if hideallexceptsax && ~isempty(nonsaxpanels)
        rvh = mycat(2, h.rvendocontour(nonsaxpanels), ...
      h.rvepicontour(nonsaxpanels),...
      h.rvendointerp(nonsaxpanels), ...
      h.rvepiinterp(nonsaxpanels) ...
      );
  else
    rvh = mycat(2, h.rvendocontour, ...
      h.rvepicontour,...
      h.rvendointerp, ...
      h.rvepiinterp, ...
      h.rvendocontourintersection, ...
      h.rvepicontourintersection);
  end
  set(rvh,'Visible', 'off');
  [do.intersections] = deal(0);
end

%nineth state corresponds to Rois
if state(9)
  set(h.roicurrent,'Visible','off');
  set([h.croi,...
    h.rroi,...
    h.mroi,...
    h.kroi,...
    h.groi,...
    h.wroi,...
    h.yroi,...
    h.broi],...
    'visible','off')
  set(h.roitext(:),'Visible','off');
  [do.roi] = deal(0);
end

%Tenth state corresponds to measures
if state(10)
  set(h.measurement,'Visible','off');
  set([h.measurementtext, h.measurementtextline],'Visible','off');
  [do.measures] = deal(0);
end

%Eleventh state corresponds to points
if state(11)
  set(h.point,'Visible','off');
  set(h.pointtext,'Visible','off');
  [do.point] = deal(0);
end

%Twelveth state corresponds to text
if state(12)
  set(h.text,'Visible','off');
  set(h.pointtext,'Visible','off');
  set([h.measurementtext, h.measurementtextline],'Visible','off');
  set(h.roitext(:),'Visible','off');
  set(h.generalpentext(:),'Visible','off');
  [do.text] = deal(0);
end

% %Thirteenth state corresponds to extent
% if state(13)
%     set(DATA.Handles.moextentcontour,'Visible', 'off')
% end

%fourteenth state corresponds to manual interaction
if state(13)
  DATA.ViewIM = cell(size(DATA.ViewIM));
  [do.manualinteraction] = deal(0);
end

%this state corresponds to general pen contours and interp
if state(14)
  generalpenhandles = [h.generalpencontour, h.generalpeninterp];
  set(generalpenhandles,'Visible', 'off');
  set(h.generalpentext(:),'Visible','off');
  [do.generalpencontour, do.generalpeninterp] = deal(0);
end

%this state corresponds to LA contours and interp
if state(15)
  lahandles = [h.lacontour, h.lainterp, h.laaxis, h.laavplane];
  set(lahandles,'Visible', 'off');
  [do.lacontour, do.lainterp] = deal(0);
end

%this state corresponds to RA contours and interp
if state(16)
  rahandles = [h.racontour, h.rainterp, h.raavplane];
  set(rahandles,'Visible', 'off');
  [do.racontour, do.rainterp] = deal(0);
end

% add to drawlist what updates need to be made.
DATA.drawlist = cell(1,length(DATA.ViewPanels));

%Drawlist is a cell of cells
for p = 1:length(DATA.ViewPanels)
  DATA.drawlist{p} = {};
end

%this only updates the draw list
for p = find(DATA.ViewPanels)
  updatedrawlist_helper(p,do(p))
end

%this resets images where we've hidden on off manual interaction
for p = find(DATA.ViewPanels)%p = find(cellfun(@isempty,DATA.ViewIM))
  drawfunctions('drawpanel',p)%drawfunctions('drawpanel',p)
end

%------------------------------------------------------------------------
function selectgeneralpen(panel) 
%------------------------------------------------------------------------
%Function that sets the current General Pen object and updates the graphics
global DATA SET

no = DATA.ViewPanels(panel);

[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,x,y);
tf = SET(no).CurrentTimeFrame;

[~, objind] = findfunctions('closestgeneralpen',no,x,y,slice,tf);

if ~isempty(objind) && objind > 0
  DATA.GeneralPenSettings.setcurrentobject(objind);
  drawfunctions('drawcontoursgeneralpen',panel);
  generalpen.generalpenfunctions('indentobjecticon',objind);
end

%------------------------------------------------------------------------
function selectroi(panel)
%------------------------------------------------------------------------
%Function that sets currentroi and updates the graphics
global DATA SET
no = DATA.ViewPanels(panel);
if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  magno = SET(no).Flow.MagnitudeNo;
else
  magno = no;
end
[y,x] = mygetcurrentpoint(DATA.Handles.imageaxes(panel));
slice = viewfunctions('clickedslice',panel,x,y);
slicestoinclude = slicesinpanel(panel);
scale = getscale(panel);

%only work with rois that exist in the included slices
roistodo = find(cellfun(@(x) ~isempty(x) && any(x==slicestoinclude),{SET(magno).Roi.Z}));
dist2 = zeros(1,length(roistodo));

%get clicked position

for loop=1:length(roistodo)
  [xl,yl] = ind2sub(DATA.ViewPanelsMatrix{panel},find(slicestoinclude==SET(magno).Roi(roistodo(loop)).Z));
  yl = (yl-1)*SET(no).YSize;
  xl = (xl-1)*SET(no).XSize;
  dist2(loop) = min((scale*(xl+SET(magno).Roi(roistodo(loop)).Y(:,SET(no).CurrentTimeFrame))-y).^2 + ...
    (scale*(yl+SET(magno).Roi(roistodo(loop)).X(:,SET(no).CurrentTimeFrame))-x).^2);
end

[val,ind] = min(dist2);
if sqrt(val)< DATA.Pref.ContourAdjustDistance%3
  currentroi = roistodo(ind);
else
  currentroi = [];
end

if ~isempty(currentroi)
  SET(magno).RoiCurrent = currentroi;
  for p = find(ismember(DATA.ViewPanels,SET(no).Linked))
    drawfunctions('drawroi',p);
  end
  segment('updateflow')
end

%-----------------------------------
function setviewtype(viewpanelstype)
%-----------------------------------
%Changes viewtype for current
global DATA NO SET

if strcmp(viewpanelstype,'montage') && SET(NO).ZSize>100
  if not(yesno('Image stack consist of many slices, this may take some seconds to display. Continue to montage view anyway?'))
    return;
  end
end

if DATA.ViewPanels(DATA.CurrentPanel) ==  0
  return
end

%If we are in orthoview mode we need to reset to single panel
if any(strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},{'orth','hla','gla','vla'}))
  viewfunctions('setview',1,1,NO,{'one'})
  return
end

DATA.ViewPanelsType{DATA.CurrentPanel} = viewpanelstype;

if (strcmp(DATA.ViewPanelsType{DATA.CurrentPanel},'one'))&&(isempty(SET(NO).CurrentSlice))
  SET(NO).CurrentSlice = 1;
  SET(NO).EndSlice = 1;
  SET(NO).StartSlice = 1;
else
  SET(NO).CurrentSlice = min(SET(NO).EndSlice,max(SET(NO).StartSlice,SET(NO).CurrentSlice));
end

if ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'montage','montagerow','montagesegmented'})
  if ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'montagesegmented'})
    slicestoshow = slicesinpanel(DATA.CurrentPanel);
    [mrows,mcols] = calcfunctions('calcrowscols',DATA.ViewPanels(DATA.CurrentPanel),length(slicestoshow));
  else
    [mrows,mcols] = calcfunctions('calcrowscols',DATA.ViewPanels(DATA.CurrentPanel));
  end
  if ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'montagerow'})
    m = mrows*mcols;
    mrows = min(mrows,2);
    mcols = ceil(m/mrows);
  end
  DATA.ViewPanelsMatrix{DATA.CurrentPanel} = [mcols mrows];
else
  DATA.ViewPanelsMatrix{DATA.CurrentPanel} = [1,1];
end

%Assert correct button is indented
DATA.setviewbuttons(0)

%this draws and updates the panel image intersections.
drawfunctions('drawplaneintersections')
drawfunctions('drawpoint3D')

%clears current viewim
DATA.ViewIM{DATA.CurrentPanel}=[];

%get the drawlist and draw panels.
updatedrawlist(DATA.CurrentPanel)
drawfunctions('drawpanel',DATA.CurrentPanel)

drawfunctions('drawselectedframe',DATA.CurrentPanel)
drawfunctions('drawselectedslice',DATA.CurrentPanel)
updatetextposition(DATA.CurrentPanel)

%---------------------
function addno2panel(panel,no,viewpanelstype,checklinked) 
%-----------------------
%Adds a no to a panel. The viewtype is one by default

global DATA SET NO

%profile on
if nargin < 4
  checklinked = true;
end
if nargin <3
  viewpanelstype = 'one';
end
if nargin < 2
  no = NO;
end

%connect interpolation in the stack that was previously selected
tools('connectallinterpolation',NO);

%Apply all settings
DATA.CurrentPanel = panel;
previousno = DATA.ViewPanels(panel);
for num = 1:numel(DATA.ViewPanels)
  DATA.Handles.point3D(num).XData(:) = nan;
  DATA.Handles.point3D(num).YData(:) = nan;
end
if previousno > 0 && length(SET(previousno).Linked) == 2 && checklinked
  % this was a stack with linked stacks
  % so check if the linked stack were also visible
  linkedpanelsind = find(ismember(DATA.ViewPanels,setdiff(SET(previousno).Linked,previousno)));
  if ~isempty(linkedpanelsind) && length(SET(no).Linked) == 2
    % look if the new no is also linked
    linkedstack = setdiff(SET(no).Linked,no);
    if numel(linkedpanelsind) > 1
      paneldist = linkedpanelsind - panel;
      paneldist(paneldist<0) = Inf;
      linkedpanelsind = linkedpanelsind(paneldist == min(paneldist));
    end
    % set linked panel first to ensure that the clicked one is set to NO
    addno2panel(linkedpanelsind(1),linkedstack,viewpanelstype, false)
  end
end

previouspaneltype = DATA.ViewPanelsType{panel};

%Delete volshow if visible
if isequal(previouspaneltype,'viewport')
  angioview('delete_Callback')
end

DATA.ViewPanels(panel) = no;
DATA.ViewPanelsType{panel} = viewpanelstype;

if ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'montage','montagerow'})
  [mrows,mcols] = calcfunctions('calcrowscols',DATA.ViewPanels(DATA.CurrentPanel));
  if ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'montagerow'})
    m = mrows*mcols;
    mrows = min(mrows,2);
    mcols = ceil(m/mrows);
  end
  DATA.ViewPanelsMatrix{DATA.CurrentPanel} = [mcols mrows];
else
  DATA.ViewPanelsMatrix{DATA.CurrentPanel} = [1,1];
end

DATA.ViewIM{panel} = [];
NO = no;

%get the drawlist and draw panels.
updatedrawlist(panel)

%update altered panels fully
panelslinked = find(ismember(DATA.ViewPanels,SET(no).Linked));
for p = panelslinked
  drawfunctions('drawpanel',p)
end

%this draws and updates the panel image intersections.
drawfunctions('drawplaneintersections')
drawfunctions('drawintersections')  %empty input udpate all viewpanels

drawfunctions('drawselectedframe',panel)
drawfunctions('drawselectedslice',panel)
drawfunctions('drawthumbnails',isempty(DATA.DATASETPREVIEW));
%drawfunctions('drawthumbnailframes')
createfunctions('addcolorbar',panel)
updatetextposition(panel)
hideallexceptsax = shouldhideallexceptsax;
if hideallexceptsax
  updatevisibility(hideallexceptsax)
end
%Assert correct button is indented
DATA.setviewbuttons(0)
if strcmp(DATA.ProgramName,'Segment') 
  segment('updatevolume');
elseif contains(DATA.ProgramName,'CMR') || contains(DATA.CurrentTheme,'lv')
  % update volume axes
  DATA.updatevolumeaxes;
end
segment('updateflow');
segment('updatemeasurement');
DATA.updatetimebaraxes

%Update toolbar if any general pen objects exists in the new stack
generalpen.generalpenfunctions('updateobjectstoolbar');

%-------------------------------------------
function switchimagestack(no,viewpanelstype) 
%-------------------------------------------
%Toggles to the first panel containing imagestack no. In the code it finds
%the first panel containing no and then uses switchpanel to switch to it

global DATA

panel = find(DATA.ViewPanels==no,1);

%if the imagestack is not out then add it to the first place
if isempty(panel)
  switchpanel(panel);
end

if nargin >1
  if ~strcmp(viewpanelstype,DATA.ViewPanelsType{panel})
    setviewtype(viewpanelstype)
  end
end

%--------------------------------------------------------
function zoomstate = getnewzoomstate(panel,no, resetflag)
%--------------------------------------------------------
%Compute initialization of zoom state. Takes panel, no, resetflag

global DATA SET

if nargin<3
  resetflag = true; %default is true
end
if nargin<2
  no = DATA.ViewPanels(panel);
end

switch DATA.ViewPanelsType{panel}
  case 'trans3DP'
    color = 'r';
    [xsize,ysize,xres,yres] = segment3dp.tools('viewsizeandres',color,no);
  case 'sag3DP'
    color = 'g';
    [xsize,ysize,xres,yres] = segment3dp.tools('viewsizeandres',color,no);
  case 'cor3DP'
    color = 'b';
    [xsize,ysize,xres,yres] = segment3dp.tools('viewsizeandres',color,no);
  case 'speedim'
    if isempty(SET(no).LevelSet) %~isfield(SET(no).LevelSet.Pen,'Color')
       color= 'r';
    else
       color = SET(no).LevelSet.Pen.Color;
    end
    [xsize,ysize,xres,yres] = segment3dp.tools('viewsizeandres',color,no);
  otherwise
  xsize = SET(no).XSize;
  ysize = SET(no).YSize;
  xres = SET(no).ResolutionX;
  yres = SET(no).ResolutionY;
end

%Find panelsize
panelpos = DATA.Handles.imageaxes(panel).Position;
figpos = DATA.fig.Position;
panelwidth = panelpos(3)*figpos(3)/(ysize*yres);
panelheight = panelpos(4)*figpos(4)/(xsize*xres);

%maxsize = max(SET(no).XSize,SET(no).YSize);
minpanelsize = min(panelwidth,panelheight);
panelwidth = panelwidth/minpanelsize; %range 1..
panelheight = panelheight/minpanelsize; %range 1..

if isnan(panelwidth) %Happens when not panel is not setup yet (i.e position and size is zeros)
  panelwidth = 1;
  panelheight = 1;
end

xcenter = xsize/2;
ycenter = ysize/2;

zoomstate = [...
  ycenter-panelwidth*ysize/2 ...
  ycenter+panelwidth*ysize/2 ...
  xcenter-panelheight*xsize/2 ...
  xcenter+panelheight*xsize/2]+0.5;

if resetflag  
    switch DATA.ViewPanelsType{panel}
    case 'trans3DP'
      SET(no).LevelSet.View.RZoomState = zoomstate;
    case 'sag3DP'
      SET(no).LevelSet.View.GZoomState = zoomstate;
    case 'cor3DP'
      SET(no).LevelSet.View.BZoomState = zoomstate;
    case 'speedim'
      %nothing
    case 'one'
      SET(no).NormalZoomState = zoomstate;
    end
end

%------------------
function zoom(incr,panel) 
%------------------
%New solution for zooming

global DATA SET

%Do not use this tool on 3D view
if strcmp(DATA.ViewPanelsType{1},'viewport')
  str = [dprintf('This feature is not available for 3D view.') newline];
  str = [str dprintf('Please use the mouse wheel instead.')];
  mywarning(str);
  return
end

if nargin < 2
  %Get panel
  panel = DATA.CurrentPanel;
end

if ~isempty(DATA.LastKey) && strcmp(DATA.LastKey,'shift')
    % all visible NOs
    nos = DATA.ViewPanels(DATA.ViewPanels > 0);
    linkedpanelsindx = find(DATA.ViewPanels > 0);
    noslinked = 0;
else
    %Get no
    no = DATA.ViewPanels(panel);
    noslinked = SET(no).Linked;
    %Get linked panels
    [nos,~,linkedpanelsindx] = intersect(noslinked,DATA.ViewPanels);
end
 
%loop over linked stacks
for loop = 1:length(nos)
  
  %extract current stack
  no = nos(loop);
  panel = linkedpanelsindx(loop);

  zoomstate = SET(no).NormalZoomState; %get zoomstate
  
  if isempty(zoomstate)
    zoomstate = getnewzoomstate(panel,no);
  end

  %Here occurs the zoom!
  
  %Compute x&y size and center
  xsize = zoomstate(2)-zoomstate(1);
  ysize = zoomstate(4)-zoomstate(3);
  xcenter = sum(zoomstate(1:2))/2;
  ycenter = sum(zoomstate(3:4))/2;
  
  %Compute factor
  f = 1-incr/5;
  
  %Compute newzoom state based on the center and how much to extend
  newzoomstate = [xcenter-xsize/2*f xcenter+xsize/2*f ycenter-ysize/2*f ycenter+ysize/2*f];

  %assign
  if all(noslinked ~= 0) && ~all(nos == noslinked)
    % set also the invisble nos to the same Zoomstate
    ind = nos ~= noslinked; %find index of nos to update
    nos2update = noslinked(ind); %find nos to update
    for nosind = 1:length(nos2update) %update linked nos
      SET(nos2update(nosind)).NormalZoomState = newzoomstate;
    end
  end
  SET(no).NormalZoomState = newzoomstate;
  if contains(DATA.ViewPanelsType{DATA.CurrentPanel},'montage')
    %clears current viewim
    DATA.ViewIM{panel}=[];
    drawfunctions('drawpanel',panel);
    drawfunctions('drawselectedslice',panel);
  else
    %update it graphically
    updatezoomandaspectratio(panel);
    drawfunctions('drawplaneintersections');
    %we want to update the text position and the frame
    drawfunctions('drawselectedframe',panel)
    updatetextposition(panel)
  end
end

%---------------------
function updatetextposition(panel)
%-----------------------
%This function updates the text position after updates such as zooming or
%panning

global DATA

%Get global graphics handles
handles = DATA.Handles;
imageaxes = handles.imageaxes;
xlim = imageaxes(panel).XLim;
ylim = imageaxes(panel).YLim;

%Get position
x = max(0.5,0.015*(xlim(2) - xlim(1)));
y = max(0.5,0.015*(ylim(2) - ylim(1)));

%Set position
set(handles.text(panel),'Position',[xlim(1) + x, ylim(1) + y]);

%--------------------------
function switchpanel(panel)
%--------------------------
%Toggles the CurrentPanel and NO to be correct aswell as does the necessary
%graphical updates. If trying to toggle to an empty panel dont change
%anything

global DATA NO

if isempty(panel)
  return
end

%connect interpolation in the stack that was previously selected
tools('connectallinterpolation',NO);

no = DATA.ViewPanels(panel);

%this is the empty panel case
if ~isequal(no,0)
  DATA.CurrentPanel = panel;
  NO = no;
  drawfunctions('drawselectedframe',panel)
  drawfunctions('drawthumbnailframes')
end

if ~DATA.issegment3dp
  %this draws and updates the panel image intersections.
  drawfunctions('drawplaneintersections')
else
  segment3dp.viewfunctions('switchpanel3dp',panel);
  return
end

DATA.setviewbuttons(0);
if strcmp(DATA.ProgramName,'Segment')
  segment('updatevolume');  
end
segment('updateflow'); 
segment('updatemeasurement');
DATA.updatetimebaraxes;
DATA.CurrentPanel = panel;

%Update toolbar is any general pen objects exists in the new stack
generalpen.generalpenfunctions('updateobjectstoolbar');

%----------------------------------------------
function placecontextmenu(panel,type,objectind) 
%----------------------------------------------
%Function that returns the clicked slice in montage and montage
global DATA SET NO

[p(1),p(2)] = mygetcurrentpoint(DATA.fig);

switch type
  case 'Measure'
    set(DATA.Handles.measurecontextmenu,'Position',p,'Visible','on');
  case 'Point'
    set(DATA.Handles.pointcontextmenu,'Position',p,'Visible','on');
  case {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'}
    set(DATA.Handles.interppointmenu,'Position',p,'Visible','on');
  case {'LAInterp','RAInterp'}
    set(DATA.Handles.larainterppointmenu,'Position',p,'Visible','on');
  case 'Roi'
    % set current roi to the clicked object in (objectind)
    SET(NO).RoiCurrent = objectind;
    drawfunctions('drawroi',panel);
    set(DATA.Handles.roicontextmenu,'Position',p,'Visible','on');
%   case 'GeneralPen'
%     % set current object to the clicked object
%     generalpen.generalpenfunctions('indentobjecticon',objectind);
%     DATA.GeneralPenSettings.setcurrentobject(objectind);
%     drawfunctions('drawcontoursgeneralpen',panel);
%     set(DATA.Handles.generalpencontextmenu,'Position',p,'Visible','on');
%   case 'GeneralPenInterp'
%     % set current object to the clicked object
%     generalpen.generalpenfunctions('indentobjecticon',objectind);
%     DATA.GeneralPenSettings.setcurrentobject(objectind);
%     drawfunctions('drawcontoursgeneralpen',panel);
%     set(DATA.Handles.interppointmenu,'Position',p,'Visible','on');
  case {'RVEndo','RVEpi'}
    set(DATA.Handles.rvcontextmenu,'Position',p,'Visible','on');
  case {'Endo','Epi'}
    set(DATA.Handles.lvcontextmenu,'Position',p,'Visible','on');
  case 'LA'
    set(DATA.Handles.lacontextmenu,'Position',p,'Visible','on');
  case 'RA'
    set(DATA.Handles.racontextmenu,'Position',p,'Visible','on');
  case {'Scar','ScarRubber'}
    set(DATA.Handles.scarcontextmenu,'Position',p,'Visible','on');
  case {'Mar','MaRRubber'}
    set(DATA.Handles.marcontextmenu,'Position',p,'Visible','on');

  otherwise
%     switch DATA.CurrentTool
%       case {...'Endo','Epi',...
%           ...'EndoInterp','EpiInterp','RVEndo','RVEpi',...
%           'Scar','MO','MaR','ScarRubber','MORubber','MaRRubber',...
%           ...'RVEndoInterp','RVEpiInterp',...
%           'EndoBalloon','Contrast'}        
%         set(DATA.Handles.lvcontextmenu,'Position',p,'Visible','on');
%       otherwise
        set(DATA.Handles.selectcontextmenu,'Position',p,'Visible','on');
%     end
end

%------------------------------------------------------------------------
function slice = clickedslice(panel,x,y) 
%------------------------------------------------------------------------
%Function that returns the clicked slice in montage and montage
global DATA SET
no = DATA.ViewPanels(panel);

if any(strcmp(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','speedim','cor3DP'}))
  switch SET(no).LevelSet.Pen.Color
    case 'r'
      slice = SET(no).LevelSet.View.RSlice;
    case 'g'
      slice = SET(no).LevelSet.View.GSlice;
    case 'b'
      slice = SET(no).LevelSet.View.BSlice;
  end
else
  scale = getscale(panel);
  imdim = zoomfunctions.getxysize(no,panel);
  x_ind = ceil(x/scale/imdim.YSize);
  y_ind = ceil(y/scale/imdim.XSize);
  
  x_ind = max(x_ind,1);
  x_ind = min(x_ind,DATA.ViewPanelsMatrix{panel}(1));
  
  y_ind = max(y_ind,0);
  y_ind = min(y_ind,DATA.ViewPanelsMatrix{panel}(2));
  
  ind = sub2ind(DATA.ViewPanelsMatrix{panel},x_ind,y_ind);
  
  slicestoinclude = slicesinpanel(panel);
  
  %clicked in black area of montage
  if ind>length(slicestoinclude)
    slice = [];
    return
  end
  slice = slicestoinclude(ind);
end

%------------------------------------------------------------
function setview(rows,cols,nos,viewpanelstype,currentnopanel)
%------------------------------------------------------------
%Creates imageaxes and boxaxes. The field boxaxes is probably unnecessary.

global DATA NO SET

if DATA.Silent || DATA.Autoloader
  return
end

h = mywaitbar([],DATA.fig); %Creates an non visible waitbar
close(h); %This clears broken waitbars

%If Segment 3DPrint then run special code
if DATA.issegment3dp
  segment3dp.viewfunctions('setview3dp');
  return
end

%selectedpaneltype = DATA.ViewPanelsType{DATA.CurrentPanel};

%Set information fields
%profile on
%this is where the refresh option finds which viewpanelbutton is indented and
%sets rows and cols
if nargin == 0 || (isempty(rows) && isempty(cols))
  rows = DATA.ViewMatrix(1);
  cols = DATA.ViewMatrix(2);
 
  panelnos = DATA.ViewPanels(DATA.ViewPanels>0);
  for loop = panelnos
    SET(loop).NormalZoomState=[];
  end

  DATA.fig.WindowButtonMotionFcn = 'segment(''toggleplaceholdermotion'')';
  DATA.fig.WindowButtonUpFcn = '';  
else
  DATA.ViewMatrix(1) = rows;
  DATA.ViewMatrix(2) = cols;
end

%Find panels to clean
paneltoclean = find(DATA.ViewPanels>0);

for panelnum = paneltoclean
  colorbar(DATA.Handles.imageaxes(panelnum),'off')
  
  if isequal(DATA.ViewPanelsType{panelnum},'viewport')
    angioview('delete_Callback')
    DATA.ViewPanelsType{panelnum} = 'one';
  end
  
  emptydata = [1,0; 0 0];
  set([DATA.Handles.scarcontour(panelnum), ...
    DATA.Handles.weightedscarcontour(panelnum),...
    DATA.Handles.marcontour(panelnum),...
    DATA.Handles.moextentcontour(panelnum), ...
    DATA.Handles.mocontour(panelnum)],...
    'ZData',emptydata,'YData',emptydata,'XData',emptydata);
  set(DATA.Handles.selectedframe(panelnum),'YData',[],'XData',[]);
end
specialview = false;
%We need to handle different input cases and the chance that some panels
%are empty. also the current panel might be out of bounds when toggling to
%a lesser amount of panels. Then we set the CurrentPanel and NO field to
%the first panel content using the switchpanel command
if nargin < 3
  %This is the case when we toggle to a panel view from orthoview
  if any(ismember(DATA.ViewPanelsType(:),{'orth','hla','gla','vla'})) && nargin~=0 %if refresh keep ortho view
    DATA.ViewPanels=zeros(1,rows*cols);
    DATA.ViewPanels(1)=NO;
    DATA.ViewPanelsType{1}='one';
    DATA.ViewPanelsType(2:end)=[];
    switchpanel(1);
  else %if we are going to a panelview with less panels than now just remove the last entries in DATA.ViewPanels NO and Current Panel will be set as first later
    if length(DATA.ViewPanels)>rows*cols
      DATA.ViewPanels = DATA.ViewPanels(1:rows*cols);
    else %if we are going to a panelview with more panels than now just add zeros in ViewPanels
      tmp = zeros(1,rows*cols);
      tmp(1:length(DATA.ViewPanels)) = DATA.ViewPanels;
      DATA.ViewPanels = tmp;
    end
  end  
  DATA.LastView = [];
else  
  if isempty(nos)
    DATA.ViewPanels = NO*ones(1,rows*cols);
  else
  	specialview = true;
    DATA.ViewPanels = nos;
    if isempty(intersect(NO,nos))
      NO = nos(1);
    end
  end
end

%this takes care of if the currentpanel field is out of bounds and places
%it within the possible integers also sets the current
if nargin > 4
  %input give the panel that should be current (NO)
  NO = currentnopanel(1);
  if ~isempty(viewpanelstype)
    % update viewpaneltype, since it is used in switchpanel function
    DATA.ViewPanelsType = viewpanelstype;
  end
  switchpanel(currentnopanel(2));
elseif DATA.CurrentPanel > rows*cols
  DATA.ViewPanels(rows*cols) = NO;
  switchpanel(rows*cols);
end

%if NO isnt among the current displayed viewpanels it is set as the stack
%in the first viewpanel
if isempty(intersect(NO,DATA.ViewPanels))
  NO = DATA.ViewPanels(1);
end

%This asserst that viewpanelstype is set even if we dont input a cell
%specifying. If there exists prior viewpanelstype it uses it otherwise set all as one.
if nargin<4
  if isempty(DATA.ViewPanelsType)
    DATA.ViewPanelsType = cell(1,length(DATA.ViewPanels));
    DATA.ViewPanelsType(:) = {'one'};
  elseif length(DATA.ViewPanelsType(:))<rows*cols
    tmp = cell(1,rows*cols);
    for panelnum = 1:length(DATA.ViewPanelsType(:))
      tmp{panelnum} = DATA.ViewPanelsType{panelnum};
    end
    DATA.ViewPanelsType = tmp;
  else
    DATA.ViewPanelsType = DATA.ViewPanelsType(1:rows*cols);
  end
else
  DATA.ViewPanelsType = viewpanelstype;
end
for typeloop = 1:length(DATA.ViewPanelsType) 
  emptytype(typeloop) = isempty(DATA.ViewPanelsType{typeloop}); %#ok<AGROW> 
end
if isempty(DATA.ViewPanelsType) || any(emptytype)
  ...DATA.ViewPanelsType = cell(1,length(DATA.ViewPanels));
  DATA.ViewPanelsType(emptytype) = {'one'};
end

% Here we set so that the viewpanelsmatrix field is correct it is [1,1] if
% not montage
DATA.ViewPanelsMatrix = cell(1,length(DATA.ViewPanels));
for panelnum = 1:numel(DATA.ViewPanels)
  if DATA.ViewPanels(panelnum) > 0 && any(strcmp(DATA.ViewPanelsType{panelnum},{'montage','montagerow','montagesegmented'}))
    if strcmp(DATA.ViewPanelsType{panelnum},'montagesegmented')
      slicestoinclude = segment('getmontagesegmentedslices',DATA.ViewPanels(panelnum));
      [mrows,mcols] = calcfunctions('calcrowscols',DATA.ViewPanels(panelnum),length(slicestoinclude));
    else
      [mrows,mcols] = calcfunctions('calcrowscols',DATA.ViewPanels(panelnum));
    end
    if ismember(DATA.ViewPanelsType{panelnum},{'montagerow'})
      m = mrows*mcols;
      mrows = min(mrows,2);
      mcols = ceil(m/mrows);
    end
    DATA.ViewPanelsMatrix{panelnum} = [mcols mrows];
  else
    DATA.ViewPanelsMatrix{panelnum} = [1,1];
  end
end

specialview = true;

if ~specialview
    %functionality to restore same zoom state as before panel settings were
  %changed
  panelstodo = find(DATA.ViewPanels(DATA.ViewPanels>0));
  incr = zeros(length(panelstodo),1);

  for actpanel = panelstodo
    no = DATA.ViewPanels(actpanel);
    if no>length(SET)
      linkednos = [];
    else
      linkednos = SET(no).Linked;
    end
    numnos = length(linkednos);
    if numnos > 1
      for cl = 2:numnos
        SET(linkednos(cl)).NormalZoomState = [];
      end
    end
    type = DATA.ViewPanelsType{actpanel};
    if any(contains({'one'},type))
      zoomstatecurrent = SET(no).NormalZoomState;
      if ~isempty(zoomstatecurrent)
        % compare current state to the original
        % where original would be zoom state without zooming
        zoomstateoriginal = getnewzoomstate(actpanel,no,false);
        szxorig = zoomstateoriginal(2)-zoomstateoriginal(1);
        szxcurrent = zoomstatecurrent(2)-zoomstatecurrent(1);
        if ~isnan(szxorig) && ~isnan(szxcurrent) && ~isapproxequal(szxorig,szxcurrent)
          % ratios not equal -> ther was zooming
          f = szxcurrent/szxorig;
          % recalculate zoom increment
          incr(actpanel) = round((1-f)*7);    
        end
      end
    end
  end  
end

%This sets the preallocated axes to fit the above settings in ViewPanels
%etc...
configaxes
if any(ismember(DATA.ViewPanelsType,{'orth'}))
  DATA.CurrentPanel = 1;
end
  
%this draws and updates the panel image intersections.
drawfunctions('drawplaneintersections')

%When updating we only want to update the panels with content. These are
%DATA.ViewPanels>0
nostodo = DATA.ViewPanels(DATA.ViewPanels>0);
panelstodo = find(DATA.ViewPanels>0);

for panelnum = panelstodo
  
  updatedrawlist(panelnum)
  DATA.ViewIM{panelnum} = [];
  shouldzoom = false;
  if ismember(DATA.ViewPanelsType{panelnum},{'montage','montagerow','montagesegmented'})
    no = nostodo(panelnum);
    zoomstatecurrent = SET(no).NormalZoomState;
    SET(no).NormalZoomState = [];
    if ~isempty(zoomstatecurrent)
      % compare current state to the original
      % where original would be zoom state without zooming
      zoomstateoriginal = getnewzoomstate(panelnum,no,false);
      szxorig = zoomstateoriginal(2)-zoomstateoriginal(1);
      szxcurrent = zoomstatecurrent(2)-zoomstatecurrent(1);
      if ~isnan(szxorig) && ~isnan(szxcurrent) && ~isapproxequal(szxorig,szxcurrent)
        % ratios not equal -> ther was zooming
        f = szxcurrent/szxorig;
        % recalculate zoom increment
        incr(panelnum) = round((1-f)*7);
        shouldzoom = true;
      end
    end
  end
  
  drawfunctions('drawpanel',panelnum)
  if shouldzoom
    inc = incr(panelnum);
    if inc > 0
      zoomfactor = 1;
    else
      zoomfactor = -1;
    end
    for n = 1:abs(inc)
      zoom(zoomfactor,panelnum);
    end
  end

  if ~specialview
    if incr(panelnum)~= 0
        % restore zoom
        SET(nostodo(panelnum)).NormalZoomState = [];
        zoom(incr(panelnum),panelnum);
    end
  end
  if ~strcmp(DATA.ViewPanelsType{panelnum},'speedim') ...
      && ~strcmp(DATA.ViewPanelsType{panelnum},'viewport')
    updatetextposition(panelnum);
  end

  createfunctions('addcolorbar',panelnum)
  
end

if DATA.issegment3dp
  %Remove slider from all others
  allpanels = 1:length(DATA.Handles.sliders);
  for panelnum = setdiff(allpanels,panelstodo)
    set(DATA.Handles.sliders(panelnum),'Visible','off');
  end
end

if specialview && DATA.DataLoaded
    %functionality to restore same zoom state as before panel settings were
  %changed
  panelstodo = find(DATA.ViewPanels>0);
  incr = zeros(length(panelstodo),1);

  for actpanel = panelstodo
    no = DATA.ViewPanels(actpanel);
    linkednos = SET(no).Linked;
    if length(linkednos) > 1
      for lno = linkednos(2:end)
        SET(lno).NormalZoomState = [];
      end
    end
    type = DATA.ViewPanelsType{actpanel};
    if any(contains({'one'},type))
      zoomstatecurrent = SET(no).NormalZoomState;
      if ~isempty(zoomstatecurrent)
        % compare current state to the original
        % where original would be zoom state without zooming
        zoomstateoriginal = getnewzoomstate(actpanel,no,false);
        szxorig = zoomstateoriginal(2)-zoomstateoriginal(1);
        szxcurrent = zoomstatecurrent(2)-zoomstatecurrent(1);
        if ~isnan(szxorig) && ~isnan(szxcurrent) && ~isapproxequal(szxorig,szxcurrent)
          % ratios not equal -> ther was zooming
          f = szxcurrent/szxorig;
          % recalculate zoom increment
          incr(actpanel) = round((1-f)*7);    
        end
      end
    end
  end
  for panelnum = numel(panelstodo)
    if incr(panelnum)~= 0
        % restore zoom
        SET(nostodo(panelnum)).NormalZoomState = [];
        zoom(incr(panelnum),panelnum);
    end
  end
end
% end of restoring zoom

hideallexceptsax = shouldhideallexceptsax;
if hideallexceptsax
  updatevisibility(hideallexceptsax)
end

drawfunctions('drawselectedframe',DATA.CurrentPanel)
drawfunctions('drawselectedslice',DATA.CurrentPanel)
buttondownfunctions('updatebuttondowns')

DATA.setviewbuttons(0);
segment('updatevolume');
segment('updateflow');
segment('updatemeasurement');
DATA.updatetimebaraxes
%DATA.showplay; %play buttons
%create thumbnails and draw frame around them
drawfunctions('drawthumbnails',isempty(DATA.DATASETPREVIEW));
%profile report

%------------------
function configaxes
%------------------
%configure imageaxes and boxaxes according to current number of panels.

global DATA

%Get global graphics handles
h = DATA.Handles;
guisettings = DATA.GUISettings;

rows = DATA.ViewMatrix(1);
cols = DATA.ViewMatrix(2);

left = guisettings.LeftGapWidth; %0.12;
right = 1 - guisettings.RightGapWidth - 0.02;  %Based on that the report panel can never be more than 220.
bottom = guisettings.BottomGapHeight; %0.013;
top = 1 - guisettings.TopGapHeight;
width = right-left;
height = top-bottom;

%bottom = bottom+0.1;
%height = height-0.1;

%Set the box axes position
set(h.boxaxes,'Position',[left bottom width height],'Visible','on');

%Set up the frame around all the panels
x = [];
y = [];
for rloop=1:(rows-1)
  x = [x nan 0 1]; %#ok<AGROW> 
  y = [y nan 1/rows*rloop 1/rows*rloop]; %#ok<AGROW> 
end

for cloop=1:(cols-1)
  x = [x nan 1/cols*cloop 1/cols*cloop]; %#ok<AGROW> 
  y = [y nan 0 1]; %#ok<AGROW> 
end

x = [x,[nan 0 1 nan 0 1 nan 0 0 nan 1 1]];
y = [y,[nan 0 0 nan 1 1 nan 0 1 nan 0 1]];

set(h.box,'XData',x,'YData',y);

%Place all panels
counter = 1;
for rloop=(rows-1):-1:0
  for cloop=0:(cols-1)
    set(h.imageaxes(counter),...
      'Position',...
      [left+cloop*width/cols bottom+rloop*height/rows width/cols height/rows]);
    counter = counter + 1;
  end
end

if ~DATA.issegment3dp
  %nan all handles and axes not shown
  drawfunctions('setxynan')
else
  %Special things for 3DPrint

  %Hide non visible panels
  for loop = (length(DATA.ViewPanelsType)+1):4
    hc = get(DATA.Handles.imageaxes(loop),'Children');
    set(hc,'Visible','off');
  end

  set(DATA.Handles.box,'visible','off'); %hide the black boxes;
  set(DATA.fig,'Color',[0 0 0]);
  
end

%------------------------------
function im = getgla(no)
%------------------------------
%get the global longitudinal axis for ortho view
global SET

[x,y] = setglacenter(no);

t0 = SET(no).CurrentTimeFrame;
z0 = SET(no).CurrentSlice;

t = 1:SET(no).TSize;
z = 1:SET(no).ZSize;

%Define image axes in xy plane and along z and t axis
xyline = [x;y;t0*ones(size(x));z0*ones(size(x))];
zline = [SET(no).HLA.slice*ones(size(z));SET(no).VLA.slice*ones(size(z));t0*ones(size(z));z];
tline = [SET(no).HLA.slice*ones(size(t));SET(no).VLA.slice*ones(size(t));t;z0*ones(size(t))];

X = meshgrid(xyline(1,:),zline(1,:),tline(1,:));
Y = meshgrid(xyline(2,:),zline(2,:),tline(2,:));
[~,~,T] = meshgrid(xyline(3,:),zline(3,:),tline(3,:));
[~,Z] = meshgrid(xyline(4,:),zline(4,:),tline(4,:));

if SET(no).TSize == 1
  im = interpn(squeeze(SET(no).IM),X,Y,Z,'nearest');
else
  im = interpn(SET(no).IM,X,Y,T,Z,'nearest');
end

%---------------------
function orthoview(no)
%---------------------
% Create orthogonal view for stack no and store necessary fields into SET
% struct
global SET NO
if nargin < 1
  no = NO;
end

SET(no).HLA.slice = round(SET(no).XSize/2);
SET(no).HLA.maxslice = SET(no).XSize;
SET(no).HLA.ZoomState = [];

SET(no).VLA.slice = round(SET(no).YSize/2);
SET(no).VLA.maxslice = SET(no).YSize;
SET(no).VLA.ZoomState = [];

ang = 2*pi/8;
SET(no).GLA.angle = ang;
SET(no).GLA.slice = 0;
SET(no).GLA.ZoomState = [];

setglacenter(no);

%------------------------------
function [x,y] = setglacenter(no)
%------------------------------
%set the center for the global longitudinal axis for the ortho view
global SET

%Define coordinates of midpoint and find image edge
x0 = SET(no).HLA.slice;
y0 = SET(no).VLA.slice;
tlim = zeros(1,4);

%Here should SET(no).ResolutionX and ResolutionY be introduced somehow
tlim(1) = SET(no).ResolutionX*(1-SET(no).HLA.slice)/sin(SET(no).GLA.angle);
tlim(2) = SET(no).ResolutionY*(1-SET(no).VLA.slice)/cos(SET(no).GLA.angle);
tlim(3) = SET(no).ResolutionX*(SET(no).XSize-SET(no).HLA.slice)/sin(SET(no).GLA.angle);
tlim(4) = SET(no).ResolutionY*(SET(no).YSize-SET(no).VLA.slice)/cos(SET(no).GLA.angle);
tlim = sort(tlim);
tmin = tlim(2);
tmax = tlim(3);

%Define extension
res = SET(no).ResolutionY*cos(SET(no).GLA.angle)+SET(no).ResolutionX*abs(sin(SET(no).GLA.angle));

if cos(SET(no).GLA.angle) >= 0
  x = x0 + sin(SET(no).GLA.angle)*(tmin:res:tmax)/SET(no).ResolutionX; %linspace(tmin,tmax,sz);
  y = y0 + cos(SET(no).GLA.angle)*(tmin:res:tmax)/SET(no).ResolutionY; %linspace(tmin,tmax,sz);
else
  x = x0 + sin(SET(no).GLA.angle)*(tmax:-res:tmin)/SET(no).ResolutionX; %linspace(tmin,tmax,sz);
  y = y0 + cos(SET(no).GLA.angle)*(tmax:-res:tmin)/SET(no).ResolutionY;
end

SET(no).GLA.x0 = x(1);
SET(no).GLA.y0 = y(1);

%------------------------------------------------------
function [measure,slice] = getmeasurecoords(panel)
%------------------------------------------------------
%Get coordinates of measurements with respect to the current view
global DATA SET
no = DATA.ViewPanels(panel);
type = DATA.ViewPanelsType{panel};

switch type
  case 'hla'
    measure = struct( ...
      'X',{SET(no).Measure.Z}, ...
      'Y',{SET(no).Measure.Y}, ...
      'Z',{SET(no).Measure.X}, ...
      'T',{SET(no).Measure.T});
    slice = SET(no).HLA.slice;
  case 'vla'
    measure = struct( ...
      'X',{SET(no).Measure.Z}, ...
      'Y',{SET(no).Measure.X}, ...
      'Z',{SET(no).Measure.Y}, ...
      'T',{SET(no).Measure.T});
    slice = SET(no).VLA.slice;
  case 'gla'
    ang = SET(no).GLA.angle;
    res = SET(no).ResolutionY*cos(ang)+SET(no).ResolutionX*abs(sin(ang));
    measure = struct( ...
      'X',{SET(no).Measure.Z}, ...
      'Y',cellfun(@(y,x)1/res*(SET(no).ResolutionY*(y-SET(no).GLA.y0)*cos(ang) + ...
      SET(no).ResolutionX*(x-SET(no).GLA.x0)*sin(ang)), ...
      {SET(no).Measure.Y},{SET(no).Measure.X},'UniformOutput',false), ...
      'Z',cellfun(@(y,x)1/res*(SET(no).ResolutionY*(y-SET(no).GLA.y0)*-sin(ang) + ...
      SET(no).ResolutionX*(x-SET(no).GLA.x0)*cos(ang)), ...
      {SET(no).Measure.Y},{SET(no).Measure.X},'UniformOutput',false), ...
      'T',{SET(no).Measure.T});
    slice = 0;
  otherwise
    measure = struct( ...
      'X',{SET(no).Measure.X}, ...
      'Y',{SET(no).Measure.Y}, ...
      'Z',{SET(no).Measure.Z}, ...
      'T',{SET(no).Measure.T});
    slice = SET(no).CurrentSlice;
end

%-----------------------------------------
function switchtimeframe(incr,synchronize) 
%-----------------------------------------
%swich time frame for all applicable image stacks

global DATA SET

panel = DATA.CurrentPanel;
no = DATA.ViewPanels(panel);

%connect any ongoing interpolation
tools('connectallinterpolation',no);

if nargin == 1
  synchronize = 0;
end
if SET(no).TSize == 1
  return
end
if abs(incr) > SET(no).TSize
  signincr = sign(incr);
  incr = signincr*mod(abs(incr),SET(no).TSize);
end
SET(no).CurrentTimeFrame = SET(no).CurrentTimeFrame+incr;
if SET(no).CurrentTimeFrame>SET(no).TSize
  tfdiff = SET(no).CurrentTimeFrame-SET(no).TSize;
  SET(no).CurrentTimeFrame = tfdiff;
end

if SET(no).CurrentTimeFrame<1
  SET(no).CurrentTimeFrame = SET(no).TSize+SET(no).CurrentTimeFrame;
end

%Calculate percentage
percent = (SET(no).CurrentTimeFrame-1)/(SET(no).TSize-1);

if DATA.issegment3dp || (findindented(DATA.Handles.hideiconholder,'synchronize')||synchronize)
  panels = find(DATA.ViewPanels>0);
  nos = unique(DATA.ViewPanels(panels));
else
  nos = no; 
  linkednos = SET(nos).Linked;
   if length(linkednos) > 1
     % look if the linked no is shown
     panels = find(ismember(DATA.ViewPanels, linkednos));
   else
      panels = panel;     
   end
end
% change timeframe also in linked nos, even if they are not shown
linkednos = cell2mat({SET(nos).Linked});
nos = unique([nos,linkednos]);

%Loop over all image stacks to update time.
for loop=nos
  SET(loop).CurrentTimeFrame = min(max(1,1+round(percent*(SET(loop).TSize-1))),SET(loop).TSize);
end

for p = panels  
  updatedrawlist(p)
  drawfunctions('drawpanel',p)
end

%also update the timebars!
updatetimebars;

%---------------------
function updatetimebars(synchronize) 
%-----------------------
%updates all timebars
global DATA SET NO

%Get global graphics handles
handles = DATA.Handles;

no = NO;

if nargin == 0
  if ~DATA.issegment3dp && findindented(DATA.Handles.hideiconholder,'synchronize')
    synchronize = 1;
  else
    synchronize = 0;
  end
end

%Set timebar for LV panel
if ~isempty(DATA.LVNO)|| strcmp(DATA.ProgramName,'Segment')
  if strcmp(DATA.ProgramName,'Segment')
    t = SET(no).TimeVector*1000;
  else
    t = SET(DATA.LVNO(1)).TimeVector*1000;
  end
  
  if  any(no == DATA.LVNO) || strcmp(DATA.ProgramName,'Segment') || (synchronize && any(ismember(DATA.ViewPanels,DATA.LVNO)))
    
    if strcmp(DATA.ProgramName,'Segment')
      xdata = t(SET(no).CurrentTimeFrame);
    else
      xdata = t(SET(DATA.LVNO(1)).CurrentTimeFrame);
    end
    set(handles.timebarlv,'XData', xdata*[1 1]);
  end
end

%Set timebar for Flow panel
if ~isempty(DATA.FlowNO)
  if strcmp(DATA.ProgramName,'Segment')
    t = SET(no).TimeVector*1000;
  else
    t = SET(DATA.FlowNO(1)).TimeVector*1000;
  end
  
  if  no == DATA.FlowNO || strcmp(DATA.ProgramName,'Segment') || (synchronize && any(DATA.ViewPanels == DATA.FlowNO))
    if strcmp(DATA.ProgramName,'Segment')
      xdata = t(SET(no).CurrentTimeFrame);
    else
      xdata = t(SET(DATA.FlowNO).CurrentTimeFrame);
    end
    set(handles.timebarflow,'XData', xdata*[1 1]);
  end
end

%Set timebar for viewing panel
t = SET(no).TimeVector*1000;
set(handles.timebar,'XData',t(SET(no).CurrentTimeFrame)*[1 1]);

%---------------------
function selectallslices 
%-----------------------
%Selects all slices in panel
global DATA SET NO
slices = slicesinpanel(DATA.CurrentPanel);

SET(NO).StartSlice = slices(1);
SET(NO).EndSlice= slices(end);
SET(NO).CurrentSlice= slices(end);

drawfunctions('drawselectedslice',DATA.CurrentPanel);
drawfunctions('drawtext',DATA.CurrentPanel);

%---------------------
function unselectallslices 
%-----------------------
%Unselects all slices in panel
global DATA SET NO
SET(NO).StartSlice = [];
SET(NO).EndSlice= [];
drawfunctions('drawselectedslice',DATA.CurrentPanel);
drawfunctions('drawtext',DATA.CurrentPanel);

%---------------------
function viewallimagestacks(viewpaneltype) 
%-----------------------
%roll outs all image stacks (up to 6) into panels in one view mode.
global SET DATA

N = length(SET);
switch N
  case 1
    rows = 1;
    cols = 1;
  case 2
    rows = 1;
    cols = 2;
  case 3
    rows = 1;
    cols = 3;
  case 4
    rows = 2;
    cols = 2;
  otherwise%case {5,6}
    rows = 2;
    cols = 3;
end

nos = 1:min(N,6);
if any(ismember(DATA.ViewPanelsType,{'orth','hla','gla','vla'}))
  % reset orthogonal view to one panel view to avoid crashing 
  setview(1,1);
end
if nargin==0
  setview(rows,cols,nos)
else
  vptcell = cell(size(nos));
  [vptcell{:}] = deal(viewpaneltype);
  setview(rows,cols,nos,vptcell)
end
%---------------------
function fastswitchslice(sliceincr,synchronize) 
%-----------------------
%This is a first attempt at a fast slice toggling using a first render with
%only one timeframe. The design  is that a timer is initiated

global DATA SET
persistent storetimer
if nargin == 1
  synchronize = 0;
end

panel = DATA.CurrentPanel;
scale = getscale(panel);
no = DATA.ViewPanels(panel);
SET(no).CurrentSlice = SET(no).CurrentSlice+sliceincr;

if SET(no).CurrentSlice <= 0
  SET(no).CurrentSlice = 1;
end

if SET(no).CurrentSlice >SET(no).ZSize
  SET(no).CurrentSlice = SET(no).ZSize;
end

if ~isempty(storetimer)
  clear('storetimer')
end

storetimer=timer('ExecutionMode','singleShot',...
  'TimerFcn',sprintf('viewfunctions(''switchslice'',%d,%d)',0,synchronize), ...
  'Name','FastSwitchSlice',...
  'StartDelay',1);
outim = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
im = calcfunctions('remapuint8viewim',outim,no);

DATA.Handles.imagehandles(panel).CData = imresize(im,scale);

for i = 2:length(DATA.drawlist{panel})
  eval(DATA.drawlist{panel}{i})
end

start(storetimer)

%------------------------------------------
function switchslice(sliceincr,synchronize) 
%------------------------------------------
%swich slice for all applicable image stacks

global DATA SET

panel = DATA.CurrentPanel;
no = DATA.ViewPanels(panel);
slice = SET(no).CurrentSlice+sliceincr;

%connect any ongoing interpolation
tools('connectallinterpolation',no);

if nargin == 1
  synchronize = 0;
end

%If one then we need to re-generate viewpanels
if strcmp(DATA.ViewPanelsType{panel},'one') || strcmp(DATA.ViewPanelsType{panel},'orth')
  DATA.ViewIM{panel} = [];
end
        
if slice > SET(no).ZSize
  slice = SET(no).ZSize;
end

if slice < 1
  slice = 1;
end

SET(no).CurrentSlice = slice;
SET(no).StartSlice = slice;

%segment('updatevolume')
SET(no).EndSlice = slice;

if ~DATA.issegment3dp && (findindented(DATA.Handles.hideiconholder,'synchronize') || synchronize)
%   pno = segment('updateparallelsets',no);
%   panels = find(ismember(DATA.ViewPanels,[pno no]));
  panels = find(DATA.ViewPanels>0);
  nos = unique(DATA.ViewPanels(panels));
else
  nos = no;
  linkednos = SET(nos).Linked;
  if linkednos > 1
    % look if the linked no is shown
    panels = find(ismember(DATA.ViewPanels, linkednos));
  else
    panels = panel;
  end
end
% change timeframe also in linked nos, even if they are not shown
linkednos = cell2mat({SET(nos).Linked});
nos = unique([nos,linkednos]);
%update Currentslice in all linked nos
invisiblenos = setdiff(nos,DATA.ViewPanels);
if not(isempty(invisiblenos))
  for n = invisiblenos
    thisslice = findcorrespondingslice(slice,no,n);
    if ~isempty(thisslice) && ~isnan(thisslice) && thisslice > 0 && thisslice <= SET(n).ZSize
      SET(n).CurrentSlice = thisslice;
      SET(n).StartSlice = thisslice;
      SET(n).EndSlice = thisslice;
    end
  end
end

%Loop over all panels to update slice
for loop = panels
  thisno = DATA.ViewPanels(loop);
  if ~isequal(no,thisno)
    thisslice = findcorrespondingslice(slice,no,thisno);
    if (~isempty(thisslice)) && (~isequal(thisslice,SET(thisno).CurrentSlice))
      SET(thisno).CurrentSlice = thisslice;
      SET(thisno).StartSlice = thisslice;
      SET(thisno).EndSlice = thisslice;
    end
  end
  %Check if need to recompute the image
  if strcmp(DATA.ViewPanelsType{loop},'one') || strcmp(DATA.ViewPanelsType{loop},'orth')
    DATA.ViewIM{loop} = [];
  end
end

%this draws and updates the panel image intersections.
drawfunctions('drawplaneintersections')

%Do graphical update. If in montage view this corresponds to changing the
%highlighted frame. If other then we need to update viewim and redraw the
%panel.

for p = panels
  updatedrawlist(p)  
end

for p = panels
  if any(strcmp(DATA.ViewPanelsType{p},{'montage','montagerow'}))
    drawfunctions('drawselectedslice',p)
    drawfunctions('drawtext',p)
  elseif strcmp(DATA.ViewPanelsType{p},'montagesegmented')
    drawfunctions('drawselectedslice',p)
    drawfunctions('drawtext',p)
    drawfunctions('drawpanel',p)
  else
    %making the viewim empty triggers generation of a new viewim
    drawfunctions('drawpanel',p)
  end
end

%----------------------------------------------------------------
function thisslice = findcorrespondingslice(refslice,refno,thisno)
%----------------------------------------------------------------
%Find correspondingslice given a refslice in a refno for the stack thisno

global SET

%Compute zaxis for ref
v1ref = SET(refno).ImageOrientation(1:3);
v2ref = SET(refno).ImageOrientation(4:6);
zref = cross(v1ref,v2ref);

%Compute zaxis for this
v1this = SET(thisno).ImageOrientation(1:3);
v2this = SET(thisno).ImageOrientation(4:6);
zthis = cross(v1this,v2this);

%Compare them
angle = acos(abs(sum(zref.*zthis)))/pi*180; %angle between them in degrees

if angle<20
  pos = calcfunctions('xyz2rlapfh',refno,SET(refno).XSize/2,SET(refno).YSize/2,refslice);
  xyz = calcfunctions('rlapfh2xyz',thisno,pos(1),pos(2),pos(3));
  
  thisslice = round(xyz(3));
  
  if (thisslice > SET(thisno).ZSize) || (thisslice < 1) || isnan(thisslice)
    thisslice = [];
  end
  
else
  thisslice = [];
  return
end

%----------------------------------
function zoomstate = updateglazoomstate(no,ysz)
%----------------------------------
%This function helps with updating the zoom for the gla view when in
%orthoview mode.
global SET
xzoomsz = max(SET(no).XSize,SET(no).YSize);
repos = (ysz-xzoomsz)/2;
zoomstate = 0.5+[repos+[0;xzoomsz];0;SET(no).ZSize];

%---------------------------------------
function viewhideallspecial_Callback(varargin) 
%---------------------------------------
%all hide buttons
global DATA
stateandicon = iconson('hideall');
state = stateandicon{1};
hidecell = getallhidecionsnames;
stateandicon = iconson(hidecell);
availableicons = find(cellfun(@(x) isa(x,'myicon'),stateandicon(:,2)))';
if state
  visiblestate = 'off';
  DATA.LastHideIconState = stateandicon;
  for i = availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay = icon.cdataIndent;
    icon.isindented = 1;
  end
else
  visiblestate = 'on';
  if ~isempty(DATA.LastHideIconState)
    stateandicon = DATA.LastHideIconState;
    for i = availableicons
      icon = stateandicon{i,2};
      if stateandicon{i,1}
        icon.cdataDisplay = icon.cdataIndent;
        icon.isindented = 1;
      else
        icon.cdataDisplay = icon.cdata;
        icon.isindented = 0;
      end
    end
    DATA.LastHideIconState = [];
  end
end

viewfunctions('updatevisibility');
DATA.Handles.configiconholder.render;
DATA.Handles.permanenticonholder.render;
viewhideframeandslice(visiblestate);

%---------------------------------------
function viewhideall_Callback(varargin) 
%---------------------------------------
%all hide buttons
global DATA
stateandicon=iconson('hideall');
if nargin > 0
  state = varargin{1};
else
  state=stateandicon{1};
end
hidecell = getallhidecionsnames;

stateandicon = iconson(hidecell);
availableicons = find(cellfun(@(x) isa(x,'myicon'),stateandicon(:,2)))';
if state
  DATA.LastHideIconState = stateandicon;
  for i = availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay=icon.cdataIndent;
    icon.isindented=1;
  end
  visiblestate = 'off';
else
  for i = availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay=icon.cdata;
    icon.isindented=0;
  end
  visiblestate = 'on';
end

viewfunctions('updatevisibility');
DATA.Handles.permanenticonholder.render;
if ~DATA.issegment3dp
  DATA.Handles.configiconholder.render;
end
viewhideframeandslice(visiblestate);

%-------------------------------------------
function viewhideframeandslice(visiblestate)
%-------------------------------------------
global DATA
set(DATA.Handles.selectedslice(:),'Visible', visiblestate);
set(DATA.Handles.selectedframe(DATA.CurrentPanel),'Visible', visiblestate);

%------------------------------------------
function hidecell = getallhidecionsnames
%-------------------------------------------
hidecell = {'hideplus','hidemar','hidescar','hidescarextent','hidescarmanual',... 'hidepins',
  'hideintersections','hideothercontour','hideinterp','hidelv','hiderv',...
  'hideroi','hidemeasure','hidepoint','hidetext','hidegeneralpen',...
  'hidela','hidera'};%,'colorbar'};

%-------------------------------------------
function state = iconson(name)
%------------------------------------------
%if given name of button return state if run with no input returns
%states of all buttons. if given cell with multiple button return icons
%in order of request.
global DATA

%Get global graphics and icons handles
handles = DATA.Handles;
dataicons = DATA.Icons;

state = 0;
if ~DATA.issegment3dp
  icons = [handles.permanenticonholder.iconCell{:},dataicons.lviconcell{:},...
    dataicons.rviconcell{:},dataicons.roiflowiconcell{:},dataicons.viabilityiconcell{:},...
    dataicons.analysisiconcell{:},dataicons.imageiconcell{:},dataicons.hidecell{:}...    
    ];
  if not(contains(DATA.ProgramName,'CT'))
    icons = [icons, dataicons.functioniconcell{:},dataicons.laraiconcell{:}, handles.playiconholder.iconCell{:}];%,icons,dataicons.generalpeniconcell{:}];
  else
    icons = [icons,dataicons.laraiconcell{:}];
  end
else
  return; %In segment3DP it is done elsewhere
end

N = length(icons);

if nargin==1
  %return icon state
  if ~iscell(name)
    ind=find(strcmp(name,{icons.name}));
    if length(ind)>1
      ind=ind(1);
    end
    
    if isempty(ind)
      state={0,nan};
    else
      state={icons(ind).isindented,icons(ind)};
    end
    
  else
    n=length(name);
    state=cell(n,2);
    %Get all icons and states for names in cell
    counter=1;
    indlist=zeros(1,n);
    for i=1:n
      ind=find(strcmp(name(i),{icons.name}));
      
      if isempty(ind)
        ind=nan;
      end
      
      if length(ind)>1
        ind=ind(1);
      end
      
      indlist(i)=ind;
    end
    
    for i=1:n
      if isnan(indlist(i))
        state{i,1} = 0;
        state{i,2} = nan;
      else
        state{i,1} = icons(indlist(i)).isindented;
        state{i,2} = icons(indlist(i));
      end
    end
  end
else
  state=cell(N,2);
  for i=1:N
    state{i,1} = icons{i}.name;
    state{i,2} = icons{i}.isindented;
  end
end


%-------------------------
function allhidden
%----------------------
%show/hide all overlayes in the image view

global DATA
stateandicon=iconson('hideall');
hideallstate=stateandicon{1};
hideallicon=stateandicon{2};

if hideallstate
  hidecell={'hideplus','hidemar','hidescar','hidescarextent','hidescarmanual',...'hidepins',
    'hideintersections','hideothercontour','hideinterp','hidelv','hiderv','hideroi','hidemeasure','hidepoint','hidetext'};%,'colorbar'};
  stateandicon = iconson(hidecell);
  availableicons = cellfun(@(x) isa(x,'myicon'),stateandicon(:,2));
  state=[stateandicon{availableicons,1}];
  
  if any(state==0)
    hideallicon.undent
    if contains(DATA.ProgramName,'3D')
      DATA.Handles.configiconholder.render;
    else
      DATA.Handles.permanenticonholder.render;
      DATA.Handles.playiconholder.render;
      DATA.Handles.approveiconholder.render;
    end
  end
end

%--------------------------
function hide(names,states) 
%--------------------------
%Gets and sets icons with specific name in the different iconholders even those not showing at the moment according to state
global DATA
stateandicon = viewfunctions('iconson',names);

for i =1:size(stateandicon,1)
  if isprop(stateandicon{i,2},'isindented')
    stateandicon{i,2}.isindented = states(i);
  end
end

DATA.Handles.configiconholder.render;
% viewfunctions('setview')
viewfunctions('updatevisibility');

%---------------------------------------
function hide_Callback(varargin) 
%---------------------------------------
%Toggle visibility for the imageoverlay
viewfunctions('updatevisibility');
allhidden

%---------------------------------
function viewhidecolorbar_Callback 
%---------------------------------
%Toggle visibility of colorbar
global DATA

stateandicon=iconson('colorbar');
state=stateandicon{1};
%stateandicon{2}.isindented=state;
if ~state
  DATA.GUISettings.ShowColorbar = false;
else
  DATA.GUISettings.ShowColorbar = true;
end
%do full refresh with colorbar
for p = find(DATA.ViewPanels)
  createfunctions('addcolorbar',p)
end

allhidden

%-----------------------------------
function viewhidepap_Callback 
%---------------------------------
%Toggle visibility of papillary overlay. This overlay is
%incorporated in the viewim so resetting them and assuring it doesnt
%trigger when button is indented

global DATA
viewfunctions('updatevisibility')
DATA.ViewIM = cell(size(DATA.ViewIM));
for p = find(DATA.ViewPanels)
  drawfunctions('drawpanel',p);
end
allhidden
%--------------------------
function viewhidemanualinteraction_Callback 
%--------------------------
%toogle visibility of manual interaction of Scar/MaR. This overlay is
%incorporated in the viewim so resetting them and assuring it doesnt
%trigger when button is indented.
global DATA

viewfunctions('updatevisibility')
DATA.ViewIM = cell(size(DATA.ViewIM));
for p = find(DATA.ViewPanels)
  drawfunctions('drawpanel',p);
end
allhidden

%--------------------------------
function viewinterp_Callback(bool) 
%--------------------------------
%Update if the view interpolated view is selected / unselected

global DATA
if nargin < 1
  bool = ~DATA.Pref.ViewInterpolated;
end
DATA.Pref.ViewInterpolated = bool;
DATA.ViewIM = cell(size(DATA.ViewIM)); %This clear ViewIM and it is reconstructed
for p = find(DATA.ViewPanels)
  drawfunctions('drawpanel',p);
end
drawfunctions('drawplaneintersections');

%------------------------------
function viewportviewfrom(type) 
%------------------------------
%Update view direction in viewport depending on direction

global DATA

%Get object handle
viewport = DATA.LevelSet.ViewPort;

if isa(viewport,'myvolshow')
  vol = true;
else
  vol = false;
end

%Get camera position
target = viewport.getcameratarget;
pos = viewport.getcameraposition;

%Check distance to center (origo) for camera
r = sqrt(sum((pos-target).^2));

if vol %For volume objects, the position and upvector depends on image orientation
  if strcmp(DATA.LevelSet.imageorientation,'Transversal')
    switch type
      case 'transversal'
        posxyz = [r 0 0];
        upvxyz = [0 0 -1];
      case 'sagittal'
        posxyz = [0 -r 0];
        upvxyz = [0 0 -1];
      case 'coronal'
        posxyz = [0 0 -r];
        upvxyz = [0 -1 0];
      case 'back'
        posxyz = [0 r 0];
        upvxyz = [0 0 -1];
      case 'left'
        posxyz = [-r 0 0];
        upvxyz = [0 0 -1];
      case 'bottom'
        posxyz = [0 0 r];
        upvxyz = [0 -1 0];
    end
  elseif strcmp(DATA.LevelSet.imageorientation,'Sagittal')
    switch type
      case 'transversal'
        posxyz = [0 0 -r];
        upvxyz = [0 -1 0];
      case 'sagittal'
        posxyz = [r 0 0];
        upvxyz = [0 -1 0];
      case 'coronal'
        posxyz = [0 -r 0];
        upvxyz = [1 0 0];
      case 'back'
        posxyz = [-r 0 0];
        upvxyz = [0 -1 0];
      case 'left'
        posxyz = [0 0 r];
        upvxyz = [0 -1 0];
      case 'bottom'
        posxyz = [0 r 0];
        upvxyz = [1 0 0];
    end

  elseif strcmp(DATA.LevelSet.imageorientation,'Coronal')
    switch type
      case 'transversal'
        posxyz = [r 0 0];
        upvxyz = [0 -1 0];
      case 'sagittal'
        posxyz = [0 0 r];
        upvxyz = [0 -1 0];
      case 'coronal'
        posxyz = [0 -r 0];
        upvxyz = [0 0 1];
      case 'back'
        posxyz = [0 0 -r];
        upvxyz = [0 -1 0];
      case 'left'
        posxyz = [-r 0 0];
        upvxyz = [0 -1 0];
      case 'bottom'
        posxyz = [0 r 0];
        upvxyz = [0 0 1];
    end
  end
else %Image orientation does not matter for surface objects
  %imageorientation = SET(NO).ImageOrientation;
  imageorientation = [1 0 0 0 1 0];

  zdir = cross(...
    imageorientation(1:3),...
    imageorientation(4:6));

  Rxyz2rlapfh = [imageorientation(4:6)' imageorientation(1:3)' -zdir'];

  Rrlapfh2xyz = inv(Rxyz2rlapfh);

  %Define in rlapfh coordinates

  switch type
    case 'transversal'
      posrlapfh = -[0 r 0];
      upvrlapfh = -[0 0 1];
    case 'sagittal'
      posrlapfh = [-r 0 0];
      upvrlapfh = -[0 0 1];
    case 'coronal'
      posrlapfh = -[0 0 r];
      upvrlapfh = [-1 0 0];
    case 'back'
      posrlapfh = [r 0 0];
      upvrlapfh = -[0 0 1];
    case 'left'
      posrlapfh = [0 r 0];
      upvrlapfh = -[0 0 1];
    case 'bottom'
      posrlapfh = [0 0 r];
      upvrlapfh = [-1 0 0];
  end
end

if ~vol 
  %Convert to xyz coordinates
  posxyz = Rrlapfh2xyz*(posrlapfh')+target(:); %#ok<MINV>
  upvxyz = Rrlapfh2xyz*(upvrlapfh'); %#ok<MINV>
end

%Update camera position
viewport.setcamera(posxyz,upvxyz);
%viewport.synchronize_volshow_with_axes();
%viewport.synchronize_axes_with_volshow();
viewport.showgraphics();

%-----------------------------------
function updatelinewidth 
%-----------------------------------
%update the line width
global DATA

%Get global graphics handles
h = DATA.Handles;

linewidth = DATA.Pref.LineWidth;
if ~isempty(linewidth)
  colortypes = 'cgbrwkym';
  type = {'Endo','Epi','RVEndo','RVEpi','LA','RA'};
  for panel = 1: length(DATA.ViewPanels)
    h.roicurrent(panel).LineWidth = linewidth + 1; %ROI contour
    h.scarcontour(panel).LineWidth = linewidth; %scar yellow line
    h.weightedscarcontour(panel).LineWidth = linewidth; %scar pink line
    h.mocontour(panel).LineWidth = linewidth; %MO dotted line
    h.moextentcontour(panel).LineWidth = linewidth; %MO plain line
    h.marcontour(panel).LineWidth = linewidth; %MaR contour
    for loop = 1: length(type)
      h.([lower(type{loop}),'contour'])(panel).LineWidth = linewidth;
    end
    for cind = 1:length(colortypes)
      h.([colortypes(cind),'roi'])(panel).LineWidth = linewidth;
    end
  end
end

%------------------------------------
function updateintersectionmarkersize
%------------------------------------
%update marker size for intersection points

global DATA

%update marker size for all intersection
for panel = 1:size(DATA.Handles.endocontourintersection,2)  
  set(DATA.Handles.endocontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.endocontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.epicontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.epicontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvendocontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvendocontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvepicontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
  set(DATA.Handles.rvepicontourintersection(panel),'markersize',DATA.Pref.MarkerSize);
end
  
if DATA.Pref.BlackWhite  
  set(DATA.Handles.endocontourintersection(panel),'color',[1 1 1]);
  set(DATA.Handles.endocontourintersection(panel),'color',[1 1 1]);
  set(DATA.Handles.epicontourintersection(panel),'color',[1 1 1]);
  set(DATA.Handles.epicontourintersection(panel),'color',[1 1 1]);
  set(DATA.Handles.rvendocontourintersection(panel),'color',[1 1 1]);
  set(DATA.Handles.rvendocontourintersection(panel),'color',[1 1 1]);
  set(DATA.Handles.rvepicontourintersection(panel),'color',[1 1 1]);
  set(DATA.Handles.rvepicontourintersection(panel),'color',[1 1 1]);
end

%update display
for panel = 1: length(DATA.ViewPanels)
  drawfunctions('drawintersections',panel)
end

%------------------------------------------------------
function z = reshape2layout(im,panel,outsideelement)
%------------------------------------------------------
%Convert a 3D array to an layout:ed image with cols, and rows
global DATA SET
no = DATA.ViewPanels(panel);
imdim = zoomfunctions.getxysize(no,panel);
if nargin < 3
  outsideelement = 0;
end
sz1 = DATA.ViewPanelsMatrix{panel}(1);
sz2 = DATA.ViewPanelsMatrix{panel}(2);
z = repmat(im(1).*outsideelement,sz2*imdim.XSize,sz1*imdim.YSize);
loop = 1;
for slice = 1:SET(no).ZSize
  c = 1+mod(loop-1,sz1);
  r = ceil(loop/sz1);
  z(...
    (1+(r-1)*imdim.XSize):(r*imdim.XSize),...
    (1+(c-1)*imdim.YSize):(c*imdim.YSize)) = im(imdim.XStart:imdim.XEnd,imdim.YStart:imdim.YEnd,slice);
  loop = loop+1;
end

%------------------------------------------------------
function updateflowresultstack  
%------------------------------------------------------
%check if flow stack in result panel should be updated

global SET NO DATA

[~,~,flowno] = findfunctions('findno');
% need to switch stack only if there is flow at all && current stack is a
% flowstack && current stack is not already presented in result panel
if ~isempty(flowno) && any(flowno == NO) && ~isempty(DATA.FlowNO) && not(any(DATA.FlowNO == SET(NO).Linked))
  % get parent stack
  if ~isempty(SET(NO).Flow) && isfield(SET(NO).Flow,'MagnitudeNo')
    magno = SET(NO).Flow.MagnitudeNo;
    if SET(magno).RoiN > 0
      DATA.FlowNO = magno;
      DATA.FlowROI = SET(magno).RoiCurrent;
      DATA.updateflow;
    end
   end
end

%---------------------------------------
function updatetoolhidestate(currenttool)  
%---------------------------------------
%function to update hide icons depending on the choosen tool, callback come
%from updatebuttondowns and pen_buttondown
%This function helps to avoid that user draw with pen but do not see the
%results
global DATA
if ~DATA.LastToolHideState
  switch currenttool
    case {'Endo','Epi'}
      iconname = 'hidelv';
    case {'RVEndo','RVEpi'}
      iconname = 'hiderv';
    case {'Roi','RoiPut'}
      iconname = 'hideroi';
    case {'Scar','MO'}
      iconname = 'hidescar';
    case 'MaR'
      iconname = 'hidemar';
    otherwise
      return
  end
  stateandicon = iconson(iconname);
  if stateandicon{1}
    DATA.Handles.configiconholder.undent(iconname,1)
  end
end

%-----------------------------
function relevantmode_Callback
%-----------------------------
global DATA SET NO
switch DATA.ProgramName
  case {'Segment 3DPrint','Segment CT'}
    return
end

updatestackmodebuttonslook;
bgroup = DATA.Handles.stacksuibuttongroup;
selectedbutton = bgroup.SelectedObject;
if contains(selectedbutton.Tag,'relevant')
  isrelevantmode = true;
else
  isrelevantmode = false;
end
DATA.ShowRelevantStacksOnly = isrelevantmode;

updatepreview = true;
drawfunctions('drawthumbnails',updatepreview);
segment('updateslider');
if DATA.ShowRelevantStacksOnly
  no = NO;
  viewstr = '';
  switch DATA.CurrentTheme
    case 'function'
      viewstr = 'lv';
    case 'roi'
      viewstr = 'flow';
    case 'strain'
      imgplane = SET(no).ImageViewPlane;
      if contains(imgplane,'short','IgnoreCase',true)
        viewstr = 'strainsax';
      elseif ismember(imgplane,{'2CH','3CH','4CH'})
        viewstr = 'strainlax';    
      end
    case 'scar'
      viewstr = 'cinescar';
    case 'txmap'
      % nothing for now
    otherwise
      set(bgroup.Children,'Enable','off');
      return
  end
  set(bgroup.Children,'Enable','on');
  isrelevantmodeempty = ~DATA.FoundRelevantStacks;
  updatepseudorelevantstackmode(isrelevantmodeempty) 
  if ~isrelevantmodeempty
    % relevant stack were found, check if all NOs in panels are also in
    if ~all(ismember(nonzeros(DATA.ViewPanels),DATA.RelevantStacks)) %check if visible NOs are part of RelevantStacks
      if isempty(viewstr)
        % update
        setview(1,1,no);
      elseif strcmpi(viewstr,'cinescar')
        setview(1,1,no);
      else
        segment('viewspecial_Callback',viewstr);
      end
    end
  end
else
  enabledthemes = {'function', 'roi', 'scar', 'strain', 'txmap'};
  enablestate = 'off';
  findfunctions('findrelevantstacks'); % this updates DATA.FoundRelevantStacks
  if ismember(DATA.CurrentTheme, enabledthemes) && DATA.FoundRelevantStacks
    enablestate = 'on';
    if isempty(bgroup.UserData)
      bgroup.UserData = 'updated';
    end
  else
    if isempty(bgroup.UserData)
      bgroup.UserData = 'init';
    end
  end
  set(bgroup.Children, 'Enable', enablestate);
 end


%--------------------------
function togglestackmode
%--------------------------
% function to toggle between relevant and all mode for the stacks
global DATA

bg = DATA.Handles.stacksuibuttongroup;
% Get all children of the button group
buttons = bg.Children;

% Find the unselected button
unselectedbutton = buttons(buttons ~= bg.SelectedObject);

% Switch to the unselected button
bg.SelectedObject = unselectedbutton;

% Manually trigger the SelectionChangedFcn
eval(bg.SelectionChangedFcn);

%----------------------------------
function updatestackmodebuttonslook
%----------------------------------
% function to update the look based on GUI main color
global DATA

% Check background theme
bgroup = DATA.Handles.stacksuibuttongroup;
if all(DATA.GUISettings.BackgroundColor == [0.2118 0.2353 0.2824])
  unselectedforegroundcolor =  [220 220 220]/255;
  selectedbackgroundcolor = [215 233 254]/255;
  selectedforegroundcolor = DATA.GUISettings.BoxAxesColor;
  javaunselectedcolor = java.awt.Color(int32(220), int32(220), int32(220));
else
  unselectedforegroundcolor = DATA.GUISettings.ForegroundColor;
  selectedbackgroundcolor = [215 233 254]/255;
  selectedforegroundcolor = DATA.GUISettings.BoxAxesColor;
  javaunselectedcolor = java.awt.Color.black;
end
% Get all button children of the button group
buttons = findobj(bgroup.Children,'style','togglebutton');

% Find the unselected button and update its appearance
unselectedbutton = buttons(buttons ~= bgroup.SelectedObject);
set(unselectedbutton,'BackgroundColor',DATA.GUISettings.BackgroundColor, ...
  'Enable','on',...
  'ForegroundColor',unselectedforegroundcolor,...
  'FontWeight','normal');

if ~isempty(bgroup.UserData) && strcmpi(bgroup.UserData,'init')
  try
    myworkon
    jButton = findjobj_fast(unselectedbutton);
    if ~isempty(jButton)
      jButton.setForeground(javaunselectedcolor); % Restore Java color
      jButton.repaint(); % Refresh GUI
    end
    bgroup.UserData = 'updated';
    myworkoff
  catch
  end
end

% Update appearance of the selected button
set(bgroup.SelectedObject,'BackgroundColor',selectedbackgroundcolor,...
  'ForegroundColor',selectedforegroundcolor,...
  'Enable','on',...
  'FontWeight','bold');

%-------------------------------------
function updatepseudorelevantstackmode(isrelevantmodeempty)
%-------------------------------------
global DATA
if isrelevantmodeempty
  bgroup = DATA.Handles.stacksuibuttongroup;

  % Get all button children of the button group
  buttons = findobj(bgroup.Children,'style','togglebutton');

  unselectedforegroundcolor = [220 220 220]/255;
  selectedbackgroundcolor = [215 233 254]/255;
  selectedforegroundcolor = DATA.GUISettings.BoxAxesColor;

  % Find the unselected button and update its appearance
  set(bgroup.SelectedObject,'BackgroundColor',selectedbackgroundcolor, ...
    'ForegroundColor',selectedforegroundcolor,...
    'Enable','off',...
    'FontWeight','normal');
  unselectedbutton = buttons(buttons ~= bgroup.SelectedObject);
  % Change appearance of the unselected button to selected button
  set(unselectedbutton,'BackgroundColor',DATA.GUISettings.BackgroundColor, ...
    'ForegroundColor',unselectedforegroundcolor,...
    'Enable','off',...
    'FontWeight','bold');
else
  updatestackmodebuttonslook;
end

%-------------------------------------
function maxnumroitext = getnumberofmaxroitext(panel)
%-------------------------------------
% Get number of possible roi text based on panel view type
global DATA
paneltype = DATA.ViewPanelsType{panel};
if contains(paneltype,'montage')
  maxnumroitext = 35; % this number is used in createfunctions
else
  maxnumroitext = 12;
end

%-------------------------------------------------
function hideallexceptsax = shouldhideallexceptsax
%-------------------------------------------------
% Function to check whether non-SAX contour should be hidden
global DATA
if DATA.issegmentcmr || DATA.issegment
  hideallexceptsax = DATA.Handles.permanenticonholder.findindented('hidenonsax'); 
else
  hideallexceptsax = false;  
end

%-------------------------------------------------
function hideshowepi
%-------------------------------------------------
% Function to toggle between show and hide LV epi
global DATA
h = DATA.Handles;
epih = [h.epicontour,h.epiinterp, h.epicontourintersection];
currentstates=(get(h.epicontour,'Visible'));
if any(strcmp(currentstates, 'on'))
    newstate = 'off';
else
    newstate = 'on';
end

set(epih,'Visible', newstate);

%----------------------------
function addtextbackground
%----------------------------
% Function to toggle between show and hide text background
global DATA

stateandicon=iconson('textbackground');
state=stateandicon{1};
if ~state
  DATA.Pref.BackgroundColor = false;
else
  DATA.Pref.BackgroundColor = true;
end

for p = find(DATA.ViewPanels)
  drawfunctions('drawpanel',p);
end
drawfunctions('drawplaneintersections');