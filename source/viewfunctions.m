function [varargout] = viewfunctions(varargin)
% Functions for querying about the view such as slices in panel and
% zoomstate
% Klas
%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%-------------------------------------
function slices = slicesinpanel(panel)
%-------------------------------------
global DATA SET

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

if ~DATA.Pref.ViewInterpolated || any(strcmp(DATA.ViewPanelsType{panel},{'montage','montagesegmented','montagerow'}))
  scale = 1;
else
  scale = 2;
end

%---------------------------------------
function updatezoomandaspectratio(panel)
%---------------------------------------
%Update zoom state

global DATA SET

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
    case {'trans3DP','sag3DP','cor3DP','speedim'}
      switch DATA.ViewPanelsType{panel}
        case 'trans3DP'
          color = 'r';
        case 'sag3DP'
          color = 'g';
        case 'cor3DP'
          color = 'b';
        case 'speedim'
          color = SET(no).LevelSet.Pen.Color;
      end    
      [xsz,ysz,xres,yres] = segment3dp.tools('viewsizeandres',color,no);
      zoomstate = SET(no).LevelSet.View.([upper(color) 'ZoomState']);
  end
end
if isempty(zoomstate) 
  zoomstate = getnewzoomstate(panel,no);
else
  % this is the previous zoom state
  if ~isempty(SET(no).NormalZoomState) && strcmp(DATA.ViewPanelsType{panel},'one')
    zoomstateold = SET(no).NormalZoomState;
    xsizeold = zoomstateold(2)-zoomstateold(1);
    ysizeold = zoomstateold(4)-zoomstateold(3);
    zoomfold = xsizeold/ysizeold;

    %recalculate zoom state without updating SET struct
    zoomstatenew = getnewzoomstate(panel,no,false);
    xsizenew = zoomstatenew(2)-zoomstatenew(1);
    ysizenew = zoomstatenew(4)-zoomstatenew(3);
    zoomfnew = xsizenew/ysizenew;  
    if ~isapproxequal(zoomfnew,zoomfold)
      % if the zoom ratio has changed, update to the newest
      SET(no).NormalZoomState = zoomstatenew;
      zoomstate = zoomstatenew;
    end
  end
end

plotboxaspectratio = [ysz/yres xsz/xres 1];
dataaspectratio = [1/yres 1/xres 1];

if ~DATA.Silent
  DATA.Handles.imageaxes(panel).DataAspectRatio = dataaspectratio;
  DATA.Handles.imageaxes(panel).PlotBoxAspectRatio = plotboxaspectratio;
  DATA.Handles.imageaxes(panel).XLim = scale*zoomstate(1:2);
  DATA.Handles.imageaxes(panel).YLim = scale*zoomstate(3:4);
end

%---------------------
function updatedrawlist(panel,isonlyannotation)
%-----------------------
%this function should be used at the end of a addition of annotation for
%a panel.
%Linked stacks are taken care of within annotationinpanel
do = viewfunctions('annotationinpanel',panel);
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
  addtodrawlist(panel,@()segment3dp.tools('helprender'));
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
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'contour'])
        addtodrawlist(panel,@()drawfunctions('drawcontours',panel,types{i}));
      end
    end
    
    %Interpolation
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'interp'])
        addtodrawlist(panel,@()drawfunctions('drawinterp',panel,[types{i},'Interp']));
      end
    end
    
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
    
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'contour'])
        addtodrawlist(panel,@()drawfunctions('drawcontours',panel,types{i}));
      end
    end
    
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'interp'])
        addtodrawlist(panel,@()drawfunctions('drawinterp',panel,[types{i},'Interp']));
      end
    end
    
    if do.measures
      addtodrawlist(panel,@()drawfunctions('drawmeasures',panel));
    end
    
    if do.text
      addtodrawlist(panel,@()drawfunctions('drawtext',panel));
    end
    
    for c = 'cgymkwrb'
      if do.([c,'roi'])
        addtodrawlist(panel,@()drawfunctions('drawroi',panel,c));
      end
    end
    
    if do.point
      addtodrawlist(panel,@()drawfunctions('drawpoint',panel));
    end
    
    if do.scar
      addtodrawlist(panel,@()drawfunctions('drawviability',panel));
    else
      do.manualinteraction = 0;
      DATA.ViewIM{panel} = [];
    end
    
    if do.manualinteraction
      addtodrawlist(panel,@()drawfunctions('showviabilityedits',panel));
    else
      DATA.ViewIM{panel} = [];
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
    
    %types = {'r','b','g','s'};
    %for i = 1:length(types)
    %  if do.([types{i} 'contour'])
    %    DATA.drawlist{panel}{end+1} = @()drawfunctions(''draw3dpoutline'',%d,''%s'')',panel,types{i});
    %  end      
    %end
    if do.threedpcontour
      addtodrawlist(panel,@()drawfunctions('draw3dpoutline',panel,'threedp'));
    end
    
    if do.point
      addtodrawlist(panel,@()drawfunctions('drawpoint',panel));
    end
    
end
%--------------------------------------------------------
function updatedrawlist_helper(panel,do,isonlyannotation)
%--------------------------------------------------------
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
  addtodrawlist(panel,sprintf('segment3dp.tools(''helprender'');'));
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
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    for i = 1:length(types)
      if do.([lower(types{i}),'contour'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawcontours'',%d,''%s'')',panel,types{i}));
      end
    end
    
    %Interpolation
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
    
    for c = 'cgymkwrb'
      if do.([c,'roi'])
        addtodrawlist(panel,sprintf('drawfunctions(''drawroi'',%d,''%s'')',panel,c));
      end
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
    
    %types = {'r','b','g','s'};
    %for i = 1:length(types)
    %  if do.([types{i} 'contour'])
    %    DATA.drawlist{panel}{end+1} = sprintf('drawfunctions(''draw3dpoutline'',%d,''%s'')',panel,types{i});
    %  end      
    %end
    if do.threedpcontour
      addtodrawlist(panel,sprintf('drawfunctions(''draw3dpoutline'',%d,''threedp'')',panel));
    end
    
    if do.point
      addtodrawlist(panel,sprintf('drawfunctions(''drawpoint'',%d)',panel));
    end
    
end

%---------------------------------------------------------
function do = annotationinpanel(panel) %#ok<DEFNU>
%-----------------------------------------------------------
%this function sets graphics objects to nan in panel if there is no
%occurence in that slice of the annotation pertaining to that object.
%Returns struct with booleans indicating occurence of annotation.
global DATA SET

%this is done on panel level so if no segmentatio found in panel then it is
%not displayed in that panel;
slices = viewfunctions('slicesinpanel',panel);
no = DATA.ViewPanels(panel);

%Assume all exists then turn off later in code. These fields are all
%lowercase and contours have the field name which matches the SET field
%name but lower case and contour after.
do.endocontour = 1;
do.epicontour = 1;
do.rvendocontour = 1;
do.rvepicontour = 1;
do.generalpencontour = 1;
do.endointerp = 1;
do.epiinterp = 1;
do.rvendointerp = 1;
do.rvepiinterp = 1;
do.generalpeninterp = 1;
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

state = iconson('click3D');
if ~state{1}
  do.point3D = 0;
end

switch DATA.ViewPanelsType{panel}
  case {'hla','gla','vla'}
    %Somethings are never shown in hla gla vla
    do.endocontour = 0;
    do.epicontour = 0;
    do.rvendocontour = 0;
    do.rvepicontour = 0;
    do.generalpencontour = 0;
    do.endointerp = 0;
    do.epiinterp = 0;
    do.rvendointerp = 0;
    do.rvepiinterp = 0;
    do.generalpeninterp = 0;
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
      set(DATA.Handles.measurement(panel),'XData', nan, 'YData', nan);
      set(DATA.Handles.measurementoutsideplane(panel),'XData', nan, 'YData', nan);   
      set(DATA.Handles.measurementtext(panel,:),'position',[nan nan]);
    end
    
  case {'trans3DP','sag3DP','cor3DP','speedim'}
    %Somethings are never shown in 3DP views
    do.endocontour = 0;
    do.epicontour = 0;
    do.rvendocontour = 0;
    do.rvepicontour = 0;
    do.generalpencontour = 0;
    do.endointerp = 0;
    do.epiinterp = 0;
    do.rvendointerp = 0;
    do.rvepiinterp = 0;
    do.generalpeninterp = 0;
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
    
    %get panels shown and do the intersection lines based on that
    panels = ismember({'trans3DP','sag3DP','cor3DP','speedim'},DATA.ViewPanelsType);
    
    %Check which combinations there are
%     if length(DATA.ViewPanelsType)>2
%       
%       do.rgline = panels(1)*panels(2);
%       do.rbline = panels(1)*panels(3);
%       do.bgline = panels(3)*panels(2);
%       do.brline = panels(3)*panels(1);
%       do.gbline = panels(2)*panels(3);
%       do.grline = panels(2)*panels(1);
%     end
    
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
        %currenttype = 'r';
        currentind = 1;
      case 'sag3DP'
        do.rgline = 0;
        do.rbline = 0;
        do.grline = 1;
        do.gbline = 1;
        do.brline = 0;
        do.bgline = 0;        
        %currenttype = 'g';
        currentind = 2;
      case 'cor3DP'
        do.rgline = 0;
        do.rbline = 0;
        do.grline = 0;
        do.gbline = 0;
        do.brline = 1;
        do.bgline = 1;        
        
        %currenttype = 'b';
        currentind = 3;
    end
    
    %for the contours we need to check the show contour button
    outlineon = findindented(DATA.Handles.configiconholder,'outline');
        
    %default is off
    do.threedpcontour = 0;
    %do.rcontour = 0;
    %do.bcontour = 0;
    %do.gcontour= 0;
    %do.scontour= 0;
    
    %do.([currenttype 'contour']) = outlineon*panels(currentind);
    do.threedpcontour = outlineon*panels(currentind);
    
    types = {'rg','gr','rb','br','bg','gb'};
    for i =  1:length(types)
      if ~do.([types{i},'line'])
        DATA.Handles.([types{i},'line'])(panel).XData = [nan nan];
        DATA.Handles.([types{i},'line'])(panel).YData = [nan nan];
      end
    end
    
    %if ~do.([currenttype,'contour'])
    if ~do.threedpcontour
      DATA.Handles.threedpcontour(panel).ZData = [127 0 ; 0 0];
    end
        
    %points in 3d
    %points---------------------------------------------------
    switch DATA.ViewPanelsType{panel}
      case 'trans3DP'
        [~,ind] = findfunctions('closestpoint3dp',panel,1,1,SET(no).LevelSet.View.RSlice);
      case 'sag3DP'
        [~,ind] = findfunctions('closestpoint3dp',panel,1,1,SET(no).LevelSet.View.GSlice);
      case 'cor3DP'
        [~,ind] = findfunctions('closestpoint3dp',panel,1,1,SET(no).LevelSet.View.BSlice);
      otherwise
        ind = [];
    end
      
    if isempty(ind)
      do.point=0;
      set(DATA.Handles.point(panel),'XData', nan, 'YData', nan);
      set(DATA.Handles.pointtext(panel,:),'position',[nan nan]);
    end
    
  case 'viewport'
    %the viewport panel has no annotations whatsoever.
    do.endocontour = 0;
    do.epicontour = 0;
    do.rvendocontour = 0;
    do.rvepicontour = 0;
    do.generalpencontour = 0;
    do.endointerp = 0;
    do.epiinterp = 0;
    do.rvendointerp = 0;
    do.rvepiinterp = 0;
    do.generalpeninterp = 0;
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
    %do.rcontour = 0;
    %do.bcontour = 0;
    %do.gcontour= 0;
    %do.scontour= 0;
    
  otherwise
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
    
    
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    %contours-------------------------------------------------
    for i = 1:length(types)
      ind = findfunctions('findframeswithsegmentation',types{i},no,slices);
      
      if ~any(ind)
        do.([lower(types{i}),'contour']) = 0;
        set(DATA.Handles.([lower(types{i}),'contour'])(panel), ...
          'XData', nan, 'YData', nan)
      end
    end
    
    types = {'Endo','Epi','RVEndo','RVEpi','GeneralPen'};
    %interpolationpoints--------------------------------------
    for i = 1:length(types)
      ind = findfunctions('findframeswithinterpolationpoints',types{i},no,slices);
      
      if ~any(ind)
        do.([lower(types{i}),'interp']) = 0;
        set(DATA.Handles.([lower(types{i}),'interp'])(panel), ...
          'XData', nan, 'YData', nan)
      end
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
      if ~any(SET(linkednos(i)).RoiCurrent == find(roisinslice))
        set(DATA.Handles.roicurrent(panel),...
          'XData', nan, 'YData', nan);
        set(DATA.Handles.roitext(panel,:),'Position',[nan nan]);
      end
    end
    
    colorlist = 'cgykmwrb';
    notincludedcolors = colorlist(~ismember(colorlist,colors));
    
    for c = notincludedcolors
      set(DATA.Handles.([c,'roi'])(panel),...
        'XData', nan, 'YData', nan);
% % %       set(DATA.Handles.([c,'roitext'])(panel,:),...
% % %         'Position', [nan,nan]);
%       set(DATA.Handles.roitext(panel,:),...
%         'Position', [nan,nan]);
      do.([c,'roi']) = 0;
    end
    
    %points---------------------------------------------------
    if ~any(ismember(SET(no).Point.Z,slices))
      do.point=0;
      set(DATA.Handles.point(panel),'XData', nan, 'YData', nan);
      set(DATA.Handles.pointtext(panel,:),'position',[nan nan]);
    end
    
    %Scar------------------------------------------------------
    %For the contour objects the equivalent of nan rendering is placing a zero matix with a one in a corner.
%     visibleslices = viewfunctions('slicesinpanel', panel)
    if isempty(SET(no).Scar)
      tmp = [1,0;0,0];
      set([DATA.Handles.scarcontour(panel),DATA.Handles.weightedscarcontour(panel),...
           DATA.Handles.moextentcontour(panel),DATA.Handles.mocontour(panel)],'ZData',tmp,'YData',tmp,'XData',tmp);
      do.scar = 0;
      do.manualinteraction = 0;
    elseif ~any(SET(no).Scar.Result(:,:,slices),'all')
      tmp = [1,0;0,0];
      set([DATA.Handles.scarcontour(panel),DATA.Handles.weightedscarcontour(panel),...
           DATA.Handles.moextentcontour(panel),DATA.Handles.mocontour(panel)],'ZData',tmp,'YData',tmp,'XData',tmp);
      do.scar = 0;
      do.manualinteraction = 0;
    elseif ~any(SET(no).Scar.NoReflow(:,:,slices),'all')
      tmp = [1,0;0,0];
      set([DATA.Handles.moextentcontour(panel),DATA.Handles.mocontour(panel)],'ZData',tmp,'YData',tmp,'XData',tmp);      
    end
    
    %Since manual interaction is incorporated into the panel
    %imagehandle we need to check the button state here
    stateandicon = iconson('hidescarmanual');
    if stateandicon{1}
      do.manualinteraction = 0;
    end
    
    %Mar--------------------------------------------------------
    if isempty(SET(no).MaR)
      %no one can see it anyway
      tmp = [1,0;0,0];
      set(DATA.Handles.marcontour(panel),'ZData',tmp,'YData',tmp,'XData',tmp);
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
      set(DATA.Handles.measurement(panel),'XData', nan, 'YData', nan);
      set(DATA.Handles.measurementoutsideplane(panel),'XData', nan, 'YData', nan);      
      set(DATA.Handles.measurementtext(panel,:),'position',[nan nan]);
    end
end


%---------------------
function updatevisibility %#ok<DEFNU>
%-----------------------
global DATA

%Turn all on first
set(DATA.Handles.marcontour,'Visible', 'on');
set(DATA.Handles.planeintersection,'Visible','on');
set(DATA.Handles.scarcontour,'Visible', 'on');
set(DATA.Handles.weightedscarcontour,'Visible', 'on'); %Weighted scar contour
set(DATA.Handles.moextentcontour,'Visible', 'on');
set(DATA.Handles.mocontour,'Visible', 'on');
set(DATA.Handles.endocontourintersection,'Visible', 'on');
set(DATA.Handles.epicontourintersection,'Visible', 'on');
set(DATA.Handles.rvendocontourintersection,'Visible', 'on');
set(DATA.Handles.rvepicontourintersection,'Visible', 'on');
set([DATA.Handles.endointerp, ...
  DATA.Handles.epiinterp, ...
  DATA.Handles.rvendointerp, ...
  DATA.Handles.rvepiinterp,...
  DATA.Handles.endocontour, ...
  DATA.Handles.epicontour,...
  DATA.Handles.rvendocontour, ...
  DATA.Handles.rvepicontour],'Visible', 'on');
set(DATA.Handles.roicurrent,'Visible','on');
set([DATA.Handles.croi,...
  DATA.Handles.rroi,...
  DATA.Handles.mroi,...
  DATA.Handles.kroi,...
  DATA.Handles.groi,...
  DATA.Handles.wroi,...
  DATA.Handles.yroi,...
  DATA.Handles.broi],...
  'visible','on')
set(DATA.Handles.roitext(:),'Visible','on');
% set([DATA.Handles.croitext(:),...
%   DATA.Handles.rroitext(:),...
%   DATA.Handles.mroitext(:),...
%   DATA.Handles.kroitext(:),...
%   DATA.Handles.groitext(:),...
%   DATA.Handles.wroitext(:),...
%   DATA.Handles.yroitext(:),...
%   DATA.Handles.broitext(:)],...
%   'Visible','on');
%set(DATA.Handles.roitext,'Visible','on');
set(DATA.Handles.measurement,'Visible','on');
set(DATA.Handles.measurementtext,'Visible','on');
set(DATA.Handles.point,'Visible','on');
set(DATA.Handles.point3D,'Visible','on');
set(DATA.Handles.pointtext,'Visible','on');
set(DATA.Handles.text,'Visible','on');
set(DATA.Handles.centercross,'Visible', 'on');

%The do structures contains information of which updates to use.
tmp = cell(1,length(DATA.ViewPanels));
for i  = 1:length(tmp)
  tmp{i} = 1;
end

do = struct('endocontour',tmp,'epicontour',tmp,'rvendocontour',tmp,'rvepicontour',tmp,'generalpencontour',tmp...
  ,'endointerp',tmp,'epiinterp',tmp,'rvendointerp',tmp,'rvepiinterp',tmp,'generalpeninterp',tmp,...
  'measures',tmp,'point',tmp,'point3D',tmp,'croi',tmp,'mroi',tmp,'groi',tmp,'rroi',tmp,...
  'broi',tmp,'wroi',tmp,'yroi',tmp,'kroi',tmp,'text',tmp,...
  'intersections',tmp,'scar',tmp,'mar',tmp,'manualinteraction',tmp,...
  'rgline',tmp,'bgline',tmp,'rbline',tmp,'gbline',tmp,'brline',tmp,...
  'grline',tmp,'threedpcontour',tmp); 

%Then we perform a check if there exists objects in any of the timeframes
%for contours interpolation points measures rois viability %this is done on panel level so if no segmentatio found in panel then it is
%not displayed in that panel;
for p = find(DATA.ViewPanels)
  do(p) = viewfunctions('annotationinpanel',p);
end

%Then hide icons are incorporated. The hide state is set for all panels so
%if a state is hidden then that is applied for all panels.
hidecell={'hideplus','hidemar','hidescar',...
  'hideintersections','hideothercontour',...
  'hideinterp','hidelv','hiderv','hideroi',...
  'hidemeasure','hidepoint','hidetext','hidescarmanual'};%'hidescarextent','hidescarmanual'};

stateandicon = iconson(hidecell);
availableicons = cellfun(@(x) isa(x,'myicon'),stateandicon(:,2));
% state=[stateandicon{availableicons,1}];
state = [stateandicon{:,1}];

% First state corresponds to centercross
if state(1)
  set(DATA.Handles.centercross,'Visible', 'off');
end

%Second state corresponds to the marcontour.
if state(2)
  set(DATA.Handles.marcontour,'Visible', 'off');
  %remake images without manual interaction
  DATA.ViewIM = cell(size(DATA.ViewIM));
end

%Third state corresponds to the scarcontour.
if state(3)
  set(DATA.Handles.scarcontour,'Visible', 'off');
  set(DATA.Handles.weightedscarcontour,'Visible', 'off'); %Weighted scar contour
  set(DATA.Handles.moextentcontour,'Visible', 'off');
  set(DATA.Handles.mocontour,'Visible', 'off');
  %remake images without manual interaction visually turned off at bottom
  %of function
  DATA.ViewIM = cell(size(DATA.ViewIM));
  [do.scar] = deal(0);
  [do.manualinteraction] = deal(0);
end

%Fourth state corresponds to the plane intersections.
if state(4)
  set(DATA.Handles.planeintersection,'Visible','off');
end

%Fifth state corresponds to the contour intersections.
if state(5)
  set(DATA.Handles.endocontourintersection,'Visible', 'off');
  set(DATA.Handles.epicontourintersection,'Visible', 'off');
  set(DATA.Handles.rvendocontourintersection,'Visible', 'off');
  set(DATA.Handles.rvepicontourintersection,'Visible', 'off');
  [do.intersections] = deal(0);
end

%Sixth state corresponds to the interpolationpoints.
if state(6)
  iph = mycat(2,DATA.Handles.endointerp, ...
    DATA.Handles.epiinterp, ...
    DATA.Handles.rvendointerp, ...
    DATA.Handles.rvepiinterp);
  set(iph,'Visible', 'off');
  [do.interp] = deal(0);
end

%Seventh state corresponds to LV contours
if state(7)
  lvh = mycat(2, DATA.Handles.endocontour, ...
    DATA.Handles.epicontour,...
    DATA.Handles.endointerp, ...
    DATA.Handles.epiinterp, ...
    DATA.Handles.endocontourintersection, ...
    DATA.Handles.epicontourintersection);
  set(lvh,'Visible', 'off');
  [do.intersections] = deal(0);
end

%Eigth state corresponds to RV contours
if state(8)
  rvh = mycat(2, DATA.Handles.rvendocontour, ...
    DATA.Handles.rvepicontour,...
    DATA.Handles.rvendointerp, ...
    DATA.Handles.rvepiinterp, ...
    DATA.Handles.rvendocontourintersection, ...
    DATA.Handles.rvendocontourintersection);
  set(rvh,'Visible', 'off');
  [do.intersections] = deal(0);
end

%nineth state corresponds to Rois
if state(9)
  set(DATA.Handles.roicurrent,'Visible','off');
  set([DATA.Handles.croi,...
    DATA.Handles.rroi,...
    DATA.Handles.mroi,...
    DATA.Handles.kroi,...
    DATA.Handles.groi,...
    DATA.Handles.wroi,...
    DATA.Handles.yroi,...
    DATA.Handles.broi],...
    'visible','off')
% % %   set([DATA.Handles.croitext(:),...
% % %     DATA.Handles.rroitext(:),...
% % %     DATA.Handles.mroitext(:),...
% % %     DATA.Handles.kroitext(:),...
% % %     DATA.Handles.groitext(:),...
% % %     DATA.Handles.wroitext(:),...
% % %     DATA.Handles.yroitext(:),...
% % %     DATA.Handles.broitext(:)],...
% % %     'Visible','off');
  set(DATA.Handles.roitext(:),'Visible','off');
  [do.roi] = deal(0);
end

%Tenth state corresponds to measures
if state(10)
  set(DATA.Handles.measurement,'Visible','off');
  set(DATA.Handles.measurementtext,'Visible','off');
  [do.measures] = deal(0);
end

%Eleventh state corresponds to points
if state(11)
  set(DATA.Handles.point,'Visible','off');
  set(DATA.Handles.pointtext,'Visible','off');
  [do.point] = deal(0);
end

%Twelveth state corresponds to text
if state(12)
  set(DATA.Handles.text,'Visible','off');
  set(DATA.Handles.pointtext,'Visible','off');
  set(DATA.Handles.measurementtext,'Visible','off');
  set(DATA.Handles.roitext(:),'Visible','off');
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
function selectroi(panel) %#ok<DEFNU>
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

oldcurrentroi = SET(magno).RoiCurrent;
if ~isempty(currentroi)
  SET(magno).RoiCurrent = currentroi;
  
  for p = find(ismember(DATA.ViewPanels,SET(no).Linked))
    drawfunctions('drawroi',p);%,SET(magno).Roi(currentroi).LineSpec(1))
  end
% % %   for p = find(ismember(DATA.ViewPanels,SET(no).Linked))
% % %     drawfunctions('drawroi',p,SET(magno).Roi(oldcurrentroi).LineSpec(1))
% % %   end
  
  segment('updateflow')
end

%---------------------
function setviewtype(viewpanelstype)
%-----------------------
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

%Assert correct button is indented
DATA.setviewbuttons(0)
if strcmp(DATA.ProgramName,'Segment') 
  segment('updatevolume');
elseif contains(DATA.ProgramName,'CMR')
  % update volume axes
  DATA.updatevolumeaxes;
end
segment('updateflow');
segment('updatemeasurement');
DATA.updatetimebaraxes


%---------------------
function switchimagestack(no,viewpanelstype) %#ok<DEFNU>
%-----------------------
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

%-----------------------------
function speedim = findspeedim 
%-----------------------------
%Returns panel for speedim otherwise empty

global DATA

speedim = find(strcmp(DATA.ViewPanelsType,'speedim'));

%---------------------------------------------
function zoomstate = getnewzoomstate(panel,no, resetflag)
%---------------------------------------------
%Compute initialization of zoom state

global DATA SET
if nargin<3
  resetflag = true;
end
switch DATA.ViewPanelsType{panel}
  case 'trans3DP'
    color = 'r';
    [xsize,ysize,~,~] = segment3dp.tools('viewsizeandres',color,no);
  case 'sag3DP'
    color = 'g';
    [xsize,ysize,~,~] = segment3dp.tools('viewsizeandres',color,no);
  case 'cor3DP'
    color = 'b';
    [xsize,ysize,~,~] = segment3dp.tools('viewsizeandres',color,no);
  case 'speedim'
    if isempty(SET(no).LevelSet) %~isfield(SET(no).LevelSet.Pen,'Color')
       color= 'r';
    else
       color = SET(no).LevelSet.Pen.Color;
    end
    [xsize,ysize,~,~] = segment3dp.tools('viewsizeandres',color,no);
  otherwise
  xsize = SET(no).XSize;
  ysize = SET(no).YSize;
end

%Find panelsize
panelpos = DATA.Handles.imageaxes(panel).Position;
figpos = DATA.fig.Position;
panelwidth = panelpos(3)*figpos(3)/ysize;
panelheight = panelpos(4)*figpos(4)/xsize;

%maxsize = max(SET(no).XSize,SET(no).YSize);
minpanelsize = min(panelwidth,panelheight);
panelwidth = panelwidth/minpanelsize; %range 1..
panelheight = panelheight/minpanelsize; %range 1..

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
  %check if 3dp
  is3dp = ismember(DATA.ViewPanelsType{panel},{'trans3DP','sag3DP','cor3DP','speedim'});
    
  if is3dp
    color = SET(no).LevelSet.Pen.Color;
    zoomstate = SET(no).LevelSet.View.([upper(color) 'ZoomState']); %get zoomstate
  else
    zoomstate = SET(no).NormalZoomState; %get zoomstate
  end
  
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
  f = 1-incr/7;
  
  %Compute newzoom state based on the center and how much to extend
  newzoomstate = [xcenter-xsize/2*f xcenter+xsize/2*f ycenter-ysize/2*f ycenter+ysize/2*f];

  if ~is3dp    
    %assign
    if all(noslinked ~= 0) && ~all(nos == noslinked)
      % set also the invisble nos to the same Zoomstate
      ind = nos ~= noslinked;
      numnos = length(ind);
      if numnos > 1
        for cl = 2:numnos
          SET(ind(cl)).NormalZoomState = [];
        end
      end
    end
    SET(no).NormalZoomState = newzoomstate;

    %update it graphically
    updatezoomandaspectratio(panel);
    drawfunctions('drawplaneintersections');
    %we want to update the text position and the frame
    drawfunctions('drawselectedframe',panel)
    updatetextposition(panel)
  else
    %--- special for 3DP
    
    %assign
    SET(no).LevelSet.View.([upper(color) 'ZoomState']) = newzoomstate;

    speedimpanel = findspeedim;
    
    %find color panel
    color = SET(no).LevelSet.Pen.Color;
    switch color
      case 'r'
        panelname = 'trans3DP';
      case 'g'        
        panelname = 'sag3DP';
      case 'b'        
        panelname = 'cor3DP';
    end
    
    colorpanel = find(strcmp(DATA.ViewPanelsType,panelname));

    %update it graphically
    if ~isempty(speedimpanel)
      updatezoomandaspectratio(speedimpanel);
      drawfunctions('drawselectedframe',speedimpanel)
      updatetextposition(speedimpanel)    
    end
    if ~isempty(colorpanel)
      updatezoomandaspectratio(colorpanel);
      drawfunctions('drawselectedframe',colorpanel)
      updatetextposition(colorpanel)    
    end
  end
end

%---------------------
function updatetextposition(panel)
%-----------------------
%This function updates the text position after updates such as zooming or
%panning

global DATA

x = max(0.5,0.015*(DATA.Handles.imageaxes(panel).XLim(2)-DATA.Handles.imageaxes(panel).XLim(1)));
y = max(0.5,0.015*(DATA.Handles.imageaxes(panel).YLim(2)-DATA.Handles.imageaxes(panel).YLim(1)));
DATA.Handles.text(panel).Position = [DATA.Handles.imageaxes(panel).XLim(1)+x,DATA.Handles.imageaxes(panel).YLim(1)+y];

%--------------------------
function switchpanel(panel)
%--------------------------
%Toggles the CurrentPanel and NO to be correct aswell as does the necessary
%graphical updates. If trying to toggle to an empty panel dont change
%anything

global DATA SET NO

no = DATA.ViewPanels(panel);

%this is the empty panel case
if no ~= 0
  DATA.CurrentPanel = panel;
  NO = no;
  drawfunctions('drawselectedframe',panel)
  drawfunctions('drawthumbnailframes')
end

%this draws and updates the panel image intersections.
drawfunctions('drawplaneintersections')

%checks if switching 3DP panel to 3DP panel if so we need to update the
%current color and speedimage
if any(strcmp(DATA.ViewPanelsType(DATA.CurrentPanel),{'trans3DP','sag3DP','cor3DP'}))
  
  %Check if switching to new color
  switch DATA.ViewPanelsType{panel}
    case 'trans3DP'
      newcolor = 'r';
    case 'sag3DP'
      newcolor = 'g';
    case 'cor3DP'
      newcolor = 'b';
  end
  
  %Store new color
  if SET(no).LevelSet.Pen.Color ~= newcolor
    SET(no).LevelSet.Pen.Color = newcolor;
    speedpanel = find(strcmp(DATA.ViewPanelsType,'speedim'));
    if ~isempty(speedpanel)
            
      DATA.ViewIM{speedpanel} = [];
      %drawfunctions('drawimages',speedpanel);    
      drawfunctions('drawpanel',speedpanel)
     
      %Remove all lines
      h = [...
        DATA.Handles.rgline(speedpanel) ...
        DATA.Handles.rbline(speedpanel) ...
        DATA.Handles.grline(speedpanel) ...
        DATA.Handles.gbline(speedpanel) ...
        DATA.Handles.brline(speedpanel) ...
        DATA.Handles.bgline(speedpanel)];
  
      set(h,'xdata',NaN,'ydata',NaN);

      %add new liness
      switch newcolor
        case 'r'
          drawfunctions('draw3dpline',speedpanel,'rg');          
          drawfunctions('draw3dpline',speedpanel,'rb');
        case 'g'
          drawfunctions('draw3dpline',speedpanel,'gr');
          drawfunctions('draw3dpline',speedpanel,'gb');
        case 'b'
          drawfunctions('draw3dpline',speedpanel,'br');          
          drawfunctions('draw3dpline',speedpanel,'bg');          
      end
      
      %Update text position on speedpanel
      viewfunctions('updatetextposition',speedpanel)
    end
  end
 
end

DATA.setviewbuttons(0);
if strcmp(DATA.ProgramName,'Segment')
  segment('updatevolume');  
end
segment('updateflow'); 
segment('updatemeasurement');
DATA.updatetimebaraxes;

%-------------------------------------
function placecontextmenu(panel,type,objectind) %#ok<DEFNU>
%-------------------------------------
%Function that returns the clicked slice in montage and montage
global DATA SET NO

[p(1),p(2)] = mygetcurrentpoint(DATA.fig);

switch type
  case 'Measure'
    set(DATA.Handles.measurecontextmenu,...
      'Position',p,'Visible','on');
  case 'Point'
    set(DATA.Handles.pointcontextmenu,...
      'Position',p,'Visible','on');
  case {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'}
    set(DATA.Handles.interppointmenu,...
      'Position',p,'Visible','on');
  case 'Roi'
    % set current roi to the clicked object in (objectind)
    SET(NO).RoiCurrent = objectind;
    drawfunctions('drawroi',panel);
    set(DATA.Handles.roicontextmenu,...
      'Position',p,'Visible','on');
  otherwise
    switch DATA.CurrentTool
      case {'Endo','Epi',...
          'EndoInterp','EpiInterp','RVEndo','RVEpi','Scar',...
          'MO','MaR','ScarRubber','MORubber','MaRRubber',...
          'RVEndoInterp','RVEpiInterp','GeneralPenInterp',...
          'GeneralPen','EndoBalloon','Contrast'}
        set(DATA.Handles.lvcontextmenu,...
          'Position',p,'Visible','on');
        
      otherwise
        set(DATA.Handles.selectcontextmenu,...
          'Position',p,'Visible','on');
    end
end

%------------------------------------------------------------------------
function slice = clickedslice(panel,x,y) %#ok<DEFNU>
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
  x_ind = ceil(x/scale/SET(no).YSize);
  y_ind = ceil(y/scale/SET(no).XSize);
  
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

%---------------------
function setview(rows,cols,nos,viewpanelstype,currentnopanel)
%-----------------------
%Creates imageaxes and boxaxes. The field boxaxes is probably unnecessary.

global DATA NO SET

if DATA.Silent 
  return
end

%Set information fields
%profile on
%this is the refresh option finds which viewpanelbutton is indented and
%sets rows and cols
if nargin == 0 || (isempty(rows) && isempty(cols))
  rows = DATA.ViewMatrix(1);
  cols = DATA.ViewMatrix(2);
 
  panelnos = DATA.ViewPanels(DATA.ViewPanels>0);
  for loop = panelnos
    SET(loop).NormalZoomState=[];
  end
  if ~isempty(SET(NO).LevelSet)  
    SET(NO).LevelSet.View.RZoomState = [];
    SET(NO).LevelSet.View.GZoomState = [];
    SET(NO).LevelSet.View.BZoomState = [];
  end
  
else
  DATA.ViewMatrix(1) = rows;
  DATA.ViewMatrix(2) = cols;
end

%Take care of outline contours when switching from 3PD to other
%if DATA.Handles.configiconholder.findindented('outline')
%  for loop = 1:length(DATA.ViewPanels)
%    if ismember(DATA.ViewPanelsType{loop},{'trans3DP','sag3DP','cor3DP','speedim'})
%      drawfunctions('setxynan',loop);
%    end
%  end
%end

%Delete colorbar
paneltoclean = find(DATA.ViewPanels>0);

for i = paneltoclean
  colorbar(DATA.Handles.imageaxes(i),'off')
end
specialview = false;
%We need to handle different input cases and the chance that some panels
%are empty. also the current panel might be out of bounds when toggling to
%a lesser amount of panels. Then we set the CurrentPanel and NO field to
%the first panel content using the switchpanel command
if nargin < 3
  %This is the case when we toggle to a panel view from orthoview
  if any(ismember(DATA.ViewPanelsType,{'orth','hla','gla','vla'})) && nargin~=0 %if refresh keep ortho view
    DATA.ViewPanels=zeros(1,rows*cols);
    DATA.ViewPanels(1)=NO;
    DATA.ViewPanelsType{1}='one';
    DATA.ViewPanelsType(2:end)=[];
    DATA.CurrentPanel = 1;
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
  DATA.CurrentPanel = currentnopanel(2);
else
  if DATA.CurrentPanel > rows*cols
    DATA.ViewPanels(rows*cols) = NO;
    DATA.CurrentPanel = rows*cols;
  end
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
    for i = 1:length(DATA.ViewPanelsType(:))
      tmp{i} = DATA.ViewPanelsType{i};
    end
    DATA.ViewPanelsType = tmp;
  else
    DATA.ViewPanelsType = DATA.ViewPanelsType(1:rows*cols);
  end
else
  DATA.ViewPanelsType = viewpanelstype;
end
for typeloop = 1:length(DATA.ViewPanelsType) 
  emptytype(typeloop) = isempty(DATA.ViewPanelsType{typeloop}); 
end
if isempty(DATA.ViewPanelsType) || any(emptytype)
  DATA.ViewPanelsType = cell(1,length(DATA.ViewPanels));
  DATA.ViewPanelsType(:) = {'one'};
end

% Here we set so that the viewpanelsmatrix field is correct it is [1,1] if
% not montage
DATA.ViewPanelsMatrix = cell(1,length(DATA.ViewPanels));
for i = 1:numel(DATA.ViewPanels)
  if DATA.ViewPanels(i) > 0 && any(strcmp(DATA.ViewPanelsType{i},{'montage','montagerow','montagesegmented'}))
    if strcmp(DATA.ViewPanelsType{i},'montagesegmented')
      slicestoinclude = segment('getmontagesegmentedslices',DATA.ViewPanels(i));
      [mrows,mcols] = calcfunctions('calcrowscols',DATA.ViewPanels(i),length(slicestoinclude));
    else
      [mrows,mcols] = calcfunctions('calcrowscols',DATA.ViewPanels(i));
    end
    if ismember(DATA.ViewPanelsType{i},{'montagerow'})
      m = mrows*mcols;
      mrows = min(mrows,2);
      mcols = ceil(m/mrows);
    end
    DATA.ViewPanelsMatrix{i} = [mcols mrows];
  else
    DATA.ViewPanelsMatrix{i} = [1,1];
  end
end

if ~specialview
    %functionality to restore same zoom state as before panel settings were
  %changed
  panelstodo = find(DATA.ViewPanels(DATA.ViewPanels>0));
  incr = zeros(length(panelstodo),1);

  for actpanel = panelstodo
    no = DATA.ViewPanels(actpanel);
    linkednos = SET(no).Linked;
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
  
%if we are dealing with 3dp views things are done slightly different.
if any(ismember(DATA.ViewPanelsType,{'trans3DP','sag3DP','cor3DP','speedim','viewport'}))
       
  %check 3dp fields asserting everything is initiated
  segment3dp.tools('check3dpfields')
  
  %Find suitable penradius
  DATA.LevelSet.maxpenradius = min(11*[SET(NO).ResolutionX SET(NO).ResolutionY ...
    SET(NO).SliceThickness+SET(NO).SliceGap]);
  SET(NO).LevelSet.Pen.Radius = min(DATA.LevelSet.maxpenradius, SET(NO).LevelSet.Pen.Radius);
  set(DATA.Handles.penradiusedit,'string',sprintf('%0.1g',SET(NO).LevelSet.Pen.Radius));

  %Set edit boxes
  segment3dp.tools('setpredefined');

  %Scale axes
  g = linspace(0,1,128)';
  DATA.LevelSet.colormap = [[g g g] ; ...
    [ones(128,1) linspace(0,1,128)' zeros(128,1)]];
  image(reshape(flipud(DATA.LevelSet.colormap),[size(DATA.LevelSet.colormap,1) 1 3]),'parent',DATA.Handles.scaleaxes);
  axis(DATA.Handles.scaleaxes,'off');
  cmap = gray(256);
  image(reshape(flipud(cmap),[size(cmap,1) 1 3]),'parent',DATA.Handles.intensityscaleaxes);
  axis(DATA.Handles.intensityscaleaxes,'off');
  drawfunctions('drawmapping');
  segment3dp.tools('updatesloperadiobuttons');
  drawfunctions('updatemapping');
  segment3dp.tools('updatemappingedit');
  segment3dp.tools('updateobjects');
  drawfunctions('drawintensitymapping');
  
  for loop = 1:length(DATA.ViewPanels)
    if ~isequal(DATA.ViewPanelsType{loop},'speedim')
      drawfunctions('drawtext',loop);
    end
  end  
  
  %Set current pen color
  switch DATA.ViewPanelsType{DATA.CurrentPanel}
    case 'trans3DP'
      SET(DATA.ViewPanels(DATA.CurrentPanel)).LevelSet.Pen.Color = 'r';
    case 'sag3DP'
      SET(DATA.ViewPanels(DATA.CurrentPanel)).LevelSet.Pen.Color = 'g';
    case 'cor3DP'
      SET(DATA.ViewPanels(DATA.CurrentPanel)).LevelSet.Pen.Color = 'b';    
  end
  
else
  %this draws and updates the panel image intersections.
  drawfunctions('drawplaneintersections')  
end
%When updating we only want to update the panels with content. These are
%DATA.ViewPanels>0
nostodo = DATA.ViewPanels(DATA.ViewPanels>0);
panelstodo = find(nostodo);

for i = panelstodo  
  
  updatedrawlist(i)
  DATA.ViewIM{i} = []; 
  
  drawfunctions('drawpanel',i)
  if ~specialview
    if incr(i)~= 0
        % restore zoom
        SET(nostodo(i)).NormalZoomState = [];
        zoom(incr(i),i);
    end
  end
  updatetextposition(i);
  createfunctions('addcolorbar',i)
  %view sliders
  if ismember(DATA.ViewPanelsType{i},{'trans3DP','sag3DP','cor3DP','viewport','speedim'})  

    %Fix position
    set(DATA.Handles.sliders(i),'units','normalized'); %needs to be set before, Matlab bug?
    p = get(DATA.Handles.imageaxes(i),'Position');
    set(DATA.Handles.sliders(i),'Position',[p(1)+p(3)-0.01 p(2) 0.01 p(4)],'Visible','on');

    %Fix min/max & value
    [rmax,gmax,bmax] = segment3dp.tools('xyz2rgb',...
      SET(NO).XSize,...
      SET(NO).YSize,...
      SET(NO).ZSize);    
    switch DATA.ViewPanelsType{i}
      case 'trans3DP'
        set(DATA.Handles.sliders(i),...
          'min',1,...
          'max',rmax,...
          'value',SET(NO).LevelSet.View.RSlice,...
          'SliderStep',...
          [1/rmax 0.1],...
          'Callback',sprintf('segment3dp.tools(''rgbslider_Callback'',''r'',%d)',i));
      case 'sag3DP'        
        set(DATA.Handles.sliders(i),...
          'min',1,....
          'max',gmax,...
          'SliderStep',[1/gmax 0.1],...
          'value',SET(NO).LevelSet.View.GSlice,...
          'Callback',sprintf('segment3dp.tools(''rgbslider_Callback'',''g'',%d)',i));
      case 'cor3DP'        
        set(DATA.Handles.sliders(i),...
          'min',1,...
          'max',bmax,...
          'SliderStep',[1/bmax 0.1],...
          'value',SET(NO).LevelSet.View.BSlice,...
          'Callback',sprintf('segment3dp.tools(''rgbslider_Callback'',''b'',%d)',i));
      case {'viewport','speedim'}
        if segment3dp.isviewportalive && DATA.Handles.configiconholder.findindented('showcutplane')
          set(DATA.Handles.sliders(i),'min',-50,'max',50,'visible','on','value',-SET(NO).LevelSet.View.CutPlaneOffset,'Callback',sprintf('segment3dp.tools(''cutplaneoffset_Callback'',%d)',i));
        end
    end

  else
    set(DATA.Handles.sliders(i),'Visible','off');    
  end
  
end

%segment3dp.tools('updateoffsetslider');

%Remove slider from all others
allpanels = 1:length(DATA.Handles.sliders);
for i = setdiff(allpanels,panelstodo)
  set(DATA.Handles.sliders(i),'Visible','off');
end
if specialview
    %functionality to restore same zoom state as before panel settings were
  %changed
  panelstodo = find(DATA.ViewPanels(DATA.ViewPanels>0));
  incr = zeros(length(panelstodo),1);

  for actpanel = panelstodo
    no = DATA.ViewPanels(actpanel);
    linkednos = SET(no).Linked;
    if length(linkednos) > 1
      SET(linkednos(2:end)).NormalZoomState = [];
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
  for i = panelstodo 
    if incr(i)~= 0
        % restore zoom
        SET(nostodo(i)).NormalZoomState = [];
        zoom(incr(i),i);
    end
  end
end
% end of restoring zoom
drawfunctions('drawselectedframe',DATA.CurrentPanel)
drawfunctions('drawselectedslice',DATA.CurrentPanel)
buttondownfunctions('updatebuttondowns')

DATA.setviewbuttons(0);
segment('updatevolume');
segment('updateflow');
segment('updatemeasurement');
DATA.updatetimebaraxes
%create thumbnails and draw frame around them
drawfunctions('drawthumbnails',isempty(DATA.DATASETPREVIEW));
%profile report

%---------------------
function configaxes
%-----------------------
%configure imageaxes and boxaxes according to current config.
global DATA

rows = DATA.ViewMatrix(1);
cols = DATA.ViewMatrix(2);

left = DATA.GUISettings.LeftGapWidth; %0.12;
right = 1-DATA.GUISettings.RightGapWidth-0.02;  %Based on that the report panel can never be more than 220.
bottom = DATA.GUISettings.BottomGapHeight; %0.013;
top = 1-DATA.GUISettings.TopGapHeight;
width = right-left;
height = top-bottom;

%bottom = bottom+0.1;
%height = height-0.1;

%Set the box axes position
DATA.Handles.boxaxes.Position = [left bottom width height];

%Set up the frame around all the panels
x = [];
y = [];
for rloop=1:(rows-1)
  x = [x nan 0 1];
  y = [y nan 1/rows*rloop 1/rows*rloop];
end

for cloop=1:(cols-1)
  x = [x nan 1/cols*cloop 1/cols*cloop];
  y = [y nan 0 1];
end

x = [x,[nan 0 1]];
y = [y,[nan 0 0]];
x = [x,[nan 0 1]];
y = [y,[nan 1 1]];
x = [x,[nan 0 0]];
y = [y,[nan 0 1]];
x = [x,[nan 1 1]];
y = [y,[nan 0 1]];

DATA.Handles.box.XData = x;
DATA.Handles.box.YData = y;

%Place all panels
counter = 1;
for rloop=(rows-1):-1:0
  for cloop=0:(cols-1)
    DATA.Handles.imageaxes(counter).Position = [left+cloop*width/cols bottom+rloop*height/rows width/cols height/rows];
    counter = counter + 1;
  end
end

%nan all handles and axes not shown
drawfunctions('setxynan')

%if viewport is up
if ~strcmp(DATA.CurrentTheme,'3dp') && ~isempty(DATA.LevelSet) && isfield(DATA.LevelSet,'ViewPort') && ~isempty(DATA.LevelSet.ViewPort) %&& DATA.Handles.configiconholder.findindented('view3d')
  try
    delete(DATA.LevelSet.ViewPort)
  catch
  end
  DATA.LevelSet.ViewPort = [];
end

%------------------------------
function im = getgla(no)
%------------------------------
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
function switchtimeframe(incr,synchronize) %#ok<DEFNU>
%-----------------------------------------
global DATA SET

panel = DATA.CurrentPanel;
no = DATA.ViewPanels(panel);
tools('connectinterpolation',no,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'});
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

if not(contains(DATA.ProgramName,'3DP')) && (findindented(DATA.Handles.hideiconholder,'synchronize')||synchronize)
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

%Loop over all image stacks to update time.
for loop=nos
  SET(loop).CurrentTimeFrame = min(max(1,1+round(percent*(SET(loop).TSize-1))),SET(loop).TSize);
end

for p = panels  
  viewfunctions('updatedrawlist',p)
  drawfunctions('drawpanel',p)
end

%also update the timebars!
viewfunctions('updatetimebars');

%---------------------
function updatetimebars(synchronize) %#ok<DEFNU>
%-----------------------
%updates all timebars
global DATA SET NO

no = NO;

if nargin == 0
  if not(contains(DATA.ProgramName,'3DP')) && findindented(DATA.Handles.hideiconholder,'synchronize')
    synchronize = 1;
  else
    synchronize = 0;
  end
end


if ~isempty(DATA.LVNO)|| strcmp(DATA.ProgramName,'Segment')
  if strcmp(DATA.ProgramName,'Segment')
    t = SET(no).TimeVector*1000;
  else
    t = SET(DATA.LVNO(1)).TimeVector*1000;
  end
  
  if  any(no == DATA.LVNO) || strcmp(DATA.ProgramName,'Segment') || (synchronize && any(DATA.ViewPanels == DATA.LVNO))
    
    if strcmp(DATA.ProgramName,'Segment')
      set(DATA.Handles.timebarlv,'xdata',...
        t(SET(no).CurrentTimeFrame)*[1 1])
    else
      set(DATA.Handles.timebarlv,'xdata',...
        t(SET(DATA.LVNO(1)).CurrentTimeFrame)*[1 1])
    end
    
  end
end




if ~isempty(DATA.FlowNO)
  if strcmp(DATA.ProgramName,'Segment')
    t = SET(no).TimeVector*1000;
  else
    t = SET(DATA.FlowNO(1)).TimeVector*1000;
  end
  
  if  no == DATA.FlowNO || strcmp(DATA.ProgramName,'Segment') || (synchronize && any(DATA.ViewPanels == DATA.FlowNO))
    if strcmp(DATA.ProgramName,'Segment')
      DATA.Handles.timebarflow.XData = t(SET(no).CurrentTimeFrame)*[1 1];
    else
      DATA.Handles.timebarflow.XData = t(SET(DATA.FlowNO).CurrentTimeFrame)*[1 1];
    end
  end
end

t = SET(no).TimeVector*1000;
set(DATA.Handles.timebar,'xdata',...
  t(SET(no).CurrentTimeFrame)*[1 1]);

%---------------------
function selectallslices %#ok<DEFNU>
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
function unselectallslices %#ok<DEFNU>
%-----------------------
%Unselects all slices in panel
global DATA SET NO
SET(NO).StartSlice = [];
SET(NO).EndSlice= [];
drawfunctions('drawselectedslice',DATA.CurrentPanel);
drawfunctions('drawtext',DATA.CurrentPanel);

%---------------------
function viewallimagestacks(viewpaneltype) %#ok<DEFNU>
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
function fastswitchslice(sliceincr,synchronize) %#ok<DEFNU>
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
  'StartDelay',1);
outim = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
im = calcfunctions('remapuint8viewim',outim,no);

DATA.Handles.imagehandles(panel).CData = imresize(im,scale);

for i = 2:length(DATA.drawlist{panel})
  eval(DATA.drawlist{panel}{i})
end

start(storetimer)

%---------------------
function switchslice(sliceincr,synchronize) %#ok<DEFNU>
%-----------------------
global DATA SET

panel = DATA.CurrentPanel;
no = DATA.ViewPanels(panel);
slice = SET(no).CurrentSlice+sliceincr;
tools('connectinterpolation',no,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})
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

if not(contains(DATA.ProgramName,'3DP')) && (findindented(DATA.Handles.hideiconholder,'synchronize') || synchronize)
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
  viewfunctions('updatedrawlist',p)
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
  
  if thisslice > SET(thisno).ZSize
    thisslice = [];
  end
  if thisslice < 1
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
function viewhideallspecial_Callback(varargin) %#ok<DEFNU>
%---------------------------------------
%all hide buttons
global DATA
stateandicon = iconson('hideall');
state = stateandicon{1};
hidecell = {'hideplus','hidemar','hidescar','hidescarextent','hidescarmanual',... 'hidepins',
  'hideintersections','hideothercontour','hideinterp','hidelv','hiderv','hideroi','hidemeasure','hidepoint','hidetext'};%,'colorbar'};
stateandicon = iconson(hidecell);
availableicons = find(cellfun(@(x) isa(x,'myicon'),stateandicon(:,2)))';
if state
  DATA.LastHideIconState = stateandicon;
  for i = availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay = icon.cdataIndent;
    icon.isindented = 1;
  end
else  
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
%---------------------------------------
function viewhideall_Callback(varargin) %#ok<DEFNU>
%---------------------------------------
%all hide buttons
global DATA
stateandicon=iconson('hideall');
state=stateandicon{1};
hidecell={'hideplus','hidemar','hidescar','hidescarextent','hidescarmanual',... 'hidepins',
  'hideintersections','hideothercontour','hideinterp','hidelv','hiderv','hideroi','hidemeasure','hidepoint','hidetext'};%,'colorbar'};
stateandicon = iconson(hidecell);
availableicons = find(cellfun(@(x) isa(x,'myicon'),stateandicon(:,2)))';
if state
  DATA.LastHideIconState = stateandicon;
  for i=availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay=icon.cdataIndent;
    icon.isindented=1;
  end
else
  for i=availableicons
    icon = stateandicon{i,2};
    icon.cdataDisplay=icon.cdata;
    icon.isindented=0;
  end
end
viewfunctions('updatevisibility');
DATA.Handles.configiconholder.render;

%-------------------------------------------
function state = iconson(name)
%------------------------------------------
%if given name of button return state if run with no input returns
%states of all buttons. if given cell with multiple button return icons
%in order of request.
global DATA

state=0;
if not(contains(DATA.ProgramName,'3D')) 
  icons=[DATA.Handles.permanenticonholder.iconCell{:},DATA.Icons.lviconcell{:},DATA.Icons.rviconcell{:},DATA.Icons.roiflowiconcell{:},DATA.Icons.viabilityiconcell{:},...
    DATA.Icons.analysisiconcell{:}, DATA.Icons.imageiconcell{:},DATA.Icons.hidecell{:}];
else
  icons=[DATA.Icons.lviconcell{:},DATA.Icons.printiconcell{:},DATA.Icons.imageiconcell{:}];
end
N=length(icons);

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
    end
  end
end

%--------------------------
function hide(names,states) %#ok<DEFNU>
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
function hide_Callback(varargin) %#ok<DEFNU>
%---------------------------------------
%Toggle visibility for the imageoverlay
viewfunctions('updatevisibility');
allhidden

%---------------------------------
function viewhidecolorbar_Callback %#ok<DEFNU>
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
function viewhidepap_Callback %#ok<DEFNU>
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
function viewhidemanualinteraction_Callback %#ok<DEFNU>
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
function viewinterp_Callback(val) %#ok<DEFNU>
%--------------------------------
%Update if the view interpolated view is selected / unselected

global DATA
if nargin < 1
  val = ~DATA.Pref.ViewInterpolated;
end
DATA.Pref.ViewInterpolated = val;
DATA.ViewIM = cell(size(DATA.ViewIM)); %This clear ViewIM and it is reconstructed
for p = find(DATA.ViewPanels)
  drawfunctions('drawpanel',p);
end
drawfunctions('drawplaneintersections');

%-----------------------------------
function viewportviewfrom(type) %#ok<DEFNU>
%-----------------------------------
%Update view direction in viewport depending on direction

global DATA SET NO

%Get object handle
viewport = DATA.LevelSet.ViewPort;

imageorientation = SET(NO).ImageOrientation;
imageorientation = [1 0 0 0 1 0];

zdir = cross(...
  imageorientation(1:3),...
  imageorientation(4:6));

Rxyz2rlapfh = [imageorientation(4:6)' imageorientation(1:3)' -zdir'];

Rrlapfh2xyz = inv(Rxyz2rlapfh);

%Get camera position
pos = viewport.getcameraposition;
  
%Check distance to center (origo) for camera
r = sqrt(sum(pos.*pos));
  
%Define in rlapfh coordinates
switch type
  case 'transversal'
    posrlapfh = [0 r 0];
    upvrlapfh = [0 0 1];
  case 'sagittal'
    posrlapfh = [-r 0 0];
    upvrlapfh = [0 0 1];
  case 'coronal'
    posrlapfh = [0 0 r];
    upvrlapfh = [-1 0 0];
end
  
%Convert to xyz coordinates
posxyz = Rrlapfh2xyz*(posrlapfh');
upvxyz = Rrlapfh2xyz*(upvrlapfh');

%Update camera position
viewport.setcamera(posxyz,upvxyz);

%-----------------------------------
function updatelinewidth %#ok<DEFNU>
%-----------------------------------
global DATA

linewidth = DATA.Pref.LineWidth;
if ~isempty(linewidth)
  colortypes = 'cgbrwkym';
  type = {'Endo','Epi','RVEndo','RVEpi'};
  for panel = 1: length(DATA.ViewPanels)
    DATA.Handles.roicurrent(panel).LineWidth = linewidth + 1;
    DATA.Handles.scarcontour(panel).LineWidth = linewidth;
    DATA.Handles.marcontour(panel).LineWidth = linewidth;
    for loop = 1: length(type)
      DATA.Handles.([lower(type{loop}),'contour'])(panel).LineWidth = linewidth;
    end
    for cind = 1:length(colortypes)
      DATA.Handles.([colortypes(cind),'roi'])(panel).LineWidth = linewidth;
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

if nargin < 3
  outsideelement = 0;
end
z = repmat(im(1).*outsideelement,DATA.ViewPanelsMatrix{panel}(2)*SET(no).XSize,DATA.ViewPanelsMatrix{panel}(1)*SET(no).YSize);
loop=1;
for slice=1:SET(no).ZSize
  c = 1+mod(loop-1,DATA.ViewPanelsMatrix{panel}(1));
  r = ceil(loop/DATA.ViewPanelsMatrix{panel}(1));
  z(...
    (1+(r-1)*SET(no).XSize):(r*SET(no).XSize),...
    (1+(c-1)*SET(no).YSize):(c*SET(no).YSize)) = im(:,:,slice);
  loop=loop+1;
end

%------------------------------------------------------
function updateflowresultstack  %#ok<DEFNU>
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
function updatetoolhidestate(currenttool)  %#ok<DEFNU>
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