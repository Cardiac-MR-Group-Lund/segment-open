function varargout = LAXselectslices(varargin)
%Select short-axis slices to perform strain analysis on

macro_helper(varargin{:});
if nargin < 1 || isempty(varargin{1})
  varargin{1} = 'init';
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%----------------
function init(taggroup) %#ok<DEFNU>
%----------------
%Initialize the GUI
global DATA SET

gui = mygui(fullfile('+straintagging','LAXselectslices.fig'));
DATA.GUI.LAXselectslices = gui;
gui.taggroup=taggroup;

twochamber=[];
threechamber=[];
fourchamber=[];
roof=max(taggroup);
xsz =zeros(1,roof);
ysz = zeros(1,roof);
xres = zeros(1,roof);
yres=zeros(1,roof);
zres=zeros(1,roof);
axh=zeros(1,roof);
box=cell(1,roof);
gui.rows=cell(1,roof);
gui.cols=cell(1,roof);

for tagno=taggroup
 
  [rows,cols] = calcfunctions('calcrowscols',tagno,SET(tagno).ZSize);
  gui.rows{tagno}=rows;
  gui.cols{tagno}=cols;
  xsz(tagno) = SET(tagno).XSize;
  ysz(tagno) = SET(tagno).YSize;
  xres(tagno) = SET(tagno).ResolutionX;
  yres(tagno) = SET(tagno).ResolutionY;
  zsz(tagno)=SET(tagno).ZSize;

  [im,slicestoplot] = calcfunctions('calcmontageviewim',tagno,[rows,cols],0,[],[],[],[],1,false);
  
  switch SET(tagno).ImageViewPlane
    case '2CH'
      twochamber = im(:,:,1);
    case '3CH'
      threechamber = im(:,:,1);
    case '4CH'
      fourchamber = im(:,:,1);
  end
end

%preallocate imagehandles
gui.handles.im2CH=[];
gui.handles.im3CH=[];
gui.handles.im4CH=[];

for tagno=taggroup
  switch SET(tagno).ImageViewPlane
    case '2CH'
      axh(tagno)=gui.handles.LAX2CH_axes;
      colormap(axh(tagno),gray(255));
      gui.handles.im2CH=image(cat(3,twochamber,twochamber,twochamber),'Parent',axh(tagno));
      axis(axh(tagno),'off','image')
      hold(axh(tagno), 'on')
      set(axh(tagno),'plotboxaspectratio',[max(zsz)*...
        ysz(tagno)*yres(tagno) xsz(tagno)*xres(tagno) 1]);
      box{tagno}=[1,1;SET(tagno).XSize,1;SET(tagno).XSize,SET(tagno).YSize;...
        1,SET(tagno).YSize;1,1];
      gui.handles.frame2CH = plot(axh(tagno),box{tagno}(:,2),box{tagno}(:,1),'y','LineWidth',3);
    case '3CH'
      axh(tagno)=gui.handles.LAX3CH_axes;
      colormap(axh(tagno),gray(255));
      gui.handles.im3CH=image(cat(3,threechamber,threechamber,threechamber),'Parent',axh(tagno));
      axis(axh(tagno),'off','image')
      hold(axh(tagno), 'on')
      set(axh(tagno),'plotboxaspectratio',[max(zsz)*...
        ysz(tagno)*yres(tagno) xsz(tagno)*xres(tagno) 1]);
      box{tagno}=[1,1;SET(tagno).XSize,1;SET(tagno).XSize,SET(tagno).YSize;...
        1,SET(tagno).YSize;1,1];
      gui.handles.frame3CH = plot(axh(tagno),box{tagno}(:,2),box{tagno}(:,1),'y','LineWidth',3);
    case '4CH'
      axh(tagno)=gui.handles.LAX4CH_axes;
      colormap(axh(tagno),gray(255));
      gui.handles.im4CH=image(cat(3,fourchamber,fourchamber,fourchamber),'Parent',axh(tagno));
      axis(axh(tagno),'off','image')
      hold(axh(tagno), 'on')
      set(axh(tagno),'plotboxaspectratio',[max(zsz)*...
        ysz(tagno)*yres(tagno) xsz(tagno)*xres(tagno) 1]);
      box{tagno}=[1,1;SET(tagno).XSize,1;SET(tagno).XSize,SET(tagno).YSize;...
        1,SET(tagno).YSize;1,1];
      gui.handles.frame4CH = plot(axh(tagno),box{tagno}(:,2),box{tagno}(:,1),'y','LineWidth',3);
  end
end

%initiates the parameter slices to use
gui.slicestouse=cell(1,roof);
for tagno=taggroup
  gui.slicestouse{tagno}=1;
end

if isempty(gui.handles.im2CH)
  gui.handles.im2CH=imshow(nan,'Parent',gui.handles.LAX2CH_axes);
else
  set(gui.handles.im2CH,'ButtonDownFcn', ...
    'straintagging.LAXselectslices(''imlight'',''2CH'')');  
end

if isempty(gui.handles.im3CH)
  gui.handles.im3CH=imshow(nan,'Parent',gui.handles.LAX3CH_axes);
else
  set(gui.handles.im3CH,'ButtonDownFcn', ...
    'straintagging.LAXselectslices(''imlight'',''3CH'')');  
end

if isempty(gui.handles.im4CH)
    gui.handles.im4CH=imshow(nan,'Parent',gui.handles.LAX4CH_axes);
else
  set(gui.handles.im4CH,'ButtonDownFcn', ...
    'straintagging.LAXselectslices(''imlight'',''4CH'')');  
end


%-------------------
function ok_callback %#ok<DEFNU>
%-------------------
%generates new stack then terminates the gui and start strain analysis

global DATA SET

gui=DATA.GUI.LAXselectslices;
taggroup=gui.taggroup;
strainno=LAXgenerateimagestack(taggroup);  %generates new image stacks

%correct the taggroup, should be in order 2CH, 3CH, 4CH as they exist

%find order for 2CH, 3CH, 4CH
chamberexist(1,:)=strcmp('2CH',{SET.ImageViewPlane});
chamberexist(2,:)=strcmp('3CH',{SET.ImageViewPlane});
chamberexist(3,:)=strcmp('4CH',{SET.ImageViewPlane});
taggroupnew = [];
for chloop = 1:3
  taggroupnew = [taggroupnew find(chamberexist(chloop,:))];
end
straingroup = strainno;
straingroup(isnan(strainno)) = taggroup(isnan(strainno));
for noloop = straingroup
  SET(noloop).StrainTagging.taggroup = taggroupnew;
end

close_callback;
straintagging.straintagging('init','cine','longaxis',straingroup(1));

%----------------------
function close_callback
%----------------------
%closes the gui without stacksplitting
global DATA

try
  DATA.GUI.LAXselectslices = close(DATA.GUI.LAXselectslices);
catch   %#ok<CTCH>
  DATA.GUI.LAXselectslices = [];
  delete(gcbf);
end

%----------------------------------
function strainnoall=LAXgenerateimagestack(taggroup)
%----------------------------------
global DATA SET
gui=DATA.GUI.LAXselectslices;
strainnoall = [];
ind = 1;

for tagno=taggroup
  if SET(tagno).ZSize>1
    %Create new SET structure
    laxset = SET(tagno);
    laxset.StartSlice = 1;
    laxset.EndSlice = 1;
    laxset.CurrentSlice = 1;
    laxset.ImageType = SET(tagno).ImageType;
    laxset.Linked = numel(SET)+1;
    %update image information
    laxset.IM = SET(tagno).IM(:,:,:,gui.slicestouse{tagno});
    laxset.ZSize = 1;
    laxset.SliceGap =SET(tagno).SliceGap;
    %update LV and RV segmentation parameters
    fields = {'Endo','Epi','RVEndo','RVEpi'};
    for k = 1:length(fields)
      segfield = fields{k};
      laxset.([segfield 'X']) = [];
      laxset.([segfield 'Y']) = [];
      if ~isempty(SET(tagno).([segfield 'X']))
        laxset.([segfield 'X'])(:,:,1) = SET(tagno).([segfield 'X'])(:,:,gui.slicestouse{tagno});
        laxset.([segfield 'Y'])(:,:,1) = SET(tagno).([segfield 'Y'])(:,:,gui.slicestouse{tagno});
      end
      laxset.([segfield 'PinX']) = [];
      laxset.([segfield 'PinY']) = [];
      laxset.([segfield 'InterpX']) = [];
      laxset.([segfield 'InterpY']) = [];
      laxset.([segfield 'PinXView']) = [];
      laxset.([segfield 'PinYView']) = [];
      if ~isempty(laxset.([segfield 'X']))
        laxset.([segfield 'XView']) = nan((DATA.NumPoints+1)*laxset.ZSize,laxset.TSize);
        laxset.([segfield 'YView']) = laxset.EndoXView;
      else
        laxset.([segfield 'XView']) = NaN;
        laxset.([segfield 'YView']) = NaN;
      end
    end
    %reset Points, Measure, Mar, Scar, Roi
    laxset.MaR = [];
    laxset.Scar = [];
    laxset.Measure = [];
    laxset.LevelSet = [];
    laxset.AtrialScar = [];
    %roi
    laxset.Roi = roi('roireset');
    laxset.RoiN = 0;
    laxset.RoiCurrent = [];
    %point
    laxset.Point.X = [];
    laxset.Point.Y = [];
    laxset.Point.T = [];
    laxset.Point.Z = [];
    laxset.Point.Label = {};
    %soom state
    laxset.NormalZoomState = [];
    laxset.MontageZoomState = [];
    laxset.MontageRowZoomState = [];
    laxset.MontageFitZoomState = [];
    
    %add new image stack
    strainno = numel(SET)+1;
    SET(strainno) = laxset;
    SET(strainno).KeptSlices = gui.slicestouse{tagno};
    strainnoall(ind) = strainno;
    
    %Change some descriptors in previous making obscured for strain
    %analysis
    %SET(tagno).ImageType='Scout';
    SET(tagno).ImageViewPlane='Unspecified';
    
    %draw new image stack and update volumes
    viewfunctions('setview',1,1,strainno,{'montage'})    
    drawfunctions('drawthumbnails');
    segment_main('updatevolume');
  else
    strainnoall(ind) = NaN;
  end
  ind = ind+1;
end

%-----------------------------------------
function imlight(view)
%-----------------------------------
global DATA SET

gui=DATA.GUI.LAXselectslices;
taggroup=gui.taggroup;

for tagno=taggroup
  if strcmp(SET(tagno).ImageViewPlane,view)
    view_ind=tagno;
    break;
  end
end

switch view
  case '2CH'
    [y,x]=mygetcurrentpoint(gui.handles.LAX2CH_axes);
  case '3CH'
    [y,x]=mygetcurrentpoint(gui.handles.LAX3CH_axes);
  case '4CH'
    [y,x]=mygetcurrentpoint(gui.handles.LAX4CH_axes);
end

X=SET(view_ind).XSize*(0:SET(view_ind).ZSize);
Y=SET(view_ind).YSize*(0:SET(view_ind).ZSize);

for i=1:length(X)-1
  if (x>X(i) && x<X(i+1))
    whichboxx=i;
  end
  if (y>Y(i) && y<Y(i+1))
    whichboxy=i;
  end
end

sliceind=(whichboxx-1)*gui.rows{view_ind}+whichboxy;
if sliceind<=SET(view_ind).ZSize
  floorbox=[1,1;SET(view_ind).XSize,1;SET(view_ind).XSize,SET(view_ind).YSize;...
    1,SET(view_ind).YSize;1,1];
  box=[floorbox(:,1)+X(whichboxx),floorbox(:,2)+Y(whichboxy)];
  switch view
    case '2CH'
      delete(gui.handles.frame2CH);
      gui.handles.frame2CH=plot(gui.handles.LAX2CH_axes,box(:,2),box(:,1),'y','LineWidth',3);
    case '3CH'
      delete(gui.handles.frame3CH);
      gui.handles.frame3CH=plot(gui.handles.LAX3CH_axes,box(:,2),box(:,1),'y','LineWidth',3);
    case '4CH'
      delete(gui.handles.frame4CH);
      gui.handles.frame4CH=plot(gui.handles.LAX4CH_axes,box(:,2),box(:,1),'y','LineWidth',3);
  end
  
  gui.slicestouse{view_ind}=sliceind;
end