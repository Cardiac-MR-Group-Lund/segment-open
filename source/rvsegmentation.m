function varargout = rvsegmentation(varargin)
%function for initialization of automatic RV segmentation

macro_helper(varargin{:});
if nargin < 1 || isempty(varargin{1})
  varargin{1} = 'init';
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init
%------------
%Initialize the GUI

global DATA SET

[saxno,~,~,ch4no,success]=lvsegmentation('croplvall',1);

if not(success) || isempty(saxno)
  myfailed('Could not find short-axis image stack for RV analysis.',DATA.GUI.Segment);
  return
end

if SET(saxno).ZSize < 5
  myfailed('Need at least 5 slices for fully automated segmentation.',DATA.GUI.Segment);
  return
end

%Open RV wizard GUI
gui = mygui('rvsegmentation.fig');
DATA.GUI.RVSegmentation = gui;

gui.saxno = saxno;
gui.ch4no = ch4no;

%set ED as current time frame
tools('enddiastole_Callback');

%Init time slider
gui.tf = SET(gui.saxno).CurrentTimeFrame;
tsz = SET(gui.saxno).TSize;
if tsz == 1
  set(gui.handles.timebaraxes,'visible','off');
else
  sliderstep = [1/(tsz-1) 5/(tsz-1)];
  set(gui.handles.timeslider,'Min',1,'Max',tsz,'Value',gui.tf, ...
    'SliderStep',sliderstep);
  set(gui.handles.timetext,'String',sprintf('%0.0f ms',1000*SET(gui.saxno).TimeVector(gui.tf)));
  
  load Icons.mat
  set(gui.handles.prevpushbutton,'CData',icon.prev);
  set(gui.handles.playtogglebutton,'CData',icon.play);
  set(gui.handles.nextpushbutton,'CData',icon.next);
  inittimebar;
end

%find lv center
if isempty(SET(gui.saxno).StartSlice) || isempty(SET(gui.saxno).EndSlice)
  SET(gui.saxno).StartSlice = round(SET(gui.saxno).ZSize);
  SET(gui.saxno).EndSlice = round(SET(gui.saxno).ZSize);
end
if SET(gui.saxno).EndSlice - SET(gui.saxno).StartSlice < 4
  SET(gui.saxno).EndSlice = min([SET(gui.saxno).ZSize ceil(SET(gui.saxno).ZSize*3/4)]);
  SET(gui.saxno).StartSlice = max([1 floor(SET(gui.saxno).ZSize/4)]);
end
gui.midslice = floor((SET(gui.saxno).StartSlice + SET(gui.saxno).EndSlice)/2);
[x,y] = findfunctions('findrvcenter',gui.saxno,SET(gui.saxno).StartSlice:SET(gui.saxno).EndSlice);
SET(gui.saxno).CenterX=x;
SET(gui.saxno).CenterY=y;

[rows,cols] = calcfunctions('calcrowscols',gui.saxno);
gui.matrix = [rows cols]; %viewmatrix for shortaxis image

gui.saxim = calcfunctions('calcmontageviewim',gui.saxno,gui.matrix);
gui.midim = calcfunctions('remapuint8',SET(gui.saxno).IM(:,:,:,gui.midslice),gui.saxno);
emptyim = uint8(zeros(1,1,SET(gui.saxno).TSize));
gui.ch4im = emptyim;
if ~isempty(gui.ch4no)
  tfs = round(linspace(1,SET(gui.ch4no).TSize,SET(gui.saxno).TSize));
  gui.ch4im = calcfunctions('remapuint8',SET(gui.ch4no).IM(:,:,tfs),gui.ch4no);
end

%plot sax image stack
gui.handles.saximh = image(gui.saxim(:,:,gui.tf),'Parent',gui.handles.shortaxisaxes);
% gui.handles.midimh = image(gui.midim(:,:,gui.tf),'Parent',gui.handles.midsliceaxes);

axh = [gui.handles.shortaxisaxes gui.handles.ch4axes]; %gui.handles.midsliceaxes 
colormap(axh(1),gray(255));
colormap(axh(2),gray(255));
hold(axh(1),'on');
hold(axh(2),'on');
set(gui.handles.saximh,'ButtonDownFcn', ...
  'rvsegmentation(''shortaxisaxes_Buttondown'')');
% set(gui.handles.midimh,'ButtonDownFcn', ...
%   'rvsegmentation(''midsliceaxes_Buttondown'')');

%plot long-axis image stacks
if ~isempty(gui.ch4no)
gui.handles.ch4imh = image(gui.ch4im(:,:,gui.tf),'Parent',gui.handles.ch4axes);
colormap(axh(2),gray(255));
hold(axh(2),'on');
set(gui.handles.ch4imh,'ButtonDownFcn', ...
  'lvsegmentation(''ch4axes_Buttondown'')');
else
  set(gui.handles.ch4axes,'visible','off');
end

axis(axh,'off');
axis(axh,'equal');

%Initiate plot handles
gui.handles.startsliceframe = plot(gui.handles.shortaxisaxes,[],[]);
gui.handles.endsliceframe = plot(gui.handles.shortaxisaxes,[],[]);
gui.handles.startslice4ch = plot(gui.handles.ch4axes,[],[]);
gui.handles.endslice4ch = plot(gui.handles.ch4axes,[],[]);
gui.handles.centerbasal = plot(gui.handles.shortaxisaxes,[],[]);
gui.handles.centerapical = plot(gui.handles.shortaxisaxes,[],[]);
gui.handles.center4ch = plot(gui.handles.ch4axes,[],[]);

set(gui.fig,'KeyPressFcn', ...
  @(hObject,eventdata)rvsegmentation('keypressed',hObject,eventdata));

%use user defined sliceselction and centerpoint if it exist
if isfield(SET(gui.saxno),'RV') && ~isempty(SET(gui.saxno).RV) && isfield(SET(gui.saxno).RV,'centerbasal') && ~isempty(SET(gui.saxno).RV.centerbasal) && isfield(SET(gui.saxno).RV,'slicebasal') && ~isempty(SET(gui.saxno).RV.slicebasal)
  SET(gui.saxno).StartSlice = SET(gui.saxno).RV.slicebasal;
  SET(gui.saxno).EndSlice = SET(gui.saxno).RV.sliceapical;
  updateselectedslices;
  %calcualte center point for montage view
  colbasal = mod(SET(gui.saxno).StartSlice,gui.matrix(2));
  if colbasal == 0, colbasal=gui.matrix(2); end
  rowbasal = ceil(SET(gui.saxno).StartSlice/gui.matrix(2));
  colapical = mod(SET(gui.saxno).EndSlice,gui.matrix(2));
  if colapical == 0, colapical=gui.matrix(2); end
  rowapical = ceil(SET(gui.saxno).EndSlice/gui.matrix(2));
  xposbasal = round((rowbasal-1)*SET(gui.saxno).XSize)+SET(gui.saxno).RV.centerbasal(1);
  yposbasal = round((colbasal-1)*SET(gui.saxno).YSize)+SET(gui.saxno).RV.centerbasal(2);
  xposapical = round((rowapical-1)*SET(gui.saxno).XSize)+SET(gui.saxno).RV.centerapical(1);
  yposapical= round((colapical-1)*SET(gui.saxno).YSize)+SET(gui.saxno).RV.centerapical(2);
  %plot center points
  updatecenterpoint(xposbasal,yposbasal,SET(gui.saxno).RV.centerbasal(1),SET(gui.saxno).RV.centerbasal(2),'basal');
  updatecenterpoint(xposapical,yposapical,SET(gui.saxno).RV.centerapical(1),SET(gui.saxno).RV.centerapical(2),'apical');  
else
  updateselectedslices;
  %default center point
  thisslicex = round(0.5*SET(gui.saxno).XSize);
  thisslicey = round(0.25*SET(gui.saxno).YSize);
  %calcualte center point for montage view
  colbasal = mod(SET(gui.saxno).StartSlice,gui.matrix(2));
  if colbasal == 0, colbasal=gui.matrix(2); end
  rowbasal = ceil(SET(gui.saxno).StartSlice/gui.matrix(2));
  colapical = mod(SET(gui.saxno).EndSlice,gui.matrix(2));
  if colapical == 0, colapical=gui.matrix(2); end
  rowapical = ceil(SET(gui.saxno).EndSlice/gui.matrix(2));
  xposbasal = round((rowbasal-1)*SET(gui.saxno).XSize)+thisslicex;
  yposbasal = round((colbasal-1)*SET(gui.saxno).YSize)+thisslicey;
  xposapical = round((rowapical-1)*SET(gui.saxno).XSize)+thisslicex;
  yposapical= round((colapical-1)*SET(gui.saxno).YSize)+thisslicey;
  
  SET(gui.saxno).RV.centerbasal = [thisslicex,thisslicey];
  SET(gui.saxno).RV.centerapical = [thisslicex,thisslicey];
  SET(gui.saxno).RV.slicebasal = SET(gui.saxno).StartSlice;
  SET(gui.saxno).RV.sliceapical = SET(gui.saxno).EndSlice;
  
  updatecenterpoint(xposbasal,yposbasal,thisslicex,thisslicey,'basal');
  updatecenterpoint(xposapical,yposapical,thisslicex,thisslicey,'apical');
end
set(gui.handles.centerbasal,'ButtonDownFcn', ...
  'rvsegmentation(''shortaxisaxes_Buttondown'')');
set(gui.handles.centerapical,'ButtonDownFcn', ...
  'rvsegmentation(''shortaxisaxes_Buttondown'')');

%-------------------
function inittimebar
%-------------------
%Initiate timebar axis
global DATA SET
gui = DATA.GUI.RVSegmentation;

h = gui.handles.timebaraxes;
no = gui.saxno;
if isempty(no)
  return
end
set(h,'ButtonDownFcn', ...
  @(hObject,eventdata)rvsegmentation('timebaraxes_ButtonDownFcn',hObject,eventdata));

delete(get(h,'Children'));

tvec = SET(no).TimeVector;

hold(h,'on');
fcn = @(hObject,eventdata)rvsegmentation('timebar_ButtonDownFcn',hObject,eventdata);

%Draw timebar (red) and set its buttondown fcn
gui.handles.timebar = plot(h,tvec(gui.tf)*[1 1],[0 1],'b','Tag','currenttime');
set(gui.handles.timebar,'ButtonDownFcn',fcn);

%Set axes options
tvec = SET(no).TimeVector;
tmin = tvec(1);
tmax = tvec(end);

marg = (tmax-tmin)/500;
axis(h,[tmin-marg tmax+marg 0 1]);
tickvec = tmin:0.2:tmax;
set(h, 'XTick', tickvec,'YTick',[],'XMinorTick','off','XGrid','on', ...
  'XTickLabel',cellfun(@num2str,num2cell(1000*tickvec), ...
  'UniformOutput',false),'TickLength',[0.005 0.025]);

temp = get(gui.handles.timebaraxes,'ylim');
ttemp = [tvec;tvec];
temp = [repmat(temp(1),size(tvec));repmat(temp(1)+0.1*(temp(2)-temp(1)),size(tvec))];
plot(gui.handles.timebaraxes,ttemp,temp,'k-');

updatetimebar;

%---------------------
function updatetimebar
%---------------------
%Update timebar axis
global DATA SET
gui = DATA.GUI.RVSegmentation;

h = gui.handles.timebaraxes;
no = gui.saxno;
if isempty(no)
  return
end

fcn = get(h,'ButtonDownFcn');
currenttime = SET(gui.saxno).TimeVector(gui.tf);

%Update timebar
set(gui.handles.timebar,'XData',currenttime*[1 1])

set(h,'ButtonDownFcn',fcn);


%---------------------------------------------
function timebaraxes_ButtonDownFcn(hObject, ~)
%---------------------------------------------
%Buttondown function for timebar axes. Changes current timeframe to the 
%one closest to position of clicked point
global DATA SET
gui = DATA.GUI.RVSegmentation;
no = gui.saxno;
x = mygetcurrentpoint(hObject);
[~,tf] = min(abs(SET(no).TimeVector-x));
gui.tf = max(min(tf,SET(no).TSize),1);
updateimages;

%-----------------------------------------
function timebar_ButtonDownFcn(hObject, ~)
%-----------------------------------------
%Buttondown function for graphical timebar object. Activates dragging of timebars.
global DATA
gui = DATA.GUI.RVSegmentation;
no = gui.saxno;
obj = hObject;
motionfcn = @(hObject,eventdata)rvsegmentation('timebaraxes_MotionFcn',hObject,eventdata,obj,no);
set(gui.fig,'WindowButtonMotionFcn',motionfcn);
buttonupfcn = @(hObject,eventdata)rvsegmentation('timebaraxes_ButtonUpFcn',hObject,eventdata);
set(gui.fig,'WindowButtonUpFcn', buttonupfcn);

%----------------------------------------------
function timebaraxes_MotionFcn(~, ~, tbobj, no)
%----------------------------------------------
%Mouse motion function for timebar axes of image specified by input 
%parameter 'field'. Used for dragging timebars to change current
%timeframe or start/end points of timeframes in which to align images.
global DATA SET
x = mygetcurrentpoint(get(tbobj,'Parent'));
[~,tf] = min(abs(SET(no).TimeVector-x));
DATA.GUI.RVSegmentation.tf = max(min(tf,SET(no).TSize),1);
updateimages;

%-------------------------------------------
function timebaraxes_ButtonUpFcn(hObject, ~)
%-------------------------------------------
%Buttonup function for timebar axes of image specified by input 
%parameter 'field'. Deactivates dragging of timebar.
set(hObject,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[]);


%---------------------------
function timeslider_Callback
%---------------------------
%Callback for timeslider
global DATA
gui = DATA.GUI.RVSegmentation;

%Update timeframe and images
gui.tf = round(mygetvalue(gui.handles.timeslider));
set(gui.handles.timeslider,'Value',gui.tf);
updateimages;

%--------------------
function updateimages
%--------------------
%Update all image stacks and also timebar
global DATA SET
gui = DATA.GUI.RVSegmentation;
set(gui.handles.saximh,'CData',gui.saxim(:,:,gui.tf));
% set(gui.handles.midimh,'CData',gui.midim(:,:,gui.tf));
if isfield(gui.handles,'ch4imh')
  set(gui.handles.ch4imh,'CData',gui.ch4im(:,:,gui.tf));
end
set(gui.handles.timetext,'String', ...
  sprintf('%0.0f ms',1000*SET(gui.saxno).TimeVector(gui.tf)));
if SET(gui.saxno).TSize > 1
  updatetimebar;
end

%---------------------
function prev_Callback
%---------------------
%Callback for previous timeframe pushbutton
global DATA SET
gui = DATA.GUI.RVSegmentation;
gui.tf = mod(gui.tf-2,SET(gui.saxno).TSize)+1;
updateimages;

%---------------------
function next_Callback
%---------------------
%Callback for next timeframe pushbutton
global DATA SET
gui = DATA.GUI.RVSegmentation;
gui.tf = mod(gui.tf,SET(gui.saxno).TSize)+1;
updateimages;

%---------------------
function play_Callback
%---------------------
%Callback for play togglebutton
global DATA SET
gui = DATA.GUI.RVSegmentation;

h = gui.handles.playtogglebutton;
no = gui.saxno;
while ishandle(h) && mygetvalue(h)
  gui.tf = mod(gui.tf,SET(no).TSize)+1;
  set(gui.handles.timeslider,'Value',gui.tf);
  updateimages;
  %Do pause to release CPU burden.
  pause(0.5*SET(no).BeatTime/SET(no).TSize);
end

%----------------------------
function updateselectedslices
%----------------------------
%Draw markers around selected slices
global DATA SET
gui = DATA.GUI.RVSegmentation;

%Clean up before drawing new frames and lines
delete([gui.handles.startsliceframe gui.handles.endsliceframe ...
  gui.handles.startslice4ch gui.handles.endslice4ch]);

if (SET(gui.saxno).XSize*gui.matrix(1)) < (SET(gui.saxno).YSize*gui.matrix(2))
  scalefactorx = 1;%(SET(gui.saxno).XSize*gui.matrix(1))/(SET(gui.saxno).YSize*gui.matrix(2));
  scalefactory = 1;
else
  scalefactorx = 1;
  scalefactory = 1;%(SET(gui.saxno).YSize*gui.matrix(2))/(SET(gui.saxno).XSize*gui.matrix(1));
end
boxx = scalefactorx*[-1 -1 1 1 -1]*SET(gui.saxno).XSize/2+0.5*[1 1 -1 -1 1];%gui.matrix(1);
boxy = scalefactory*[-1 1 1 -1 -1]*SET(gui.saxno).YSize/2+0.5*[1 -1 -1 1 1];%gui.matrix(2);
for slice = unique([SET(gui.saxno).StartSlice SET(gui.saxno).EndSlice])
  row = ceil(slice/gui.matrix(2));
  col = slice-(row-1)*gui.matrix(2);
  centerx = (row-0.5)*SET(gui.saxno).XSize+0.5;%;%+1;%
  centery = (col-0.5)*SET(gui.saxno).YSize+0.5;%;%+1;%
  if slice == SET(gui.saxno).StartSlice
    gui.handles.startsliceframe = plot( ...
      gui.handles.shortaxisaxes,centery+boxy,centerx+boxx,'b','LineWidth',2);
  end
  if slice == SET(gui.saxno).EndSlice
    gui.handles.endsliceframe = plot( ...
      gui.handles.shortaxisaxes,centery+boxy,centerx+boxx,'r','LineWidth',2);
  end
end

if ~isempty(gui.ch4no)
  [ssx,ssy] = calcfunctions('calcplaneintersections', ...
    gui.ch4no,gui.saxno,'one','one',SET(gui.saxno).StartSlice);
  [esx,esy] = calcfunctions('calcplaneintersections', ...
    gui.ch4no,gui.saxno,'one','one',SET(gui.saxno).EndSlice);
  gui.handles.startslice4ch = plot( ...
    gui.handles.ch4axes,ssy,ssx,'b');
  gui.handles.endslice4ch = plot( ...
    gui.handles.ch4axes,esy,esx,'r');
end

segment('updateselectedslices');


%----------------------------------------------------------------
function updatecenterpoint(x,y,thisslicex,thisslicey,centerslice)
%----------------------------------------------------------------
%Update center point
global DATA SET
gui = DATA.GUI.RVSegmentation;

if nargin < 1
  x = SET(gui.saxno).CenterX;
  y = SET(gui.saxno).CenterY;
  thisslicex = SET(gui.saxno).CenterX;
  thisslicey = SET(gui.saxno).CenterY;
  centerslice = 'basal';
end

crosscolor = [1 0.6 0.2]; %DATA.centercrossdef
[tempx,tempy] = ndgrid(0:(gui.matrix(1)-1),0:(gui.matrix(2)-1));
switch centerslice
  case 'basal'
    %Clean up before drawing
    delete(gui.handles.centerbasal);
    gui.handles.centerbasal = plot(gui.handles.shortaxisaxes,y,x, ...
      '+','color',crosscolor,'LineWidth',3,'MarkerSize',12);
    startxyz = [thisslicex;thisslicey;SET(gui.saxno).StartSlice];
    endxyz = [SET(gui.saxno).RV.centerapical(1); SET(gui.saxno).RV.centerapical(2); SET(gui.saxno).EndSlice];
    set(gui.handles.centerbasal,'ButtonDownFcn', ...
      'rvsegmentation(''shortaxisaxes_Buttondown'')');
  case 'apical'
    %Clean up before drawing
    delete(gui.handles.centerapical);
    gui.handles.centerapical = plot(gui.handles.shortaxisaxes,y,x, ...
      '+','color',crosscolor,'LineWidth',3,'MarkerSize',12);
    startxyz = [SET(gui.saxno).RV.centerbasal(1); SET(gui.saxno).RV.centerbasal(2); SET(gui.saxno).StartSlice];
    endxyz = [thisslicex;thisslicey; SET(gui.saxno).EndSlice];
    set(gui.handles.centerapical,'ButtonDownFcn', ...
      'rvsegmentation(''shortaxisaxes_Buttondown'')');
end

%Get center point coordinates in longaxis images
[saxpos,saxmat] = calcfunctions('calcormat',gui.saxno);
truecenter = [saxpos saxpos] + saxmat*([startxyz endxyz]-1);

if ~isempty(gui.ch4no)
  %Clean up before drawing
  delete(gui.handles.center4ch);
  [ch4pos,ch4mat] = calcfunctions('calcormat',gui.ch4no);
  ch4center = 1 + ch4mat\(truecenter - [ch4pos ch4pos]);
  gui.handles.center4ch = plot( ...
    gui.handles.ch4axes,ch4center(2,:),ch4center(1,:),'Color',crosscolor);
end

%--------------------------------
function shortaxisaxes_Buttondown %#ok<DEFNU>
%--------------------------------
global DATA SET
gui = DATA.GUI.RVSegmentation;

%Extract coordinates clicked
[y,x] = mygetcurrentpoint(gui.handles.shortaxisaxes);
%Find slice
col = 1+floor((y-0.5)/SET(gui.saxno).YSize);
row = 1+floor((x-0.5)/SET(gui.saxno).XSize);
slice = col+(row-1)*gui.matrix(2);

if slice > SET(gui.saxno).ZSize
  return;
end
slicediff = [abs(SET(gui.saxno).StartSlice-slice) abs(SET(gui.saxno).EndSlice-slice)];
[~,closestslice] = min(slicediff);
switch closestslice
  case 1
    SET(gui.saxno).StartSlice = slice;
    SET(gui.saxno).RV.slicebasal = slice;
    %storing points for usage in RV segmentation
    thisslicex = max(1,min(SET(gui.saxno).XSize,round(x-(row-1)*SET(gui.saxno).XSize)));
    SET(gui.saxno).RV.centerbasal(1) = thisslicex;
    thisslicey = max(1,min(SET(gui.saxno).YSize,round(y-(col-1)*SET(gui.saxno).YSize)));
    SET(gui.saxno).RV.centerbasal(2) = thisslicey;
    %plot centerpoint
    updatecenterpoint(x,y,thisslicex,thisslicey,'basal');
  case 2
    SET(gui.saxno).EndSlice = slice;
    SET(gui.saxno).RV.sliceapical= slice;
    %storing points for usage in RV segmentation
    thisslicex = max(1,min(SET(gui.saxno).XSize,round(x-(row-1)*SET(gui.saxno).XSize)));
    SET(gui.saxno).RV.centerapical(1) = thisslicex;
    thisslicey = max(1,min(SET(gui.saxno).YSize,round(y-(col-1)*SET(gui.saxno).YSize)));
    SET(gui.saxno).RV.centerapical(2) = thisslicey;
    %plot centerpoint
    updatecenterpoint(x,y,thisslicex,thisslicey,'apical');
end
updateselectedslices;


%-------------------------------------
function keypressed(hObject,eventdata)
%-------------------------------------
%Keypress function for GUI
global DATA SET
gui = DATA.GUI.RVSegmentation;

switch eventdata.Key
%   case 'uparrow'
%     if isempty(eventdata.Modifier)
%       SET(gui.saxno).StartSlice = max(1,SET(gui.saxno).StartSlice-1);
%     else
%       SET(gui.saxno).EndSlice = max(SET(gui.saxno).StartSlice,SET(gui.saxno).EndSlice-1);
%     end
%     updateselectedslices;
%   case 'downarrow'
%     if isempty(eventdata.Modifier)
%       SET(gui.saxno).StartSlice = min(SET(gui.saxno).EndSlice,SET(gui.saxno).StartSlice+1);
%     else
%       SET(gui.saxno).EndSlice = min(SET(gui.saxno).ZSize,SET(gui.saxno).EndSlice+1);
%     end
%     updateselectedslices;
  case 'leftarrow'
    prev_Callback;
  case 'rightarrow'
    next_Callback;
  case 'p'
    set(gui.handles.playtogglebutton,'value',1-mygetvalue(gui.handles.playtogglebutton));
    play_Callback;
  case 'd'    
    if get(gui.handles.playtogglebutton,'value')
      set(gui.handles.playtogglebutton,'value',0)
      play_Callback;
    end
    gui.tf = SET(gui.saxno).EDT;
    updateimages;
  case 's'
    if get(gui.handles.playtogglebutton,'value')
      set(gui.handles.playtogglebutton,'value',0)
      play_Callback;
    end
    gui.tf = SET(gui.saxno).EST;
    updateimages;
end

%---------------------
function dorv_Callback %#ok<DEFNU>
%---------------------
%Do RV segmentation in selected slices
global DATA
gui = DATA.GUI.RVSegmentation;
segment('switchtoimagestack',gui.saxno);
oldpol=DATA.ThisFrameOnly;
DATA.ThisFrameOnly=0;
rv('segmentrvendo_Callback',0);
close_Callback
drawfunctions('drawall',DATA.ViewMatrix(1),DATA.ViewMatrix(2));
figure(DATA.GUI.Segment.fig);
DATA.ThisFrameOnly=oldpol;

%----------------------
function close_Callback  
%----------------------
%Close RV segmentation GUI
global DATA

try
  DATA.GUI.RVSegmentation = close(DATA.GUI.RVSegmentation);
catch   %#ok<CTCH>
  DATA.GUI.RVSegmentation=[];
  delete(gcbf);
end
