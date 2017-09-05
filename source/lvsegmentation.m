function varargout = lvsegmentation(varargin)
% LVSEGMENTATION MATLAB code for lvsegmentation.fig

macro_helper(varargin{:});
if nargin < 1 || isempty(varargin{1})
  varargin{1} = 'init';
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init
%------------
%Initialize the GUI

global DATA SET NO

cineno = findfunctions('findno');
saxno = find(strcmp({SET.ImageViewPlane},'Short-axis'));
zno = find([SET.ZSize] > 1);
saxno = intersect(cineno,union(saxno,zno));

if isempty(saxno)
  if ismember(NO,zno) && ...
      yesno('Could not find image stack defined as Short-axis Cine. Use current stack?');
    saxno = NO;
    SET(saxno).ImageViewPlane = 'Short-axis';
    SET(saxno).ImageType = 'Cine';
    drawfunctions('updatenopanels',saxno);
  else
    myfailed('Could not find short-axis image stack');
    return
  end
end

if ismember(NO,saxno)
  saxno = NO;
else
  nos = cellfun(@(x)dprintf('Short-axis stack %d',x), ...
    num2cell(saxno),'UniformOutput',false);
  m = mymenu('Do segmentation on short-axis stack?', ...
    nos{:});
  if m
    saxno = saxno(m);
  else
    return;
  end
end

if SET(saxno).ZSize < 5
  myfailed('Need at least 5 slices for fully automated segmentation (endo&epi).',DATA.GUI.Segment);
  return
end

nostocrop = [];
if SET(saxno).XSize * SET(saxno).ResolutionX > 210 || ...
    SET(saxno).YSize * SET(saxno).ResolutionY > 210
  nostocrop = saxno;
end

ch2no = [];
ch3no = [];
ch4no = [];
for cham = 2:4
  chno = find(strcmp({SET(:).ImageViewPlane},sprintf('%dCH',cham)));
  if ~isempty(chno)
    if length(chno) > 1
      cinechno = find(strcmp({SET(cineno).ImageViewPlane},sprintf('%dCH',cham)));
      if ~isempty(cinechno) && length(cinechno) == 1
        chno = cineno(cinechno);
      else
        %ask user
        nostri = '(';
        for loop = 1:length(chno)
          if loop == length(chno)
            nostri = [nostri sprintf('%d)',chno(loop))];
          else
            nostri = [nostri sprintf('%d, ',chno(loop))];
          end
        end
        uniquechno = inputdlg({dprintf('Select %dCH image stack out of %s.',cham,nostri)},'Image Stack',1,{sprintf('%d',chno(1))});
        if ~ismember(str2num(uniquechno{1}),chno)
          uniquechno = [];
        end
        if isempty(uniquechno)
          myfailed('Invalid image stack.',DATA.GUI.Segment);
          return;
        else
          [chno,ok] = str2num(uniquechno{1}); %#ok<ST2NM>
          if not(ok)
            myfailed('Invalid image stack.',DATA.GUI.Segment);
            return;
          end;
        end;
      end
    end
    if cham == 2
      ch2no = chno;
    elseif cham == 3
      ch3no = chno;
    elseif cham == 4
      ch4no = chno;
    end
    %Check if image needs cropping
    if SET(chno).XSize * SET(chno).ResolutionX > 260 || ...
    SET(chno).YSize * SET(chno).ResolutionY > 260
      nostocrop = [nostocrop chno]; %#ok<AGROW>
    end
  end
end

if ~isempty(nostocrop)
  if ~autocropall(true,nostocrop)
    if ismember(saxno,nostocrop)
      myfailed('Need to crop short-axis stack in order to do LV segmentation');
      return
    end
  elseif SET(saxno).XSize * SET(saxno).ResolutionX > 250 || ...
      SET(saxno).YSize * SET(saxno).ResolutionY > 250
    myfailed('Need to crop short-axis stack more in order to do LV segmentation');
    autocropall(true,nostocrop)
    if SET(saxno).XSize * SET(saxno).ResolutionX > 250 || ...
        SET(saxno).YSize * SET(saxno).ResolutionY > 250
      myfailed('Need to crop short-axis stack more in order to do LV segmentation');
      return
    end
  end
end

%Open LV wizard GUI
gui = mygui('lvsegmentation.fig');
DATA.GUI.LVSegmentation = gui;

gui.saxno = saxno;
gui.ch4no = ch4no;
gui.ch3no = ch3no;
gui.ch2no = ch2no;

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
[x,y] = findfunctions('findlvcenter',gui.saxno,SET(gui.saxno).StartSlice:SET(gui.saxno).EndSlice);
SET(gui.saxno).CenterX=x;
SET(gui.saxno).CenterY=y;

[rows,cols] = calcfunctions('calcrowscols',gui.saxno);
gui.matrix = [rows cols]; %viewmatrix for shortaxis image

gui.saxim = calcfunctions('calcmontageviewim',gui.saxno,gui.matrix);
gui.midim = calcfunctions('remapuint8',SET(gui.saxno).IM(:,:,:,gui.midslice),gui.saxno);
emptyim = uint8(zeros(1,1,SET(gui.saxno).TSize));
gui.ch4im = emptyim;
gui.ch3im = emptyim;
gui.ch2im = emptyim;
if ~isempty(gui.ch4no)
  tfs = round(linspace(1,SET(gui.ch4no).TSize,SET(gui.saxno).TSize));
  gui.ch4im = calcfunctions('remapuint8',SET(gui.ch4no).IM(:,:,tfs),gui.ch4no);
end
if ~isempty(gui.ch3no)
  tfs = round(linspace(1,SET(gui.ch3no).TSize,SET(gui.saxno).TSize));
  gui.ch3im = calcfunctions('remapuint8',SET(gui.ch3no).IM(:,:,tfs),gui.ch3no);
end
if ~isempty(gui.ch2no)
  tfs = round(linspace(1,SET(gui.ch2no).TSize,SET(gui.saxno).TSize));
  gui.ch2im = calcfunctions('remapuint8',SET(gui.ch2no).IM(:,:,tfs),gui.ch2no);
end

%plot sax image stack
gui.handles.saximh = image(gui.saxim(:,:,gui.tf),'Parent',gui.handles.shortaxisaxes);
gui.handles.midimh = image(gui.midim(:,:,gui.tf),'Parent',gui.handles.midsliceaxes);

axh = [gui.handles.shortaxisaxes gui.handles.midsliceaxes ...
  gui.handles.ch4axes gui.handles.ch3axes gui.handles.ch2axes];
colormap(axh(1),gray(255));
colormap(axh(2),gray(255));
hold(axh(1),'on');
hold(axh(2),'on');
set(gui.handles.saximh,'ButtonDownFcn', ...
  'lvsegmentation(''shortaxisaxes_Buttondown'')');
set(gui.handles.midimh,'ButtonDownFcn', ...
  'lvsegmentation(''midsliceaxes_Buttondown'')');
 
%plot long-axis image stacks
if ~isempty(gui.ch4no)
gui.handles.ch4imh = image(gui.ch4im(:,:,gui.tf),'Parent',gui.handles.ch4axes);
colormap(axh(3),gray(255));
hold(axh(3),'on');
set(gui.handles.ch4imh,'ButtonDownFcn', ...
  'lvsegmentation(''ch4axes_Buttondown'')');
else
  set(gui.handles.ch4axes,'visible','off');
end
if ~isempty(gui.ch3no)
gui.handles.ch3imh = image(gui.ch3im(:,:,gui.tf),'Parent',gui.handles.ch3axes);
colormap(axh(4),gray(255));
hold(axh(4),'on');
set(gui.handles.ch3imh,'ButtonDownFcn', ...
  'lvsegmentation(''ch3axes_Buttondown'')');
else
  set(gui.handles.ch3axes,'visible','off');
end
if ~isempty(gui.ch2no)
gui.handles.ch2imh = image(gui.ch2im(:,:,gui.tf),'Parent',gui.handles.ch2axes);
colormap(axh(5),gray(255));
hold(axh(5),'on');
set(gui.handles.ch2imh,'ButtonDownFcn', ...
  'lvsegmentation(''ch2axes_Buttondown'')');
else
  set(gui.handles.ch2axes,'visible','off');
end

axis(axh,'off');
axis(axh,'equal');

%Initiate plot handles
gui.handles.startsliceframe = plot(gui.handles.shortaxisaxes,[],[]);
gui.handles.endsliceframe = plot(gui.handles.shortaxisaxes,[],[]);
gui.handles.startslice4ch = plot(gui.handles.ch4axes,[],[]);
gui.handles.endslice4ch = plot(gui.handles.ch4axes,[],[]);
gui.handles.startslice3ch = plot(gui.handles.ch3axes,[],[]);
gui.handles.endslice3ch = plot(gui.handles.ch3axes,[],[]);
gui.handles.startslice2ch = plot(gui.handles.ch2axes,[],[]);
gui.handles.endslice2ch = plot(gui.handles.ch2axes,[],[]);
gui.handles.center = plot(gui.handles.shortaxisaxes,[],[]);
gui.handles.centermid = plot(gui.handles.midsliceaxes,[],[]);
gui.handles.center4ch = plot(gui.handles.ch4axes,[],[]);
gui.handles.center3ch = plot(gui.handles.ch3axes,[],[]);
gui.handles.center2ch = plot(gui.handles.ch2axes,[],[]);

set(gui.fig,'KeyPressFcn', ...
  @(hObject,eventdata)lvsegmentation('keypressed',hObject,eventdata));

updateselectedslices;
updatecenterpoint;

%-------------------
function inittimebar
%-------------------
%Initiate timebar axis
global DATA SET
gui = DATA.GUI.LVSegmentation;

h = gui.handles.timebaraxes;
no = gui.saxno;
if isempty(no)
  return
end
set(h,'ButtonDownFcn', ...
  @(hObject,eventdata)lvsegmentation('timebaraxes_ButtonDownFcn',hObject,eventdata));

delete(get(h,'Children'));

tvec = SET(no).TimeVector;

hold(h,'on');
fcn = @(hObject,eventdata)lvsegmentation('timebar_ButtonDownFcn',hObject,eventdata);

%Draw timebar (red) and set its buttondown fcn
gui.handles.timebar = plot(h,tvec(gui.tf)*[1 1],[0 1],'b','Tag','currenttime');
set(gui.handles.timebar,'ButtonDownFcn',fcn);

%Set axes options
tvec = SET(no).TimeVector;
tmin = tvec(1);
tmax = tvec(end);

marg = (tmax-tmin)/500;
axis(h,[tmin-marg tmax+marg 0 1]);
%tstep = 5*(tmax-tmin)/(numel(tvec)-1);
%tickvec = [ceil(tmin/tstep)*tstep:tstep:floor(tmax/tstep)*tstep];
%tickvec = tvec(1:5:end);
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
gui = DATA.GUI.LVSegmentation;

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
gui = DATA.GUI.LVSegmentation;
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
gui = DATA.GUI.LVSegmentation;
no = gui.saxno;
obj = hObject;
motionfcn = @(hObject,eventdata)lvsegmentation('timebaraxes_MotionFcn',hObject,eventdata,obj,no);
set(gui.fig,'WindowButtonMotionFcn',motionfcn);
buttonupfcn = @(hObject,eventdata)lvsegmentation('timebaraxes_ButtonUpFcn',hObject,eventdata);
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
DATA.GUI.LVSegmentation.tf = max(min(tf,SET(no).TSize),1);
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
gui = DATA.GUI.LVSegmentation;

%Update timeframe and images
gui.tf = round(mygetvalue(gui.handles.timeslider));
set(gui.handles.timeslider,'Value',gui.tf);
updateimages;

%--------------------
function updateimages
%--------------------
%Update all image stacks and also timebar
global DATA SET
gui = DATA.GUI.LVSegmentation;
set(gui.handles.saximh,'CData',gui.saxim(:,:,gui.tf));
set(gui.handles.midimh,'CData',gui.midim(:,:,gui.tf));
if isfield(gui.handles,'ch4imh')
  set(gui.handles.ch4imh,'CData',gui.ch4im(:,:,gui.tf));
end
if isfield(gui.handles,'ch3imh')
  set(gui.handles.ch3imh,'CData',gui.ch3im(:,:,gui.tf));
end
if isfield(gui.handles,'ch2imh')
  set(gui.handles.ch2imh,'CData',gui.ch2im(:,:,gui.tf));
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
gui = DATA.GUI.LVSegmentation;
gui.tf = mod(gui.tf-2,SET(gui.saxno).TSize)+1;
updateimages;

%---------------------
function next_Callback
%---------------------
%Callback for next timeframe pushbutton
global DATA SET
gui = DATA.GUI.LVSegmentation;
gui.tf = mod(gui.tf,SET(gui.saxno).TSize)+1;
updateimages;

%---------------------
function play_Callback
%---------------------
%Callback for play togglebutton
global DATA SET
gui = DATA.GUI.LVSegmentation;

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
gui = DATA.GUI.LVSegmentation;

%Clean up before drawing new frames and lines
delete([gui.handles.startsliceframe gui.handles.endsliceframe ...
  gui.handles.startslice4ch gui.handles.endslice4ch ...
  gui.handles.startslice3ch gui.handles.endslice3ch ...
  gui.handles.startslice2ch gui.handles.endslice2ch]);

%circlex = cos(linspace(0,2*pi,100))*2*SET(gui.saxno).XSize/gui.matrix(1);
%circley = sin(linspace(0,2*pi,100))*2*SET(gui.saxno).YSize/gui.matrix(2);
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
newmid = floor((SET(gui.saxno).StartSlice + SET(gui.saxno).EndSlice)/2);
if newmid ~= gui.midslice
  gui.midslice = newmid;
  gui.midim = calcfunctions('remapuint8',SET(gui.saxno).IM(:,:,:,gui.midslice));
  set(gui.handles.midimh,'Cdata',gui.midim(:,:,gui.tf));
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

if ~isempty(gui.ch3no)
  [ssx,ssy] = calcfunctions('calcplaneintersections', ...
    gui.ch3no,gui.saxno,'one','one',SET(gui.saxno).StartSlice);
  [esx,esy] = calcfunctions('calcplaneintersections', ...
    gui.ch3no,gui.saxno,'one','one',SET(gui.saxno).EndSlice);
  gui.handles.startslice3ch = plot( ...
    gui.handles.ch3axes,ssy,ssx,'b');
  gui.handles.endslice3ch = plot( ...
    gui.handles.ch3axes,esy,esx,'r');
end

if ~isempty(gui.ch2no)
  [ssx,ssy] = calcfunctions('calcplaneintersections', ...
    gui.ch2no,gui.saxno,'one','one',SET(gui.saxno).StartSlice);
  [esx,esy] = calcfunctions('calcplaneintersections', ...
    gui.ch2no,gui.saxno,'one','one',SET(gui.saxno).EndSlice);
  gui.handles.startslice2ch = plot( ...
    gui.handles.ch2axes,ssy,ssx,'b');
  gui.handles.endslice2ch = plot( ...
    gui.handles.ch2axes,esy,esx,'r');
end

segment('updateselectedslices');


%------------------------------
function recalculatecenterpoint
%------------------------------
%Recalculate center ponit based on new LV slices

global SET DATA
gui = DATA.GUI.LVSegmentation;
[x,y] = findfunctions('findlvcenter',gui.saxno,SET(gui.saxno).StartSlice:SET(gui.saxno).EndSlice);
SET(gui.saxno).CenterX=x;
SET(gui.saxno).CenterY=y;


%-------------------------
function updatecenterpoint
%-------------------------
%Update center point
global DATA SET
gui = DATA.GUI.LVSegmentation;

%Clean up before drawing
delete([gui.handles.center gui.handles.centermid ...
  gui.handles.center4ch gui.handles.center3ch gui.handles.center2ch]);

startxyz = [SET(gui.saxno).CenterX; ...
  SET(gui.saxno).CenterY; ...
  SET(gui.saxno).StartSlice];
endxyz = [startxyz(1:2); SET(gui.saxno).EndSlice];

crosscolor = [1 0.6 0.2]; %DATA.centercrossdef
[tempx,tempy] = ndgrid(0:(gui.matrix(1)-1),0:(gui.matrix(2)-1));
gui.handles.center = plot(...
  gui.handles.shortaxisaxes,...
  startxyz(2)+tempy(:)*SET(gui.saxno).YSize,...
  startxyz(1)+tempx(:)*SET(gui.saxno).XSize,'+',...
  'color',crosscolor);
gui.handles.centermid = plot( ...
  gui.handles.midsliceaxes,startxyz(2),startxyz(1),'+','Color',crosscolor);

%Get center point coordinates in longaxis images
[saxpos,saxmat] = calcfunctions('calcormat',gui.saxno);
truecenter = [saxpos saxpos] + saxmat*([startxyz endxyz]-1);

if ~isempty(gui.ch4no)
  [ch4pos,ch4mat] = calcfunctions('calcormat',gui.ch4no);
  ch4center = 1 + ch4mat\(truecenter - [ch4pos ch4pos]);
  gui.handles.center4ch = plot( ...
    gui.handles.ch4axes,ch4center(2,:),ch4center(1,:),'Color',crosscolor);
end

if ~isempty(gui.ch3no)
  [ch3pos,ch3mat] = calcfunctions('calcormat',gui.ch3no);
  ch3center = 1 + ch3mat\(truecenter - [ch3pos ch3pos]);
  gui.handles.center3ch = plot( ...
    gui.handles.ch3axes,ch3center(2,:),ch3center(1,:),'Color',crosscolor);
end

if ~isempty(gui.ch2no)
  [ch2pos,ch2mat] = calcfunctions('calcormat',gui.ch2no);
  ch2center = 1 + ch2mat\(truecenter - [ch2pos ch2pos]);
  gui.handles.center2ch = plot( ...
    gui.handles.ch2axes,ch2center(2,:),ch2center(1,:),'Color',crosscolor);
end

%--------------------------------
function shortaxisaxes_Buttondown %#ok<DEFNU>
%--------------------------------
global DATA SET
gui = DATA.GUI.LVSegmentation;

%Extract coordinates clicked
[y,x] = mygetcurrentpoint(gui.handles.shortaxisaxes);
%Find slice
col = 1+floor((y-0.5)/SET(gui.saxno).YSize);
row = 1+floor((x-0.5)/SET(gui.saxno).XSize);
slice = col+(row-1)*gui.matrix(2);

if slice > SET(gui.saxno).ZSize
  return;
end
switch get(gui.fig,'SelectionType')
  case 'normal'
    if slice <= SET(gui.saxno).EndSlice
      SET(gui.saxno).StartSlice = slice;
    end
  case 'alt'
    if slice >= SET(gui.saxno).StartSlice
      SET(gui.saxno).EndSlice = slice;
    end
end
updateselectedslices;
recalculatecenterpoint;
updatecenterpoint;

%-------------------------------
function midsliceaxes_Buttondown
%-------------------------------
%Button down function for mid slice axes. Move center point.
global DATA SET
gui = DATA.GUI.LVSegmentation;

[y,x] = mygetcurrentpoint(gui.handles.midsliceaxes);
SET(gui.saxno).CenterX = x;
SET(gui.saxno).CenterY = y;
updatecenterpoint;

%--------------------------
function ch4axes_Buttondown
%--------------------------



%--------------------------
function ch3axes_Buttondown
%--------------------------



%--------------------------
function ch2axes_Buttondown
%--------------------------

%-------------------------------------
function keypressed(hObject,eventdata)
%-------------------------------------
%Keypress function for GUI
global DATA SET
gui = DATA.GUI.LVSegmentation;

switch eventdata.Key
  case 'uparrow'
    if isempty(eventdata.Modifier)
      SET(gui.saxno).StartSlice = max(1,SET(gui.saxno).StartSlice-1);
    else
      SET(gui.saxno).EndSlice = max(SET(gui.saxno).StartSlice,SET(gui.saxno).EndSlice-1);
    end
    updateselectedslices;
    recalculatecenterpoint;
    updatecenterpoint;
  case 'downarrow'
    if isempty(eventdata.Modifier)
      SET(gui.saxno).StartSlice = min(SET(gui.saxno).EndSlice,SET(gui.saxno).StartSlice+1);
    else
      SET(gui.saxno).EndSlice = min(SET(gui.saxno).ZSize,SET(gui.saxno).EndSlice+1);
    end
    updateselectedslices;
    recalculatecenterpoint;
    updatecenterpoint;
  case 'leftarrow'
    prev_Callback;
  case 'rightarrow'
    next_Callback;
  case 'p'
    set(gui.handles.playtogglebutton,'value',1-mygetvalue(gui.handles.playtogglebutton));
    play_Callback;
end

%---------------------
function dolv_Callback %#ok<DEFNU>
%---------------------
%Do LV segmentation in selected slices
global DATA
gui = DATA.GUI.LVSegmentation;
segment('switchtoimagestack',gui.saxno);
oldpol=DATA.ThisFrameOnly;
segment('framemode_Callback',2);%segment('thisframeonly_Callback',false);
lvpeter('segmentfullyautomatic_Callback');
close_Callback
drawfunctions('drawall');
figure(DATA.GUI.Segment.fig);
switch oldpol
  case 1
    segment('framemode_Callback',1)
  case 0
    segment('framemode_Callback',2)
end

%----------------------
function close_Callback  
%----------------------
%Close LV segmentation GUI
global DATA

try
  DATA.GUI.LVSegmentation = close(DATA.GUI.LVSegmentation);
catch   %#ok<CTCH>
  DATA.GUI.LVSegmentation=[];
  delete(gcbf);
end