function slicereport(type)
%GUI for graphical display of slice based regional function.
%Adapted for use with mygui class by Nils Lundahl
if nargin == 0
  type = 'init';
end
macro_helper(type);
feval(type); % FEVAL switchyard

%------------
function init
%------------
%Initiate GUI
global DATA SET NO

if SET(NO).TSize<2
  myfailed('Data needs to be time resolved.',DATA.GUI.Segment);
  return;
end

if isempty(SET(NO).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end

if sum(findfunctions('findslicewithendo',NO))==0
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

%Generate figure and GUI
gui = mygui('slicereport.fig');
DATA.GUI.SliceReport = gui;

temp = get(gui.handles.sectorslistbox,'String');
gui.numsectors = str2num(temp{mygetlistbox(gui.handles.sectorslistbox)}); %#ok<ST2NM>
gui.slice = SET(NO).CurrentSlice;
if SET(NO).ZSize > 1
  set(gui.handles.sliceslider,...
    'min',1,...
    'max',SET(NO).ZSize,...
    'value',SET(NO).ZSize-SET(NO).CurrentSlice+1,...
    'SliderStep',1/(SET(NO).ZSize-1)*[1 1]);  
  set(gui.handles.slicetext,'String',dprintf('Slice %d',gui.slice));
end
set(gui.handles.rotationslider,...
  'value',SET(NO).SectorRotation,...
  'min',-180,...
  'max',180,...
  'SliderStep',[3/360 3/360]);

if ~DATA.GUISettings.PointsEnabled
  set(gui.handles.rotationfromannotationcheckbox,'visible','off','value',0);
end

update; %Update more
try set(gui.figure1,'Visible','on'); catch, end

%----------------------
function close_Callback
%----------------------
%Close GUI
global DATA

try
  DATA.GUI.SliceReport=close(DATA.GUI.SliceReport);  %close the slice gui
catch %#ok<CTCH>
  delete(gcf)
end

DATA.GUI.SliceReport= [];

%--------------
function update
%--------------
%Update GUI
global DATA
gui = DATA.GUI.SliceReport;

if DATA.Record
  drawnow;
  DATA.MovieFrame = mygetframe(gui.fig);
  export('exportmovierecorder_Callback','newframe');
end
parameter;

%--------------
function sector
%--------------
%Update plots upon changed sector rotation

%Take from slider just in case
global DATA SET NO
gui = DATA.GUI.SliceReport;

SET(NO).SectorRotation = mygetvalue(gui.handles.rotationslider);

if get(gui.handles.rotationfromannotationcheckbox,'value')
  %Take sector rotation from insertion point
  px = NaN;
  py = NaN;
  pz = NaN;
  for loop=1:length(SET(NO).Point.Label)
    if isequal(SET(NO).Point.Label{loop},'Sector start')
      if ~isnan(px)
        mywarning('More than one sector start found, take ''latest''.',DATA.GUI.Segment);
      end
      px = SET(NO).Point.X(loop);
      py = SET(NO).Point.Y(loop);
      pz = SET(NO).Point.Z(loop);
    end
  end
  if isnan(px)
    myfailed('Could not find an annotation point labelled ''Sector start''',DATA.GUI.Segment);
    set(gui.handles.rotationfromannotationcheckbox,'value',0);
  else
    %--- Adjust slice rotation
    mx = mean(SET(NO).EndoX(:,SET(NO).EDT,pz));
    my = mean(SET(NO).EndoY(:,SET(NO).EDT,pz));
    SET(NO).SectorRotation = angle(complex(my-py,mx-px))*180/pi;
  end
  set(gui.handles.rotationslider,'value',...
    SET(NO).SectorRotation);
end

if gui.numsectors>DATA.Pref.RadialProfiles
  mywarning('Reporting in more sectors than evaluating. Reverting back to same number of sectors as evaluation the profile along. For details see preferences. ',DATA.GUI.Segment);
  gui.numsectors = DATA.Pref.RadialProfiles;
end

[gui.meanx,gui.meany,gui.sectors] = ...
  calcfunctions('findmeaninsectorslice',gui.type,DATA.NumPoints,SET(NO).CurrentTimeFrame,gui.slice,gui.numsectors);

%--- Plot play image
tf = SET(NO).CurrentTimeFrame;
axes(gui.handles.playaxes);
gui.handles.playimage = image(calcfunctions('remapuint8',SET(NO).IM(:,:,tf,gui.slice),NO));
colormap(gray(256));
axis image off;
hold on;
gui.handles.endocontour = line('parent',gui.handles.playaxes,'XData',SET(NO).EndoY(:,tf,gui.slice),'YData',SET(NO).EndoX(:,tf,gui.slice),'Color','r','LineStyle','-');

if ~isempty(SET(NO).EpiX)
  gui.handles.epicontour = line('parent',gui.handles.playaxes,'XData',SET(NO).EpiY(:,tf,gui.slice),'YData',SET(NO).EpiX(:,tf,gui.slice),'Color','g','LineStyle','-');
else
  gui.handles.epicontour = line('parent',gui.handles.playaxes,'XData',NaN,'YData',NaN);
end
hold off;

%--- Plot overview image
axes(gui.handles.imageaxes);
gui.handles.image = image(calcfunctions('remapuint8',SET(NO).IM(:,:,tf,gui.slice),NO));
axis image off;

hold on;
h = line('parent',gui.handles.imageaxes,'XData',SET(NO).EndoY(:,tf,gui.slice),'YData',SET(NO).EndoX(:,tf,gui.slice),'Color','r','LineStyle','-');
set(h,'linewidth',3);
if ~isempty(SET(NO).EpiX)
  h = line('parent',gui.handles.imageaxes,'XData',SET(NO).EpiY(:,tf,gui.slice),'YData',SET(NO).EpiX(:,tf,gui.slice),'Color','g','LineStyle','-');
else
  h = line('parent',gui.handles.imageaxes,'XData',NaN,'YData',NaN);
end

set(h,'linewidth',3);
if not(isnan(gui.meanx))
  for loop=1:gui.numsectors
    if isequal(gui.type,'epi')
      xposm = 0.5*(...
        SET(NO).EpiY(gui.sectors(loop),tf,gui.slice)+...
        SET(NO).EpiY(gui.sectors(loop+1),tf,gui.slice));
      yposm = 0.5*(...
        SET(NO).EpiX(gui.sectors(loop),tf,gui.slice)+...
        SET(NO).EpiX(gui.sectors(loop+1),tf,gui.slice));
      xpos = SET(NO).EpiY(gui.sectors(loop),tf,gui.slice);
      ypos = SET(NO).EpiX(gui.sectors(loop),tf,gui.slice);
    else
      xposm = 0.5*(...
        SET(NO).EndoY(gui.sectors(loop),tf,gui.slice)+...
        SET(NO).EndoY(gui.sectors(loop+1),tf,gui.slice));
      yposm = 0.5*(...
        SET(NO).EndoX(gui.sectors(loop),tf,gui.slice)+...
        SET(NO).EndoX(gui.sectors(loop+1),tf,gui.slice));
      xpos = SET(NO).EndoY(gui.sectors(loop),tf,gui.slice);
      ypos = SET(NO).EndoX(gui.sectors(loop),tf,gui.slice);
    end
    line('parent',gui.handles.imageaxes,'XData',[gui.meany xpos],'YData',[gui.meanx ypos],'Color','y','LineStyle','-');
    if (loop==1)
      line('parent',gui.handles.imageaxes,'XData',[xpos xpos-(gui.meany-xpos)],'YData',[ypos ypos-(gui.meanx-ypos)],'Color','y','LineStyle','-');
    end
    if (loop==1)||(gui.numsectors<20)
      h = text(xposm,yposm,sprintf('%d',loop));
      set(h,'color','y','fontsize',14);
    end
  end %loop
end %if isnan
hold off;

%-------------
function slice
%-------------
%Different slice, update a lot of things
global DATA SET NO
gui = DATA.GUI.SliceReport;

gui.slice = round(mygetvalue(gui.handles.sliceslider));
gui.slice = SET(NO).ZSize-gui.slice+1;
SET(NO).CurrentSlice = gui.slice;
gui.slice = min(max(gui.slice,1),SET(NO).ZSize);
set(gui.handles.slicetext,'String',dprintf('Slice %d',gui.slice));
set(gui.handles.sliceslider,'value',SET(NO).ZSize-gui.slice+1);
sector;
parameter;

%-----------------
function parameter
%-----------------
%Update things upon changed selection in parameter listbox

global DATA SET NO
gui = DATA.GUI.SliceReport;

t = (0:(SET(NO).TSize-1))*SET(NO).TIncr;
switch mygetlistbox(gui.handles.parameterlistbox)
  case 1
    %Wallthickness
    if isempty(SET(NO).EpiX)
      myfailed('No LV epicardium available.',DATA.GUI.Segment);
      set(gui.handles.parameterlistbox,'value',3);
      parameter;
    else
      gui.type = 'epi';
      wallthickness = calcfunctions('calcwallthickness',gui.numsectors,NO);
      gui.outdata = squeeze(wallthickness(:,gui.slice,:));
      gui.title = dprintf('Wall thickness');
      gui.outunit = 'mm';
    end
  case 2
    %Fractional thickening
    if isequal(SET(NO).EDT,SET(NO).EST)
      mywarning('Warning, end-diastole occurs at the same time as end-systole. Use autodetect under edit menu.',DATA.GUI.Segment);
    end
    if isempty(SET(NO).EpiX)
      myfailed('No LV epicardium available.',DATA.GUI.Segment);
      set(gui.handles.parameterlistbox,'value',3);
      parameter;
    else
      gui.type = 'epi';
      wallthickness = calcfunctions('calcwallthickness',gui.numsectors,NO);
      gui.outdata = squeeze(wallthickness(:,gui.slice,:));
      minthick = repmat(gui.outdata(:,SET(NO).EDT),[1 SET(NO).TSize]);
      gui.outdata = (gui.outdata-minthick)./minthick;
      gui.outdata = gui.outdata*100;
      gui.title = dprintf('Fractional wallthickening');
      gui.outunit = '%';
    end
  case 3
    %Radial velocity
    gui.type = 'endo';
    radvel = calcfunctions('calcradialvelocity',NO);
    if not(isnan(radvel(1,1,gui.slice)))
      gui.outdata = calcfunctions('findmeaninsector',....
        'endo',radvel,gui.slice,gui.numsectors);
    else
      gui.outdata = [];
    end
    gui.title = dprintf('Radial velocity');
    gui.outunit = 'cm/s';
  case 4
    %Radius
    gui.type = 'endo';
    rad = calcfunctions('calcendoradius',NO);
    if not(isnan(rad(1,1,gui.slice)))
      gui.outdata = calcfunctions('findmeaninsector',....
        'endo',rad,gui.slice,gui.numsectors);
    else
      gui.outdata = [];
    end
    gui.title = dprintf('Radius');
    gui.outunit = 'mm';
  otherwise
    myfailed('Unknown option to reportslice.',DATA.GUI.Segment);
end
gui.outdata = squeeze(gui.outdata);
axes(gui.handles.plotaxes);
if not(isempty(gui.outdata))  
  h = plot(repmat(t,size(gui.outdata,1),1)',gui.outdata');
  set(h,'linewidth',2,'marker','.');
  %hold on;
  %plot([t(1),t(end)],[0 0],'k:');
  %hold off;
  
  title(sprintf('%s [%s]',gui.title,gui.outunit),'color',DATA.GUISettings.ForegroundColor);
  
  %Fix legend
  if gui.numsectors<20
    l = cell(1,gui.numsectors);
    for loop=1:gui.numsectors
      l{loop} = dprintf('Sector %d',loop);
    end
    legend(l{:});
  end
  
  ylabel(sprintf('[%s]',gui.outunit),'color',DATA.GUISettings.ForegroundColor);
  set(gui.handles.plotaxes,'xcolor',DATA.GUISettings.ForegroundColor);
  set(gui.handles.plotaxes,'ycolor',DATA.GUISettings.ForegroundColor);
else
  plot('parent',gui.handles.plotaxes,'XData',t,'YData',zeros(size(t)));
  title(dprintf('%s not available.',gui.title),'color',DATA.GUISettings.ForegroundColor);
end
xlabel(dprintf('Time'));
sector;

%-------------------
function frameupdate
%-------------------
%Different timeframe, update less.
global DATA SET NO
gui = DATA.GUI.SliceReport;

set(gui.handles.playimage,'cdata',calcfunctions('remapuint8',SET(NO).IM(:,:,SET(NO).CurrentTimeFrame,gui.slice),NO));
set(gui.handles.endocontour,...
  'xdata',SET(NO).EndoY(:,SET(NO).CurrentTimeFrame,gui.slice),...
  'ydata',SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,gui.slice));
if ~isempty(SET(NO).EpiX)
  set(gui.handles.epicontour,...
    'xdata',SET(NO).EpiY(:,SET(NO).CurrentTimeFrame,gui.slice),...
    'ydata',SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,gui.slice));
end
if DATA.Record
  drawnow;
  DATA.MovieFrame = mygetframe(gui.fig);
  export('exportmovierecorder_Callback','newframe');
end

%------------
function play
%------------
%Callback for video playback

global DATA SET NO
gui = DATA.GUI.SliceReport;

DATA.StartFrame = SET(NO).CurrentTimeFrame;
DATA.StartTime = now;
try
  while get(gui.handles.playtogglebutton,'value')
    SET(NO).CurrentTimeFrame = segment('getframenumber');
    pause(0.5*SET(NO).BeatTime/SET(NO).TSize);
    %pause(0.01);
    frameupdate;
  end
catch %#ok<CTCH>
  %Do nothing
end

%-----------------------
function export_Callback
%-----------------------
%Export data to clipboard. Renamed to avoid confusion with export.m
global DATA SET NO
gui = DATA.GUI.SliceReport;

t = (0:(SET(NO).TSize-1))*SET(NO).TIncr;
stri = sprintf('%s\tSlice:%d\n',SET(NO).PatientInfo.Name,gui.slice);
stri = [stri sprintf('%s\t%s\n',gui.title,gui.outunit)];
stri = [stri sprintf('Time\t')];
for sector=1:gui.numsectors
  stri = [stri sprintf('Sector%d\t',sector)]; %#ok<AGROW>
end
stri = [stri sprintf('\n')];
for tloop=1:SET(NO).TSize
  stri = [stri sprintf('%0.5g\t',t(tloop))]; %#ok<AGROW>
  for sector=1:gui.numsectors
    stri = [stri sprintf('%0.5g\t',gui.outdata(sector,tloop))]; %#ok<AGROW>
  end
  stri = [stri sprintf('\n')]; %#ok<AGROW>
end
clipboard('copy',stri);
mymsgbox('Data Exported.','Done!',DATA.GUI.Segment);

%------------
function next
%------------
%Callback to go to next frame
global SET NO
SET(NO).CurrentTimeFrame=SET(NO).CurrentTimeFrame+1;
if SET(NO).CurrentTimeFrame>SET(NO).TSize
  SET(NO).CurrentTimeFrame=1;
end
frameupdate;

%------------
function prev
%------------
%Callback to go to previous frame
global SET NO
SET(NO).CurrentTimeFrame=SET(NO).CurrentTimeFrame-1;
if SET(NO).CurrentTimeFrame<1
  SET(NO).CurrentTimeFrame=SET(NO).TSize;
end
frameupdate;
      
%------------
function wall
%------------
%Callback upon changing frac/wall radiobutton selection to wall
global DATA
gui = DATA.GUI.SliceReport;

set(gui.handles.wallradiobutton,'value',1);
set(gui.handles.fracradiobutton,'value',0);
slice;

%------------
function frac
%------------
%Callback upon changing frac/wall radiobutton selection to frac
global DATA
gui = DATA.GUI.SliceReport;

set(gui.handles.wallradiobutton,'value',0);
set(gui.handles.fracradiobutton,'value',1);
slice;

%---------------
function sectors
%---------------
%Update upon changing number of sectors
global DATA
gui = DATA.GUI.SliceReport;

temp = get(gui.handles.sectorslistbox,'String');
gui.numsectors = str2num(temp{mygetlistbox(gui.handles.sectorslistbox)}); %#ok<ST2NM>
slice;
