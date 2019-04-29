function [varargout] = reportroi(varargin)
%GUI for ROI analysis.
%Adapted for use with mygui class by Nils Lundahl

macro_helper(varargin{:});
if nargin == 0 || isempty(varargin{1})
  varargin{1} = 'init';
end
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init
%------------
global DATA SET NO

if SET(NO).RoiN<1
  myfailed('No ROIs available. Use blue pen to draw ROIs.',DATA.GUI.Segment);
  return;
end;

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

if isempty(SET(NO).EndAnalysis)
  SET(NO).EndAnalysis = SET(NO).TSize;
end;

if isempty(SET(NO).StartAnalysis)
  SET(NO).StartAnalysis = 1;
end;

%if (SET(NO).EndAnalysis-SET(NO).StartAnalysis)<1
%  myfailed('No timespan, adjust Start/End analysis time (under stress menu or drag red bar).',DATA.GUI.Segment);
%  return;
%
%end;

%Open gui
%check if open
if isopengui('roi.fig')
  close_Callback;
end
gui = mygui('roi.fig');
DATA.GUI.ROI = gui;

if ~isrectilinear(SET(NO).TimeVector)
  gui.rectilinear = false;
end;

%  mywarning(['Non uniform time steps detected. Smoothing and display of '...
%    'min/max slope disabled.']);
%  set([gui.handles.smoothcheckbox ...
%    gui.handles.showslopescheckbox],'enable','off');
%end

%Ask for what ROI's to include
[gui.rois,~,gui.normalized] = roi('roiselector','',...
  get(gui.handles.thissliceonlycheckbox,'value'),...
  '',get(gui.handles.normalizedcheckbox,'value'));

%Initalize variables
gui.no = NO;
gui.t = SET(NO).TimeVector;
gui.t = gui.t*1000; %ms
gui.roiint = nan([SET(NO).TSize length(gui.rois)]);
gui.roistd = gui.roiint;
gui.roisizepix = gui.roiint;
gui.roisize = gui.roiint;
gui.roimin = gui.roiint;
gui.roimax = gui.roiint;
gui.sigma = mygetvalue(gui.handles.smoothslider);
gui.offset = 0;
gui.detrendstartms = 0;
gui.detrendendvalue = 0;
set(gui.handles.sigmaedit,'string',sprintf('%0.5g',gui.sigma));

recalc; %do plot

% datacursormode(gui.fig)

% sprintf('fig nr: %0.3g \n', gui.fig)
%set pointer, eriks
set(gui.fig,'Pointer','arrow');

%------------------
function smoothedit
%------------------
global DATA
gui = DATA.GUI.ROI;

s = mygetedit(gui.handles.sigmaedit);
n = str2double(s);
if not(isnan(n))
  gui.sigma = n;
end;
gui.sigma = min(max(gui.sigma,get(gui.handles.smoothslider,'min')),get(gui.handles.smoothslider,'max'));
set(gui.handles.sigmaedit,'string',sprintf('%0.5g',gui.sigma));
set(gui.handles.smoothslider,'value',gui.sigma);

%--------------------
function smoothslider
%--------------------
global DATA
gui = DATA.GUI.ROI;

gui.sigma = mygetvalue(gui.handles.smoothslider);
set(gui.handles.sigmaedit,'string',sprintf('%0.5g',gui.sigma));
doplot;
    
%------------------
function offsetedit %#ok<DEFNU>
%------------------
global DATA
gui = DATA.GUI.ROI;

temp = mygetedit(gui.handles.offsetedit);
[num,ok] = str2num(temp); %#ok<ST2NM>
if ok
  gui.offset = num;
  set(gui.handles.offsetslider,'value',num);
else
  myfailed('Could not interpret number.',DATA.GUI.Segment);
  return;
end;
doplot;
    
%--------------------
function offsetslider %#ok<DEFNU>
%--------------------
%User adjusted offsetslider

global DATA
gui = DATA.GUI.ROI;

gui.offset = mygetvalue(gui.handles.offsetslider);
set(gui.handles.offsetedit,'string',sprintf('%0.5g',gui.offset));
doplot;

%---------------------------
function setdetrend_Callback %#ok<DEFNU>
%---------------------------
%User sets detrend information

global DATA
gui = DATA.GUI.ROI;

s.DetrendStart_ms = gui.detrendstartms;
s.DetrendEndValue = gui.detrendendvalue;
[s,ok] = inputstruct(s);
if ~ok
  myfailed('Aborted.');
  return;
end;

gui.detrendstartms = s.DetrendStart_ms;
gui.detrendendvalue = s.DetrendEndValue;

%Call plotting routine
doplot;

%---------------------------------
function parameterlistbox_Callback %#ok<DEFNU>
%---------------------------------
%Callback from parameter listbox
doplot;

%---------------------------
function showminmax_Callback %#ok<DEFNU>
%---------------------------
%Callback for show min/max checkbox
doplot;

%---------------------------
function showslopes_Callback %#ok<DEFNU>
%---------------------------
%Callback for show min/max slopes checkbox
doplot;

%-------------------------
function showfwhm_Callback %#ok<DEFNU>
%-------------------------
%Callback for show FWHM checkbox
doplot;

%----------------------------------
function showcentergravity_Callback %#ok<DEFNU>
%----------------------------------
%Callback for show min/max checkbox
doplot;

%-----------------------
function smooth_Callback %#ok<DEFNU>
%-----------------------
%Callback for smooth checkbox
doplot;


%-------------------------
function cropzero_Callback %#ok<DEFNU>
%-------------------------
%Callback for crop at zero checkbox
doplot;


%------------------------
function detrend_Callback %#ok<DEFNU>
%------------------------
%Callback for detrend checkbox
doplot;


%-------------------------------
function roibar_Buttondown(type) %#ok<DEFNU>
%-------------------------------
%Button down function for pressing timebar in plot flow gui. Sets button 
%motion and button up function for changing start time and end time for 
%calculation of flow.
global SET NO

SET(NO).Stress.temphandle = gcbo;
set(gcf,'WindowButtonMotionFcn','reportroi(''roibar_Motion'')');
set(gcf,'WindowButtonUpFcn',sprintf(...
  'reportroi(''roibar_Buttonup'',''%s'')',type));

%-----------------------------
function roibar_Buttonup(type) %#ok<DEFNU>
%-----------------------------
%Button up function called after moving timebar in plot flow gui. 
%Changes start time and end time for calculation of flow.
global SET NO

set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');

%Get Convert to timeframe
x = get(SET(NO).Stress.temphandle,'xdata');
x = x(1);
x = round(1+x/(1000*SET(NO).TIncr));

switch type
  case 'startbar'
    SET(NO).StartAnalysis = x;
  case 'endbar'
    SET(NO).EndAnalysis = x;    
end;

if SET(NO).StartAnalysis>SET(NO).EndAnalysis
  temp = SET(NO).StartAnalysis;
  SET(NO).StartAnalysis=SET(NO).EndAnalysis;
  SET(NO).EndAnalysis=temp;
end;

%Do graphical update
doplot;

%---------------------
function roibar_Motion %#ok<DEFNU>
%---------------------
%Button motion function called when moving timebar in plot flow gui. 
%Changes start time and end time for calculation of flow.
global SET NO

[x,y] = mygetcurrentpoint(gca);

%Convert to timeframes and round
x = round(1+x/(1000*SET(NO).TIncr));
x = min(max(1,x),SET(NO).TSize);

%Convert back to ms
x = (x-1)*SET(NO).TIncr*1000;

%Update display
set(SET(NO).Stress.temphandle,'xdata',[x x]);

%------------------
function plotfilter %#ok<DEFNU>
%------------------
%get sigma from slider
global DATA
gui = DATA.GUI.ROI;

gui.sigma = mygetvalue(gui.handles.smoothslider);

%create filter
x = linspace(-60,60,121);
f = exp(-x.^2/gui.sigma.^2);
figure(22);
plot(x,f);
xlabel('Timeframes');
title('Smoothing applicability.');

%--------------
function recalc
%--------------
global DATA SET
gui = DATA.GUI.ROI;

%--- Ask for what ROI's to include
[gui.rois,~,gui.normalized] = roi('roiselector','',...
  get(gui.handles.thissliceonlycheckbox,'value'),...
  '',get(gui.handles.normalizedcheckbox,'value'));

%--- Extract data
h = mywaitbarstart(length(gui.rois),'Please wait, calculating data',[],DATA.GUI.Segment);
for rloop=1:length(gui.rois)
  z = SET(gui.no).Roi(gui.rois(rloop)).Z;
  for tloop=SET(gui.no).Roi(gui.rois(rloop)).T
    if not(isnan(SET(gui.no).Roi(gui.rois(rloop)).X(1,tloop)))
      if gui.normalized
        temp = SET(gui.no).IM(:,:,tloop,z);
      else
        temp = calcfunctions('calctruedata',SET(gui.no).IM(:,:,tloop,z),gui.no);
      end;
      roimask = segment('createmask',...
        [SET(gui.no).XSize SET(gui.no).YSize],...
        SET(gui.no).Roi(gui.rois(rloop)).Y(:,tloop),...
        SET(gui.no).Roi(gui.rois(rloop)).X(:,tloop));
      ind = find(roimask);
      if ~isempty(ind)
        gui.roimin(tloop,rloop) = min(temp(ind));
        gui.roimax(tloop,rloop) = max(temp(ind));
        gui.roiint(tloop,rloop) = mean(temp(ind));
        gui.roistd(tloop,rloop) = std(temp(ind));
      end;
      gui.roisize(tloop,rloop) = (1/100)*polyarea(...
        SET(gui.no).ResolutionY*SET(gui.no).Roi(gui.rois(rloop)).Y(:,tloop),...
        SET(gui.no).ResolutionX*SET(gui.no).Roi(gui.rois(rloop)).X(:,tloop));
      gui.roisizepix(tloop,rloop) = sum(roimask(:))*SET(gui.no).ResolutionX*SET(gui.no).ResolutionY/100;
    end;
  end;
  h = mywaitbarupdate(h);
end;
mywaitbarclose(h);
doplot;
set(gui.fig,'Pointer','Arrow');

%--------------
function doplot
%--------------
global DATA SET
gui = DATA.GUI.ROI;

%--- Check what to plot
switch mygetlistbox(gui.handles.parameterlistbox)
  case 1
    gui.outdata = gui.roiint;
    gui.outname = translation.dictionary('Mean Intensity');
    if (SET(gui.no).TSize==1) %sbt20160516, mod EH
        auc = gui.outdata;%sbt20160516
    elseif (numel(gui.t)>1)%sbt20160516
      for lloop=1:length(gui.rois)
        auc(lloop) = trapz(gui.t/1000,gui.outdata(:,lloop));
      end
    else%sbt20160516
        myfailed('Error: unknown number of timeframes found.')%sbt20160516
        auc = 0;%sbt20160516
    end%sbt20160516
    set(gui.handles.auctext,'Visible','on');
    aucstring = sprintf('Area under curve:\n');
    for lloop=1:length(gui.rois)
      aucstring = [aucstring dprintf('%s:  %0.2f\n',SET(gui.no).Roi(gui.rois(lloop)).Name,auc(lloop))];
    end
    set(gui.handles.auctext,'String',aucstring);
  case 2
    gui.outdata = gui.roistd;
    gui.outname = translation.dictionary('Standard deviation');
    set(gui.handles.auctext,'Visible','off');
  case 3
    gui.outdata = gui.roisize;
    gui.outname = translation.dictionary('Area [cm^2]');
    set(gui.handles.auctext,'Visible','off');
  case 4
    gui.outdata = gui.roisizepix;
    gui.outname = translation.dictionary('Area based on pixels [cm^2] ');
    set(gui.handles.auctext,'Visible','off');
  case 5
    gui.outdata = gui.roimax;
    gui.outname = translation.dictionary('Maximum intensity');
    set(gui.handles.auctext,'Visible','off');
  case 6
    gui.outdata = gui.roimin;
    gui.outname = translation.dictionary('Minimum intensity');
    set(gui.handles.auctext,'Visible','off');
  otherwise
    myfailed('Not yet implemented.',DATA.GUI.Segment);
    gui.outdata = [];
    gui.outname = '';
end;

%Add offset
gui.outdata = gui.outdata+gui.offset;

%Crop data
if get(gui.handles.cropzerocheckbox,'value')
  gui.outdata(gui.outdata<0) = 0;
end;

%Detrend data
if get(gui.handles.detrendcheckbox,'value')
  t = SET(gui.no).TimeVector*1000;
  t(t < gui.detrendstartms) = NaN;
  t = t - min(t);
  t = t./max(t); %Normalize
  t(isnan(t)) = 0;
  detrend = t*gui.detrendendvalue;
  detrend = detrend(:); %reshape
  detrend = repmat(detrend,[1 size(gui.outdata,2)]); %reshape
  gui.outdata = gui.outdata-detrend;
end;

%--- Do smoothing
if dosmooth
  %get sigma from slider
  gui.sigma = mygetvalue(gui.handles.smoothslider);
  
  %create filter
  x = linspace(-60,60,121);
  f = exp(-x.^2/gui.sigma.^2);
  
  gui.smoothoutdata = gui.outdata;
  for rloop=1:size(gui.outdata,2)
    
    %Extract signal
    signal = gui.outdata(:,rloop);
        
    if ~gui.rectilinear
      signalfiltered = smooth_helper(gui.sigma,SET(gui.no).TimeVector,signal);
      %signalfiltered = conv2(signal,f','same')./(eps+conv2(ones(size(signal)),f','same'));
    else
      %Normalized averaging, rectilinear code
      signalfiltered = conv2(signal,f','same')./(eps+conv2(ones(size(signal)),f','same'));      
    end;
    
    %Store
    gui.smoothoutdata(:,rloop) = signalfiltered;        
    
  end;
  
else
  gui.smoothoutdata = gui.outdata;
  gui.t = SET(gui.no).TimeVector*1000;
end;

%--- Plot it
  
if SET(gui.no).TSize>1
  
  %Init
  gui.maxv = zeros(1,length(gui.rois));
  gui.minv = gui.maxv;
  gui.maxvd = gui.maxv;
  gui.minvd = gui.maxv;
  gui.maxind = gui.maxv;
  gui.minind = gui.maxv;
  gui.maxindd = gui.maxv;
  gui.minindd = gui.maxv;
  gui.fwhm = gui.maxv;
  gui.fwhmstart = gui.maxv;
  gui.fwhmend = gui.maxv;
  gui.centergravity = gui.maxv;

  %Start plotting
  hold(gui.handles.plotaxes,'off');
  istimeresolved = (sum(~isnan(gui.outdata(:,1)))>1);
  if istimeresolved
    plot(gui.handles.plotaxes,gui.t,gui.outdata(:,1),SET(gui.no).Roi(gui.rois(1)).LineSpec);
  else
    plot(gui.handles.plotaxes,gui.t,gui.outdata(:,1),sprintf('%so',SET(gui.no).Roi(gui.rois(1)).LineSpec(1)));
  end
  hold(gui.handles.plotaxes,'on');
  for rloop=2:length(gui.rois)    
    istimeresolved = (sum(~isnan(gui.outdata(:,rloop)))>1);
    if istimeresolved
      plot(gui.handles.plotaxes,gui.t,gui.outdata(:,rloop),SET(gui.no).Roi(gui.rois(rloop)).LineSpec);
    else
      plot(gui.handles.plotaxes,gui.t,gui.outdata(:,rloop),sprintf('%so',SET(gui.no).Roi(gui.rois(rloop)).LineSpec(1)));
    end
  end;
  %--- Calc min/max/slope  
  
  %Loop over rois
  for rloop=1:size(gui.outdata,2)
    
    %Max/Min
    temp = gui.smoothoutdata(:,rloop);
    temp(1:SET(gui.no).StartAnalysis-1) = NaN;
    temp(SET(gui.no).EndAnalysis+1:end) = NaN;
    [gui.minv(rloop),gui.minind(rloop)] = min(temp);
    [gui.maxv(rloop),gui.maxind(rloop)] = max(temp);
    
    %Slopes
    if dosmooth
      tempd = conv2(gui.smoothoutdata(:,rloop),[1;0;-1]/(2*SET(gui.no).TIncr*1000),'same');
      tempd(1:SET(gui.no).StartAnalysis-1) = NaN;
      tempd(SET(gui.no).EndAnalysis+1:end) = NaN;
      [gui.minvd(rloop),gui.minindd(rloop)] = min(tempd);
      [gui.maxvd(rloop),gui.maxindd(rloop)] = max(tempd);
    end;
    
    %FWHM
    gui.fwhm(rloop) = 0.5*(gui.maxv(rloop)+gui.minv(rloop));
    pos = find(temp>=gui.fwhm(rloop));
    if (length(pos)<1)
      %Constant value? For whatever reason, skip FWHM
      gui.fwhmstart(rloop) = NaN;
      gui.fwhmend(rloop) = NaN;
      %myfailed('FWHM not found, adjust start/end bars.',DATA.GUI.Segment);
      %return;
    else
      gui.fwhmstart(rloop) = pos(1);
      gui.fwhmend(rloop) = pos(end);
    end;
    
    %Centerofgravity
    temp = gui.smoothoutdata(SET(gui.no).StartAnalysis:SET(gui.no).EndAnalysis,rloop);
    tempt = (SET(gui.no).StartAnalysis:SET(gui.no).EndAnalysis);
    temp_nom = sum(temp'.*tempt);
    temp_den = sum(temp);
    if ~isequal(temp_den,0)
      gui.centergravity(rloop) = temp_nom/temp_den;
    else
      gui.centergravity(rloop) = NaN;
    end;
    
    %Plot max/min
    if mygetvalue(gui.handles.showminmaxcheckbox)
      plot(gui.handles.plotaxes,gui.t(gui.minind(rloop)),gui.minv(rloop),'k*');
      plot(gui.handles.plotaxes,gui.t(gui.maxind(rloop)),gui.maxv(rloop),'k*');
      plot(gui.handles.plotaxes,[gui.t(1) gui.t(end)],[gui.minv(rloop) gui.minv(rloop)],'k:');
      plot(gui.handles.plotaxes,[gui.t(1) gui.t(end)],[gui.maxv(rloop) gui.maxv(rloop)],'k:');
    end;
    
    %Plot min/max slopes
    if mygetvalue(gui.handles.showslopescheckbox)
      temp = gui.smoothoutdata(:,rloop);
      plot(gui.handles.plotaxes,gui.t(gui.maxindd(rloop)),temp(gui.maxindd(rloop)),'ko');
      plot(gui.handles.plotaxes,gui.t(gui.minindd(rloop)),temp(gui.minindd(rloop)),'ko');
      
      ylim = get(gui.handles.plotaxes,'ylim');
      
      deltat = gui.t(end)*0.03;
      left = temp(gui.minindd(rloop))-gui.minvd(rloop)*deltat;
      right = temp(gui.minindd(rloop))+gui.minvd(rloop)*deltat;
      plot(gui.handles.plotaxes,gui.t(gui.minindd(rloop))+[-deltat deltat],[left right],'k-');
      plot(gui.handles.plotaxes,[gui.t(gui.minindd(rloop)) gui.t(gui.minindd(rloop))],ylim,'k:');
      
      deltat = gui.t(end)*0.03;
      left = temp(gui.maxindd(rloop))-gui.maxvd(rloop)*deltat;
      right = temp(gui.maxindd(rloop))+gui.maxvd(rloop)*deltat;
      plot(gui.handles.plotaxes,gui.t(gui.maxindd(rloop))+[-deltat deltat],[left right],'k-');
      plot(gui.handles.plotaxes,[gui.t(gui.maxindd(rloop)) gui.t(gui.maxindd(rloop))],ylim,'k:');
    end;
    
    %Plot FWHM
    if mygetvalue(gui.handles.showfwhmcheckbox)&&~isnan(gui.fwhmstart(rloop))
      ylim = get(gui.handles.plotaxes,'ylim');
      plot(gui.handles.plotaxes,[gui.t(1) gui.t(end)],[gui.fwhm(rloop) gui.fwhm(rloop)],'k:');
      plot(gui.handles.plotaxes,[gui.t(gui.fwhmstart(rloop)) gui.t(gui.fwhmstart(rloop))],ylim,'k:');
      plot(gui.handles.plotaxes,[gui.t(gui.fwhmend(rloop)) gui.t(gui.fwhmend(rloop))],ylim,'k:');
    end;
    
    %Plot CenterGravity
    if mygetvalue(gui.handles.showcentergravitycheckbox)
      ylim = get(gui.handles.plotaxes,'ylim');
      pos = (gui.centergravity(rloop)-1)*SET(gui.no).TIncr*1000;
      plot(gui.handles.plotaxes,[pos pos],ylim,'k:');
    end;
    
  end; %rloop
  
  %Plot smoothed curves
  if dosmooth
    for rloop=1:size(gui.outdata,2)
      h = plot(gui.handles.plotaxes,gui.t,gui.smoothoutdata(:,rloop),'k:'); %SET(gui.no).RoiLineSpec{gui.rois(rloop)});
      set(h,'linewidth',2);
    end;
  end;
  
  %Add red bars
  ylim = get(gui.handles.plotaxes,'ylim');
  h = plot([gui.t(SET(gui.no).StartAnalysis) gui.t(SET(gui.no).StartAnalysis)],ylim,'r-');
  set(h,'linewidth',2,'ButtonDownFcn','reportroi(''roibar_Buttondown'',''startbar'')');
  h = plot([gui.t(SET(gui.no).EndAnalysis)   gui.t(SET(gui.no).EndAnalysis)],ylim,'r-');
  set(h,'linewidth',2,'ButtonDownFcn','reportroi(''roibar_Buttondown'',''endbar'')');
  hold(gui.handles.plotaxes,'off');
  
  %Add legend and labels
  legendstring=cell(1,length(gui.rois));
  for lloop=1:length(gui.rois)
    legendstring{lloop}=SET(gui.no).Roi(gui.rois(lloop)).Name;
  end
  legend(gui.handles.plotaxes,legendstring);
  xlabel(gui.handles.plotaxes,translation.dictionary('Time [ms]'));
  ylabel(gui.handles.plotaxes,gui.outname);
  
else %not time resolved
  
  legendstring = cell(1,length(gui.rois));
  for lloop=1:length(gui.rois)
    legendstring{lloop} = SET(gui.no).Roi(gui.rois(lloop)).Name;
  end
  bar(gui.handles.plotaxes,gui.outdata(gui.rois)');
  set(gui.handles.plotaxes,'xtick',1:1:length(legendstring));
  set(gui.handles.plotaxes,'xticklabel',legendstring);

end;

%--------------
function export %#ok<DEFNU>
%--------------
global DATA SET
gui = DATA.GUI.ROI;

doplot;
c = cell(4+size(gui.outdata,1),1+length(gui.rois)*7);
c{1,1} = SET(gui.no).PatientInfo.Name;
if dosmooth
  c{1,2} = 'Smoothing sigma [timeframes]';
  c{1,3} = gui.sigma;
end;
c{4,1} = 'Time[ms]';
for rloop=1:length(gui.rois)
  if dosmooth
    c{4,(rloop-1)*7+2} = sprintf('Smoothed:%s',gui.outname);
  else
    c{4,(rloop-1)*7+2} = sprintf('%s',gui.outname);
  end;
  c{4,(rloop-1)*7+3} = 'Intensity';
  c{4,(rloop-1)*7+4} = 'SD';
  c{4,(rloop-1)*7+5} = 'Size[cm2]';
  c{4,(rloop-1)*7+6} = 'Pixsize[cm2]';
  c{4,(rloop-1)*7+7} = 'Max';
  c{4,(rloop-1)*7+8} = 'Min';
end;

for rloop=1:length(gui.rois)
  roiname{rloop} = SET(gui.no).Roi(gui.rois(rloop)).Name;
  isroinumber(rloop) = (sum(ismember(['ROI-'],roiname{rloop}))==4); %check if ROI names are ROI-1, ROI-2, etc
end
if sum(isroinumber)==length(isroinumber) %sort based on numbers after "ROI-"
  for rloop=1:length(gui.rois)
    thisroiname = roiname{rloop};
    roinumber(rloop) = str2num(thisroiname(5:end));
  end
  [~,exportindex] = sort(roinumber);  
else %sort in alphabetic order
  [~,exportindex] = sort(roiname);
end
for tloop=1:SET(gui.no).TSize
  c{tloop+4,1} = gui.t(tloop);
  for rloop=1:length(gui.rois)%Print values
    c{3,(rloop-1)*7+2} = roiname{exportindex(rloop)};
    c{tloop+4,(rloop-1)*7+2} = gui.smoothoutdata(tloop,exportindex(rloop));
    c{tloop+4,(rloop-1)*7+3} = gui.roiint(tloop,exportindex(rloop));
    c{tloop+4,(rloop-1)*7+4} = gui.roistd(tloop,exportindex(rloop));
    c{tloop+4,(rloop-1)*7+5} = gui.roisize(tloop,exportindex(rloop));
    c{tloop+4,(rloop-1)*7+6} = gui.roisizepix(tloop,exportindex(rloop));
    c{tloop+4,(rloop-1)*7+7} = gui.roimax(tloop,exportindex(rloop));
    c{tloop+4,(rloop-1)*7+8} = gui.roimin(tloop,exportindex(rloop));
  end;
end;

segment('cell2clipboard',c);

%--------------------
function minmaxexport
%--------------------
global DATA SET
gui = DATA.GUI.ROI;

if SET(gui.no).TSize <= 1
  myfailed('The image is not time resolved.',DATA.GUI.Segment);
  return;
end

if any(gui.maxindd == 0)
  myfailed('Need to smooth data first.');
  return;
end

doplot;
c = cell(12,1+2*length(gui.rois));
c{1,1} = SET(gui.no).PatientInfo.Name;
c{4,1} = 'Type';
c{5,1} = 'StartAnalysis';
c{6,1} = 'EndAnalysis';
c{7,1} = 'Offset';
c{8,1} = 'MaxValue';
c{9,1} = 'MinValue';
c{10,1} = 'MaxSlope[unit/ms]';
c{11,1} = 'MinSlope[unit/ms]';
c{12,1} = 'FWHMstart';
c{13,1} = 'FWHMend';
c{14,1} = 'CenterGravity';
c{5,3} = 1000*(SET(gui.no).StartAnalysis-1)*SET(gui.no).TIncr;
c{6,3} = 1000*(SET(gui.no).EndAnalysis-1)*SET(gui.no).TIncr;
if dosmooth
  c{16,1} = 'Sigma:';
  c{16,2} = gui.sigma;
else
  c{16,1} = 'No smoothing';
end;
if get(gui.handles.cropzerocheckbox,'value')
  c{17,1} = 'Cropping';
else
  c{17,1} = 'No cropping';
end;
for rloop=1:length(gui.rois)
  c{3,2+(rloop-1)*2} = SET(gui.no).Roi(gui.rois(rloop)).Name;
  c{4,2+(rloop-1)*2} = 'Value';
  c{4,3+(rloop-1)*2} = 'Time[ms]';
  
  c{7,2+(rloop-1)*2} = gui.offset;
  
  c{8,2+(rloop-1)*2} = gui.maxv(rloop);
  c{8,3+(rloop-1)*2} = gui.t(gui.maxind(rloop));
  
  c{9,2+(rloop-1)*2} = gui.minv(rloop);
  c{9,3+(rloop-1)*2} = gui.t(gui.minind(rloop));
  
  c{10,2+(rloop-1)*2} = gui.maxvd(rloop);
  c{10,3+(rloop-1)*2} = gui.t(gui.maxindd(rloop));
  
  c{11,2+(rloop-1)*2} = gui.minvd(rloop);
  c{11,3+(rloop-1)*2} = gui.t(gui.minindd(rloop));
  
  c{12,2+(rloop-1)*2} = gui.fwhm(rloop);
  c{12,3+(rloop-1)*2} = gui.t(gui.fwhmstart(rloop));
  
  c{13,2+(rloop-1)*2} = gui.fwhm(rloop);
  c{13,3+(rloop-1)*2} = gui.t(gui.fwhmend(rloop));
  
  c{14,3+(rloop-1)*2} = (gui.centergravity(rloop)-1)*1000*SET(gui.no).TIncr;
end;
segment('cell2clipboard',c);
  
%----------------------
function close_Callback
%-----------------------
global DATA
datacursormode off;
try
  DATA.GUI.ROI=close(DATA.GUI.ROI);  %close the flow gui
catch %#ok<CTCH>
  close(gcf)
end

DATA.GUI.ROI = [];

%---------------------------
function smoothen = dosmooth
%---------------------------
global DATA
gui = DATA.GUI.ROI;

%--- checkboxes
if mygetvalue(gui.handles.smoothcheckbox)
  smoothen = true;
else
  smoothen = false;
end;

%---------------------------
function datacursor_Callback %#ok<DEFNU>
%---------------------------
%Function which enables datacursor probing on roi-analysis plots.
%/SB
global DATA
gui = DATA.GUI.ROI;

datacursormode(gui.fig);

%-----------------------------------------
function [xf] = smooth_helper(sigma,t,x)
%-----------------------------------------
%Helper function to interpolate non regular signals

%Interpolate non-rectilinear to rectilinear, use twice as many points as in
%original signal
ti = linspace(t(1),t(end),length(t)*2);
xi = interp1(t,x,ti,'linear');

%create filter
fx = linspace(-60,60,189); %was 121, approximate how much more points needed to get same smoothing as for rectilinear.
f = exp(-fx.^2/sigma.^2); %double time resolution
  
%Filter using normalized averaging same algorithm as for rectilinear
%signals.
xif = conv2(xi,f,'same')./(eps+conv2(ones(size(xi)),f,'same'));

%Resample back to non rectilinear signal
xf = interp1(ti,xif,t);