function [varargout] = transittimetool(varargin)
%GUI for transit time analysis.
%Nils Lundahl

macro_helper(varargin{:});
if nargin == 0 || isempty(varargin{1})
  varargin{1} = 'init';
end
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init
%------------
global DATA SET NO

if SET(NO).RoiN<2
  myfailed('Need 2 ROIs for transit time calculations.',DATA.GUI.Segment);
  return;
end

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

DATA.investigationaluselabel; %In Segment CMR: Warning Only for investigational use

if isempty(SET(NO).EndAnalysis)
  SET(NO).EndAnalysis = SET(NO).TSize;
end

if isempty(SET(NO).StartAnalysis)
  SET(NO).StartAnalysis = 1;
end

%if (SET(NO).EndAnalysis-SET(NO).StartAnalysis)<1
%  myfailed('No timespan, adjust Start/End analysis time (under stress menu or drag red bar).',DATA.GUI.Segment);
%  return;
%
%end;

%Show gui
gui = mygui('transittimetool.fig');
DATA.GUI.TransitTimeTool = gui;

% if ~isrectilinear(SET(NO).TimeVector)
%   mywarning(['Non uniform time steps detected. Display of '...
%     'min/max slope disabled.']);
% end

%Initalize variables
gui.no = NO;
gui.t = SET(NO).TimeVector;
gui.roiint = nan([SET(NO).TSize 2]);
gui.roistd = gui.roiint;
gui.roisizepix = gui.roiint;
gui.roisize = gui.roiint;
gui.roimin = gui.roiint;
gui.roimax = gui.roiint;
gui.offset = 0;
gui.detrendstartms = 0;
gui.detrendendvalue = 0;
gui.baseline = [0 0];

gui.timeframes = 1:SET(gui.no).TSize;
gui.normalized = 1;
gui.limits = [SET(gui.no).StartAnalysis SET(gui.no).EndAnalysis; ...
  SET(gui.no).StartAnalysis SET(gui.no).EndAnalysis];

recalc; %do plot

%set pointer, eriks
set(gui.fig,'Pointer','arrow');


%------------------
function plotfilter %#ok<DEFNU>
%------------------
%get sigma from slider
global DATA
gui = DATA.GUI.TransitTimeTool;

gui.sigma = mygetvalue(gui.handles.smoothslider);

%create filter
x = linspace(-60,60,121);
f = exp(-x.^2/gui.sigma.^2);
figure(22);
plot(x,f);
xlabel(dprintf('Timeframes'));
title(dprintf('Smoothing applicability.'));

%--------------
function recalc
%--------------
global DATA SET
gui = DATA.GUI.TransitTimeTool;

%--- Ask for what ROI's to include
roi1 = mymenu('Onset ROI',{SET(gui.no).Roi.Name});
roi2inds = setdiff(1:SET(gui.no).RoiN,roi1);
if numel(roi2inds) > 1
  roi2 = mymenu('Offset ROI',{SET(gui.no).Roi(roi2inds).Name});
else
  roi2 = 1;
end
gui.rois = [roi1 roi2inds(roi2)];
gui.timeframes = 1:SET(gui.no).TSize;
gui.normalized = 1;

%--- Extract data
h = mywaitbarstart(length(gui.rois)*length(gui.timeframes),'Please wait, calculating data',[],DATA.GUI.Segment);
for rloop=1:length(gui.rois)
  z = SET(gui.no).Roi(gui.rois(rloop)).Z;
  for tloop=gui.timeframes
    if not(isnan(SET(gui.no).Roi(gui.rois(rloop)).X(1,tloop)))
      if gui.normalized
        temp = SET(gui.no).IM(:,:,tloop,z);
      else
        temp = calcfunctions('calctruedata',SET(gui.no).IM(:,:,tloop,z),gui.no);
      end
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
      end
      gui.roisize(tloop,rloop) = (1/100)*stablepolyarea(...
        SET(gui.no).ResolutionY*SET(gui.no).Roi(gui.rois(rloop)).Y(:,tloop),...
        SET(gui.no).ResolutionX*SET(gui.no).Roi(gui.rois(rloop)).X(:,tloop));
      gui.roisizepix(tloop,rloop) = sum(roimask(:))*SET(gui.no).ResolutionX*SET(gui.no).ResolutionY/100;
      h = mywaitbarupdate(h);
    end
  end
end
mywaitbarclose(h);
doplot;
set(gui.fig,'Pointer','Arrow');

%--------------
function doplot
%--------------
global DATA SET
gui = DATA.GUI.TransitTimeTool;

%--- Check what to plot
parameter = 1;
switch parameter
  case 1
    gui.outdata = gui.roiint;
    gui.outname = 'Mean intensity';
  case 2
    gui.outdata = gui.roistd;
    gui.outname = 'Standard deviation';
  case 3
    gui.outdata = gui.roisize;
    gui.outname = 'Area [cm^2]';
  case 4
    gui.outdata = gui.roisizepix;
    gui.outname = 'Area [cm^2] based on pixels';
  case 5
    gui.outdata = gui.roimax;
    gui.outname = 'Maximum intensity.';
  case 6
    gui.outdata = gui.roimin;
    gui.outname = 'Minimum intensity.';
  otherwise
    myfailed('Not yet implemented.',DATA.GUI.Segment);
    gui.outdata = [];
    gui.outname = '';
end

bl1 = str2double(get(gui.handles.inbaseedit,'String'));
if ~isnan(bl1)
  gui.baseline(1) = bl1;
else
  set(gui.handles.inbaseedit,'String',sprintf('%0.3f',gui.baseline(1)));
end
bl2 = str2double(get(gui.handles.outbaseedit,'String'));
if ~isnan(bl2)
  gui.baseline(2) = bl2;
else
  set(gui.handles.outbaseedit,'String',sprintf('%0.3f',gui.baseline(2)))
end

%Add offset
gui.outdata = gui.outdata+gui.offset;

%Detrend data
if false && get(gui.handles.detrendcheckbox,'value')
  t = SET(gui.no).TimeVector;
  t(t < gui.detrendstartms) = NaN;
  t = t - min(t);
  t = t./max(t); %Normalize
  t(isnan(t)) = 0;
  detrend = t*gui.detrendendvalue;
  detrend = detrend(:); %reshape
  detrend = repmat(detrend,[1 size(gui.outdata,2)]); %reshape
  gui.outdata = gui.outdata-detrend;
end


%--- Do smoothing
dosmooth = false;
if dosmooth
  %get sigma from slider
  gui.sigma = mygetvalue(gui.handles.smoothslider);
  
  %create filter
  x = linspace(-60,60,121);
  f = exp(-x.^2/gui.sigma.^2);
  
  gui.smoothoutdata = gui.outdata;
  for rloop=1:size(gui.outdata,2)
    temp = gui.outdata(:,rloop);
    temp = conv2(temp,f','same')./(eps+conv2(ones(size(temp)),f','same'));
    gui.smoothoutdata(:,rloop) = temp;
  end
else
  gui.smoothoutdata = gui.outdata;
  gui.t = SET(gui.no).TimeVector;
end

basecheckbox = [gui.handles.inbasecheckbox gui.handles.outbasecheckbox];

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
	hold(gui.handles.inflowaxes,'off');
  plot(gui.handles.inflowaxes,gui.t,gui.outdata(:,1),SET(gui.no).Roi(gui.rois(1)).LineSpec);
  hold(gui.handles.inflowaxes,'on');
  plot(gui.handles.outflowaxes,gui.t,gui.outdata(:,2),SET(gui.no).Roi(gui.rois(2)).LineSpec);
  hold(gui.handles.outflowaxes,'on');
  
  %--- Calc min/max/slope  
  
  %Loop over rois
  handles = [gui.handles.inflowaxes gui.handles.outflowaxes];
  centeredithandles = [gui.handles.oncenteredit gui.handles.offcenteredit];
  for rloop=1:size(gui.outdata,2)
    
    %Max/Min
    temp = gui.smoothoutdata(:,rloop);
    temp(1:gui.limits(rloop,1)-1) = NaN;
    temp(gui.limits(rloop,2)+1:end) = NaN;
    [gui.minv(rloop),gui.minind(rloop)] = min(temp);
    [gui.maxv(rloop),gui.maxind(rloop)] = max(temp);
    
    %Slopes
    tempd = conv2(gui.smoothoutdata(:,rloop),[1;0;-1]/(2*SET(gui.no).TIncr),'same');
    tempd(1:gui.limits(rloop,1)-1) = NaN;
    tempd(gui.limits(rloop,2)+1:end) = NaN;
    [gui.minvd(rloop),gui.minindd(rloop)] = min(tempd);
    [gui.maxvd(rloop),gui.maxindd(rloop)] = max(tempd);
    
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
    end
    
    %Centerofgravity
    temp = gui.smoothoutdata(gui.limits(rloop,1):gui.limits(rloop,2),rloop);
    
    if mygetvalue(basecheckbox(rloop))
      temp = temp-gui.baseline(rloop);
      temp(temp < 0) = 0;
    end
    tempt = (gui.limits(rloop,1):gui.limits(rloop,2));
    temp_nom = sum(temp'.*tempt);
    temp_den = sum(temp);
    if ~isequal(temp_den,0)
      gui.centergravity(rloop) = temp_nom/temp_den;
    else
      gui.centergravity(rloop) = NaN;
    end
    
    %Plot max/min
    if get(gui.handles.showminmaxcheckbox,'value')
      plot(handles(rloop),gui.t(gui.minind(rloop)),gui.minv(rloop),'k*');
      plot(handles(rloop),gui.t(gui.maxind(rloop)),gui.maxv(rloop),'k*');
      plot(handles(rloop),[gui.t(1) gui.t(end)],[gui.minv(rloop) gui.minv(rloop)],'k:');
      plot(handles(rloop),[gui.t(1) gui.t(end)],[gui.maxv(rloop) gui.maxv(rloop)],'k:');
    end
    
    %Plot min/max slopes
    if get(gui.handles.showslopescheckbox,'value')
      temp = gui.smoothoutdata(:,rloop);
      plot(handles(rloop),gui.t(gui.maxindd(rloop)),temp(gui.maxindd(rloop)),'ko');
      plot(handles(rloop),gui.t(gui.minindd(rloop)),temp(gui.minindd(rloop)),'ko');
      
      ylim = get(handles(rloop),'ylim');
      
      deltat = gui.t(end)*0.03;
      left = temp(gui.minindd(rloop))-gui.minvd(rloop)*deltat;
      right = temp(gui.minindd(rloop))+gui.minvd(rloop)*deltat;
      plot(handles(rloop),gui.t(gui.minindd(rloop))+[-deltat deltat],[left right],'k-');
      plot(handles(rloop),[gui.t(gui.minindd(rloop)) gui.t(gui.minindd(rloop))],ylim,'k:');
      
      deltat = gui.t(end)*0.03;
      left = temp(gui.maxindd(rloop))-gui.maxvd(rloop)*deltat;
      right = temp(gui.maxindd(rloop))+gui.maxvd(rloop)*deltat;
      plot(handles(rloop),gui.t(gui.maxindd(rloop))+[-deltat deltat],[left right],'k-');
      plot(handles(rloop),[gui.t(gui.maxindd(rloop)) gui.t(gui.maxindd(rloop))],ylim,'k:');
    end
    
    %Plot FWHM
    if get(gui.handles.showfwhmcheckbox,'value')&&~isnan(gui.fwhmstart(rloop))
      ylim = get(handles(rloop),'ylim');
      plot(handles(rloop),[gui.t(1) gui.t(end)],[gui.fwhm(rloop) gui.fwhm(rloop)],'k:');
      plot(handles(rloop),[gui.t(gui.fwhmstart(rloop)) gui.t(gui.fwhmstart(rloop))],ylim,'k:');
      plot(handles(rloop),[gui.t(gui.fwhmend(rloop)) gui.t(gui.fwhmend(rloop))],ylim,'k:');
    end
    
    %Plot CenterGravity
    
    centergravityt = (gui.centergravity(rloop)-1)*SET(gui.no).TIncr;
    set(centeredithandles(rloop),'String',sprintf('%0.2f',centergravityt));
    if get(gui.handles.showcentergravitycheckbox,'value')
      ylim = get(handles(rloop),'ylim');
      plot(handles(rloop),centergravityt*[1 1],ylim,'k:');
    end
    
  end %rloop
  
  onradiobuttons = [...
    gui.handles.onminradiobutton ...
    gui.handles.onmaxupradiobutton ...
    gui.handles.onupradiobutton ...
    gui.handles.onmaxradiobutton ...
    gui.handles.ondownradiobutton ...
    gui.handles.oncenterradiobutton];
  onval = find(cell2mat(get(onradiobuttons,'value')));
  switch onval
    case 6
      set(gui.handles.transittimeedit,'String',sprintf('%0.3f',...
        diff(gui.centergravity)*SET(gui.no).TIncr));
  end
  
  %Plot smoothed curves
  if dosmooth
    plot(gui.handles.inflowaxes,gui.t,gui.smoothoutdata(:,1),'k:');
    plot(gui.handles.outflowaxes,gui.t,gui.smoothoutdata(:,2),'k:');
%     for rloop=1:size(gui.outdata,2)
%       h = plot(handles(rloop),gui.t,gui.smoothoutdata(:,rloop),'k:'); %SET(gui.no).RoiLineSpec{gui.rois(rloop)});
%       set(h,'linewidth',2);
%       
%     end;
  end
    
  %Add baseline crop
  if mygetvalue(gui.handles.inbasecheckbox)
    xlim = get(gui.handles.inflowaxes,'xlim');
    h = plot(gui.handles.inflowaxes,xlim,gui.baseline(1)*[1 1],'r-');
    set(h,'linewidth',2,'ButtonDownFcn', ...
      'transittimetool(''baseline_Buttondown'',''inflow'')');
    gui.handles.inbaseline = h;
  end
    
  if mygetvalue(gui.handles.outbasecheckbox)
  xlim = get(gui.handles.outflowaxes,'xlim');
  h = plot(gui.handles.outflowaxes,xlim,gui.baseline(2)*[1 1],'r-');
  set(h,'linewidth',2,'ButtonDownFcn', ...
    'transittimetool(''baseline_Buttondown'',''outflow'')');
  gui.handles.outbaseline = h;
  end
    
  %Add red time bars
  ylim = get(gui.handles.inflowaxes,'ylim');
  h = plot(gui.handles.inflowaxes,[gui.t(gui.limits(1,1)) gui.t(gui.limits(1,1))],ylim,'r-');
  set(h,'linewidth',2,'ButtonDownFcn','transittimetool(''roibar_Buttondown'',''startbar'',''onset'')');
  h = plot(gui.handles.inflowaxes,[gui.t(gui.limits(1,2))   gui.t(gui.limits(1,2))],ylim,'r-');
  set(h,'linewidth',2,'ButtonDownFcn','transittimetool(''roibar_Buttondown'',''endbar'',''onset'')');
  hold(gui.handles.inflowaxes,'off');
  
  ylim = get(gui.handles.outflowaxes,'ylim');
  h = plot(gui.handles.outflowaxes,[gui.t(gui.limits(2,1)) gui.t(gui.limits(2,1))],ylim,'r-');
  set(h,'linewidth',2,'ButtonDownFcn','transittimetool(''roibar_Buttondown'',''startbar'',''offset'')');
  h = plot(gui.handles.outflowaxes,[gui.t(gui.limits(2,2))   gui.t(gui.limits(2,2))],ylim,'r-');
  set(h,'linewidth',2,'ButtonDownFcn','transittimetool(''roibar_Buttondown'',''endbar'',''offset'')');
  hold(gui.handles.outflowaxes,'off');
  
  %Add legend and labels
  legendstring=cell(1,length(gui.rois));
  for lloop=1:length(gui.rois)
   legendstring{lloop}=SET(gui.no).Roi(gui.rois(lloop)).Name;
  end
  legend(gui.handles.inflowaxes,legendstring{1});
  legend(gui.handles.outflowaxes,legendstring{2});
  
  xlabel(gui.handles.inflowaxes,dprintf('Time [s]'));
  ylabel(gui.handles.inflowaxes,dprintf('Signal Intensity [a.u.]'));
  xlabel(gui.handles.outflowaxes,dprintf('Time [s]'));
  ylabel(gui.handles.outflowaxes,dprintf('Signal Intensity [a.u.]'));
  
  yticks = get(gui.handles.inflowaxes,'YTick');
  yticksspaced = yticks(1):0.01:yticks(end);
  if numel(yticks) > numel(yticksspaced)
    set(gui.handles.inflowaxes,'YTick',yticksspaced);
  end
  
  yticks = get(gui.handles.outflowaxes,'YTick');
  yticksspaced = yticks(1):0.01:yticks(end);
  if numel(yticks) > numel(yticksspaced)
    set(gui.handles.outflowaxes,'YTick',yticksspaced);
  end
  
else %not time resolved
  
  %legendstring = cell(1,length(gui.rois));
  %for lloop=1:length(gui.rois)
  %  legendstring{lloop} = SET(gui.no).Roi(gui.rois(lloop)).Name;
  %end
  bar(gui.handles.inflowaxes,gui.outdata(gui.rois)');
  %set(gui.handles.inflowaxes,'xticklabel',legendstring);

end

%------------------------------------
function baseline_Buttondown(inorout)
%------------------------------------
global DATA
gui = DATA.GUI.TransitTimeTool;
set(gui.fig,'WindowButtonMotionFcn',...
  sprintf('transittimetool(''baseline_Motion'',''%s'')',inorout));
set(gui.fig,'WindowButtonUpFcn',...
  sprintf('transittimetool(''baseline_Buttonup'',''%s'')',inorout));

%------------------------------------
function baseline_Motion(inorout)
%------------------------------------
global DATA
gui = DATA.GUI.TransitTimeTool;
switch inorout
  case 'inflow'
    [x,y] = mygetcurrentpoint(gui.handles.inflowaxes);
    set(gui.handles.inbaseline,'YData',[y y]);
    set(gui.handles.inbaseedit,'String',sprintf('%0.3f',y));
  case 'outflow'
    [x,y] = mygetcurrentpoint(gui.handles.outflowaxes);
    set(gui.handles.outbaseline,'YData',[y y]);
    set(gui.handles.outbaseedit,'String',sprintf('%0.3f',y));   
end

%----------------------------------
function baseline_Buttonup(inorout)
%----------------------------------
global DATA
gui = DATA.GUI.TransitTimeTool;
set(gui.fig,'WindowButtonMotionFcn','');
set(gui.fig,'WindowButtonUpFcn','');


%---------------------------------------
function roibar_Buttondown(type,inorout) %#ok<DEFNU>
%---------------------------------------
%Button down function for pressing timebar in transit time tool gui. 
%Sets button motion and button up function for changing start time and end 
%time for calculation.
global SET NO

SET(NO).Stress.temphandle = gcbo;
set(gcf,'WindowButtonMotionFcn','transittimetool(''roibar_Motion'')');
set(gcf,'WindowButtonUpFcn',sprintf(...
  'transittimetool(''roibar_Buttonup'',''%s'',''%s'')',type,inorout));


%---------------------
function roibar_Motion %#ok<DEFNU>
%---------------------
%Button motion function called when moving timebar in plot flow gui. 
%Changes start time and end time for calculation of flow.
global SET NO

[x,y] = mygetcurrentpoint(gca);

%Convert to timeframes and round
x = round(1+x/SET(NO).TIncr);
x = min(max(1,x),SET(NO).TSize);

%Convert back to ms
x = (x-1)*SET(NO).TIncr;

%Update display
set(SET(NO).Stress.temphandle,'xdata',[x x]);

%-------------------------------------
function roibar_Buttonup(type,inorout) %#ok<DEFNU>
%-------------------------------------
%Button up function called after moving timebar in transit time tool gui. 
%Changes start time and end time for calculation.
global DATA SET NO
gui = DATA.GUI.TransitTimeTool;

set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');

%Get Convert to timeframe
x = get(SET(NO).Stress.temphandle,'xdata');
x = x(1);
x = round(1+x/(SET(NO).TIncr));

switch [inorout '-' type]
  case 'onset-startbar'
    gui.limits(1,1) = x;
  case 'onset-endbar'
    gui.limits(1,2) = x;   
  case 'offset-startbar'
    gui.limits(2,1) = x;
  case 'offset-endbar'
    gui.limits(2,2) = x;    
end

for rloop = 1:2
  if gui.limits(rloop,1) > gui.limits(rloop,2)
    temp = gui.limits(rloop,1);
    gui.limits(rloop,1)=gui.limits(rloop,2);
    gui.limits(rloop,2)=temp;
  end
end

%------------------
function smoothedit
%------------------
global DATA
gui = DATA.GUI.TransitTimeTool;

s = mygetedit(gui.handles.sigmaedit);
n = str2double(s);
if not(isnan(n))
  gui.sigma = n;
end
gui.sigma = min(max(gui.sigma,get(gui.handles.smoothslider,'min')),get(gui.handles.smoothslider,'max'));
set(gui.handles.sigmaedit,'string',sprintf('%0.5g',gui.sigma));
set(gui.handles.smoothslider,'value',gui.sigma);

%--------------------
function smoothslider
%--------------------
global DATA
gui = DATA.GUI.TransitTimeTool;

gui.sigma = mygetvalue(gui.handles.smoothslider);
set(gui.handles.sigmaedit,'string',sprintf('%0.5g',gui.sigma));

%--------------
function export
%--------------
global DATA
gui = DATA.GUI.TransitTimeTool;
c = {get(gui.handles.transittimeedit,'String')};

segment('cell2clipboard',c);


%--------------------
function minmaxexport
%--------------------
global DATA SET
gui = DATA.GUI.TransitTimeTool;

if SET(gui.no).TSize <= 1
  myfailed('The image is not time resolved.',DATA.GUI.Segment);
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
c{5,3} = (SET(gui.no).StartAnalysis-1)*SET(gui.no).TIncr;
c{6,3} = (SET(gui.no).EndAnalysis-1)*SET(gui.no).TIncr;
if dosmooth
  c{16,1} = 'Sigma:';
  c{16,2} = gui.sigma;
else
  c{16,1} = 'No smoothing';
end
if get(gui.handles.cropzerocheckbox,'value')
  c{17,1} = 'Cropping';
else
  c{17,1} = 'No cropping';
end
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
  
  c{14,3+(rloop-1)*2} = (gui.centergravity(rloop)-1)*SET(gui.no).TIncr;
end
segment('cell2clipboard',c);

%----------------------
function close_Callback %#ok<DEFNU>
%-----------------------
%Close transit time tool GUI
global DATA

try
  DATA.GUI.TransitTimeTool=close(DATA.GUI.TransitTimeTool);  %close the gui
catch %#ok<CTCH>
  delete(gcbf)
end

DATA.GUI.TransitTimeTool = [];