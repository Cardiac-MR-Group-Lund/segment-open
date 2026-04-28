function [varargout] = reportflow(varargin)
%GUI for 2D flow analysis

%Adapted for use with mygui class by Nils Lundahl

%#ok<*GVMIS>

if nargin == 0 || isempty(varargin{1})
  varargin{1} = 'init';
elseif strncmp(varargin{1},'get',3)
  [varargout{1:nargout}] = getdata(varargin{:});
  return
end
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%---------------------------
function resampedit_Callback
%---------------------------
global DATA  

gui = DATA.GUI.Flow;

[temp,ok] = str2num(mygetedit(gui.handles.resampedit));
if not(ok)
	disp('Illegal value');
else
	gui.resampsize = temp;
end

%----------------------------
function resampreset_Callback 
%----------------------------
global NO

%Recalculate to get correct flow
init(NO);

%-------------------------
function resamp_Callback 
%-------------------------
global DATA SET NO
no = NO;
gui = DATA.GUI.Flow;

gui.resamped = 0;
init(no); %get real values from recalc everytime!
numrois = gui.numrois;
t = gui.t;
size(t)
if ~isempty(gui.resampsize)
	TSize2 = gui.resampsize;
else
	disp('TSize2 error');
	TSize2 = 5;
end

%create new data structures
newnetflow = NaN(TSize2, numrois);
newvelmean = newnetflow;
newvelstd = newnetflow;
newvelmax = newnetflow;
newvelmin = newnetflow;
newkenergy = newnetflow;
newarea = newnetflow;
newposflow = newnetflow;
newnegflow = newnetflow;

for loop=1:numrois
	[newt, newy] = calcfunctions('resamplemodel',t',gui.netflow(:,loop),TSize2);
	newnetflow(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.velmean(:,loop),TSize2);
	newvelmean(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.velstd(:,loop),TSize2);
	newvelstd(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.velmax(:,loop),TSize2);
	newvelmax(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.velmin(:,loop),TSize2);
	newvelmin(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.kenergy(:,loop),TSize2);
	newkenergy(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.area(:,loop),TSize2);
	newarea(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.posflow(:,loop),TSize2);
	newposflow(:,loop) = newy;
	
	[~, newy] = calcfunctions('resamplemodel',t',gui.negflow(:,loop),TSize2);
	newnegflow(:,loop) = newy;
end

%assignment time:
gui.t = newt';
gui.velmean = newvelmean;
gui.velstd = newvelstd;
gui.velmax = newvelmax;
gui.velmin = newvelmin;
gui.kenergy = newkenergy;
gui.area = newarea;
gui.netflow = newnetflow;
gui.posflow = newposflow;
gui.negflow = newnegflow;
%Calculate diameter equivalent
gui.diameter = 2*10*sqrt(gui.area/pi); %mm*10 => cm

nom = SET(no).Flow.MagnitudeNo;

%backup variables:
gui.TIncr = ((SET(nom).TSize-1)*SET(nom).TIncr)/(TSize2-1);
gui.tsize = TSize2;

gui.resamped = 1;
recalculate(no); %Also do an update
% gui.resamped = 0;

segment('updateflow'); %DATA.updateaxestables('flow',nom,nop); %maybe unnecessary?

%---------------------------------
function ok = init(no,eddycheck, invisible)
%---------------------------------
%Init flow report GUI

global DATA SET NO

ok = false; %inform if unable to launch GUI
if nargin<1
  no = NO;
end

if nargin<2
  eddycheck = true;
end
if nargin<3
  invisible = false;
end

if isempty(SET(no).Flow)
  myfailed('No flow data in current image stack.',DATA.GUI.Segment);
  return;
end

if ~isempty(SET(no).Flow)
  nom = SET(no).Flow.MagnitudeNo;
else
  return;
end

if isempty(SET(nom).Flow.PhaseNo)
  myfailed('No through plane velocity data found.',DATA.GUI.Segment);
  return;
end

if SET(nom).RoiN==0
  myfailed('No ROIs available.',DATA.GUI.Segment);
  return;
end

tempnos=[SET(no).Flow.MagnitudeNo SET(no).Flow.PhaseNo];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  disp('non single image')
  return;
end

if ~DATA.Silent
  viewfunctions('switchimagestack',nom)
end

if isempty(SET(nom).EndAnalysis)
  SET(nom).EndAnalysis = SET(nom).TSize;
end

if isempty(SET(nom).StartAnalysis)
  SET(nom).StartAnalysis = 1;
end

if (SET(nom).EndAnalysis-SET(nom).StartAnalysis) ==0
  if ~isopengui('flowonetimephase.fig')
    if invisible
      gui = mygui('flowonetimephase.fig','invisible');
    else
      gui = mygui('flowonetimephase.fig');
    end
    DATA.GUI.Flow = gui;
  else
    gui = DATA.GUI.Flow;
  end

  %Calculate number of ROI's
  rois2take = [];

  for loop=1:SET(nom).RoiN
    %   if sum(~isnan(SET(nom).Roi(loop).X(1,:)))>1
    %More than one timeframe
    rois2take = [rois2take loop]; %#ok<AGROW>
    %   end;
  end
  gui.rois2take = rois2take;
  gui.numrois = length(rois2take);
  gui.nom = nom;
  %gui.nop = nop;

  showflowroidata(no,0);
elseif (SET(nom).EndAnalysis-SET(nom).StartAnalysis)<0 %1
  myfailed('No timespan, adjust Start/End analysis time (under stress menu or drag red bar).',DATA.GUI.Segment);
  return;
else
  nop = SET(nom).Flow.PhaseNo;
  if not(isrectilinear(SET(nom).TimeVector))
    switch mymenu(['Non uniform time steps detected. This is not '...
        'supported for flow calculations. Will use mean time step.'],...
        'Accept for now',...
        'Accept and apply to stack',...
        'Cancel')
      case 2
        SET(nom).TimeVector = SET(nom).TIncr*(0:SET(nom).TSize-1);
        SET(nop).TimeVector = SET(nom).TimeVector;
      case 3
        return
    end
  end

  %Automatically perform Eddy current. Added Einar Heiberg
  if eddycheck
    if isfield(SET(nom).Flow,'PhaseCorrAsk') && SET(nom).Flow.PhaseCorrAsk == false
      %do not ask if perform eddy current
    else
      if ~isequal(SET(nom).Scanner,'Philips')
        %Only perform for other scanners than Philips
        try
          if isempty(SET(nop).Flow.PhaseCorr)
            plotonok = true;
            floweddycurrent('initsmall',plotonok); %When specifying plotonok, then this function is called again when user press ok.
            return;
          end
        catch %Do point of crasch if we can't perform
        end
      end
    end
  end

  %Calculate number of ROI's
  rois2take = [];
  for loop=1:SET(nom).RoiN
    if sum(~isnan(SET(nom).Roi(loop).X(1,:)))>1
      %More than one timeframe
      rois2take = [rois2take loop]; %#ok<AGROW>
    else
      if matches(SET(nom).Roi(loop).Name,'Static tissue')
        % ensure time resolved static tissue
        roi('copyrpoi_helper',nom,loop);
        rois2take = [rois2take loop]; %#ok<AGROW>
      end
    end
  end
  if isempty(rois2take)
    return
  end

  %Open GUI
  if ~isopengui('flow.fig')
    if invisible
      gui = mygui('flow.fig','invisible');
    else
      gui = mygui('flow.fig');
    end
    DATA.GUI.Flow = gui;
  else
    gui = DATA.GUI.Flow;
  end
  set(gui.fig,'Name',dprintf('Flow report'));
  if ~invisible

    if ~isfield(gui.handles,'settingsmenu')
      % add menu if there is no menu
      gui.handles.settingsmenu = uimenu(...
        'Label', dprintf('Settings'), ...
        'Tag','settingsmenu');
      g.Handles.manualheartbeatsmenu = uimenu(gui.handles.settingsmenu,...
        'Label',dprintf('Manually Define Number of Heart Beats (for real time flow)'), ...
        'Tag','manualheartbeatsmenu',...
        'Callback','reportflow(''defineheartrateingui'')');
    end
  end

  gui.outsize = [SET(nom).XSize SET(nom).YSize];
  if not(isfield(SET(nom).Flow,'HeartRate')) || isempty(SET(nom).Flow.HeartRate)
    temphr = (60/(SET(nom).TSize*SET(nom).TIncr));
    SET(nom).Flow.HeartRate = temphr;
    segment('updateflow');  %update result panel
  end

  gui.hr = SET(nom).Flow.HeartRate;
  txtstr = sprintf('%s: %0.3g bpm\n%s: %d cm/s',dprintf('Heart rate'),gui.hr,'VENC',SET(nom).VENC);
  set(gui.handles.heartratetext,'String',txtstr);

  if not(isfield(SET(nop).Flow,'PhaseCorr'))
    SET(nop).Flow.PhaseCorr = [];
  end

  gui.t = SET(nom).TIncr*(-1+(1:SET(nom).TSize));
  gui.d = '.';
  gui.dt = 9;

  gui.nom = nom;
  gui.nop = nop;

  gui.rois2take = rois2take;
  gui.numrois = length(rois2take);

  %Add instance regarding flow
  gui.show = true(1,gui.numrois);

  %new variables for internal sampling of data: only for use in resampling
  gui.TIncr = [];
  gui.tsize = [];
  gui.resamped = 0;
  resampedit_Callback;
  recalculate(no);
  ok = true;

  set(gui.fig,'pointer','arrow');

  if SET(nom).Flow.HeartRate < 35
    force = false;
    defineheartrate(gui.nom,force);
  end

end

%----------------------------
function defineheartrateingui
%----------------------------
%Define heartrate called from gui

global NO

% update heart rate
force = true;
defineheartrate(NO,force);

% update all graphicals elements based on the new heart rate
init(NO);

%---------------------------------
function defineheartrate(no,force)
%---------------------------------
%Ask user to manually define the heart reate by giving the number of heart
%beats

% persistent asked 

global SET DATA NO  

% if isempty(asked)
%   asked = false;
% end

if nargin < 2
  force = false; 
end

% %Check if we should ignore
% if asked && ~force
%   return
% end

logdisp('defineheartrate')

if nargin < 1
  no = NO;
end

if isempty(SET(no).Flow)
  myfailed(dprintf('No flow data in current image stack.'));
  return
end

% asked = true;

if isfield(SET(no), 'OrigHeartRate') && ~isempty(SET(no).OrigHeartRate) && (SET(no).OrigHeartRate > 0)
  if SET(no).Flow.HeartRate < SET(no).OrigHeartRate
    SET(no).OrigHeartRate = SET(no).Flow.HeartRate;
  end
  temphr = SET(no).OrigHeartRate;
else
  temphr = (60/(SET(no).TSize*SET(no).TIncr));
end

if ~force
  yes = yesno('Low heart rate detected, possible real time images. Do you want to set number of heart beats?');
  if ~yes
    return
  end
end

%assume real time flow, ask user to define number of hearts beat
[nbrheartbeats,ok] = mygetnumber('Enter number of heart beats (decimal with .)', 'Heart beats',1.0,0,[]);
if not(ok)
  return;
end

if ~isempty(nbrheartbeats) && nbrheartbeats > 0 && isnumeric(nbrheartbeats)
  newheartrate = nbrheartbeats*temphr;
else
  myfailed('Incorrect input value for number of heart beats, using 1',DATA.GUI.Segment);
  newheartrate = temphr;
end

%Store it
tools('sethearrateok_Callback',newheartrate, no);

%Notify
mymsgbox(dprintf('The heart rate is now set to %0.3g',newheartrate),'Heart rate',DATA.GUI.Segment);

%-----------------------
function showflowroidata(no,flag)
%-----------------------
%To display flow (cm/s) from  ROI's  for non-time resolved Contrast Phase
%Images
global DATA SET

gui = DATA.GUI.Flow;
nom = gui.nom;
%nop = gui.nop;
rois2take = gui.rois2take;
numrois = gui.numrois;

%SET(no).Flow.Result(roinbr)
if isempty(SET(no).Flow.Result)
  close_Callback
  return
end

if flag == 0
  %Calculate Flow
  gui.netflow = zeros(1,numrois);
  gui.posflow = zeros(1,numrois);
  gui.negflow = zeros(1,numrois);
  gui.velmean = zeros(1,numrois);
  gui.velstd = zeros(1,numrois);
  gui.velmax = zeros(1,numrois);
  gui.velmin = zeros(1,numrois);
  gui.roiname = cell(1,numrois);
  gui.roinbr = zeros(1,numrois);
  gui.kenergy=0;
  gui.area=zeros(1,numrois);
  gui.volstri = sprintf('%s\t%s\t%s\t%s\n',dprintf('ROI Name'),dprintf('Total vol'),dprintf('Forward vol'),dprintf('Backward vol'));

  line = 3;
  for rloop=1:numrois
    %--- Interpolate up

    %Sum
    gui.netflow(rloop) = SET(no).Flow.Result(rloop).netflow;
    gui.posflow(rloop) = SET(no).Flow.Result(rloop).posflow;
    gui.negflow(rloop) = SET(no).Flow.Result(rloop).negflow;
    gui.velmean(rloop)  = SET(no).Flow.Result(rloop).velmean;
    gui.velstd(rloop)  = SET(no).Flow.Result(rloop).velstd;

    gui.velmax(rloop) = SET(no).Flow.Result(rloop).velmax;
    gui.velmin(rloop) = SET(no).Flow.Result(rloop).velmean;
    gui.area(rloop) =SET(no).Flow.Result(rloop).area;

    %Add to stri
    gui.volstri = [gui.volstri sprintf('%s\t%0.4g\t%0.4g\t%0.4g\n',...
      SET(nom).Roi(rois2take(rloop)).Name,...
      gui.netflow(rloop),...
      gui.posflow(rloop),...
      gui.negflow(rloop))];
  % 	gui.voltext = [gui.voltext dprintf('ROI: %s:\nNet: %0.3g ml\nNet: %0.3g l/min\nForward: %0.3g ml\nBackward: %0.3g ml\nRegurgitant fraction: %0.3g%%\n\n',...
  % 		SET(nom).Roi(rois2take(rloop)).Name,...
  % 		gui.nettotvol(rloop),...
  % 		gui.nettotvol(rloop)*gui.hr/1000,...
  % 		gui.netforwardvol(rloop),...
  % 		gui.netbackwardvol(rloop),...
  %     abs(100*gui.netbackwardvol(rloop)/gui.netforwardvol(rloop)))];
    gui.voltext{line} = sprintf('%s: %s:',dprintf('ROI'),SET(nom).Roi(rois2take(rloop)).Name);
    gui.voltext{line+1} = sprintf('%s: %0.3g ml/s',dprintf('Net flow'),gui.netflow(rloop));
    gui.voltext{line+2} = sprintf('%s: %0.3g ml/s',dprintf('Pos Flow'),gui.posflow(rloop));
    gui.voltext{line+3} = sprintf('%s: %0.3g ml/s',dprintf('Neg Flow'),gui.negflow(rloop));
    gui.voltext{line+4} = sprintf('%s: %0.3g cm/s',dprintf('Mean Vel'),gui.velmean(rloop));
    gui.voltext{line+5} = sprintf('%s: %0.3g cm/s',dprintf('Std Vel'),gui.velstd(rloop));
    gui.voltext{line+6} = sprintf('%s: %0.3g cm/s',dprintf('Max Vel'),gui.velmax(rloop));
    gui.voltext{line+7} = sprintf('%s: %0.3g cm/s',dprintf('Min Vel'),gui.velmin(rloop));
    gui.voltext{line+8} = sprintf('%s: %0.3g cm^2',dprintf('ROI Area'),gui.area(rloop));
    line = line+10;
  end
else
  line = 3;
  for rloop=1:numrois
	%--- Interpolate up


	%Add to stri
	gui.volstri = [gui.volstri sprintf('%s\t%0.4g\t%0.4g\t%0.4g\n',...
		SET(nom).Roi(rois2take(rloop)).Name,...
		gui.netflow(rloop),...
		gui.posflow(rloop),...
		gui.negflow(rloop))];
% 	gui.voltext = [gui.voltext dprintf('ROI: %s:\nNet: %0.3g ml\nNet: %0.3g l/min\nForward: %0.3g ml\nBackward: %0.3g ml\nRegurgitant fraction: %0.3g%%\n\n',...
% 		SET(nom).Roi(rois2take(rloop)).Name,...
% 		gui.nettotvol(rloop),...
% 		gui.nettotvol(rloop)*gui.hr/1000,...
% 		gui.netforwardvol(rloop),...
% 		gui.netbackwardvol(rloop),...
%     abs(100*gui.netbackwardvol(rloop)/gui.netforwardvol(rloop)))];
  gui.voltext{line} = sprintf('%s: %s:',dprintf('ROI'),SET(nom).Roi(rois2take(rloop)).Name);
  gui.voltext{line+1} = sprintf('%s: %0.3g cm/s',dprintf('Net flow'),gui.netflow(rloop));
	gui.voltext{line+2} = sprintf('%s: %0.3g cm/s',dprintf('Pos Flow'),gui.posflow(rloop));
	gui.voltext{line+3} = sprintf('%s: %0.3g cm/s',dprintf('Neg Flow'),gui.negflow(rloop));
  gui.voltext{line+4} = sprintf('%s: %0.3g cm/s',dprintf('Mean Vel'),gui.velmean(rloop));
  gui.voltext{line+5} = sprintf('%s: %0.3g cm/s',dprintf('Std Vel'),gui.velstd(rloop));
  gui.voltext{line+6} = sprintf('%s: %0.3g cm/s',dprintf('Max Vel'),gui.velmax(rloop));
  gui.voltext{line+7} = sprintf('%s: %0.3g cm/s',dprintf('Min Vel'),gui.velmin(rloop));
  gui.voltext{line+8} = sprintf('%s: %0.3g cm^2',dprintf('ROI Area'),gui.area(rloop));

  line = line+10;
  end
end
  
set(gui.handles.volumelistbox,'String',gui.voltext);


%-----------------------
function recalculate(no)
%-----------------------
%Calculate flow data
global DATA SET
gui = DATA.GUI.Flow;
nom = gui.nom;
nop = gui.nop;
if gui.rois2take ~= SET(nom).RoiN
  % Check if the number of time resolvedof ROIs has been updated outside the figure
  rois2take = [];
  for loop = 1:SET(nom).RoiN
    if sum(~isnan(SET(nom).Roi(loop).X(1,:)))>1
      %More than one timeframe
      rois2take = [rois2take loop]; %#ok<AGROW>
    end
  end
  if ~isempty(rois2take)
    gui.rois2take = rois2take;
    gui.numrois = length(rois2take);
  else
    return
  end
end
rois2take = gui.rois2take;
numrois = gui.numrois;

calcfunctions('calcflow',no);

if (~gui.resamped)  
  
	gui.velmean = NaN(SET(nom).TSize,numrois);
	gui.velstd = gui.velmean;
	gui.velmax = gui.velmean;
	gui.velmin = gui.velmean;
	gui.kenergy = gui.velmean;
	gui.area = gui.velmean;
	gui.netflow = gui.velmean;
	gui.posflow = gui.velmean;
	gui.negflow = gui.velmean;
  gui.phasecorrstring = '';
  if ~isempty(SET(nop).Flow.PhaseCorr)
    if ~isfield(SET(nop).Flow,'PhaseCorrTimeResolved')
      mywarning('Incompatible eddy current correction. Correction reset.',DATA.GUI.Segment);
      SET(nop).Flow.PhaseCorr = [];
      gui.phasecorrstring = '';
    else
      if SET(nop).Flow.PhaseCorrTimeResolved
        gui.phasecorrstring = dprintf('[Time-resolved eddy current compensation applied]');
      else
        gui.phasecorrstring = dprintf('[Stationary eddy current compensation applied]');
      end
    end
  end
  warnedempty = false;
  if DATA.Silent
    numiter = numrois+1;
  else
    numiter = 1;
  end
  h = mywaitbarstart(numrois,'Please wait, calculating flow.',numiter);
  for rloop = 1:numrois
    for tloop = SET(nom).Roi(rois2take(rloop)).T
      %Create mask
      mask = logical(segment('createmask',...
        gui.outsize,...
        SET(nom).Roi(rois2take(rloop)).Y(:,tloop),...
        SET(nom).Roi(rois2take(rloop)).X(:,tloop)));
      
      %Extract phase image
      temp = SET(nop).IM(:,:,tloop,SET(nom).Roi(rois2take(rloop)).Z);
      
      %If empty phasecorr, the do not add phase correction.
      if isempty(SET(nop).Flow.PhaseCorr)
        veldata = SET(nom).Roi(rois2take(rloop)).Sign*...
          (temp-0.5)*2*SET(nop).VENC;
      else
        %Phase correction
        if SET(nop).Flow.PhaseCorrTimeResolved
          %Time resolved phase correction
          veldata = SET(nom).Roi(rois2take(rloop)).Sign*...
            (temp-0.5-SET(nop).Flow.PhaseCorr(:,:,tloop,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(nop).VENC;
        else
          %Stationary phase correction
          veldata = SET(nom).Roi(rois2take(rloop)).Sign*...
            (temp-0.5-SET(nop).Flow.PhaseCorr(:,:,1,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(nop).VENC;
        end
      end
      
      veldata = veldata(mask);
      if isempty(veldata)
        if not(warnedempty)
          mywarning('Empty ROI detected. Should not occur.',DATA.GUI.Segment);
        end
        warnedempty = true;
      else
        %Okej to go
        posveldata = veldata(veldata>0);
        negveldata = veldata(veldata<0);
        gui.velmean(tloop,rloop) = mean(veldata);
        gui.velstd(tloop,rloop) = std(veldata);
        gui.velmax(tloop,rloop) = max(veldata);
        gui.velmin(tloop,rloop) = min(veldata);
        gui.kenergy(tloop,rloop) = sum((veldata/100).^3/2*(SET(nom).ResolutionX*SET(nom).ResolutionY/1e6)*1060); %kg/m^3
        gui.area(tloop,rloop) = SET(nom).ResolutionX*SET(nom).ResolutionY*sum(mask(:))/100; %cm^2
        gui.netflow(tloop,rloop) = (10/1000)*SET(nom).ResolutionX*SET(nom).ResolutionY*sum(veldata); %cm^3
        gui.posflow(tloop,rloop) = (10/1000)*SET(nom).ResolutionX*SET(nom).ResolutionY*sum(posveldata); %cm^3
        gui.negflow(tloop,rloop) = (10/1000)*SET(nom).ResolutionX*SET(nom).ResolutionY*sum(negveldata); %cm^3
      end
      
    end
    h = mywaitbarupdate(h);
  end
  mywaitbarclose(h);
  
  
  timeframes = SET(nom).StartAnalysis:SET(nom).EndAnalysis;
  %Calculate diameter equivalent
  gui.diameter = 2*10*sqrt(gui.area/pi); %mm*10 => cm
  TIncr = SET(nom).TIncr;
  TSize = SET(nom).TSize;
else
  gui.phasecorrstring = '';
  % 	StartAnalysis = 1;
  % 	EndAnalysis = length(gui.t);
  timeframes = 1:gui.tsize;
  TIncr = gui.TIncr;
  TSize = gui.tsize;
end


if (isempty(TIncr) || isempty(TSize))
  disp('TIncr or TSize structures are empty due to resampling error or bad SET!');
end

gui.volstri = sprintf('%s\t%s\t%s\t%s\n',dprintf('ROI Name'),dprintf('Total vol'),dprintf('Forward vol'),dprintf('Backward vol'));
gui.voltext{1} = dprintf('Time between green bars:%d ms\n\n',...
  round(1000*TIncr*(length(timeframes)-1)));
%--------------------------------------------------------------------------

%Calculate stroke volume
gui.nettotvol = zeros(1,numrois);
gui.nettotposvol = zeros(1,numrois);
gui.nettotnegvol = zeros(1,numrois);
gui.netforwardvol = zeros(1,numrois);
gui.netbackwardvol = zeros(1,numrois);
gui.meanflow = zeros(1,numrois);
gui.maxflow = zeros(1,numrois);
gui.maxmaxvel = zeros(1,numrois);
gui.roiname = cell(1,numrois);
gui.roinbr = zeros(1,numrois);
gui.hr = (60/(TSize*TIncr));

line = 3;
for rloop=1:numrois
	%--- Interpolate up
		
	%Sum
	gui.nettotvol(rloop) = sum(gui.netflow(timeframes,rloop), 'omitnan')*TIncr;
	gui.nettotposvol(rloop) = sum(gui.posflow(timeframes,rloop), 'omitnan')*TIncr;
	gui.nettotnegvol(rloop) = sum(gui.negflow(timeframes,rloop), 'omitnan')*TIncr;
	gui.netforwardvol(rloop) = sum(...
		gui.netflow(timeframes,rloop).*(gui.netflow(timeframes,rloop)>0), 'omitnan')*TIncr;
	gui.netbackwardvol(rloop) = sum(...
		gui.netflow(timeframes,rloop).*(gui.netflow(timeframes,rloop)<0), 'omitnan')*TIncr;
	gui.roiname{rloop}=SET(nom).Roi(rois2take(rloop)).Name;
	gui.roinbr(rloop)=rois2take(rloop);
	gui.sv(rloop)=gui.nettotvol(rloop)/gui.hr*60;
	gui.meanflow(rloop)=sum(gui.netflow(timeframes,rloop), 'omitnan')/length(timeframes)*60/1000;
	gui.maxflow(rloop)=max(gui.netflow(:,rloop))*60/1000;
	gui.maxmaxvel(rloop)=max(gui.velmax(:,rloop));
	gui.meanmeanvel(rloop)=mynanmean(gui.velmean(:,rloop));
	%Add to stri
	gui.volstri = [gui.volstri sprintf('%s\t%0.4g\t%0.4g\t%0.4g\n',...
		SET(nom).Roi(rois2take(rloop)).Name,...
		gui.nettotvol(rloop),...
		gui.netforwardvol(rloop),...
    gui.netbackwardvol(rloop))];
  % 	gui.voltext = [gui.voltext dprintf('ROI: %s:\nNet: %0.3g ml\nNet: %0.3g l/min\nForward: %0.3g ml\nBackward: %0.3g ml\nRegurgitant fraction: %0.3g%%\n\n',...
  % 		SET(nom).Roi(rois2take(rloop)).Name,...
  % 		gui.nettotvol(rloop),...
  % 		gui.nettotvol(rloop)*gui.hr/1000,...
  % 		gui.netforwardvol(rloop),...
  % 		gui.netbackwardvol(rloop),...
  %     abs(100*gui.netbackwardvol(rloop)/gui.netforwardvol(rloop)))];
  gui.voltext{line} = sprintf('%s: %s:',dprintf('ROI'),SET(nom).Roi(rois2take(rloop)).Name);
  gui.voltext{line+1} = sprintf('%s: %0.3g ml',dprintf('Net vol'),gui.nettotvol(rloop)); 
  gui.voltext{line+2} = sprintf('%s: %0.3g ml',dprintf('Forward vol'),gui.netforwardvol(rloop));
  gui.voltext{line+3} = sprintf('%s: %0.3g ml',dprintf('Backward vol'),gui.netbackwardvol(rloop)); 
  gui.voltext{line+4} = sprintf('%s: %0.3g%%',dprintf('Regurgitant fraction'),abs(100*gui.netbackwardvol(rloop)/gui.netforwardvol(rloop)));
  gui.voltext{line+5} = sprintf('%s: %0.3g L/min',dprintf('FlowCO'),gui.nettotvol(rloop)*gui.hr/1000);
  
  line = line+7;
end
set(gui.handles.volumelistbox,'String',gui.voltext);

set(gui.fig,'pointer','arrow');
update(nom); %Also do an update

exporttosetstruct;
for loop = 1:SET(nom).RoiN
  SET(nom).RoiCurrent = loop;
  calcfunctions('calcflow',nom);
end
segment('updateflow'); %DATA.updateaxestables('flow',nom,nop);

grid(gui.handles.plotaxes,'on');
gui.handles.plotaxes.GridColor = 'k';
set(gui.handles.plotaxes,'Color',[0.94 0.94 0.94]);

%---------------------
function zoom_Callback 
%---------------------
%Callback for checkbox toggling zoom on/off
global DATA
gui = DATA.GUI.Flow;
if isequal(get(gui.handles.zoomcheckbox,'value'),1)
	h=zoom(gui.fig);
	set(h,'enable','on');
else
	h=zoom(gui.fig);
	set(h,'enable','off');
end

%-------------------------------
function s = getvisiblestring(n)
%-------------------------------
%Get if a ROI is visible or not

global DATA

gui = DATA.GUI.Flow;

v = gui.show(n); %Get if it is shown or not
if v %return as a string suitable for visible
  s = 'on';
else
  s = 'off';
end

%--------------------------------
function varargout = getdata(arg)
%--------------------------------
%Output requested data from flow calculations
global DATA
gui = DATA.GUI.Flow;
varargout = cell(1,nargout); %initialise as empty
value = [];
switch arg
  case 'getalldata'
    value = gui;
  case 'gettotal'
    value = gui.nettotvol;
  case 'getfwdflow'
    value = gui.netforwardvol;
  case 'getbwdflow'
    value = gui.netbackwardvol;
  case 'getmeanflow'
    value = mean(gui.netflow);
  case 'getpeakflow'
    value = max(gui.netflow);
  case 'getvelocity'
    value = gui.velmean;
  case 'getmeanvelocity'
    value = mean(gui.velmean,1);
  case 'getpeakvelocity'
    value = max(gui.velmax);
    varargout{2} = min(gui.velmin);
  case 'getstrokevolume'
    value = gui.nettotvol.*gui.hr/1000;
  case 'getdistensibility'
    value = gui.diameter/10;%cm/10=mm
  case 'getcurve'
    value = gui.netflow;
  case 'getroiname'
    value = gui.roiname;
  case 'getpeakvelocitycurve'
    value = gui.velmax;
end
varargout{1} = value;

%-------------------------
function exporttosetstruct
%-------------------------
global DATA SET
gui = DATA.GUI.Flow;

no = gui.nom;
for roinbr = 1:gui.numrois%gui.rois2take
  flowdata = [];
  flowdata.velmean = gui.velmean(:,roinbr);
  flowdata.velstd = gui.velstd(:,roinbr);
  flowdata.velmax = gui.velmax(:,roinbr);
  flowdata.velmin = gui.velmin(:,roinbr);
  flowdata.kenergy = gui.kenergy(:,roinbr);
  
  flowdata.netflow = gui.netflow(:,roinbr);
  flowdata.posflow = gui.posflow(:,roinbr);
  flowdata.negflow = gui.negflow(:,roinbr);
  flowdata.diameter = gui.diameter(:,roinbr);
  
  flowdata.nettotvol = gui.nettotvol(roinbr);
  flowdata.nettotposvol = gui.nettotposvol(roinbr);
  flowdata.nettotnegvol = gui.nettotnegvol(roinbr);
  flowdata.netforwardvol = gui.netforwardvol(roinbr);
  flowdata.netbackwardvol = gui.netbackwardvol(roinbr);
  
  flowdata.meanflow = gui.meanflow(roinbr);
  flowdata.maxflow = gui.maxflow(roinbr);
  flowdata.maxmaxvel = gui.maxmaxvel(roinbr);
  flowdata.sv = gui.sv(roinbr);
  flowdata.meanmeanvel = gui.meanmeanvel(roinbr);
  
  SET(no).Roi(gui.rois2take(roinbr)).Flow = flowdata;
end 

%----------------------
function close_Callback 
%-----------------------
%Close flow report GUI

global DATA

datacursormode off;

try
  DATA.GUI.Flow = close(DATA.GUI.Flow);  %close the flow gui
catch %#ok<CTCH>
  delete(gcbf)
end

DATA.GUI.Flow = [];
segment('updateflow'); %DATA.updateaxestables('flow');

%---------------------------
function hidestatic_Callback
%---------------------------
%Hides / shows static ROIs

global DATA

gui = DATA.GUI.Flow;

v = get(gui.handles.hidestaticcheckbox,'Value');

%Set all to true (i.e visible)
gui.show = true(1,gui.numrois);

if v
  %Show static
  for loop = 1:gui.numrois
    if contains(lower(gui.roiname{loop}),'static')
      gui.show(loop) = false;
    end
  end
end

update();

%---------------------------
function selectrois_Callback
%---------------------------
%Select what ROIs to plot

global DATA

gui = DATA.GUI.Flow;

s = [];
s(1).Field = 'Text';
s(1).Label = dprintf('Select ROIs to show');
s(1).Default = dprintf('Select ROIs to show');
s(1).Type = 'text';

for loop = 1:gui.numrois
  s(loop+1).Field = sprintf('Field%d',loop); %#ok<AGROW> 
  s(loop+1).Label = gui.roiname{loop}; %#ok<AGROW> 
  s(loop+1).Default = gui.show(loop); %#ok<AGROW>
end

[outs,ok] = myinputstruct(s,dprintf('Select ROIs'));
if ~ok
  return
end

%Here we go
gui.show = true(1,gui.numrois); %Set all to true (i.e visible)
for loop = 1:gui.numrois
  if ~outs.(sprintf('Field%d',loop))
    gui.show(loop) = false;
  end
end

set(gui.handles.hidestaticcheckbox,'Value',false); %Turn it off
update();

%--------------
function update(no)
%--------------
%Make graphical update of plot
global DATA SET NO

gui = DATA.GUI.Flow;
t = gui.t;
d = gui.d;
dt = gui.dt;
nom = gui.nom;
rois2take = gui.rois2take;
numrois = gui.numrois;
if nargin < 1
  no = NO;
end
  
%check if there has been roi updates in the main gui
if ~(gui.nop==no) && ~(gui.nom==no)
  reportflow('init');
elseif ~(numrois==length(SET(no).Roi))
  %reportflow('init');
end

v = mygetlistbox(gui.handles.parameterlistbox);
timestr = makeunitstring(dprintf('Time'),'ms');
switch v
  case 1 %'Flow'
    %--- Plot flow
    h = plot(gui.handles.plotaxes,t*1000,gui.netflow(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);    
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on')
    for loop=2:size(gui.netflow,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.netflow(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    ttlstr = sprintf('%s. %s',dprintf('Net flow'),gui.phasecorrstring);
    title(gui.handles.plotaxes,ttlstr,'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(dprintf('Flow'),'ml/s'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
        
  case 2 %'Positive/Negative Flow'
    %--- plot forward/backward flow
    h = plot(gui.handles.plotaxes,t*1000,gui.posflow(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on');
    for loop=2:size(gui.netflow,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.posflow(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    for loop=1:size(gui.negflow,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.negflow(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'linestyle',':','Visible',getvisiblestring(loop));
    end
    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    title(gui.handles.plotaxes,dprintf('Positive (-) / Negative (:) flow. %s',gui.phasecorrstring),'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(dprintf('Flow'),'ml/s'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
  case 3 %'Velocity'
    %--- Plot velocity
    h = plot(gui.handles.plotaxes,t*1000,gui.velmean(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on')
    for loop=2:size(gui.velmean,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.velmean(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    
    %--- Plot std bars
    h1 = plot(...
      gui.handles.plotaxes,...
      [t'*1000 t'*1000]',...
      [gui.velmean(:,1)-gui.velstd(:,1) gui.velmean(:,1)+gui.velstd(:,1)]',...
      SET(nom).Roi(1).LineSpec);
    h2 = plot(...
      gui.handles.plotaxes,...
      [t'*1000-dt t'*1000+dt]',...
      [gui.velmean(:,1)-gui.velstd(:,1) gui.velmean(:,1)-gui.velstd(:,1)]',...
      SET(nom).Roi(1).LineSpec);
    h3 = plot(...
      gui.handles.plotaxes,...
      [t'*1000-dt t'*1000+dt]',...
      [gui.velmean(:,1)-gui.velstd(:,1) gui.velmean(:,1)-gui.velstd(:,1)]',...
      SET(nom).Roi(1).LineSpec);
    h4 = plot(...
      gui.handles.plotaxes,...
      [t'*1000-dt t'*1000+dt]',...
      [gui.velmean(:,1)+gui.velstd(:,1) gui.velmean(:,1)+gui.velstd(:,1)]',...
      SET(nom).Roi(1).LineSpec);
    set([h1 h2 h3 h4],'Visible',getvisiblestring(1));

    for loop=2:size(gui.velmean,2)
      h1 = plot(...
        gui.handles.plotaxes,...
        [t'*1000 t'*1000]',...
        [gui.velmean(:,loop)-gui.velstd(:,loop) gui.velmean(:,loop)+gui.velstd(:,loop)]',...
        SET(nom).Roi(loop).LineSpec);
      h2 = plot(...
        gui.handles.plotaxes,...
        [t'*1000-dt t'*1000+dt]',...
        [gui.velmean(:,loop)-gui.velstd(:,loop) gui.velmean(:,loop)-gui.velstd(:,loop)]',...
        SET(nom).Roi(loop).LineSpec);
      h3 = plot(...
        gui.handles.plotaxes,...
        [t'*1000-dt t'*1000+dt]',...
        [gui.velmean(:,loop)-gui.velstd(:,loop) gui.velmean(:,loop)-gui.velstd(:,loop)]',...
        SET(nom).Roi(loop).LineSpec);
      h4 = plot(...
        gui.handles.plotaxes,...
        [t'*1000-dt t'*1000+dt]',...
        [gui.velmean(:,loop)+gui.velstd(:,loop) gui.velmean(:,loop)+gui.velstd(:,loop)]',...
        SET(nom).Roi(1).LineSpec);
      set([h1 h2 h3 h4],'Visible',getvisiblestring(loop));
    end

    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    paramstr = dprintf('Mean velocity');
    title(gui.handles.plotaxes,sprintf('%s. %s',paramstr,gui.phasecorrstring),'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(paramstr,'cm/s'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
  case 4 %'Max Velocity'
    %--- Plot max
    h = plot(gui.handles.plotaxes,t*1000,gui.velmax(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on');
    for loop=2:size(gui.area,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.velmax(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    paramstr = dprintf('Max velocity');
    title(gui.handles.plotaxes,sprintf('%s. %s',paramstr,gui.phasecorrstring),'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(paramstr, 'cm/s'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
  case 5 %'Min Velocity'
    %--- Plot min
    h = plot(gui.handles.plotaxes,t*1000,gui.velmin(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on');
    for loop=2:size(gui.area,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.velmin(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    paramstr = dprintf('Min velocity');
    title(gui.handles.plotaxes,sprintf('%s. %s',paramstr,gui.phasecorrstring),'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(paramstr, 'cm/s'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
  case 6 %'Signed Kinetic Energy'
    %--- Plot KE
    h = plot(gui.handles.plotaxes,t*1000,gui.kenergy(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on');
    for loop=2:size(gui.area,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.kenergy(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    paramstr = dprintf('Signed Kinetic Energy');
    title(gui.handles.plotaxes,sprintf('%s. %s',paramstr,gui.phasecorrstring),'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(paramstr,'Nm/s'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
  case 7 %'Area'
    %--- Plot area
    h = plot(gui.handles.plotaxes,t*1000,gui.area(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on');
    for loop=2:size(gui.area,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.area(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    set(gui.handles.plotaxes,'ylim',[0 1.2*max(max(gui.area))]);
    paramstr = dprintf('Area');    
    title(gui.handles.plotaxes,sprintf('%s. %s',paramstr,gui.phasecorrstring),'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(paramstr,'cm^2'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
  case 8 %'Diameter'
    %--- Plot diameter
    h = plot(gui.handles.plotaxes,t*1000,gui.diameter(:,1),[d SET(nom).Roi(rois2take(1)).LineSpec]);
    set(h,'linewidth',2,'Visible',getvisiblestring(1));
    hold(gui.handles.plotaxes,'on');
    for loop=2:size(gui.diameter,2)
      h = plot(gui.handles.plotaxes,t*1000,gui.diameter(:,loop),[d SET(nom).Roi(rois2take(loop)).LineSpec]);
      set(h,'linewidth',2,'Visible',getvisiblestring(loop));
    end
    set(gui.handles.plotaxes,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
    hold(gui.handles.plotaxes,'off');
    set(gui.handles.plotaxes,'xlim',[0 (SET(nom).TSize-1)*1000*SET(nom).TIncr]);
    set(gui.handles.plotaxes,'ylim',[0 1.2*max(max(gui.diameter))]);
    paramstr = dprintf('Diameter');
    title(gui.handles.plotaxes,sprintf('%s. %s',paramstr,gui.phasecorrstring),'FontSize',12,'Color',DATA.GUISettings.ForegroundColor);
    xlabel(gui.handles.plotaxes,timestr,'Color',DATA.GUISettings.ForegroundColor);
    ylabel(gui.handles.plotaxes,makeunitstring(paramstr,'mm'),'Color',DATA.GUISettings.ForegroundColor);
    grid(gui.handles.plotaxes,'on');
  case {9,10,11,12} %{'3D plot','3D plot one frame','Pixel Export'}
    %--- 3D plot
    if v == 10 %strcmp(s{v},'3D plot one frame')
      rois=SET(nom).RoiCurrent;
    else
      rois=1:numrois;
    end
    for rloop=rois
      
      timeframes = SET(nom).StartAnalysis:SET(nom).EndAnalysis;
      
      %Create mask
      mask = false([SET(nom).XSize SET(nom).YSize length(timeframes) 1]);
      left = SET(nom).YSize;
      right = 1;
      up = SET(nom).XSize;
      down = 1;
      for tloop=timeframes
        if ~isempty(find(SET(nom).Roi(rois2take(rloop)).T==tloop,1))
          mask(:,:,tloop-timeframes(1)+1) = segment('createmask',...
            [SET(nom).XSize SET(nom).YSize],...
            SET(nom).Roi(rois2take(rloop)).Y(:,tloop),...
            SET(nom).Roi(rois2take(rloop)).X(:,tloop));
          
          %Find left/right extent
          temp = sum(mask(:,:,tloop-timeframes(1)+1),1);
          pos = find(temp);
          if ~isempty(pos)
            left = min(left,pos(1));
            right = max(right,pos(end));
          end
          
          %Find up/down extent
          temp = sum(mask(:,:,tloop-timeframes(1)+1),2);
          pos = find(temp);
          if ~isempty(pos)
            up = min(up,pos(1));
            down = max(down,pos(end));
          end
        end
      end
      
      %Add extra space around
      left = max(left-1,1);
      right = min(right+1,SET(nom).YSize);
      up = max(up-1,1);
      down = min(down+1,SET(nom).XSize);
      
      x = (up:down)-up+1;
      y = (left:right)-left+1;
      xs = round(sqrt(length(timeframes)));
      
      %Find scaling
      scaling = double((down-up)*xs/SET(nom).VENC);
      
      nop = SET(nom).Flow.PhaseNo;
      
      if v == 10 %strcmp(s{v},'3D plot one frame')
        %3D plot option
        flow('flow3dmovie_Callback',mask,double(scaling),nop,up,down,left,right);
      end
      
      if v == 9 %strcmp(s{v},'3D plot')
        %Plot all 3d plots in one figure.
        fig = figure(18+rloop);
        setupicon(fig);
        set(fig,...
          'Name',dprintf('3D velocity profile %s',SET(nom).Roi(rois2take(rloop)).Name),...
          'Numbertitle','off');
        clf;
        
        hold on;
        for tloop=1:length(timeframes)
          
          if isempty(SET(nop).Flow.PhaseCorr)
            %No phase correction
            temp = SET(nom).Roi(rois2take(rloop)).Sign*mask(:,:,tloop).*...
              (SET(nop).IM(:,:,timeframes(tloop),SET(nom).Roi(rois2take(rloop)).Z)-0.5)*2*SET(SET(nom).Flow.PhaseNo).VENC;
          else
            %Phase correction
            if SET(nop).Flow.PhaseCorrTimeResolved
              %Time resolved phase correction
              temp = SET(nom).Roi(rois2take(rloop)).Sign*mask(:,:,tloop).*...
                (SET(nop).IM(:,:,timeframes(tloop),SET(nom).Roi(rois2take(rloop)).Z)-0.5-...
                SET(nop).Flow.PhaseCorr(:,:,tloop,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(SET(nom).Flow.PhaseNo).VENC;
            else
              %Stationary phase correction
              temp = SET(nom).Roi(rois2take(rloop)).Sign*mask(:,:,tloop).*...
                (SET(nop).IM(:,:,timeframes(tloop),SET(nom).Roi(rois2take(rloop)).Z)-0.5-...
                SET(nop).Flow.PhaseCorr(:,:,1,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(SET(nom).Flow.PhaseNo).VENC;
            end
          end
          
          temp = temp(up:down,left:right);
          temp = double(temp);
          dx = mod(tloop-1,xs);
          dy = ceil(tloop/xs)-1;
          h = surf(y+dy*y(end),x+dx*x(end),0.25*temp*scaling); %#ok<NASGU>
        end
        hold off;
        axis off image vis3d;
        rotate3d on;
        cameratoolbar(fig);
      end
      
      if v == 11 %strcmp(s{v},'Pixel Export')
        %Export all pixeldata.
        
        %Reserve memory
        ysize = length(up:down);
        xsize = length(left:right);
        outdata = cell((ysize+2)*length(timeframes)+11,xsize+2);
        rowoffset = 12;
        outdata{1,1} = 'RoiName:';
        outdata{1,2} = SET(nom).Roi(rois2take(rloop)).Name;
        outdata{3,2} = 'RL[mm]';
        outdata{3,3} = 'AP[mm]';
        outdata{3,4} = 'FH[mm]';
        outdata{4,1} = 'TopLeftCorner';
        outdata{5,1} = 'TopRightCorner';
        outdata{6,1} = 'BottomRightCorner';
        outdata{7,1} = 'BottomLeftCorner';
        outdata{8,1} = 'ResolutionRows[mm]:';
        outdata{8,2} = SET(nom).ResolutionY;
        outdata{9,1} = 'ResolutionCols[mm]:';
        outdata{9,2} = SET(nom).ResolutionX;
        outdata{10,1} = 'VelocityScale[cm/s]';
        
        %top left corner
        pos = calcfunctions('xyz2rlapfh',...
          nom,...
          up,...
          left,...
          SET(nom).Roi(rois2take(rloop)).Z);
        outdata{4,2} = pos(1);
        outdata{4,3} = pos(2);
        outdata{4,4} = pos(3);
        
        %top right corner
        pos = calcfunctions('xyz2rlapfh',...
          nom,...
          up,...
          right,...
          SET(nom).Roi(rois2take(rloop)).Z);
        outdata{5,2} = pos(1);
        outdata{5,3} = pos(2);
        outdata{5,4} = pos(3);
        
        %Bottom right corner
        pos = calcfunctions('xyz2rlapfh',...
          nom,...
          down,...
          right,...
          SET(nom).Roi(rois2take(rloop)).Z);
        outdata{6,2} = pos(1);
        outdata{6,3} = pos(2);
        outdata{6,4} = pos(3);
        
        %Bottom left corner
        pos = calcfunctions('xyz2rlapfh',...
          nom,...
          down,...
          left,...
          SET(nom).Roi(rois2take(rloop)).Z);
        outdata{7,2} = pos(1);
        outdata{7,3} = pos(2);
        outdata{7,4} = pos(3);
        
        for tloop=1:length(timeframes)
          
          if isempty(SET(nop).Flow.PhaseCorr)
            %No phase correction
            temp = SET(nom).Roi(rois2take(rloop)).Sign*mask(:,:,tloop).*...
              (SET(nop).IM(:,:,timeframes(tloop),SET(nom).Roi(rois2take(rloop)).Z)-0.5)*2*SET(SET(nom).Flow.PhaseNo).VENC;
          else
            %Phase correction
            if SET(nop).Flow.PhaseCorrTimeResolved
              %Time resolved phase correction
              temp = SET(nom).Roi(rois2take(rloop)).Sign*mask(:,:,tloop).*...
                (SET(nop).IM(:,:,timeframes(tloop),SET(nom).Roi(rois2take(rloop)).Z)-0.5-...
                SET(nop).Flow.PhaseCorr(:,:,tloop,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(SET(nom).Flow.PhaseNo).VENC;
            else
              %Stationary phase correction
              temp = SET(nom).Roi(rois2take(rloop)).Sign*mask(:,:,tloop).*...
                (SET(nop).IM(:,:,timeframes(tloop),SET(nom).Roi(rois2take(rloop)).Z)-0.5-...
                SET(nop).Flow.PhaseCorr(:,:,1,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(SET(nom).Flow.PhaseNo).VENC;
            end
          end
          
          temp = temp(up:down,left:right);
          temp = double(temp);
          outdata(rowoffset+(tloop-1)*(ysize+2)+(1:ysize),1:xsize) = num2cell(temp);
          outdata{rowoffset+(tloop-1)*(ysize+2),1} = 'Time [ms]:';
          outdata{rowoffset+(tloop-1)*(ysize+2),2} = (tloop-1)*SET(nom).TIncr*1000;
        end
        segment('cell2clipboard',outdata);
      end %v=11
      
    end
end %Case clause

%Add the gren side bars and the red baseline bar, legend, not 3D
xstart = SET(nom).TIncr*(SET(nom).StartAnalysis-1)*1000;
xend = SET(nom).TIncr*(SET(nom).EndAnalysis-1)*1000;
if mygetlistbox(gui.handles.parameterlistbox)<=(length(get(gui.handles.parameterlistbox,'String'))-2)

  displaybaselinebar = get(gui.handles.baselinebarcheckbox,'Value');

  %Add zero line
  hzeroline = yline(0, 'k-', 'LineWidth', 2,'parent',gui.handles.plotaxes); % Horizontal zero line
  set(hzeroline, 'HandleVisibility', 'off'); % Exclude from legend
  if displaybaselinebar
    set(hzeroline, 'ButtonDownFcn','reportflow(''baselinebar_Buttondown'')');
  end
  gui.handles.hzeroline = hzeroline;

  ylim = get(gui.handles.plotaxes,'ylim');
  
  hold(gui.handles.plotaxes,'on');
  h = plot(gui.handles.plotaxes,[xstart xstart],ylim,'g-');
  set(h,'linewidth',2,'ButtonDownFcn','reportflow(''sidelinebar_Buttondown'',''startbar'')');
  gui.handles.startbar = h;
  h = plot(gui.handles.plotaxes,[xend xend],ylim,'g-');
  set(h,'linewidth',2,'ButtonDownFcn','reportflow(''sidelinebar_Buttondown'',''endbar'')');
  gui.handles.endbar = h;

  if displaybaselinebar
    %--- Code for baselinebar
    baselinelevel = computebaselinelevel(gui.netflow); %Compute baselinelevel

    if isequal(mygetlistbox(gui.handles.parameterlistbox),1)
      h = plot(gui.handles.plotaxes,1000*[t(1) t(end)],[baselinelevel baselinelevel],'r-');
      set(h,'linewidth',2,'ButtonDownFcn','reportflow(''baselinebar_Buttondown'')');
      gui.handles.baselinebar = h;
    end
  end

  hold(gui.handles.plotaxes,'off');
  
  %--- Legend
  temp = cell(1,size(gui.diameter,2));
  for loop=1:size(gui.diameter,2)
    if SET(nom).Roi(rois2take(loop)).Sign<0
      temp{loop} = [SET(nom).Roi(rois2take(loop)).Name ' (Sign switched)'];
    else
      temp{loop} = SET(nom).Roi(rois2take(loop)).Name;
    end
  end
  legend(gui.handles.plotaxes,temp{:});
end

gui.handles.plotaxes.GridColor = 'k';
set(gui.handles.plotaxes,'Color',[0.94 0.94 0.94]);

%Work around to change yellow lines to orange
helperfunctions('changelinecolorfromyellowtoorange',gui.handles.plotaxes);

%--------------
function export 
%--------------
%Export flow data to clipboard
global DATA SET
gui = DATA.GUI.Flow;
nom = gui.nom;
numrois = gui.numrois;
rois2take = gui.rois2take;
tincr = SET(nom).TIncr;

if (SET(nom).EndAnalysis - SET(nom).StartAnalysis) ==0 %For non-time resolved Phase Contrast Images
  %Create output matrix
  outdata = cell(numrois,9);
  outdata{1, 1} = 'ROI no';
  outdata{1, 2} = 'ROI name';
  outdata{1, 3} = 'Total flow [ml/s]';
  outdata{1, 4} = 'Positive flow [ml/s]';
  outdata{1, 5} = 'Negative flow [ml/s]';
  outdata{1, 6} = 'Mean vel [cm/s]';
  outdata{1, 7} = 'SD of vel [cm/s]';
  outdata{1, 8} = 'Min velocity [cm/s]';
  outdata{1, 9} = 'Max velocity [cm/s]';
  outdata{1, 10} = 'Area[cm^2]';

  for loop=1:numrois
    outdata{1+loop, 1} = loop;
    outdata{1+loop, 2} = SET(nom).Roi(rois2take(loop)).Name;
    outdata{1+loop, 3} = gui.netflow(loop);
    outdata{1+loop, 4} = gui.posflow(loop);
    outdata{1+loop, 5} = gui.negflow(loop);
    outdata{1+loop, 6} = gui.velmean(loop);
    outdata{1+loop, 7} = gui.velstd(loop);
    outdata{1+loop, 8} = gui.velmin(loop);
    outdata{1+loop, 9} = gui.velmax(loop);
    outdata{1+loop, 10} = gui.area(loop); 
  end
  
else %For time resolved Phase Contrast Images
if gui.resamped
	timeframes = 1:gui.tsize;
	tincr = gui.TIncr;
else
	timeframes = SET(nom).StartAnalysis:SET(nom).EndAnalysis;
end
	outdata = cell(1+length(timeframes)+5,1+9*SET(nom).RoiN);

%--- First header
outdata{2,1} = 'Time[ms]';
for loop=1:numrois
  outdata{1,2+(loop-1)*9} = SET(nom).Roi(rois2take(loop)).Name;
end

%--- Second header
for loop=1:numrois
  outdata{2,2+(loop-1)*9} = 'Area[cm^2]';
  outdata{2,3+(loop-1)*9} = 'Mean vel [cm/s]';
  outdata{2,4+(loop-1)*9} = 'SD of vel [cm/s]';
  outdata{2,5+(loop-1)*9} = 'Min velocity [cm/s]';
  outdata{2,6+(loop-1)*9} = 'Max velocity [cm/s]';
  outdata{2,7+(loop-1)*9} = 'Signed Kinetic Energy [Nm/s]';
  outdata{2,8+(loop-1)*9} = 'Total flow [ml/s]';
  outdata{2,9+(loop-1)*9} = 'Positive flow [ml/s]';
  outdata{2,10+(loop-1)*9} = 'Negative flow [ml/s]';
end

%--- Data
for tloop=timeframes
  row = 3+tloop-timeframes(1);
  outdata{row,1} = (tloop-1)*tincr*1000;
  for loop=1:numrois
    outdata{row,2+(loop-1)*9} = gui.area(tloop,loop);
    outdata{row,3+(loop-1)*9} = gui.velmean(tloop,loop);
    outdata{row,4+(loop-1)*9} = gui.velstd(tloop,loop);
    outdata{row,5+(loop-1)*9} = gui.velmin(tloop,loop);
    outdata{row,6+(loop-1)*9} = gui.velmax(tloop,loop);
    outdata{row,7+(loop-1)*9} = gui.kenergy(tloop,loop);
    outdata{row,8+(loop-1)*9} = gui.netflow(tloop,loop);
    outdata{row,9+(loop-1)*9} = gui.posflow(tloop,loop);
    outdata{row,10+(loop-1)*9} = gui.negflow(tloop,loop);
  end
end

%--- Integral flow
row = row+2;
outdata{row,1} = 'Name';
outdata{row,2} = 'Net volume [ml]';
outdata{row,3} = 'Forward volume [ml]';
outdata{row,4} = 'Backward volume [ml]';
outdata{row,5} = 'Peak min velocity [cm/s]';
outdata{row,6} = 'Peak max velocity [cm/s]';
outdata{row,7} = 'Heart rate [bpm]';
outdata{row,8} = 'Time between bars [ms]';
outdata{row,9} = 'Cardiac output [L/min]';
% outdata{row,10} = 'Mean flow [L/min]';

row = row+1;
for loop=1:numrois
  outdata{row,1} = SET(nom).Roi(rois2take(loop)).Name;
  outdata{row,2} = sum(gui.netflow(timeframes,loop))*tincr;
  outdata{row,3} = tincr*sum(...
    gui.netflow(timeframes,loop).*(gui.netflow(timeframes,loop)>0));
  outdata{row,4} = tincr*sum(...
    gui.netflow(timeframes,loop).*(gui.netflow(timeframes,loop)<0));
  outdata{row,5} = min(gui.velmin(timeframes,loop));
  outdata{row,6} = max(gui.velmax(timeframes,loop));
  outdata{row,7} = gui.hr;
  outdata{row,8} = round(1000*tincr*(length(timeframes)-1));
  outdata{row,9} = gui.hr*sum(gui.netflow(timeframes,loop))*tincr/1000;
%   outdata{row,10} = (outdata{row,2}/1000)*(1000/outdata{row,8})*60; %SV [L] * bpm => L/min was incorrect :(outdata{row,8}/1000)*(outdata{row,2}/1000*60);
  %(net ml)/1000 => L
  %time ms)/1000*60 => [min]
  row = row+1;
end
end
segment('cell2clipboard',outdata);

%------------------------------------
function sidelinebar_Buttondown(type) 
%------------------------------------
%Button down function for baselinebar in flow report GUI.

set(gcf,'WindowButtonMotionFcn',sprintf('reportflow(''sidelinebar_Motion'',''%s'')',type));
set(gcf,'WindowButtonUpFcn',sprintf('reportflow(''sidelinebar_Buttonup'',''%s'')',type));

%----------------------------------
function sidelinebar_Buttonup(type) 
%----------------------------------
%Button up function for baselinebar in flow report GUI.

global SET DATA

gui = DATA.GUI.Flow;
nom = gui.nom;
set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');

%Get Convert to timeframe
switch type
  case 'startbar'
    h = gui.handles.startbar;
  case 'endbar'
    h = gui.handles.endbar;
end

x = get(h,'xdata');
x = x(1);
x = round(1+x/(1000*SET(nom).TIncr));
x = max(min(x,SET(nom).TSize),1);

switch type
  case 'startbar'
    SET(nom).StartAnalysis = x;
  case 'endbar'
    SET(nom).EndAnalysis = x;    
end

if SET(nom).StartAnalysis > SET(nom).EndAnalysis
  temp = SET(nom).StartAnalysis;
  SET(nom).StartAnalysis = SET(nom).EndAnalysis;
  SET(nom).EndAnalysis = temp;
end

%Recalculate to get correct flow
recalculate(nom);

%--------------------------------
function sidelinebar_Motion(type)
%--------------------------------
%Motion function for sidelinebar.

global DATA SET NO

[x,~] = mygetcurrentpoint(gca);

gui = DATA.GUI.Flow;

%Convert to timeframes and round
x = round(1+x/(1000*SET(NO).TIncr));
x = min(max(0,x),SET(NO).TSize);

%Convert back to ms
x = (x-1)*SET(NO).TIncr*1000;

%Update display
switch type
  case 'startbar'
    set(gui.handles.startbar,'xdata',[x x]);
  case 'endbar'
    set(gui.handles.endbar,'xdata',[x x]);
end

%------------------------------
function baselinebar_Buttondown 
%------------------------------
%Button down function for baselinebar in flow report GUI.

set(gcf,'WindowButtonMotionFcn','reportflow(''baselinebar_Motion'')');
set(gcf,'WindowButtonUpFcn','reportflow(''baselinebar_Buttonup'')');

%--------------------------
function baselinebar_Motion 
%--------------------------
%Motion function for baselinebar.

global DATA

[~,y] = mygetcurrentpoint(gca);

%Update display
gui = DATA.GUI.Flow;
x = get(gui.handles.baselinebar,'XData');
set(gui.handles.baselinebar,'XData',[x(1) x(end)],'YData',[y y]);

%-----------------------------------------
function y = computebaselinelevel(netflow)
%-----------------------------------------
%Compute baseline level

extractedflow = netflow(:);
sortedflow = sort(abs(extractedflow));

%Take the threshold as the 2/3 largest flow
flowthreshold = sortedflow(round(0.67*length(extractedflow))); 
ind = abs(extractedflow) < flowthreshold;

%Take the mean of those
lowflow = extractedflow(ind);
y = mean(lowflow);

%----------------------------
function baselinebar_Buttonup 
%----------------------------
%Button up function for flow bar in flow report GUI.

global DATA SET NO

set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');

gui = DATA.GUI.Flow;

phno = SET(NO).Flow.PhaseNo;

%Compute baselinelevel
oldbaselinelevel = computebaselinelevel(gui.netflow);
y = get(gui.handles.baselinebar,'ydata');
y = y(1);
d = oldbaselinelevel-y; %ml/s

area = sum(mean(gui.area)); %sum of the mean area
dvel = d/area;
dvel = dvel/SET(phno).VENC/2;

%Check if no phase correction computed
if isempty(SET(phno).Flow.PhaseCorr)
  SET(phno).Flow.PhaseCorr = zeros(SET(phno).XSize,SET(phno).YSize,SET(phno).ZSize);
end

SET(phno).Flow.PhaseCorr = SET(phno).Flow.PhaseCorr+dvel;

%Recalculate to get correct flow
recalculate(NO);

%---------------------------
function datacursor_Callback 
%---------------------------
%Function which enables datacursor probing on roi-analysis plots.
%/SB
global DATA
gui = DATA.GUI.Flow;

datacursormode(gui.fig);


%---------------------------
function switchsign_Callback 
%---------------------------
%switch flow sign of current roi and recalucate and update
global NO DATA SET

gui = DATA.GUI.Flow;
roi('roiswitchsign_Callback');



if SET(NO).TSize <= 1   %For non-time resolced Phase Contrast Images
  %switch all the signs where it has impact
  gui.velmean = - gui.velmean;
  tmpmax = -gui.velmax;
  gui.velmax = -gui.velmin;
  gui.velmin = tmpmax;
  gui.kenergy = -gui.kenergy; %kg/m^3
  gui.area = gui.area; %cm^2  %to keep positive sign in rOI area for non-time resolved phase contrast data
  gui.netflow = -gui.netflow;  %cm^3
  tmpposflow = -gui.posflow;  %cm^3
  gui.posflow = abs(gui.negflow);%cm^3
  gui.negflow = tmpposflow;
  
  showflowroidata(NO,1);
else
  %switch all the signs where it has impact
  gui.velmean = - gui.velmean;
  tmpmax = gui.velmax;
  gui.velmax = gui.velmin;
  gui.velmin = tmpmax;
  gui.kenergy = -gui.kenergy; %kg/m^3
  gui.area = -gui.area; %cm^2
  gui.netflow = -gui.netflow;  %cm^3
  tmpposflow = gui.posflow;  %cm^3
  gui.posflow = gui.negflow;%cm^3
  gui.negflow = tmpposflow;
  recalculate(NO);
end

%----------------------------------
function clearcompensation_Callback 
%----------------------------------
%Clear eddycurrent compensation

global SET NO 

phno = SET(NO).Flow.PhaseNo;

SET(phno).Flow.PhaseCorr = [];

recalculate(NO);

update;

%-----------------------------
function compensation_Callback 
%-----------------------------
%open the eddy current compensation gui and recalculate and update

global SET NO 

phno = SET(NO).Flow.PhaseNo;

if isempty(SET(phno).Flow.PhaseCorr)
  floweddycurrent('initsmall');
else
  floweddycurrent('init');
end

update;