function varargout = floweddycurrent(varargin)
%Eddy current code.

%Original code by Einar Heiberg.
%Refactured summer 2015 by Einar Heiberg
if nargin==0
  init;
  return;
end;

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%----------------------
function init(plotonok)
%----------------------
%Initialize the main eddy current GUI

global DATA

%Check for errors
[ok,no] = safetychecks;
if ~ok
  return;
end;

%Check if we should plot after user press ok
if nargin<1
  plotonok = false;
end;

%Initalize
DATA.GUI.EddyCurrent = mygui('floweddycurrent.fig');
gui = DATA.GUI.EddyCurrent;

initvariables(no);
gui.plotonok = plotonok;

%Initialize graphics
graphicalinit;

%Find static tissue
findstatic;

%Calculate
recalculate_Callback;

%---------------------------
function initsmall(plotonok) %#ok<DEFNU>
%---------------------------
%Initialize the small automated eddy current GUI

global DATA

%Check for errors
[ok,no] = safetychecks;
if ~ok
  return;
end;

%Check if we should plot after user press ok
if nargin<1
  plotonok = false;
end;

%Initalize
DATA.GUI.EddyCurrent = mygui('floweddycurrentsmall.fig');
gui = DATA.GUI.EddyCurrent;

%Hide gui
set(gui.fig,'visible','off');

initvariables(no);
gui.biggui = false; %Mark that the smaller GUI
gui.plotonok = plotonok;

%Find static tissue
findstatic;

%Calculate
recalculate_Callback; %This also calls update

%Show GUI
set(gui.fig,'visible','on');

%-------------------------
function initvariables(no)
%-------------------------
global DATA SET

gui = DATA.GUI.EddyCurrent;

gui.phaseno = no;
gui.magno = SET(no).Flow.MagnitudeNo;
gui.slice = SET(no).CurrentSlice;
gui.timeframe = SET(no).CurrentTimeFrame;
gui.maxvel = 1;
gui.phasecorr = [];
gui.logim = [];
gui.badfit = [];
gui.biggui = true; %When true it is the main gui which is initialized.
gui.ignorerois = false; %If true then we do not care about stationary tissue or non stationary tissue

if isfield(SET(gui.phaseno).Flow, 'PhaseCorrRemoveBadFit');
  gui.usebadfit = SET(gui.phaseno).Flow.PhaseCorrRemoveBadFit;
else
  gui.usebadfit = true;
end

if ~isfield(SET(gui.phaseno).Flow,'PhaseCorr')
  SET(gui.no).Flow.PhaseCorr = [];
end;

if ~isfield(SET(gui.phaseno).Flow,'PhaseCorrPercentiles')
  SET(gui.phaseno).Flow.PhaseCorrPercentiles = 0.6;
end;
if ~isfield(SET(gui.phaseno).Flow,'PhaseCorrMagnitudeThreshold')
  SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold = 0.1;
end;
if ~isfield(SET(gui.phaseno).Flow,'PhaseCorrMethod')
  SET(gui.phaseno).Flow.PhaseCorrMethod = 'linear';
end;
if ~isfield(SET(gui.phaseno).Flow,'PhaseCorrTimeResolved')
  SET(gui.phaseno).Flow.PhaseCorrTimeResolved = false;
end;
if ~isfield(SET(gui.phaseno).Flow,'PhaseCorrStaticTissueRois')
  SET(gui.phaseno).Flow.PhaseCorrStaticTissueRois = false;
end;

%---------------------
function graphicalinit
%---------------------
%Initialize graphics
global DATA SET

gui = DATA.GUI.EddyCurrent;

%This function only updates the big gui
if ~gui.biggui
  return;
end;

%Update boxes, sliders etc.
set(gui.handles.phaseslider,'Value',100*SET(gui.phaseno).Flow.PhaseCorrPercentiles);
set(gui.handles.phaseedit,'Value',100*SET(gui.phaseno).Flow.PhaseCorrPercentiles,'String',round(100*SET(gui.phaseno).Flow.PhaseCorrPercentiles));
imvals = SET(gui.magno).IM(:);
set(gui.handles.magnitudeslider, ...
  'value',SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold, ...
  'min',min(imvals(imvals > 0)));
set(gui.handles.magnitudeedit, ...
  'Value',SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold, ...
  'String',round(100*SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold)/100);

switch SET(gui.phaseno).Flow.PhaseCorrMethod
  case 'linear'
    v = 1;
  case 'quadratic'
    v = 2;
  case 'static'
    v = 1; %Backwards compability
    SET(gui.phaseno).Flow.PhaseCorrMethod = 'linear';
  case 'ge'
    v = 3;
  otherwise
    myfailed(dprintf('Unknown method %s.',SET(gui.phaseno).Flow.PhaseCorrMethod));
    return;
end;
set(gui.handles.methodlistbox,'value',v);
set(gui.handles.timeresolvedcheckbox,'value',double(SET(gui.phaseno).Flow.PhaseCorrTimeResolved));
set(gui.handles.statictissueroischeckbox,'value',double(SET(gui.phaseno).Flow.PhaseCorrStaticTissueRois));
set(gui.handles.usebadfitcheckbox,'value',double(gui.usebadfit));

if SET(gui.phaseno).ZSize>1
  set(gui.handles.sliceslider,'min',1,'max',SET(gui.phaseno).ZSize,'value',gui.slice);
else
  set([gui.handles.slicetext gui.handles.sliceslider],'visible','off');
end;

set(gui.handles.timeslider,'min',1,'max',SET(gui.phaseno).TSize,'value',gui.timeframe);
if ~SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  set([gui.handles.timetext gui.handles.timeslider],'visible','off');
end;

%-------------------------------
function [res,no] = safetychecks
%-------------------------------
global SET NO

res = false;
no = NO;

%Check so that is single as the code requires it.
imissingle = classcheckim(NO);
if not(imissingle)
  return;
end

%Error checking
if isempty(SET(NO).Flow)
  myfailed('No flow data.');
  return;
end;

if isequal(NO,SET(NO).Flow.MagnitudeNo)
  %Check if only one phase is available.
  temp = [...
    SET(NO).Flow.PhaseNo
    SET(NO).Flow.PhaseX
    SET(NO).Flow.PhaseY];
  if length(temp)==1
    no = temp;
  else
    myfailed('Can not do eddy current compensation on magnitude image with multiple phase images. Select a phase image stack.');
    return;
  end;
end;

%If this statement is reached then ok.
res = true;

%---------------------------
function ignorerois_Callback %#ok<DEFNU>
%---------------------------
%Turn on ignore ROIs mode. This will disable adding static and non static
%tissue ROIS

global DATA

gui = DATA.GUI.EddyCurrent;
gui.ignorerois = true;

%Find static tissue
findstatic;

%Calculate
recalculate_Callback;

%---------------------------
function timeslider_Callback %#ok<DEFNU>
%---------------------------
%Timeslider Callback

global DATA SET

gui = DATA.GUI.EddyCurrent;

%Find temporal
SET(gui.phaseno).Flow.PhaseCorrTimeResolved = (1==get(gui.handles.timeresolvedcheckbox,'value'));
if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  set([gui.handles.timetext gui.handles.timeslider],'visible','on');
  gui.timeframe = mygetvalue(gui.handles.timeslider);
  gui.timeframe = min(max(round(gui.timeframe),1),SET(gui.phaseno).TSize);
else
  set([gui.handles.timetext gui.handles.timeslider],'visible','off');
  gui.timeframe = 1;
end;

update;

%-----------------------
function usebadfit_Callback %#ok<DEFNU>
%-----------------------
% Remove bad fit pixels checkbox callback

global DATA SET
gui = DATA.GUI.EddyCurrent;

chkboxvalue = get(gui.handles.usebadfitcheckbox, 'Value');
gui.usebadfit = logical(chkboxvalue);

SET(gui.phaseno).Flow.PhaseCorrRemoveBadFit = gui.usebadfit;

if ~gui.usebadfit
  gui.badfit = [];
end

findstatic;
recalculate_Callback;
update;

%-----------------------
function static_Callback %#ok<DEFNU>
%-----------------------
%Static tissue checkbox callback

global DATA SET

gui = DATA.GUI.EddyCurrent;

SET(gui.phaseno).Flow.PhaseCorrStaticTissueRois = get(gui.handles.statictissueroischeckbox,'value');

findstatic;
recalculate_Callback;
update;

%-----------------------------
function sliceslider_Callback %#ok<DEFNU>
%-----------------------------
%Slice slider callback

global DATA SET

gui = DATA.GUI.EddyCurrent;

%Find slice
if SET(gui.phaseno).ZSize>1
  gui.slice = mygetvalue(gui.handles.sliceslider);
  gui.slice = min(max(round(gui.handles.slice),1),SET(gui.phaseno).ZSize);
else
  gui.slice = 1;
end;

update;

%------------------------------
function methodlistbox_Callback %#ok<DEFNU>
%------------------------------
%Update the method listbox

global SET DATA

gui = DATA.GUI.EddyCurrent;

%Find method
v = mygetvalue(gui.handles.methodlistbox);
switch v
  case 1
    SET(gui.phaseno).Flow.PhaseCorrMethod = 'linear';
  case 2
    SET(gui.phaseno).Flow.PhaseCorrMethod = 'quadratic';
  case 3
    SET(gui.phaseno).Flow.PhaseCorrMethod = 'ge';
end;

recalculate_Callback;

%----------------------------------
function magnitudeslider_Callback %#ok<DEFNU>
%----------------------------------
%Magnitude threshold slider

global SET DATA

gui = DATA.GUI.EddyCurrent;

SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold = mygetvalue(gui.handles.magnitudeslider);

%Update
set(gui.handles.magnitudeslider,'value',SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold);
set(gui.handles.magnitudeedit,'Value',SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold,'String',round(100*SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold)/100);

findstatic;
recalculate_Callback;
update;

%------------------------------
function magnitudeedit_Callback %#ok<DEFNU>
%------------------------------
%Magnitude edit

global SET DATA

gui = DATA.GUI.EddyCurrent;

SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold = max(min(str2num(get(gui.handles.magnitudeedit,'String')), ...
  get(gui.handles.magnitudeslider,'max')),get(gui.handles.magnitudeslider,'min'));

%Update
set(gui.handles.magnitudeslider,'value',SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold);
set(gui.handles.magnitudeedit,'Value',SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold,'String',round(100*SET(gui.phaseno).Flow.PhaseCorrMagnitudeThreshold)/100);

%----------------------------
function phaseslider_Callback %#ok<DEFNU>
%----------------------------
%Phase threshold slider

global SET DATA

gui = DATA.GUI.EddyCurrent;

SET(gui.phaseno).Flow.PhaseCorrPercentiles = min(1,max(0,mygetvalue(gui.handles.phaseslider)/100));

%Update
set(gui.handles.phaseslider,'value',100*SET(gui.phaseno).Flow.PhaseCorrPercentiles);
set(gui.handles.phaseedit,'Value',100*SET(gui.phaseno).Flow.PhaseCorrPercentiles,'String',round(100*SET(gui.phaseno).Flow.PhaseCorrPercentiles));

findstatic;
recalculate_Callback;
update;


%--------------------------
function phaseedit_Callback %#ok<DEFNU>
%--------------------------
%Phase threshold editor

global SET DATA
gui = DATA.GUI.EddyCurrent;

SET(gui.phaseno).Flow.PhaseCorrPercentiles = min(1,max(0,str2num(get(gui.handles.phaseedit,'String'))/100));

%Update
set(gui.handles.phaseslider,'value',100*SET(gui.phaseno).Flow.PhaseCorrPercentiles);
set(gui.handles.phaseedit,'Value',100*SET(gui.phaseno).Flow.PhaseCorrPercentiles,'String',round(100*SET(gui.phaseno).Flow.PhaseCorrPercentiles));

findstatic;
recalculate_Callback;
update;


%------------------------------
function timedependent_Callback %#ok<DEFNU>
%------------------------------
%User clicks on checkbox time dependent 

findstatic;
timeslider_Callback;
recalculate_Callback;
update;

%--------------------------------
function [varargout] = findstatic
%--------------------------------
%Find static pixels store in gui variable

global SET DATA

gui = DATA.GUI.EddyCurrent;

myworkon;

varargout = cell(1,nargout);

%To simplify the code
no = gui.phaseno;

%Depending on method do different stuff.
switch SET(no).Flow.PhaseCorrMethod
  case 'ge'
    %Phase correction is given directly corresponding phantom
    %experiment.
    
    type = '';
    if isequal(SET(no).Flow.PhaseX,no);
      type = 'x';
    end;
    
    if isequal(SET(no).Flow.PhaseY,no);
      type = 'y';
    end;
    
    if isequal(SET(no).Flow.PhaseNo,no);
      type = 'z';
    end;
    
    if isempty(type)
      myfailed(dprintf('Could not determine direction of flow (x,y,z).'));
      myworkoff;
      return;
    end;
    
    %--- Find corresponding phase image
    point = SET(no).ImagePosition;
    match = [];
    for loop=1:length(SET)
      if not(isempty(SET(loop).Flow))
        %Flow data set at least...
        if ~isequal(SET(loop).Flow.MagnitudeNo,SET(no).Flow.MagnitudeNo)
          %They do not point to the same magnitude image i.e
          %candidates.
          
          candidate = false;
          switch type
            case 'x'
              if isequal(SET(loop).Flow.PhaseX,loop)
                candidate = true;
              end;
            case 'y'
              if isequal(SET(loop).Flow.PhaseY,loop)
                candidate = true;
              end;
            case 'z'
              if isequal(SET(loop).Flow.PhaseNo,loop)
                candidate = true;
              end;
          end;
          
          if candidate
            tempdir = cross(...
              SET(loop).ImageOrientation(1:3),...
              SET(loop).ImageOrientation(4:6));
            temppoint = SET(loop).ImagePosition;
            d = sum(tempdir(:).*temppoint(:));
            
            %Check point
            if abs(sum(point(:).*tempdir(:))-d)<1e-4
              %The point fits in this plane, add it!.
              match = [match loop]; %#ok<AGROW>
            end;
          end;
          
        end; %not same magnitude
      end; %Flow data set
    end; %Loop over image stacks
    
    if isempty(match)
      myfailed(dprintf('No matches found.'));
      SET(no).Flow.PhaseCorrMethod = 'linear';
      try
        set(gui.handles.methodlistbox,'value',1);
      catch
      end;
      myworkoff;
      return;
    end;
    
    %Create cell array
    matchcell = cell(1,length(match));
    for loop=1:length(match)
      matchcell{loop} = match(loop);
    end;
    
    if length(match)>1
      mywarning(dprintf('More than one plane candidate found.'));
      m = mymenu('Choose image stack',mathcell);
    else
      m = 1;
    end;
    
    if isequal(m,0)
      myfailed(dprintf('Aborted.'));
      myworkoff;
      return;
    end;
    
    match = match(m);
    
    gui.phasecorr = 0.5-SET(match).IM;
    
    gui.logim = SET(SET(match).Flow.MagnitudeNo).IM(:,:,1)>0.05; %ge limit
    temp = gui.phasecorr;
    %templog = repmat(~gui.logim,[1 1 SET(match).TSize]);
    %temp(templog) = 0;
    gui.phasecorr = temp;
    clear templog;
    gui.maxvel = max(temp(:));
    gui.minvel = min(temp(:));
  case {'linear','quadratic','fourth'}
    
    %--- Find static pixels
    if SET(no).Flow.PhaseCorrStaticTissueRois
      
      if ~gui.ignorerois
        %--- Use static tissue ROI's
        rois = [];
        for rloop=1:SET(gui.magno).RoiN
          if isequal(SET(gui.magno).Roi(rloop).Name,'Static tissue')
            rois = [rois rloop]; %#ok<AGROW>
          end;
        end;
      
        if isempty(rois)
          myfailed(dprintf('No static tissue regions marked.'));
          SET(gui.phaseno).Flow.PhaseCorrStaticTissueRois = false;
          try
            set(gui.handles.statictissueroischeckbox,'value',0);
          catch
          end
          return;
        end;
      
        gui.logim = false([SET(no).XSize SET(no).YSize SET(no).ZSize]);
      
        for rloop=1:length(rois)
          frames = SET(gui.magno).Roi(rois(rloop)).T(1); %handles.NO is the magnitude image stack
          slice = SET(gui.magno).Roi(rois(rloop)).Z(1);
          for tloop = 1:length(frames)
            gui.logim(:,:,slice) = gui.logim(:,:,slice) | segment('createmask',...
              [SET(no).XSize SET(no).YSize],...
              SET(gui.magno).Roi(rois(rloop)).Y(:,frames(tloop)),...
              SET(gui.magno).Roi(rois(rloop)).X(:,frames(tloop)));
          end;
        end;
      end;
      
    else
      
      %Find ROI's to exclude
      if ~gui.ignorerois        
        rois = SET(gui.magno).Roi;
        rois = rois(strcmp({rois.Name},'Non-static tissue'));
        excludeim = false([SET(no).XSize SET(no).YSize SET(no).ZSize]);
      
        for rloop=1:length(rois)
          frames = rois(rloop).T(1); %handles.NO is the magnitude image stack
          slice = rois(rloop).Z(1);
          for tloop = 1:length(frames)
            excludeim(:,:,slice) = segment('createmask',...
              [SET(no).XSize SET(no).YSize],...
              rois(rloop).Y(:,frames(tloop)),...
              rois(rloop).X(:,frames(tloop)));
          end;
        end;
      else
        excludeim = false([SET(no).XSize SET(no).YSize SET(no).ZSize]);
      end;
      
      %Calculate the temporal standard deviation
      stdmap = std(SET(no).IM,0,3);
      magmap = mean(SET(gui.magno).IM,3);
      stdmap(magmap<SET(no).Flow.PhaseCorrMagnitudeThreshold) = 0; %Exact zero
      stdmap(excludeim) = 0;
      
      %Divide into quadrants and find cut off value
      xdiv = round(SET(no).XSize/2);
      ydiv = round(SET(no).YSize/2);
      
      q1 = stdmap(1:xdiv,1:ydiv,:,:); qv1 = sort(q1(:));
      q2 = stdmap(1:xdiv,ydiv:end,:,:); qv2 = sort(q2(:));
      q3 = stdmap(xdiv:end,1:ydiv,:,:); qv3 = sort(q3(:));
      q4 = stdmap(xdiv:end,ydiv:end,:,:); qv4 = sort(q4(:));
      
      %Remove zeros
      qv1 = qv1(qv1~=0);
      qv2 = qv2(qv2~=0);
      qv3 = qv3(qv3~=0);
      qv4 = qv4(qv4~=0);
      
      %Find 10% lower percentile in each quadrant
      qq = {qv1 qv2 qv3 qv4};
      qp = zeros(1,4);
      for qi = 1:4
        qv = qq{qi};
        lfpcp = round(length(qv)*SET(no).Flow.PhaseCorrPercentiles);
        if ~isempty(qv) && lfpcp > 0
          qp(qi) = qv(lfpcp);
        else
          qp(qi) = 0;
        end
      end
      
      q1p = qp(1); q2p = qp(2); q3p = qp(3); q4p = qp(4);
      
      %Use threshold in the four different quadrants
      logim = false(size(stdmap));
      logim(1:xdiv,1:ydiv,:,:) = (q1<q1p);
      logim(1:xdiv,ydiv:end,:,:) = (q2<q2p);
      logim(xdiv:end,1:ydiv,:,:) = (q3<q3p);
      logim(xdiv:end,ydiv:end,:,:) = (q4<q4p);
      
      %Remove zeros
      logim(stdmap==0) = false;
      gui.logim = logim;
    end;
    %Recalculate
    %floweddycurrent_Callback('recalculate');
  otherwise
    myfailed(dprintf('Method %s not implemented.',SET(no).Flow.PhaseCorrMethod));
    myworkoff;
    return;
end;

%--- Exclude non static tissue
if ~gui.ignorerois

  %Find ROI's
  rois = [];
  for rloop=1:SET(gui.magno).RoiN
    if isequal(SET(gui.magno).Roi(rloop).Name,'Non-static tissue')
      rois = [rois rloop]; %#ok<AGROW>
    end;
  end;
      
  %Loop over ROI's and exclude them
  for rloop=1:length(rois)
    frames = SET(gui.magno).Roi(rois(rloop)).T(1);
    slice = SET(gui.magno).Roi(rois(rloop)).Z(1);
    for tloop = 1:length(frames)
      gui.logim(:,:,slice) = gui.logim(:,:,slice) & ~segment('createmask',...
        [SET(no).XSize SET(no).YSize],...
        SET(gui.magno).Roi(rois(rloop)).Y(:,frames(tloop)),...
        SET(gui.magno).Roi(rois(rloop)).X(:,frames(tloop)));
    end;
  end;
  
  %Later fix for multiple slices
  gui.logim = staticfix(gui.logim);
  
  %--- Include static tissue
  
  %Find ROI's
  rois = [];
  for rloop=1:SET(gui.magno).RoiN
    if isequal(SET(gui.magno).Roi(rloop).Name,'Static tissue')
      rois = [rois rloop]; %#ok<AGROW>
    end;
  end;
  
  %Loop over ROI's and exclude them
  for rloop=1:length(rois)
    frames = SET(gui.magno).Roi(rois(rloop)).T(1);
    slice = SET(gui.magno).Roi(rois(rloop)).Z(1);
    for tloop = 1:length(frames)
      gui.logim(:,:,slice) = gui.logim(:,:,slice) | segment('createmask',...
        [SET(no).XSize SET(no).YSize],...
        SET(gui.magno).Roi(rois(rloop)).Y(:,frames(tloop)),...
        SET(gui.magno).Roi(rois(rloop)).X(:,frames(tloop)));
    end;
  end;
  
  %Later fix for multiple slices
  gui.logim = staticfix(gui.logim);
end; %do not do this if ignorerois are true

%Calculate number of static pixels, this is only for research purposes and
%can later be deleted.
if nargout>1
  rows = size(gui.logim,1);
  half = round(rows/2);
  numupper = sum(sum(gui.logim(1:half,:)));
  numlower = sum(sum(gui.logim((half+1):end,:)));  
  varargout{1} = numupper;
  varargout{2} = numlower;
end;

myworkoff;

%----------------------------
function recalculate_Callback 
%----------------------------
%Recalculate, called on top level

global SET DATA

gui = DATA.GUI.EddyCurrent;

%Find method
if isequal(SET(gui.phaseno).Flow.PhaseCorrMethod,'ge')
  %Just to recalculate as we assume no bad fit  
  stop = recalculate;
  if stop
    return;
  end
else
  %First recalculation
  stop = recalculate;
  
  numtimes = 10;
  if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
    numtimes = 2;
  end;
  if stop
    return;
  end
  
  if gui.usebadfit
    for loop = 1:numtimes
      %Find pixels with bad fit
      findbadfit;
      
      %Repeat recalculation with bad fit pixels removed
      stop = recalculate;
      if stop
        return
      end
    end
  end
  
end;

%Graphical update
update;

%-------------------
function stop = recalculate
%-------------------
%Workhorse in recalculation

global SET DATA

stop = false;
gui = DATA.GUI.EddyCurrent;
no = gui.phaseno;

myworkon;

%Calculate phase correction map;
switch SET(no).Flow.PhaseCorrMethod
  case 'linear'
    order = 1;
  case 'quadratic'
    order = 2;
  case 'fourth'
    myfailed('Not yet implemented.');
    return;
  case 'ge'
    findstatic; %Check that this works
    return;
end;

%Creat grid
[x,y,z] = ndgrid(1:SET(no).XSize,1:SET(no).YSize,1:SET(no).ZSize);

%Extract static pixels
if isempty(gui.badfit)
  ind = find(gui.logim);
else
  ind = find(gui.logim & (~gui.badfit)); %Remove badfit pixels
end;

if length(ind)<50
  myfailed('Too few points, aborting.');
  stop = true;
  myworkoff;
  return;
end;

%Cut away coordinates
xi = x(:); yi=y(:); zi=z(:);
xi = xi(ind);
yi = yi(ind);
zi = zi(ind);

%Set up equation system
switch order
  case 1
    if SET(no).ZSize>1
      A = [ones(size(xi)) xi yi zi];
    else
      A = [ones(size(xi)) xi yi];
    end;
  case 2
    if SET(no).ZSize>1
      A = [ones(size(xi)) xi yi zi xi.*yi xi.*zi yi.*zi xi.^2 yi.^2 zi.^2];
    else
      A = [ones(size(xi)) xi yi xi.*yi xi.^2 yi.^2];
    end;
    %case 3
    %  if SET(handles.no).ZSize>1
    %    myfailed('Not yet implemented.');
    %    return;
    %  else
    %    A = zeros(length(xi),1+9);
    %    A(:,1) = 1;
    %    A(:,2) = xi;          %x
    %    A(:,3) = A(:,2).*xi;  %x^2
    %    A(:,4) = A(:,3).*xi;  %x^3
    %    A(:,5) = yi;          %y
    %    A(:,6) = A(:,5).*yi;  %y^2
    %    A(:,7) = A(:,6).*yi;  %y^3
    %  end;
    
end;

%Build matrix
switch order
  case 1
    if SET(no).ZSize>1
      Abuild = [ones([size(x,1)*size(x,2)*size(x,3) 1]) x(:) y(:) z(:)];
    else
      Abuild = [ones([size(x,1)*size(x,2) 1]) x(:) y(:)];
    end;
  case 2
    if SET(no).ZSize>1
      Abuild = [ones([size(x,1)*size(x,2)*size(x,3) 1]) x(:) y(:) z(:) x(:).*y(:) x(:).*z(:) y(:).*z(:) (x(:)).^2 (y(:)).^2 (z(:)).^2];
    else
      Abuild = [ones([size(x,1)*size(x,2) 1]) x(:) y(:) x(:).*y(:) (x(:)).^2 (y(:)).^2];
    end;
  case 3
    
end;

%Check if timeresolved
if SET(no).Flow.PhaseCorrTimeResolved
  timeframes = SET(no).TSize;
  handles.phasecorr = repmat(single(0),size(SET(no).IM));
else
  timeframes = 1;
end;

%--- Loop over timeframes if timeresolved other wise loop just once.
handles.maxvel = 0;
h = mywaitbarstart(timeframes,'Please wait.',1);
gui.maxvel = 0;
for tloop=1:timeframes
  
  if timeframes>1
    %Timeresolved
    phase = SET(no).IM(:,:,tloop,:)-0.5;
  else
    phase = mean(SET(no).IM(:,:,:,:),3)-0.5;
  end;
  
  %Make as vector
  phase = phase(:);
  
  %Take only stationary points
  phase = phase(ind);
  
  B = phase;
  tempstate = warning;
  warning('off'); %#ok<WNOFF>
  c = A\double(B); %Solve!
  warning(tempstate);
  
  %Build back and store
  if timeframes>1
    gui.phasecorr(:,:,tloop,:) = reshape(Abuild*c,[size(x,1) size(x,2) 1 SET(no).ZSize]);
    temp = gui.phasecorr(:,:,tloop,:);
  else
    gui.phasecorr = reshape(Abuild*c,[size(x,1) size(x,2) 1 SET(no).ZSize]);
    temp = gui.phasecorr;
  end;
  
  temp(gui.logim) = 0;
  gui.maxvel = max(gui.maxvel,max(temp(:)));
  
  h = mywaitbarupdate(h);
  
end; %End of tloop
mywaitbarclose(h);
myworkoff;

%--------------------------------
function [varargout] = findbadfit
%--------------------------------
%Find pixels with bad fit, optionally returns how many unconnected badfits
%there were.

global DATA SET

gui = DATA.GUI.EddyCurrent;

varargout = cell(1,nargout);

%Calculate phase image, take mean phase
vel = 2*SET(gui.phaseno).VENC*(mean(SET(gui.phaseno).IM,3)-0.5);

%Calculate difference
if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  diffvel = vel-mean(gui.phasecorr,3)*2*SET(gui.phaseno).VENC;
else
  diffvel = vel-gui.phasecorr*2*SET(gui.phaseno).VENC;
end;

%Take only static pixels
diffvel(~gui.logim) = 0;

%--- Do this also in quadrants (later)

%Sort the fitting error by size
diffthreshold = sort(abs(diffvel(:)));
diffthreshold = diffthreshold(diffthreshold~=0); %Remove zeros

%Remove the 25% worst fits
diffthreshold = diffthreshold(round(length(diffthreshold)*0.75));

mask = (diffvel>diffthreshold) | (diffvel<-diffthreshold);

%Store
gui.badfit = mask;

if nargout>0
  %--- warning messages if not enough (bad balanced support).
  [bwim] = bwlabel(mask,4); %4 connected regions
  varargout{1} = max(bwim(:));
end;

%----------------------
function clear_Callback %#ok<DEFNU>
%----------------------
%Clear phase background data
global SET DATA

gui = DATA.GUI.EddyCurrent;
no = gui.phaseno;

SET(no).Flow.PhaseCorr = [];

update;

%----------------------
function apply_Callback %#ok<DEFNU>
%----------------------
%Apply phase correction to data

global SET DATA

gui = DATA.GUI.EddyCurrent;
no = gui.phaseno;

SET(no).Flow.PhaseCorr = gui.phasecorr;
SET(gui.phaseno).Flow.PhaseCorrRemoveBadFit = gui.usebadfit;

update;

%---------------------
function done_Callback %#ok<DEFNU>
%---------------------
%Done button

global DATA

gui = DATA.GUI.EddyCurrent;

magno = gui.magno;
plotonok = gui.plotonok;

close(gui);
drawfunctions('drawimageno');

if isopengui('flow.fig')
  reportflow('recalculate', magno);
  return;
end;

if plotonok
  %Call plot function
  eddycheck = false;
  reportflow('init',magno,eddycheck);
end;

%----------------------------
function smallcancel_Callback %#ok<DEFNU>
%----------------------------
%Cancel button from small interface

global SET DATA

gui = DATA.GUI.EddyCurrent;

magno = gui.magno;
plotonok = gui.plotonok;
SET(magno).Flow.PhaseCorrAsk = false;

close(gui);
drawfunctions('drawimageno');

if isopengui('flow.fig')
  reportflow('recalculate', magno);
  return;
end;

if plotonok
  %Call plot function
  eddycheck = false;
  reportflow('init',magno,eddycheck);
end;


%----------------
function cancel_Callback %#ok<DEFNU>
%-----------------------
%callback for cancel Eddy current compensation
 
if ~isopengui('flow.fig')
  drawfunctions('drawimageno'); %update in main GUI
end;
close_Callback;

%----------------------
function close_Callback
%----------------------
%Close the main Eddy current compensation interface

global DATA SET

gui = DATA.GUI.EddyCurrent;
magno = gui.magno;
SET(magno).Flow.PhaseCorrAsk = false; %do not ask anymore for Eddy current compansation

try
  DATA.GUI.EddyCurrent = close(DATA.GUI.EddyCurrent);
catch   %#ok<CTCH>
  DATA.GUI.EddyCurrent=[];
  delete(gcbf);
end

  
%--------------------------------
function fixim = staticfix(logim)
%-------------------------------
%Test function to create better static regions

numslices = size(logim,3);

%Create space
fixim = false(size(logim));

%Loop over slices to remove small isolated regions
for slice = 1:numslices
  
  %Label objects in this slice
  [im] = bwlabel(logim(:,:,slice),4); %4 connected regions
  rowim = im(:);

  %Find number of objects
  maxclasses = max(rowim);
  
  %Create space for sizes
  sizevec = zeros(1,maxclasses);
  
  %Loop over objects and find number of pixels
  for loop = 1:maxclasses
    sizevec(loop) = sum(rowim==loop);
  end;
  
  %Find suitable threshold of object sizes
  maxsize = max(sizevec);
  sizethreshold = maxsize*0.1;
  
  %Merge the objects that are large enough
  mergeim = false(size(im));
  for loop = 1:maxclasses
    if sizevec(loop)>sizethreshold
      mergeim(im==loop) = true;
    end;
  end;

  %Store it
  fixim(:,:,slice) = mergeim;
end;

%----------------------------
function smalladjust_Callback %#ok<DEFNU>
%----------------------------
%Adjust in the small GUI. Take small GUI away and show big gui.

global DATA

gui = DATA.GUI.EddyCurrent;

plotonok = gui.plotonok;
close(gui);

%Start the big gui
init(plotonok);

%------------------------
function smallok_Callback %#ok<DEFNU>
%------------------------
%Apply and close

global DATA SET

gui = DATA.GUI.EddyCurrent;

%Apply
SET(gui.phaseno).Flow.PhaseCorr = gui.phasecorr;

magno = gui.magno;
plotonok = gui.plotonok;

%Close
close(gui);

drawfunctions('drawimageno');

%updated the flow interface if it is open
if isopengui('flow.fig')
  reportflow('recalculate', magno);
  return;
end;

if plotonok
  eddycheck = false;
  reportflow('init',magno,eddycheck);
end

%--------------
function update
%--------------
%Graphical update

global DATA SET

gui = DATA.GUI.EddyCurrent;

%If not big gui then call another function
if ~gui.biggui
  updatesmall;
  return;
end;

%Update timeresolved or not
if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  set([gui.handles.timetext gui.handles.timeslider],'visible','on');  
else
  set([gui.handles.timetext gui.handles.timeslider],'visible','off');
end;

%Extract maximum velocity
gui.maxvel = max(abs(gui.phasecorr(:)));

if isempty(gui.maxvel)
  gui.maxvel = 1;
end;

%--- Show velocity
if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  %timeresolved
  temp = 2*SET(gui.phaseno).VENC*(SET(gui.phaseno).IM(:,:,gui.timeframe,gui.slice)-0.5);
  temp(not(gui.logim(:,:,gui.slice))) = 0; %just show the static pixels
else
  temp = 2*SET(gui.phaseno).VENC*(mean(SET(gui.phaseno).IM(:,:,:,gui.slice),3)-0.5);
  temp(not(gui.logim(:,:,gui.slice))) = 0; %just show the static pixels
end;

if ~isempty(gui.badfit)
  temp(gui.badfit(:,:,gui.slice)) = 0; %clear also badfit
end;
  
imagesc(temp,'parent',gui.handles.phaseaxes);
axis(gui.handles.phaseaxes,'image','off');
colormap(gui.handles.phaseaxes,jet);
clim = 2*SET(gui.phaseno).VENC*[-gui.maxvel gui.maxvel];

%Check if valid range
if isequal(clim(1),0)
  clim = [-1 1];
end;
set(gui.handles.phaseaxes,'clim',clim);
%colorbar('peer',handles.phaseaxes,'south');

%--- Show calculated phaseerror
if ~isempty(gui.phasecorr)
  if gui.timeframe<size(gui.phasecorr,3)
    imagesc(gui.phasecorr(:,:,gui.timeframe,gui.slice)*2*SET(gui.phaseno).VENC,'parent',gui.handles.phaseerroraxes);
  else
    imagesc(gui.phasecorr(:,:,1,gui.slice)*2*SET(gui.phaseno).VENC,'parent',gui.handles.phaseerroraxes);
  end;
  axis(gui.handles.phaseerroraxes,'image','off');
  %clim = get(handles.phaseaxes,'clim');
  %clim = max(abs(clim));
  clim = 2*SET(gui.phaseno).VENC*[-gui.maxvel gui.maxvel];
  set(gui.handles.phaseerroraxes,'clim',clim);
  colorbar('peer',gui.handles.phaseerroraxes,'East');
else
  imagesc(zeros(SET(gui.phaseno).XSize,SET(gui.phaseno).YSize),'parent',gui.handles.phaseerroraxes);
  colorbar('peer',gui.handles.phaseerroraxes,'East');
  axis(gui.handles.phaseerroraxes,'image','off');
end;

%--- Show magnitude image
magno = SET(gui.phaseno).Flow.MagnitudeNo;
graymap = (0:DATA.GUISettings.ColorMapSize-1)'/(DATA.GUISettings.ColorMapSize-1);
if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  tempim = calcfunctions('remapuint8',SET(magno).IM(:,:,gui.timeframe,gui.slice),magno,graymap);
else
  %Not timeresolved
  tempim = calcfunctions('remapuint8',mean(SET(magno).IM(:,:,:,gui.slice),3),magno,graymap);
end;

%Fix the red image
tempr = tempim;
tempr(gui.logim(:,:,gui.slice)) = uint8(255);

%Fix the green image
if ~isempty(gui.badfit)
  tempg = tempim;
  tempg(gui.badfit(:,:,gui.slice)) = uint8(255);
else
  tempg = tempim;
end;

%Combine to an RGB image
tempim = cat(3,tempr,tempg,tempim);
image(tempim,'parent',gui.handles.previewaxes);
axis(gui.handles.previewaxes,'image','off');

%--- Show applied phase correction
if ~isempty(SET(gui.phaseno).Flow.PhaseCorr)
  if size(SET(gui.phaseno).Flow.PhaseCorr,3)>1
    temp = SET(gui.phaseno).Flow.PhaseCorr(:,:,gui.timeframe,gui.slice);
  else
    temp = SET(gui.phaseno).Flow.PhaseCorr(:,:,1,gui.slice);
  end;
else
  temp = zeros(SET(gui.phaseno).XSize,SET(gui.phaseno).YSize);
end;
temp = 2*SET(gui.phaseno).VENC*temp;
imagesc(temp,'parent',gui.handles.appliedphasecorraxes);
axis(gui.handles.appliedphasecorraxes,'image','off');
colormap(gui.handles.appliedphasecorraxes,jet);
if ~exist('clim','var')
  clim = [-gui.maxvel gui.maxvel];
end;
set(gui.handles.appliedphasecorraxes,'clim',clim);
colorbar('peer',gui.handles.appliedphasecorraxes,'East');

%Find ROI's to graphically display
nrois = [];
srois = [];
for rloop=1:SET(gui.magno).RoiN
  if isequal(SET(gui.magno).Roi(rloop).Name,'Non-static tissue')
    nrois = [nrois rloop]; %#ok<AGROW>
  end;
  if isequal(SET(gui.magno).Roi(rloop).Name,'Static tissue')
    srois = [srois rloop]; %#ok<AGROW>
  end;
end;            
      
%Loop over ROI's and exclude them
for rloop=1:length(nrois)  
  slice = SET(gui.magno).Roi(nrois(rloop)).Z(1);
  if isequal(slice,gui.slice)
    hold(gui.handles.phaseaxes,'on');
    plot(gui.handles.phaseaxes,...
      SET(gui.magno).Roi(nrois(rloop)).Y(:,1),...
      SET(gui.magno).Roi(nrois(rloop)).X(:,1),'r-');
    hold(gui.handles.phaseaxes,'off');

    hold(gui.handles.previewaxes,'on');
    plot(gui.handles.previewaxes,...
      SET(gui.magno).Roi(nrois(rloop)).Y(:,1),...
      SET(gui.magno).Roi(nrois(rloop)).X(:,1),'r-');
    hold(gui.handles.previewaxes,'off');

  end;
end;

for rloop=1:length(srois)  
  slice = SET(gui.magno).Roi(srois(rloop)).Z(1);
  if isequal(slice,gui.slice)
    hold(gui.handles.phaseaxes,'on');
    plot(gui.handles.phaseaxes,...
      SET(gui.magno).Roi(srois(rloop)).Y(:,1),...
      SET(gui.magno).Roi(srois(rloop)).X(:,1),'b-');
    hold(gui.handles.phaseaxes,'off');
  
    hold(gui.handles.previewaxes,'on');
    plot(gui.handles.previewaxes,...
      SET(gui.magno).Roi(srois(rloop)).Y(:,1),...
      SET(gui.magno).Roi(srois(rloop)).X(:,1),'b-');
    hold(gui.handles.previewaxes,'off');
  end;  
end;

%Enable/disable relevant handles.
h = [...
  gui.handles.magnitudeslider ...
  gui.handles.magnitudeedit ...
  gui.handles.phaseslider ...
  gui.handles.phaseedit];

static = SET(gui.phaseno).Flow.PhaseCorrStaticTissueRois;
if static 
  set(h,'enable','off');
else
  set(h,'enable','on');
end;


%-------------------
function updatesmall
%-------------------
%Graphical update for small GUI

global DATA SET

gui = DATA.GUI.EddyCurrent;

%Extract maximum velocity
gui.maxvel = max(abs(gui.phasecorr(:)));

if isempty(gui.maxvel)
  gui.maxvel = 1;
end;

%--- Show magnitude with overlayed velocity

%Extract
if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  %timeresolved
  vel = 2*SET(gui.phaseno).VENC*(SET(gui.phaseno).IM(:,:,gui.timeframe,gui.slice)-0.5);
  vel(not(gui.logim(:,:,gui.slice))) = 0; %just show the static pixels
else
  vel = 2*SET(gui.phaseno).VENC*(mean(SET(gui.phaseno).IM(:,:,:,gui.slice),3)-0.5);
  vel(not(gui.logim(:,:,gui.slice))) = 0; %just show the static pixels
end;

%Create magnitude image
magno = SET(gui.phaseno).Flow.MagnitudeNo;
graymap = (0:DATA.GUISettings.ColorMapSize-1)'/(DATA.GUISettings.ColorMapSize-1);
if SET(gui.phaseno).Flow.PhaseCorrTimeResolved
  im = calcfunctions('remapuint8',SET(magno).IM(:,:,gui.timeframe,gui.slice),magno,graymap);
else
  %Not timeresolved
  im = calcfunctions('remapuint8',mean(SET(magno).IM(:,:,:,gui.slice),3),magno,graymap);
end;

%Create the RGB image
rim = im;
gim = im;
bim = im;

%Create mask
if gui.usebadfit
  mask = gui.logim & ~gui.badfit;
else
  mask = gui.logim;
end
maskind = find(mask);

%Create colormap
cmap = jet(255);
cmap = uint8(255*cmap);

scale = 1/gui.maxvel;

ind = max(1,min(255,round(128+vel(mask)*scale)));
ind = ind(:);

%interpolate to get color
rim(maskind) = cmap(ind,1);
gim(maskind) = cmap(ind,2);
bim(maskind) = cmap(ind,3);

%Compose the RGB image
im = cat(3,rim,gim,bim);
image(im,'parent',gui.handles.magaxes);
axis(gui.handles.magaxes,'image','off');

%--- Applied phase correction
if ~isempty(gui.phasecorr)
  if gui.timeframe<size(gui.phasecorr,3)
    imagesc(gui.phasecorr(:,:,gui.timeframe,gui.slice)*2*SET(gui.phaseno).VENC,'parent',gui.handles.corraxes);
  else
    imagesc(gui.phasecorr(:,:,1,gui.slice)*2*SET(gui.phaseno).VENC,'parent',gui.handles.corraxes);
  end;
  axis(gui.handles.corraxes,'image','off');
  clim = 2*SET(gui.phaseno).VENC*[-gui.maxvel gui.maxvel];
  set(gui.handles.corraxes,'clim',clim);
  colorbar('peer',gui.handles.corraxes,'East');
else
  imagesc(zeros(SET(gui.phaseno).XSize,SET(gui.phaseno).YSize),'parent',gui.handles.phaseerroraxes);
  colorbar('peer',gui.handles.corraxes,'East');
  axis(gui.handles.corraxes,'image','off');
end;

%Find ROI's to graphically display
nrois = [];
srois = [];
for rloop=1:SET(gui.magno).RoiN
  if isequal(SET(gui.magno).Roi(rloop).Name,'Non-static tissue')
    nrois = [nrois rloop]; %#ok<AGROW>
  end;
  if isequal(SET(gui.magno).Roi(rloop).Name,'Static tissue')
    srois = [srois rloop]; %#ok<AGROW>
  end;
end;            
      
%Loop over ROI's and exclude them
for rloop=1:length(nrois)  
  slice = SET(gui.magno).Roi(nrois(rloop)).Z(1);
  if isequal(slice,gui.slice)
    hold(gui.handles.magaxes,'on');
    plot(gui.handles.magaxes,...
      SET(gui.magno).Roi(nrois(rloop)).Y(:,1),...
      SET(gui.magno).Roi(nrois(rloop)).X(:,1),'r-');
    hold(gui.handles.magaxes,'off');
    
  end;
end;

for rloop=1:length(srois)  
  slice = SET(gui.magno).Roi(srois(rloop)).Z(1);
  if isequal(slice,gui.slice)
    hold(gui.handles.magaxes,'on');
    plot(gui.handles.magaxes,...
      SET(gui.magno).Roi(srois(rloop)).Y(:,1),...
      SET(gui.magno).Roi(srois(rloop)).X(:,1),'b-');
    hold(gui.handles.magaxes,'off');
  
  end;  
end;

