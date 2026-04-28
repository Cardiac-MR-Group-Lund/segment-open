function varargout = segmentation(varargin)
% SEGMENTATION
% Functions for doing operations on segmentation

% Moved out from segment_main by Nisse Lundahl
%#ok<*GVMIS>
%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%----------------------------------------
function rvcopyfromlvendo_Callback
%----------------------------------------
%Copy RV from LV segmentation (endocardium).
global DATA SET NO

tools('enableundo');
if isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX = NaN(size(SET(NO).EndoX));
  SET(NO).RVEndoY = NaN(size(SET(NO).EndoY));
end
SET(NO).RVEndoX = SET(NO).EndoX;
SET(NO).RVEndoY = SET(NO).EndoY;
mymsgbox('LV endo copied to RV endo.','Done!',DATA.GUI.Segment);
drawfunctions('drawno',NO)
segment('updatevolume');


%----------------------------------------
function epicopyfromlvendo_Callback
%----------------------------------------
%Copy LV Epi from LV Endo in one slice one tf
global SET NO

tools('enableundo');
if isempty(SET(NO).EpiX)
  SET(NO).EpiX = NaN(size(SET(NO).EndoX));
  SET(NO).EpiY = NaN(size(SET(NO).EndoY));
end
SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice) = SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice);
SET(NO).EpiY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice) = SET(NO).EndoY(:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice);
lv('segmentexpandcontract_Callback',7,'epi');
logdisp('LV endo copied to LV epi.');
drawfunctions('drawno',NO)
segment('updatevolume');


%---------------------------------------
function rvcopyfromlvepi_Callback
%---------------------------------------
%Copy RV from LV segmentation (epicardium).
global DATA SET NO

tools('enableundo')
if isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX = NaN(size(SET(NO).EpiX));
  SET(NO).RVEpiY = NaN(size(SET(NO).EpiY));
end
SET(NO).RVEpiX = SET(NO).EpiX;
SET(NO).RVEpiY = SET(NO).EpiY;
mymsgbox('LV epi copied to RV epi.','Done!',DATA.GUI.Segment);

drawfunctions('drawno',NO)
segment('updatevolume');

%-------------------------------------
function ratransferfromrvendo_Callback
%-------------------------------------
%Transfer RV Endo to RA segmentation (endocardium).
type = 'RA';
transferfromrvendo_helper(type);

%-------------------------------------
function latransferfromrvendo_Callback
%-------------------------------------
%Transfer RV Endo to LA segmentation (endocardium).
type = 'LA';
transferfromrvendo_helper(type);

%----------------------------------------
function transferfromrvendo_helper(type,no)
%----------------------------------------

global SET NO
if nargin < 2
  no = NO;
end

if isempty(SET(no).RVEndoX) || isempty(SET(no).RVEndoY) || all(isnan(SET(no).RVEndoX(:))) || all(isnan(SET(no).RVEndoY(:)))
  logdisp('No RV is available');
  return
end

tools('enableundo');

opencontour = true;


%generate LA object
if isempty(SET(no).(type)) %object is empty, create a new object
  generalpen.atriumpenfunctions('createnewobject',lower(type),no);
end
numpoints = SET(no).(type).getnumpoints;
if isempty(SET(no).(type).X)
  SET(no).(type).Y = nan(numpoints,SET(no).TSize,SET(no).ZSize);
  SET(no).(type).X = nan(numpoints,SET(no).TSize,SET(no).ZSize);
end

hasfoundrv = false;

%for slices including points "TV plane", cut LA segmentation
if not(isempty(SET(no).Point.X))
  tvind = find(matches(SET(no).Point.Label,'TV plane'));
  if length(tvind) == 2 && isequal(SET(no).Point.T(tvind(1)), SET(no).Point.T(tvind(2))) && isequal(SET(no).Point.Z(tvind(1)), SET(no).Point.Z(tvind(2)))
    [rvcon, ~] = strainmitt.strainmitt('detectopenrvcontour',no,SET(no).Point.T(tvind(1)));
    [x,y] = tools('distributepointsalongcontour',rvcon(:,1),rvcon(:,2),numpoints,opencontour);
    tf = SET(no).Point.T(tvind(1));
    z = SET(no).Point.Z(tvind(1));
    SET(no).(type).X(:,tf,z) = x;
    SET(no).(type).Y(:,tf,z) = y;
    hasfoundrv = true;
  end
end

if hasfoundrv
  clearallrv_helper(no); %clear RV segmentation
  callbackfunctions('pointclearallname','TV plane',no); %remove TV-plane points

  %Update Strain flag
  if ~isempty(SET(no).StrainMitt)
    updateredo(SET(no).StrainMitt,type);
  end
  msg = sprintf('RV endo transferred to %s endo.',type);
  logdisp(msg);
  drawfunctions('drawno',no)
  segment('updatevolume');
else
  failstr = dprintf('Place two %s plane points.','TV');
  myfailed(failstr);
end

%---------------------------------
function resetedgedetection
%---------------------------------
%Reset the edge detection.
global DATA
DATA.EndoEdgeDetected = false;
DATA.EpiEdgeDetected = false;


%-----------------------------------------
function clearalllv_Callback(silent, onlyepi)
%-----------------------------------------
%Clear all lv segmentation, both endo and epi
global DATA SET NO

if nargin < 1
  silent = false;
end

if nargin < 2
  onlyepi = false;
end

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

if ~onlyepi
  SET(no).EndoX = [];
  SET(no).EndoY = [];
  SET(no).EndoDraged = false([SET(no).TSize SET(no).ZSize]);
  SET(no).EndoInterpX = [];
  SET(no).EndoInterpY = [];
end
SET(no).EpiX = [];
SET(no).EpiY = [];
SET(no).EpiDraged = false([SET(no).TSize SET(no).ZSize]);
SET(no).EpiInterpX = [];
SET(no).EpiInterpY = [];

if ~isempty(SET(no).StrainMitt)
  if ~isempty(DATA.GUI.StrainMitt)
    strainmitt.strainmitt('close_Callback');
  end
  deleteresult(SET(no).StrainMitt,'lv');
end

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  logdisp('Clear LV segmentation also clears the Strain tagging quantification');
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  if isfield(SET(no).StrainTagging,'runningregistration')
    runningregistration = SET(no).StrainTagging.runningregistration;
  else
    runningregistration=0;
  end
  if isfield(SET(no).StrainTagging,'transformparameters')
    transformparameters = SET(no).StrainTagging.transformparameters;
  end

  if isfield(SET(no).StrainTagging,'bwdtransformparameters')
    bwdtransformparameters = SET(no).StrainTagging.bwdtransformparameters;
  end

  if isfield(SET(no).StrainTagging,'mhdframe1')
    mhdframe1 = SET(no).StrainTagging.mhdframe1;
    mhdframe2 = SET(no).StrainTagging.mhdframe2;
  end

  if isfield(SET(no).StrainTagging,'mhdframe3')
    mhdframe3 = SET(no).StrainTagging.mhdframe3;
    mhdframe4 = SET(no).StrainTagging.mhdframe4;
  end

  SET(no).StrainTagging = [];
  SET(no).StrainTagging.runningregistration = runningregistration;
  if ismember('transformparameters',who)
    SET(no).StrainTagging.transformparameters = transformparameters;
  end
  if ismember('bwdtransformparameters',who)
    SET(no).StrainTagging.bwdtransformparameters = bwdtransformparameters;
  end
  if ismember('mhdframe1',who)
    SET(no).StrainTagging.mhdframe1 = mhdframe1;
    SET(no).StrainTagging.mhdframe2 = mhdframe2;
  end

  if ismember('mhdframe3',who)
    SET(no).StrainTagging.mhdframe3 = mhdframe3;
    SET(no).StrainTagging.mhdframe4 = mhdframe4;
  end

end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

if ~isempty(SET(no).Scar)
    viability('viabilityclear_Callback'); %no scar if LV is deleted
end

if ~isempty(SET(no).MaR)
    mar('clearall_Callback'); %no MaR if LV is deleted
end

%safety check so we dont go out of bounds when doing the montage LV
ind=find(findfunctions('findslicewithrvendo',no)+...
  findfunctions('findslicewithrvepi',no));
msind = find(strcmp({DATA.ViewPanelsType{DATA.ViewPanels==no}},'montagesegmented'));
if ~isempty(msind) && isempty(ind)
  DATA.ViewPanelsType{msind}='montage';
end

DATA.updatevolumeaxes
segment('updatevolume');
if ~silent
  drawfunctions('drawno',no)
end

%-----------------------------------------
function clearallrv_Callback(silent)
%-----------------------------------------
%Clear all rv segmentation, both endo and epi
arguments
  silent = false;
end
global SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

clearallrv_helper(no,silent);

%-----------------------------------------
function clearallrv_helper(no,silent)
%-----------------------------------------
%Clear all rv segmentation, both endo and epi in specified stack
arguments
  no
  silent = false;
end
global DATA SET

tools('enableundo')

SET(no).RVEndoX = [];
SET(no).RVEndoY = [];
SET(no).RVEpiX = [];
SET(no).RVEpiY = [];

SET(no).RVEndoInterpX = [];
SET(no).RVEndoInterpY = [];
SET(no).RVEpiInterpX = [];
SET(no).RVEpiInterpY = [];

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

if ~isempty(SET(no).StrainMitt)
  if ~isempty(DATA.GUI.StrainMitt)
    strainmitt.strainmitt('close_Callback');
  end
  deleteresult(SET(no).StrainMitt,'rv');
end

%safety check so we dont go out of bounds when doing the montage LV
ind = find(findfunctions('findslicewithendoall',no)+...
  findfunctions('findslicewithepiall',no), 1);

msind = find(strcmp(DATA.ViewPanelsType(DATA.ViewPanels==no),'montagesegmented'));
if ~isempty(msind) && isempty(ind)
  DATA.ViewPanelsType{msind}='montage';
end

segment('updatevolume');
if ~silent
  drawfunctions('drawno',no)
end

%---------------------------------------
function clearall_Callback(no)
%---------------------------------------
%Clear all segmentation, both endo and epi, lv and rv
global DATA SET NO

if nargin < 1
  no = NO;
end

if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

removeallinterp_Callback(true);
SET(no).EndoX = [];
SET(no).EndoY = [];
SET(no).EpiX = [];
SET(no).EpiY = [];
SET(no).RVEndoX = [];
SET(no).RVEndoY = [];
SET(no).RVEpiX = [];
SET(no).RVEpiY = [];
SET(no).EndoDraged = false([SET(no).TSize SET(no).ZSize]);
SET(no).EpiDraged = false([SET(no).TSize SET(no).ZSize]);

% SET(no).EST = 1;
% SET(no).EDT = 1;

if not(isempty(SET(no).Scar))
  viability('viabilityclear_Callback');
end

if ~isempty(SET(no).StrainMitt)
  if ~isempty(DATA.GUI.StrainMitt)
    strainmitt.strainmitt('close_Callback');
  end
  deleteresult(SET(no).StrainMitt,'lv');
  deleteresult(SET(no).StrainMitt,'rv');
end

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  logdisp('Clear LV segmentation also clears the Strain tagging quantification');
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  if isfield(SET(no).StrainTagging,'runningregistration')
    runningregistration = SET(no).StrainTagging.runningregistration;
  else
    runningregistration =  false;
  end
  if isfield(SET(no).StrainTagging,'transformparameters')
    transformparameters = SET(no).StrainTagging.transformparameters;
  end
  if isfield(SET(no).StrainTagging,'bwdtransformparameters')
    bwdtransformparameters = SET(no).StrainTagging.bwdtransformparameters;
  end
  if isfield(SET(no).StrainTagging,'mhdframe1')
    mhdframe1 = SET(no).StrainTagging.mhdframe1;
    mhdframe2 = SET(no).StrainTagging.mhdframe2;
    mhdframe3 = SET(no).StrainTagging.mhdframe3;
    mhdframe4 = SET(no).StrainTagging.mhdframe4;
  end
  SET(no).StrainTagging = [];
  SET(no).StrainTagging.runningregistration = runningregistration;
  if ismember('transformparameters',who)
    SET(no).StrainTagging.transformparameters = transformparameters;
  end
  if ismember('bwdtransformparameters',who)
    SET(no).StrainTagging.bwdtransformparameters = bwdtransformparameters;
  end
  if ismember('mhdframe1',who)
    SET(no).StrainTagging.mhdframe1 = mhdframe1;
    SET(no).StrainTagging.mhdframe2 = mhdframe2;
    SET(no).StrainTagging.mhdframe3 = mhdframe3;
    SET(no).StrainTagging.mhdframe4 = mhdframe4;
  end
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

msind = find(strcmp({DATA.ViewPanelsType{DATA.ViewPanels==no}},'montagesegmented'));
if ~isempty(msind) && isempty(ind)
  DATA.ViewPanelsType{msind}='montage';
end

DATA.updatevolumeaxes
segment('updatevolume');
drawfunctions('drawno',no)

%--------------------------------------------------------
function clearthis_Callback(endo,epi,rvendo,rvepi)
%--------------------------------------------------------
%Helper function to clear segmentation in this (current slice) slice.
global DATA SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

if nargin<4
  rvepi = true;
end

if nargin<3
  rvendo = true;
end

if nargin<2
  epi = true;
end

if nargin<1
  endo = true;
end

if DATA.ThisFrameOnly
  timeframes = SET(no).CurrentTimeFrame;
else
  timeframes = 1:SET(no).TSize;
end
ind = true(SET(no).ZSize,1);
ind(SET(no).CurrentSlice) = false;
clearslices(no,ind,timeframes,endo,epi,rvendo,rvepi);

%----------------------------------------------------------
function clearslices_Callback(endo,epi,rvendo,rvepi)
%----------------------------------------------------------
%Helper function to clear segmentation in selected slices.
global DATA SET NO

if nargin==0
  endo = true;
  epi = true;
  rvendo = true;
  rvepi = true;
end

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

%Generate index
ind = true(SET(no).ZSize,1);
ind(SET(no).StartSlice:SET(no).EndSlice)=false;

if nargin<4
  rvepi = true;
end

if nargin<3
  rvendo = true;
end

if nargin<2
  epi = true;
end

if nargin<1
  endo = true;
end

if sum(ind)==length(ind)
  myfailed('No slices selected.',DATA.GUI.Segment);
  return;
end

clearslices(no,ind,1:SET(no).TSize,endo,epi,rvendo,rvepi);

%--------------------------------------------------------------
function clearslicesthis_Callback(endo,epi,rvendo,rvepi)
%--------------------------------------------------------------
%Clear segmentation for selected slices according to mode.
global DATA SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

if nargin<4
  rvepi = true;
end

if nargin<3
  rvendo = true;
end

if nargin<2
  epi = true;
end

if nargin<1
  endo = true;
end

ind = true(SET(no).ZSize,1);
ind(SET(no).StartSlice:SET(no).EndSlice)=false;

if sum(ind)==length(ind)
  myfailed('No slices selected.',DATA.GUI.Segment);
  return;
end
% if DATA.ThisFrameOnly
tf = SET(no).CurrentTimeFrame;
% else
%   tf=1:SET(no).TSize;
% end
clearslices(no,ind,tf,endo,epi,rvendo,rvepi);

%----------------------------------------------------------------
function clearslices(no,ind,timeframes,endo,epi,rvendo,rvepi)
%----------------------------------------------------------------
%Workhorse in clearing slices.
global DATA SET

if nargin<7
  rvepi = true;
end

if nargin<6
  rvendo = true;
end

if nargin<5
  epi = true;
end

if nargin<4
  endo = true;
end

tools('enableundo')

for zloop=1:SET(no).ZSize
  if ~ind(zloop)

    %LV Endo
    if endo&&~isempty(SET(no).EndoX)
      SET(no).EndoX(:,timeframes,zloop) = NaN;
      SET(no).EndoY(:,timeframes,zloop) = NaN;
    end

    %LV Epi
    if epi&&~isempty(SET(no).EpiX)
      SET(no).EpiX(:,timeframes,zloop) = NaN;
      SET(no).EpiY(:,timeframes,zloop) = NaN;
    end

    %RV Endo
    if rvendo&&~isempty(SET(no).RVEndoX)
      SET(no).RVEndoX(:,timeframes,zloop) = NaN;
      SET(no).RVEndoY(:,timeframes,zloop) = NaN;
    end

    %RV Epi
    if rvepi&&~isempty(SET(no).RVEpiX)
      SET(no).RVEpiX(:,timeframes,zloop) = NaN;
      SET(no).RVEpiY(:,timeframes,zloop) = NaN;
    end

    %Scar
    if endo&&epi&&~isempty(SET(no).Scar)
      SET(no).Scar.Auto(:,:,zloop) = false;
      SET(no).Scar.Result(:,:,zloop) = false;
      SET(no).Scar.Manual(:,:,zloop) = int8(0);
      SET(no).Scar.NoReflow(:,:,zloop) = false;
      SET(no).Scar.MyocardMask(:,:,zloop) = false;
      viability('viabilitycalc');
    end

    %MaR
    if endo&&epi&&~isempty(SET(no).MaR)
      SET(no).MaR.Auto(:,:,timeframes,zloop) = false;
      SET(no).MaR.Result(:,:,timeframes,zloop) = false;
      SET(no).MaR.Manual(:,:,timeframes,zloop) = int8(0);
      SET(no).MaR.NoReflow(:,:,timeframes,zloop) = false;
      SET(no).MaR.MyocardMask(:,:,timeframes,zloop) = false;
      mar('update');
    end
    %Interp Points
    for tloop=1:length(timeframes)
      if endo
        if ~isempty(SET(no).EndoInterpX)
          SET(no).EndoInterpX{timeframes(tloop),zloop} = [];
          SET(no).EndoInterpY{timeframes(tloop),zloop} = [];
        end
      end
      if rvendo
        if ~isempty(SET(no).RVEndoInterpX)
          SET(no).RVEndoInterpX{timeframes(tloop),zloop} = [];
          SET(no).RVEndoInterpY{timeframes(tloop),zloop} = [];
        end
      end
      if epi
        if ~isempty(SET(no).EpiInterpX)
          SET(no).EpiInterpX{timeframes(tloop),zloop} = [];
          SET(no).EpiInterpY{timeframes(tloop),zloop} = [];
        end
      end
      if rvepi
        if ~isempty(SET(no).RVEpiInterpX)
          SET(no).RVEpiInterpX{timeframes(tloop),zloop} = [];
          SET(no).RVEpiInterpY{timeframes(tloop),zloop} = [];
        end
      end
    end
    if endo
      SET(no).EndoDraged(timeframes,zloop) = false;
    end
    if epi
      SET(no).EpiDraged(timeframes,zloop) = false;
    end
  end
end

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  logdisp('Clear LV segmentation also clears the Strain tagging quantification');
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  if isfield(SET(no).StrainTagging,'runningregistration')
    runningregistration = SET(no).StrainTagging.runningregistration;
  else
    runningregistration = 0;
  end
  if isfield(SET(no).StrainTagging,'transformparameters')
    transformparameters = SET(no).StrainTagging.transformparameters;
  end
  if isfield(SET(no).StrainTagging,'bwdtransformparameters')
    bwdtransformparameters = SET(no).StrainTagging.bwdtransformparameters;
  end
  if isfield(SET(no).StrainTagging,'mhdframe1')
    mhdframe1 = SET(no).StrainTagging.mhdframe1;
    mhdframe2 = SET(no).StrainTagging.mhdframe2;
    mhdframe3 = SET(no).StrainTagging.mhdframe3;
    mhdframe4 = SET(no).StrainTagging.mhdframe4;
  end
  SET(no).StrainTagging = [];
  SET(no).StrainTagging.runningregistration = runningregistration;
  if ismember('transformparameters',who)
    SET(no).StrainTagging.transformparameters = transformparameters;
  end
  if ismember('bwdtransformparameters',who)
    SET(no).StrainTagging.bwdtransformparameters = bwdtransformparameters;
  end
  if ismember('mhdframe1',who)
    SET(no).StrainTagging.mhdframe1 = mhdframe1;
    SET(no).StrainTagging.mhdframe2 = mhdframe2;
    SET(no).StrainTagging.mhdframe3 = mhdframe3;
    SET(no).StrainTagging.mhdframe4 = mhdframe4;
  end
end

%If straintagging initiated adjust LVupdated
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

if not(isempty(SET(no).StrainMitt)) && (ismember(1,timeframes) || ismember(SET(no).EDT,timeframes))
  logdisp('Clear segmentation also clears the Strain MITT quantification');
  if ~isempty(DATA.GUI.StrainTagging)
    strainmitt.strainmitt('close_Callback');
  end
  if rvendo
    contourtype = 'rvendo';
    deleteresult(SET(no).StrainMitt,contourtype);
  end

  if endo || epi
    contourtype = 'endo';
    deleteresult(SET(no).StrainMitt,contourtype);
  end
  strainmitt.strainmitt('deletenofromlaxgroup',no);
end

%safety check so we dont go out of bounds when doing the montage LV
ind = find(findfunctions('findslicewithendo',no,timeframes)+...
  findfunctions('findslicewithepi',no,timeframes)+...
  findfunctions('findslicewithrvendo',no,timeframes)+...
  findfunctions('findslicewithrvepi',no,timeframes));

msind = find(strcmp({DATA.ViewPanelsType{DATA.ViewPanels==no}},'montagesegmented'));
if ~isempty(msind) && isempty(ind)
  DATA.ViewPanelsType{msind} = 'montage';
end

DATA.updatevolumeaxes
segment('updatevolume');
DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawno',no)

%---------------------------
function clear_helper
%---------------------------
%Helper fcn to clear segmentation

global DATA SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

if DATA.ThisFrameOnly
  timeframes=SET(no).CurrentTimeFrame;
else
  timeframes=1:SET(no).TSize;
end
ind = true(SET(no).ZSize,1);
ind(SET(no).StartSlice:SET(no).EndSlice)=false;
clearslices(no,ind,timeframes);

%--------------------------------------------------------------------------------------------------
function [destx,desty,desttimeframes] = importcoordinatehelper(destno,desttime,sourceno,sourcex,sourcey,sourcetime,sourceslice)
%--------------------------------------------------------------------------------------------------
%Converts coordinates xsource,ysource to coordinated xdest,ydest. Assumes slices fixed separately
%
%destime = vector of destination times
%xsource = matrix (numpoints x sourceframes)
%xdest = matrix(numpoints x destframes)

if (length(desttime)==1)
  %There is just take one timeframe in destination, then take only one from source
  desttimeframes = 1;
  [pos] = calcfunctions('xyz2rlapfh',sourceno,sourcex(:,sourcetime),sourcey(:,sourcetime),repmat(sourceslice,[size(sourcex,1) 1]));
  pos = calcfunctions('rlapfh2xyz',destno,pos(:,1),pos(:,2),pos(:,3));
  destx = pos(1,:)';
  desty =pos(2,:)';
else
  if (size(sourcex,2)>1) && (sum(~isnan(sourcex(1,:)))==size(sourcex,2))
    %--- Time resolved segmentation => interpolate is method to use

    x = interp1(sourcetime,sourcex',desttime)'; %note flip it and flip it back
    y = interp1(sourcetime,sourcey',desttime)'; %note flip it and flip it back

    %initialize output
    destx = zeros(size(sourcex,1),size(x,2));
    desty = destx;

    %Loop over destination timeframes to do coordinate conversion
    for tloop=1:size(x,2)
      [pos] = calcfunctions('xyz2rlapfh',sourceno,x(:,tloop),y(:,tloop),repmat(sourceslice,size(sourcex,1),1));
      pos = calcfunctions('rlapfh2xyz',destno,pos(:,1),pos(:,2),pos(:,3));
      destx(:,tloop) = pos(1,:)';
      desty(:,tloop) =pos(2,:)';
    end
    desttimeframes = 1:size(x,2);
    %--- end of timeresolved (interpolation)
  else

    %--- Not time resolved = loop over source to find what to import
    destx = zeros(size(sourcex,1),length(desttime));
    desty = destx;
    desttimeframes = zeros(1,size(sourcex,2));
    for tloop=1:size(sourcex,2) %frames in source
      if ~isnan(sourcex(1,tloop))
        %source frame contains segmentation
        [~,nt] = min(abs(desttime-sourcetime(tloop))); %find best matching timeframe in destination
        if isnan(nt)
          nt = 1;
        end

        [pos] = calcfunctions('xyz2rlapfh',sourceno,sourcex(:,tloop),sourcey(:,tloop),repmat(sourceslice,size(sourcex,1),1));
        pos = calcfunctions('rlapfh2xyz',destno,pos(:,1),pos(:,2),pos(:,3));
        destx(:,tloop) = pos(1,:)';
        desty(:,tloop) =pos(2,:)';
        desttimeframes(tloop) = nt;
      end %source frame contains segmentation
    end %loop over souce timeframes

    ind = find(desttimeframes~=0);
    desttimeframes = desttimeframes(ind);

    if ~isempty(ind)
      destx = destx(:,ind);
      desty = desty(:,ind);
    end

  end
  %--- End not time resolved
end


%----------------------------------------------------------------------------------
function [desttimeframes,destslices,sourceslices,failedstr] = importsegmentationhelper ...
  (tono,fromno,doendo,doepi,dorvendo,dorvepi,importtf,txmapimport)
%----------------------------------------------------------------------------------
%Helper function to segmentimportsegmention Callback.
%- tono is destination of segmentation.
%- fromno is source.
%
%The function is capable of handeling slice offsets and different
%pixelssizes as well as situations when number of timeframes differ. When
%destination is not timeresolved and source is timeresolved then user is
%asked from what timeframe to take the segmentation from.

global SET

if nargin < 8
  txmapimport = false;
end

if nargin > 6
  [sourceslices,destslices,sourcetime,desttime,~,~,failedstr] = ...
    findmatchingslices(tono,fromno,doendo,doepi,dorvendo,dorvepi,false,importtf);
else
  [sourceslices,destslices,sourcetime,desttime,~,~,failedstr] = ...
    findmatchingslices(tono,fromno,doendo,doepi,dorvendo,dorvepi);
end

if isempty(sourceslices)
  desttimeframes = [];
  destslices = [];
  return;
end

%Convert desttime to destframes
desttimeframes = ones(1,length(desttime));
for tloop = 1:length(desttime)
  [~,ind] = min(abs(SET(tono).TimeVector-desttime(tloop)));
  desttimeframes(tloop) = ind;
end

%Loop over the number of slices in destination images
for zloop=1:SET(tono).ZSize

  %Match slices, see above.
  sourceslice = sourceslices(zloop);
  destslice = destslices(zloop);

  %Endo
  if doendo
    if txmapimport
      for destloop = 1:length(desttimeframes)
        if ~isempty(desttimeframes)
          [tempx,tempy,desttimeframes(destloop)] = importcoordinatehelper(tono,desttime(destloop),fromno,SET(fromno).EndoX(:,:,sourceslice),SET(fromno).EndoY(:,:,sourceslice),importtf,sourceslice);
          SET(tono).EndoX(:,destloop,destslice) = tempx;
          SET(tono).EndoY(:,destloop,destslice) = tempy;
        end
      end
    else
      [tempx,tempy,desttimeframes] = importcoordinatehelper(tono,desttime,fromno,SET(fromno).EndoX(:,:,sourceslice),SET(fromno).EndoY(:,:,sourceslice),sourcetime,sourceslice);
      if ~isempty(desttimeframes)
        SET(tono).EndoX(:,desttimeframes,destslice) = tempx;
        SET(tono).EndoY(:,desttimeframes,destslice) = tempy;
      end
    end
  end

  %Epi
  if doepi
    if txmapimport
      for destloop = 1:length(desttimeframes)
        if ~isempty(desttimeframes)
          [tempx,tempy,desttimeframes(destloop)] = importcoordinatehelper(tono,desttime(destloop),fromno,SET(fromno).EpiX(:,:,sourceslice),SET(fromno).EpiY(:,:,sourceslice),importtf,sourceslice);
          SET(tono).EpiX(:,destloop,destslice) = tempx;
          SET(tono).EpiY(:,destloop,destslice) = tempy;
        end
      end
    else
      [tempx,tempy,desttimeframes] = importcoordinatehelper(tono,desttime,fromno,SET(fromno).EpiX(:,:,sourceslice),SET(fromno).EpiY(:,:,sourceslice),sourcetime,sourceslice);
      if ~isempty(desttimeframes)
        SET(tono).EpiX(:,desttimeframes,destslice) = tempx;
        SET(tono).EpiY(:,desttimeframes,destslice) = tempy;
      end
    end
  end

  %RVEndo
  if dorvendo
    [tempx,tempy,desttimeframes] = importcoordinatehelper(tono,desttime,fromno,SET(fromno).RVEndoX(:,:,sourceslice),SET(fromno).RVEndoY(:,:,sourceslice),sourcetime,sourceslice);
    if ~isempty(desttimeframes)
      SET(tono).RVEndoX(:,desttimeframes,destslice) = tempx;
      SET(tono).RVEndoY(:,desttimeframes,destslice) = tempy;
    end
  end

  %RVEpi
  if dorvepi
    [tempx,tempy,desttimeframes] = importcoordinatehelper(tono,desttime,fromno,SET(fromno).RVEpiX(:,:,sourceslice),SET(fromno).RVEpiY(:,:,sourceslice),sourcetime,sourceslice);
    if ~isempty(desttimeframes)
      SET(tono).RVEpiX(:,desttimeframes,destslice) = tempx;
      SET(tono).RVEpiY(:,desttimeframes,destslice) = tempy;
    end
  end

end

%--------------------------------------------------------------------------
function [sourceslice,destslice,sourcetime,desttime,zdirsource,zdirdest, failedstr] = findmatchingslices ...
  (tono,fromno,doendo,doepi,dorvendo,dorvepi,takefromclosestseg,importtf)
%--------------------------------------------------------------------------
%Find matching slices between source image stack and destionation image stack
%The matching is based on camera position

global DATA SET

if nargin < 7
  takefromclosestseg = false;
end

zdirsource = cross(...
  SET(fromno).ImageOrientation(1:3),...
  SET(fromno).ImageOrientation(4:6));

zdirdest = cross(...
  SET(tono).ImageOrientation(1:3),...
  SET(tono).ImageOrientation(4:6));

anglediff = acos(zdirsource*(zdirdest'))*180/pi;
degthreshold = 10;
failedstr = '';
if anglediff > degthreshold && any([doendo,doepi,dorvendo,dorvepi])
  %myfailed('Angle difference between the two image stacks is greater than 10 degrees.',DATA.GUI.Segment);
  sourcestr = dprintf('Stack #%d',fromno);
  deststackstr = dprintf('Stack #%d',tono);
  failedstr1 = dprintf('Angle difference between source (%s) and destination (%s) is greater than %d degrees.',sourcestr,deststackstr,degthreshold);
  failedstr2 = dprintf('Import is not possible.');
  failedstr = sprintf('%s %s',failedstr1,failedstr2);
  logdisp(sprintf('Angle difference between stack %d and %d is greater than %d degrees.',sourcestr,deststackstr,degthreshold));
  sourceslice = [];
  destslice = [];
  sourcetime = [];
  desttime = [];
  zdirsource = [];
  zdirdest = [];
  %Later offer to do anyway
  return;
end

%Set up structure on positions for source slices in the
ind = (0:(SET(fromno).ZSize-1))';

%--- Loop over slices in destination image
matches = zeros(SET(tono).ZSize,1);
for loop=1:SET(tono).ZSize

  %Calculate 3D position source
  slicepos = SET(tono).ImagePosition-zdirdest*(loop-1)*(SET(tono).SliceThickness+SET(tono).SliceGap);

  %Calculate position
  pos = repmat(slicepos,SET(fromno).ZSize,1)-(repmat(SET(fromno).ImagePosition,SET(fromno).ZSize,1)-...
    repmat(zdirsource,SET(fromno).ZSize,1).*repmat(ind,1,3).*(SET(fromno).SliceThickness+SET(fromno).SliceGap));

  %Project onto zdir and take minimum
  [~,minind] = sort(abs(pos*(zdirsource')));

  if takefromclosestseg
    %take the first minind that includes segmentation.
    haveseg = false(SET(fromno).ZSize,1);
    if doendo
      haveendo = squeeze(sum(~isnan(SET(fromno).EndoX(1,:,minind)),2));
      haveseg = haveseg | haveendo>0;
    end
    if doepi
      haveepi = squeeze(sum(~isnan(SET(fromno).EpiX(1,:,minind)),2));
      haveseg = haveseg | haveepi>0;
    end
    if dorvendo
      havervendo = squeeze(sum(~isnan(SET(fromno).RVEndoX(1,:,minind)),2));
      haveseg = haveseg | havervendo>0;
    end
    if dorvepi
      havervepi = squeeze(sum(~isnan(SET(fromno).RVEpiX(1,:,minind)),2));
      haveseg = haveseg | havervepi>0;
    end
    if sum(haveseg) == 0
      myfailed('Could not find a matching slice');
      return;
    end
    matches(loop) = minind(find(haveseg,1,'first'));
  else
    matches(loop) = minind(1);
  end
end

temp = unique(matches);
if length(temp)<length(matches) && any([doendo,doepi,dorvendo,dorvepi])
  mywarning('Slices in source image stack are denser or non overlapping compared to destination slices. Result may be corrupted.',DATA.GUI.Segment);
end

%Number of timeframes
if (SET(tono).TSize==1)&&(SET(fromno).TSize>1)
  if nargin > 7
    %automatically define time frame to import from
    sourcetime = importtf;
    desttime = 1;
  else
    %Manually select time frame to import from.
    figure(16);
    set(16,'Name','Volume vs time','Numbertitle','off');
    myadjust(16,DATA.GUI.Segment);
    plot(SET(fromno).LVV,'r.-');
    grid on;

    %[~,indt] = max(SET(fromno).LVV);
    indt = round(SET(fromno).TSize*2/3); %This is approximately diastasis
    indt = max(indt,1);
    desttime = 1;
    s = myinputdlg({'Enter timeframe'},'Time frame',1,{sprintf('%d',indt)});
    close(16);
    if isempty(s)
      % user pressed cancel
      sourceslice = [];
      destslice = [];
      sourcetime = [];
      desttime = [];
      zdirsource = [];
      zdirdest = [];
      return;
    else
      [sourcetime,ok] = str2num(s{1});
      if not(ok)
        myfailed('Invalid timeframe.',DATA.GUI.Segment);
        return;
      end
      if (sourcetime<1)||(sourcetime>length(SET(fromno).LVV))
        rangestr = sprintf('[1..%d]',length(SET(fromno).LVV));
        errmsg = dprintf('Needs to be %s',rangestr);
        myfailed(errmsg,DATA.GUI.Segment);
        return;
      end
    end
  end
else
  if SET(tono).TSize==1 && SET(fromno).TSize==1
    desttime = 1;
    sourcetime = 1;
  else
    desttime = SET(tono).TimeVector;
    sourcetime = SET(fromno).TimeVector;
  end
end

%Ensure that destination data exists
if doendo
  if isempty(SET(tono).EndoX)
    numpoints = size(SET(fromno).EndoX,1);
    SET(tono).EndoX = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).EndoY = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
  end
end
if doepi
  if isempty(SET(tono).EpiX)
    numpoints = size(SET(fromno).EpiX,1);
    SET(tono).EpiX = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).EpiY = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
  end
end
if dorvendo
  if isempty(SET(tono).RVEndoX)
    numpoints = size(SET(fromno).RVEndoX,1);
    SET(tono).RVEndoX = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).RVEndoY = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
  end
end
if dorvepi
  if isempty(SET(tono).RVEpiX)
    numpoints = size(SET(fromno).RVEpiX,1);
    SET(tono).RVEpiX = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).RVEpiY = nan([numpoints SET(tono).TSize SET(tono).ZSize]);
  end
end

%Loop over the number of slices in destination images
sourceslice = nan(SET(tono).ZSize,1);
destslice = sourceslice;
for zloop = 1:SET(tono).ZSize
  %Match slices, see above.
  sourceslice(zloop) = matches(zloop);
  destslice(zloop) = zloop;
end

%----------------------------------------------
function importfromcine2scar_Callback
%----------------------------------------------
%Import segmentation from cine to scar image stack.
global DATA SET NO

if length(DATA.LVNO)==1 &&  DATA.LVNO ~= NO%(isempty(SET(DATA.LVNO).EndoX) || isempty(SET(DATA.LVNO).EpiX))%doesnt make sense to import to same
  cineshortaxisno = DATA.LVNO;
else
  cineshortaxisno = findfunctions('findcineshortaxisno');
end
if isempty(cineshortaxisno)
  myfailed('Could not find cine short-axis to import from.',DATA.GUI.Segment);
  return;
end

if isempty(SET(cineshortaxisno).EndoX) || isempty(SET(cineshortaxisno).EpiX)
  myfailed('No LV segmentation available.');
  return;
end

%Check what timeframe to take from
importtf = round(2/3*SET(cineshortaxisno).TSize); %approximately diastasis

%If no in diastasis, then try ESV (tested this seems better than EDT)
if isnan(SET(cineshortaxisno).LVV(importtf))
  importtf = SET(cineshortaxisno).EST;
end

%If no in EST try EDT
if isnan(SET(cineshortaxisno).LVV(importtf))
  importtf = SET(cineshortaxisno).EDT;
end

importsegmentationwithsnap_Callback(cineshortaxisno,importtf,[1,1,0,0]);


%----------------------------------------------
function importfromcine2txmap_Callback(doadjust)
%----------------------------------------------
%Import segmentation from cine to Tx maps image stack.
global DATA SET NO

if nargin < 1
  doadjust = false; %a flag when the icon with import and adjust postiion is clicked
end

%find cine stack to import from
if ~isempty(DATA.LVNO)
  cineshortaxisno = DATA.LVNO;
else
  cineshortaxisno = findfunctions('findcineshortaxisno');
end

%error checks
if isempty(cineshortaxisno)
  myfailed('Could not find cine short-axis to import from.',DATA.GUI.Segment);
  return;
end
if isempty(SET(cineshortaxisno).EndoX) || isempty(SET(cineshortaxisno).EpiX)
  myfailed('No LV segmentation available.');
  return;
end

%Check what timeframe to take from
importtf = round(2/3*SET(cineshortaxisno).TSize); %approximately diastasis

%If no in diastasis, then try ESV (tested this seems better than EDT)
if isnan(SET(cineshortaxisno).LVV(importtf))
  importtf = SET(cineshortaxisno).EST;
end

%If no in EST try EDT
if isnan(SET(cineshortaxisno).LVV(importtf))
  importtf = SET(cineshortaxisno).EDT;
end

%import contour
txmapimport = true;
if doadjust
  importsegmentationwithsnap_Callback(cineshortaxisno,importtf,[1,1,0,0],txmapimport);
else
  %do the just import of semgnetaion
  [destno,sourceno] = importsegmentation_Callback(cineshortaxisno,importtf,[1,1,0,0],txmapimport);
  
  if destno==0 && sourceno==0 %Added condition, to prevent bug when importsegmentation is cancelled by clicking 'Cancel' button
    return;
  end
  lvsegchanged = true; segment('updatevolume',lvsegchanged);
  drawfunctions('drawno',NO)
end

%--------------------------------------------------------------------------
function importsegmentationwithsnap_Callback(no,importtf,doseg,txmapimport)
%--------------------------------------------------------------------------
%Same as importsegmentation but also snaps contour (rigid registration)

global NO SET

switch nargin
  case 0
    [destno,sourceno] = importsegmentation_Callback;
  case 1
    [destno,sourceno] = importsegmentation_Callback(no);
  case 2
    [destno,sourceno] = importsegmentation_Callback(no,importtf);
  case 3
    [destno,sourceno] = importsegmentation_Callback(no,importtf,doseg);
  case 4
    [destno,sourceno] = importsegmentation_Callback(no,importtf,doseg,txmapimport);
end

if destno==0 && sourceno==0 %Added condition, to prevent bug when importsegmentation is cancelled by clicking 'Cancel' button
  return;
end

%Adjust to contours
if ~isequal(SET(destno).AcquisitionTime,SET(sourceno).AcquisitionTime)
  importadjust(NO);
end

lvsegchanged = true; segment('updatevolume',lvsegchanged);
drawfunctions('drawno',NO)

%---------------------------------------------------------------------------
function [destno,sourceno] = importsegmentation_Callback(no,importtf,doseg,txmapimport)
%---------------------------------------------------------------------------
%Import segmentation from another image stack.
%Imports to current image stack NO from no or if called with no input
%arguments user is asked.

global DATA SET NO

if nargin < 4
  txmapimport = false;
end

tools('enableundo');

if length(SET)<2
  myfailed('Only one image stack in memory, import from file instead (under File menu).',DATA.GUI.Segment);
  return;
end

if nargin == 0
  %Find what imagestack
  numstacks = numel(SET); 
  menuitems = cell(1,numstacks-1);
  imgplane = SET(NO).ImageViewPlane;
  imgtype = SET(NO).ImageType;
  impstacks = setdiff(1:numstacks,NO);

  if strcmpi(imgplane,'short-axis') && contains(imgtype,'T1')
    heartpart = 'LV';
    nos = findfunctions('findsaxnowithheartpart',heartpart,impstacks);
    if ~isempty(nos)
      impstacks = nos;
    end
  end

  for n = 1:numel(impstacks)
    nn = impstacks(n);
    menuitems{n} = sprintf('%d. %s',nn,[SET(nn).ImageType ' / ' SET(nn).ImageViewPlane]);
  end
  s = mymenu('Select which stack to import from',menuitems,DATA.GUI.Segment);

  if s == 0 %operation cancelled
    sourceno=0;
    destno=0;
    return;
  else
    no = impstacks(s);
  end
end

if no==NO
  myfailed('Cannot import from same image stack.',DATA.GUI.Segment);
  return;
end
if (no>length(SET))||(no<1)
  myfailed('Invalid image stack selected.',DATA.GUI.Segment);
  return;
end

sourceno = no;
destno = NO;

if nargin<3
  %Check what to do - later optionally take from input arguments.
  if ~isempty(SET(sourceno).EndoX) && ~all(isnan(SET(sourceno).EndoX(:)))
    doendo = true;
  else
    doendo = false;
  end
  if ~isempty(SET(sourceno).EpiX) && ~all(isnan(SET(sourceno).EpiX(:)))
    doepi = true;
  else
    doepi = false;
  end
  if ~isempty(SET(sourceno).RVEndoX)&& ~all(isnan(SET(sourceno).RVEndoX(:)))
    dorvendo = true;
  else
    dorvendo = false;
  end
  if ~isempty(SET(sourceno).RVEpiX)&& ~all(isnan(SET(sourceno).RVEpiX(:)))
    dorvepi = true;
  else
    dorvepi = false;
  end
else
  %assure it is a logical;
  doseg=logical(doseg);
  doendo=doseg(1);
  doepi=doseg(2);
  dorvendo=doseg(3);
  dorvepi=doseg(4);
end

if ~any([doendo, doepi, dorvendo, dorvepi])
  myfailed('No segmentation in image to import from.',DATA.GUI.Segment);
  return;
end
tools('connectinterpolation',no,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'});
if nargin > 1
  [~,~,~,failedstr] = importsegmentationhelper(destno,sourceno,doendo,doepi,dorvendo,dorvepi,importtf,txmapimport);
else
  [~,~,~,failedstr] = importsegmentationhelper(destno,sourceno,doendo,doepi,dorvendo,dorvepi);
end

if ~isempty(failedstr)
  myfailed(failedstr)
  sourceno = 0;
  destno = 0;
  return
end

if ~txmapimport
  viewfunctions('setview')
end

%-----------------------------
function e = edgehelper(im)
%-----------------------------
%Calculates an edge image

%matlab
%e = double(edge(im,'sobel'));

%matlab
%e = double(edge(im,'canny'));

%--- simple
f1 = [-1 0 1];
f2 = f1';

e1 = conv2(im,f1,'same');
e2 = conv2(im,f2,'same');

e =sqrt(e1.^2+e2.^2);

%---------------------------
function importadjust(no)
%---------------------------
%Snap segmentation to another frame

%Einar Heiberg

global SET NO

if nargin<1
  no = NO;
end

tools('enableundo');

m = 7; %Size of search region (unit is pixels)

%Check if there is segmentation
if isempty(SET(no).EndoX) && isempty(SET(no).EpiX)
  return;
end

%Loop over slices
otable = zeros(m*2+1,m*2+1,SET(no).ZSize);
h = mywaitbarstart(SET(no).ZSize,dprintf('Please wait registering contours.'));
for z = 1:SET(no).ZSize

  %Loop over timeframes and update otable
  for t = 1:SET(no).TSize

    doimage = false;
    if ~isempty(SET(no).EndoX) && ~isnan(SET(no).EndoX(1,t,z))
      doimage = true;
    end
    if ~isempty(SET(no).EpiX) && ~isnan(SET(no).EpiX(1,t,z))
      doimage = true;
    end

    if doimage
      %Extract image
      im = SET(no).IM(:,:,t,z);

      edgeim = edgehelper(im);

      %Loop over displacement
      for dx = -m:m
        for dy = -m:m

          %Add score for endo
          if ~isempty(SET(no).EndoX) && ~isnan(SET(no).EndoX(1,t,z))
            endox = SET(no).EndoX(:,t,z);
            endoy = SET(no).EndoY(:,t,z);
            otable(dx+m+1,dy+m+1,z)  = otable(dx+m+1,dy+m+1,z)+score(edgeim,endoy+dy,endox+dx);
          end

          if ~isempty(SET(no).EpiX) && ~isnan(SET(no).EpiX(1,t,z))
            epix = SET(no).EpiX(:,t,z);
            epiy = SET(no).EpiY(:,t,z);
            otable(dx+m+1,dy+m+1,z)  = otable(dx+m+1,dy+m+1,z)+score(edgeim,epiy+dy,epix+dx);
          end

        end  %dy
      end %dx
    end %doimage
  end %t

  h = mywaitbarupdate(h);
end %Z

%Loop over slices again
for z = 1:SET(no).ZSize
  %Find optimal, not elegant code, but does not take long to run
  maxvalue = 0;
  ox = 0;
  oy = 0;
  for dx = -m:m
    for dy = -m:m
      if otable(dx+m+1,dy+m+1,z)>maxvalue
        ox = dx;
        oy = dy;
        maxvalue = otable(dx+m+1,dy+m+1,z);
      end
    end
  end

  %[ox oy]

  %Update LV
  SET(no).EndoX(:,:,z) = SET(no).EndoX(:,:,z)+ox;
  SET(no).EndoY(:,:,z) = SET(no).EndoY(:,:,z)+oy;
  SET(no).EpiX(:,:,z) = SET(no).EpiX(:,:,z)+ox;
  SET(no).EpiY(:,:,z) = SET(no).EpiY(:,:,z)+oy;
end %Z

mywaitbarclose(h);

% figure(29);
% slice = SET(no).CurrentSlice;
% im = SET(no).IM(:,:,1,slice);
%
% e = edgehelper(im);
%
% endox = SET(no).EndoX(:,1,slice);
% endoy = SET(no).EndoY(:,1,slice);
% epix = SET(no).EpiX(:,1,slice);
% epiy = SET(no).EpiY(:,1,slice);
% imagesc(e);
% hold on;
% plot(endoy,endox,'w-');
% plot(epiy,epix,'w-');
% hold off;

%-------------------------------
function importadjust_Callback
%-------------------------------
%Snaps imported segmentation to image by trying to look at edge detection

%Einar Heiberg

global NO

if nargin<1
  no = NO;
end

importadjust(no);

segment('updatevolume');
drawfunctions('drawno',no)
%----------------------
function z = score(im,x,y)
%----------------------
%score for endo segmentation
z = sum(interp2(im,x,y,'linear'));

%------------------------------------------------------
function removeallinterp_Callback(silent,no,arg,indarg)
%------------------------------------------------------
%Remove all interp points.
global SET NO

if nargin == 0
  silent=false;
end
if nargin < 2 || isempty(no)
  no = NO;
end

if nargin<4
  indarg.endoind = true(1,SET(no).TSize);
  indarg.epiind = true(1,SET(no).TSize);
  indarg.rvendoind = true(1,SET(no).TSize);
  indarg.rvepiind = true(1,SET(no).TSize);
  indarg.laind = true(1,SET(no).TSize);
  indarg.raind = true(1,SET(no).TSize);
  if nargin<3
    arg.endo=true;
    arg.epi=true;
    arg.rvendo=true;
    arg.rvepi=true;
    arg.la=true;
    arg.ra=true;
  end
end

%First we assure that the sizes of interp fields are Tsize * Zsize
if size(SET(no).EndoInterpX,1)~=SET(no).TSize  && ~isempty(SET(no).EndoInterpX)
  tmp=cell(SET(no).TSize,SET(no).ZSize);
  tmpEndoInterpX=tmp;
  tmpEndoInterpY=tmp;

  for t=1:SET(no).TSize
    for z=1:SET(no).ZSize
      tmpEndoInterpX=SET(no).EndoInterpX{t,z};
      tmpEndoInterpY=SET(no).EndoInterpY{t,z};
    end
  end
  SET(no).EndoInterpX=tmpEndoInterpX;
  SET(no).EndoInterpY=tmpEndoInterpY;
end

if size(SET(no).EpiInterpX,1)~=SET(no).TSize && ~isempty(SET(no).EpiInterpX)
  tmp=cell(SET(no).TSize,SET(no).ZSize);
  tmpEpiInterpX=tmp;
  tmpEpiInterpY=tmp;

  for t=1:SET(no).TSize
    for z=1:SET(no).ZSize
      tmpEpiInterpX=SET(no).EndoInterpX{t,z};
      tmpEpiInterpY=SET(no).EndoInterpY{t,z};
    end
  end
  SET(no).EpiInterpX=tmpEpiInterpX;
  SET(no).EpiInterpY=tmpEpiInterpY;
end

if size(SET(no).RVEndoInterpX,1)~=SET(no).TSize && ~isempty(SET(no).RVEndoInterpX)
  tmp=cell(SET(no).TSize,SET(no).ZSize);
  tmpRVEndoInterpX=tmp;
  tmpRVEndoInterpY=tmp;

  for t=1:SET(no).TSize
    for z=1:SET(no).ZSize
      tmpRVEndoInterpX=SET(no).RVEndoInterpX{t,z};
      tmpRVEndoInterpY=SET(no).RVEndoInterpY{t,z};
    end
  end
  SET(no).RVEndoInterpX=tmpRVEndoInterpX;
  SET(no).RVEndoInterpY=tmpRVEndoInterpY;
end

% Then the timeframes to be removed are are removed

if arg.endo && ~isempty(SET(no).EndoInterpX)
  SET(no).EndoInterpX(indarg.endoind,:) = cell(1,1);
  SET(no).EndoInterpY(indarg.endoind,:) = cell(1,1);
end
if arg.epi && ~isempty(SET(no).EpiInterpX)
  SET(no).EpiInterpX(indarg.epiind,:) = cell(1,1);
  SET(no).EpiInterpY(indarg.epiind,:) = cell(1,1);
end

if arg.rvendo && ~isempty(SET(no).RVEndoInterpX)
  SET(no).RVEndoInterpX(indarg.rvendoind,:) = cell(1,1);
  SET(no).RVEndoInterpY(indarg.rvendoind,:) = cell(1,1);
end
if arg.rvepi && ~isempty(SET(no).RVEpiInterpX)
  SET(no).RVEpiInterpX(indarg.rvepiind,:) = cell(1,1);
  SET(no).RVEpiInterpY(indarg.rvepiind,:) = cell(1,1);
end
if arg.la && ~isempty(SET(no).LA)
  SET(no).LA.InterpX(indarg.laind,:) = cell(1,1);
  SET(no).LA.InterpY(indarg.laind,:) = cell(1,1);
end
if arg.ra && ~isempty(SET(no).RA)
  SET(no).RA.InterpX(indarg.raind,:) = cell(1,1);
  SET(no).RA.InterpY(indarg.raind,:) = cell(1,1);
end


if not(silent)
  drawfunctions('drawno',no)
end

%-------------------------------------------------
function interpolatedelineationovertime_Callback % icon callback
%-------------------------------------------------
% Interpolate LV or RV delineation over time from existing delineations

global SET NO DATA

no = NO;
Ntf = SET(no).TSize;
if Ntf  == 1
  str = dprintf('Data needs to be time resolved.');
  myfailed(str,DATA.GUI.Segment); 
  return
end

segtodo = mymenu('Select delineation to interpolate', 'LV Endocardium', 'LV Epicardium', 'RV Endocardium', 'RV Epicardium');

if isequal(segtodo,0)
  return;
end
interptypes = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
tools('connectinterpolation',no,interptypes(segtodo));
switch segtodo
  case 1 %LV endocardium
    X = SET(no).EndoX;
    Y = SET(no).EndoY;
    str = dprintf('No LV endocardium available.');
  case 2 %LV epicardium
    X = SET(no).EpiX;
    Y = SET(no).EpiY;
    str = dprintf('No LV epicardium available.');
  case 3 %RV endocardium
    X = SET(no).RVEndoX;
    Y = SET(no).RVEndoY;
    str = dprintf('No RV endocardium available.');
  case 4 %RV epicardium
    X = SET(no).RVEpiX;
    Y = SET(no).RVEpiY;
    str = dprintf('No RV epicardium available.');
end

if isempty(X) || isempty (Y) %check if there is any delineation present
  myfailed(str,DATA.GUI.Segment);
  return;
end

interpolationtype = mymenu('Linear or sinusoidal interpolation?', 'Sinusoidal', 'Linear');

selectedslices = SET(no).StartSlice:SET(no).EndSlice;

hasseg = squeeze(~isnan(X(1,:,:)));%squeeze(~isnan(X(1,:,selectedslices))); % hasseg(tf, sl)
if isrow(hasseg)
  hasseg = hasseg';
end
%tf_hasseg = sum(hasseg,2) > 0; % timeframes with segmentation

%Nsl = SET(no).ZSize;


myworkon;

%New scheme does only for selected slices % Perform interpolation for all slices
for sloop = selectedslices%1:Nsl

  tf_hasseg = hasseg(:,sloop);
  if any(tf_hasseg)
    for tloop = 1:Ntf
      % no interpolation needed if we already have one
      if tf_hasseg(tloop) %tf_hasseg(tloop)
        continue
      end
      % Find lower and upper timeframes for interpolation
      tf_hasseg_ind = find(tf_hasseg)';
      tf_hasseg_ind_periodic = [tf_hasseg_ind-Ntf tf_hasseg_ind tf_hasseg_ind+Ntf];
      tf_hasseg_ind_periodic_wrap = [tf_hasseg_ind tf_hasseg_ind tf_hasseg_ind];

      lowerindpos = find(tf_hasseg_ind_periodic < tloop, 1, 'last');
      upperindpos = find(tf_hasseg_ind_periodic > tloop, 1, 'first');

      lowerind = tf_hasseg_ind_periodic_wrap(lowerindpos);
      upperind = tf_hasseg_ind_periodic_wrap(upperindpos);

      % Interpolation weight

      iw = (tloop-tf_hasseg_ind_periodic(lowerindpos))/(tf_hasseg_ind_periodic(upperindpos)-tf_hasseg_ind_periodic(lowerindpos));
      if interpolationtype==1 %if sinusoidal interpolation is chosen
        iw = 0.5*(1-cos(pi*iw));
      end
      %
      %   New scheme does only for selected slices % Perform interpolation for all slices
      %   for sloop = selectedslices%1:Nsl
      % Convert segmentations to polar coordinates (R,Theta)
      lowerX = X(1:end-1,lowerind,sloop);
      lowerY = Y(1:end-1,lowerind,sloop);
      lower_mx = mean(lowerX);
      lower_my = mean(lowerY);
      lowerR = sqrt((lowerX - lower_mx).^2 + (lowerY - lower_my).^2);
      lowerT = atan2(lowerY-lower_my, lowerX-lower_mx);

      [~,I] = sort(lowerT);
      lowerRsort = lowerR(I);
      lowerTsort = lowerT(I);

      upperX = X(1:end-1,upperind,sloop);
      upperY = Y(1:end-1,upperind,sloop);
      upper_mx = mean(upperX);
      upper_my = mean(upperY);
      upperR = sqrt((upperX - upper_mx).^2 + (upperY - upper_my).^2);
      upperT = atan2(upperY-upper_my, upperX-upper_mx);

      [~,I] = sort(upperT);
      upperRsort = upperR(I);
      upperTsort = upperT(I);

      % Skip slices where we don't have stuff to interpolate from
      if any(isnan([lowerR' upperR' lowerT' upperT']))
        continue
      end

      % Interpolate to normalized Theta range
      Npts = size(X, 1);
      Tnorm = linspace(-pi,pi, Npts);
      Tnorm = Tnorm(1:end-1);

      x = [lowerTsort-2*pi; lowerTsort; lowerTsort+2*pi];
      y = [lowerRsort; lowerRsort; lowerRsort];
      len = sqrt(...
        conv2(x',[1 -1],'valid').^2+...
        conv2(y',[1 -1],'valid').^2);
      len = [0;len(:)]; %Add zero first
      len = cumsum(len);
      tempind = find(conv2(len,[1;-1],'valid')~=0); %Remove doublets
      %     len = [len(1);len(tempind+1)]; %used in interpolation later
      x = [x(1); x(tempind+1)];
      y = [y(1); y(tempind+1)];
      %     totallength = len(end);%used in interpolation later

      lowerRint = interp1(x,y, Tnorm, 'linear');

      x = [upperTsort-2*pi; upperTsort; upperTsort+2*pi];
      y = [upperRsort; upperRsort; upperRsort];
      len = sqrt(...
        conv2(x',[1 -1],'valid').^2+...
        conv2(y',[1 -1],'valid').^2);
      len = [0;len(:)]; %Add zero first
      len = cumsum(len);
      tempind = find(conv2(len,[1;-1],'valid')~=0); %Remove doublets
      %     len = [len(1);len(tempind+1)]; %used in interpolation later
      x = [x(1); x(tempind+1)];
      y = [y(1); y(tempind+1)];


      upperRint = interp1(x, y, Tnorm, 'linear');

      % Perform interpolation
      midR = (1-iw)*lowerRint + iw*upperRint;
      mid_mx = (1-iw)*lower_mx + iw*upper_mx;
      mid_my = (1-iw)*lower_my + iw*upper_my;

      % Go back to cartesian
      segoutX = mid_mx + midR.*cos(Tnorm);
      segoutY = mid_my + midR.*sin(Tnorm);

      % A little smoothing
      segoutX = conv([segoutX segoutX segoutX], [1 1 1]/3, 'same');
      segoutX = segoutX(Npts:(2*Npts-1));

      segoutY = conv([segoutY segoutY segoutY], [1 1 1]/3, 'same');
      segoutY = segoutY(Npts:(2*Npts-1));

      % Store as periodic
      segoutX(Npts) = segoutX(1);
      segoutY(Npts) = segoutY(1);

      % Store
      switch segtodo
        case 1 %LV endocardium
          SET(no).EndoX(:,tloop,sloop) = segoutX;
          SET(no).EndoY(:,tloop,sloop) = segoutY;
        case 2 %LV epicardium
          SET(no).EpiX(:,tloop,sloop) = segoutX;
          SET(no).EpiY(:,tloop,sloop) = segoutY;
        case 3 %RV endocardium
          SET(no).RVEndoX(:,tloop,sloop) = segoutX;
          SET(no).RVEndoY(:,tloop,sloop) = segoutY;
        case 4 %RV epicardium
          SET(no).RVEpiX(:,tloop,sloop) = segoutX;
          SET(no).RVEpiY(:,tloop,sloop) = segoutY;
      end
    end
  end
end

myworkoff;

segment('updatevolume');
drawfunctions('drawno',no)

