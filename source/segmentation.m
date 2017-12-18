function varargout = segmentation(varargin)
% SEGMENTATION
% Functions for doing operations on segmentation
%
% Moved out from segment_main by Nisse Lundahl

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%----------------------------------------
function rvcopyfromlvendo_Callback %#ok<DEFNU>
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
drawfunctions('drawallslices');
segment('updatevolume');

%---------------------------------------
function rvcopyfromlvepi_Callback %#ok<DEFNU>
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
drawfunctions('drawallslices');
segment('updatevolume');

%---------------------------------
function resetedgedetection %#ok<DEFNU>
%---------------------------------
%Reset the edge detection.
global DATA 
DATA.EndoEdgeDetected = false;
DATA.EpiEdgeDetected = false;


%-----------------------------------------
function clearalllv_Callback(silent) %#ok<DEFNU>
%-----------------------------------------
%Clear all lv segmentation, both endo and epi 
global DATA SET NO

if nargin < 1
 silent = false;
end

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

removeallpins_Callback(true); %side effect calls enableundo

SET(no).EndoX = []; 
SET(no).EndoY = []; 
SET(no).EpiX = []; 
SET(no).EpiY = []; 
SET(no).EndoDraged = false([SET(no).TSize SET(no).ZSize]);
SET(no).EpiDraged = false([SET(no).TSize SET(no).ZSize]);

SET(no).EndoInterpX = []; 
SET(no).EndoInterpY = []; 
SET(no).EpiInterpX = []; 
SET(no).EpiInterpY = []; 

% SET(no).EST = 1;
% SET(no).EDT = 1;

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  disp('Clear LV segmentation also clears the Strain tagging quantification');
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

DATA.updatevolumeaxes

segment('updatemodeldisplay');
segment('updatevolume');
if ~silent
  drawfunctions('drawimageno');
  drawfunctions('drawallslices');
end

%-----------------------------------------
function clearallrv_Callback %#ok<DEFNU>
%-----------------------------------------
%Clear all rv segmentation, both endo and epi 
global SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

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
segment('updatemodeldisplay');
segment('updatevolume');
drawfunctions('drawimageno');
drawfunctions('drawallslices');

%---------------------------------------
function clearall_Callback %#ok<DEFNU>
%---------------------------------------
%Clear all segmentation, both endo and epi, lv and rv
global DATA SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

removeallpins_Callback(true); %side effect calls enableundo
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

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  disp('Clear LV segmentation also clears the Strain tagging quantification');
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  runningregistration = SET(no).StrainTagging.runningregistration;
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

DATA.updatevolumeaxes

segment('updatemodeldisplay');
segment('updatevolume');
drawfunctions('drawimageno');
drawfunctions('drawallslices');

%--------------------------------------------------------
function clearthis_Callback(endo,epi,rvendo,rvepi) %#ok<DEFNU>
%--------------------------------------------------------
%Helper function to clear segmentation in this (current slice) slice.
global DATA SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

if nargin<4
  rvepi = true;
end;

if nargin<3
  rvendo = true;
end;

if nargin<2
  epi = true;
end;

if nargin<1
  endo = true;
end;

if DATA.ThisFrameOnly 
  timeframes = SET(no).CurrentTimeFrame;
else
  timeframes = 1:SET(no).TSize;
end;
ind = true(SET(no).ZSize,1);
ind(SET(no).CurrentSlice) = false;
clearslices(no,ind,timeframes,endo,epi,rvendo,rvepi);

%----------------------------------------------------------
function clearslices_Callback(endo,epi,rvendo,rvepi) %#ok<DEFNU>
%----------------------------------------------------------
%Helper function to clear segmentation in selected slices.
global DATA SET NO

if nargin==0
  endo = true;
  epi = true;
  rvendo = true;
  rvepi = true;    
end;

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

%Generate index
ind = true(SET(no).ZSize,1);
ind(SET(no).StartSlice:SET(no).EndSlice)=false;

if nargin<4
  rvepi = true;
end;

if nargin<3
  rvendo = true;
end;

if nargin<2
  epi = true;
end;

if nargin<1
  endo = true;
end;

if sum(ind)==length(ind)
  myfailed('No slices selected.',DATA.GUI.Segment);
  return;
end;

clearslices(no,ind,1:SET(no).TSize,endo,epi,rvendo,rvepi);

%--------------------------------------------------------------
function clearslicesthis_Callback(endo,epi,rvendo,rvepi) %#ok<DEFNU>
%--------------------------------------------------------------
%Clear segmentation for selected slices according to mode.
global DATA SET NO

no = NO;
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

if nargin<4
  rvepi = true;
end;

if nargin<3
  rvendo = true;
end;

if nargin<2
  epi = true;
end;

if nargin<1
  endo = true;
end;

ind = true(SET(no).ZSize,1);
ind(SET(no).StartSlice:SET(no).EndSlice)=false;

if sum(ind)==length(ind)
  myfailed('No slices selected.',DATA.GUI.Segment);
  return;
end;
if DATA.ThisFrameOnly
  tf=SET(no).CurrentTimeFrame;
else
  tf=1:SET(no).TSize;
end
clearslices(no,ind,tf,endo,epi,rvendo,rvepi);

%----------------------------------------------------------------
function clearslices(no,ind,timeframes,endo,epi,rvendo,rvepi)
%----------------------------------------------------------------
%Workhorse in clearing slices.
global DATA SET

if nargin<7
  rvepi = true;
end;

if nargin<6
  rvendo = true;
end;

if nargin<5
  epi = true;
end;

if nargin<4
  endo = true;
end;

tools('enableundo')

for zloop=1:SET(no).ZSize
  if ~ind(zloop)

    %LV Endo
    if endo&&~isempty(SET(no).EndoX)
      SET(no).EndoX(:,timeframes,zloop) = NaN;
      SET(no).EndoY(:,timeframes,zloop) = NaN;
    end;
    
    %LV Epi
    if epi&&~isempty(SET(no).EpiX)
      SET(no).EpiX(:,timeframes,zloop) = NaN;
      SET(no).EpiY(:,timeframes,zloop) = NaN;
    end;
    
    %RV Endo
    if rvendo&&~isempty(SET(no).RVEndoX)
      SET(no).RVEndoX(:,timeframes,zloop) = NaN;
      SET(no).RVEndoY(:,timeframes,zloop) = NaN;
    end;    
    
    %RV Epi
    if rvepi&&~isempty(SET(no).RVEpiX)
      SET(no).RVEpiX(:,timeframes,zloop) = NaN;
      SET(no).RVEpiY(:,timeframes,zloop) = NaN;
    end;    
    
    %Scar
    if endo&&epi&&~isempty(SET(no).Scar)
      SET(no).Scar.Auto(:,:,zloop) = false;
      SET(no).Scar.Result(:,:,zloop) = false;
      SET(no).Scar.Manual(:,:,zloop) = int8(0);
      SET(no).Scar.NoReflow(:,:,zloop) = false;
      SET(no).Scar.MyocardMask(:,:,zloop) = false;
      viability('viabilitycalc');
    end;
    
    %MaR
    if endo&&epi&&~isempty(SET(no).MaR)
      SET(no).MaR.Auto(:,:,timeframes,zloop) = false;
      SET(no).MaR.Result(:,:,timeframes,zloop) = false;
      SET(no).MaR.Manual(:,:,timeframes,zloop) = int8(0);
      SET(no).MaR.NoReflow(:,:,timeframes,zloop) = false;
      SET(no).MaR.MyocardMask(:,:,timeframes,zloop) = false;
      mar('update');
    end;
    
    %Pins
    for tloop=1:length(timeframes)
      if endo
        if ~isempty(SET(no).EndoPinX)
          SET(no).EndoPinX{timeframes(tloop),zloop} = [];
          SET(no).EndoPinY{timeframes(tloop),zloop} = [];
        end;
      end
      if rvendo
        if ~isempty(SET(no).RVEndoPinX)
          SET(no).RVEndoPinX{timeframes(tloop),zloop} = [];
          SET(no).RVEndoPinY{timeframes(tloop),zloop} = [];
        end;
      end;        
      if epi
        if ~isempty(SET(no).EpiPinX)        
          SET(no).EpiPinX{timeframes(tloop),zloop} = [];
          SET(no).EpiPinY{timeframes(tloop),zloop} = [];
        end;
      end
      if rvepi
        if ~isempty(SET(no).RVEpiPinX)        
          SET(no).RVEpiPinX{timeframes(tloop),zloop} = [];
          SET(no).RVEpiPinY{timeframes(tloop),zloop} = [];
        end;
      end;
    end;
    %Interp Points
    for tloop=1:length(timeframes)
      if endo
        if ~isempty(SET(no).EndoInterpX)
          SET(no).EndoInterpX{timeframes(tloop),zloop} = [];
          SET(no).EndoInterpY{timeframes(tloop),zloop} = [];
        end;
      end
      if rvendo
        if ~isempty(SET(no).RVEndoInterpX)
          SET(no).RVEndoInterpX{timeframes(tloop),zloop} = [];
          SET(no).RVEndoInterpY{timeframes(tloop),zloop} = [];
        end;
      end;        
      if epi
        if ~isempty(SET(no).EpiInterpX)        
          SET(no).EpiInterpX{timeframes(tloop),zloop} = [];
          SET(no).EpiInterpY{timeframes(tloop),zloop} = [];
        end;
      end
      if rvepi
        if ~isempty(SET(no).RVEpiInterpX)        
          SET(no).RVEpiInterpX{timeframes(tloop),zloop} = [];
          SET(no).RVEpiInterpY{timeframes(tloop),zloop} = [];
        end;
      end;
    end;
    if endo
      SET(no).EndoDraged(timeframes,zloop) = false;
    end;
    if epi
      SET(no).EpiDraged(timeframes,zloop) = false;
    end;
  end;
end;

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  disp('Clear LV segmentation also clears the Strain tagging quantification');
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

DATA.updatevolumeaxes
segment('updatevolume');
segment('updatemodeldisplay');
drawfunctions('drawimageno');
drawfunctions('drawallslices');

%---------------------------
function clear_helper %#ok<DEFNU>
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
end;
ind = true(SET(no).ZSize,1);
ind(SET(no).StartSlice:SET(no).EndSlice)=false;
clearslices(no,ind,timeframes);


%----------------------------------------------------------------------------------
function [desttimeframes,destslices,sourceslices] = importsegmentationhelper ...
  (tono,fromno,doendo,doepi,dorvendo,dorvepi,importtf)
%----------------------------------------------------------------------------------
%Helper function to segmentimportsegmention Callback. 
%- tono is destination of segmentation.
%- fromno is source.
%
%The function is capable of handeling slice offsets and different
%pixelssizes as well as situations when number of timeframes differ. When
%destination is not timeresolved and source is timeresolved then user is
%asked from what timeframe to take the segmentation from. 
%
%Work horse in importing. This function would benefit from anti-cut 
%and paste treatment.

global SET

if nargin > 6
  [sourceslices,destslices,sourcetime,desttime,zdirsource,zdirdest] = ...
    findmatchingslices(tono,fromno,doendo,doepi,dorvendo,dorvepi,false,importtf);
else
  [sourceslices,destslices,sourcetime,desttime,zdirsource,zdirdest] = ...
    findmatchingslices(tono,fromno,doendo,doepi,dorvendo,dorvepi);
end

desttimeframes = [];
if isempty(sourceslices)
  return;
end

%Loop over the number of slices in destination images
for zloop=1:SET(tono).ZSize
  
  %Match slices, see above.
  sourceslice = sourceslices(zloop);
  destslice = destslices(zloop);
  
  %Match positions
  cornerpossource = SET(fromno).ImagePosition-zdirsource*(zloop-1)*(SET(fromno).SliceThickness+SET(fromno).SliceGap);
  cornerposdest = SET(tono).ImagePosition-zdirdest*(destslice-1)*(SET(tono).SliceThickness+SET(tono).SliceGap);
  
  xofs = (cornerpossource-cornerposdest)*(SET(tono).ImageOrientation(4:6)'); %Project on one axis in mm
  yofs = (cornerpossource-cornerposdest)*(SET(tono).ImageOrientation(1:3)'); %Project on the other axis in mm

  xofs = xofs / SET(tono).ResolutionX; %Now in pixels in destination coordinates
  yofs = yofs / SET(tono).ResolutionY; %Now in pixels in destination coordinates  
  
  %factor between
  fx = SET(fromno).ResolutionX/SET(tono).ResolutionX;
  fy = SET(fromno).ResolutionY/SET(tono).ResolutionY;  
  
  %--- Store this slice
  if (length(desttime)==1)
    desttimeframes = 1;
    if doendo
      SET(tono).EndoX(:,desttime,destslice) = SET(fromno).EndoX(:,sourcetime,sourceslice)*fx+xofs;
      SET(tono).EndoY(:,desttime,destslice) = SET(fromno).EndoY(:,sourcetime,sourceslice)*fy+yofs;
    end;
    if doepi
      SET(tono).EpiX(:,desttime,destslice) = SET(fromno).EpiX(:,sourcetime,sourceslice)*fx+xofs;
      SET(tono).EpiY(:,desttime,destslice) = SET(fromno).EpiY(:,sourcetime,sourceslice)*fy+yofs;
    end;
    if dorvendo
      SET(tono).RVEndoX(:,desttime,destslice) = SET(fromno).RVEndoX(:,sourcetime,sourceslice)*fx+xofs;
      SET(tono).RVEndoY(:,desttime,destslice) = SET(fromno).RVEndoY(:,sourcetime,sourceslice)*fy+yofs;
    end;
    if dorvepi
      SET(tono).RVEpiX(:,desttime,destslice) = SET(fromno).RVEpiX(:,sourcetime,sourceslice)*fx+xofs;
      SET(tono).RVEpiY(:,desttime,destslice) = SET(fromno).RVEpiY(:,sourcetime,sourceslice)*fy+yofs;
    end;
    
  else
    %not equal, needs to interpolate
    
    %--- Do endo
    if doendo
      if (SET(fromno).TSize>1) && (sum(~isnan(SET(fromno).EndoX(1,:,sourceslice)))==SET(fromno).TSize)
        %Time resolved segmentation
        SET(tono).EndoX(:,:,destslice) = xofs+fx*interp1(...
          sourcetime,...
          SET(fromno).EndoX(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        SET(tono).EndoY(:,:,destslice) = yofs+fy*interp1(...
          sourcetime,...
          SET(fromno).EndoY(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        desttimeframes = 1:SET(tono).TSize;
      else
        %Non timeresolved segmentation
        for tloop=1:SET(fromno).TSize
          if ~isnan(SET(fromno).EndoX(1,tloop,sourceslice))
%             nt = round(1+(tloop-1)/(SET(fromno).TSize-1)*(SET(tono).TSize-1));
            [~,nt] = min(abs(desttime-sourcetime(tloop)));
            if isnan(nt)
              nt = 1;
            end;
            SET(tono).EndoX(:,nt,destslice) = xofs+fx*SET(fromno).EndoX(:,tloop,sourceslice);
            SET(tono).EndoY(:,nt,destslice) = yofs+fy*SET(fromno).EndoY(:,tloop,sourceslice);
            desttimeframes = [desttimeframes nt];
          end;
        end;
      end;
    end; %Do endo
    
    %--- Do epi
    if doepi
      if (SET(fromno).TSize>1) && (sum(~isnan(SET(fromno).EpiX(1,:,sourceslice)))==SET(fromno).TSize)
        %Time resolved segmentation
        SET(tono).EpiX(:,:,destslice) = xofs+fx*interp1(...
          sourcetime,...
          SET(fromno).EpiX(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        SET(tono).EpiY(:,:,destslice) = yofs+fy*interp1(...
          sourcetime,...
          SET(fromno).EpiY(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        desttimeframes = 1:SET(tono).TSize;
      else
        %Non timeresolved segmentation
        desttimeframes = [];
        for tloop=1:SET(fromno).TSize
          if ~isnan(SET(fromno).EpiX(1,tloop,sourceslice))
%             nt = round(1+(tloop-1)/(SET(fromno).TSize-1)*(SET(tono).TSize-1));
            [~,nt] = min(abs(desttime-sourcetime(tloop)));
            if isnan(nt)
              nt = 1;
            end;            
            SET(tono).EpiX(:,nt,destslice) = xofs+fx*SET(fromno).EpiX(:,tloop,sourceslice);
            SET(tono).EpiY(:,nt,destslice) = yofs+fy*SET(fromno).EpiY(:,tloop,sourceslice);
            desttimeframes = [desttimeframes nt];
          end;
        end;
      end;
    end; %Do epi
    
    %--- Do rvendo
    if dorvendo
      if sum(isnan(SET(fromno).RVEndoX(1,:,sourceslice)))==0
        %Time resolved segmentation
        SET(tono).RVEndoX(:,:,destslice) = xofs+fx*interp1(...
          sourcetime,...
          SET(fromno).RVEndoX(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        SET(tono).RVEndoY(:,:,destslice) = yofs+fy*interp1(...
          sourcetime,...
          SET(fromno).RVEndoY(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        desttimeframes = 1:SET(tono).TSize;
      else
        %Non timeresolved segmentation
        desttimeframes = [];
        for tloop=1:SET(fromno).TSize
          if ~isnan(SET(fromno).RVEndoX(1,tloop,sourceslice))
%             nt = round(1+(tloop-1)/(SET(fromno).TSize-1)*(SET(tono).TSize-1));
            [~,nt] = min(abs(desttime-sourcetime(tloop)));
            SET(tono).RVEndoX(:,nt,destslice) = xofs+fx*SET(fromno).RVEndoX(:,tloop,sourceslice);
            SET(tono).RVEndoY(:,nt,destslice) = yofs+fy*SET(fromno).RVEndoY(:,tloop,sourceslice);
            desttimeframes = [desttimeframes nt];
          end;
        end;
      end;
    end; %Do rvendo
        
    %--- Do rvepi    
    if dorvepi
      if sum(isnan(SET(fromno).RVEpiX(1,:,sourceslice)))==0
        %Time resolved segmentation
        SET(tono).RVEpiX(:,:,destslice) = xofs+fx*interp1(...
          sourcetime,...
          SET(fromno).RVEpiX(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        SET(tono).RVEpiY(:,:,destslice) = yofs+fy*interp1(...
          sourcetime,...
          SET(fromno).RVEpiY(:,:,sourceslice)',... %y flip it
          desttime)'; %flip it again
        desttimeframes = 1:SET(tono).TSize;
      else
        %Non timeresolved segmentation
        desttimeframes = [];
        for tloop=1:SET(fromno).TSize
          if ~isnan(SET(fromno).RVEpiX(1,tloop,sourceslice))
%             nt = round(1+(tloop-1)/(SET(fromno).TSize-1)*(SET(tono).TSize-1));
            [~,nt] = min(abs(desttime-sourcetime(tloop)));
            SET(tono).RVEpiX(:,nt,destslice) = xofs+fx*SET(fromno).RVEpiX(:,tloop,sourceslice);
            SET(tono).RVEpiY(:,nt,destslice) = yofs+fy*SET(fromno).RVEpiY(:,tloop,sourceslice);
            desttimeframes = [desttimeframes nt];
          end;
        end;
      end;
    end; %Do rvepi
    
  end; %Need to interpolate
  
end;


%--------------------------------------------------------------------------
function [sourceslice,destslice,sourcetime,desttime,zdirsource,zdirdest] = findmatchingslices ...
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

if anglediff>10
  myfailed('Angle difference between the two image stacks is greater than 10 degrees.',DATA.GUI.Segment);
  sourceslice = [];
  destslice = [];
  sourcetime = [];
  desttime = [];
  zdirsource = [];
  zdirdest = [];
  %Later offer to do anyway
  return;
end;

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
end;
  
temp = unique(matches);
if length(temp)<length(matches)
  mywarning('Slices in source image stack are denser or non overlapping compared to destination slices. Result may be corrupted.',DATA.GUI.Segment);
end;

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
    
    s = inputdlg({'Enter timeframe'},'TimeFrame',1,{sprintf('%d',indt)});
    close(16);
    if isempty(s)
      myfailed('Invalid timeframe.',DATA.GUI.Segment);
      return;
    else
      [sourcetime,ok] = str2num(s{1}); %#ok<ST2NM>
      if not(ok)
        myfailed('Invalid timeframe.',DATA.GUI.Segment);
        return;
      end;
      if (sourcetime<1)||(sourcetime>length(SET(fromno).LVV))
        myfailed(dprintf('Needs to be [1..%d',length(SET(fromno).LVV)),DATA.GUI.Segment);
        return;
      end;
    end;
  end
else
  if SET(tono).TSize==1 && SET(fromno).TSize==1
    desttime=1;
    sourcetime=1;
  elseif not(SET(fromno).Cyclic) || not(SET(tono).Cyclic) && ...
      SET(tono).HeartRate>0 && SET(fromno).HeartRate>0
    % interpolation based on heart rate
    desttime = SET(tono).TimeVector*SET(tono).HeartRate/60;
    sourcetime = SET(fromno).TimeVector*SET(fromno).HeartRate/60;
  else
    desttime = linspace(0,1,SET(tono).TSize);
    sourcetime = linspace(0,1,SET(fromno).TSize);
  end
end;

%Ensure that destination data exists
if doendo
  if isempty(SET(tono).EndoX)
    SET(tono).EndoX = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).EndoY = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
  end;    
end;
if doepi
  if isempty(SET(tono).EpiX)
    SET(tono).EpiX = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).EpiY = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
  end;      
end;
if dorvendo
  if isempty(SET(tono).RVEndoX)
    SET(tono).RVEndoX = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).RVEndoY = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
  end;
end;
if dorvepi
  if isempty(SET(tono).RVEpiX)
    SET(tono).RVEpiX = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
    SET(tono).RVEpiY = nan([DATA.NumPoints SET(tono).TSize SET(tono).ZSize]);
  end;
end;

%Loop over the number of slices in destination images
sourceslice = nan(SET(tono).ZSize,1);
destslice = sourceslice;
for zloop=1:SET(tono).ZSize  
  %Match slices, see above.
  sourceslice(zloop) = matches(zloop);
  destslice(zloop) = zloop;
end

%----------------------------------------------
function importfromcine2scar_Callback %#ok<DEFNU>
%----------------------------------------------
%Import segmentation from cine to scar image stack.
global DATA SET

cineshortaxisno = findfunctions('findcineshortaxisno');
if isempty(cineshortaxisno)
  myfailed('Could not find cine short-axis to import from.',DATA.GUI.Segment);
  return;
end
importtf = round(2/3*SET(cineshortaxisno).TSize); %approximately diastasis
importsegmentationwithsnap_Callback(cineshortaxisno,importtf);

%---------------------------------------------------------------
function importsegmentationwithsnap_Callback(no,importtf) 
%---------------------------------------------------------------
%Same as importsegmentation but also snaps contour (rigid registration)

global NO

switch nargin
  case 0
    importsegmentation_Callback;
  case 1
    importsegmentation_Callback(no);
  case 2
    importsegmentation_Callback(no,importtf);
end;

%Adjust to contours
importadjust(NO);
  
segment('updatemodeldisplay');
lvsegchanged = true; segment('updatevolume',lvsegchanged);
drawfunctions('drawallslices');

%-----------------------------------------------------
function importsegmentation_Callback(no,importtf) 
%-----------------------------------------------------
%Import segmentation from another image stack.
%Imports to current image stack NO from no or if called with no input
%arguments user is asked.

global DATA SET NO

tools('enableundo');

if length(SET)<2
  myfailed('Only one image stack in memory, import from file instead (under File menu).',DATA.GUI.Segment);
  return;
end;

if nargin==0
  %Find what imagestack
  menuitems = cell(1,numel(SET)-1);
  impstacks = setdiff(1:numel(SET),NO);
  for n = 1:numel(impstacks)
    nn = impstacks(n);
    menuitems{n} = sprintf('%d. %s',nn,[SET(nn).ImageType ' / ' SET(nn).ImageViewPlane]);
  end
  s = mymenu('Select which stack to import from',menuitems,DATA.GUI.Segment);

  if s == 0 %operation cancelled
    return;
  else
    no = impstacks(s);
  end;
end;

if no==NO
  myfailed('Cannot import from same image stack.',DATA.GUI.Segment);
  return;
end;
if (no>length(SET))||(no<1)
  myfailed('Invalid image stack selected.',DATA.GUI.Segment);
  return;
end;

%Check what to do - later optionally take from input arguments.
if ~isempty(SET(no).EndoX) && ~all(isnan(SET(no).EndoX(:)))
  doendo = true;
else
  doendo = false;
end;
if ~isempty(SET(no).EpiX) && ~all(isnan(SET(no).EpiX(:)))
  doepi = true;
else
  doepi = false;
end;
if ~isempty(SET(no).RVEndoX)&& ~all(isnan(SET(no).RVEndoX(:)))
  dorvendo = true;
else
  dorvendo = false;
end;
if ~isempty(SET(no).RVEpiX)&& ~all(isnan(SET(no).RVEpiX(:)))
  dorvepi = true;
else
  dorvepi = false;
end;

if ~any([doendo, doepi, dorvendo, dorvepi])
  myfailed('No segmentation in image to import from.',DATA.GUI.Segment);
  return;
end

if nargin > 1
  [~,~,sourceslices]=importsegmentationhelper(NO,no,doendo,doepi,dorvendo,dorvepi,importtf);
else
  [~,~,sourceslices]=importsegmentationhelper(NO,no,doendo,doepi,dorvendo,dorvepi);
end

%---------------------------
function importadjust(no)
%---------------------------
%Snap segmentation to another frame 

%Einar Heiberg

global SET NO

if nargin<1
  no = NO;
end;

tools('enableundo');

m = 7; %Size of search region

%Check if there is segmentation
if isempty(SET(no).EndoX) && isempty(SET(no).EpiX)
  return;
end;

%Loop over slices
otable = zeros(m*2+1,m*2+1,SET(no).ZSize);   
h = waitbar(0,'Please wait registering contours.');
for z = 1:SET(no).ZSize          
  
  %Loop over timeframes and update otable
  for t = 1:SET(no).TSize
    
    doimage = false;
    if ~isempty(SET(no).EndoX) && ~isnan(SET(no).EndoX(1,t,z))
      doimage = true;
    end;
    if ~isempty(SET(no).EpiX) && ~isnan(SET(no).EpiX(1,t,z))
      doimage = true;
    end;
    
    if doimage
      %Extract image
      im = SET(no).IM(:,:,t,z);
      e = double(edge(im,'sobel'));      
      %e = double(edge(im,'canny'));
      
      %Loop over displacement
      for dx = -m:m
        for dy = -m:m
          
          %Add score for endo
          if ~isempty(SET(no).EndoX) && ~isnan(SET(no).EndoX(1,t,z))
            endox = SET(no).EndoX(:,t,z);
            endoy = SET(no).EndoY(:,t,z);
            otable(dx+m+1,dy+m+1,z)  = otable(dx+m+1,dy+m+1,z)+score(e,endoy+dy,endox+dx);
          end;
          
          if ~isempty(SET(no).EndoX) && ~isnan(SET(no).EpiX(1,t,z))
            epix = SET(no).EpiX(:,t,z);
            epiy = SET(no).EpiY(:,t,z);
            otable(dx+m+1,dy+m+1,z)  = otable(dx+m+1,dy+m+1,z)+score(e,epiy+dy,epix+dx);
          end;
          
        end;  %dy
      end; %dx
    end; %doimage
  end; %t
    
  waitbar(z/SET(no).ZSize,h);
end; %Z

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
      end;
    end;
  end;
  
  %[ox oy]
  
  %Update LV
  SET(no).EndoX(:,:,z) = SET(no).EndoX(:,:,z)+ox;
  SET(no).EndoY(:,:,z) = SET(no).EndoY(:,:,z)+oy;
  SET(no).EpiX(:,:,z) = SET(no).EpiX(:,:,z)+ox;
  SET(no).EpiY(:,:,z) = SET(no).EpiY(:,:,z)+oy;
end; %Z

close(h);

%-------------------------------
function importadjust_Callback
%-------------------------------
%Snaps imported segmentation to image by trying to look at edge detection

%Einar Heiberg

global NO

if nargin<1
  no = NO;
end;

importadjust(no);

segment('updatemodeldisplay');
segment('updatevolume');
drawfunctions('updatenopanels',no);

%----------------------
function z = score(im,x,y)
%----------------------
z = sum(interp2(im,x,y,'linear'));

%--------------------------------------
function removepin_Callback(type,m) %#ok<INUSD,DEFNU>
%--------------------------------------
%Removes pins
global SET NO

tools('enableundo');
if nargin~=1
  mydisp('Remove pin by right mouse click on the pin.');
end;

[x,y,slice] = segment('getclickedcoords');

%Find what pin.
switch type
  case  'endo'
    if isempty(SET(NO).EndoPinX)
      return;
    end;
    xm = SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice};
    ym = SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice};
  case 'epi'
    if isempty(SET(NO).EpiPinX)
      return;
    end;
    xm = SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice};
    ym = SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice};
  case 'rvendo'
    if isempty(SET(NO).RVEndoPinX)
      return;
    end;
    xm = SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice};
    ym = SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice};    
end;

%Calculate distance
[~,tempind] = min(sqrt((xm-y).^2+(ym-x).^2));

%remove from endo and put back
ind = true(length(xm),1);
ind(tempind)=false;

switch type
  case 'endo'
    SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice} = xm(ind);
    SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice} = ym(ind);
  case 'epi'
    SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice} = xm(ind);
    SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice} = ym(ind);
  case 'rvendo'
    SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice} = xm(ind);
    SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice} = ym(ind);    
end;

segment('updatemodeldisplay');
drawfunctions('drawsliceno');


%--------------------------------------
function removethispins_Callback(endo,epi,rvendo,rvepi) %#ok<DEFNU>
%--------------------------------------
%Remove pins in this slice and timeframe
global SET NO


if nargin<4
  endo=true;
  epi=true;
  rvendo=true;
  rvepi=true;
end

tools('enableundo');
[~,~,slice] = segment('getclickedcoords');

if endo && ~isempty(SET(NO).EndoPinX)
  SET(NO).EndoPinX{SET(NO).CurrentTimeFrame,slice} = [];
  SET(NO).EndoPinY{SET(NO).CurrentTimeFrame,slice} = [];
end;

if epi && ~isempty(SET(NO).EpiPinX)
  SET(NO).EpiPinX{SET(NO).CurrentTimeFrame,slice} = [];
  SET(NO).EpiPinY{SET(NO).CurrentTimeFrame,slice} = [];
end;

if rvendo && ~isempty(SET(NO).RVEndoPinX)
  SET(NO).RVEndoPinX{SET(NO).CurrentTimeFrame,slice} = [];
  SET(NO).RVEndoPinY{SET(NO).CurrentTimeFrame,slice} = [];
end;

if rvepi && ~isempty(SET(NO).RVEpiPinX)
  SET(NO).RVEpiPinX{SET(NO).CurrentTimeFrame,slice} = [];
  SET(NO).RVEpiPinY{SET(NO).CurrentTimeFrame,slice} = [];
end;

segment('updatemodeldisplay');
drawfunctions('drawsliceno');

%-----------------------------------------
function removeallpins_Callback(silent,endo,epi,rvendo,rvepi)
%-----------------------------------------
%Remove all pins
global SET NO

if nargin == 0
  silent=false;
end

if nargin<5
  endo=true;
  epi=true;
  rvendo=true;
  rvepi=true;
end

tools('enableundo');

if endo && ~isempty(SET(NO).EndoPinX)
  SET(NO).EndoPinX = [];
  SET(NO).EndoPinY = [];
end
if epi && ~isempty(SET(NO).EpiPinX)
  SET(NO).EpiPinX = [];
  SET(NO).EpiPinY = [];
end

if rvendo && ~isempty(SET(NO).RVEndoPinX)
  SET(NO).RVEndoPinX = [];
  SET(NO).RVEndoPinY = [];
end;
if rvepi && ~isempty(SET(NO).RVEpiPinX)
  SET(NO).RVEpiPinX = [];
  SET(NO).RVEpiPinY = [];
end;  

if not(silent)
  segment('updatemodeldisplay');
  drawfunctions('drawallslices');
end;

%---------------------------------------------------
function removeallpinsthisslice_Callback(current,endo,epi,rvendo,rvepi) %#ok<DEFNU>
%---------------------------------------------------
%Remove all pins in this slice
global SET NO

tools('enableundo');

if (nargin>0) && (current)
  slice = SET(NO).CurrentSlice;
else
  % never happens, or?? /JU
  [~,~,slice] = segment('getclickedcoords');
end;

if nargin<5
  endo=true;
  epi=true;
  rvendo=true;
  rvepi=true;
end

for tloop=1:SET(NO).TSize
  if endo && ~isempty(SET(NO).EndoPinX)
    SET(NO).EndoPinX{tloop,slice} = [];
    SET(NO).EndoPinY{tloop,slice} = [];
  end;
  if epi && ~isempty(SET(NO).EpiPinX)  
    SET(NO).EpiPinX{tloop,slice} = [];
    SET(NO).EpiPinY{tloop,slice} = [];
  end;
  if rvendo && ~isempty(SET(NO).RVEndoPinX)
    SET(NO).RVEndoPinX{tloop,slice} = [];
    SET(NO).RVEndoPinY{tloop,slice} = [];
  end;
  if rvepi && ~isempty(SET(NO).RVEpiPinX)  
    SET(NO).RVEpiPinX{tloop,slice} = [];
    SET(NO).RVEpiPinY{tloop,slice} = [];
  end;  
end;

segment('updatemodeldisplay');
drawfunctions('drawsliceno');

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
  if nargin<3
    arg.endo=true;
    arg.epi=true;
    arg.rvendo=true;
    arg.rvepi=true;
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
end;
if arg.rvepi && ~isempty(SET(no).RVEpiInterpX)
  SET(no).RVEpiInterpX(indarg.rvepiind,:) = cell(1,1);
  SET(no).RVEpiInterpY(indarg.rvepiind,:) = cell(1,1);
end;  

if not(silent)
  segment('updatemodeldisplay');
  drawfunctions('drawallslices');
end;

%-------------------------------------------------
function interpolatedelineationovertime_Callback
%-------------------------------------------------
% Interpolate LV or RV delineation over time from existing delineations

global SET NO DATA

no = NO;


segtodo=mymenu('Select delineation to interpolate', 'LV endocardium', 'LV epicardium', 'RV endocardium', 'RV epicardium');

if isequal(segtodo,0)
  myfailed('Aborted.');
  return;
end;

switch segtodo
  case 1 %LV endocardium
    X=SET(no).EndoX;
    Y=SET(no).EndoY;
    str='No LV endocardium to interpolate.'; 
  case 2 %LV epicardium
    X=SET(no).EpiX;
    Y=SET(no).EpiY;
    str = 'No LV epicardium to interpolate.'; 
  case 3 %RV endocardium
    X=SET(no).RVEndoX;
    Y=SET(no).RVEndoY;
    str= 'No RV endocardium to interpolate.'; 
  case 4 %RV epicardium
    X=SET(no).RVEpiX;
    Y=SET(no).RVEpiY;
    str = 'No RV epicardium to interpolate.';
end

if isempty(X) || isempty (Y) %check if there is any delineation present
  myfailed(str,DATA.GUI.Segment);
  return;
end

interpolationtype=mymenu('Linear or sinusoidal interpolation?', 'Sinusoidal', 'Linear')

selectedslices=SET(no).StartSlice:SET(no).EndSlice;

hasseg = squeeze(~isnan(X(1,:,:)));%squeeze(~isnan(X(1,:,selectedslices))); % hasseg(tf, sl)
if isrow(hasseg)
  hasseg=hasseg';
end
%tf_hasseg = sum(hasseg,2) > 0; % timeframes with segmentation

%Nsl = SET(no).ZSize;
Ntf = SET(no).TSize;

myworkon;

%New scheme does only for selected slices % Perform interpolation for all slices
for sloop = selectedslices%1:Nsl
  
tf_hasseg=hasseg(:,sloop);
  for tloop = 1:Ntf
  % no interpolation needed if we already have one  
  if tf_hasseg(tloop);%tf_hasseg(tloop)
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
    
    lowerRint = interp1([lowerTsort-2*pi; lowerTsort; lowerTsort+2*pi], ...
      [lowerRsort; lowerRsort; lowerRsort], Tnorm, 'linear');
    upperRint = interp1([upperTsort-2*pi; upperTsort; upperTsort+2*pi], ...
      [upperRsort; upperRsort; upperRsort], Tnorm, 'linear');
    
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



% for sloop = selectedslices%1:Nsl
% for tloop = 1:Ntf
%   % no interpolation needed if we already have one  
%   if hasseg(tloop,sloop-SET(no).StartSlice+1);%tf_hasseg(tloop)
%     continue
%   end
% 
%   % Find lower and upper timeframes for interpolation
%   tf_hasseg_ind = find(tf_hasseg)';
%   tf_hasseg_ind_periodic = [tf_hasseg_ind-Ntf tf_hasseg_ind tf_hasseg_ind+Ntf];
%   tf_hasseg_ind_periodic_wrap = [tf_hasseg_ind tf_hasseg_ind tf_hasseg_ind];
% 
%   lowerindpos = find(tf_hasseg_ind_periodic < tloop, 1, 'last');
%   upperindpos = find(tf_hasseg_ind_periodic > tloop, 1, 'first');
% 
%   lowerind = tf_hasseg_ind_periodic_wrap(lowerindpos);
%   upperind = tf_hasseg_ind_periodic_wrap(upperindpos);
% 
%   % Interpolation weight
%   
%   iw = (tloop-tf_hasseg_ind_periodic(lowerindpos))/(tf_hasseg_ind_periodic(upperindpos)-tf_hasseg_ind_periodic(lowerindpos));
%   if interpolationtype==1 %if sinusoidal interpolation is chosen
%     iw = 0.5*(1-cos(pi*iw));
%   end
% %   
% %   New scheme does only for selected slices % Perform interpolation for all slices
% %   for sloop = selectedslices%1:Nsl
%     % Convert segmentations to polar coordinates (R,Theta)
%     lowerX = X(1:end-1,lowerind,sloop);
%     lowerY = Y(1:end-1,lowerind,sloop);
%     lower_mx = mean(lowerX);
%     lower_my = mean(lowerY);
%     lowerR = sqrt((lowerX - lower_mx).^2 + (lowerY - lower_my).^2);
%     lowerT = atan2(lowerY-lower_my, lowerX-lower_mx);
%     
%     [~,I] = sort(lowerT);
%     lowerRsort = lowerR(I);
%     lowerTsort = lowerT(I);
% 
%     upperX = X(1:end-1,upperind,sloop);
%     upperY = Y(1:end-1,upperind,sloop);
%     upper_mx = mean(upperX);
%     upper_my = mean(upperY);
%     upperR = sqrt((upperX - upper_mx).^2 + (upperY - upper_my).^2);
%     upperT = atan2(upperY-upper_my, upperX-upper_mx);
%     
%     [~,I] = sort(upperT);
%     upperRsort = upperR(I);
%     upperTsort = upperT(I);
%     
%     % Skip slices where we don't have stuff to interpolate from
%     if any(isnan([lowerR' upperR' lowerT' upperT']))
%       continue
%     end
%     
%     % Interpolate to normalized Theta range
%     Npts = size(X, 1);
%     Tnorm = linspace(-pi,pi, Npts);
%     Tnorm = Tnorm(1:end-1);
%     
%     lowerRint = interp1([lowerTsort-2*pi; lowerTsort; lowerTsort+2*pi], ...
%       [lowerRsort; lowerRsort; lowerRsort], Tnorm, 'linear');
%     upperRint = interp1([upperTsort-2*pi; upperTsort; upperTsort+2*pi], ...
%       [upperRsort; upperRsort; upperRsort], Tnorm, 'linear');
%     
%     % Perform interpolation
%     midR = (1-iw)*lowerRint + iw*upperRint;
%     mid_mx = (1-iw)*lower_mx + iw*upper_mx;
%     mid_my = (1-iw)*lower_my + iw*upper_my;
% 
%     % Go back to cartesian
%     segoutX = mid_mx + midR.*cos(Tnorm);
%     segoutY = mid_my + midR.*sin(Tnorm);
% 
%     % A little smoothing
%     segoutX = conv([segoutX segoutX segoutX], [1 1 1]/3, 'same');
%     segoutX = segoutX(Npts:(2*Npts-1));
% 
%     segoutY = conv([segoutY segoutY segoutY], [1 1 1]/3, 'same');
%     segoutY = segoutY(Npts:(2*Npts-1));
% 
%     % Store as periodic
%     segoutX(Npts) = segoutX(1);
%     segoutY(Npts) = segoutY(1);
% 
%     % Store
% 
%     switch segtodo
%       case 1 %LV endocardium
%         SET(no).EndoX(:,tloop,sloop) = segoutX;
%         SET(no).EndoY(:,tloop,sloop) = segoutY;
%       case 2 %LV epicardium
%         SET(no).EpiX(:,tloop,sloop) = segoutX;
%         SET(no).EpiY(:,tloop,sloop) = segoutY;
%       case 3 %RV endocardium
%         SET(no).RVEndoX(:,tloop,sloop) = segoutX;
%         SET(no).RVEndoY(:,tloop,sloop) = segoutY;
%       case 4 %RV epicardium
%         SET(no).RVEpiX(:,tloop,sloop) = segoutX;
%         SET(no).RVEpiY(:,tloop,sloop) = segoutY;
%     end
%     
%     
%   end
% end

myworkoff;

segment('updatevolume');
segment('updatemodeldisplay');
drawfunctions('drawimageno');
drawfunctions('drawallslices');

