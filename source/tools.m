function [varargout] = tools(varargin)
%This file contains functions for image tools.
%Examples are imageresampling, rotation, manipulation, undo tools, copying
%tools.

%Einar Heiberg

%#ok<*GVMIS>

%TODO:
%- make so that flow stacks are correctly flipped and avoid usageof NO
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------------------------------------------
function setcurrenttimeframeasfirst_Callback(applyto,showmessage)
%------------------------------------------------
%Sets current timeframe as first time frame by cycling data in time.
%Works with existing contours, but not general segmentation module, or
%time resolved annotation points.
arguments
  applyto = 'selected'
  showmessage = true
end

global DATA SET NO

if SET(NO).TSize < 2
  myfailed('Data needs to be time resolved.',DATA.GUI.Segment);
  return
end

if showmessage
  str1 = dprintf('Setting current timeframe as first timeframe is irreversible.');
  str2 = dprintf('Are you sure?');
  str = sprintf('%s %s',str1,str2);
  if ~yesno(str,[],DATA.GUI.Segment)
    return
  end
end
if ~isempty(SET(NO).StrainTagging)
  if ~yesno('Reorder time frames clears the Strain tagging analysis. Are you sure?',[],DATA.GUI.Segment)
    return
  end
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  SET(NO).StrainTagging = [];
end

if ~isempty(SET(NO).StrainMitt)
  logdisp(sprintf('Reordering time frames clears the Strain MITT analysis in stack %d.',NO));
  if ~isempty(DATA.GUI.StrainMitt)
    strainmitt.strainmitt('close_Callback');
  end
  SET(NO).StrainMitt = [];
end

% if user still sure to do that, then proceed

switch applyto
  case 'selected'
    slices = SET(NO).StartSlice:SET(NO).EndSlice;
  case 'all'
    slices = 1:SET(NO).ZSize;
end

shiftvec1 = zeros(1,2);
shiftvec2 = zeros(1,2);
shiftvec3 = zeros(1,3);
shiftvec3(3) = -(SET(NO).CurrentTimeFrame-1); %1 gives 0, 2 gives -1,...
shiftvec2(2) = shiftvec3(3);
shiftvec1(1) = shiftvec3(3);

shiftby = shiftvec1(1);

%checks if flow data
nos = SET(NO).Linked;

for no = nos
  for zloop = slices
    %Image
    SET(no).IM(:,:,:,zloop) = circshift(SET(no).IM(:,:,:,zloop),shiftvec3);

    %Segmentation
    types = {'Endo','Epi','RVEndo','RVEpi'};
    for loop = 1:length(types)
      if ~isempty(SET(no).([types{loop},'X']))
        SET(no).([types{loop},'X'])(:,:,zloop) = circshift(SET(no).([types{loop},'X'])(:,:,zloop),shiftvec2);
        SET(no).([types{loop},'Y'])(:,:,zloop) = circshift(SET(no).([types{loop},'Y'])(:,:,zloop),shiftvec2);
      end
    end
    types = {'LA','RA'};
    for loop = 1:length(types)
      if (~isempty(SET(no).(types{loop})) && ~isempty(SET(no).(types{loop}).X))
        SET(no).(types{loop}).X(:,:,zloop) = circshift(SET(no).(types{loop}).X(:,:,zloop),shiftvec2);
        SET(no).(types{loop}).Y(:,:,zloop) = circshift(SET(no).(types{loop}).Y(:,:,zloop),shiftvec2);
      end
    end

    %Pins
    types = {'EndoPin','EpiPin','RVEndoPin','RVEpiPin'};
    for loop = 1:length(types)
      if ~isempty(SET(no).([types{loop},'X']))
        SET(no).([types{loop},'X'])(:,zloop) = circshift(SET(no).([types{loop},'X'])(:,zloop),shiftvec1);
        SET(no).([types{loop},'Y'])(:,zloop) = circshift(SET(no).([types{loop},'Y'])(:,zloop),shiftvec1);
      end
    end

    %Interp Points
    types = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
    for loop = 1:length(types)
      if ~isempty(SET(no).([types{loop},'X']))
        SET(no).([types{loop},'X'])(:,zloop) = circshift(SET(no).([types{loop},'X'])(:,zloop),shiftvec1);
        SET(no).([types{loop},'Y'])(:,zloop) = circshift(SET(no).([types{loop},'Y'])(:,zloop),shiftvec1);
      end
    end
    types = {'LA','RA'};
    for loop = 1:length(types)
      if (~isempty(SET(no).(types{loop})) && ~isempty(SET(no).(types{loop}).InterpX))
        SET(no).(types{loop}).InterpX(:,zloop) = circshift(SET(no).(types{loop}).InterpX(:,zloop),shiftvec1);
        SET(no).(types{loop}).InterpY(:,zloop) = circshift(SET(no).(types{loop}).InterpY(:,zloop),shiftvec1);
      end
    end

    %Viability (ignored)

    %ROI's
    for rloop = 1:SET(no).RoiN
      if SET(no).Roi(rloop).Z == zloop
        SET(no).Roi(rloop).X = circshift(SET(no).Roi(rloop).X,shiftvec2);
        SET(no).Roi(rloop).Y = circshift(SET(no).Roi(rloop).Y,shiftvec2);
        SET(no).Roi(rloop).T = circshift(SET(no).Roi(rloop).T,shiftvec1);
      end
    end

    %Levelset
    %Ignore, see error check above.

    %MAR
    if ~isempty(SET(no).MaR)
      SET(no).MaR.Auto(:,:,:,zloop) = circshift(SET(no).MaR.Auto(:,:,:,zloop),shiftvec3);
      SET(no).MaR.Result(:,:,:,zloop) = circshift(SET(no).MaR.Result(:,:,:,zloop),shiftvec3);
      SET(no).MaR.Manual(:,:,:,zloop) = circshift(SET(no).MaR.Manual(:,:,:,zloop),shiftvec3);
      SET(no).MaR.MyocardMask(:,:,:,zloop) = circshift(SET(no).MaR.MyocardMask(:,:,:,zloop),shiftvec3);
    end

    %Measurements
    % Points
    %get points in the current slice
    if ~isempty(SET(no).Point)
      ind = find(SET(no).Point.Z == zloop);
      if ~isempty(ind)
        newtf = mod(SET(no).Point.T(ind)+shiftby, SET(no).TSize);
        newtf(newtf == 0) = SET(no).TSize;
        SET(no).Point.T(ind) = newtf;
      end
    end
    %get points in the current slice
    if ~isempty(SET(no).Measure)
      for nm = 1: length(SET(no).Measure)
        ind = find(SET(no).Measure(nm).Z == zloop, 1);
        if ~isempty(ind)
          newtf = mod(SET(no).Measure(nm).T+shiftby, SET(no).TSize);
          newtf(newtf==0) = SET(no).TSize;
          SET(no).Measure(nm).T = newtf;
        end
      end
    end
  end %end loop over slices

  if SET(no).EDT ~= 1 || SET(no).EST ~= 1
    %shift ED and ES definitions
    SET(no).EDT = mod(SET(no).EDT+shiftvec1(1),SET(no).TSize);
    SET(no).EST = mod(SET(no).EST+shiftvec1(1),SET(no).TSize);
  end
  SET(no).CurrentTimeFrame = mod(SET(no).CurrentTimeFrame+shiftvec1(1),SET(no).TSize);

end %end loop over nos

DATA.ViewIM{DATA.CurrentPanel} = [];
calcfunctions('calcvolume',NO);

if ~DATA.Batch && ~DATA.Autoloader
  panels = find(ismember(DATA.ViewPanels, nos));
  for p = panels
    viewfunctions('updatedrawlist',p)
    drawfunctions('drawpanel',p)
  end
end

DATA.updatevolumeaxes;
DATA.updatetimebaraxes;
segment('updateflow');

%------------------------------------------------------
function removeduplicatetimeframes_Callback(duplicates)
%------------------------------------------------------
%Removes duplicate timeframes, askes for the number of duplicates. Two
%means keep every second frame.
global SET NO

if nargin==0
  [duplicates,ok] = mygetnumber('Enter number of duplicate images','Duplicates',...
    2,... %default
    1,... %min
    10); %max

  if not(ok)
    myfailed('Invalid number of duplicates or aborted.');
    return;
  end
end

%Remove it
ind = false(1,SET(NO).TSize);
ind(1:(duplicates+1):end) = true;
removetimeframes(ind);

%------------------------------
function addnoise_Callback(f)
%------------------------------
%Adds noise to current image stack
global DATA SET NO

if nargin==0
  str1 = dprintf('Adding noise is irreversible.');
  str2 = dprintf('Are you sure?');
  str = sprintf('%s %s',str1,str2);
  if not(yesno(str,[],DATA.GUI.Segment))
    return;
  end

  [f,ok] = mygetnumber('Enter std of noise','STD',0.01,0,inf);
  if ~ok
    myfailed('Invalid STD or aborted.',DATA.GUI.Segment);
    return;
  end
end

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

SET(NO).IM = SET(NO).IM+f*randn(size(SET(NO).IM));

%segment('viewrefresh_Callback');
viewfunctions('setview')

%-----------------------------------------
function temporalmeanvalue_Callback(silent)
%-----------------------------------------
%Calculate temporal mean image by averaging over time.

%Written by Helen Soneson
global DATA SET NO

if nargin == 0
  silent = 0;
end

if isequal(SET(NO).TSize,1)
  if not(silent)
    myfailed('Already non time resolved. Aborted.',DATA.GUI.Segment);
    return;
  end
end

meanovertime = zeros(SET(NO).XSize,SET(NO).YSize,1,SET(NO).ZSize);   %mean temporal image
if silent == 0
  h = mywaitbarstart(SET(NO).ZSize,'Calculating mean over time',1,DATA.GUI.Segment);
end
for sliceloop = 1:SET(NO).ZSize
  for xloop = 1:SET(NO).XSize
    for yloop = 1:SET(NO).YSize
      meanovertime(xloop,yloop,1,sliceloop) = mean(SET(NO).IM(xloop,yloop,:,sliceloop));
    end
  end
  if silent == 0
    h = mywaitbarupdate(h);
  end
end
error = zeros(1,SET(NO).TSize);  %the difference between the mean temporal image and all time frames
for timeloop = 1:SET(NO).TSize
  error(timeloop) = sum(sum(sum(abs(squeeze(meanovertime(:,:,1,:))-squeeze(SET(NO).IM(:,:,timeloop,:))))));
end
if silent == 0
  mywaitbarclose(h);
end
nbr = SET(NO).XSize*SET(NO).YSize*SET(NO).ZSize;  %nbr of pixels in the image
error = error/nbr;
[~, timeframe] = min(error);  %find the time frame with the smallest difference to the mean temporal image
if silent == 0
  %plot the difference between the mean temporal image and all time frames
  figure,plot(timeframe,error(timeframe),'ro',1:SET(NO).TSize,error,'b'), grid on,xlabel(dprintf('Time frame')),ylabel(dprintf('Intensity error compared to mean image'))
  output{1,1} = 'timeframe';
  output{1,2} = 'error';
  for loop = 1:length(error)
    output{1+loop,1} = loop;
    output{1+loop,2} = error(loop);
  end
  %   segment('cell2clipboard',output);
end

oldNO = NO;
number = length(SET)+1;
NO = number;

for no = oldNO %nos
  NO = number;
  SET(NO) = SET(no);
  if no > 1
    for sliceloop = 1:SET(NO).ZSize
      for xloop = 1:SET(NO).XSize
        for yloop = 1:SET(NO).YSize
          meanovertime(xloop,yloop,1,sliceloop) = mean(SET(NO).IM(xloop,yloop,:,sliceloop));
        end
      end
    end
  end
  SET(NO).IM = meanovertime;

  %LV segmentation
  if not(isempty(SET(NO).EndoX))
    tempEndoX = SET(oldNO).EndoX(:,timeframe,:);
    tempEndoY = SET(oldNO).EndoY(:,timeframe,:);
    %     SET(NO).EndoX = nan(size(SET(NO).EndoX,1),1,size(SET(NO).EndoX,3));
    %     SET(NO).EndoY = nan(size(SET(NO).EndoY,1),1,size(SET(NO).EndoY,3));
    if ~isempty(tempEndoX)
      SET(NO).EndoX(:,1,:) = tempEndoX;
      SET(NO).EndoY(:,1,:) = tempEndoY;
    end
  end
  if not(isempty(SET(NO).EpiX))
    tempEpiX = SET(oldNO).EpiX(:,timeframe,:);
    tempEpiY = SET(oldNO).EpiY(:,timeframe,:);
    %     SET(NO).EpiX = nan(size(SET(NO).EpiX,1),1,size(SET(NO).EpiX,3));
    %     SET(NO).EpiY= nan(size(SET(NO).EpiY,1),1,size(SET(NO).EpiY,3));
    if ~isempty(tempEpiX)
      SET(NO).EpiX(:,1,:) = tempEpiX;
      SET(NO).EpiY(:,1,:) = tempEpiY;
    end
  end

  %RV segmentation
  if not(isempty(SET(NO).RVEndoX))
    tempEndoX = SET(oldNO).RVEndoX(:,timeframe,:);
    tempEndoY = SET(oldNO).RVEndoY(:,timeframe,:);
    %     SET(NO).RVEndoX = nan(size(SET(NO).RVEndoX,1),1,size(SET(NO).RVEndoX,3));
    %     SET(NO).RVEndoY = nan(size(SET(NO).RVEndoY,1),1,size(SET(NO).RVEndoY,3));
    if ~isempty(tempEndoX)
      SET(NO).RVEndoX(:,1,:) = tempEndoX;
      SET(NO).RVEndoY(:,1,:) = tempEndoY;
    end
  end
  if not(isempty(SET(NO).RVEpiX))
    tempEpiX = SET(oldNO).RVEpiX(:,timeframe,:);
    tempEpiY = SET(oldNO).RVEpiY(:,timeframe,:);
    %     SET(NO).RVEpiX = nan(size(SET(NO).RVEpiX,1),1,size(SET(NO).RVEpiX,3));
    %     SET(NO).RVEpiY = nan(size(SET(NO).RVEpiY,1),1,size(SET(NO).RVEpiY,3));
    if ~isempty(tempEpiX)
      SET(NO).RVEpiX(:,1,:) = tempEpiX;
      SET(NO).RVEpiY(:,1,:) = tempEpiY;
    end
  end

  %Set the parameter for a non-gated image stack
  SET(NO).Strain = [];
  SET(NO).Scar = [];
  SET(NO).Perfusion = [];
  SET(NO).PerfusionScoring = [];
  SET(NO).GadgetronScoring = [];
  SET(NO).Stress = [];
  SET(NO).TSize = 1;
  SET(NO).EndAnalysis = 1;
  SET(NO).CurrentTimeFrame =1;
  SET(NO).OrgTSize = 1;
  SET(NO).TIncr = 0;
  SET(NO).TimeVector = 0;
  SET(NO).TDelay = 0;
  SET(NO).EndoDraged = zeros(1,SET(NO).ZSize);
  SET(NO).EpiDraged = zeros(1,SET(NO).ZSize);
  SET(NO).EST = 1;
  SET(NO).EDT = 1;
  %LV volumes
  SET(NO).LVV = SET(oldNO).LVV(timeframe);
  SET(NO).EPV = SET(oldNO).EPV(timeframe);
  SET(NO).PV = SET(oldNO).PV(timeframe);
  SET(NO).LVM = SET(oldNO).LVM(timeframe);
  SET(NO).EDV = SET(NO).LVV(SET(NO).EDT)-SET(NO).PV(SET(NO).EDT);
  SET(NO).ESV = SET(NO).LVV(SET(NO).EST)-SET(NO).PV(SET(NO).EST);
  SET(NO).EF = 0;
  SET(NO).SV = 0;
  %RV volumes
  SET(NO).RVV = SET(oldNO).RVV(timeframe);
  SET(NO).RVEPV = SET(oldNO).RVEPV(timeframe);
  SET(NO).RVM = SET(oldNO).RVM(timeframe);
  SET(NO).RVEDV = SET(NO).RVV(SET(NO).EDT);
  SET(NO).RVESV = SET(NO).RVV(SET(NO).EST);
  SET(NO).RVSV = 0;
  SET(NO).RVEF = 0;
  %LV segmentation
  SET(NO).EndoXView = nan(size(SET(NO).EndoXView,1),1);
  SET(NO).EndoYView = nan(size(SET(NO).EndoYView,1),1);
  SET(NO).EpiXView = nan(size(SET(NO).EpiXView,1),1);
  SET(NO).EpiYView = nan(size(SET(NO).EpiYView,1),1);
  %ROIs
  if ~isempty(SET(NO).Roi)
    roistf=[];
    for roiloop = 1:SET(NO).RoiN
      if length(SET(NO).Roi(roiloop).T) == 1
        roistf(roiloop) = SET(NO).Roi(roiloop).T;
        SET(NO).Roi(roiloop).X = SET(NO).Roi(roiloop).X(:,timeframe);
        SET(NO).Roi(roiloop).Y = SET(NO).Roi(roiloop).Y(:,timeframe);
      else
        %erase time resolved rois
        roistf(roiloop) = NaN;
      end
    end
    eraserois = find(roistf ~= timeframe);
    if ~isempty(eraserois)
      %erase ROIs for other time frames
      roi('roidelete_Callback',eraserois);
    end
    roistf=[];
    for roiloop = 1:SET(NO).RoiN
      roistf(roiloop) = SET(NO).Roi(roiloop).T;
    end
    keeprois = find(roistf == timeframe);
    if ~isempty(keeprois)
      for roiloop = 1:length(keeprois)
        SET(NO).Roi(roiloop).T = 1;
      end
    end
  end
  %LevelSet
  if ~isempty(SET(NO).LevelSet)
    logdisp('3D segmentation data is erased in the mean temporal image stack.',true);
    SET(NO).LevelSet = [];
  end

  %Add stack links
  SET(oldNO).Children = [SET(oldNO).Children NO];
  SET(NO).Parent = []; %oldNO
  SET(NO).Children = [];
  SET(NO).Linked = NO;
  %   newlinkies = [SET(oldNO).Linked NO];
  %   for sloop = newlinkies
  %     SET(sloop).Linked = newlinkies;
  %   end

  %remove Strain tagging analysis
  if not(isempty(SET(NO).StrainTagging))
    logdisp('Merging to a mean image stack clears the Strain tagging analysis',true);
    if ~isempty(DATA.GUI.StrainTagging)
      straintagging.straintagging('close_Callback');
    end
    SET(NO).StrainTagging = [];
  end
  %remove Strain MITT analysis
  logmsg = 'Merging to a mean image stack clears the Strain MITT analysis';
  strainmitt.strainmitt('removestrainanalysis',NO,logmsg);

  number = number+1;
end
if silent == 0
  viewfunctions('setview');
end

%--------------------------------------
function extraslice_Callback(type)
%--------------------------------------
%Add an extra slice. Type is either 'basal' or 'apical'
%Takes also care of delineations.

global DATA SET NO

if nargin==0
  type = 'basal';
end

%--- Extra basal
switch type
  case 'basal'
    %Duplicate image slice
    SET(NO).IM = cat(4,SET(NO).IM(:,:,:,1),SET(NO).IM);
    SET(NO).ZSize=SET(NO).ZSize+1;

    if ~isempty(SET(NO).EndoPinX)
      SET(NO).EndoPinX = cat(2,SET(NO).EndoPinX(:,1),SET(NO).EndoPinX);
      SET(NO).EndoPinY = cat(2,SET(NO).EndoPinY(:,1),SET(NO).EndoPinY);
    end
    if ~isempty(SET(NO).EpiPinX)
      SET(NO).EpiPinX = cat(2,SET(NO).EpiPinX(:,1),SET(NO).EpiPinX);
      SET(NO).EpiPinY = cat(2,SET(NO).EpiPinY(:,1),SET(NO).EpiPinY);
    end
    if ~isempty(SET(NO).RVEndoPinX)
      SET(NO).RVEndoPinX = cat(2,SET(NO).RVEndoPinX(:,1),SET(NO).RVEndoPinX);
      SET(NO).RVEndoPinY = cat(2,SET(NO).RVEndoPinY(:,1),SET(NO).RVEndoPinY);
    end
    if ~isempty(SET(NO).RVEpiPinX)
      SET(NO).RVEpiPinY = cat(2,SET(NO).RVEpiPinY(:,1),SET(NO).RVEpiPinY);
      SET(NO).RVEpiPinX = cat(2,SET(NO).RVEpiPinX(:,1),SET(NO).RVEpiPinX);
    end

    if ~isempty(SET(NO).EndoInterpX)
      SET(NO).EndoInterpX = cat(2,SET(NO).EndoInterpX(:,1),SET(NO).EndoInterpX);
      SET(NO).EndoInterpY = cat(2,SET(NO).EndoInterpY(:,1),SET(NO).EndoInterpY);
    end
    if ~isempty(SET(NO).EpiInterpX)
      SET(NO).EpiInterpX = cat(2,SET(NO).EpiInterpX(:,1),SET(NO).EpiInterpX);
      SET(NO).EpiInterpY = cat(2,SET(NO).EpiInterpY(:,1),SET(NO).EpiInterpY);
    end
    if ~isempty(SET(NO).RVEndoInterpX)
      SET(NO).RVEndoInterpX = cat(2,SET(NO).RVEndoInterpX(:,1),SET(NO).RVEndoInterpX);
      SET(NO).RVEndoInterpY = cat(2,SET(NO).RVEndoInterpY(:,1),SET(NO).RVEndoInterpY);
    end
    if ~isempty(SET(NO).RVEpiInterpX)
      SET(NO).RVEpiInterpX = cat(2,SET(NO).RVEpiInterpX(:,1),SET(NO).RVEpiInterpX);
      SET(NO).RVEpiInterpY = cat(2,SET(NO).RVEpiInterpY(:,1),SET(NO).RVEpiInterpY);
    end

    DATA.ViewIM{DATA.CurrentPanel} = []; %Clear data panel

    %Copy from systole
    tf = SET(NO).EST;

    if ~isempty(SET(NO).EndoX)
      SET(NO).EndoX = cat(3,SET(NO).EndoX(:,:,1),SET(NO).EndoX);
      SET(NO).EndoY = cat(3,SET(NO).EndoY(:,:,1),SET(NO).EndoY);
      SET(NO).EndoX(:,1,1) = SET(NO).EndoX(:,tf,2);
      SET(NO).EndoY(:,1,1) = SET(NO).EndoY(:,tf,2);
      %Remove in systole basal slice
      SET(NO).EndoX(:,tf,1) = NaN;
      SET(NO).EndoY(:,tf,1) = NaN;
    end

    if ~isempty(SET(NO).EpiX)
      SET(NO).EpiX = cat(3,SET(NO).EpiX(:,:,1),SET(NO).EpiX);
      SET(NO).EpiY = cat(3,SET(NO).EpiY(:,:,1),SET(NO).EpiY);
      SET(NO).EpiX(:,1,1) = SET(NO).EpiX(:,tf,2);
      SET(NO).EpiY(:,1,1) = SET(NO).EpiY(:,tf,2);
      %Remove in systole basal slice
      SET(NO).EpiX(:,tf,1) = NaN;
      SET(NO).EpiY(:,tf,1) = NaN;
    end

    if ~isempty(SET(NO).RVEndoX)
      SET(NO).RVEndoX = cat(3,SET(NO).RVEndoX(:,:,1),SET(NO).RVEndoX);
      SET(NO).RVEndoY = cat(3,SET(NO).RVEndoY(:,:,1),SET(NO).RVEndoY);
      SET(NO).RVEndoX(:,1,1) = SET(NO).RVEndoX(:,tf,2);
      SET(NO).RVEndoY(:,1,1) = SET(NO).RVEndoY(:,tf,2);
      SET(NO).RVEndoX(:,tf,1) = NaN;
      SET(NO).RVEndoY(:,tf,1) = NaN;
    end

    if ~isempty(SET(NO).RVEpiX)
      SET(NO).RVEpiX = cat(3,SET(NO).RVEpiX(:,:,1),SET(NO).RVEpiX);
      SET(NO).RVEpiY = cat(3,SET(NO).RVEpiY(:,:,1),SET(NO).RVEpiY);
      SET(NO).RVEpiX(:,1,1) = SET(NO).RVEpiX(:,tf,2);
      SET(NO).RVEpiY(:,1,1) = SET(NO).RVEpiY(:,tf,2);
      %Remove in systole basal slice
      SET(NO).RVEpiX(:,tf,1) = NaN;
      SET(NO).RVEpiY(:,tf,1) = NaN;
    end

    %Draged (t*z)
    SET(NO).EndoDraged = cat(2,SET(NO).EndoDraged(:,1),SET(NO).EndoDraged);
    SET(NO).EpiDraged = cat(2,SET(NO).EpiDraged(:,1),SET(NO).EpiDraged);

    %Scar
    if not(isempty(SET(NO).Scar))
      SET(NO).Scar.IM = squeeze(SET(NO).IM(:,:,1,:));
      SET(NO).Scar.Auto = cat(3,SET(NO).Scar.Auto(:,:,1),SET(NO).Scar.Auto);
      SET(NO).Scar.Result = cat(3,SET(NO).Scar.Result(:,:,1),SET(NO).Scar.Result);
      SET(NO).Scar.Manual = cat(3,SET(NO).Scar.Manual(:,:,1),SET(NO).Scar.Manual);
      SET(NO).Scar.MyocardMask = cat(3,SET(NO).Scar.MyocardMask(:,:,1),SET(NO).Scar.MyocardMask);
      SET(NO).Scar.NoReflow = cat(3,SET(NO).Scar.NoReflow(:,:,1),SET(NO).Scar.NoReflow);
      if numel(SET(NO).Scar.mthreshold) == SET(NO).ZSize
        SET(NO).Scar.mthreshold = cat(3,SET(NO).Scar.mthreshold(1),SET(NO).Scar.mthreshold);
      end
    end

    %MaR
    if not(isempty(SET(NO).MaR))
      SET(NO).MaR.Auto = cat(4,SET(NO).MaR.Auto(:,:,:,1),SET(NO).MaR.Auto);
      SET(NO).MaR.Result = cat(4,SET(NO).MaR.Result(:,:,:,1),SET(NO).MaR.Result);
      SET(NO).MaR.Manual = cat(4,SET(NO).MaR.Manual(:,:,:,1),SET(NO).MaR.Manual);
      SET(NO).MaR.MyocardMask = cat(4,SET(NO).MaR.MyocardMask(:,:,:,1),SET(NO).MaR.MyocardMask);
    end

    %ROI
    if SET(NO).RoiN > 0
      roiz = {SET(NO).Roi.Z};
      roiz = cellfun(@(x)(x+1),roiz,'UniformOutput',false);
      [SET(NO).Roi.Z] = deal(roiz{:});
    end

    %Measure
    if ~isempty(SET(NO).Measure)
      measz = {SET(NO).Measure.Z};
      measz = cellfun(@(x)(x+1),measz,'UniformOutput',false);
      [SET(NO).Measure.Z] = deal(measz{:});
    end

    %Annotation points
    if isfield(SET(NO).Point,'Z')
      SET(NO).Point.Z = SET(NO).Point.Z + 1;
    end

    zdir = cross(...
      SET(NO).ImageOrientation(1:3),...
      SET(NO).ImageOrientation(4:6));
    zdir = zdir(:)';
    SET(NO).ImagePosition = SET(NO).ImagePosition+zdir*(SET(NO).SliceThickness+SET(NO).SliceGap);

    %LA, RA
    heartparts = {'LA','RA'};
    for loop = 1:numel(heartparts)
      hp = heartparts{loop};
      if ~isempty(SET(NO).(hp)) && ~isempty(SET(NO).(hp).X)
        %remap contour's coordinates
        numpoints = SET(NO).(hp).getnumpoints;
        tempx = nan(numpoints,SET(NO).TSize,SET(NO).ZSize);
        tempy = nan(numpoints,SET(NO).TSize,SET(NO).ZSize);
        tempx(:,:,2:end) = SET(NO).(hp).X;
        tempy(:,:,2:end) = SET(NO).(hp).Y;
        SET(NO).(hp).X = tempx;
        SET(NO).(hp).Y = tempy;

        %reset interpolation points
        SET(NO).(hp).InterpX = [];
        SET(NO).(hp).InterpY = [];
      end
    end

  case 'apical'
    %Duplicate image slice
    SET(NO).IM = cat(4,SET(NO).IM,SET(NO).IM(:,:,:,end));
    SET(NO).ZSize=SET(NO).ZSize+1;

    if ~isempty(SET(NO).EndoPinX)
      SET(NO).EndoPinX = cat(2,SET(NO).EndoPinX,SET(NO).EndoPinX(:,end));
      SET(NO).EndoPinY = cat(2,SET(NO).EndoPinY,SET(NO).EndoPinY(:,end));
    end
    if ~isempty(SET(NO).EpiPinX)
      SET(NO).EpiPinY = cat(2,SET(NO).EpiPinY,SET(NO).EpiPinY(:,end));
      SET(NO).EpiPinX = cat(2,SET(NO).EpiPinX,SET(NO).EpiPinX(:,end));
    end
    if ~isempty(SET(NO).RVEndoPinX)
      SET(NO).RVEndoPinX = cat(2,SET(NO).RVEndoPinX,SET(NO).RVEndoPinX(:,end));
      SET(NO).RVEndoPinY = cat(2,SET(NO).RVEndoPinY,SET(NO).RVEndoPinY(:,end));
    end
    if ~isempty(SET(NO).RVEpiPinX)
      SET(NO).RVEpiPinY = cat(2,SET(NO).RVEpiPinY,SET(NO).RVEpiPinY(:,end));
      SET(NO).RVEpiPinX = cat(2,SET(NO).RVEpiPinX,SET(NO).RVEpiPinX(:,end));
    end

    if ~isempty(SET(NO).EndoInterpX)
      SET(NO).EndoInterpX = cat(2,SET(NO).EndoInterpX,SET(NO).EndoInterpX(:,end));
      SET(NO).EndoInterpY = cat(2,SET(NO).EndoInterpY,SET(NO).EndoInterpY(:,end));
    end
    if ~isempty(SET(NO).EpiInterpX)
      SET(NO).EpiInterpY = cat(2,SET(NO).EpiInterpY,SET(NO).EpiInterpY(:,end));
      SET(NO).EpiInterpX = cat(2,SET(NO).EpiInterpX,SET(NO).EpiInterpX(:,end));
    end
    if ~isempty(SET(NO).RVEndoInterpX)
      SET(NO).RVEndoInterpX = cat(2,SET(NO).RVEndoInterpX,SET(NO).RVEndoInterpX(:,end));
      SET(NO).RVEndoInterpY = cat(2,SET(NO).RVEndoInterpY,SET(NO).RVEndoInterpY(:,end));
    end
    if ~isempty(SET(NO).RVEpiInterpX)
      SET(NO).RVEpiInterpY = cat(2,SET(NO).RVEpiInterpY,SET(NO).RVInterpInterpY(:,end));
      SET(NO).RVEpiInterpX = cat(2,SET(NO).RVEpiInterpX,SET(NO).RVInterpInterpX(:,end));
    end

    DATA.ViewIM{DATA.CurrentPanel} = [];

    if ~isempty(SET(NO).EndoX)
      SET(NO).EndoX = cat(3,SET(NO).EndoX,SET(NO).EndoX(:,:,end));
      SET(NO).EndoY = cat(3,SET(NO).EndoY,SET(NO).EndoY(:,:,end));
    end

    if ~isempty(SET(NO).EpiX)
      SET(NO).EpiY = cat(3,SET(NO).EpiY,SET(NO).EpiY(:,:,end));
      SET(NO).EpiX = cat(3,SET(NO).EpiX,SET(NO).EpiX(:,:,end));
    end

    if ~isempty(SET(NO).RVEndoX)
      SET(NO).RVEndoX = cat(3,SET(NO).RVEndoX,SET(NO).RVEndoX(:,:,end));
      SET(NO).RVEndoY = cat(3,SET(NO).RVEndoY,SET(NO).RVEndoY(:,:,end));
    end

    if ~isempty(SET(NO).RVEpiX)
      SET(NO).RVEpiY = cat(3,SET(NO).RVEpiY,SET(NO).RVEpiY(:,:,end));
      SET(NO).RVEpiX = cat(3,SET(NO).RVEpiX,SET(NO).RVEpiX(:,:,end));
    end

    %Draged (t*z)
    SET(NO).EndoDraged = cat(2,SET(NO).EndoDraged,SET(NO).EndoDraged(:,1));
    SET(NO).EpiDraged = cat(2,SET(NO).EpiDraged,SET(NO).EpiDraged(:,1));

    %Scar
    if not(isempty(SET(NO).Scar))
      SET(NO).Scar.IM = squeeze(SET(NO).IM(:,:,1,:));
      SET(NO).Scar.Auto = cat(3,SET(NO).Scar.Auto,SET(NO).Scar.Auto);
      SET(NO).Scar.Result = cat(3,SET(NO).Scar.Result,SET(NO).Scar.Result(:,:,end));
      SET(NO).Scar.Manual = cat(3,SET(NO).Scar.Manual,SET(NO).Scar.Manual(:,:,end));
      SET(NO).Scar.MyocardMask = cat(3,SET(NO).Scar.MyocardMask,SET(NO).Scar.MyocardMask(:,:,end));
      SET(NO).Scar.NoReflow = cat(3,SET(NO).Scar.NoReflow,SET(NO).Scar.NoReflow(:,:,end));
      if numel(SET(NO).Scar.mthreshold) == SET(NO).ZSize
        SET(NO).Scar.mthreshold = cat(3,SET(NO).Scar.mthreshold,SET(NO).Scar.mthreshold(end));
      end
    end

    %Scar
    if not(isempty(SET(NO).MaR))
      SET(NO).MaR.Auto = cat(4,SET(NO).MaR.Auto,SET(NO).MaR.Auto);
      SET(NO).MaR.Result = cat(4,SET(NO).MaR.Result,SET(NO).MaR.Result(:,:,:,end));
      SET(NO).MaR.Manual = cat(4,SET(NO).MaR.Manual,SET(NO).MaR.Manual(:,:,:,end));
      SET(NO).MaR.MyocardMask = cat(4,SET(NO).MaR.MyocardMask,SET(NO).MaR.MyocardMask(:,:,:,end));
    end

    %LA, RA
    heartparts = {'LA','RA'};
    for loop = 1:numel(heartparts)
      hp = heartparts{loop};
      if ~isempty(SET(NO).(hp)) && ~isempty(SET(NO).(hp).X)
        %remap contour's coordinates 
        numpoints = SET(NO).(hp).getnumpoints;
        tempx = nan(numpoints,SET(NO).TSize,SET(NO).ZSize);
        tempy = nan(numpoints,SET(NO).TSize,SET(NO).ZSize);
        tempx(:,:,1:end-1) = SET(NO).(hp).X;
        tempy(:,:,1:end-1) = SET(NO).(hp).Y;
        SET(NO).(hp).X = tempx;
        SET(NO).(hp).Y = tempy;

        %reset interpolation points
        SET(NO).(hp).InterpX = [];
        SET(NO).(hp).InterpY = [];
      end
    end
end

%remove Strain tagging analysis
if not(isempty(SET(NO).StrainTagging))
  logdisp('Adding a slices clears the Strain tagging analysis',true);
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  SET(NO).StrainTagging = [];
end

%remove Strain MITT analysis
logmsg = 'Adding a slice clears the Strain MITT analysis';
strainmitt.strainmitt('removestrainanalysis',NO,logmsg);

%segment('viewrefresh_Callback');
viewfunctions('setview')

%----------------------------
function kspace_Callback
%----------------------------
%Shows KSPACE of image (Fourier Transform). Displays the log of the
%Fouries Transform.

global SET NO

tempnos = NO;
imissingle = classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

f = fftshift(fft2(SET(NO).IM(:,:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice)));

figure(123);
set(123,'NumberTitle','off','Name','K-space');
imagesc([-1 1],[-1 1],log(abs(f)));
lbl = 'Normalized freq';
xlabel(lbl);
ylabel(lbl);
ttl = 'log(abs(fft(im)))';
title(ttl);

%--------------------------------------
function intensitymapping_Callback
%--------------------------------------
%Plots intensity mapping function.
global SET NO

c = SET(NO).IntensityMapping.Contrast;
b = SET(NO).IntensityMapping.Brightness-0.5;

figure(5);
temp = c*calcfunctions('returnmapping',NO)+b;
plot(...
  linspace(0,1,length(temp)),...
  temp,'k-');
hold on;
plot([0 1],[1 1],'k--');
plot([0 1],[0 0],'k--');
title(dprintf('Image intensity mapping'));
axis([0 1 min(0,min(temp(:))) max(1,max(temp(:)))])

%---------------------------------------
function viewtrueintensity_Callback
%---------------------------------------
%View true image intensities in current timeframe and slice.
global SET NO

figure(5);
set(5,'numbertitle','off','name','True image intensity');

temp = SET(NO).IM(:,:,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice);
if isempty(SET(NO).Flow)
  imagesc(calcfunctions('calctruedata',temp,NO));
else
  imagesc(2*SET(NO).VENC*(-0.5+temp));
  set(gca,'clim',[-SET(NO).VENC SET(NO).VENC]);
end

colormap(jet);
axis image off;
colorbar;

%--------------------------------------
function duplicateimagestack_Callback(frompreviewmenu)
%--------------------------------------
%Duplicates current image stack. This function breaks links to
%linked datasets, parent and children.

global DATA SET NO

if nargin > 0 && frompreviewmenu %from context menu
  dupno = segment_main('thumbnailno');
else
  dupno = NO;
end
if DATA.ShowRelevantStacksOnly && strcmpi(DATA.CurrentTheme,'roi')
  mywarning('This is not implemented for flow images.',DATA.GUI.Segment)
  return
end

%Duplicate current dataset
no = length(SET) + 1;
SET(no) = SET(dupno);

if DATA.issegment3dp
  if ~isempty(SET(no).LevelSet)
    SET(no).LevelSet.Object = copy(SET(dupno).LevelSet.Object); %Ensure that it is not copied directly
    SET(no).LevelSet.Object.NO = no;
  end
end

% functions to ensure correct duplication/creation of a new instance of the LA or RA object
generalpen.atriumpenfunctions('duplicate','LA',no);
generalpen.atriumpenfunctions('duplicate','RA',no);

% Delete Strain MITT
SET(no).StrainMitt = [];

SET(no).Linked = no;
SET(no).Parent = [];
SET(no).Children = [];

if ~isempty(SET(no).Flow) || length(SET(no).Linked) > 1
  mywarning('Image stack contains flow data and/or other stack links, not duplicating internal references.',DATA.GUI.Segment);
  SET(no).Flow = [];
end

if ~contains(DATA.CurrentTheme,'3dp')% dont do switch if 3dp mode
  viewfunctions('setview',1,1,no,{'one'})
end
calcpreview = false;
drawfunctions('drawthumbnails',calcpreview);

%-----------------------------------------------
function upsampleimage_Callback(f,inputNO,silent)
%-----------------------------------------------
%Upsamples current image stack (in slice only). Takes care of
%segmentation in the upsampling process.

global DATA SET NO

if nargin==0
  str1 = dprintf('Upsampling/downsampling image stack is irreversible.');
  str2 = dprintf('Are you sure?');
  str = sprintf('%s %s',str1,str2);
  if not(yesno(str,[],DATA.GUI.Segment))
    return;
  end

  %Get factor
  [f,ok] = mygetnumber('Enter upsampling/downsampling factor (<1 => downsample)','Resampling Factor',2,0.001,10);
  if ~ok
    myfailed('Invalid factor or aborted.',DATA.GUI.Segment);
    return;
  end
end

if f == 1
  return
end

if nargin<2
  inputNO = NO;
end

if nargin<3
  silent = DATA.Silent;
end

oldsize = [SET(inputNO).XSize SET(inputNO).YSize SET(inputNO).ZSize];

%Find image stacks to crop
nos = SET(inputNO).Linked; %Crop current image stack

for no = nos

  %Restore image enhancement
  if isfield(SET(no),'IMBackup') && (~isempty(SET(no).IMBackup))
    SET(no).IM = SET(no).IMBackup;
    SET(no).IMBackup = [];
  end

  %Resample
  SET(no).IM = upsampleinplane(f,SET(no).IM,silent);

  %--- Things to change

  %Sizes
  SET(no).XSize = size(SET(no).IM,1);
  SET(no).YSize = size(SET(no).IM,2);

  %Delineation
  if ~isempty(SET(no).EndoX)
    SET(no).EndoX = resamplehelper(f,SET(no).EndoX);
    SET(no).EndoY = resamplehelper(f,SET(no).EndoY);
  end

  if ~isempty(SET(no).EpiX)
    SET(no).EpiX = resamplehelper(f,SET(no).EpiX);
    SET(no).EpiY = resamplehelper(f,SET(no).EpiY);
  end

  if ~isempty(SET(no).RVEndoX)
    SET(no).RVEndoX = resamplehelper(f,SET(no).RVEndoX);
    SET(no).RVEndoY = resamplehelper(f,SET(no).RVEndoY);
  end

  if ~isempty(SET(no).RVEpiX)
    SET(no).RVEpiX = resamplehelper(f,SET(no).RVEpiX);
    SET(no).RVEpiY = resamplehelper(f,SET(no).RVEpiY);
  end

  if ~isempty(SET(no).EndoPinX)
    SET(no).EndoPinX = resamplehelper(f,SET(no).EndoPinX);
    SET(no).EndoPinY = resamplehelper(f,SET(no).EndoPinY);
  end
  if ~isempty(SET(no).EpiPinX)
    SET(no).EpiPinX = resamplehelper(f,SET(no).EpiPinX);
    SET(no).EpiPinY = resamplehelper(f,SET(no).EpiPinY);
  end
  if ~isempty(SET(no).RVEndoPinX)
    SET(no).RVEndoPinX = resamplehelper(f,SET(no).RVEndoPinX);
    SET(no).RVEndoPinY = resamplehelper(f,SET(no).RVEndoPinY);
  end
  if ~isempty(SET(no).RVEpiPinX)
    SET(no).RVEpiPinX = resamplehelper(f,SET(no).RVEpiPinX);
    SET(no).RVEpiPinY = resamplehelper(f,SET(no).RVEpiPinY);
  end

  if ~isempty(SET(no).EndoInterpX)
    SET(no).EndoInterpX = resamplehelper(f,SET(no).EndoInterpX);
    SET(no).EndoInterpY = resamplehelper(f,SET(no).EndoInterpY);
  end
  if ~isempty(SET(no).EpiInterpX)
    SET(no).EpiInterpX = resamplehelper(f,SET(no).EpiInterpX);
    SET(no).EpiInterpY = resamplehelper(f,SET(no).EpiInterpY);
  end
  if ~isempty(SET(no).RVEndoInterpX)
    SET(no).RVEndoInterpX = resamplehelper(f,SET(no).RVEndoInterpX);
    SET(no).RVEndoInterpY = resamplehelper(f,SET(no).RVEndoInterpY);
  end
  if ~isempty(SET(no).RVEpiInterpX)
    SET(no).RVEpiInterpX = resamplehelper(f,SET(no).RVEpiInterpX);
    SET(no).RVEpiInterpY = resamplehelper(f,SET(no).RVEpiInterpY);
  end

  if ~isempty(SET(no).LA) && ~isempty(SET(no).LA.X)
    SET(no).LA.X = resamplehelper(f,SET(no).LA.X);
    SET(no).LA.Y = resamplehelper(f,SET(no).LA.Y);
  end

  if ~isempty(SET(no).LA) && ~isempty(SET(no).LA.InterpX)
    SET(no).LA.InterpX = resamplehelper(f,SET(no).LA.InterpX);
    SET(no).LA.InterpY = resamplehelper(f,SET(no).LA.InterpY);
  end

  if ~isempty(SET(no).RA) && ~isempty(SET(no).RA.X)
    SET(no).RA.X = resamplehelper(f,SET(no).RA.X);
    SET(no).RA.Y = resamplehelper(f,SET(no).RA.Y);
  end

  if ~isempty(SET(no).RA) && ~isempty(SET(no).RA.InterpX)
    SET(no).RA.InterpX = resamplehelper(f,SET(no).RA.InterpX);
    SET(no).RA.InterpY = resamplehelper(f,SET(no).RA.InterpY);
  end

  %Point
  if ~isempty(SET(no).Point)
    SET(no).Point.X = resamplehelper(f,SET(no).Point.X);
    SET(no).Point.Y = resamplehelper(f,SET(no).Point.Y);
  end

  %Measurement
  if ~isempty(SET(no).Measure)
    for mloop = 1:length(SET(no).Measure)
      SET(no).Measure(mloop).X = resamplehelper(f,SET(no).Measure(mloop).X);
      SET(no).Measure(mloop).Y = resamplehelper(f,SET(no).Measure(mloop).Y);
    end
  end

  %Resolution
  SET(no).ResolutionX = SET(no).ResolutionX/f;
  SET(no).ResolutionY = SET(no).ResolutionY/f;

  %Mmode stuff
  SET(no).Mmode.X = resamplehelper(f,SET(no).Mmode.X);
  SET(no).Mmode.Y = resamplehelper(f,SET(no).Mmode.Y);
  SET(no).Mmode.M1 = resamplehelper(f,SET(no).Mmode.M1);
  SET(no).Mmode.M2 = resamplehelper(f,SET(no).Mmode.M2);

  %Center
  SET(no).CenterX = resamplehelper(f,SET(no).CenterX);
  SET(no).CenterY = resamplehelper(f,SET(no).CenterY);

  %ROI's
  for rloop=1:SET(no).RoiN
    SET(no).Roi(rloop).X = resamplehelper(f,SET(no).Roi(rloop).X);
    SET(no).Roi(rloop).Y = resamplehelper(f,SET(no).Roi(rloop).Y);
  end

  %Scar
  if not(isempty(SET(no).Scar))
    SET(no).Scar.Auto = upsampleinplane(f,SET(no).Scar.Auto)>0.5;
    SET(no).Scar.Result = upsampleinplane(f,SET(no).Scar.Result)>0.5;
    SET(no).Scar.Manual = int8(upsampleinplane(f,SET(no).Scar.Manual));
    SET(no).Scar.Manual(SET(no).Scar.Manual<0) = int8(-1);
    SET(no).Scar.MyocardMask = upsampleinplane(f,SET(no).Scar.MyocardMask)>0.5;
    SET(no).Scar.NoReflow = upsampleinplane(f,SET(no).Scar.NoReflow)>0.5;
    SET(no).Scar.IM = upsampleinplane(f,SET(no).Scar.IM);
    if ~isempty(SET(no).Scar.GreyZone)
      SET(no).Scar.GreyZone.map = uint8(upsampleinplane(f,SET(no).Scar.GreyZone.map));
    end
  end

  %MaR
  if not(isempty(SET(no).MaR))
    SET(no).MaR.Auto = upsampleinplane(f,SET(no).MaR.Auto)>0.5;
    SET(no).MaR.Result = upsampleinplane(f,SET(no).MaR.Result)>0.5;
    SET(no).MaR.Manual = int8(upsampleinplane(f,SET(no).MaR.Manual));
    SET(no).MaR.Manual(SET(no).MaR.Manual<0) = int8(-1);
    SET(no).MaR.MyocardMask = upsampleinplane(f,SET(no).MaR.MyocardMask)>0.5;
  end

  %Phase correction
  if not(isempty(SET(no).Flow))
    if not(isempty(SET(no).Flow.PhaseCorr))
      SET(no).Flow.PhaseCorr = upsampleinplane(f,SET(no).Flow.PhaseCorr);
    end
  end

  %Zoom state
  SET(no).NormalZoomState = [];
  SET(no).MontageZoomState = [];
  SET(no).MontageRowZoomState = [];

  %Current slice
  SET(no).CurrentSlice = min(max(1,round(SET(no).CurrentSlice*f)),SET(no).ZSize);

  %Resample .LevelSet
  if DATA.issegment3dp && ~isempty(SET(no).LevelSet)
    O = SET(no).LevelSet.Object;

    %resample stored objects, important to resample stored object before
    for obj = 1:O.getnumobjects()

      if O.iscompressed(obj)
        O.uncompress(obj,oldsize);
      end

      temp = O.getbw(obj);

      if ~isempty(temp)
        temp = upsampleinplane(f,temp,silent);
        O.setbw(obj,temp);
        segment3dp.tools('setobjectcut',obj,[],no);
        segment3dp.tools('setobjectdeform',obj,[],no);
        O.compress(obj);
      end

    end

    %---things to change

    %Zoomstate
    SET(no).LevelSet.View.RZoomState = [];
    SET(no).LevelSet.View.GZoomState = [];
    SET(no).LevelSet.View.BZoomState = [];

    %slice
    %slice
    switch lower(DATA.LevelSet.imageorientation)
      case 'transversal'
        fieldy = 'GSlice';
        fieldx = 'BSlice';
      case 'sagittal'
        fieldy = 'BSlice';
        fieldx = 'RSlice';
      case 'coronal'
        fieldy = 'GSlice';
        fieldx = 'RSlice';
    end

    SET(no).LevelSet.View.(fieldx) = max(1,min(round(resamplehelper(f,SET(no).LevelSet.View.(fieldx))),SET(no).XSize));
    SET(no).LevelSet.View.(fieldy) = max(1,min(round(resamplehelper(f,SET(no).LevelSet.View.(fieldy))),SET(no).YSize));

    %resets zoom
    silentmessage = true;
    silentupdate = true;
    segment3dp.tools('init3DP',silentmessage,silentupdate);
  end

  %remove Strain tagging analysis
  if not(isempty(SET(no).StrainTagging))
    logdisp('Upsample image stack clears the Strain tagging analysis',true);
    if ~isempty(DATA.GUI.StrainTagging)
      straintagging.straintagging('close_Callback');
    end
    SET(no).StrainTagging = [];
  end

  %remove Strain MITT analysis
  logmsg = 'Upsample image stack clears the Strain MITT analysis';
  strainmitt.strainmitt('removestrainanalysis',no,logmsg);


  if isfield(SET(no),'RV') && ~isempty(SET(no).RV) && ~isempty(SET(no).RV.centerbasal)
    SET(no).RV.centerbasal = resamplehelper(f,SET(no).RV.centerbasal);
    SET(no).RV.centerapical = resamplehelper(f,SET(no).RV.centerapical);
  end

  %CT LV reconstructions points
  try
    if ~isempty(SET(no).CT)
      SET(no).CT.Apex(1)=min(size(SET(no).IM,1),max(1,round(resamplehelper(f,SET(no).CT.Apex(1)))));
      SET(no).CT.Apex(2)=min(size(SET(no).IM,2),max(1,round(resamplehelper(f,SET(no).CT.Apex(2)))));
      SET(no).CT.AV(1)=min(size(SET(no).IM,1),max(1,round(resamplehelper(f,SET(no).CT.AV(1)))));
      SET(no).CT.AV(2)=min(size(SET(no).IM,2),max(1,round(resamplehelper(f,SET(no).CT.AV(2)))));
      SET(no).CT.RVPulmi(1)=min(size(SET(no).IM,1),max(1,round(resamplehelper(f,SET(no).CT.RVPulmi(1)))));
      SET(no).CT.RVPulmi(2)=min(size(SET(no).IM,2),max(1,round(resamplehelper(f,SET(no).CT.RVPulmi(2)))));
    end
  catch
    SET(no).CT = [];
  end

end

DATA.EndoEDGE0 = [];
DATA.EndoEDGE1 = [];
DATA.EndoEDGE2 = [];
DATA.EndoEDGE3 = [];
DATA.EpiEDGE0 = [];
DATA.EpiEDGE1 = [];
DATA.EpiEDGE2 = [];
DATA.EpiEDGE3 = [];
DATA.EndoEdgeDetected = false;
DATA.EpiEdgeDetected = false;
disableundo;

[tf,loc] = ismember(nos,DATA.ViewPanels);
[DATA.ViewIM{loc(tf)}] = deal([]);
drawfunctions('drawno',no);
for p = loc(tf)
  drawfunctions('drawselectedslice',p)
  drawfunctions('drawselectedframe',p)
end

%--------------------------------------
function upsampletemporal_Callback(f)
%--------------------------------------
%Usamples current image stack in time. f is upsampling factor.
%Takes care of segmentation.

global DATA SET NO

if nargin==0
  %Get factor
  [f,ok] = mygetnumber('Enter upsampling/downsampling factor (<1 => downsample)','Resampling Factor',2,0,10);
  if ~ok
    myfailed('Invalid upsampling/downsampling factor or aborted.',DATA.GUI.Segment);
    return;
  end
end

% check for still image
forcecopy = false;
if size(SET(NO).IM,3)==1
  if ~yesno('Image has no temporal dimension. Do you want to proceed?')
    return
  end
  forcecopy = true;
  SET(NO).TimeVector = linspace(0,1,f);
end

if not(isrectilinear(SET(NO).TimeVector))
  if ~yesno(['Non uniform time steps detected. This is not '...
      'supported for upsampling. Mean time step will be used ' ...
      'and time step information will be lost. Proceed?'])
    return
  end
end

% Ask about Fourier upsampling
immethod = 'image';
if f>1
  if yesno('Use Fourier temporal upsampling?')
    immethod = 'image-fourier';
  end
end

%Find image stacks to upsample temporal
nos = SET(NO).Linked; %Upsample temporal
if ~isempty(SET(NO).Flow)
  nos = [nos SET(NO).Flow.PhaseNo SET(NO).Flow.PhaseX SET(NO).Flow.PhaseY SET(NO).Flow.Angio SET(NO).Flow.VelMag SET(NO).Flow.MagnitudeNo];
  SET(SET(NO).Flow.MagnitudeNo).Flow.Result = []; %This needs to be cleared to avoid future errors.
end
nos = union(nos,nos); %Remove duplicates

%Loop over image stacks.
for no = nos
  if forcecopy || (SET(no).TSize>1)
    %Only upsample if timeresolved. Note linke stacks can have different
    %Number of timeframes (one if master have many as for instance in T2*)

    %--- Upsample volume
    SET(no).IM = upsampletemporal(f,SET(no).IM,immethod);
    newtsz = size(SET(no).IM,3);
    oldtsize = SET(no).TSize;
    SET(no).TSize = newtsz;

    %SET(no).TIncr = SET(no).TIncr/f; %Incorrect!
    if forcecopy || isequal(SET(no).TIncr,0)
      SET(no).TIncr = 1/newtsz;
      %TimeVector is already fixed above
    else
      SET(no).TIncr = SET(no).TIncr*(oldtsize)/(newtsz);
      SET(no).TimeVector = SET(no).TIncr*(0:SET(no).TSize-1);
    end

    %--- Upsample volume
    %--- Upsample contour
    if ~isempty(SET(NO).EndoX)
      SET(no).EndoX = upsampletemporal(f,SET(no).EndoX);
      SET(no).EndoY = upsampletemporal(f,SET(no).EndoY);
    end

    if ~isempty(SET(NO).EpiX)
      SET(no).EpiX = upsampletemporal(f,SET(no).EpiX);
      SET(no).EpiY = upsampletemporal(f,SET(no).EpiY);
    end

    if ~isempty(SET(NO).RVEndoX)
      SET(no).RVEndoX = upsampletemporal(f,SET(no).RVEndoX);
      SET(no).RVEndoY = upsampletemporal(f,SET(no).RVEndoY);
    end

    if ~isempty(SET(NO).RVEpiX)
      SET(no).RVEpiX = upsampletemporal(f,SET(no).RVEpiX);
      SET(no).RVEpiY = upsampletemporal(f,SET(no).RVEpiY);
    end

    % The condition "no EndoXView/EpiXView" may be represented by a single NaN.
    if ~isempty(SET(NO).EndoXView) ...
        && ~(length(SET(NO).EndoXView) == 1 && all(isnan(SET(NO).EndoXView)))
      SET(no).EndoXView = upsampletemporal(f,SET(no).EndoXView);
      SET(no).EndoYView = upsampletemporal(f,SET(no).EndoYView);
    end

    if ~isempty(SET(NO).EpiXView) ...
        && ~(length(SET(NO).EpiXView) == 1 && all(isnan(SET(NO).EpiXView)))
      SET(no).EpiXView = upsampletemporal(f,SET(no).EpiXView);
      SET(no).EpiYView = upsampletemporal(f,SET(no).EpiYView);
    end

    if ~isempty(SET(no).LA) && ~isempty(SET(no).LA.X)
      SET(no).LA.X = upsampletemporal(f,SET(no).LA.X);
      SET(no).LA.Y = upsampletemporal(f,SET(no).LA.Y);
    end

    if ~isempty(SET(no).RA) && ~isempty(SET(no).RA.X)
      SET(no).RA.X = upsampletemporal(f,SET(no).RA.X);
      SET(no).RA.Y = upsampletemporal(f,SET(no).RA.Y);
    end

    SET(no).LVV = upsampletemporal(f,SET(no).LVV,'vector');
    SET(no).PV = upsampletemporal(f,SET(no).PV,'vector');
    SET(no).EPV = upsampletemporal(f,SET(no).EPV,'vector');
    SET(no).LVM = upsampletemporal(f,SET(no).LVM,'vector');
    SET(no).RVV = upsampletemporal(f,SET(no).RVV,'vector');
    SET(no).RVEPV = upsampletemporal(f,SET(no).RVEPV,'vector');
    SET(no).RVM = upsampletemporal(f,SET(no).RVM,'vector');

    %SET(no).LVV = interp1(0:1/(length(SET(no).LVV)-1):1,SET(no).LVV,...
    %  0:1/((SET(no).TSize)-1):1);

    %Reset pins
    SET(no).EndoPinX = [];
    SET(no).EndoPinY = [];
    SET(no).EpiPinX = [];
    SET(no).EpiPinY = [];
    SET(no).RVEndoPinX = [];
    SET(no).RVEndoPinY = [];
    SET(no).RVEpiPinX = [];
    SET(no).RVEpiPinY = [];

    %Reset Interp Pts
    SET(no).EndoInterpX = [];
    SET(no).EndoInterpY = [];
    SET(no).EpiInterpX = [];
    SET(no).EpiInterpY = [];
    SET(no).RVEndoInterpX = [];
    SET(no).RVEndoInterpY = [];
    SET(no).RVEpiInterpX = [];
    SET(no).RVEpiInterpY = [];
    if ~isempty(SET(no).LA)
      SET(no).LA.InterpX = [];
      SET(no).LA.InterpY = [];
    end
    if ~isempty(SET(no).RA)
      SET(no).RA.InterpX = [];
      SET(no).RA.InterpY = [];
    end

    SET(no).EndoDraged = zeros(SET(no).TSize,SET(no).ZSize);
    SET(no).EpiDraged = zeros(SET(no).TSize,SET(no).ZSize);

    %--- Reset edge detection
    DATA.EndoEDGE0 = [];
    DATA.EndoEDGE1 = [];
    DATA.EndoEDGE2 = [];
    DATA.EndoEDGE3 = [];
    DATA.EpiEDGE0 = [];
    DATA.EpiEDGE1 = [];
    DATA.EpiEDGE2 = [];
    DATA.EpiEDGE3 = [];
    DATA.EndoEdgeDetected = false;
    DATA.EpiEdgeDetected = false;

    %Scar
    if not(isempty(SET(no).Scar))
      SET(no).Scar.Auto = upsampletemporal(f,single(SET(no).Scar.Auto)) > 0.5;
      SET(no).Scar.Result = upsampletemporal(f,single(SET(no).Scar.Result)) > 0.5;
      SET(no).Scar.Manual = int8(upsampletemporal(f,single(SET(no).Scar.Manual)));
      SET(no).Scar.Manual(SET(no).Scar.Manual < 0) = int8(-1);
      SET(no).Scar.MyocardMask = upsampletemporal(f,SET(no).Scar.MyocardMask);
      SET(no).Scar.NoReflow = upsampletemporal(f,SET(no).Scar.NoReflow);
      SET(no).Scar.IM = upsampletemporal(f,SET(no).Scar.IM);
    end

    %MaR
    if not(isempty(SET(no).MaR))
      SET(no).MaR.Auto = upsampletemporal(f,single(SET(no).MaR.Auto)) > 0.5;
      SET(no).MaR.Result = upsampletemporal(f,single(SET(no).MaR.Result)) > 0.5;
      SET(no).MaR.Manual = int8(upsampletemporal(f,single(SET(no).MaR.Manual)));
      SET(no).MaR.Manual(SET(no).MaR.Manual<0) = int8(-1);
      SET(no).MaR.MyocardMask = upsampletemporal(f,SET(no).MaR.MyocardMask);
    end

    %Reset zoomstate
    SET(no).NormalZoomState = [];
    SET(no).MontageZoomState = [];
    SET(no).MontageRowZoomState = [];

    %ROI's
    for rloop=1:SET(no).RoiN
      SET(no).Roi(rloop).X = upsampletemporal(f,SET(no).Roi(rloop).X);
      SET(no).Roi(rloop).Y = upsampletemporal(f,SET(no).Roi(rloop).Y);
      SET(no).Roi(rloop).T = 1:SET(no).TSize;
    end

    %make sure not out of bounds.
    %SET(no).CurrentSlice = 1;
    %SET(no).StartSlice = 1;
    %SET(no).EndSlice = 1;
    SET(no).CurrentTimeFrame = 1;
    SET(no).EndAnalysis = SET(no).TSize;
    SET(no).EST = round(SET(no).EST * f);
    SET(no).EDT = round(SET(no).EDT * f);

    if SET(no).EST > SET(no).TSize
      SET(no).EST = SET(no).TSize;
    end

    if SET(no).EDT > SET(no).TSize
      SET(no).EDT = SET(no).TSize;
    end

    if SET(no).EST == 0
      SET(no).EST = 1;
    end

    if SET(no).EDT == 0
      SET(no).EDT = 1;
    end
    %SET(no).EST = min([1 max([SET(no).TSize round(SET(no).EST*f)])]);
    %SET(no).EDT = min([1 max([SET(no).TSize round(SET(no).EDT*f)])]);
    SET(no).CurrentTimeFrame = min([1 max([SET(no).TSize round(SET(no).CurrentTimeFrame*f)])]);

    %remove Strain tagging analysis
    if not(isempty(SET(no).StrainTagging))
      logdisp('Upsample image stack clears the Strain tagging analysis',true);
      if ~isempty(DATA.GUI.StrainTagging)
        straintagging.straintagging('close_Callback');
      end
      SET(no).StrainTagging = [];
    end

    %remove Strain MITT analysis
    logmsg = 'Upsample image stack clears the Strain MITT analysis';
    strainmitt.strainmitt('removestrainanalysis',no,logmsg);

  end %If Tsize(loop)>1
end %Loop over image stacks

segment('updatevolume');
segment('updatetimebaraxes');
DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel)

%----------------------------------------
function newvol = upsampletemporal(f, vol, varargin)
%----------------------------------------
%Helper function to upsample in time. f is factor and vol is volume to
%upsample. Third argument is optional and is type;
% (image, segmentation (default), and vector

% Erik Bergvall
% Mainly rewritten Einar Heiberg

if abs(f - 1) < sqrt(eps) % resample factor is equal to one
  newvol = vol;
  return
end

% Check type of input. Default to 'segmentation' mode. (type 2)
if nargin <= 2
  type = 'segmentation';
else
  type = varargin{1};
end

switch type
  case 'image'
    newvol= upsamplevolume2([1,1,f,1],vol);
    %     N = size(vol,3);
    %
    %     if N>1
    %       %Volume is already time resolved
    %       tnew = linspace(0,N,round(f*N)+1); tnew(end) = [];
    %
    %       newvol = zeros(size(vol,1),size(vol,2),length(tnew),size(vol,4));
    %       hw = waitbar(0,'Please wait.');
    %       for i=1:size(newvol,3)
    %         g = 0;
    %         kfac = 6;
    %         for j=-kfac*size(vol,3):((kfac+1)*size(vol,3)+1)
    %           tj = j-1;
    %           if abs(tnew(i)-tj)<sqrt(eps)
    %             q = 1;
    %           else
    %             q = sin(pi*(tnew(i)-tj))/(pi*(tnew(i)-tj));
    %           end
    %           jmod = mod(j-1,size(vol,3))+1;
    %           g = g + q*vol(:,:,jmod,:);
    %         end
    %         newvol(:,:,i,:) = g;
    %
    %         waitbar(i/size(newvol,3),hw);
    %       end
    %       close(hw);
    %
    %       if not(isempty(vol))
    %         newvol(:,:,1,:) = vol(:,:,1,:);
    %       else
    %         newvol = [];
    %       end
    %     else
    %       %Volume is not timeresolved duplicate
    %       newvol = repmat(vol,[1 1 round(f) 1]);
    %     end;
    %
  case 'image-fourier'
    % Fourier upsampling of image - useful for flow upsampling

    % Allocate new data
    N = size(vol,3);
    Nnew = round(f*N);
    newvolfft = zeros(size(vol,1),size(vol,2),Nnew,size(vol,4));

    % do FFT and zero-pad (insert in larger allocated block)
    volfft = fft(vol,[],3);
    newvolfft(:,:,1:(N/2+1),:) = volfft(:,:,1:(N/2+1),:);
    newvolfft(:,:,(Nnew-N/2+2):Nnew,:) = volfft(:,:,(N/2+2):N,:);

    % Transform back
    newvol = real(ifft(newvolfft,[],3));
    newvol = newvol*f; % Fourier gain

  case 'segmentation'
    %curve
    % resampling of curves, uses nearest interp as we appearently
    % have no temporal point-2-point correspondence

    [X,Y] = meshgrid(0:1/(size(vol,2)-1):1,0:1/(size(vol,1)-1):1);
    [XI,YI] = meshgrid(0:1/(size(vol,2)*f-1):1,0:1/(size(vol,1)-1):1);
    newvol = zeros(size(XI,1),size(XI,2),size(vol,3));
    if (numel(varargin)>1) && (nnz(strcmp(varargin{2}, {'nearest','spline','pchip'}))>0)
      interptype = varargin{2};
    else
      %default
      interptype = 'linear';
    end
    if ~isequal(size(vol,2),1)
      %more than one temporal dimension => interpolate
      for loop = 1:size(vol,3)
        ZI = interp2(X,Y,vol(:,:,loop),XI,YI,interptype);
        newvol(:,:,loop) = ZI;
      end
    else
      %Simply duplicate
      newvol = repmat(vol,[1 round(f) 1]);
    end
  case 'vector'
    if isequal(length(vol),1)
      newvol = repmat(vol,1,round(f));
    else
      newvol = interp1(1:length(vol),vol,linspace(1,length(vol),round(f*length(vol))));
    end
end

%-------------------------------------------
function upsampleslices_Callback(f,doupdate)
%-------------------------------------------
%Upsample current image stack in slice direction with factor f. Takes care of
%segmentation. doupdate is default on, and if on then graphical update is
%performed

global DATA SET NO

if nargin < 2
  doupdate = true;
end

if nargin == 0
  str1 = dprintf('Upsampling/downsampling slices is irreversible.');
  str2 = dprintf('Are you sure?');
  str = sprintf('%s %s',str1,str2);
  if not(yesno(str,[],DATA.GUI.Segment))
    return;
  end

  %Get factor
  [f,ok] = mygetnumber('Enter upsampling/downsampling factor (<1 => downsample)','Resampling Factor',2,0,5);
  if ~ok
    myfailed('Invalid upsampling/downsampling factor or aborted.',DATA.GUI.Segment);
    return;
  end
end

if f == 1
  return
end

%Find image stacks to crop
nos = SET(NO).Linked; %Crop current image stack

oldsize = [SET(NO).XSize SET(NO).YSize SET(NO).ZSize];

for no = nos

  %Restore image enhancement
  if isfield(SET(no),'IMBackup') && (~isempty(SET(no).IMBackup))
    SET(no).IM = SET(no).IMBackup;
    SET(no).IMBackup = [];
  end

  %--- Upsample volume
  %SET(no).IM = upsamplevolume2([1,1,1,f],SET(no).IM,true);
  SET(no).IM = upsampleslices(f,SET(no).IM,false);

  SET(no).ZSize = size(SET(no).IM,4);

  %--- Upsample contour
  if ~isempty(SET(NO).EndoX)
    SET(no).EndoX = upsampleslices(f,SET(no).EndoX);
    SET(no).EndoY = upsampleslices(f,SET(no).EndoY);
  end

  if ~isempty(SET(NO).EpiX)
    SET(no).EpiX = upsampleslices(f,SET(no).EpiX);
    SET(no).EpiY = upsampleslices(f,SET(no).EpiY);
  end

  if ~isempty(SET(NO).RVEndoX)
    SET(no).RVEndoX = upsampleslices(f,SET(no).RVEndoX);
    SET(no).RVEndoY = upsampleslices(f,SET(no).RVEndoY);
  end

  if ~isempty(SET(NO).RVEpiX)
    SET(no).RVEpiX = upsampleslices(f,SET(no).RVEpiX);
    SET(no).RVEpiY = upsampleslices(f,SET(no).RVEpiY);
  end

  %Point
  if ~isempty(SET(no).Point)
    SET(no).Point.Z=round(resamplehelper(f,SET(no).Point.Z));
  end

  %CurrentSlice
  SET(no).CurrentSlice=min(max(1,round(resamplehelper(f,SET(no).CurrentSlice))),SET(no).ZSize);

  %Reset pins
  SET(no).EndoPinX = [];
  SET(no).EndoPinY = [];
  SET(no).EpiPinX = [];
  SET(no).EpiPinY = [];
  SET(no).RVEndoPinX = [];
  SET(no).RVEndoPinY = [];
  SET(no).RVEpiPinX = [];
  SET(no).RVEpiPinY = [];
  %Reset Interp Pts
  SET(no).EndoInterpX = [];
  SET(no).EndoInterpY = [];
  SET(no).EpiInterpX = [];
  SET(no).EpiInterpY = [];
  SET(no).RVEndoInterpX = [];
  SET(no).RVEndoInterpY = [];
  SET(no).RVEpiInterpX = [];
  SET(no).RVEpiInterpY = [];

  SET(no).EndoDraged = zeros(SET(no).TSize,SET(no).ZSize);
  SET(no).EpiDraged = zeros(SET(no).TSize,SET(no).ZSize);

  SET(no).SliceThickness = SET(no).SliceThickness/f;
  SET(no).SliceGap = SET(no).SliceGap/f;

  %--- Reset edge detection
  DATA.EndoEDGE0 = [];
  DATA.EndoEDGE1 = [];
  DATA.EndoEDGE2 = [];
  DATA.EndoEDGE3 = [];
  DATA.EpiEDGE0 = [];
  DATA.EpiEDGE1 = [];
  DATA.EpiEDGE2 = [];
  DATA.EpiEDGE3 = [];
  DATA.EndoEdgeDetected = false;
  DATA.EpiEdgeDetected = false;

  %Scar
  if not(isempty(SET(no).Scar))
    SET(no).Scar.Auto = upsampleslices(f,single(SET(no).Scar.Auto))>0.5;
    SET(no).Scar.Result = upsampleslices(f,single(SET(no).Scar.Result))>0.5;
    SET(no).Scar.Manual = int8(upsampleslices(f,single(SET(no).Scar.Manual)));
    SET(no).Scar.Manual(SET(no).Scar.Manual<0) = int8(-1);
    SET(no).Scar.MyocardMask = upsampleslices(f,SET(no).Scar.MyocardMask);
    SET(no).Scar.NoReflow = upsampleslices(f,SET(no).Scar.NoReflow);
    SET(no).Scar.IM = upsampleslices(f,SET(no).Scar.IM);
  end

  %MaR
  if not(isempty(SET(no).MaR))
    SET(no).MaR.Auto = upsampleslices(f,single(SET(no).MaR.Auto))>0.5;
    SET(no).MaR.Result = upsampleslices(f,single(SET(no).MaR.Result))>0.5;
    SET(no).MaR.Manual = int8(upsampleslices(f,single(SET(no).MaR.Manual)));
    SET(no).MaR.Manual(SET(no).MaR.Manual<0) = int8(-1);
    SET(no).MaR.MyocardMask = upsampleslices(f,SET(no).MaR.MyocardMask);
  end

  %Reset zoomstate
  SET(no).NormalZoomState = [];
  SET(no).MontageZoomState = [];
  SET(no).MontageRowZoomState = [];

  %ROI's
  for rloop = 1:SET(no).RoiN
    SET(no).Roi(rloop).Z = min(max(1,round(f*SET(no).Roi(rloop).Z)),SET(no).ZSize);
  end

  %Resample .LevelSet
  if ~isempty(SET(no).LevelSet)

    O = SET(no).LevelSet.Object;

    %resample stored objects, important to resample stored object    
    for obj = 1:O.getnumobjects()

      if O.iscompressed(obj)
        O.uncompress(obj,oldsize);
      end

      temp = O.getbw(obj);

      if ~isempty(temp)
        temp = upsampleslices(f,temp);
        O.setbw(obj,temp);
        O.setcut(obj,[]);
        O.setdeform(obj,[]);
        if ~O.getselected(obj)
          O.compress(obj); %Compress as we want to save memory in this step if there are many objects
        end
      end
    end

    %---things to change

    %Zoomstate
    SET(no).LevelSet.View.RZoomState=[];
    SET(no).LevelSet.View.GZoomState=[];
    SET(no).LevelSet.View.BZoomState=[];

    %slice
    if ~isfield(DATA.LevelSet,'imageorientation')
      DATA.LevelSet.imageorientation = 'transversal';
    end

    switch lower(DATA.LevelSet.imageorientation)
      case 'transversal'
        field = 'RSlice';
        slice = SET(no).LevelSet.View.RSlice;
      case 'sagittal'
        field = 'GSlice';
        slice = SET(no).LevelSet.View.GSlice;
      case 'coronal'
        field = 'BSlice';
        slice = SET(no).LevelSet.View.BSlice;
    end

    SET(no).LevelSet.View.(field)= max(1,min(round(resamplehelper(f,slice)),SET(no).ZSize));
    if doupdate
      segment3dp.tools('init3DP',1,1);
    end

  end

  %remove Strain tagging analysis
  if not(isempty(SET(no).StrainTagging))
    logdisp('Upsample image stack clears the Strain tagging analysis',true);
    if ~isempty(DATA.GUI.StrainTagging)
      straintagging.straintagging('close_Callback');
    end
    SET(no).StrainTagging = [];
  end
  %remove Strain MITT analysis
  logmsg = 'Upsample image stack clears the Strain MITT analysis';
  strainmitt.strainmitt('removestrainanalysis',no,logmsg);

  if isfield(SET(no),'RV') && ~isempty(SET(no).RV) && ~isempty(SET(no).RV.slicebasal)
    SET(no).RV.slicebasal = max(1,min(round(resamplehelper(f,SET(no).RV.slicebasal)),SET(no).ZSize));
    SET(no).RV.sliceapical = max(1,min(round(resamplehelper(f,SET(no).RV.sliceapical)),SET(no).ZSize));
  end

  %CT LV reconstructions points
  if ~isempty(SET(no).CT)&&isfield(SET(no).CT,'Apex')
    SET(no).CT.Apex(3)=min(size(SET(no).IM,4),max(1,round(resamplehelper(f,SET(no).CT.Apex(3)))));
    SET(no).CT.AV(3)=min(size(SET(no).IM,4),max(1,round(resamplehelper(f,SET(no).CT.AV(3)))));
    SET(no).CT.RVPulmi(3)=min(size(SET(no).IM,4),max(1,round(resamplehelper(f,SET(no).CT.RVPulmi(3)))));
  end

end %Loop over image stacks

%make sure not out of bounds.
SET(no).CurrentSlice = 1;
SET(no).StartSlice = 1;
SET(no).EndSlice = 1;

DATA.ViewIM{DATA.CurrentPanel} = [];

% DATA.ViewIM{DATA.CurrentPanel} = [];
% drawfunctions('drawimageno');
if doupdate
  viewfunctions('setview');
end

%----------------------------------
function setimageinfo_Callback(arg)
%----------------------------------
%Asks for and sets image details of current image stack

global DATA SET NO


if ~isa(DATA,'maingui')
  delete(gcbf)
  return
end

if not(DATA.DataLoaded)
  return;
end

gui = DATA.GUI.SetImageInfo;

if nargin < 2
  no = NO;
  if nargin < 1
    arg = 'init';
  end
end

switch arg
  case 'init'
    gui = mygui('setimageinfo.fig');
    DATA.GUI.SetImageInfo = gui;
    gui.no = no;
    setimageinfo_Callback('update');

  case 'update'
    %Populate edit boxes
    set(gui.handles.slicethicknessedit,'String',SET(gui.no).SliceThickness);
    set(gui.handles.slicegapedit,'String',SET(gui.no).SliceGap);
    set(gui.handles.resolutionxedit,'String',SET(gui.no).ResolutionX);
    set(gui.handles.resolutionyedit,'String',SET(gui.no).ResolutionY);
    set(gui.handles.timeincrementedit,'String',SET(gui.no).TIncr*1000);
    set(gui.handles.heartrateedit,'String',SET(gui.no).HeartRate);
    if SET(gui.no).Rotated
      set(gui.handles.rotatedyesradiobutton,'Value',1);
      set(gui.handles.rotatednoradiobutton,'Value',0);
    else
      set(gui.handles.rotatedyesradiobutton,'Value',0);
      set(gui.handles.rotatednoradiobutton,'Value',1);
    end


  case 'timeincrementedit'
    tincr = mygetedit(gui.handles.timeincrementedit);
    tincr = str2double(tincr)/1000;
    if not(isnan(tincr))
      heartrate = round(60/(tincr*SET(gui.no).TSize));
      set(gui.handles.heartrateedit,'String',heartrate);
    end
  case 'heartrateedit'
    hr = mygetedit(gui.handles.heartrateedit);
    hr = str2double(hr);
    if not(isnan(hr))
      tincr = (60/hr)/SET(gui.no).TSize;
      set(gui.handles.timeincrementedit,'String',tincr*1000);
    end
  case 'rotatedyes'
    set(gui.handles.rotatedyesradiobutton,'Value',1);
    set(gui.handles.rotatednoradiobutton,'Value',0);
  case 'rotatedno'
    set(gui.handles.rotatedyesradiobutton,'Value',0);
    set(gui.handles.rotatednoradiobutton,'Value',1);
  case 'apply'
    %Check input
    ok = true;
    if str2double(get(gui.handles.slicethicknessedit,'String'))<=0
      myfailed('Slice thickness must be positive.',DATA.GUI.Segment);
      ok = false;
    end
    if str2double(get(gui.handles.resolutionxedit,'String'))<=0
      myfailed('X Resolution must be positive.',DATA.GUI.Segment);
      ok = false;
    end
    if str2double(get(gui.handles.resolutionyedit,'String'))<=0
      myfailed('Y Resolution must be positive.',DATA.GUI.Segment);
      ok = false;
    end
    if str2double(get(gui.handles.timeincrementedit,'String'))<0
      myfailed('Time increment must be positive.',DATA.GUI.Segment);
      ok = false;
    end
    if str2double(get(gui.handles.heartrateedit,'String'))<0
      myfailed('Heart rate must be positive.',DATA.GUI.Segment);
      ok = false;
    end

    if ok
      if ~isempty(SET(gui.no).Flow) && (SET(gui.no).Flow.MagnitudeNo==gui.no || SET(gui.no).Flow.PhaseNo==gui.no)
        nos = [SET(gui.no).Flow.MagnitudeNo, SET(gui.no).Flow.PhaseNo];
      else
        nos = gui.no;
      end

      %For Flow images these should match therefore if you update one the
      %other is also adjusted.---------------
      for no2 = nos
        SET(no2).SliceThickness = str2double(get(gui.handles.slicethicknessedit,'String'));
        SET(no2).SliceGap =str2double(get(gui.handles.slicegapedit,'String'));
        SET(no2).ResolutionX = str2double(get(gui.handles.resolutionxedit,'String'));
        SET(no2).ResolutionY = str2double(get(gui.handles.resolutionyedit,'String'));
        SET(no2).TIncr = str2double(get(gui.handles.timeincrementedit,'String'))/1000;
        SET(no2).TimeVector = linspace(0,(SET(gui.no).TSize-1)*SET(gui.no).TIncr,SET(gui.no).TSize);
        SET(no2).BeatTime = SET(gui.no).TIncr*SET(gui.no).TSize;
        hr = str2double(get(gui.handles.heartrateedit,'String'));
        SET(no2).HeartRate = hr;
        if not(isempty(SET(no2).Flow)) && isfield(SET(no2).Flow,'HeartRate')
          SET(no2).Flow.HeartRate = hr;
        end
        if get(gui.handles.rotatedyesradiobutton,'Value')
          SET(no2).Rotated = true;
        else
          SET(no2).Rotated = false;
        end
      end
      %--------------

      SET(gui.no).NormalZoomState = [];
      SET(gui.no).MontageZoomState = [];
      SET(gui.no).MontageRowZoomState = [];
      viewfunctions('setview')

      if length(nos)>1
        mymsgbox(dprintf('Information applied to image stacks %d, %d.',nos(1),nos(2)),'Done!',DATA.GUI.Segment);
      else
        mymsgbox('Information applied to current image stack.','Done!',DATA.GUI.Segment);
      end
    else
      return;
    end
  case 'cancel'
    mymsgbox('No image info updated.');
    try
      DATA.GUI.SetImageInfo = close(DATA.GUI.SetImageInfo);
    catch %#ok<CTCH>
      delete(gui.fig);
    end
  otherwise
    myfailed(dprintf('Unknown argument %s to Set image info',arg),DATA.GUI.Segment);
    return;
end

%--------------------------------------
function removeallbutedes_Callback
%--------------------------------------
%Remove all timeframes in current image stack except diastole and systole.
global DATA SET NO

if ~DATA.Silent
  str1 = dprintf('Deleting all timeframes but diastole and systole (irreversible).');
  str2 = dprintf('Are you sure?');
  str = sprintf('%s %s',str1,str2);
  if ~yesno(str,[],DATA.GUI.Segment)
    return;
  end
end

ind = false(1,SET(NO).TSize);
ind(SET(NO).EDT) = true;
ind(SET(NO).EST) = true;
removetimeframes(ind);

%----------------------------------------------
function removeallbutthistimeframe_Callback
%----------------------------------------------
%Remove all timeframes except current timeframe.
global DATA SET NO

if ~DATA.Silent
  str1 = dprintf('Deleting all but this timeframe (irreversible).');
  str2 = dprintf('Are you sure?');
  str = sprintf('%s %s',str1,str2);
  if ~yesno(str,[],DATA.GUI.Segment)
    return;
  end
end

ind = ((1:SET(NO).TSize)==SET(NO).CurrentTimeFrame);
removetimeframes(ind);

%--------------------------------------------
function removecurrenttimeframe_Callback
%--------------------------------------------
%Remove current timeframe
global DATA SET NO
no = NO;
currenttf = SET(NO).CurrentTimeFrame;
str1 = dprintf('Deleting current timeframe (irreversible).');
str2 = sprintf('%s %d, %s #%d',dprintf('Time frame'),currenttf, dprintf('Stack'),no);
str3 = dprintf('Are you sure?');
str = sprintf('%s\n%s\n%s',str1,str2,str3);
if ~yesno(str,[],DATA.GUI.Segment)
  return;
end

ind = ((1:SET(no).TSize)~= currenttf);
removetimeframes(ind);


%------------------------------------------
function removenexttimeframes_Callback
%------------------------------------------
%Remove next timeframe and upto and including last timeframe
global DATA SET NO

if SET(NO).TSize==1
  myfailed('Only one timeframe exists. Can not remove current timeframe.',DATA.GUI.Segment);
  return;
end

if SET(NO).CurrentTimeFrame==SET(NO).TSize
  myfailed('Last timeframe. No more timeframes to delete.',DATA.GUI.Segment);
  return;
end

if ~yesno('Deleting next timeframe and all timeframes to last timeframe. Are you sure?',[],DATA.GUI.Segment)
  return;
end

ind = true(1,SET(NO).TSize);
ind((SET(NO).CurrentTimeFrame+1):end) = false;
removetimeframes(ind);

%----------------------------------------------
function removeprevioustimeframes_Callback
%----------------------------------------------
%Removes previous timeframe and all frames to first timeframe.
global DATA SET NO

if SET(NO).TSize==1
  myfailed('Only one timeframe exists. Can not remove current timeframe.',DATA.GUI.Segment);
  return;
end

if ~yesno('Deleting timeframes from start up to previous timeframe. Are you sure?',[],DATA.GUI.Segment)
  return;
end

ind = true(1,SET(NO).TSize);
ind(1:(SET(NO).CurrentTimeFrame-1)) = false;
removetimeframes(ind);

%---------------------------------
function removetimeframes(ind)
%---------------------------------
%Helper function to remove timeframes. This is the workhorse
%when removing timeframes

global DATA SET NO

disableundo;

%Generic deletions
if DATA.EndoEdgeDetected
  DATA.EndoEDGE0 = DATA.EndoEDGE0(:,:,ind,:);
  DATA.EndoEDGE1 = DATA.EndoEDGE1(:,:,ind,:);
  DATA.EndoEDGE2 = DATA.EndoEDGE2(:,:,ind,:);
  DATA.EndoEDGE3 = DATA.EndoEDGE3(:,:,ind,:);
end
if DATA.EpiEdgeDetected
  DATA.EpiEDGE0 = DATA.EpiEDGE0(:,:,ind,:);
  DATA.EpiEDGE1 = DATA.EpiEDGE1(:,:,ind,:);
  DATA.EpiEDGE2 = DATA.EpiEDGE2(:,:,ind,:);
  DATA.EpiEDGE3 = DATA.EpiEDGE3(:,:,ind,:);
end

nop = SET(NO).Linked;

for loop = 1:length(DATA.ViewIM)
  DATA.ViewIM{loop} = [];
end

%Loop over imagestacks
for nloop=1:length(nop)

  no = nop(nloop);

  if SET(no).TSize>1
    %Remove from variables

    %IM
    SET(no).IM = SET(no).IM(:,:,ind,:);

    %Timevector
    indexes=find(ind);
    if indexes(1)>1
      SET(no).TDelay=SET(no).TimeVector(indexes(1));
      SET(no).TimeVector = SET(no).TimeVector(ind)-SET(no).TDelay;
    else
      SET(no).TimeVector = SET(no).TimeVector(ind);
    end

    if isfield(SET(no),'TriggerTime') && numel(SET(no).TriggerTime)>1
      %TriggerTime
      SET(no).TriggerTime=SET(no).TriggerTime(:,ind);
    end

    if numel(SET(no).InversionTime)>1
      %InversionTime
      SET(no).InversionTime = SET(no).InversionTime(:,ind);
    end

    if numel(SET(no).EchoTime)>1
      %EchoTime
      SET(no).EchoTime = SET(no).EchoTime(:,ind);
    end

    if isfield(SET(no),'T2preptime') &&numel(SET(no).T2preptime)>1
      %T2-Preparation Time:
      SET(no).T2preptime = SET(no).T2preptime(:,ind);
    end

    %PapillaryIM
    if not(isempty(SET(no).PapillaryIM))
      SET(no).PapillaryIM = SET(no).PapillaryIM(:,:,ind,:);
    end

  end

  if SET(no).TSize>1
    newEST=find(find(ind)==SET(no).EST);
    newEDT=find(find(ind)==SET(no).EDT);
    if isempty(newEST)||isempty(newEDT)
      logdisp('Systolic or diastolic timeframe deleted. Autodetecting.',true);
      [~,SET(no).EST] = min(SET(no).LVV(ind));
      [~,SET(no).EDT] = max(SET(no).LVV(ind));
    else
      SET(no).EST=newEST;
      SET(no).EDT=newEDT;
    end

    SET(no).LVV = SET(no).LVV(ind);
    SET(no).EPV = SET(no).EPV(ind);
    SET(no).PV = SET(no).PV(ind);
    SET(no).RVV = SET(no).RVV(ind);
    SET(no).RVEPV = SET(no).RVEPV(ind);

    if ~isempty(SET(no).EndoX)
      SET(no).EndoX = SET(no).EndoX(:,ind,:);
      SET(no).EndoY = SET(no).EndoY(:,ind,:);
    end
    if ~isempty(SET(no).EpiX)
      SET(no).EpiX = SET(no).EpiX(:,ind,:);
      SET(no).EpiY = SET(no).EpiY(:,ind,:);
    end
    if ~isempty(SET(no).RVEndoX)
      SET(no).RVEndoX = SET(no).RVEndoX(:,ind,:);
      SET(no).RVEndoY = SET(no).RVEndoY(:,ind,:);
    end
    if ~isempty(SET(no).RVEpiX)
      SET(no).RVEpiX = SET(no).RVEpiX(:,ind,:);
      SET(no).RVEpiY = SET(no).RVEpiY(:,ind,:);
    end

    if ~isempty(SET(no).EndoPinX)
      SET(no).EndoPinX = SET(no).EndoPinX(ind,:);
      SET(no).EndoPinY = SET(no).EndoPinY(ind,:);
    end
    if ~isempty(SET(no).EpiPinX)
      SET(no).EpiPinX = SET(no).EpiPinX(ind,:);
      SET(no).EpiPinY = SET(no).EpiPinY(ind,:);
    end
    if ~isempty(SET(no).RVEndoPinX)
      SET(no).RVEndoPinX = SET(no).RVEndoPinX(ind,:);
      SET(no).RVEndoPinY = SET(no).RVEndoPinY(ind,:);
    end
    if ~isempty(SET(no).RVEpiPinX)
      SET(no).RVEpiPinX = SET(no).RVEpiPinX(ind,:);
      SET(no).RVEpiPinY = SET(no).RVEpiPinY(ind,:);
    end

    if ~isempty(SET(no).EndoInterpX)
      SET(no).EndoInterpX = SET(no).EndoInterpX(ind,:);
      SET(no).EndoInterpY = SET(no).EndoInterpY(ind,:);
    end
    if ~isempty(SET(no).EpiInterpX)
      SET(no).EpiInterpX = SET(no).EpiInterpX(ind,:);
      SET(no).EpiInterpY = SET(no).EpiInterpY(ind,:);
    end
    if ~isempty(SET(no).RVEndoInterpX)
      SET(no).RVEndoInterpX = SET(no).RVEndoInterpX(ind,:);
      SET(no).RVEndoInterpY = SET(no).RVEndoInterpY(ind,:);
    end
    if ~isempty(SET(no).RVEpiInterpX)
      SET(no).RVEpiInterpX = SET(no).RVEpiInterpX(ind,:);
      SET(no).RVEpiInterpY = SET(no).RVEpiInterpY(ind,:);
    end
    if ~isempty(SET(no).EndoDraged)
      SET(no).EndoDraged = SET(no).EndoDraged(ind,:);
    end
    if ~isempty(SET(no).EpiDraged)
      SET(no).EpiDraged = SET(no).EpiDraged(ind,:);
    end

    if ~isempty(SET(no).LA)
      % calls clpen class method
      SET(no).LA = SET(no).LA.deletetimeframes(ind);
    end

    if ~isempty(SET(no).RA)
      % calls clpen class method
      SET(no).RA = SET(no).RA.deletetimeframes(ind);
    end

    if SET(no).RoiN>0
      for rloop = 1:SET(no).RoiN
        SET(no).Roi(rloop).X = SET(no).Roi(rloop).X(:,ind);
        SET(no).Roi(rloop).Y = SET(no).Roi(rloop).Y(:,ind);
        SET(no).Roi(rloop).T = find(not(isnan(SET(no).Roi(rloop).X(1,:))));
      end
    end

    SET(no).TSize = size(SET(no).IM,3);

  else
    mywarning('Cannot delete current timeframe when there is only one timeframe.',DATA.GUI.Segment);
  end

  if SET(no).EndAnalysis>SET(no).TSize
    SET(no).EndAnalysis = SET(no).TSize;
  end

  if SET(no).StartAnalysis>SET(no).TSize
    SET(no).StartAnalysis = SET(no).TSize;
  end

  SET(no).CurrentTimeFrame = min(SET(no).CurrentTimeFrame,SET(no).TSize);

  %remove timeframes in .LevelSet
  if ~isempty(SET(no).LevelSet)

    %Zoomstate
    SET(no).LevelSet.View.RZoomState=[];
    SET(no).LevelSet.View.GZoomState=[];
    SET(no).LevelSet.View.BZoomState=[];

    %time slice
    SET(no).LevelSet.View.TSlice=round(SET(no).TSize/2);
  end

  %MaR
  if not(isempty(SET(no).MaR))
    SET(no).MaR.Auto = SET(no).MaR.Auto(:,:,ind,:);
    SET(no).MaR.Result = SET(no).MaR.Result(:,:,ind,:);
    SET(no).MaR.Manual = SET(no).MaR.Manual(:,:,ind,:);
    SET(no).MaR.MyocardMask = SET(no).MaR.MyocardMask(:,:,ind,:);
    SET(no).MaR.Percentage = SET(no).MaR.Percentage(ind);
  end

  % remove measures
  indi=find(ind);
  if ~isempty(SET(no).Measure)
    mind=[];
    for mloop=1:length(SET(no).Measure)
      if ismember(SET(no).Measure(mloop).T,indi)
        SET(no).Measure(mloop).T=find(SET(no).Measure(mloop).T==indi);
        mind = [mind mloop]; %#ok<AGROW>
      end
    end
    SET(no).Measure = SET(no).Measure(mind);
  end

  %remove points
  pind=[];
  if ~isempty(SET(no).Point.X)
    for ploop=1:length(SET(no).Point.X)
      if ismember(SET(no).Point.T(ploop),indi)
        SET(no).Point.T(ploop)=find(SET(no).Point.T(ploop)==indi);
        pind = [pind ploop]; %#ok<AGROW>
      elseif isnan(SET(no).Point.T(ploop))
        pind = [pind ploop]; %#ok<AGROW>
      end
    end
    SET(no).Point.X = SET(no).Point.X(pind);
    SET(no).Point.Y = SET(no).Point.Y(pind);
    SET(no).Point.T = SET(no).Point.T(pind);
    SET(no).Point.Z = SET(no).Point.Z(pind);
    SET(no).Point.Label = SET(no).Point.Label(pind);
  end

  %remove Strain tagging analysis
  if not(isempty(SET(no).StrainTagging))
    logdisp('Removing time frames clears the Strain tagging analysis',true);
    if ~isempty(DATA.GUI.StrainTagging)
      straintagging.straintagging('close_Callback');
    end
    SET(no).StrainTagging = [];
  end

  %remove Strain MITT analysis
  logmsg = 'Removing time frames clears the Strain MITT analysis';
  strainmitt.strainmitt('removestrainanalysis',no,logmsg);

end %Loop over all linked image stacks.

drawfunctions('drawpanel',DATA.CurrentPanel);
segment('updatevolume');
DATA.updatetimebaraxes;

if ~DATA.Silent
  mymsgbox(dprintf('Removed timeframe(s) %s.',num2str(find(not(ind)))),'',DATA.GUI.Segment);
end

%   With disabletimethings commented out, this does nothing,
%   so I commented out the whole if statement /JU
%if SET(NO).TSize==1
%  %disabletimethings;
%  drawimageno;
%end;

%-------------------------------------
function removeunselected_Callback
%-------------------------------------
%Remove unselected slices from current image stack
global DATA SET NO

%Remove unselected slices, i.e keep selected
ind = false(SET(NO).ZSize,1);
ind(SET(NO).StartSlice:SET(NO).EndSlice)=true;

if all(ind)
  myfailed('No slices selected.',DATA.GUI.Segment)
  return;
end

removeslices_Callback(ind);

%-----------------------------------
function removeselected_Callback
%-----------------------------------
%Remove selected image stacks from current image stack.

global DATA SET NO

%Remove selected slices
ind = true(SET(NO).ZSize,1);
ind(SET(NO).StartSlice:SET(NO).EndSlice) = false;

if all(ind)
  myfailed('No slices selected.',DATA.GUI.Segment)
  return;
end

removeslices_Callback(ind);


%-------------------------------------
function removesliceshelper(ind,newno)
%-------------------------------------
%This is the workhorse when removing slices ind is a logical input of which
%slices to keep.
global DATA SET NO

if nargin<2
  newno = NO;
end

oldsize = [SET(newno).XSize SET(newno).YSize SET(newno).ZSize];

%Used for ROI removal
oldzsize = SET(newno).ZSize;
oldslices = 1:oldzsize;

slicemap = oldslices;
slicemap(not(ind)) = NaN;
slicemap(ind) = 1:sum(ind);

if DATA.EndoEdgeDetected
  DATA.EndoEDGE0 = DATA.EndoEDGE0(:,:,:,ind);
  DATA.EndoEDGE1 = DATA.EndoEDGE1(:,:,:,ind);
  DATA.EndoEDGE2 = DATA.EndoEDGE2(:,:,:,ind);
  DATA.EndoEDGE3 = DATA.EndoEDGE3(:,:,:,ind);
end

if not(isempty(DATA.BALLOON))
  DATA.BALLOON = DATA.BALLOON(:,:,:,ind);
end

if DATA.EpiEdgeDetected
  DATA.EpiEDGE0 = DATA.EpiEDGE0(:,:,:,ind);
  DATA.EpiEDGE1 = DATA.EpiEDGE1(:,:,:,ind);
  DATA.EpiEDGE2 = DATA.EpiEDGE2(:,:,:,ind);
  DATA.EpiEDGE3 = DATA.EpiEDGE3(:,:,:,ind);
end

%Find image stacks to crop
nop = SET(newno).Linked; %Crop current and linked image stacks
if ~isempty(SET(newno).Parent)
  no = SET(newno).Parent;
else
  no = newno;
end

orignumslices = SET(no).ZSize;

if size(SET(no).InversionTime,1) == orignumslices
  % InversionTime
  SET(no).InversionTime = SET(no).InversionTime(ind,:);
end

if size(SET(no).EchoTime,1) == orignumslices
  % EchoTime
  SET(no).EchoTime = SET(no).EchoTime(ind,:);
end

if isfield(SET(no),'T2preptime') && size(SET(no).T2preptime,1) == orignumslices
  % T2-Preparation Time
  SET(no).T2preptime = SET(no).T2preptime(ind,:);
end

%---- Loop over image stacks for Flow-y data,
for nloop = 1:length(nop)

  nol=nop(nloop);

  %Restore image enhancement
  if isfield(SET(nol),'IMBackup') && (~isempty(SET(nol).IMBackup))
    SET(nol).IM = SET(nol).IMBackup;
    SET(nol).IMBackup = [];
  end

  %Remove from image
  SET(nol).IM = SET(nol).IM(:,:,:,ind);

  %Update image position
  if ~ind(1)
    % Removing the first slice.
    % If we remove N slices at the top of the stack (low indices), we need to shift
    % ImagePosition by N*slicedistance.
    firstremoved = find(ind, 1, 'first');
    slicestoshift = firstremoved - 1; % N

    zdir = cross(...
      SET(nol).ImageOrientation(1:3),...
      SET(nol).ImageOrientation(4:6));
    zdir = zdir(:)';

    slicedistance = SET(nol).SliceThickness + SET(nol).SliceGap;
    SET(nol).ImagePosition = SET(nol).ImagePosition - slicestoshift*zdir*slicedistance;
    SET(nol).RotatedImagePosition = SET(nol).RotatedImagePosition - slicestoshift*zdir*slicedistance;

  end

  %PapillaryIM
  if not(isempty(SET(nol).PapillaryIM))
    SET(nol).PapillaryIM = SET(nol).PapillaryIM(:,:,:,ind);
  end

end

%Scar
if not(isempty(SET(no).Scar))
  SET(no).Scar.IM   = SET(no).Scar.IM(:,:,ind);
  SET(no).Scar.Auto = SET(no).Scar.Auto(:,:,ind);
  SET(no).Scar.Result = SET(no).Scar.Result(:,:,ind);
  SET(no).Scar.Manual = SET(no).Scar.Manual(:,:,ind);
  %SET(no).Scar.Undo = SET(no).Scar.Undo(:,:,ind);
  SET(no).Scar.MyocardMask = SET(no).Scar.MyocardMask(:,:,ind);
  SET(no).Scar.NoReflow = SET(no).Scar.NoReflow(:,:,ind);
  if numel(SET(no).Scar.mthreshold) == SET(no).ZSize
    SET(no).Scar.mthreshold = SET(no).Scar.mthreshold(ind);
  end
end

%MaR
if not(isempty(SET(no).MaR))
  SET(no).MaR.Auto = SET(no).MaR.Auto(:,:,:,ind);
  SET(no).MaR.Result = SET(no).MaR.Result(:,:,:,ind);
  SET(no).MaR.Manual = SET(no).MaR.Manual(:,:,:,ind);
  SET(no).MaR.MyocardMask = SET(no).MaR.MyocardMask(:,:,:,ind);
end

%CT definition of slices included in LV
if isfield(SET(no),'CT') && isfield(SET(no).CT,'Apex')
  newstartslice = find(ind==1,1,'first');
  SET(no).CT.AV(3) = min(length(ind),max(1,SET(no).CT.AV(3)-newstartslice+1));
  SET(no).CT.Apex(3)= min(length(ind),max(1,SET(no).CT.Apex(3)-newstartslice+1));
end

%Update variables
SET(no).ZSize = sum(ind);
SET(no).StartSlice = 1;
SET(no).EndSlice = SET(no).ZSize;
SET(no).CurrentSlice = min(SET(no).CurrentSlice,SET(no).ZSize);

%Remove from contour
if ~isempty(SET(no).EndoX)
  SET(no).EndoX = SET(no).EndoX(:,:,ind);
  SET(no).EndoY = SET(no).EndoY(:,:,ind);
end

if ~isempty(SET(no).EpiX)
  SET(no).EpiX = SET(no).EpiX(:,:,ind);
  SET(no).EpiY = SET(no).EpiY(:,:,ind);
end

if ~isempty(SET(no).RVEndoX)
  SET(no).RVEndoX = SET(no).RVEndoX(:,:,ind);
  SET(no).RVEndoY = SET(no).RVEndoY(:,:,ind);
end

if ~isempty(SET(no).RVEpiX)
  SET(no).RVEpiX = SET(no).RVEpiX(:,:,ind);
  SET(no).RVEpiY = SET(no).RVEpiY(:,:,ind);
end

if ~isempty(SET(no).EndoPinX)
  SET(no).EndoPinX = SET(no).EndoPinX(:,ind);
  SET(no).EndoPinY = SET(no).EndoPinY(:,ind);
end
if ~isempty(SET(no).EpiPinX)
  SET(no).EpiPinX = SET(no).EpiPinX(:,ind);
  SET(no).EpiPinY = SET(no).EpiPinY(:,ind);
end
if ~isempty(SET(no).RVEndoPinX)
  SET(no).RVEndoPinX = SET(no).RVEndoPinX(:,ind);
  SET(no).RVEndoPinY = SET(no).RVEndoPinY(:,ind);
end
if ~isempty(SET(no).RVEpiPinX)
  SET(no).RVEpiPinX = SET(no).RVEpiPinX(:,ind);
  SET(no).RVEpiPinY = SET(no).RVEpiPinY(:,ind);
end

if ~isempty(SET(no).EndoInterpX)
  SET(no).EndoInterpX = SET(no).EndoInterpX(:,ind);
  SET(no).EndoInterpY = SET(no).EndoInterpY(:,ind);
end
if ~isempty(SET(no).EpiInterpX)
  SET(no).EpiInterpX = SET(no).EpiInterpX(:,ind);
  SET(no).EpiInterpY = SET(no).EpiInterpY(:,ind);
end
if ~isempty(SET(no).RVEndoInterpX)
  SET(no).RVEndoInterpX = SET(no).RVEndoInterpX(:,ind);
  SET(no).RVEndoInterpY = SET(no).RVEndoInterpY(:,ind);
end
if ~isempty(SET(no).RVEpiInterpX)
  SET(no).RVEpiInterpX = SET(no).RVEpiInterpX(:,ind);
  SET(no).RVEpiInterpY = SET(no).RVEpiInterpY(:,ind);
end

if ~isempty(SET(no).LA)
  % calls clpen class method
  SET(no).LA = SET(no).LA.deleteslices(ind);
end

if ~isempty(SET(no).RA)
  % calls clpen class method
  SET(no).RA = SET(no).RA.deleteslices(ind);
end

SET(no).EndoDraged = SET(no).EndoDraged(:,ind);
SET(no).EpiDraged = SET(no).EpiDraged(:,ind);

%Move ROI'S (if necessary)
roiind = zeros(1,SET(no).RoiN);
for loop=1:SET(no).RoiN
  if not(isnan(slicemap(SET(no).Roi(loop).Z)))
    SET(no).Roi(loop).Z = slicemap(SET(no).Roi(loop).Z);
    roiind(loop)=1;
  end
end
roiind=logical(roiind);

%Remove ROI's (if necessary)
SET(no).Roi = SET(no).Roi(roiind);
SET(no).RoiN = sum(roiind(:));
if ~all(roiind)
  SET(no).RoiCurrent = [];
end

if not(isempty(SET(no).Scar))
  segment('updatevolume');
  viability('viabilitycalcvolume');
end

%remove slices in .LevelSet
if ~isempty(SET(no).LevelSet) && isfield(SET(no).LevelSet,'Object')

  O = SET(no).LevelSet.Object;

  %remove slices in stored objects
  for obj=1:O.getnumobjects()

    %Uncompress 
    if O.iscompressed(obj)
      O.uncompress(obj,oldsize); %Special taking size into account
    end

    temp = O.getbw(obj);

    if ~isempty(temp)
      temp = temp(:,:,ind);
      O.setbw(obj,temp);
      O.setcut(obj,[]);
      O.setdeform(obj,[]);    
      if ~O.getselected(obj)
        O.compress(obj); %Compress
      end
    end

  end
  
  %---things to change

  %Zoomstate
  SET(no).LevelSet.View.RZoomState = [];
  SET(no).LevelSet.View.GZoomState = [];
  SET(no).LevelSet.View.BZoomState = [];

  %slice
  if SET(no).ZSize>1
    SET(no).LevelSet.View.RSlice=round(SET(no).ZSize/2);
  else
    SET(no).LevelSet.View.RSlice=round(SET(no).TSize/2);
  end
end

if ~isempty(SET(no).Point.X)
  pointind = not(isnan(slicemap(round(SET(no).Point.Z))));
  for loop=1:length(SET(no).Point.X)
    if not(isnan(slicemap(round(SET(no).Point.Z(loop)))))
      SET(no).Point.Z(loop) = slicemap(round(SET(no).Point.Z(loop)));
    end
  end

  %Remove points if necessary
  SET(no).Point.X = SET(no).Point.X(pointind);
  SET(no).Point.Y = SET(no).Point.Y(pointind);
  SET(no).Point.T = SET(no).Point.T(pointind);
  SET(no).Point.Z = SET(no).Point.Z(pointind);
  SET(no).Point.Label = SET(no).Point.Label(pointind);
end

% remove measures
if ~isempty(SET(no).Measure)
  mind=[];
  %   for mloop=1:length(SET(no).Measure)
  %     if sum(unique(slicemap(SET(no).Measure(mloop).Z))) >= 1 %ismember(SET(no).Measure(mloop).Z,slicemap)
  %       SET(no).Measure(mloop).Z=slicemap(SET(no).Measure(mloop).Z);
  %       mind = [mind mloop]; %#ok<AGROW>
  %     end
  %   end
  slicestoinclude = find(~isnan(slicemap));
  lowerslice =  slicestoinclude(1);
  upperslice = slicestoinclude(end);

  for mloop=1:length(SET(no).Measure)
    if all(SET(no).Measure(mloop).Z>=lowerslice) && all(SET(no).Measure(mloop).Z<=upperslice) ...
        %&& SET(no).Measure(mloop).Z(2)>=lowerslice && SET(no).Measure(mloop).Z(2)<=upperslice%sum(unique(slicemap(SET(no).Measure(mloop).Z))) >= 1 %ismember(SET(no).Measure(mloop).Z,slicemap)
      SET(no).Measure(mloop).Z=slicemap(SET(no).Measure(mloop).Z);
      mind = [mind mloop]; %#ok<AGROW>
    end
  end

  if ~isempty(mind)
    SET(no).Measure = SET(no).Measure(mind);
  else
    SET(no).Measure = [];
  end
end

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  logdisp('Removing slices clears the Strain tagging analysis',true);
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  SET(no).StrainTagging = [];
end
%remove Strain MITT analysis
logmsg = 'Removing slices clears the Strain MITT analysis';
strainmitt.strainmitt('removestrainanalysis',no,logmsg);


%------------------------------------
function removeslices_Callback(ind)
%------------------------------------
%Remove slices from current image stack.

global DATA SET NO

disableundo;

if nargin==0
  %Generate index
  ind = true(SET(NO).ZSize,1);
  ind(SET(NO).StartSlice:SET(NO).EndSlice)=false;
end

if ~any(ind)
  myfailed('If you do this there will be no slices left!',DATA.GUI.Segment);
  return;
end

if nargin~=2
  str1 = dprintf('Do you really want to remove the slices?');
  str2 = dprintf('This operation is irreversible.');
  str = sprintf('%s %s',str1,str2);
  if ~yesno(str,[],DATA.GUI.Segment)
    %DATA.failedaborted;
    return;
  end
end

if ~isempty(SET(NO).LevelSet) && length(ind)==1 && SET(NO).TSize==1
  if ~yesno('Image stack contains general segmentation tool data. This will be lost since only 2D image will be left after this operation . Are you sure?',[],DATA.GUI.Segment)
    %DATA.failedaborted;
    return;
  end
  SET(NO).LevelSet = [];
end

if ind(1)&&ind(end)&&(sum(~ind)>0)
  %take first and last and remove som in between => image position will
  %fail.
  if ~yesno('This operation will make the 3D data inconsistent. Do you want to proceed?',[],DATA.GUI.Segment)
    %DATA.failedaborted;
    return;
  end
end

removesliceshelper(ind); %Call to remove

segment('updatevolume');

DATA.ViewIM{DATA.CurrentPanel} = [];
if SET(NO).ZSize<20
  viewfunctions('setviewtype','montage')
else
  viewfunctions('setviewtype','one');
end

%--------------------------------------
function nos = crophelper(no0,xind,yind)
%--------------------------------------
%This function is the enginge in cropping

global SET DATA

oldsize = [SET(no0).XSize SET(no0).YSize SET(no0).ZSize];

%Find image stacks to crop
nos = SET(no0).Linked; %Crop current and linked image stacks

%---- Loop over image stacks to crop
for no = nos

  if SET(no).Rotated
    SET(no).RotationCenter = SET(no).RotationCenter-yind(1)-1;
  end

  %Restore image enhancement
  if isfield(SET(no),'IMBackup') && (~isempty(SET(no).IMBackup))
    SET(no).IM = SET(no).IMBackup;
    SET(no).IMBackup = [];
  end

  %ImageOrientation(1:3) is y in Segment

  %y direction
  SET(no).ImagePosition = SET(no).ImagePosition+(yind(1)-1)*SET(no).ResolutionY*SET(no).ImageOrientation(1:3);
  SET(no).RotatedImagePosition = SET(no).RotatedImagePosition+(yind(1)-1)*SET(no).ResolutionY*SET(no).RotatedImageOrientation(1:3);

  %x direction
  SET(no).ImagePosition = SET(no).ImagePosition+(xind(1)-1)*SET(no).ResolutionX*SET(no).ImageOrientation(4:6);
  SET(no).RotatedImagePosition = SET(no).RotatedImagePosition+(xind(1)-1)*SET(no).ResolutionX*SET(no).RotatedImageOrientation(4:6);

  %Remove from image
  SET(no).IM = SET(no).IM(xind,yind,:,:);
  oldxsize = SET(no).XSize; %Remember old size
  oldysize = SET(no).YSize;

  %Set new size
  SET(no).XSize = size(SET(no).IM,1);
  SET(no).YSize = size(SET(no).IM,2);

  if not(isempty(SET(no).PapillaryIM))
    SET(no).PapillaryIM = SET(no).PapillaryIM(xind,yind,:,:);
  end

  if (DATA.EndoEdgeDetected) && isequal(oldxsize,size(DATA.EndoEDGE0,1)) && isequal(oldysize,size(DATA.EndoEDGE0,2))
    DATA.EndoEDGE0 = DATA.EndoEDGE0(xind,yind,:,:);
    DATA.EndoEDGE1 = DATA.EndoEDGE1(xind,yind,:,:);
    DATA.EndoEDGE2 = DATA.EndoEDGE2(xind,yind,:,:);
    DATA.EndoEDGE3 = DATA.EndoEDGE3(xind,yind,:,:);
  end

  if not(isempty(DATA.BALLOON)) && isequal(oldxsize,size(DATA.BALLOON,1)) && isequal(oldysize,size(DATA.BALLOON,2))
    DATA.BALLOON = DATA.BALLOON(xind,yind,:,:);
  end

  if (DATA.EpiEdgeDetected) && isequal(oldxsize,size(DATA.EpiEDGE0,1)) && isequal(oldysize,size(DATA.EpiEDGE0,2))
    DATA.EpiEDGE0 = DATA.EpiEDGE0(xind,yind,:,:);
    DATA.EpiEDGE1 = DATA.EpiEDGE1(xind,yind,:,:);
    DATA.EpiEDGE2 = DATA.EpiEDGE2(xind,yind,:,:);
    DATA.EpiEDGE3 = DATA.EpiEDGE3(xind,yind,:,:);
  end

  if ~isempty(SET(no).Scar)
    SET(no).Scar.IM   = SET(no).Scar.IM(xind,yind,:);
    SET(no).Scar.Auto = SET(no).Scar.Auto(xind,yind,:);
    SET(no).Scar.Result = SET(no).Scar.Result(xind,yind,:);
    SET(no).Scar.Manual = SET(no).Scar.Manual(xind,yind,:);
    %SET(no).Scar.Undo = SET(no).Scar.Undo(xind,yind,:);
    SET(no).Scar.MyocardMask = SET(no).Scar.MyocardMask(xind,yind,:);
    SET(no).Scar.NoReflow = SET(no).Scar.NoReflow(xind,yind,:);
    if ~isempty(SET(no).Scar.GreyZone)
      if ~isempty(SET(no).Scar.GreyZone.map)
        SET(no).Scar.GreyZone.map = SET(no).Scar.GreyZone.map(xind,yind,:,:,:);
      end
    end
  end

  if ~isempty(SET(no).MaR)
    SET(no).MaR.Auto = SET(no).MaR.Auto(xind,yind,:,:);
    SET(no).MaR.Result = SET(no).MaR.Result(xind,yind,:,:);
    SET(no).MaR.Manual = SET(no).MaR.Manual(xind,yind,:,:);
    SET(no).MaR.MyocardMask = SET(no).MaR.MyocardMask(xind,yind,:,:);
  end

  %Translate segmentation
  xofs = -xind(1)+1;
  yofs = -yind(1)+1;

  if ~isempty(SET(no).CT )
    if isfield(SET(no).CT,'Apex')
      SET(no).CT.AV(1) = SET(no).CT.AV(1)+xofs;
      SET(no).CT.AV(2) = SET(no).CT.AV(2)+yofs;
      SET(no).CT.Apex(1)= SET(no).CT.Apex(1)+xofs;
      SET(no).CT.Apex(2)= SET(no).CT.Apex(2)+yofs;
    end
    if isfield(SET(no).CT,'avx')
      SET(no).CT.avx = SET(no).CT.avx+xofs;
      SET(no).CT.avy = SET(no).CT.avy+yofs;
      SET(no).CT.apexx = SET(no).CT.apexx+xofs;
      SET(no).CT.apexy = SET(no).CT.apexy+yofs;
      SET(no).CT.orthavx = [nan nan];
      SET(no).CT.orthavy = [nan nan];
      SET(no).CT.orthapexx = nan;
      SET(no).CT.orthapexy = nan;
    end
  end

  SET(no).CenterX = SET(no).CenterX+xofs;
  SET(no).CenterY = SET(no).CenterY+yofs;
  % If outside or near edge, recenter
  if (SET(no).CenterX<0.5+2) || (SET(no).CenterY<0.5+2)...
      || (SET(no).CenterX>SET(no).XSize-2) || (SET(no).CenterY>SET(no).YSize-2)
    SET(no).CenterX = SET(no).XSize/2;
    SET(no).CenterY = SET(no).YSize/2;
  end

  SET(no).Mmode.X = SET(no).Mmode.X+xofs;
  SET(no).Mmode.Y = SET(no).Mmode.Y+yofs;

  warn = false;
  types = {'Endo','Epi','RVEndo','RVEpi'};
  for loop = 1:numel(types)
    if ~isempty(SET(no).([types{loop},'X']))
      SET(no).([types{loop},'X']) = SET(no).([types{loop},'X']) + xofs;
      SET(no).([types{loop},'Y']) = SET(no).([types{loop},'Y']) + yofs;

      tblr = (SET(no).([types{loop},'X']) < 0.5) | ...
        (SET(no).([types{loop},'X']) > (SET(no0).XSize + 0.5)) | ...
        (SET(no).([types{loop},'Y']) < 0.5) | ...
        (SET(no).([types{loop},'Y']) > (SET(no0).YSize + 0.5));
      if any(tblr(:))
        warn = true;
        [SET(no).([types{loop},'X']), SET(no).([types{loop},'Y'])] = ...
          cropcontour(SET(no).([types{loop},'X']),SET(no).([types{loop},'Y']));
      end
    end
  end

  types = {'LA','RA'};
  for loop = 1:numel(types)
    if (~isempty(SET(no).(types{loop})) && ~isempty(SET(no).(types{loop}).X))
      SET(no).(types{loop}).X = SET(no).(types{loop}).X + xofs;
      SET(no).(types{loop}).Y = SET(no).(types{loop}).Y + yofs;

      tblr = (SET(no).(types{loop}).X < 0.5) | ...
        (SET(no).(types{loop}).X > (SET(no0).XSize + 0.5)) | ...
        (SET(no).(types{loop}).Y < 0.5) | ...
        (SET(no).(types{loop}).Y > (SET(no0).YSize + 0.5));
      if any(tblr(:))
        warn = true;
        [SET(no).(types{loop}).X, SET(no).(types{loop}).Y] = ...
          cropcontour(SET(no).(types{loop}).X, SET(no).(types{loop}).Y);
      end
    end
  end

  % Interpolation points are lost, no big deal.
  segmentation('removeallinterp_Callback',true,no);

  if not(isempty(SET(no).Flow))
    if not(isempty(SET(no).Flow.PhaseCorr))
      switch ndims(SET(no).Flow.PhaseCorr)
        case 2
          SET(no).Flow.PhaseCorr = SET(no).Flow.PhaseCorr(xind,yind);
        case 3
          SET(no).Flow.PhaseCorr = SET(no).Flow.PhaseCorr(xind,yind,:);
        case 4
          SET(no).Flow.PhaseCorr = SET(no).Flow.PhaseCorr(xind,yind,:,:);
      end
    end
  end

  SET(no).NormalZoomState = [];
  SET(no).MontageRowZoomState = [];
  SET(no).MontageZoomState = [];

  %crop .LevelSet
  if ~isempty(SET(no).LevelSet) && isfield(SET(no).LevelSet,'Object')

    O = SET(no).LevelSet.Object;

    %crop stored objects
    for obj = 1:O.getnumobjects()

      %Special uncompress
      if O.iscompressed(obj)
        O.uncompress(obj,oldsize);
      end

      temp = O.getbw(obj);

      if ~isempty(temp)
        temp = temp(xind,yind,:);
        O.setbw(obj,temp);
        O.setcut(obj,[]);
        O.setdeform(obj,[]);
        O.compress(obj);
      end

    end
  
    %---things to change

    %Zoomstate
    SET(no).LevelSet.View.RZoomState = [];
    SET(no).LevelSet.View.GZoomState = [];
    SET(no).LevelSet.View.BZoomState = [];

    %slice
    SET(no).LevelSet.View.GSlice = max(1,min(round(SET(no).LevelSet.View.GSlice+yofs),SET(no).YSize));
    SET(no).LevelSet.View.BSlice = max(1,min(round(SET(no).LevelSet.View.BSlice+xofs),SET(no).XSize));

    segment3dp.graphics('update3D');
  end

  %remove Strain tagging analysis  
  if not(isempty(SET(no).StrainTagging))
    logdisp('Cropping image stacks clears the Strain tagging analysis',true);
    if ~isempty(DATA.GUI.StrainTagging)
      straintagging.straintagging('close_Callback');
    end
    SET(no).StrainTagging = [];
  end
  %remove Strain MITT analysis
  logmsg = 'Cropping image stacks clears the Strain MITT analysis';
  strainmitt.strainmitt('removestrainanalysis',no,logmsg);


end %loop over image stacks

%--- Some things that are only done for Mag in flow case.
no=nos(1);

if SET(no).RoiN>0
  roidels=[];
  for rloop=1:SET(no).RoiN
    SET(no).Roi(rloop).X = SET(no).Roi(rloop).X+xofs;
    SET(no).Roi(rloop).Y = SET(no).Roi(rloop).Y+yofs;
    tblr=   (SET(no).Roi(rloop).X<0.5) ...
      | (SET(no).Roi(rloop).X>(SET(no0).XSize+0.5)) ...
      | (SET(no).Roi(rloop).Y<0.5) ...
      | (SET(no).Roi(rloop).Y>(SET(no0).YSize+0.5));
    if any(tblr(:))
      [SET(no).Roi(rloop).X,SET(no).Roi(rloop).Y]=...
        cropcontour(SET(no).Roi(rloop).X,SET(no).Roi(rloop).Y,true);
      if all(all(isnan(SET(no).Roi(rloop).X)))
        roidels=[roidels rloop]; %#ok<AGROW>
      else
        warn=true;
      end
    end
  end
  for loop=1:length(roidels)
    roi('roidelete_Callback',roidels(loop));
  end
end

%Point
if ~isempty(SET(no).Point)
  SET(no).Point.X=SET(no).Point.X+xofs;
  SET(no).Point.Y=SET(no).Point.Y+yofs;

  ind=((SET(no).Point.X>0.5)&(SET(no).Point.X<(SET(no).XSize+0.5))...
    &(SET(no).Point.Y>0.5)&(SET(no).Point.Y<(SET(no).YSize+0.5)));
  SET(no).Point.X = SET(no).Point.X(ind);
  SET(no).Point.Y = SET(no).Point.Y(ind);
  SET(no).Point.T = SET(no).Point.T(ind);
  SET(no).Point.Z = SET(no).Point.Z(ind);
  SET(no).Point.Label = SET(no).Point.Label(ind);
end

%Measure
if ~isempty(SET(no).Measure)
  ind=false(length(SET(no).Measure));
  for mloop=1:length(SET(no).Measure)
    SET(no).Measure(mloop).X=SET(no).Measure(mloop).X+xofs;
    SET(no).Measure(mloop).Y=SET(no).Measure(mloop).Y+yofs;
    ind(mloop)=all(((SET(no).Measure(mloop).X>0.5)&(SET(no).Measure(mloop).X<(SET(no).XSize+.5))...
      &(SET(no).Measure(mloop).Y>0.5)&(SET(no).Measure(mloop).Y<(SET(no).YSize+.5))));
  end
  SET(no).Measure = SET(no).Measure(ind);

end
%---

%Max & min value
if isfield(SET(no),'maxValue')
  SET(no).minValue = single(min(SET(no).IM(:)));
  SET(no).maxValue = single(max(SET(no).IM(:)));

  %Check if IM is a single then need to convert to true data.
  if isa(SET(no).IM,'single')
    SET(no).minValue = calcfunctions('calctruedata',SET(no).minValue,no);
    SET(no).maxValue = calcfunctions('calctruedata',SET(no).maxValue,no);
  end

end

if (warn)
  mywarning('Cropping contours.',DATA.GUI.Segment);
  segment('updatevolume');
end

%----------------------------------------------------------
function [outx,outy] = cropcontour(inx,iny,isroi)
%----------------------------------------------------------
%Used to crop and reinterpolate contours.

global DATA SET NO

if nargin==2
  isroi=false;
end

if isroi
  tloops=1:SET(NO).TSize;
  zloops=1;
else
  tloops=1:SET(NO).TSize;
  zloops=1:SET(NO).ZSize;
end

% using NO is ok, since flow is always same X/Y/T/Z size.

% Consider rewriting to store things in SET/DATA, so that not so much data
% is passed.

%prebuild
outx=nan(DATA.NumPoints,tloops(end),zloops(end));
outy=nan(DATA.NumPoints,tloops(end),zloops(end));
for tloop=tloops
  for zloop=zloops
    top= inx(:,tloop,zloop)<0.5;
    left=iny(:,tloop,zloop)<0.5;
    bottom=inx(:,tloop,zloop)>(SET(NO).XSize+0.5);
    right=iny(:,tloop,zloop)>(SET(NO).YSize+0.5);

    ind = ~(top|bottom|left|right);
    if any(ind)&&~all(isnan(inx(:,tloop,zloop)))
      tmpx=inx(ind,tloop,zloop)';
      tmpy=iny(ind,tloop,zloop)';
      % in case 1/end closure has been cropped
      if (tmpx(1)~=tmpx(end))
        tmpx=[tmpx tmpx(1)]; %#ok<AGROW>
        tmpy=[tmpy tmpy(1)]; %#ok<AGROW>
      end
      % calculate length
      len = sqrt(...
        conv2(tmpx,[1 -1],'valid').^2+...
        conv2(tmpy,[1 -1],'valid').^2);
      len = [0;len(:)]; %Add zero first
      len = cumsum(len);
      tempind = find(conv2(len,[1;-1],'valid')~=0); %Remove doublets
      len = [len(1);len(tempind+1)];
      tmpx = [tmpx(1) tmpx(tempind+1)];
      tmpy = [tmpy(1) tmpy(tempind+1)];
      totallength = len(end);
      %Resample
      xr = interp1(len,tmpx,linspace(0,totallength,DATA.NumPoints),'linear');
      yr = interp1(len,tmpy,linspace(0,totallength,DATA.NumPoints),'linear');
      % orientation
      %Code did not seem to do anything commented out EH:
      %mx = mean(xr);
      %my = mean(yr);
      %       if sum(unwrap(conv2(angle(complex(xr-mx,yr-my)),[1 -1],'valid')))<0
      %         tmpx = fliplr(xr);
      %         tmpy = fliplr(yr);
      %         %disp('counterclockwise');
      %       else
      %         %disp('clockwise');
      %       end;
      % commit
      outx(:,tloop,zloop)=xr(:);
      outy(:,tloop,zloop)=yr(:);
    end
  end
end

%-------------------------------
function removethis_Callback
%-------------------------------
%Remove current slice.
global SET NO

%Removes current slice
ind = true(SET(NO).ZSize,1);
ind(SET(NO).CurrentSlice) = false;
removeslices_Callback(ind);
disableundo;

%-------------------------------
function x = resamplehelper(f,x)
%-------------------------------
%Helper function to resample image stacks.

%factor 2 => x' = 2*x-0.5
%factor 3 => x' = 3*x-1
%factor 4 => x' = 4*x-1.5
%factor 5 => x' = 5*x-2
%factor 2.5 => x' = 2.5*x-1
%factor 3.5 => x' = 3.5*x-1.5
%factor 0.5 => x' = x'*0.5 (a odd)
%factor 0.5 => x' = x'*0.5+0.5 (a even)

if f>0
  d = (ceil(f)-1)/2;
else
  d = 0;
end

if isa(x,'double') || isa(x,'single')
  x = x*f-d;
else
  for xloop=1:size(x,1)
    for yloop=1:size(x,2)
      x{xloop,yloop} = x{xloop,yloop}*f-d;
    end
  end
end

%------------------------------
function makeisotropic_Callback
%------------------------------
%Make an image stack isotropic

makeisotropic;
%call the root file, not very elegant, this function should be moved.

%----------------------------------------------
function newvol = upsampleinplane(f,vol,silent)
%----------------------------------------------
%Helper function to upsample a volume vol.

global DATA

if nargin < 3
  silent = DATA.Silent;
end

if DATA.Silent
  silent = 1;
end

%Find new size
newsize = size(imresize(vol(:,:,1,1),f,'bicubic'));

%Reserve memory
newvol = repmat(0*vol(1),[...
  newsize(1) ...
  newsize(2) ...
  size(vol,3) ...
  size(vol,4)]);

%Loop over image volume
if size(vol,3)>1 && ~silent
  h = mywaitbarstart(size(vol,3),'Upsampling/downsampling volume',1,DATA.GUI.Segment);
end

for tloop=1:size(vol,3)
  for zloop=1:size(vol,4)
    newvol(:,:,tloop,zloop) = single(imresize(vol(:,:,tloop,zloop),f,'bicubic'));
  end
  if size(vol,3)>1 && ~silent
    h = mywaitbarupdate(h);
  end
end

if size(vol,3)>1 && ~silent
  mywaitbarclose(h);
end

%--------------------------------------
function newvol = upsamplevolume2(f,vol,silent)
%--------------------------------------
%Helper function to upsample a volume vol. in this code f is a vector  with
%input [inplane x scalefactor, inplane y scale factor, temporal scale factor, slice scale factor]
global DATA

%This function is not recommended. Instead use upsampleinplane

if nargin < 3
  silent = DATA.Silent;
end

%Find new size
oldsize = size(vol);
f = f(1:length(oldsize)); %safety check for example if there is no slice dimension
newsize = round(size(vol).*f);

oldtypefunc = str2func(class(vol));

if  ~silent
  h = mywaitbarstart(3,'Please wait.',1,DATA.GUI.Segment);
  h = mywaitbarupdate(h);
end

vectors = cell(1,length(oldsize));

for i=1:length(vectors)
  vectors{i} = linspace(1,oldsize(i),newsize(i));
end
grid = cell(1,sum(oldsize~=1));
[grid{:}] = ndgrid(vectors{oldsize~=1});

if  ~silent
  h = mywaitbarupdate(h);
end

%GPU does not work nicely
if false %DATA.Pref.GPU.Use
  try
    newvol = interpn(gpuArray(squeeze(vol)),grid{:},'linear');

    newvol = gather(newvol);
  catch
    %newvol = interpn(single(squeeze(vol)),grid{:},'nearest');
    %newvol = interpn(single(squeeze(vol)),grid{:},'linear');
    resetgpu;

    %resort to CPU
    logdisp('Resorting to CPU.',true);
    newvol = interpn(single(squeeze(vol)),grid{:},'linear');

    %newvol = interpn(gpuArray(single(squeeze(vol))),grid{:},'linear');
    %newvol = gather(newvol);
  end
else
  newvol = interpn(single(squeeze(vol)),grid{:},'linear');
end
if  ~silent
  h = mywaitbarupdate(h);
end

newvol = oldtypefunc(reshape(newvol,newsize));

if  ~silent
  mywaitbarclose(h);
end
myworkoff;

%---------------------------------------------
function newvol = upsampleslices(f,vol,silent)
%---------------------------------------------
%Upsamples along last dimension.

%This is a major rewrite compared to the old one that was slow.

%Einar Heiberg

global DATA

if nargin < 3
  silent = DATA.Silent;
end

%Find new size
switch ndims(vol)
  case 4
    temp = squeeze(vol(1,1,1,:));
  case 3
    temp = squeeze(vol(1,1,:));
  case 2
    temp = squeeze(vol(1,:));
  case 1
    temp = vol;
end

xi = linspace(1,length(temp),round(length(temp)*f));
newlastdimsize = length(interp1(single(temp),xi,'nearest'));

if newlastdimsize<1
  myfailed('Too few resulting slices. Aborting.',DATA.GUI.Segment);
  return;
end

%Reserve memory
type = class(vol);
zeroel = feval(type,0); %Casts zero to the class

%Compute the new size
newsize = size(vol);
newsize(ndims(vol)) = newlastdimsize;

method = 'bilinear'; %affects case 4,3
switch ndims(vol)
  case 4
    goterror = false;
    if usegpu()
      try
        %Loop over image volume
        if ~silent
          h = mywaitbarstart(size(vol,3)*size(vol,2),'Upsampling/downsampling image stack',1,DATA.GUI.Segment);
        end
        %---Use GPU acceleration
        logdisp('GPU accelerated upsampleslice 4dim',true);
        gvol = single(gpuArray(vol)); %send to gpu and convert to single
        gnewvol = gpuArray(repmat(zeroel,newsize)); %reserve memory on gpu
        for tloop = 1:size(vol,3)
          for yloop = 1:size(vol,2)
            gim = imresize(squeeze(gvol(:,yloop,tloop,:)),[newsize(1) newsize(4)],'cubic');
            gnewvol(:,yloop,tloop,:) = gim;
            if (~silent) && isequal(rem(yloop,20),0)
              h = mywaitbarset(yloop/size(vol,2),h);
            end
          end
        end
        newvol = feval(type,gather(gnewvol)); %get and convert to type
      catch e %in case got error mark it
        goterror = true;
        logdisp('--- Error in GPU upsample ---',true);
        mydispexception(e);
        logdisp('--- Error in GPU upsample ---',true);
      end
      if ~silent
        mywaitbarclose(h);
      end
    end

    if ~usegpu() || goterror
      logdisp('non GPU accelerated upsampleslice 4dim',true);
      if ~silent
        h = mywaitbarstart(size(vol,3)*size(vol,2),'Upsampling/downsampling image stack',1,DATA.GUI.Segment);
      end
      newvol = repmat(zeroel,newsize);
      for tloop = 1:size(vol,3)
        for yloop = 1:size(vol,2)
          im = imresize(single(squeeze(vol(:,yloop,tloop,:))),[newsize(1) newsize(4)],method);
          newvol(:,yloop,tloop,:) = feval(type,im);
          if (~silent) && isequal(rem(yloop,20),0)
            h = mywaitbarset(yloop/size(vol,2),h);
          end
        end
      end
      if ~silent
        mywaitbarclose(h);
      end
    end

  case 3 %three dimensional
    goterror = false;
    if usegpu()
      try
        %Loop over image volume
        if ~silent
          h = mywaitbarstart(size(vol,2),'Upsampling/downsampling image stack',1,DATA.GUI.Segment);
        end
        %---Use GPU acceleration
        logdisp('GPU accelerated upsampleslice 3dim',true);
        gvol = single(gpuArray(vol)); %send to gpu and convert to single
        gnewvol = gpuArray(repmat(zeroel,newsize)); %reserve memory
        for yloop = 1:size(vol,2)
          gim = imresize(squeeze(gvol(:,yloop,:)),[newsize(1) newsize(3)],'cubic');
          gnewvol(:,yloop,:) = gim;
          if (~silent) && isequal(rem(yloop,20),0)
            h = mywaitbarset(yloop/size(vol,2),h);
          end
        end
        newvol = feval(type,gather(gnewvol)); %get and convert to type
      catch e %in case got error mark it
        goterror = true;
        logdisp('--- Error in GPU upsample ---',true);
        mydispexception(e);
        logdisp('--- Error in GPU upsample ---',true);
      end
      if ~silent
        mywaitbarclose(h);
      end
    end

    if ~usegpu() || goterror
      logdisp('non GPU accelerated upsampleslice 3dim',true);
      if ~silent
        h = mywaitbarstart(size(vol,2),'Upsampling/downsampling image stack',1,DATA.GUI.Segment);
      end
      newvol = repmat(zeroel,newsize); %reserve memory
      for yloop = 1:size(vol,2)
        im = imresize(single(squeeze(vol(:,yloop,:))),[newsize(1) newsize(3)],method);
        newvol(:,yloop,:) = feval(type,im);
        if (~silent) && isequal(rem(yloop,20),0)
          h = mywaitbarset(yloop/size(vol,2),h);
        end
      end
      if ~silent
        mywaitbarclose(h);
      end
    end
  case 2
    %a*z
    switch class(vol)
      case {'single','double'}
        newvol = interp1(vol',xi,'linear')';
      otherwise
        newvol = interp1(single(vol'),xi,'linear')';
    end
  case 1
    %z
    newvol = interp1(vol(:),xi,'linear');
end

%-----------------------------
function precomp_Callback
%------------------------------
%Function to make intensity precompensation of MR
%gradient echo images. This function is somewhat obsoleted.

global SET NO

tempnos = NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

if (isequal(SET(NO).ImagingTechnique,'MRGE')&&(SET(NO).TSize>1))
  %gradientecho
  %Precompensate for inlow artifacts
  h = mywaitbarstart(SET(NO).ZSize,'Please wait, compensating inflow artifacts.');
  m = zeros(1,SET(NO).TSize);
  s = zeros(1,SET(NO).TSize);

  for zloop=1:SET(NO).ZSize
    for tloop=1:SET(NO).TSize
      temp = SET(NO).IM(:,:,tloop,zloop);
      m(tloop)=mean(temp(:));
      s(tloop)=std(temp(:));
    end
    mm = mean(m);
    ss = mean(s);
    for tloop=1:SET(NO).TSize
      if not(ss==0)
        SET(NO).IM(:,:,tloop,zloop)=(SET(NO).IM(:,:,tloop,zloop)-m(tloop))*(ss/s(tloop))+mm;
      end
    end
    h = mywaitbarupdate(h);
  end
  mywaitbarclose(h);
end

%-----------------------------
function mirrorx_Callback
%-----------------------------
%Mirros in x direction.
%Need to flip both x and z to maintain a righthand system.
flipx_helper;
segment('updatevolume');

%full refresh needed
viewfunctions('setview')

%----------------------------------
function invertcolors_Callback
%----------------------------------
%Invert colors in current image stack (essentially 1-x).
%Need to fix more to get intensity offset correct.

global DATA SET NO

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

SET(NO).IM = 1-SET(NO).IM;
DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimage',DATA.CurrentPanel);

%Find overall max
tempmin = min(SET(NO).IM(:));
tempmax = max(SET(NO).IM(:));

%Normalize %only if bigger than one.
if tempmax>1
  if not(tempmax==0)
    SET(NO).IM = SET(NO).IM-tempmin;
    SET(NO).IM = SET(NO).IM./(tempmax-tempmin);
    SET(NO).IntensityScaling = tempmax-tempmin;
    SET(NO).IntensityOffset = tempmin;
  else
    SET(NO).IntensityScaling = 0;
    SET(NO).IntensityOffset = tempmin;
  end
else
  SET(NO).IntensityScaling = 1;
end

%-------------------------------
function normalize_Callback(no)
%-------------------------------
%Normalize image data of current image stack. Make sure image
%intensities are in the range [0..1]. Also stores in so that
%true image itensities can be retrieved, see calctruedata

global DATA SET NO

if nargin == 0
  no = NO;
end

if not(isempty(SET(no).Flow))
  if ~isempty(SET(no).Flow.PhaseNo)&&(no==SET(no).Flow.PhaseNo)
    myfailed('You can not normalize phase image',DATA.GUI.Segment);
    return;
  end
  if ~isempty(SET(no).Flow.PhaseX)&&(no==SET(no).Flow.PhaseX)
    myfailed('You can not normalize phase image',DATA.GUI.Segment);
    return;
  end
  if ~isempty(SET(no).Flow.PhaseY)&&(no==SET(no).Flow.PhaseY)
    myfailed('You can not normalize phase image',DATA.GUI.Segment);
    return;
  end
elseif isfield(SET(no),'Linked') && ~isempty(SET(no).Linked)
  %Normalize all linked stacks
  logdisp('Normalizing data',true);
  for linked_no = SET(no).Linked
    if isa(SET(linked_no).IM,'single')
      [tempmin,tempmax] = mynormalize(SET(linked_no).IM);

      SET(linked_no).IntensityScaling = tempmax-tempmin;
      SET(linked_no).IntensityOffset = tempmin;
    else
      SET(linked_no).IntensityScaling = 1;
      SET(linked_no).IntensityOffset = 0;
    end
  end

  DATA.EndoEdgeDetected = false;
  DATA.EpiEdgeDetected = false;
  return
end

%Normalize %only if bigger than one.

logdisp('Normalizing data',true);

% This fcn makes changes to SET(NO).IM directly
% Very uggly but saves a lot of memory for large datasets.
if isa(SET(no).IM,'single')
  [tempmin,tempmax] = mynormalize(SET(no).IM);

  SET(no).IntensityScaling = tempmax-tempmin;
  SET(no).IntensityOffset = tempmin;
else
  SET(no).IntensityScaling = 1;
  SET(no).IntensityOffset = 0;
end
%The above call is equavivalent to:
%tempmin = min(SET(NO).IM(:));
%tempmax = max(SET(NO).IM(:));
%SET(NO).IM = SET(NO).IM-tempmin;
%SET(NO).IM = SET(NO).IM./(tempmax-tempmin);

DATA.EndoEdgeDetected = false;
DATA.EpiEdgeDetected = false;

%------------------------------------------
function undosegmentation_Callback(no)
%------------------------------------------
%Revert segmentation from undo history

global DATA SET NO

if nargin==0
  no = NO;
  if ~isempty(SET(NO).Parent)
    no = SET(NO).Parent;
  end
end

if isempty(DATA.Undo)
  mywarning('Nothing to undo.',DATA.GUI.Segment);
  disableundo;
  return;
end

if (DATA.UndoN(no)<1)
  mywarning('Nothing to undo.',DATA.GUI.Segment);
  disableundo;
  return;
end

SET(no).EndoX =    DATA.Undo(no,DATA.UndoN(no)).EndoXBackup;
SET(no).EndoY =    DATA.Undo(no,DATA.UndoN(no)).EndoYBackup;
SET(no).EpiX =     DATA.Undo(no,DATA.UndoN(no)).EpiXBackup;
SET(no).EpiY =     DATA.Undo(no,DATA.UndoN(no)).EpiYBackup;
SET(no).RVEndoX =  DATA.Undo(no,DATA.UndoN(no)).RVEndoXBackup;
SET(no).RVEndoY =  DATA.Undo(no,DATA.UndoN(no)).RVEndoYBackup;
SET(no).RVEpiX =   DATA.Undo(no,DATA.UndoN(no)).RVEpiXBackup;
SET(no).RVEpiY =   DATA.Undo(no,DATA.UndoN(no)).RVEpiYBackup;
SET(no).EndoInterpX = DATA.Undo(no,DATA.UndoN(no)).EndoInterpXBackup;
SET(no).EndoInterpY = DATA.Undo(no,DATA.UndoN(no)).EndoInterpYBackup;
SET(no).EpiInterpX =  DATA.Undo(no,DATA.UndoN(no)).EpiInterpXBackup;
SET(no).EpiInterpY =  DATA.Undo(no,DATA.UndoN(no)).EpiInterpYBackup;
SET(no).RVEndoInterpX = DATA.Undo(no,DATA.UndoN(no)).RVEndoInterpXBackup;
SET(no).RVEndoInterpY = DATA.Undo(no,DATA.UndoN(no)).RVEndoInterpYBackup;
SET(no).RVEpiInterpX =  DATA.Undo(no,DATA.UndoN(no)).RVEpiInterpXBackup;
SET(no).RVEpiInterpY =  DATA.Undo(no,DATA.UndoN(no)).RVEpiInterpYBackup;
SET(no).RoiN =    DATA.Undo(no,DATA.UndoN(no)).RoiN;
SET(no).RoiCurrent = DATA.Undo(no,DATA.UndoN(no)).RoiCurrent;
SET(no).Roi = DATA.Undo(no,DATA.UndoN(no)).Roi;
SET(no).EndoDraged = DATA.Undo(no,DATA.UndoN(no)).EndoDragedBackup;
SET(no).EpiDraged = DATA.Undo(no,DATA.UndoN(no)).EpiDragedBackup;
SET(no).CenterX = DATA.Undo(no,DATA.UndoN(no)).CenterX;
SET(no).CenterY = DATA.Undo(no,DATA.UndoN(no)).CenterY;
SET(no).Point = DATA.Undo(no,DATA.UndoN(no)).Point;
SET(no).Measure = DATA.Undo(no,DATA.UndoN(no)).Measure;
undopos = DATA.UndoN(no);
SET(no).LA = DATA.Undo(no,undopos).LA;
SET(no).RA = DATA.Undo(no,undopos).RA;

%OBS: For flow data these are still SET-specific /JU
% Need a better solution for undoing contrast and zoom.
% SET(NO).IntensityMapping = DATA.Undo(no,DATA.UndoN(NO)).IntensityMapping;
% SET(NO).IntensityOffset = DATA.Undo(no,DATA.UndoN(NO)).IntensityOffset;
% SET(NO).IntensityScaling = DATA.Undo(no,DATA.UndoN(NO)).IntensityScaling;
% SET(NO).NormalZoomState = DATA.Undo(NO,DATA.UndoN(NO)).NormaZoomState;
% SET(NO).MontageZoomState = DATA.Undo(NO,DATA.UndoN(NO)).MontageZoomState;
% SET(NO).MontageRowZoomState = DATA.Undo(NO,DATA.UndoN(NO)).MontageRowZoomState;
% % Must be regenerated for contrast updates
% DATA.ViewIM{DATA.CurrentPanel} = [];
% makeviewim;
% if DATA.UseLight
%   DATA.BalloonLevel = -1; %Force update of ballonimage
% end;

if ~isempty(SET(NO).Scar)
  SET(NO).Scar.Manual(:)=int8(0);
  SET(NO).Scar.Manual(DATA.Undo(no,DATA.UndoN(no)).ScarManualP)=int8(1);
  SET(NO).Scar.Manual(DATA.Undo(no,DATA.UndoN(no)).ScarManualN)=int8(-1);
  SET(NO).Scar.Result(:)=false;
  SET(NO).Scar.Result(DATA.Undo(no,DATA.UndoN(no)).ScarNoReflow)=true;
  SET(NO).Scar.NoReflow(:)=false;
  SET(NO).Scar.NoReflow(DATA.Undo(no,DATA.UndoN(no)).ScarNoReflow)=true;
  viability('viabilitycalc');
end

if ~isempty(SET(NO).MaR)
  SET(NO).MaR.Manual(:)=int8(0);
  SET(NO).MaR.Manual(DATA.Undo(no,DATA.UndoN(no)).MaRManualP)=int8(1);
  SET(NO).MaR.Manual(DATA.Undo(no,DATA.UndoN(no)).MaRManualN)=int8(-1);
  SET(NO).MaR.Result(:)=false;
  mar('update');
end

%Clear some data to store memory, do not need to clear all, just big ones
DATA.Undo(no,DATA.UndoN(no)).EndoXBackup = [];
DATA.Undo(no,DATA.UndoN(no)).EndoYBackup = [];
DATA.Undo(no,DATA.UndoN(no)).EpiXBackup = [];
DATA.Undo(no,DATA.UndoN(no)).EpiYBackup = [];
DATA.Undo(no,DATA.UndoN(no)).LA = [];
DATA.Undo(no,DATA.UndoN(no)).RA = [];
% DATA.Undo(no,DATA.UndoN(no)).RoiX = [];
% DATA.Undo(no,DATA.UndoN(no)).RoiY = [];how should this be changed after
% Roi turned to struct and only saved as .Roi in DATA.Undo

%Decrease pointer
DATA.UndoN(no) = DATA.UndoN(no)-1;
if DATA.UndoN(no)<1
  DATA.UndoN(no) = 0;
  disableundo;
end

segment('updatevolume');
segment('updateflow');
segment('updatemeasurement');
drawfunctions('drawno',no,true)
% for p = find(DATA.ViewPanels)
%   DATA.ViewIM{p} = [];
%   drawfunctions('drawno',p,true)
% end

%-----------------------------------------------
function success = smoothsegmentation_helperold(no)
%-----------------------------------------------
%This function is no longer used and should be removed when closing ticket
%#2907. /Einar

global DATA SET

success = true;
%smooth segmentation in current slice and current time frame
switch DATA.CurrentTool
  case {'Endo','EndoInterp','EndoBalloon'}
    %Double check so that there is a cell with the correct size for the used
    %type
    %delete interpolation points if they exist
    tools('connectinterpolation',no,{'EndoInterp'})
    if ~isempty(SET(no).EndoInterpX)
      if ~all([size(SET(no).EndoInterpX,1)==SET(no).TSize,size(SET(no).EndoInterpX,2)==SET(no).ZSize])
        tmpX = cell(SET(no).TSize,SET(no).ZSize);
        tmpY = cell(SET(no).TSize,SET(no).ZSize);
        sz = size(SET(no).EndoInterpX);
        tmpX(1:sz(1),1:sz(2)) = SET(no).EndoInterpX(1:sz(1),1:sz(2));
        tmpY(1:sz(1),1:sz(2)) = SET(no).EndoInterpX(1:sz(1),1:sz(2));
        SET(no).EndoInterpX = tmpX;
        SET(no).EndoInterpY = tmpY;
      end
      SET(no).EndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
      SET(no).EndoInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
    end

    %   callbackfunctions('interpdeletepointall','EndoInterp');
    %smooth LV endo
    if not(isnan(SET(no).EndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      endox = SET(no).EndoX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      endoy = SET(no).EndoY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      sigma = 1;
      %re-parameterize to radial components
      centerx = mean(endox);
      centery = mean(endoy);
      dist2center = sqrt((endox-centerx).^2+(endoy-centery).^2);
      tempsmoothedline = lvrvtools('smooth1D',[dist2center ; dist2center ; dist2center],sigma);
      tempsmoothedline = lvrvtools('smooth1D',tempsmoothedline,sigma);
      smoothedline = tempsmoothedline(length(dist2center)+1:2*length(dist2center));
      %re-paramterize to x-y-components again
      tempangle = asin((endoy-centery)./dist2center);
      angle = zeros(length(endox),1);
      for angleloop = 1:length(endoy)
        if (endoy(angleloop)-centery) > 0 && (endox(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop);
        elseif (endoy(angleloop)-centery) < 0 && (endox(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop)+2*pi;
        elseif (endoy(angleloop)-centery) < 0 && (endox(angleloop)-centerx) > 0
          angle(angleloop) = pi-tempangle(angleloop);
        else
          angle(angleloop) = pi-tempangle(angleloop);
        end
      end
      SET(no).EndoX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = centerx+smoothedline.*sin(angle-pi/2);
      SET(no).EndoY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = centery+smoothedline.*cos(angle-pi/2);
      SET(no).EndoX(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).EndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      SET(no).EndoY(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).EndoY(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    else
      mywarning('No LV endo segmentation to smooth.',DATA.GUI.Segment);
      success = false;
      return;
    end

  case {'Epi','EpiInterp','EpiBalloon'}
    %Double check so that there is a cell with the correct size for the used
    %type
    tools('connectinterpolation',no,{'EpiInterp'})
    %delete interpolation points if they exist
    if ~isempty(SET(no).EpiInterpX)
      if ~all([size(SET(no).EpiInterpX,1)==SET(no).TSize,size(SET(no).EpiInterpX,2)==SET(no).ZSize])
        tmpX = cell(SET(no).TSize,SET(no).ZSize);
        tmpY = cell(SET(no).TSize,SET(no).ZSize);
        sz = size(SET(no).EpiInterpX);
        tmpX(1:sz(1),1:sz(2)) = SET(no).EpiInterpX(1:sz(1),1:sz(2));
        tmpY(1:sz(1),1:sz(2)) = SET(no).EpiInterpX(1:sz(1),1:sz(2));
        SET(no).EpiInterpX = tmpX;
        SET(no).EpiInterpY = tmpY;
      end

      SET(no).EpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
      SET(no).EpiInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
    end

    %   callbackfunctions('interpdeletepointall','EpiInterp');
    %smooth LV epi
    if not(isnan(SET(no).EpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      epix = SET(no).EpiX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      epiy = SET(no).EpiY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      sigma = 1;
      %re-parameterize to radial components
      centerx = mean(epix);
      centery = mean(epiy);
      dist2center = sqrt((epix-centerx).^2+(epiy-centery).^2);
      tempsmoothedline = lvrvtools('smooth1D',[dist2center ; dist2center ; dist2center],sigma);
      tempsmoothedline = lvrvtools('smooth1D',tempsmoothedline,sigma);
      smoothedline = tempsmoothedline(length(dist2center)+1:2*length(dist2center));
      %re-paramterize to x-y-components again
      tempangle = asin((epiy-centery)./dist2center);
      angle = zeros(length(epix),1);
      for angleloop = 1:length(epiy)
        if (epiy(angleloop)-centery) > 0 && (epix(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop);
        elseif (epiy(angleloop)-centery) < 0 && (epix(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop)+2*pi;
        elseif (epiy(angleloop)-centery) < 0 && (epix(angleloop)-centerx) > 0
          angle(angleloop) = pi-tempangle(angleloop);
        else
          angle(angleloop) = pi-tempangle(angleloop);
        end
      end
      SET(no).EpiX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = centerx+smoothedline.*sin(angle-pi/2);
      SET(no).EpiY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = centery+smoothedline.*cos(angle-pi/2);
      SET(no).EpiX(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).EpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      SET(no).EpiY(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).EpiY(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    else
      success = false;
      mywarning('No LV epi segmentation to smooth.',DATA.GUI.Segment);
      return;
    end
  case {'RVEndo','RVEndoInterp','RVEndoBalloon'}
    %Double check so that there is a cell with the correct size for the used
    %type
    %delete interpolation points if they exist
    tools('connectinterpolation',no,{'RVEndoInterp'})
    if ~isempty(SET(no).RVEndoInterpX)
      if ~all([size(SET(no).RVEndoInterpX,1)==SET(no).TSize,size(SET(no).RVEndoInterpX,2)==SET(no).ZSize])
        tmpX = cell(SET(no).TSize,SET(no).ZSize);
        tmpY = cell(SET(no).TSize,SET(no).ZSize);
        sz = size(SET(no).RVEndoInterpX);
        tmpX(1:sz(1),1:sz(2)) = SET(no).RVEndoInterpX(1:sz(1),1:sz(2));
        tmpY(1:sz(1),1:sz(2)) = SET(no).RVEndoInterpX(1:sz(1),1:sz(2));
        SET(no).RVEndoInterpX = tmpX;
        SET(no).RVEndoInterpY = tmpY;
      end
      SET(no).RVEndoInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
      SET(no).RVEndoInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
    end
    %smooth RV endo
    if not(isnan(SET(no).RVEndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      endox = SET(no).RVEndoX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      endoy = SET(no).RVEndoY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      sigma = 1;
      %re-parameterize to radial components
      centerx = mean(endox);
      centery = mean(endoy);
      dist2center = sqrt((endox-centerx).^2+(endoy-centery).^2);
      tempsmoothedline = lvrvtools('smooth1D',[dist2center ; dist2center ; dist2center],sigma);
      tempsmoothedline = lvrvtools('smooth1D',tempsmoothedline,sigma);
      smoothedline = tempsmoothedline(length(dist2center)+1:2*length(dist2center));
      %re-paramterize to x-y-components again
      tempangle = asin((endoy-centery)./dist2center);
      angle = zeros(length(endox),1);
      for angleloop = 1:length(endoy)
        if (endoy(angleloop)-centery) > 0 && (endox(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop);
        elseif (endoy(angleloop)-centery) < 0 && (endox(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop)+2*pi;
        elseif (endoy(angleloop)-centery) < 0 && (endox(angleloop)-centerx) > 0
          angle(angleloop) = pi-tempangle(angleloop);
        else
          angle(angleloop) = pi-tempangle(angleloop);
        end
      end
      x = centerx+smoothedline.*sin(angle-pi/2);
      y = centery+smoothedline.*cos(angle-pi/2);
      [x,y] = calcfunctions('calcpointsoutsideLV',x,y,no);
      SET(no).RVEndoX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = x;
      SET(no).RVEndoY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = y;
      SET(no).RVEndoX(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).RVEndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      SET(no).RVEndoY(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).RVEndoY(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    else
      success = false;
      mywarning('No RV endo segmentation to smooth.',DATA.GUI.Segment);
      return;
    end

  case {'RVEpi','RVEpiInterp'}
    %Double check so that there is a cell with the correct size for the used
    %type
    %delete interpolation points if they exist
    tools('connectinterpolation',no,{'RVEpiInterp'})
    if ~isempty(SET(no).RVEpiInterpX)
      if ~all([size(SET(no).RVEpiInterpX,1)==SET(no).TSize,size(SET(no).RVEpiInterpX,2)==SET(no).ZSize])
        tmpX = cell(SET(no).TSize,SET(no).ZSize);
        tmpY = cell(SET(no).TSize,SET(no).ZSize);
        sz = size(SET(no).RVEpiInterpX);
        tmpX(1:sz(1),1:sz(2)) = SET(no).RVEpiInterpX(1:sz(1),1:sz(2));
        tmpY(1:sz(1),1:sz(2)) = SET(no).RVEpiInterpX(1:sz(1),1:sz(2));
        SET(no).RVEpiInterpX = tmpX;
        SET(no).RVEpiInterpY = tmpY;
      end
      SET(no).RVEpiInterpX{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
      SET(no).RVEpiInterpY{SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
    end
    %smooth RV epi
    if not(isnan(SET(no).RVEpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      epix = SET(no).RVEpiX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      epiy = SET(no).RVEpiY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      sigma = 1;
      %re-parameterize to radial components
      centerx = mean(epix);
      centery = mean(epiy);
      dist2center = sqrt((epix-centerx).^2+(epiy-centery).^2);
      tempsmoothedline = lvrvtools('smooth1D',[dist2center ; dist2center ; dist2center],sigma);
      tempsmoothedline = lvrvtools('smooth1D',tempsmoothedline,sigma);
      smoothedline = tempsmoothedline(length(dist2center)+1:2*length(dist2center));
      %re-paramterize to x-y-components again
      tempangle = asin((epiy-centery)./dist2center);
      angle = zeros(length(epix),1);
      for angleloop = 1:length(epiy)
        if (epiy(angleloop)-centery) > 0 && (epix(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop);
        elseif (epiy(angleloop)-centery) < 0 && (epix(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop)+2*pi;
        elseif (epiy(angleloop)-centery) < 0 && (epix(angleloop)-centerx) > 0
          angle(angleloop) = pi-tempangle(angleloop);
        else
          angle(angleloop) = pi-tempangle(angleloop);
        end
      end
      SET(no).RVEpiX(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = centerx+smoothedline.*sin(angle-pi/2);
      SET(no).RVEpiY(1:end-1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = centery+smoothedline.*cos(angle-pi/2);
      SET(no).RVEpiX(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).RVEpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      SET(no).RVEpiY(end,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = SET(no).RVEpiY(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    else
      success = false;
      mywarning('No RV epi segmentation to smooth.',DATA.GUI.Segment);
      return;
    end
  case {'Roi','drawroi','putroi','RoiBalloon'}
    %smooth ROI
    if not(isempty(SET(no).RoiCurrent))
      x = SET(no).Roi(SET(no).RoiCurrent).X(1:end-1,SET(no).CurrentTimeFrame);
      y = SET(no).Roi(SET(no).RoiCurrent).Y(1:end-1,SET(no).CurrentTimeFrame);
      sigma = 1;
      %re-parameterize to radial components
      centerx = mean(x);
      centery = mean(y);
      dist2center = sqrt((x-centerx).^2+(y-centery).^2);
      tempsmoothedline = lvrvtools('smooth1D',[dist2center ; dist2center ; dist2center],sigma);
      tempsmoothedline = lvrvtools('smooth1D',tempsmoothedline,sigma);
      smoothedline = tempsmoothedline(length(dist2center)+1:2*length(dist2center));
      %re-paramterize to x-y-components again
      tempangle = asin((y-centery)./dist2center);
      angle = zeros(length(x),1);
      for angleloop = 1:length(y)
        if (y(angleloop)-centery) > 0 && (x(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop);
        elseif (y(angleloop)-centery) < 0 && (x(angleloop)-centerx) < 0
          angle(angleloop) = tempangle(angleloop)+2*pi;
        elseif (y(angleloop)-centery) < 0 && (x(angleloop)-centerx) > 0
          angle(angleloop) = pi-tempangle(angleloop);
        else
          angle(angleloop) = pi-tempangle(angleloop);
        end
      end
      roiXtemp = centerx+smoothedline.*sin(angle-pi/2);
      roiYtemp = centery+smoothedline.*cos(angle-pi/2);
      roiXtemp(end+1) = roiXtemp(1);
      roiYtemp(end+1) = roiYtemp(1);
      SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame) = roiXtemp;
      SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame) = roiYtemp;
    else
      success = false;
      mywarning('No ROI segmentation to smooth.',DATA.GUI.Segment);
      return;
    end
end

%ensure strain calculation is re-done if lv segmentation is changed
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

if ~isempty(SET(no).StrainMitt)
  updateredo(SET(no).StrainMitt,DATA.CurrentTool);
end
calcfunctions('updatemarandscar',no);

%---------------------------
function smoothlvall(no, sigma)
%---------------------------
% Function to smooth LV contours in all slices and timeframes
global SET

if nargin < 4
  sigma = 1;
end
zsize = SET(no).ZSize;
tsize = SET(no).TSize;

%removes interpolation points for the slice and timeframe
tools('connectinterpolation',no,{'EndoInterp'})
removeinterp(no, 'EndoInterp');
tools('connectinterpolation',no,{'EpiInterp'})
removeinterp(no, 'EpiInterp');

existendo = squeeze(any(~isnan(SET(no).EndoX), 1));
existepi = squeeze(any(~isnan(SET(no).EpiX), 1));

for t = 1:tsize %loop through time
  for z = 1:zsize %loop through slices
    if existendo(t, z)
      smoothcontour('Endo', no, t, z, sigma);
    end
    if existepi(t, z)
      smoothcontour('Epi', no, t, z, sigma);
    end
  end
end

segment('updatevolume');
drawfunctions('drawno', no);

%---------------------------
function smoothcontour(type, no, tf, slice, sigma)
%---------------------------
global SET

if strcmp(type, 'Endo')
  x = SET(no).EndoX(1:end-1, tf, slice);
  y = SET(no).EndoY(1:end-1, tf, slice);
elseif strcmp(type, 'Epi')
  x = SET(no).EpiX(1:end-1, tf, slice);
  y = SET(no).EpiY(1:end-1, tf, slice);
else
  return;
end
if any(isnan(x))
  return;
end

% Center and radial distance
cx = mean(x);
cy = mean(y);
r = sqrt((x - cx).^2 + (y - cy).^2);

% Smooth radial distance
r_smooth = lvrvtools('smooth1D', [r; r; r], sigma);
r_smooth = lvrvtools('smooth1D', r_smooth, sigma);
r_smooth = r_smooth(length(r)+1:2*length(r));

% calculate angle
theta = atan2(y - cy, x - cx);

% Reconstruct smoothed contour
x_new = cx + r_smooth .* cos(theta);
y_new = cy + r_smooth .* sin(theta);

% Store back
x_new(end+1) = x_new(1);
y_new(end+1) = y_new(1);

if strcmp(type, 'Endo')
  SET(no).EndoX(:, tf, slice) = x_new;
  SET(no).EndoY(:, tf, slice) = y_new;
elseif strcmp(type, 'Epi')
  SET(no).EpiX(:, tf, slice) = x_new;
  SET(no).EpiY(:, tf, slice) = y_new;
end


%-----------------------------------------------
function success = smoothsegmentation_helper(no)
%-----------------------------------------------
%Smooth the segmentation, uses SET struct and checks DATA.CurrentTool for
%which contour to smooth. 
%
%Major rewrite of old function by Einar. Old function ensured that the
%segmentation starts in the same orientation as (i.e at 3 o clock). New
%code does not do that. The reason that for the rewrite was that the old
%code had some bug making the contour not smooth.

global DATA SET

success = true;

%smooth segmentation in current slice and current time frame
switch DATA.CurrentTool
  case {'Endo','EndoInterp','EndoBalloon'}
    %Work on EndoX,EndoY
    if ~isempty(SET(no).EndoX) && not(isnan(SET(no).EndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      X = SET(no).EndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      Y = SET(no).EndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      fieldstring = 'Endo';
      interpstring = 'EndoInterp';
    else
      mywarning('No LV endo segmentation to smooth.',DATA.GUI.Segment);
      success = false;
      return;
    end

  case {'Epi','EpiInterp','EpiBalloon'}
    %Double check so that there is a cell with the correct size for the used
    %type
    tools('connectinterpolation',no,{'EpiInterp'})
    if not(isnan(SET(no).EpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      X = SET(no).EpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      Y = SET(no).EpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      interpstring = 'EpiInterp';
      fieldstring = 'Epi';
    else
      success = false;
      mywarning('No LV epi segmentation to smooth.',DATA.GUI.Segment);
      return;
    end
  case {'RVEndo','RVEndoInterp','RVEndoBalloon'}
    tools('connectinterpolation',no,{'RVEndoInterp'});    
    %smooth RV endo
    if not(isnan(SET(no).RVEndoX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      X = SET(no).RVEndoX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      Y = SET(no).RVEndoY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      interpstring = 'RVEndoInterp';
      fieldstring = 'RVEndo';
    else
      success = false;
      mywarning('No RV endo segmentation to smooth.',DATA.GUI.Segment);
      return;
    end

  case {'RVEpi','RVEpiInterp'}
    if not(isnan(SET(no).RVEpiX(1,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)))
      X = SET(no).RVEpiX(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      Y = SET(no).RVEpiY(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      interpstring = 'RVEpiInterp';
      fieldstring = 'RVEpi';
    else
      success = false;
      mywarning('No RV epi segmentation to smooth.',DATA.GUI.Segment);
      return;
    end
  case {'Roi','drawroi','putroi','RoiBalloon'}
    %smooth ROI
    interpstring = '';
    fieldstring = 'Roi';
    if not(isempty(SET(no).RoiCurrent))
      X = SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame);
      Y = SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame);
    else
      success = false;
      mywarning('No ROI segmentation to smooth.',DATA.GUI.Segment);
      return;
    end
  otherwise
    return
end

%Take cae of interpolation points, this code just ensures correct size and
%removes interpolation points for the slice and timeframe
tools('connectinterpolation',no,{interpstring})
if ~isempty(interpstring)
  removeinterp(no, interpstring);
end

%--- Interpolation of the field

x = X(1:end-1);
y = Y(1:end-1);
sigma = 1;

np = length(x);
x = lvrvtools('smooth1D',[x ; x; x],sigma);
y = lvrvtools('smooth1D',[y ; y; y],sigma);
x = x((np+1):(2*np)); %cut back
y = y((np+1):(2*np));

%Store
if ~isequal(fieldstring,'Roi')
  %Normal case
  SET(no).([fieldstring 'X'])(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = [x;x(1)];
  SET(no).([fieldstring 'Y'])(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice) = [y;y(1)];
else
  SET(no).Roi(SET(no).RoiCurrent).X(:,SET(no).CurrentTimeFrame) = [x;x(1)];
  SET(no).Roi(SET(no).RoiCurrent).Y(:,SET(no).CurrentTimeFrame) = [y;y(1)];
end

%ensure strain calculation is re-done if lv segmentation is changed
if ~isempty(SET(no).StrainTagging) && isfield(SET(no).StrainTagging, 'LVupdated')
  SET(no).StrainTagging.LVupdated = 1;
end

if ~isempty(SET(no).StrainMitt)
  updateredo(SET(no).StrainMitt,DATA.CurrentTool);
end
calcfunctions('updatemarandscar',no);

%------------------------------------------
function removeinterp(no, interpstring)
%------------------------------------------
global SET
 if ~isempty(SET(no).([interpstring 'X']))
    if ~all([size(SET(no).([interpstring 'X']),1)==SET(no).TSize,size(SET(no).([interpstring 'X']),2)==SET(no).ZSize])
      tmpX = cell(SET(no).TSize,SET(no).ZSize);
      tmpY = cell(SET(no).TSize,SET(no).ZSize);
      sz = size(SET(no).([interpstring 'X']));
      tmpX(1:sz(1),1:sz(2)) = SET(no).([interpstring 'X'])(1:sz(1),1:sz(2));
      tmpY(1:sz(1),1:sz(2)) = SET(no).([interpstring 'Y'])(1:sz(1),1:sz(2));
      SET(no).([interpstring 'X']) = tmpX;
      SET(no).([interpstring 'Y']) = tmpY;
    end

    %Store
    SET(no).([interpstring 'X']){SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
    SET(no).([interpstring 'Y']){SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
  end
%------------------------------------------
function smoothsegmentation_Callback(no)
%------------------------------------------
%Smooth latest segmentation (LV, RV and ROI)
global SET NO

if nargin==0
  no = NO;
  if ~isempty(SET(NO).Parent)
    no = SET(NO).Parent;
  end
end

%smooth segmentation in current slice and current time frame
success = smoothsegmentation_helper(no);
if ~success
  return
end
segment('updatevolume');

drawfunctions('drawno',no);


%--------------------------------
function endsystole_Callback(noupdate,single,saxes)
%--------------------------------
%Change current time frame to end systole
global DATA SET NO

if (DATA.Run) && (findindented(DATA.Handles.playiconholder,'play'))
  % stop playing
  DATA.Run = 0;
  undent(DATA.Handles.playiconholder,'play',0)
end

if nargin==0
  noupdate=false;
end

if nargin<2
  single=1;
end

if nargin<3
  saxes=false; %define if using SAX ES tf, otherwise ES tf for current stack
end

nos = DATA.ViewPanels(DATA.ViewPanels~=0);%SET(NO).Linked;

%also check for imageviewplane
imvpcell={'Short-axis','Transversal','2CH','3CH','4CH','Flow'};
imtcell={'Flow (magnitude)','Flow (through plane)'};
include=zeros(1,length(nos));
for i = 1 : length(nos)
  include(i) = ismember(SET(nos(i)).ImageViewPlane,imvpcell)||ismember(SET(nos(i)).ImageType,imtcell);
end

%if no image stack fulfilling the criteria above, still change for current
%image stack
if sum(include(:)) == 0
  [~,NOindex] = find(nos==NO);
  include(NOindex) = 1;
end

nos=nos(logical(include));
if length(DATA.LVNO)==1
  no = DATA.LVNO;
else
  no = findfunctions('findcineshortaxisno');
end
useother=1;

if isempty(no) || SET(no).EDT==SET(no).EST
  %not correct use NO instead
  no=NO;
  useother=0;
end

%Adjustments have been made to NO use this instead
if ~SET(NO).EDT==SET(NO).EST
  no=NO;
  useother=0;
end
%nos = unique([nos no]); %ensure that the current image stack is changed

if useother
  for loop = nos
    if SET(loop).EDT==SET(loop).EST
      [~,SET(loop).EST]=min(abs(SET(no).EST/SET(no).TSize-(1:SET(loop).TSize)/SET(loop).TSize));
    end
  end
end

if single
  SET(NO).CurrentTimeFrame = SET(NO).EST;
  drawfunctions('drawpanel',DATA.CurrentPanel);
else
  for loop=1:length(nos)
    if saxes %go to tf that corresponds to SAX ES tf
      %find tf by using the same code as calcfunctions('calcsegmentationintersections')
      tfsaxes = SET(no).EST;
      if tfsaxes ~= 1
        nophase = SET(nos(loop)).TimeVector(tfsaxes)/SET(nos(loop)).TimeVector(end);
        TV_tmp = nan(length(no),max([SET(no).TSize]));
        for i = 1:length(no)
          TV_tmp(i,1:SET(no(i)).TSize) = SET(no(i)).TimeVector;
        end
        denomntr = [SET(no).TSize]-1;
        if any(denomntr == 0)
          % if denominator can be 0, then set time frames to 1
          tfs = ones(1,length(no));
        else
          %Then we find the matching phases in the other nos
          [~,tfs] = min((TV_tmp./([SET(no).TIncr].*(denomntr))'-nophase)'.^2);
        end
      else
        tfs = 1;
      end
      %[~,tfs] = min(abs(SET(no).TimeVector(SET(no).EST)-SET(nos(loop)).TimeVector));
      SET(nos(loop)).CurrentTimeFrame = tfs;
    else %go to ES tf
      SET(nos(loop)).CurrentTimeFrame = SET(nos(loop)).EST;
    end
  end
  for p = find(ismember(DATA.ViewPanels,nos))
    drawfunctions('drawpanel',p);
  end
end
if ~noupdate
  segment('updatevolume');
end
viewfunctions('updatetimebars');

%---------------------------------
function enddiastole_Callback(noupdate,single,saxed)
%---------------------------------
%change current time frame to end diastole if there exists SAX with EDT and
%
global DATA SET NO

if (DATA.Run) && (findindented(DATA.Handles.playiconholder,'play'))
  % stop playing
  DATA.Run = 0;
  undent(DATA.Handles.playiconholder,'play',0)
end


if nargin==0
  noupdate=false;
end

if nargin<2
  single=1;
end

if nargin<3
  saxed=false;
end

nos = DATA.ViewPanels(DATA.ViewPanels~=0);%SET(NO).Linked;

%also check for imageviewplane
imvpcell={'Short-axis','Transversal','2CH','3CH','4CH','Flow'};
imtcell={'Flow (magnitude)','Flow (through plane)'};
include=zeros(1,length(nos));
for i = 1 : length(nos)
  include(i) = ismember(SET(nos(i)).ImageViewPlane,imvpcell)||ismember(SET(nos(i)).ImageType,imtcell);
end

%if no image stack fulfilling the criteria above, still change for current
%image stack
if sum(include(:)) == 0
  [~,NOindex] = find(nos==NO);
  include(NOindex) = 1;
end

nos=nos(logical(include));
if length(DATA.LVNO)==1%~isempty(DATA.LVNO)
  no = DATA.LVNO;
else
  no = findfunctions('findcineshortaxisno');
end
useother=1;

if isempty(no) || SET(no).EDT==SET(no).EST
  %not correct use NO instead
  no=NO;
  useother=0;
end
%nos = unique([nos no]); %ensure that the current image stack is changed

%Adjustments have been made to NO use this instead
if ~SET(NO).EDT==SET(NO).EST
  no=NO;
  useother=0;
end

if useother
  for loop = nos
    if SET(loop).EDT==SET(loop).EST
      [~,SET(loop).EDT]=min(abs(SET(no).EDT/SET(no).TSize-(1:SET(loop).TSize)/SET(loop).TSize));
    end
  end
end


if single
  SET(NO).CurrentTimeFrame = SET(NO).EDT;
  drawfunctions('drawpanel',DATA.CurrentPanel);
else
  for loop=1:length(nos)

    if saxed %go to tf that corresponds to SAX ED tf
      %find tf by using the same code as calcfunctions('calcsegmentationintersections')
      tfsaxed = SET(no).EDT;
      if tfsaxed ~= 1
        nophase = SET(nos(loop)).TimeVector(tfsaxed)/SET(nos(loop)).TimeVector(end);
        TV_tmp = nan(length(no),max([SET(no).TSize]));
        for i = 1:length(no)
          TV_tmp(i,1:SET(no(i)).TSize) = SET(no(i)).TimeVector;
        end
        denomntr = [SET(no).TSize]-1;
        if any(denomntr == 0)
          % if denominator can be 0, then set time frames to 1
          tfs = ones(1,length(no));
        else
          %Then we find the matching phases in the other nos
          [~,tfs] = min((TV_tmp./([SET(no).TIncr].*(denomntr))'-nophase)'.^2);
        end
      else
        tfs = 1;
      end
      SET(nos(loop)).CurrentTimeFrame = tfs;

    else
      SET(nos(loop)).CurrentTimeFrame = SET(nos(loop)).EDT;
    end
  end
  for p = find(ismember(DATA.ViewPanels,nos))
    drawfunctions('drawpanel',p);
  end
end
if ~noupdate
  segment('updatevolume');
end
viewfunctions('updatetimebars');

%--------------------------------
function endsystoleall_Callback
%--------------------------------
%Change current timeframe of all image stacks to end systole.
endsystole_Callback(0,0)
%global DATA SET NO
%
% oldNO=NO;
%
% nos=setdiff(unique(DATA.ViewPanels),0);
% noupdate=true;
% for loop=1:length(nos)
%   NO=nos(loop);
%   if SET(NO).EDT == SET(NO).EST
%     SET(NO).CurrentTimeFrame = round(1+(SET(oldNO).EST-1)/(SET(oldNO).TSize-1)*(SET(NO).TSize-1));
%     drawfunctions('drawsliceno',NO);
%   else
%     endsystole_Callback(noupdate);
%   end
% end
%
% NO=oldNO;
% segment('updatevolume');

%--------------------------------
function enddiastoleall_Callback
%--------------------------------
%Change current timeframe of all image stacks to end diastole.
enddiastole_Callback(0,0)

% global DATA SET NO
%
% oldNO=NO;
%
% nos=setdiff(unique(DATA.ViewPanels),0);
% noupdate=true;
% for loop=1:length(nos)
%   NO=nos(loop);
%   if SET(NO).EDT == SET(NO).EST
%     SET(NO).CurrentTimeFrame = round(1+(SET(oldNO).EDT-1)/(SET(oldNO).TSize-1)*(SET(NO).TSize-1));
%     drawfunctions('drawsliceno',NO);
%   else
%     enddiastole_Callback(noupdate);
%   end
% end
%
% NO=oldNO;
% segment('updatevolume');

%-----------------------------------------
function copyforwardselected_Callback(lv)
%-----------------------------------------
%Copy segmentation of selected slices forward in time.
global DATA SET NO

enableundo;

if SET(NO).CurrentTimeFrame==SET(NO).TSize
  tf = 1;
else
  tf = SET(NO).CurrentTimeFrame+1;
end

if nargin > 0 && isequal(lv,'lv')
  currenttool = 'Endo';
else
  currenttool = DATA.CurrentTool;
end

switch currenttool
  case {'Roi','RoiPut'}
    for n = SET(NO).RoiCurrent
      SET(NO).Roi(n).X(:,tf) = SET(NO).Roi(n).X(:,SET(NO).CurrentTimeFrame);
      SET(NO).Roi(n).Y(:,tf) = SET(NO).Roi(n).Y(:,SET(NO).CurrentTimeFrame);
      SET(NO).Roi(n).T = find(not(isnan(SET(NO).Roi(n).X(1,:))));
    end
  case {'Endo','Epi'}
    ind = SET(NO).StartSlice:SET(NO).EndSlice;
    if isempty(ind)
      myfailed('No slices selected.',DATA.GUI.Segment);
      return;
    end

    if ~isempty(SET(NO).EndoX)
      SET(NO).EndoX(:,tf,ind) = SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,ind);
      SET(NO).EndoY(:,tf,ind) = SET(NO).EndoY(:,SET(NO).CurrentTimeFrame,ind);
    end

    if ~isempty(SET(NO).EpiX)
      SET(NO).EpiX(:,tf,ind) = SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,ind);
      SET(NO).EpiY(:,tf,ind) = SET(NO).EpiY(:,SET(NO).CurrentTimeFrame,ind);
    end
end

SET(NO).CurrentTimeFrame = SET(NO).CurrentTimeFrame+1;
drawfunctions('drawpanel',DATA.CurrentPanel);


%--------------------------------------------------
function copyupward_Callback(type,silent,dolv,lvalg)
%--------------------------------------------------
%Copies upward and refine.
%type is either endo or epi.
%silent is true if to avoid graphic update.
%if dolv false then copy rv.

global DATA SET NO

if nargin<2
  silent = false;
end

if nargin<3
  if isequal(DATA.CurrentTheme,'function')
    dolv = true;
  else
    dolv = false;
  end
end

if nargin<4
  lvalg = 'new';
end

if nargin<1
  endo = 1;
  epi = 1;
  rvendo = 1;
  rvepi = 1;
else
  endo = isequal(type,'endo');
  epi = isequal(type,'epi');
  rvendo = isequal(type,'rvendo');
  rvepi = isequal(type,'rvepi');
  if isequal(type,'lv')
    endo =1;
    epi = 1;
    dolv = true;
  end  
end

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = false;
%   if not(existfunctions('existrvendoinselected',NO))
%     DATA.ThisFrameOnly = true;
%   end
% end

if DATA.ThisFrameOnly
  timeframes = SET(NO).CurrentTimeFrame;
else
  timeframes = 1:SET(NO).TSize;
end

%check if contour exist in selected all slices and timeframes
docontinue=true;
if dolv
  if endo %lvendo
    tools('connectinterpolation',NO,{'EndoInterp'})
    if not(existfunctions('existendoinselected',NO))%check so that lvendo exist in all slices and frames
      docontinue=false;
      errmsg='Endo must exist in all timeframes or ';
    end
  else %lvepi
    tools('connectinterpolation',NO,{'EpiInterp'})
    if not(existfunctions('existepiinselected',NO))%check so that lvepi exist in all slices and frames
      docontinue=false;
      errmsg='Epi must exist in all timeframes or ';
    end
  end
else%rvendo
  tools('connectinterpolation',NO,{'RVEndoInterp'})
  if not(existfunctions('existrvendoinselected',NO))%check so that rvendo exist in all slices and frames
    docontinue=false;
    errmsg='RV Endo must exist in all timeframes or ';
  end
end

%do not continue calculations if contour does not exist in all selectedtimeframes
%ans slices
if not(docontinue)
  tip = DATA.singleframetip;
  myfailed(sprintf('%s%s',errmsg,tip));
  return;
end


if SET(NO).CurrentSlice>1

  if dolv
    if (endo)&&(~isempty(SET(NO).EndoX))
      SET(NO).EndoX(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).EndoX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).EndoY(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).EndoY(:,timeframes,SET(NO).CurrentSlice);
    end
  else
    if (rvendo)&&(~isempty(SET(NO).RVEndoX))
      SET(NO).RVEndoX(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).RVEndoX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).RVEndoY(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).RVEndoY(:,timeframes,SET(NO).CurrentSlice);
    end
  end

  if dolv
    if (epi)&&(~isempty(SET(NO).EpiY))
      SET(NO).EpiX(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).EpiX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).EpiY(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).EpiY(:,timeframes,SET(NO).CurrentSlice);
    end
  else
    if (rvepi)&&(~isempty(SET(NO).RVEpiY))
      SET(NO).RVEpiX(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).RVEpiX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).RVEpiY(:,timeframes,SET(NO).CurrentSlice-1) = SET(NO).RVEpiY(:,timeframes,SET(NO).CurrentSlice);
    end
  end

  SET(NO).CurrentSlice = SET(NO).CurrentSlice-1;
  SET(NO).StartSlice = SET(NO).CurrentSlice;
  SET(NO).EndSlice = SET(NO).StartSlice;

  if dolv
    if endo
      switch lvalg
        case 'old'
          lv('segmentrefineendo_Callback',silent);
        case 'new'
          lvrvtools('segmentrefineendo_Callback',silent);
      end
    end
  elseif rvendo
    rv('segmentrefinervendo_Callback',silent,false);
  end
  if epi && dolv
    switch lvalg
      case 'old'
        lv('segmentrefineepi_Callback',silent);
      case 'new'
        lvrvtools('segmentrefineepi_Callback',silent);
    end
  end

  if not(silent)
    DATA.ViewIM{DATA.CurrentPanel} = [];
    drawfunctions('drawno',NO)
  end
end

%----------------------------------------------------
function copydownward_Callback(type,silent,dolv,lvalg)
%----------------------------------------------------
%Copy segmentation in current slice downwards and refine
global DATA SET NO

if nargin<2
  silent = false;
end

if nargin<3
  % Be extra sure
  if isequal(DATA.CurrentTheme,'function')
    dolv = true;
  else
    dolv = false;
  end
end

if nargin<4
  lvalg = 'new';
end

if nargin<1
  endo = 1;
  epi = 1;
  rvendo = 1;
  rvepi = 1;
else
  endo = isequal(type,'endo');
  epi = isequal(type,'epi');
  rvendo = isequal(type,'rvendo');
  rvepi = isequal(type,'rvepi');
  if isequal(type,'lv')
    endo =1;
    epi = 1;
    dolv = true;
  end
end

%if strcmp(DATA.ProgramName,'Segment CMR')
% for Segment CMR, do it for all time frames if endo exist in all,
% otherwise only in current time frame
%   DATA.ThisFrameOnly = false;
%   if not(existfunctions('existrvendoinselected',NO))
%     DATA.ThisFrameOnly = true;
%   end
%end
if DATA.ThisFrameOnly
  timeframes = SET(NO).CurrentTimeFrame;
else
  timeframes = 1:SET(NO).TSize;
end

%check if contour exist in selected all slices and timeframes
docontinue=true;
if dolv
  tools('connectinterpolation',NO,{'EndoInterp'})
  if endo %lvendo
    if not(existfunctions('existendoinselected',NO))%check so that lvendo exist in all slices and frames
      docontinue=false;
      errmsg='Endo must exist in all timeframes or ';
    end
  end
else%rvendo
  tools('connectinterpolation',NO,{'RVEndoInterp'})
  if not(existfunctions('existrvendoinselected',NO))%check so that rvendo exist in all slices and frames
    docontinue=false;
    errmsg='RV Endo must exist in all timeframes or ';
  end
end
if dolv && epi %lvepi
  tools('connectinterpolation',NO,{'EpiInterp'})
  if not(existfunctions('existepiinselected',NO))%check so that lvepi exist in all slices and frames
    docontinue=false;
    errmsg='Epi must exist in all timeframes or ';
  end
end
%do not continue calculations if contour does not exist in all selectedtimeframes
%ans slices
if not(docontinue)
  tip = DATA.singleframetip;
  myfailed(sprintf('%s%s',errmsg,tip));
  return;
end


if SET(NO).CurrentSlice<SET(NO).ZSize

  if dolv
    if (endo)&&(~isempty(SET(NO).EndoX))
      SET(NO).EndoX(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).EndoX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).EndoY(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).EndoY(:,timeframes,SET(NO).CurrentSlice);
    end
  else
    if (rvendo)&&(~isempty(SET(NO).RVEndoX))
      SET(NO).RVEndoX(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).RVEndoX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).RVEndoY(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).RVEndoY(:,timeframes,SET(NO).CurrentSlice);
    end
  end
  if dolv
    if (epi)&&(~isempty(SET(NO).EpiX))
      SET(NO).EpiX(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).EpiX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).EpiY(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).EpiY(:,timeframes,SET(NO).CurrentSlice);
    end
  else
    if (rvepi)&&(~isempty(SET(NO).RVEpiX))
      SET(NO).RVEpiX(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).RVEpiX(:,timeframes,SET(NO).CurrentSlice);
      SET(NO).RVEpiY(:,timeframes,SET(NO).CurrentSlice+1) = SET(NO).RVEpiY(:,timeframes,SET(NO).CurrentSlice);
    end
  end

  SET(NO).CurrentSlice = SET(NO).CurrentSlice+1;
  SET(NO).StartSlice = SET(NO).CurrentSlice;
  SET(NO).EndSlice = SET(NO).StartSlice;
  if dolv
    if endo
      switch lvalg
        case 'old'
          lv('segmentrefineendo_Callback',silent);
        case 'new'
          lvrvtools('segmentrefineendo_Callback',silent);
      end
    end
  else
    rv('segmentrefinervendo_Callback',silent,false);
  end
  if dolv && epi
    switch lvalg
      case 'old'
        lv('segmentrefineepi_Callback',silent);
      case 'new'
        lvrvtools('segmentrefineepi_Callback',silent);
    end
  end

  if not(silent)
    DATA.ViewIM{DATA.CurrentPanel} = [];
    drawfunctions('drawno',NO)
  end
end

%--------------------------------------------------
function copytoalltimeframes_Callback(type,silent,dolv)
%--------------------------------------------------
%Copies segmentation in the selected timeframe and all slices to all timeframes.
%type is either endo, epi, or lv.
%silent is true if to avoid graphic update.
%if dolv false then copy RV.

global DATA SET NO


if nargin<2
  silent = false;
end

if nargin<3
  if isequal(DATA.CurrentTheme,'lv')
    dolv = true;
  elseif isequal(DATA.CurrentTheme,'rv')
    dolv = false;
  else
    %select manually if no or several 2ch found.
    s=cell(1,2);
    s(1,1)=cellstr('LV');
    s(1,2)=cellstr('RV');
    m=mymenu('Copy LV or RV segmentation?',s);
    if m==1
      dolv = true;
    else
      dolv = false;
    end
  end
end

if nargin<1
  endo = 1;
  epi = 1;
else
  endo = isequal(type,'endo');
  epi = isequal(type,'epi');
  if isequal(type,'lv')
    endo = 1;
    epi = 1;
  end
end

%check if contour exist in selected timeframe
docontinue=true;
haveendo = true;
haveepi = true;
if dolv
  if endo
    anyseg=existfunctions('existendoinselected',NO,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice);
    if anyseg==0%check so that lvendo exist in any slice and current time frame
      haveendo = false;
      errmsg='LV must exist in this timeframe.';
    end
  end
  if epi
    anyseg=existfunctions('existepiinselected',NO,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice);
    if anyseg==0 %check so that lvepi exist in any slice and current time frame
      haveepi = false;
      errmsg='LV must exist in this timeframe.';
    end
    if not(haveendo) && not(haveepi)
      docontinue=false;
    end
  end
else %RV endo
  anyseg=existfunctions('existrvendoinselected',NO,SET(NO).CurrentTimeFrame,SET(NO).CurrentSlice);
  if anyseg==0%check so that rvendo exist in all slices and frames
    docontinue=false;
    errmsg='RV Endo must exist in all timeframe.';
  end
end

%do not continue calculations if contour does not exist in the selected timeframe
if not(docontinue)
  myfailed(sprintf('%s',errmsg));
  return;
end

for t=1:SET(NO).TSize
  if dolv
    if (endo)&&(~isempty(SET(NO).EndoX))
      SET(NO).EndoX(:,t,:) = SET(NO).EndoX(:,SET(NO).CurrentTimeFrame,:);
      SET(NO).EndoY(:,t,:) = SET(NO).EndoY(:,SET(NO).CurrentTimeFrame,:);
    end
  else
    if (endo)&&(~isempty(SET(NO).RVEndoX))
      SET(NO).RVEndoX(:,t,:) = SET(NO).RVEndoX(:,SET(NO).CurrentTimeFrame,:);
      SET(NO).RVEndoY(:,t,:) = SET(NO).RVEndoY(:,SET(NO).CurrentTimeFrame,:);
    end
  end

  if dolv
    if (epi)&&(~isempty(SET(NO).EpiY))
      SET(NO).EpiX(:,t,:) = SET(NO).EpiX(:,SET(NO).CurrentTimeFrame,:);
      SET(NO).EpiY(:,t,:) = SET(NO).EpiY(:,SET(NO).CurrentTimeFrame,:);
    end
  else
    if (epi)&&(~isempty(SET(NO).RVEpiY))
      SET(NO).RVEpiX(:,t,:) = SET(NO).RVEpiX(:,SET(NO).CurrentTimeFrame,:);
      SET(NO).RVEpiY(:,t,:) = SET(NO).RVEpiY(:,SET(NO).CurrentTimeFrame,:);
    end
  end
end


if not(silent)
  drawfunctions('drawpanel',DATA.CurrentPanel);
end


%----------------------------------------
function opensetimagedescription_Callback(frompreviewmenu)
%----------------------------------------
%Open the gui for setting the image type, image view plane and
%imaging technique for current image stack
global DATA SET NO

if nargin>0 && frompreviewmenu %from context menu
  no = segment_main('thumbnailno');
else
  no = NO;
end

%Initalize
DATA.GUI.SetImageDescription = mygui('setimagedescription.fig');
gui = DATA.GUI.SetImageDescription;
gui.NO = no;
h = gui.handles;
% update figure name with stack number
stackstr = dprintf('Stack #%d',no);

h.figure1.Name = sprintf('%s (%s)',h.figure1.Name,stackstr);
[type,viewplane,technique] = segment('imagedescription');
set(h.listboximagetype,'String',type);
set(h.listboximageviewplane,'String',viewplane);
set(h.listboximagingtechnique,'String',technique);

%set the parameters
imagetype = SET(no).ImageType;
allimagetypes = get(h.listboximagetype,'String');
imagetypenbr = [];
for loop = 1:length(allimagetypes)
  if strcmp(imagetype,allimagetypes{loop})
    imagetypenbr = loop;
  end
end
if isempty(imagetypenbr)
  type{length(type)} = imagetype;
  type{length(type)+1} = 'User defined';
  set(h.listboximagetype,'String',type,'value',length(type)-1);
else
  set(h.listboximagetype,'value',imagetypenbr);
end

imageviewplane = SET(no).ImageViewPlane;
allimageviewplanes = get(h.listboximageviewplane,'String');
imageviewplanenbr = [];
for loop = 1:length(allimageviewplanes)
  if strcmp(imageviewplane,allimageviewplanes{loop})
    imageviewplanenbr = loop;
  end
end
if isempty(imageviewplanenbr)
  viewplane{length(viewplane)} = imageviewplane;
  viewplane{length(viewplane)+1} = 'User defined';
  set(h.listboximageviewplane,'String',viewplane,'value',length(viewplane)-1);
else
  set(h.listboximageviewplane,'value',imageviewplanenbr);
end

imagingtechnique = SET(no).ImagingTechnique;
allimagingtechniques = get(h.listboximagingtechnique,'String');
imagingtechniquenbr = [];
for loop = 1:length(allimagingtechniques)
  if strcmp(imagingtechnique,allimagingtechniques{loop})
    imagingtechniquenbr = loop;
  end
end
if isempty(imagingtechniquenbr)
  technique{length(technique)} = imagingtechnique;
  technique{length(technique)+1} = 'User defined';
  set(h.listboximagingtechnique,'String',technique,'value',length(technique)-1);
else
  set(h.listboximagingtechnique,'value',imagingtechniquenbr);
end
DATA.setsectorrotation(no);

%------------------------------------
function setimagedescription_Callback
%------------------------------------
%Set the image type, image view plane and
%imaging technique for current image stack
global DATA SET

gui = DATA.GUI.SetImageDescription;
no = gui.NO;
h = gui.handles;

imagetypenbr = get(h.listboximagetype,'value');
imageviewplanenbr = get(h.listboximageviewplane,'value');
imagingtechniquenbr = get(h.listboximagingtechnique,'value');
s = [];
if imagetypenbr == length(get(h.listboximagetype,'String'))
  s.Image_Type = '';
end
if imageviewplanenbr == length(get(h.listboximageviewplane,'String'))
  s.Image_View_Plane = '';
end
if imagingtechniquenbr == length(get(h.listboximagingtechnique,'String'))
  s.Imaging_Technique = '';
end

if ~isempty(s)
  [res,ok] = inputstruct(s,'User defined');
  if ~ok
    myfailed('Invalid user defined image description',DATA.GUI.Segment);
    return;
  end
  if isfield(res,'Image_Type')
    SET(no).ImageType = res.Image_Type;
  else
    temp = get(h.listboximagetype,'String');
    SET(no).ImageType = temp{imagetypenbr};
  end
  if isfield(res,'Image_View_Plane')
    SET(no).ImageViewPlane = res.Image_View_Plane;
  else
    temp = get(h.listboximageviewplane,'String');
    SET(no).ImageViewPlane = temp{imageviewplanenbr};
  end
  if isfield(res,'Imaging_Technique')
    SET(no).ImagingTechnique = res.Imaging_Technique;
  else
    temp = get(h.listboximagingtechnique,'String');
    SET(no).ImagingTechnique = temp{imagingtechniquenbr};
  end
else
  temp = get(h.listboximagetype,'String');
  SET(no).ImageType = temp{imagetypenbr};
  temp = get(h.listboximageviewplane,'String');
  SET(no).ImageViewPlane = temp{imageviewplanenbr};
  temp = get(h.listboximagingtechnique,'String');
  SET(no).ImagingTechnique = temp{imagingtechniquenbr};
end

%set modality based on imaging tecknique
ImTe = SET(no).ImagingTechnique(1:2);
if strcmp(ImTe,'NM')
  SET(no).Modality = 'NM';
elseif strcmp(ImTe,'MR')
  SET(no).Modality = 'MR';
elseif strcmp(ImTe,'CT')
  SET(no).Modality = 'CT';
elseif strcmp(ImTe,'PT')
  SET(no).Modality = 'PT';
elseif strcmp(ImTe,'US')
  SET(no).Modality = 'US';
elseif strcmp(ImTe,'XR')
  SET(no).Modality = 'XR';
elseif strcmp(ImTe,'CR')
  SET(no).Modality = 'CR';
else
  SET(no).Modality = 'OT';
end

%only for tagged imnage stacks
if strcmp(SET(no).ImageType,'Strain from tagging')
  SET(no).Cyclic = false;
end

%close the gui
closesetimagedescription_Callback

%Update stack name
calcpreview = false;
drawfunctions('drawthumbnails',calcpreview);

%update results
segment('updatemeasurement');
panel = find(DATA.ViewPanels == no,1);
if ~isempty(panel) && ~DATA.issegment3dp
  drawfunctions('drawlaparameters',no,SET(no).CurrentTimeFrame,DATA.CurrentPanel);
end

%-----------------------------
function setheartrate_Callback(frompreviewmenu)
%-----------------------------
%Set heart rate of current image stack
global SET NO DATA

if nargin>0 && frompreviewmenu %from context menu
  no = segment_main('thumbnailno');
else
  no = NO;
end

%Init gui
DATA.GUI.HeartRate = mygui('adjustheartrate.fig');
gui = DATA.GUI.HeartRate;
gui.NO = no;
set([gui.handles.textrrinterval, gui.handles.editrrinterval], 'Enable','off');
set(gui.handles.stacknumbertext, 'String',dprintf('Stack #%d',no));

curhr = num2str(SET(no).HeartRate);  % Current Heart Rate
set(gui.handles.editheartrate,'String', curhr);

currr = num2str(SET(no).RRInterval); % Current R-R Interval
set(gui.handles.editrrinterval,'String', currr);

% newhr = str2double(myinputdlg('Heart Rate (bpm)','Enter new value',1,{curhr}));
% if isempty(newhr)
%   return
% end
% if ~isnan(newhr) && newhr >= 0
%   SET(no).HeartRate = newhr;
% 
%   %update Time vector and Time Increment
%   tsize = SET(no).TSize;
%   [tvector, tinc] = updatetimevector_helper(newhr, tsize);
% 
%   if ~isempty(tinc) %new HR was not set to zero
%     SET(no).TimeVector = tvector;
%     SET(no).TIncr = tinc;
%   end
%   
%   %update Time Bar
%   DATA.updatetimebaraxes
% 
%   %Update result panels
%   segment('updatevolume');
%   segment('updateflow');
% 
% else
%   myfailed('Invalid heart rate input');
% end

%-----------------------------------------
function setrrinterval_Callback
%-----------------------------------------
% Helper function to handle R-R interval checkbox in adjustheartrate.fig
 global DATA
 
 gui = DATA.GUI.HeartRate;
 
 if get(gui.handles.checkboxuserr,'Value') == 1
   set([gui.handles.textrrinterval, gui.handles.editrrinterval], 'Enable','on');
   set([gui.handles.textheartrate, gui.handles.editheartrate], 'Enable','off');
 else
   set([gui.handles.textrrinterval, gui.handles.editrrinterval], 'Enable','off');
   set([gui.handles.textheartrate, gui.handles.editheartrate], 'Enable','on');
 end
 
%-----------------------------------------
function sethearrateok_Callback(newhr, no)
%-----------------------------------------
% Helper function to handle OK button in adjustheartrate.fig

global DATA SET

if nargin < 1
  gui = DATA.GUI.HeartRate;
  no = gui.NO;

  if get(gui.handles.checkboxuserr,'Value') == 1 % use R-R interval to calculate Heart Rate
    rr =str2double(get(gui.handles.editrrinterval,'String')); % R-R (ms)
    newhr = gethrfromrrinterval(rr,no); % Calculate HR (bpm) from R-R interval (ms)
  else
    newhr = str2double(get(gui.handles.editheartrate, 'String'));
  end

  %close the gui
  closeadjustheartrate_Callback

  if isempty(newhr)
    return
  end
end

if ~isnan(newhr) && newhr >= 0
  nos = SET(no).Linked;
  for n = nos
    SET(n).HeartRate = newhr;
    SET(n).RRInterval = 60000/newhr; %update RR Intreval
    if ~isempty(SET(n).Flow)
      SET(n).Flow.HeartRate = newhr;
    end
  end

  %update Time vector and Time Increment
  tsize = SET(no).TSize;
  [tvector, tinc] = updatetimevector_helper(newhr, tsize);

  if ~isempty(tinc) %new HR was not set to zero
    for n = nos
      SET(n).TimeVector = tvector;
      SET(n).TIncr = tinc;
    end
  end

  %update Time Bar
  DATA.updatetimebaraxes

  %Update result panels
  segment('updatevolume');
  segment('updateflow');

else
  myfailed('Invalid heart rate input');
end

%-----------------------------------------
function hr = gethrfromrrinterval(newrr,no)
%-----------------------------------------
% Function to calculate Heart Rate (bpm) from R-R Interval (ms)
% Helper function used in adjustheartrate.fig
global SET;

if newrr > 0
  hr = 60000/newrr;
  SET(no).RRInterval = newrr;
elseif newrr == 0
  hr = 0;
  myfailed('R-R interval set to 0 then Heart Rate set to 0');
else
  myfailed('Invalid R-R interval input');
  hr = 0;
end

%-----------------------------------------
function closeadjustheartrate_Callback
%-----------------------------------------
%Close the gui for setting heart rate (adjustheartrate.fig)
global DATA

try
  close(DATA.GUI.HeartRate);
catch
  close(gcbf);
end
DATA.GUI.HeartRate = [];

%-----------------------------------------
function [tvector, tinc] = updatetimevector_helper(newhr, tsize)
%-----------------------------------------
% Function to update time vector after changing Heart Rate (HR)
% Helper funtion to setheartrate_Callback

if newhr > 0
  tinc = 1/(newhr/60)/tsize; %new time increment
  tvector = (0:tsize-1)*tinc; % new time vector
else
  tvector = []; %invalid input of new HR
  tinc = []; %invalid  input of new HR
end
 
%-----------------------------------------
function closesetimagedescription_Callback
%-----------------------------------------
%Close the gui for setting the image description
global DATA

try
  close(DATA.GUI.SetImageDescription);
catch
  close(gcbf);
end
DATA.GUI.SetImageDescription = [];

%--------------------------------
function outstring = maskpatientstrings(instring,patientname,patientid)
%--------------------------------
% function to mask patient identity in strings
global SET
try
  if nargin < 2
    % no orginal name was provided
    if ~isempty(SET) && isfield(SET(1), 'PatientInfo') && isfield(SET(1).PatientInfo, 'Name')
      if ~isempty(SET(1).PatientInfo.Name)
        patientname = SET(1).PatientInfo.Name;
      else
        patientname = 'ThereIsNoPatientNameInThisDataSet';
      end
    end
  elseif isempty(SET) && nargin < 3
    % return when SET is empty and no extra values were provided
    outstring = instring;
    return
  end
  if nargin < 3
    % no orginal ID was provided
    if ~isempty(SET) && isfield(SET(1), 'PatientInfo') && isfield(SET(1).PatientInfo, 'Name')
      if ~isempty(SET(1).PatientInfo.Name)
        patientid = SET(1).PatientInfo.ID;
      else
        patientid = 'ThereIsNoPatientIDInThisDataSet';
      end
    end
  end
  if ~isempty(patientname)
    % replace space with '_'
    patientname = strrep(patientname,' ','_');
    %check for patientname and patient ID in the provided instring
    if contains(instring,patientname,'IgnoreCase',true)
      instring = regexprep(instring,patientname,'patientname','ignorecase');
    end
  end
  if ~isempty(patientid)
    patientid = removeinvalidchars(patientid);
    patientid = strrep(patientid,':','-');
    if contains(instring,patientid,'IgnoreCase',true)
      instring = regexprep(instring,patientid,'patientID','ignorecase');
    end
  end
catch
end
outstring = instring;

%--------------------------------
function anonymous_Callback(silent,newname)
%--------------------------------
%Makes a data set anonymous by removing
%- PatientInfo.Name
%- PatientInfo.ID
%- PatientInfo.BirthDate
%- FileName
%- OrigFileName
%- PathName

global DATA SET
if nargin<1
  silent=false;
end

if not(silent)
  str1 = dprintf('Making data set anonymous is irreversible.');
  str2 = dprintf('This feature will remove subject identity (patient ID and name) on the current .mat files.');
  str3 = dprintf('Are you sure?');
  str = sprintf('%s %s %s',str1,str2,str3);
  if not(yesno(str,[],DATA.GUI.Segment))
    return;
  end
end

if nargin<2
  %Ask for new name
  newname = myinputdlg({'Enter new name for dataset'},'New name',1,{'NoName'});
  if isempty(newname)
    myfailed('No name entered. Aborted.',DATA.GUI.Segment);
    return;
  else
    newname = newname{1};
  end
end

%Remove info
for loop=1:length(SET)
  %calculate age since birthdate is removed
  if (not(isfield(SET(loop).PatientInfo,'Age'))||isempty(SET(loop).PatientInfo.Age))&&...
      (not(isempty(SET(loop).PatientInfo.BirthDate)) &&not(isempty(SET(loop).PatientInfo.AcquisitionDate)))
    age = datenum(SET(loop).PatientInfo.AcquisitionDate,'yyyymmdd')-datenum(SET(loop).PatientInfo.BirthDate,'yyyymmdd');
    SET(loop).PatientInfo.Age = datestr(age,'yy');
  end
  %remove patientinfo
  SET(loop).PatientInfo.Name = newname;
  SET(loop).PatientInfo.ID = '';
  SET(loop).PatientInfo.BirthDate = '';
  SET(loop).FileName = newname;
  SET(loop).OrigFileName = newname;
  SET(loop).PathName = '';
  if isfield(SET(loop),'Report')
    % remove report information
    SET(loop).Report = [];
  end
end

if not(silent)
  mymsgbox('Patient identity removed. Store the file to save the changes.','Done!',DATA.GUI.Segment);
  viewfunctions('setview');
end
segment('updatetitle');


%--------------------------------
function anonymoustotal_Callback(silent,newname)
%--------------------------------
%Makes a data set completely anonymous by removing
%- PatientInfo.Name
%- PatientInfo.ID
%- PatientInfo.BirthDate
%- PatientInfo.Sex
%- PatientInfo.Age
%- PatientInfo.AcquisitionDate
%- PatientInfo.Length
%- PatientInfo.Weight
%- PatientInfo.BSA
%- PatientInfo.Institution
%- FileName
%- OrigFileName
%- PathName

global DATA SET
if nargin<1
  silent=false;
end

if not(silent)
  str1 = dprintf('Making data set anonymous is irreversible.');
  str2 = dprintf('This feature will remove all subject identity and info (patient ID, name, sex, age, length, weight and acquisition date) on the current .mat files.');
  str3 = dprintf('Are you sure?');
  str = sprintf('%s %s %s',str1,str2,str3);
  if not(yesno(str,[],DATA.GUI.Segment))
    return;
  end
end

if nargin<2
  %Ask for new name
  newname = myinputdlg({'Enter new name for dataset'},'New name',1,{'NoName'});
  if isempty(newname)
    myfailed('No name entered. Aborted.',DATA.GUI.Segment);
    return;
  else
    newname = newname{1};
  end
end

%Remove info
for loop=1:length(SET)
  %remove patient info
  SET(loop).PatientInfo.Name = newname;
  SET(loop).PatientInfo.ID = '';
  SET(loop).PatientInfo.BirthDate = '';
  SET(loop).FileName = newname;
  SET(loop).OrigFileName = newname;
  SET(loop).PathName = '';
  SET(loop).PatientInfo.Sex = '';
  SET(loop).PatientInfo.Age = '';
  SET(loop).PatientInfo.AcquisitionDate = '';
  SET(loop).PatientInfo.Length = '';
  SET(loop).PatientInfo.Weight = '';
  SET(loop).PatientInfo.BSA = '';
  SET(loop).PatientInfo.Institution = '';
  if isfield(SET(loop),'Report')
    % remove report information
    SET(loop).Report = [];
  end
end

if not(silent)
  mymsgbox('Patient identity and info removed. Store the file to save the changes.','Done!',DATA.GUI.Segment);
  viewfunctions('setview');
end
segment('updatetitle');

%--------------------------------
function anonymouscorrective_Callback(silent,newname)
%--------------------------------
%Correction of pseudonymization by removing
%- OrigFileName
%- PathName

global DATA SET
if nargin<1
  silent=false;
end

if not(silent)
  if not(yesno('This feature will make a correction of previous pseudonymization by removing original filename and pathname on the current .mat files. Are you sure?',[],DATA.GUI.Segment))
    return;
  end
end

%Remove info
for loop =1:length(SET)
  SET(loop).OrigFileName = newname;
  SET(loop).PathName = '';
  if isfield(SET(loop),'Report')
    % remove report information
    SET(loop).Report = [];
  end
end

if not(silent)
  mymsgbox('Original filename and pathname removed. Store the file to save the changes.','Done!',DATA.GUI.Segment);
  segment('viewrefresh_Callback');
end


%--------------------------------
function applylight_Callback
%--------------------------------
%Apply current light setting to current image stack. Makes
%permanent changes in the IM field.

global DATA SET NO

str1 = dprintf('This feature is irreversible.');
str2 = dprintf('Are you sure?');
str = sprintf('%s %s', str1,str2);
if yesno(str,[],DATA.GUI.Segment)

  tempnos=NO;
  imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
  if not(imissingle)
    return;
  end

  SET(NO).IM=SET(NO).IntensityMapping.Contrast*SET(NO).IM+...
    (SET(NO).IntensityMapping.Brightness-0.5);
  SET(NO).IM = min(SET(NO).IM,1);
  SET(NO).IM = max(SET(NO).IM,0);
  DATA.EndoEdgeDetected = false;
  DATA.EpiEdgeDetected = false;
  DATA.ViewIM{DATA.CurrentPanel} = [];
  segment('resetlight_Callback');
  DATA.ViewIM{DATA.CurrentPanel} = [];
  drawfunctions('drawimage',DATA.CurrentPanel);
end

%-------------------------------
function flip_Callback(dim)
%-------------------------------
%Helper function to image stack flipping tools
global DATA SET NO

oldsize = [SET(NO).XSize SET(NO).YSize SET(NO).ZSize];

%--- Find image stacks to flip
nos = NO;
if length(SET(NO).Linked) > 1
  nos = SET(NO).Linked;

  tempnos=NO;
  imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
  if not(imissingle)
    return;
  end
end

disableundo;

%--- Loop over image stacks
for no = nos

  %Restore image enhancement
  if isfield(SET(no),'IMBackup') && (~isempty(SET(no).IMBackup))
    SET(no).IM = SET(no).IMBackup;
    SET(no).IMBackup = [];
  end

  %IM
  SET(no).IM = flip(SET(no).IM,dim);

  %PapillaryIM
  if not(isempty(SET(no).PapillaryIM))
    SET(no).PapillaryIM = flip(SET(no).PapillaryIM,dim);
  end


  %Scar
  if not(isempty(SET(no).Scar))
    if dim>2
      tempdim = 3;
    else
      tempdim = dim;
    end
    SET(no).Scar.IM = flip(SET(no).Scar.IM,tempdim);
    SET(no).Scar.Auto = flip(SET(no).Scar.Auto,tempdim);
    SET(no).Scar.Result = flip(SET(no).Scar.Result,tempdim);
    SET(no).Scar.Manual = flip(SET(no).Scar.Manual,tempdim);
    SET(no).Scar.NoReflow = flip(SET(no).Scar.NoReflow,tempdim);
    SET(no).Scar.MyocardMask = flip(SET(no).Scar.MyocardMask,tempdim);
    %SET(no).Scar.Undo = flip(SET(no).Scar.Undo,tempdim);
  end

  %MaR
  if not(isempty(SET(no).MaR))
    SET(no).MaR.Auto = flip(SET(no).MaR.Auto,dim);
    SET(no).MaR.Result = flip(SET(no).MaR.Result,dim);
    SET(no).MaR.Manual = flip(SET(no).MaR.Manual,dim);
    SET(no).MaR.MyocardMask = flip(SET(no).MaR.MyocardMask,dim);
  end

  %LevelSet
  if not(isempty(SET(no).LevelSet)) && isfield(SET(no).LevelSet,'Object')
    mywarning('Close general segmentation view before continuing.',DATA.GUI.Segment);

    O = SET(no).LevelSet.Object;

    % stored objects, important to flip stored objects
    for obj = 1:O.getnumobjects()

      if O.iscompressed(obj)
        O.uncompress(obj,oldsize);
      end

      temp = O.getbw(obj);

      if ~isempty(temp)
        temp = flip(temp,dim);
        O.setbw(obj,temp);
        O.setcut(obj,[]);
        O.setdeform(obj,[]);
        O.compress(obj);
      end

    end 

  end

  %Fix with segmentation
  for zloop=1:SET(no).ZSize
    for tloop=1:SET(no).TSize
      segment('checkconsistency',SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
    end
  end

end

%Flip image coordinates
switch dim
  case 1
    %x
    SET(no).ImagePosition = SET(no).ImagePosition+(SET(no).XSize-1)*SET(no).ResolutionX*SET(no).ImageOrientation(4:6);
    SET(no).ImageOrientation(4:6) = -SET(no).ImageOrientation(4:6); %Negate direction;
  case 2
    %y
    SET(no).ImagePosition = SET(no).ImagePosition+(SET(no).YSize-1)*SET(no).ResolutionY*SET(no).ImageOrientation(1:3);
    SET(no).ImageOrientation(1:3) = -SET(no).ImageOrientation(1:3); %Negate direction;
  case 3
    %t
    %Do nothing!
  case 4
    %z
    zdir = cross(...
      SET(no).ImageOrientation(1:3),...
      SET(no).ImageOrientation(4:6));
    %Update position, zdir will be automatically negated since x is always
    %flipped at the same time as z is flipped.
    SET(no).ImagePosition = SET(no).ImagePosition+zdir*(SET(no).ZSize-1)*...
      (SET(no).SliceThickness+SET(no).SliceGap);
end

%Negate velocity data if necessary
if not(isempty(SET(no).Flow))
  switch dim
    case 1
      if isequal(SET(no).Flow.PhaseX,no)
        SET(no).IM = 1-SET(no).IM; %Negate
      end
    case 2
      if isequal(SET(no).Flow.PhaseY,no)
        SET(no).IM = 1-SET(no).IM; %Negate
      end
    case 3
      if isequal(SET(no).Flow.PhaseNo,no)
        SET(no).IM = 1-SET(no).IM; %Negate
      end
    case 4
      %Negate all since reversing time
      if isequal(SET(no).Flow.PhaseX,no)
        SET(no).IM = 1-SET(no).IM; %Negate
      end
      if isequal(SET(no).Flow.PhaseY,no)
        SET(no).IM = 1-SET(no).IM; %Negate
      end
      if isequal(SET(no).Flow.PhaseNo,no)
        SET(no).IM = 1-SET(no).IM; %Negate
      end
  end

end

%remove Strain tagging analysis
if not(isempty(SET(no).StrainTagging))
  logdisp('Flip image stack clears the Strain tagging analysis',true);
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  SET(no).StrainTagging = [];
end
%remove Strain MITT analysis
logmsg = 'Flip image stack clears the Strain MITT analysis';
strainmitt.strainmitt('removestrainanalysis',no,logmsg);

%Graphical update
DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel)

DATA.EndoEdgeDetected = false;
DATA.EpiEdgeDetected = false;

%---------------------------
function flipx_Callback
%---------------------------
%Flip x direction of current image stack.
%Need to flip both x and z to maintain a righthand system.
flipx_helper;
flipz_helper;
segment('updatevolume');
viewfunctions('setview') %this is the refresh call

%Create thumbnails.
recreate = true;
sliderupdated = false;
drawfunctions('drawthumbnails',recreate,sliderupdated);

%---------------------------
function flipy_Callback
%---------------------------
%Flip y directon of current image stack.
%Need to flip both y and z to maintain a righthand system.
flipy_helper;
flipz_helper;
segment('updatevolume');
viewfunctions('setview') %this is the refresh call

%Create thumbnails.
recreate = true;
sliderupdated = false;
drawfunctions('drawthumbnails',recreate,sliderupdated);

%---------------------------
function flipz_Callback
%---------------------------
%Flip z direction of current image stack.
%Need to flip both z and x to maintain a righthand system.
flipx_helper;
flipz_helper;
segment('updatevolume');
viewfunctions('setview') %this is the refresh call

%-----------------------------------
function rotate90right_Callback
%-----------------------------------
%Rotate current image stack 90 degrees right. Currently not working
%properly.

flipxy_helper;
flipy_helper;
segment('updatevolume');
viewfunctions('setview') %this is the refresh call

%Create thumbnails.
recreate = true;
sliderupdated = false;
drawfunctions('drawthumbnails',recreate,sliderupdated);

%--------------------------
function flipx_helper
%--------------------------
%Helper function to flip in x direction. Takes care of segmentation.
global SET NO

%Make sure magnitude image is current.
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
  viewfunctions('switchimagestack',no);
end

if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX = (SET(NO).XSize+1)-SET(NO).EndoX;
  %maintain rotation dir
  SET(NO).EndoX = flip(SET(NO).EndoX,1);
  SET(NO).EndoY = flip(SET(NO).EndoY,1);
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX = (SET(NO).XSize+1)-SET(NO).EpiX;
  %maintain rotation dir
  SET(NO).EpiX = flip(SET(NO).EpiX,1);
  SET(NO).EpiY = flip(SET(NO).EpiY,1);
end

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX = (SET(NO).XSize+1)-SET(NO).RVEndoX;
  %maintain rotation dir
  SET(NO).RVEndoX = flip(SET(NO).RVEndoX,1);
  SET(NO).RVEndoY = flip(SET(NO).RVEndoY,1);
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX = (SET(NO).XSize+1)-SET(NO).RVEpiX;
  %maintain rotation dir
  SET(NO).RVEpiX = flip(SET(NO).RVEpiX,1);
  SET(NO).RVEpiY = flip(SET(NO).RVEpiY,1);
end

if ~isempty(SET(NO).LA)
  SET(NO).LA.X = (SET(NO).XSize+1)-SET(NO).LA.X;
  %maintain rotation dir
  SET(NO).LA.X = flip(SET(NO).LA.X,1);
  SET(NO).LA.Y = flip(SET(NO).LA.Y,1);
end

if ~isempty(SET(NO).RA)
  SET(NO).RA.X = (SET(NO).XSize+1)-SET(NO).RA.X;
  %maintain rotation dir
  SET(NO).RA.X = flip(SET(NO).RA.X,1);
  SET(NO).RA.Y = flip(SET(NO).RA.Y,1);
end


for zloop=1:SET(NO).ZSize
  for tloop=1:SET(NO).TSize
    if ~isempty(SET(NO).EndoPinX)
      SET(NO).EndoPinX{tloop,zloop} = (SET(NO).XSize+1)-SET(NO).EndoPinX{tloop,zloop};
    end
    if ~isempty(SET(NO).EpiPinX)
      SET(NO).EpiPinX{tloop,zloop} = (SET(NO).XSize+1)-SET(NO).EpiPinX{tloop,zloop};
    end
  end
end

%Forget about trying to move interp points
segmentation('removeallinterp_Callback',true);

%roi
if SET(NO).RoiN>0
  for rloop=1:SET(NO).RoiN
    SET(NO).Roi(rloop).X = (SET(NO).XSize+1)-SET(NO).Roi(rloop).X;
  end
end

%point
if not(isempty(SET(NO).Point))
  SET(NO).Point.X = SET(NO).XSize-SET(NO).Point.X+1;
end

%measure
for loop=1:length(SET(NO).Measure)
  SET(NO).Measure(loop).X = SET(NO).XSize-SET(NO).Measure(loop).X+1;
end

%Flip image data and update
flip_Callback(1);

%-------------------------
function flipy_helper
%-------------------------
%Helper function to flip in y direction. Takes care of segmentation.
global SET NO

%Make sure magnitude image is current.
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
  viewfunctions('switchimagestack',no)
end

if ~isempty(SET(NO).EndoX)
  SET(NO).EndoY = (SET(NO).YSize+1)-SET(NO).EndoY;
  %maintain rotation dir
  SET(NO).EndoX = flip(SET(NO).EndoX,1);
  SET(NO).EndoY = flip(SET(NO).EndoY,1);
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiY = (SET(NO).YSize+1)-SET(NO).EpiY;
  %maintain rotation dir
  SET(NO).EpiX = flip(SET(NO).EpiX,1);
  SET(NO).EpiY = flip(SET(NO).EpiY,1);
end

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoY = (SET(NO).YSize+1)-SET(NO).RVEndoY;
  %maintain rotation dir
  SET(NO).RVEndoX = flip(SET(NO).RVEndoX,1);
  SET(NO).RVEndoY = flip(SET(NO).RVEndoY,1);
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiY = (SET(NO).YSize+1)-SET(NO).RVEpiY;
  %maintain rotation dir
  SET(NO).RVEpiX = flip(SET(NO).RVEpiX,1);
  SET(NO).RVEpiY = flip(SET(NO).RVEpiY,1);
end

if ~isempty(SET(NO).LA)
  SET(NO).LA.Y = (SET(NO).YSize+1)-SET(NO).LA.Y;
  %maintain rotation dir
  SET(NO).LA.X = flip(SET(NO).LA.X,1);
  SET(NO).LA.Y = flip(SET(NO).LA.Y,1);
end
if ~isempty(SET(NO).RA)
  SET(NO).RA.Y = (SET(NO).YSize+1)-SET(NO).RA.Y;
  %maintain rotation dir
  SET(NO).RA.X = flip(SET(NO).RA.X,1);
  SET(NO).RA.Y = flip(SET(NO).RA.Y,1);
end

for zloop=1:SET(NO).ZSize
  for tloop=1:SET(NO).TSize
    if ~isempty(SET(NO).EndoPinX)
      SET(NO).EndoPinY{tloop,zloop} = (SET(NO).YSize+1)-SET(NO).EndoPinY{tloop,zloop};
    end
    if ~isempty(SET(NO).EndoPinY)
      SET(NO).EpiPinY{tloop,zloop} = (SET(NO).YSize+1)-SET(NO).EpiPinY{tloop,zloop};
    end
  end
end

%Forget about trying to move interp points
segmentation('removeallinterp_Callback',true);

%roi
if SET(NO).RoiN>0
  for rloop=1:SET(NO).RoiN
    SET(NO).Roi(rloop).Y = (SET(NO).YSize+1)-SET(NO).Roi(rloop).Y;
  end
end

%point
if not(isempty(SET(NO).Point))
  SET(NO).Point.Y = SET(NO).YSize-SET(NO).Point.Y+1;
end

%measure
for loop=1:length(SET(NO).Measure)
  SET(NO).Measure(loop).Y = SET(NO).YSize-SET(NO).Measure(loop).Y+1;
end

%rotation center
if SET(NO).Rotated
  SET(NO).RotationCenter = SET(NO).YSize+1-SET(NO).RotationCenter;
end

%Flip image data and update
flip_Callback(2);

%-------------------------
function flipz_helper
%-------------------------
%Helper function to flip in z direction. Takes care of segmentation.

global SET NO

%Make sure magnitude image is current.
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
  viewfunctions('switchimagestack',no)
end

SET(NO).CurrentSlice = SET(NO).ZSize-SET(NO).CurrentSlice+1;

ind = SET(NO).ZSize:-1:1;

if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX = SET(NO).EndoX(:,:,ind);
  SET(NO).EndoY = SET(NO).EndoY(:,:,ind);
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX = SET(NO).EpiX(:,:,ind);
  SET(NO).EpiY = SET(NO).EpiY(:,:,ind);
end

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX = SET(NO).RVEndoX(:,:,ind);
  SET(NO).RVEndoY = SET(NO).RVEndoY(:,:,ind);
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX = SET(NO).RVEpiX(:,:,ind);
  SET(NO).RVEpiY = SET(NO).RVEpiY(:,:,ind);
end

if ~isempty(SET(NO).LA)
  SET(NO).LA.X = SET(NO).LA.X(:,:,ind);
  SET(NO).LA.Y = SET(NO).LA.Y(:,:,ind);
end
if ~isempty(SET(NO).RA)
  SET(NO).RA.X = SET(NO).RA.X(:,:,ind);
  SET(NO).RA.Y = SET(NO).RA.Y(:,:,ind);
end

SET(NO).EndoDraged = SET(NO).EndoDraged(:,ind);
SET(NO).EpiDraged = SET(NO).EpiDraged(:,ind);

if ~isempty(SET(NO).EndoPinX)
  SET(NO).EndoPinX = SET(NO).EndoPinX(:,ind);
  SET(NO).EndoPinY = SET(NO).EndoPinY(:,ind);
end
if ~isempty(SET(NO).EpiPinY)
  SET(NO).EpiPinX = SET(NO).EpiPinX(:,ind);
  SET(NO).EpiPinY = SET(NO).EpiPinY(:,ind);
end
if ~isempty(SET(NO).RVEndoPinX)
  SET(NO).RVEndoPinX = SET(NO).RVEndoPinX(:,ind);
  SET(NO).RVEndoPinY = SET(NO).RVEndoPinY(:,ind);
end
if ~isempty(SET(NO).RVEpiPinY)
  SET(NO).RVEpiPinX = SET(NO).RVEpiPinX(:,ind);
  SET(NO).RVEpiPinY = SET(NO).RVEpiPinY(:,ind);
end

%Forget about trying to move interp points
segmentation('removeallinterp_Callback',true);

%roi
if SET(NO).RoiN>0
  for rloop=1:SET(NO).RoiN
    SET(NO).Roi(rloop).Z = SET(NO).ZSize-SET(NO).Roi(rloop).Z+1;
  end
end

%point
if not(isempty(SET(NO).Point))
  SET(NO).Point.Z = SET(NO).ZSize-SET(NO).Point.Z+1;
end

%measure
for loop=1:length(SET(NO).Measure)
  SET(NO).Measure(loop).Z = SET(NO).ZSize-SET(NO).Measure(loop).Z+1;
end

%Flip image data and update
flip_Callback(4);

%---------------------------
function flipt_Callback
%---------------------------
%Helper function to flip in t direction. Takes care of segmentation.

global SET NO

%Make sure magnitude image is current.
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
  viewfunctions('switchimagestack',no)

end

ind = SET(NO).TSize:-1:1;

if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX = SET(NO).EndoX(:,ind,:);
  SET(NO).EndoY = SET(NO).EndoY(:,ind,:);
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX = SET(NO).EpiX(:,ind,:);
  SET(NO).EpiY = SET(NO).EpiY(:,ind,:);
end

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX = SET(NO).RVEndoX(:,ind,:);
  SET(NO).RVEndoY = SET(NO).RVEndoY(:,ind,:);
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX = SET(NO).RVEpiX(:,ind,:);
  SET(NO).RVEpiY = SET(NO).RVEpiY(:,ind,:);
end

if ~isempty(SET(NO).LA.X)
  SET(NO).LA.X = SET(NO).LA.X(:,ind,:);
  SET(NO).LA.Y = SET(NO).LA.Y(:,ind,:);
end

if ~isempty(SET(NO).RA.X)
  SET(NO).RA.X = SET(NO).RA.X(:,ind,:);
  SET(NO).RA.Y = SET(NO).RA.Y(:,ind,:);
end

SET(NO).EndoDraged = SET(NO).EndoDraged(ind,:);
SET(NO).EpiDraged = SET(NO).EpiDraged(ind,:);
SET(NO).EDT = SET(NO).TSize+1-SET(NO).EDT;
SET(NO).EST = SET(NO).TSize+1-SET(NO).EST;

if ~isempty(SET(NO).EndoPinX)
  SET(NO).EndoPinX = SET(NO).EndoPinX(ind,:);
  SET(NO).EndoPinY = SET(NO).EndoPinY(ind,:);
end
if ~isempty(SET(NO).EpiPinX)
  SET(NO).EpiPinX = SET(NO).EpiPinX(ind,:);
  SET(NO).EpiPinY = SET(NO).EpiPinY(ind,:);
end
if ~isempty(SET(NO).RVEndoPinX)
  SET(NO).RVEndoPinX = SET(NO).RVEndoPinX(ind,:);
  SET(NO).RVEndoPinY = SET(NO).RVEndoPinY(ind,:);
end
if ~isempty(SET(NO).EpiPinX)
  SET(NO).RVEpiPinX = SET(NO).RVEpiPinX(ind,:);
  SET(NO).RVEpiPinY = SET(NO).RVEpiPinY(ind,:);
end

%Forget about trying to move interp points
segmentation('removeallinterp_Callback',true);

if SET(NO).RoiN>0
  for rloop=1:SET(NO).RoiN
    SET(NO).Roi(rloop).X = SET(NO).Roi(rloop).X(:,ind);
    SET(NO).Roi(rloop).Y = SET(NO).Roi(rloop).Y(:,ind);
    SET(NO).Roi(rloop).T = find(not(isnan(SET(NO).Roi(rloop).X(1,:))));
  end
end

%point
if not(isempty(SET(NO).Point))
  SET(NO).Point.T = SET(NO).TSize-SET(NO).Point.T+1;
end

%measure, now timeresolved
if not(isempty(SET(NO).Measure))
  SET(NO).Measure.T = SET(NO).TSize-SET(NO).Measure.T+1; %wors for NaN
end

%Flip image data and update
flip_Callback(3);
viewfunctions('setview');

%----------------------------
function flipzt_Callback
%----------------------------
%Interchange z & t for current image stack. Takes care of segmentation.
global DATA SET NO

if not(isempty(SET(NO).Flow))
  mywarning('Flip z-t not implemented for flow image stacks. Select view one and flip all connected image stacks manually.',DATA.GUI.Segment);
end

oldsize = [SET(NO).XSize SET(NO).YSize SET(NO).OrgZSize];

%Restore image enhancement
if isfield(SET(no),'IMBackup') && (~isempty(SET(no).IMBackup))
  SET(no).IM = SET(no).IMBackup;
  SET(no).IMBackup = [];
end

permorder = [1 2 4 3];
SET(NO).IM = permute(SET(NO).IM,permorder);
SET(NO).TSize = size(SET(NO).IM,3);
SET(NO).ZSize = size(SET(NO).IM,4);
SET(NO).LVV = nan(1,SET(NO).TSize);
SET(NO).PV = zeros(1,SET(NO).TSize);
SET(NO).EPV = SET(NO).LVV;
SET(NO).LVM = SET(NO).LVV;
SET(NO).PFR = 0;
SET(NO).PER = 0;
SET(NO).PFRT = 1;
SET(NO).PERT = 1;
SET(NO).ESV = 0;
SET(NO).EDV = 0;
SET(NO).EF = 0;
SET(NO).SV = 0;
SET(NO).RVV = SET(NO).LVV;
SET(NO).EPV = SET(NO).RVV;
SET(NO).RVM = 0;
SET(NO).RVEDV = 0;
SET(NO).RVESV = 0;
SET(NO).RVEF = 0;

if isequal(SET(NO).TIncr,0)
  SET(NO).TIncr = 1/SET(NO).TSize;
  mywarning('Time increment was 0, now set so that R-R is one second.',DATA.GUI.Segment);
end

SET(NO).TimeVector = linspace(0,(SET(NO).TSize-1)*SET(NO).TIncr,SET(NO).TSize); % 180502: Added Update of timevector when switching z and t dimensions
SET(NO).CurrentTimeFrame = 1;
SET(NO).CurrentSlice = 1;
temp = SET(NO).OrgTSize;
SET(NO).OrgTSize = SET(NO).OrgZSize;
SET(NO).OrgZSize = temp;

%endox,endoy,epix,epix
permorder = [1 3 2];  %orig [n t z]
if ~isempty(SET(NO).EndoX)
  SET(NO).EndoX = permute(SET(NO).EndoX,permorder);
  SET(NO).EndoY = permute(SET(NO).EndoY,permorder);
end

if ~isempty(SET(NO).EpiX)
  SET(NO).EpiX = permute(SET(NO).EpiX,permorder);
  SET(NO).EpiY = permute(SET(NO).EpiY,permorder);
end

if ~isempty(SET(NO).RVEndoX)
  SET(NO).RVEndoX = permute(SET(NO).RVEndoX,permorder);
  SET(NO).RVEndoY = permute(SET(NO).RVEndoY,permorder);
end

if ~isempty(SET(NO).RVEpiX)
  SET(NO).RVEpiX = permute(SET(NO).RVEpiX,permorder);
  SET(NO).RVEpiY = permute(SET(NO).RVEpiY,permorder);
end

if ~isempty(SET(NO).LA.X)
  SET(NO).LA.X = permute(SET(NO).LA.X,permorder);
  SET(NO).LA.Y = permute(SET(NO).LA.Y,permorder);
end

if ~isempty(SET(NO).RA.X)
  SET(NO).RA.X = permute(SET(NO).RA.X,permorder);
  SET(NO).RA.Y = permute(SET(NO).RA.Y,permorder);
end

SET(NO).EndoDraged = SET(NO).EndoDraged';
SET(NO).EpiDraged = SET(NO).EpiDraged';

%pins
permorder = [2 1]; %orig [t z]

if ~isempty(SET(NO).EndoPinX)
  SET(NO).EndoPinX = permute(SET(NO).EndoPinX,permorder);
  SET(NO).EndoPinY = permute(SET(NO).EndoPinY,permorder);
end

if ~isempty(SET(NO).EpiPinX)
  SET(NO).EpiPinX = permute(SET(NO).EpiPinX,permorder);
  SET(NO).EpiPinY = permute(SET(NO).EpiPinY,permorder);
end

if ~isempty(SET(NO).RVEndoPinX)
  SET(NO).RVEndoPinX = permute(SET(NO).RVEndoPinX,permorder);
  SET(NO).RVEndoPinY = permute(SET(NO).RVEndoPinY,permorder);
end

if ~isempty(SET(NO).RVEpiPinX)
  SET(NO).RVEpiPinX = permute(SET(NO).RVEpiPinX,permorder);
  SET(NO).RVEpiPinY = permute(SET(NO).RVEpiPinY,permorder);
end

%Forget about trying to move interp points
segmentation('removeallinterp_Callback',true);

if not(isempty(SET(NO).LevelSet)) && isfield(SET(NO).LevelSet,'Object')
  permorder = [1 2 4 3];
  mywarning('Close general segmentation view before continuing.',DATA.GUI.Segment);
  % stored objects, important to flip stored objects

  O = SET(NO).LevelSet.Object;

  for obj = 1:O.getnumobjects()

    if O.iscompressed(obj)
      O.uncompress(obj,oldsize);
    end

    temp = O.getbw(obj);

    if ~isempty(temp)
      temp = permute(temp,permorder);
      O.setbw(obj,temp);
      O.setcut(obj,[]);
      O.setdeform(obj,[]);
      O.compress(obj);
    end
  end

  [r,g,b] = segment3dp.tools('xyz2rgb',SET(NO).XSize,SET(NO).YSize,SET(NO).ZSize);

  %slice
  SET(NO).LevelSet.View.RSlice = round(r/2);
  SET(NO).LevelSet.View.GSlice = round(g/2);
  SET(NO).LevelSet.View.BSlice = round(b/2);

  %ZoomState
  SET(NO).LevelSet.View.RZoomState = [];
  SET(NO).LevelSet.View.BZoomState = [];
  SET(NO).LevelSet.View.GZoomState = [];
end

%roi
if SET(NO).RoiN>0
  mywarning('Flip in z-t of ROI''s not implemented.',DATA.GUI.Segment);
  roi('roiclearall_Callback');
end

%point
if not(isempty(SET(NO).Point))
  for loop=1:length(SET(NO).Point.X)
    if isnan(SET(NO).Point.T(loop))
      mywarning('Non time resolved point detected. Ignored.',DATA.GUI.Segment);
    else
      [SET(NO).Point.T(loop),SET(NO).Point.Z(loop)] = flipvars(...
        SET(NO).Point.T(loop),SET(NO).Point.Z(loop));
    end
  end
end

%measure now time resolved JU
if not(isempty(SET(NO).Measure))
  for loop=1:length(SET(NO).Measure)
    if isnan(SET(NO).Measure(loop).T)
      mywarning('Non time resolved measure detected. Ignored.',DATA.GUI.Segment);
    else
      [SET(NO).Measure(loop).T,SET(NO).Measure(loop).Z] = flipvars(...
        SET(NO).Measure(loop).T,SET(NO).Measure(loop).Z);
    end
  end
  %  measureclearall_Callback;
end

if not(isempty(SET(NO).Scar))
  myfailed('Flip z-t is not implemented with scar data present. Scar data removed.',DATA.GUI.Segment);
  SET(NO).Scar = [];
end

segment('updatevolume');
DATA.ViewIM{DATA.CurrentPanel} = [];
viewfunctions('setview');

%--------------------------
function flipxy_helper
%--------------------------
%Interchange z & t of current image stack. Takes care of segmentation.
global DATA SET NO

%--- Find image stacks to flip
nos = SET(NO).Linked;

oldsize = [SET(NO).XSize SET(NO).YSize SET(NO).ZSize];

%--- Loop over image stacks
for no = nos

  %Restore image enhancement
  if isfield(SET(no),'IMBackup') && (~isempty(SET(no).IMBackup))
    SET(no).IM = SET(no).IMBackup;
    SET(no).IMBackup = [];
  end

  SET(no).IM = permute(SET(no).IM,[2 1 3 4]);
  SET(no).XSize = size(SET(no).IM,1);
  SET(no).YSize = size(SET(no).IM,2);
  temp = SET(no).ResolutionX;
  SET(no).ResolutionX = SET(no).ResolutionY;
  SET(no).ResolutionY = temp;
  temp = SET(no).OrgXSize;
  SET(no).OrgXSize = SET(no).OrgYSize;
  SET(no).OrgYSize = temp;
  temp = SET(no).XMin;
  SET(no).XMin = SET(no).YMin;
  SET(no).YMin = temp;

  %endo,epi,roi
  if ~isempty(SET(NO).EndoX)
    [SET(no).EndoX,SET(no).EndoY] = flipvars(SET(no).EndoX,SET(no).EndoY);
  end

  if ~isempty(SET(NO).EpiX)
    [SET(no).EpiX,SET(no).EpiY] = flipvars(SET(no).EpiX,SET(no).EpiY);
  end

  if ~isempty(SET(NO).RVEndoX)
    [SET(no).RVEndoX,SET(no).RVEndoY] = flipvars(SET(no).RVEndoX,SET(no).RVEndoY);
  end

  if ~isempty(SET(NO).RVEpiX)
    [SET(no).RVEpiX,SET(no).RVEpiY] = flipvars(SET(no).RVEpiX,SET(no).RVEpiY);
  end

  if ~isempty(SET(NO).LA)
    [SET(NO).LA.X,SET(NO).LA.Y] = flipvars(SET(NO).LA.X,SET(NO).LA.Y);
  end

  if ~isempty(SET(NO).RA)
    [SET(NO).RA.X,SET(NO).RA.Y] = flipvars(SET(NO).RA.X,SET(NO).RA.Y);
  end


  [SET(no).EndoPinX,SET(no).EndoPinY] = flipvars(SET(no).EndoPinX,SET(no).EndoPinY);
  [SET(no).EpiPinX,SET(no).EpiPinY] = flipvars(SET(no).EpiPinX,SET(no).EpiPinY);
  [SET(no).RVEndoPinX,SET(no).RVEndoPinY] = flipvars(SET(no).RVEndoPinX,SET(no).RVEndoPinY);
  [SET(no).RVEpiPinX,SET(no).RVEpiPinY] = flipvars(SET(no).RVEpiPinX,SET(no).RVEpiPinY);

  %Forget about trying to move interp points
  segmentation('removeallinterp_Callback',true);

  for rloop=1:SET(no).RoiN
    [SET(no).Roi(rloop).X,SET(no).Roi(rloop).Y] = flipvars(SET(no).Roi(rloop).X,SET(no).Roi(rloop).Y);
  end

  %Scar
  if not(isempty(SET(no).Scar))
    permorder = [2 1 3];
    SET(no).Scar.IM = permute(SET(no).Scar.IM,permorder);
    SET(no).Scar.Auto = permute(SET(no).Scar.Auto,permorder);
    SET(no).Scar.Result = permute(SET(no).Scar.Result,permorder);
    SET(no).Scar.Manual = permute(SET(no).Scar.Manual,permorder);
    %SET(no).Scar.Undo = permute(SET(no).Scar.Undo,permorder);
    SET(no).Scar.MyocardMask = permute(SET(no).Scar.MyocardMask,permorder);
    SET(no).Scar.NoReflow = permute(SET(no).Scar.NoReflow,permorder);
  end

  %MaR
  if not(isempty(SET(no).MaR))
    permorder = [2 1 3 4];
    SET(no).MaR.Auto = permute(SET(no).MaR.Auto,permorder);
    SET(no).MaR.Result = permute(SET(no).MaR.Result,permorder);
    SET(no).MaR.Manual = permute(SET(no).MaR.Manual,permorder);
    SET(no).MaR.MyocardMask = permute(SET(no).MaR.MyocardMask,permorder);
  end

  %LevelSet
  if not(isempty(SET(no).LevelSet))
    permorder = [2 1 3 4];
    mywarning('Close general segmentation view before continuing.',DATA.GUI.Segment);

    O = SET(no).LevelSet.Object;

    % stored objects, important to flip stored objects
    for obj = 1:O.getnumobjects()

      if O.iscompressed(obj)
        O.uncompress(obj,oldsize);
      end

      temp = O.getbw(obj);

      if ~isempty(temp)
        temp = permute(temp,permorder);
        O.setbw(obj,temp);
        O.setcut(obj,[]);
        O.setdeform(obj,[]);
        O.compress(obj);
      end

    end

    [r,g,b] = segment3dp.tools('xyz2rgb',SET(no).XSize,SET(no).YSize,SET(no).ZSize);

    %slice
    SET(no).LevelSet.View.RSlice = round(r/2);
    SET(no).LevelSet.View.GSlice = round(g/2);
    SET(no).LevelSet.View.BSlice = round(b/2);

    %ZoomState
    SET(no).LevelSet.View.RZoomState = [];
    SET(no).LevelSet.View.BZoomState = [];
    SET(no).LevelSet.View.GZoomState = [];
  end

  %point,measure
  if not(isempty(SET(no).Point))
    [SET(no).Point.X,SET(no).Point.Y] = flipvars(SET(no).Point.X,SET(no).Point.Y);
  end
  for loop=1:length(SET(no).Measure)
    [SET(no).Measure(loop).X,SET(no).Measure(loop).Y] = flipvars(SET(no).Measure(loop).X,SET(no).Measure(loop).Y);
  end

  %Flip phase directions
  if not(isempty(SET(no).Flow))
    temp = SET(no).Flow.PhaseX;
    SET(no).Flow.PhaseX = SET(no).Flow.PhaseY;
    SET(no).Flow.PhaseY = temp;
  end

  %Coordinates
  %SET(no).ImagePosition = SET(no).ImagePosition+...
  %  SET(no).ImageOrientation(4:6)*(SET(no).XSize-1)*SET(no).ResolutionX+...
  %  SET(no).ImageOrientation(1:3)*(SET(no).YSize-1)*SET(no).ResolutionY;
  orient = SET(no).ImageOrientation;
  SET(no).ImageOrientation(1:3) = orient(4:6);
  SET(no).ImageOrientation(4:6) = orient(1:3);

end %loop over image stacks

DATA.ViewIM{DATA.CurrentPanel} = [];

%----------------------------
function flipxz_Callback
%----------------------------
%Interchange z & x of current image stack. Takes care of segmentation.
global DATA SET NO

if ~yesno('Warning, this will clear measurement and segmentation if present, (except general segmentation). Are you sure?',[],DATA.GUI.Segment)
  return;
end

if ~isempty(SET(NO).Flow)
  mywarning('This is not implemented for flow images, you need to run it on each flow component',DATA.GUI.Segment);
end

if (SET(NO).SliceThickness+SET(NO).SliceGap)<0
  myfailed('Zero or negative slice thickness.',DATA.GUI.Segment);
  return;
end

oldsize = [SET(NO).XSize SET(NO).YSize SET(NO).ZSize];

%Restore image enhancement
if isfield(SET(NO),'IMBackup') && (~isempty(SET(NO).IMBackup))
  SET(NO).IM = SET(NO).IMBackup;
  SET(NO).IMBackup = [];
end

SET(NO).IM = permute(SET(NO).IM,[4 2 3 1]);
SET(NO).IM = flip(SET(NO).IM,1);
SET(NO).XSize = size(SET(NO).IM,1);
SET(NO).ZSize = size(SET(NO).IM,4);

SET(NO).CurrentTimeFrame = 1;
SET(NO).CurrentSlice = 1;
temp = SET(NO).ResolutionX;
SET(NO).ResolutionX = SET(NO).SliceThickness+SET(NO).SliceGap;
SET(NO).SliceThickness = temp;
SET(NO).SliceGap = 0;
temp = SET(NO).OrgXSize;
SET(NO).OrgXSize = SET(NO).OrgZSize;
SET(NO).OrgZSize = temp;

%Clear segmentation
callbackfunctions('segmentclearall_Callback',true);
SET(NO).Scar = [];
%LevelSet is not cleared and thuse need to be flipped
if not(isempty(SET(NO).LevelSet)) && isfield(SET(no).LevelSet,'Object')
  permorder = [4 2 3 1];
  mywarning('Close general segmentation view before continuing.',DATA.GUI.Segment);

  O = SET(NO).LevelSet.Object;

  % stored objects, important to flip stored objects
  for obj = 1:O.getnumobjects()

    if O.iscompressed(obj)
      O.uncompress(obj,oldsize);
    end

    temp = O.getbw(obj);

    if ~isempty(temp)
      temp = permute(temp,permorder);
      temp = flip(temp,1);
      O.setbw(obj,temp);
      O.setcut(obj,[]);
      O.setdeform(obj,[]);
      O.compress(obj);
    end

  end

  [r,g,b] = segment3dp.tools('xyz2rgb',SET(NO).XSize,SET(NO).YSize,SET(NO).ZSize);

  %slice
  SET(NO).LevelSet.View.RSlice = round(r/2);
  SET(NO).LevelSet.View.GSlice = round(g/2);
  SET(NO).LevelSet.View.BSlice = round(b/2);

  %ZoomState
  SET(NO).LevelSet.View.RZoomState = [];
  SET(NO).LevelSet.View.BZoomState = [];
  SET(NO).LevelSet.View.GZoomState = [];
  
  segment3dp.graphics('updatergbsliders');
end

%point
if not(isempty(SET(NO).Point))
  [SET(NO).Point.X,SET(NO).Point.Z] = flipvars(SET(NO).Point.X,SET(NO).Point.Z);
  SET(NO).Point.Z = round(SET(NO).Point.Z);
  SET(NO).Point.X = SET(NO).XSize - SET(NO).Point.X + 1;
end
%delete measure
SET(NO).Measure=[];

%Coordinates
ydir = SET(NO).ImageOrientation(1:3);
zdir = cross(...
  SET(NO).ImageOrientation(1:3),...
  SET(NO).ImageOrientation(4:6));
oldpos = SET(NO).ImagePosition;
%x^~ = zdir;
%y^~ = ydir;
%z^~ = xdir;
SET(NO).ImageOrientation(4:6) = zdir;
SET(NO).ImageOrientation(1:3) = ydir;

%z^~ = xdir;
SET(NO).ImagePosition = oldpos....
  -zdir*(SET(NO).XSize-1)*SET(NO).ResolutionX;

SET(NO).NormalZoomState = [];
SET(NO).MontageZoomState = [];
SET(NO).MontageRowZoomState = [];

DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimage',DATA.CurrentPanel);

%Create thumbnails.
recreate = true;
sliderupdated = false;
drawfunctions('drawthumbnails',recreate,sliderupdated);

%----------------------
function enableundo(no)
%----------------------
%Enables undo function in menus and icons. Also
%copies old segmentation data so undo i possible
%This should allways be called before routines that
%change the segmentation
global DATA SET NO

segment('stopmovie_Callback');

if (nargin==0) || isempty(no)
  no = NO;
end

% write
if ~isempty(SET(no).Parent)
  no = SET(no).Parent;
end

DATA.NeedToSave = true;

%Increase pointer
try
  if ~isempty(DATA.Undo)
    DATA.UndoN(no) = DATA.UndoN(no)+1;
  else
    DATA.UndoN = zeros(1,length(SET));
    DATA.UndoN(no) = 1;
  end
catch %#ok<CTCH>
  DATA.UndoN = zeros(1,length(SET));
  DATA.UndoN(no) = 1;
end

%Check if larger than UndoHistory
if DATA.UndoN(no)>DATA.Pref.UndoHistory
  DATA.Undo(no,1:(DATA.Pref.UndoHistory-1)) = DATA.Undo(no,2:end);
  DATA.UndoN(no) = DATA.Pref.UndoHistory;
end

%Store
DATA.Undo(no,DATA.UndoN(no)).EndoXBackup = SET(no).EndoX;
DATA.Undo(no,DATA.UndoN(no)).EndoYBackup = SET(no).EndoY;
DATA.Undo(no,DATA.UndoN(no)).RVEndoXBackup = SET(no).RVEndoX;
DATA.Undo(no,DATA.UndoN(no)).RVEndoYBackup = SET(no).RVEndoY;
DATA.Undo(no,DATA.UndoN(no)).EpiXBackup = SET(no).EpiX;
DATA.Undo(no,DATA.UndoN(no)).EpiYBackup = SET(no).EpiY;
DATA.Undo(no,DATA.UndoN(no)).RVEpiXBackup = SET(no).RVEpiX;
DATA.Undo(no,DATA.UndoN(no)).RVEpiYBackup = SET(no).RVEpiY;
DATA.Undo(no,DATA.UndoN(no)).EndoInterpXBackup = SET(no).EndoInterpX;
DATA.Undo(no,DATA.UndoN(no)).EndoInterpYBackup = SET(no).EndoInterpY;
DATA.Undo(no,DATA.UndoN(no)).EpiInterpXBackup = SET(no).EpiInterpX;
DATA.Undo(no,DATA.UndoN(no)).EpiInterpYBackup = SET(no).EpiInterpY;
DATA.Undo(no,DATA.UndoN(no)).RVEndoInterpXBackup = SET(no).RVEndoInterpX;
DATA.Undo(no,DATA.UndoN(no)).RVEndoInterpYBackup = SET(no).RVEndoInterpY;
DATA.Undo(no,DATA.UndoN(no)).RVEpiInterpXBackup = SET(no).RVEpiInterpX;
DATA.Undo(no,DATA.UndoN(no)).RVEpiInterpYBackup = SET(no).RVEpiInterpY;
DATA.Undo(no,DATA.UndoN(no)).EndoDragedBackup = SET(no).EndoDraged;
DATA.Undo(no,DATA.UndoN(no)).EpiDragedBackup = SET(no).EpiDraged;
DATA.Undo(no,DATA.UndoN(no)).RoiCurrent = SET(no).RoiCurrent;
DATA.Undo(no,DATA.UndoN(no)).RoiN = SET(no).RoiN;
DATA.Undo(no,DATA.UndoN(no)).Roi = SET(no).Roi;
DATA.Undo(no,DATA.UndoN(no)).CenterX = SET(no).CenterX;
DATA.Undo(no,DATA.UndoN(no)).CenterY = SET(no).CenterY;
DATA.Undo(no,DATA.UndoN(no)).Point = SET(no).Point;
DATA.Undo(no,DATA.UndoN(no)).Measure = SET(no).Measure;
% LA/RA
undopos = DATA.UndoN(no);
raobj = SET(no).RA;
if ~isempty(raobj)
  % copy function from generalpen.clatriumpen to copy corrdinates to a new
  % instance of the object
  DATA.Undo(no,undopos).RA = copy(raobj);
else
  DATA.Undo(no,undopos).RA = raobj;
end

laobj = SET(no).LA;
if ~isempty(laobj)
  % copy function from generalpen.clatriumpen to copy corrdinates to a new
  % instance of the object
  DATA.Undo(no,undopos).LA = copy(laobj);
else
  DATA.Undo(no,undopos).LA = laobj;
end

% TO DO General pen

if ~isempty(SET(NO).Scar)
  DATA.Undo(no,DATA.UndoN(no)).ScarResult  = find(SET(NO).Scar.Result);
  DATA.Undo(no,DATA.UndoN(no)).ScarManualP = find(SET(NO).Scar.Manual==1);
  DATA.Undo(no,DATA.UndoN(no)).ScarManualN = find(SET(NO).Scar.Manual==-1);
  DATA.Undo(no,DATA.UndoN(no)).ScarNoReflow = find(SET(NO).Scar.NoReflow);
else
  DATA.Undo(no,DATA.UndoN(no)).ScarResult  = [];
  DATA.Undo(no,DATA.UndoN(no)).ScarManualP = [];
  DATA.Undo(no,DATA.UndoN(no)).ScarManualN = [];
  DATA.Undo(no,DATA.UndoN(no)).ScarNoReflow = [];
end

if ~isempty(SET(NO).MaR)
  DATA.Undo(no,DATA.UndoN(no)).MaRResult  = find(SET(NO).MaR.Result);
  DATA.Undo(no,DATA.UndoN(no)).MaRManualP = find(SET(NO).MaR.Manual==1);
  DATA.Undo(no,DATA.UndoN(no)).MaRManualN = find(SET(NO).MaR.Manual==-1);
else
  DATA.Undo(no,DATA.UndoN(no)).MaRResult  = [];
  DATA.Undo(no,DATA.UndoN(no)).MaRManualP = [];
  DATA.Undo(no,DATA.UndoN(no)).MaRManualN = [];
end

% OBS! For flow these are still SET specific. Need a better undo solution
% for contrast and zoom.
%DATA.Undo(NO,DATA.UndoN(NO)).IntensityMapping = SET(NO).IntensityMapping;
%DATA.Undo(NO,DATA.UndoN(NO)).IntensityOffset = SET(NO).IntensityOffset;
%DATA.Undo(NO,DATA.UndoN(NO)).IntensityScaling = SET(NO).IntensityScaling;
%DATA.Undo(NO,DATA.UndoN(NO)).NormaZoomState = SET(NO).NormalZoomState;
%DATA.Undo(NO,DATA.UndoN(NO)).MontageZoomState = SET(NO).MontageZoomState;
%DATA.Undo(NO,DATA.UndoN(NO)).MontageRowZoomState = SET(NO).MontageRowZoomState;

%-----------------------
function disableundo(no)
%-----------------------
%Disable undo function. Also copies segmentation data
%to make sure data is consistent
global DATA NO

if nargin==0
  no = NO;
end

if ~isempty(DATA.Undo)
  DATA.UndoN(no) = 0;
else
  % Implies that a non-undoable (irreversible) fcn was called before any irreversible
  % functions have been called, and then the NeedToSave bit should be
  % flipped, so we know there are unsaved things that have happened.
  DATA.NeedToSave=true;
end

%--------------------------
function seted_Callback
%--------------------------
%Set diastole to be current timeframe
global SET NO DATA

%Loop over all linked images to update EDT
nos = SET(NO).Linked;
for no = nos
  SET(no).EDT = SET(NO).CurrentTimeFrame;
end

DATA.updatetimebaraxes;
segment('updatevolume');
segment('updatemeasurement');

%--------------------------
function setes_Callback
%--------------------------
%Set systole to be current timeframe
global SET NO DATA

%Loop over all linked images to update EDT
nos = SET(NO).Linked;
for no = nos
  SET(no).EST = SET(NO).CurrentTimeFrame;
end

DATA.updatetimebaraxes;
segment('updatevolume');
segment('updatemeasurement');

%-------------------------------------
function autoesed_Callback(silent,no)
%-------------------------------------
%Autodetect and store ED, & ES.
global SET NO DATA

if nargin==0
  silent = false;
end

if nargin<2
  no = NO;
end

calcfunctions('calcvolume',no);

%Store
%LV
%set EDT
[SET(no).EDV,SET(no).EDT] = max(SET(no).LVV);
%set EST
[SET(no).ESV,SET(no).EST] = min(SET(no).LVV);
if isequal(SET(no).EDT,SET(no).EST)
  SET(no).ESV=0;
  SET(no).EDV=0;
end
if isnan(SET(no).EDV)
  SET(no).EDV=0;
  SET(no).ESV=0;
  SET(no).EDT=1;
  SET(no).EST=1;
end

SET(no).SV = SET(no).EDV-SET(no).ESV; %Stroke volume
if SET(no).EDV>0
  SET(no).EF = SET(no).SV/SET(no).EDV; %Ejection fraction
else
  SET(no).EF = 0;
end

%RV
if SET(no).EDT==1 && SET(no).EST==1
  [SET(no).RVEDV,SET(no).EDT] = max(SET(no).RVV);
  [SET(no).RVESV,SET(no).EST] = min(SET(no).RVV);
  if isequal(SET(no).EDT,SET(no).EST)
    SET(no).RVESV=0;
    SET(no).RVEDV=0;
  end
  if isnan(SET(no).RVEDV)
    SET(no).RVEDV=0;
    SET(no).RVESV=0;
    SET(no).EDT=1;
    SET(no).EST=1;
  end
end

SET(no).RVSV = SET(no).RVEDV-SET(no).RVESV; %Stroke volume
if SET(no).RVEDV>0
  SET(no).RVEF = SET(no).RVSV/SET(no).RVEDV; %Ejection fraction
else
  SET(no).RVEF = 0;
end

%also check for imageviewplane
imvpcell={'Short-axis','Transversal','2CH','3CH','4CH','Flow'};
imtcell={'Flow (magnitude)','Flow (through plane)'};
include=zeros(1,length(SET));
for i = 1 : length(SET)
  include(i) = ismember(SET(i).ImageViewPlane,imvpcell) || ismember(SET(i).ImageType,imtcell);
end

nos=find(include);
%apply estimated
for tmpno=nos
  [~,SET(tmpno).EST] = min(abs(SET(no).EST/SET(no).TSize-(1:SET(tmpno).TSize)/SET(tmpno).TSize));
  [~,SET(tmpno).EDT] = min(abs(SET(no).EDT/SET(no).TSize-(1:SET(tmpno).TSize)/SET(tmpno).TSize));
end

if not(silent)
  segment('updatevolume');
end
DATA.updatetimebaraxes;

%---------------------------
function hideid_Callback
%---------------------------
global DATA
% check status of anonymous mode

DATA.HideID = DATA.Handles.configiconholder.findindented('hideid');

DATA.updatetitle;

%----------------------------------------
function viewpatientinfo_Callback(arg,no)
%----------------------------------------
%GUI to view and adjust patient information.

global DATA SET NO

if ~isa(DATA,'maingui')
  delete(gcbf)
  return
end

if not(DATA.DataLoaded)
  return;
end

gui = DATA.GUI.PatientInfo;

if nargin < 2
  no = NO;
  if nargin < 1
    arg = 'init';
  end
end

switch arg
  case 'init'
    gui = mygui('patientinfo.fig');
    myadjust(gui.fig)
    set(gui.handles.bsaunittxt,'String', ['m', char(178)])
    set(gui.handles.bmiunittxt,'String', ['kg/m', char(178)])
    DATA.GUI.PatientInfo = gui;
    gui.no = no;
   
    viewpatientinfo_Callback('update');

  case 'update'
    %Calculate bsa
    SET(gui.no).PatientInfo.BSA = calcfunctions('calcbsa',...
      SET(gui.no).PatientInfo.Weight,...
      SET(gui.no).PatientInfo.Length);

    %Calculate bmi
    SET(gui.no).PatientInfo.BMI = calcfunctions('calcbmi',...
      SET(gui.no).PatientInfo.Weight,...
      SET(gui.no).PatientInfo.Length);

    %Populate edit boxes
    set(gui.handles.nameedit,'String',SET(gui.no).PatientInfo.Name);
    set(gui.handles.idedit,'String',SET(gui.no).PatientInfo.ID);
    set(gui.handles.birthdateedit,'String',SET(gui.no).PatientInfo.BirthDate);
    set(gui.handles.acquisitiondateedit,'String',SET(gui.no).PatientInfo.AcquisitionDate);
    [agedigit,ageunit] = calcfunctions('calcagewithunits',SET(gui.no).PatientInfo.Age);
    set(gui.handles.patientagetext,'String',ageunit);
    set(gui.handles.ageedit,'String',agedigit);
    set(gui.handles.lengthedit,'String',SET(gui.no).PatientInfo.Length);
    set(gui.handles.weightedit,'String',SET(gui.no).PatientInfo.Weight);
    set(gui.handles.bsaedit,'String',SET(gui.no).PatientInfo.BSA);
    set(gui.handles.bmiedit,'String',SET(gui.no).PatientInfo.BMI);
    if isempty(SET(gui.no).PatientInfo.Sex)
      SET(gui.no).PatientInfo.Sex = '-';
    end
    switch lower(SET(gui.no).PatientInfo.Sex(1))
      case 'm'
        set(gui.handles.sexlistbox,'value',1);
        SET(gui.no).PatientInfo.Sex = 'M';
      case 'f'
        set(gui.handles.sexlistbox,'value',2);
        SET(gui.no).PatientInfo.Sex = 'F';
      otherwise
        set(gui.handles.sexlistbox,'value',3);
        SET(gui.no).PatientInfo.Sex = '-';
    end
    
    %Call to update title
    segment('updatetitle');
    
  case 'nameedit'
    s = mygetedit(gui.handles.nameedit);
    SET(gui.no).PatientInfo.Name = s;
    viewpatientinfo_Callback('update');
  case 'idedit'
    s = mygetedit(gui.handles.idedit);
    SET(gui.no).PatientInfo.ID = s;
    viewpatientinfo_Callback('update');
  case 'birthdateedit'
    s = mygetedit(gui.handles.birthdateedit);
    SET(gui.no).PatientInfo.BirthDate = s;
    %check birthdate, age and acquisition date consistency
    if not(isempty(SET(gui.no).PatientInfo.AcquisitionDate)) && not(isempty(SET(gui.no).PatientInfo.Age))
      [agedigitnew,ageunitnew] = calcfunctions('calcagefrombirthdate',SET(gui.no).PatientInfo.BirthDate, SET(gui.no).PatientInfo.AcquisitionDate);
      [agedigitold,ageunitold] = calcfunctions('calcagewithunits',SET(gui.no).PatientInfo.Age);
      if not(isequal(ageunitnew,ageunitold(1))) && isequal(agedigitnew,agedigitold)
        % something went wrong in the setting
        mywarning(dprintf('Birth date and age are not consitent.\nSetting age to the value corresponding Acqusition Date - Birth Date.'));
        SET(gui.no).PatientInfo.Age = calcfunctions('calcageyearfraction',agedigitnew,ageunitnew);
      end
    end
    viewpatientinfo_Callback('update');
  case 'ageedit'
    s = mygetedit(gui.handles.ageedit);
    n = str2double(s);
    if not(isnan(n))
      SET(gui.no).PatientInfo.Age = n;
      %check birthdate, age and acquisition date consistency
      if not(isempty(SET(gui.no).PatientInfo.AcquisitionDate)) && not(isempty(SET(gui.no).PatientInfo.BirthDate))
        [agedigitnew,ageunitnew] = calcfunctions('calcagefrombirthdate',SET(gui.no).PatientInfo.BirthDate, SET(gui.no).PatientInfo.AcquisitionDate);
        [agedigitold,ageunitold] = calcfunctions('calcagewithunits',SET(gui.no).PatientInfo.Age);
        if not((ageunitnew == ageunitold(1)) && (agedigitnew == agedigitold))
          % something went wrong in the setting
          mywarning(dprintf('Birth date and age are not consitent.\nSetting age to the value corresponding Acqusition Date - Birth Date.'));
          SET(gui.no).PatientInfo.Age = calcfunctions('calcageyearfraction',agedigitnew,ageunitnew);
        end
      end
    end
    viewpatientinfo_Callback('update');
  case 'weightedit'
    s = mygetedit(gui.handles.weightedit);
    n = str2double(s);
    if not(isnan(n))
      SET(gui.no).PatientInfo.Weight = n;
    end
    viewpatientinfo_Callback('update');
  case 'lengthedit'
    s = mygetedit(gui.handles.lengthedit);
    n = str2double(s);
    if not(isnan(n))
      SET(gui.no).PatientInfo.Length = n;
    end
    viewpatientinfo_Callback('update');
  case 'sexlistbox'
    temp = 'MF-';
    SET(gui.no).PatientInfo.Sex = temp(mygetlistbox(gui.handles.sexlistbox));
    viewpatientinfo_Callback('update');
  case 'apply'
    temp = SET(gui.no).PatientInfo;
    for no=1:length(SET)
      SET(no).PatientInfo = temp;
    end
    mymsgbox('Information applied to all image stacks.','Done!',DATA.GUI.Segment);
  case 'close'
    %mymsgbox('Save dataset in order to keep changes to patient info.');
    try
      if ~isempty(DATA.GUI.ReportSheetGenerator)
        if strcmp(get(DATA.GUI.ReportSheetGenerator.fig,'Name'),'Segment CT Report')
          reporter.reportsheet('initreferencelistbox','CT');
        else
          reporter.reportsheet('setfiguretitle');
          reporter.reportsheet('initpreferenceselection');
        end
      end
      DATA.GUI.PatientInfo = close(DATA.GUI.PatientInfo);
    catch %#ok<CTCH>
      delete(gui.fig);
    end
  otherwise
    myfailed(dprintf('Unknown argument %s to patientinfo',arg),DATA.GUI.Segment);
    return;
end

%-----------------------------
function [a,b] = flipvars(b,a)
%-----------------------------
%Flip variables. Simple and elegant


%--------------------------------------
function x = translatecontours_helper(x)
%--------------------------------------
%Fix out of bounds indices

x(x==0) = 1;
x(x==length(x)+1) = length(x);

%--------------------------------
function translatecontours(dx,dy)
%--------------------------------
%Translate all contours. Called by hotkeys.

%Einar Heiberg

%Call the helper function
translatecontour_helper(dx,dy,false); %false = do not translate image

%---------------------------------------
function translatecontoursandimage(dx,dy)
%---------------------------------------
%Translate all contours and the image. Called by hotkeys.

%Einar Heiberg

%Call the helper function
translatecontour_helper(dx,dy,true);

%-----------------------------------------------------
function translatecontour_helper(dx,dy,translateimage)
%-----------------------------------------------------
%Helper function to translate contours and image.
%If translateimage is true then image is also translated.

%Not yet supported:
%- Scar
%- MaR
%- Pins
%- Levelset

%Einar Heiberg

global DATA SET NO

no = NO;
slice = SET(no).StartSlice:SET(no).EndSlice;
if isempty(slice)
  slice = SET(no).CurrentSlice;
end

%Endo
if ~isempty(SET(no).EndoX)
  SET(no).EndoX(:,:,slice) = SET(no).EndoX(:,:,slice)+dx;
  SET(no).EndoY(:,:,slice) = SET(no).EndoY(:,:,slice)+dy;
end

%Epi
if ~isempty(SET(no).EpiX)
  SET(no).EpiX(:,:,slice) = SET(no).EpiX(:,:,slice)+dx;
  SET(no).EpiY(:,:,slice) = SET(no).EpiY(:,:,slice)+dy;
end

%RVEndo
if ~isempty(SET(no).RVEndoX)
  SET(no).RVEndoX(:,:,slice) = SET(no).RVEndoX(:,:,slice)+dx;
  SET(no).RVEndoY(:,:,slice) = SET(no).RVEndoY(:,:,slice)+dy;
end

%RVEpi
if ~isempty(SET(no).RVEpiX)
  SET(no).RVEpiX(:,:,slice) = SET(no).RVEpiX(:,:,slice)+dx;
  SET(no).RVEpiY(:,:,slice) = SET(no).RVEpiY(:,:,slice)+dy;
end

if (~isempty(SET(no).EndoInterpX))
  for tloop = 1:SET(no).TSize
    for zloop = slice
      %Endo-Interp
      if ~isempty(SET(no).EndoInterpX{tloop,zloop})
        SET(no).EndoInterpX{tloop,zloop} = SET(no).EndoInterpX{tloop,zloop}+dx;
        SET(no).EndoInterpY{tloop,zloop} = SET(no).EndoInterpY{tloop,zloop}+dy;
      end
    end
  end
end
if (~isempty(SET(no).EpiInterpX))
  for tloop = 1:SET(no).TSize
    for zloop = slice
      %Endo-Interp
      if ~isempty(SET(no).EpiInterpX{tloop,zloop})
        SET(no).EpiInterpX{tloop,zloop} = SET(no).EpiInterpX{tloop,zloop}+dx;
        SET(no).EpiInterpY{tloop,zloop} = SET(no).EpiInterpY{tloop,zloop}+dy;
      end
    end
  end
end
if (~isempty(SET(no).RVEndoInterpX))
  for tloop = 1:SET(no).TSize
    for zloop = slice
      %RVEndo-Interp
      if ~isempty(SET(no).RVEndoInterpX{tloop,zloop})
        SET(no).RVEndoInterpX{tloop,zloop} = SET(no).RVEndoInterpX{tloop,zloop}+dx;
        SET(no).RVEndoInterpY{tloop,zloop} = SET(no).RVEndoInterpY{tloop,zloop}+dy;
      end
    end
  end
end
if (~isempty(SET(no).RVEpiInterpX))
  for tloop = 1:SET(no).TSize
    for zloop = slice
      %RVEndo-Interp
      if ~isempty(SET(no).RVEpiInterpX{tloop,zloop})
        SET(no).RVEpiInterpX{tloop,zloop} = SET(no).RVEpiInterpX{tloop,zloop}+dx;
        SET(no).RVEpiInterpY{tloop,zloop} = SET(no).RVEpiInterpY{tloop,zloop}+dy;
      end
    end
  end
end

%ROI
if ~isempty(SET(no).Roi) && ~isempty(SET(no).Roi.X)
  for loop=1:length(SET(no).Roi)
    if (SET(no).Roi(loop).Z>=slice(1)) && (SET(no).Roi(loop).Z<=slice(end))
      SET(no).Roi(loop).X = SET(no).Roi(loop).X+dx;
      SET(no).Roi(loop).Y = SET(no).Roi(loop).Y+dy;
    end
  end
end

%Measurements
if ~isempty(SET(no).Measure)
  for loop = 1:length(SET(no).Measure)
    if (SET(no).Measure(loop).Z(1)>=slice(1)) && (SET(no).Measure(loop).Z(2)<=slice(end))
      SET(no).Measure(loop).X = SET(no).Measure(loop).X+dx;
      SET(no).Measure(loop).Y = SET(no).Measure(loop).Y+dy;
    end
  end
end

%Point
if ~isempty(SET(no).Point)
  for loop = 1:length(SET(no).Point.X)
    if (SET(no).Point.Z(loop)>=slice(1)) && (SET(no).Point.Z(loop)<=slice(end))
      SET(no).Point.X(loop) = SET(no).Point.X(loop)+dx;
      SET(no).Point.Y(loop) = SET(no).Point.Y(loop)+dy;
    end
  end
end

if translateimage
  %--- Image
  x = (1:SET(no).XSize)-dx;
  y = (1:SET(no).YSize)-dy;

  %Get rid of out of bands
  x = translatecontours_helper(x);
  y = translatecontours_helper(y);

  SET(no).IM(:,:,:,slice) = SET(no).IM(x,y,:,slice);
end

DATA.ViewIM{DATA.CurrentPanel} = [];

%Update screen
drawfunctions('drawpanel',DATA.CurrentPanel);

viewpanels = DATA.ViewPanels;
for loop = 1:length(viewpanels)
  if viewpanels(loop)>0
    drawfunctions('drawintersections',loop);
  end
end

%-------------------------------------
function setcolormap_Callback(type,no)
%-------------------------------------
%Set colormap for current image stack.
global DATA SET NO

if nargin < 2
  no = NO;
end

n = DATA.GUISettings.ColorMapSize;
set([DATA.Handles.colormapgraymenu ...
  DATA.Handles.hsvcolormapmenu ...
  DATA.Handles.jetcolormapmenu ...
  DATA.Handles.hotcolormapmenu ...
  DATA.Handles.spectcolormapmenu ...
  DATA.Handles.originalcolormapmenu],'checked','off');

switch type
  case 'original'
    %Apply original colormap (retrieved from DICOM upon loading) if exists
    if isempty(SET(no).ColormapOriginal)
      myfailed('Unknown colormap',DATA.GUI.Segment);
      return
    end
    tempmap = SET(no).ColormapOriginal;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
  case 'gadgetron'
    %Apply Gadgetron colormap
    load('colormapgadgetron');
    tempmap = colormapgadgetron;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
  case 'siemens'
    %Apply Gadgetron colormap
    load('colormapsiemens');
    tempmap = colormapsiemens;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
  case 'gray'
    SET(no).Colormap = []; %gray(n);
    %set(DATA.Handles.colormapgraymenu,'checked','on');
  case 'hsv'
    SET(no).Colormap = hsv(n);
    %set(DATA.Handles.hsvcolormapmenu,'checked','on');
  case 'jet'
    SET(no).Colormap = jet(n);
    %set(DATA.Handles.jetcolormapmenu,'checked','on');
  case 'hot'
    SET(no).Colormap = hot(n);
    %set(DATA.Handles.hotcolormapmenu,'checked','on');
  case 'spect'
    SET(no).Colormap = spect(n);
    %set(DATA.Handles.spectcolormapmenu,'checked','on');
  case 't1'
    %Try to figure out one to use and then call the
    %correct one "recursively"

    stri = t1checktypehelper(no);
    setcolormap_Callback(stri,no);
    return; %we should be done

  case 't1pre'
    %set colormap for T1 images
    load('colormapt1');
    tempmap = flipud(colormapt1)./256;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
    %set contrast and brightness
    window =3000;  %1400, widht
    level =1500; %1300, center
    [contrast, brightness] = calcfunctions('win2con',window,level);
    SET(no).IntensityMapping.Contrast = contrast;
    SET(no).IntensityMapping.Brightness = brightness;
    %update image
    %     drawfunctions('drawcontrastimage',no);  %done by drawall below
  case 't1post'
    %set colormap for T1 images
    load('colormapt1');
    tempmap = flipud(colormapt1)./256;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
    %set contrast and brightness
    window = 600; %widht
    level = 500; %center
    [contrast, brightness] = calcfunctions('win2con',window,level);
    SET(no).IntensityMapping.Contrast = contrast;
    SET(no).IntensityMapping.Brightness = brightness;
    %update image
  case 't2'
    %set colormap for T2 images
    load('colormapt1');
    tempmap = flipud(colormapt1)./256;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
    %set contrast and brightness
    window = 120; %widht
    level = 60; %center
    [contrast, brightness] = calcfunctions('win2con',window,level);
    SET(no).IntensityMapping.Contrast = contrast;
    SET(no).IntensityMapping.Brightness = brightness;
    %update image
    %     drawfunctions('drawcontrastimage',no);  %done by drawall below
  case 't2star'
    %set colormap for T2*
    load('colormapt1');
    tempmap = flipud(colormapt1)./256;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
    %set contrast and brightness
    window = 80; %widht
    level = 45; %center
    [contrast, brightness] = calcfunctions('win2con',window,level);
    SET(no).IntensityMapping.Contrast = contrast;
    SET(no).IntensityMapping.Brightness = brightness;
    %update image
    %     drawfunctions('drawcontrastimage',no); %done by drawall below
  case 'ecv'
    %set colormap for ECV
    load('colormapecv');
    tempmap = flipud(colormapecv)./256;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(tempmap),tempmap(:,1),linspace(1,length(tempmap),n))';
    cmap(:,2) = interp1(1:length(tempmap),tempmap(:,2),linspace(1,length(tempmap),n))';
    cmap(:,3) = interp1(1:length(tempmap),tempmap(:,3),linspace(1,length(tempmap),n))';
    SET(no).Colormap = cmap;
    %set contrast and brightness
    window = 100; %widht
    level = 50; %center
    [contrast, brightness] = calcfunctions('win2con',window,level);
    SET(no).IntensityMapping.Contrast = contrast;
    SET(no).IntensityMapping.Brightness = brightness;
    %update image
    %     drawfunctions('drawcontrastimage',no);  % done by drawall below
  otherwise
    myfailed('Unknown colormap',DATA.GUI.Segment);
    return;
end

panelstodo = find(DATA.ViewPanels==no);

for panel=panelstodo
  if DATA.GUISettings.ShowColorbar
    createfunctions('addcolorbar',panel)
  end
  DATA.ViewIM{panel} = [];
end

drawfunctions('drawno',no)

recalc=1;
updatedslider=0;
drawfunctions('drawthumbnails',recalc,updatedslider);
% drawfunctions('drawcontrastimage',no);  %drawfunctions('drawall',DATA.ViewMatrix);

%------------------------------------
function stri = t1checktypehelper(no)
%------------------------------------
%Check the type of image, if it is either 't1pre' or 't1post'

global SET 

stri = 't1pre'; %Guess for now

%Compute intensity in bloodpool
si = calcfunctions('calcbloodpool',no);

if false %~isnan(si)
  si = calcfunctions('calctruedata',si,no);
  if si>1200
    stri = 't1pre';
  else
    stri = 't1post';
  end
  return
end

%Try without endocardium
im = SET(no).IM(:,:,SET(no).CurrentTimeFrame,max(1,round(SET(no).ZSize/2)));
im = calcfunctions('calctruedata',im,no);
threshold = 1200;
percentbloodpool = sum(im(:)> threshold)/numel(im);

if percentbloodpool>0.08 %in tested is seems like it should be around 12%, ones that are not are in the order of 0.6%
  stri = 't1pre';
else
  stri = 't1post';
end

%-------------------------------
function slidingaverage_Callback
%-------------------------------
%Performs sliding average of an image

global SET NO;

%Define filter
n = 10;
f = single(repmat(1/n,[1 1 n]));
h = waitbar(0,'Please wait.');
for zloop = 1:SET(NO).ZSize
  SET(NO).IM(:,:,:,zloop) = econv3(SET(NO).IM(:,:,:,zloop),f);
  waitbar(zloop/SET(NO).ZSize,h);
end
close(h);

%-----------------------------------------
function imageenhancement_Callback
%-----------------------------------------
%Image enhancement using adapthisteq (CLAE), please see adapthisteq for details.

%Einar Heiberg

global SET NO DATA
str1 = dprintf('Image enhancement is irreversible.');
str2 = dprintf('Are you sure?');
str = sprintf('%s %s',str1,str2);
if not(yesno(str,[],DATA.GUI.Segment))
  return;
end

totalnumsteps =SET(NO).TSize*SET(NO).ZSize ;

no = NO;

numsteps = 0;
h = waitbar(0,'Please wait.');
for tloop = 1:SET(NO).TSize
  for zloop = 1:SET(NO).ZSize
    numsteps = numsteps+1;
    waitbar(numsteps/totalnumsteps,h);
    SET(no).IM(:,:,tloop,zloop) = adapthisteq(SET(no).IM(:,:,tloop,zloop));
  end
end
close(h);
DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimage',DATA.CurrentPanel);

%--------------------------------
function unlinkimages_Callback
%--------------------------------
%Unlink images
global SET NO
mywarning('There is no built-in undo function for this. Are you sure?');

no = NO;
if ~isempty(SET(NO).Parent)
  no = SET(NO).Parent;
end

linked = SET(no).Linked;
for loop = 1:length(linked)
  SET(linked(loop)).Parent = [];
  SET(linked(loop)).Children = []; % Just to be sure
  SET(linked(loop)).Linked = linked(loop);
end
SET(no).Parent = [];
SET(no).Children = [];
SET(no).Linked = no;

viewfunctions('setview');

%------------------------------------------------------------
function orthogonalresamplehelper(no,P,flipv,pos,newdimorder)
%------------------------------------------------------------
%Helper function to transversal2sagittal (and more in the future), this
%function should be genereric for all coordinate transformation if I have
%not been thinking completly crazy.

%Einar Heiberg
global SET

%Calculate 3D coordinates, used to calculate image position and image
%orientation.

[pos1] = calcfunctions('xyz2rlapfh',no,0,0,0); %"origo"
[pos2] = calcfunctions('xyz2rlapfh',no,P(1,1),P(1,2),P(1,3));
[pos3] = calcfunctions('xyz2rlapfh',no,P(2,1),P(2,2),P(2,3));
[pos4] = calcfunctions('xyz2rlapfh',no,P(3,1),P(3,2),P(3,3));

%Calculate directions and resolution
dir1 = pos2-pos1;
step1 = norm(dir1,2);
dir1 = dir1./step1;

dir2 = pos3-pos1;
step2 = norm(dir2,2);
dir2 = dir2./step2;

dir3 = pos4-pos1;
step3 = norm(dir3,2);

%Orientation and imageposition
orientation = [dir1 dir2];
imageposition = calcfunctions('xyz2rlapfh',no,pos(1),pos(2),pos(3)); %"origo"

SET(no).ImageOrientation = orientation;
SET(no).ImagePosition = imageposition;
SET(no).ResolutionX = step2;
SET(no).ResolutionY = step1;
SET(no).SliceThickness = step3;
SET(no).SliceGap = 0;
SET(no).IM = permute(SET(no).IM,newdimorder);
SET(no).XSize = size(SET(no).IM,1);
SET(no).YSize = size(SET(no).IM,2);
SET(no).TSize = size(SET(no).IM,3);
SET(no).ZSize = size(SET(no).IM,4);
SET(no).NormalZoomState = [];
SET(no).MontageZoomState = [];

%Flip if necessary
for loop = 1:4
  if flipv(loop)
    SET(no).IM = flip(SET(no).IM,loop);
  end
end

%Clear segmentation and other stuff
callbackfunctions('segmentclearall_Callback',true);
SET(no).Scar = []; %Clear scar
SET(no).LevelSet = []; %Clear levelset
% callbackfunctions('pointclearall_helper',no); %Clear points
callbackfunctions('pointclearall_Callback');
SET(no).PapillaryIM = [];
SET(no).Measure = [];
roi('roireset',no); %reset ROI's
currentslice = round(SET(no).ZSize/2);
SET(no).StartSlice = currentslice;
SET(no).EndSlice = currentslice;
SET(no).CurrentSlice = currentslice;
SET(no).Strain = [];
SET(no).StrainTagging = [];
SET(no).Stress = [];
SET(no).MaR = [];
SET(no).T2 = [];
SET(no).ECV = [];

%Calculate the preview
sliderupdated = false;
calculatepreview = true;
viewfunctions('switchimagestack',no);
viewfunctions('switchslice',0);
drawfunctions('drawthumbnails',calculatepreview,sliderupdated);

%------------------------------------------------------
function ok = orthogonalresamplecheck(no,dir1,dir2,str)
%------------------------------------------------------
%Check if can create orthogonal stack, give error message if not.

%Einar Heiberg
global SET

ok = false;

if ~yesno('This will clear all segmentation and measurements. Are you sure?')
  return;
end

%Check it is correct direction
if (acos(abs(dir1.*SET(no).ImageOrientation(1:3)))*180/pi)>10
  myfailed(sprintf('Current image stack is not a %s stack.',str));
  return;
end

if (acos((dir2.*SET(no).ImageOrientation(4:6)))*180/pi)>10
  myfailed(sprintf('Current image stack is not a %s stack.',str));
  return;
end

if ~isempty(SET(no).Flow)
  myfailed('Can not resample flow image stack.');
  return;
end

if (SET(no).XSize<10) || (SET(no).YSize<10) || (SET(no).ZSize<10)
  myfailed('Can resample slices with slices.');
  return;
end

%So far good!
ok = true;

%------------------------------------
function transversal2coronal_Callback
%------------------------------------
%Takes a transversal image stack and converts it to a coronal image stack

%Einar Heiberg
global SET NO

%1 0 0   0 1 0
if ~orthogonalresamplecheck(NO,[1 0 0],[0 1 0],'transversal')
  return;
end

%Define three points in the image; (x,y,z)
P = [...
  0 1 0 ; ...
  0 0 1 ; ...
  1 0 0];

flipv = [0 0 0 1]; %flip z

pos = [size(SET(NO).IM,1) 1 1];

%xs,ys,ts,zs => zs,ys,ts,xs
newdimorder = [4 2 3 1];

orthogonalresamplehelper(NO,P,flipv,pos,newdimorder);

SET(NO).ImageViewPlane = 'Coronal';

%-------------------------------------
function transversal2sagittal_Callback
%-------------------------------------
%Takes a transversal image stack and converts it to a sagittal image stack

%Einar Heiberg

global SET NO


if ~orthogonalresamplecheck(NO,[1 0 0],[0 1 0],'transversal')
  return;
end

%Define three points in the image; (x,y,z)
P = [...
  1 0 0 ; ... %dir1
  0 0 1 ; ... %dir2
  0 1 0];     %dir3

flipv = [0 0 0 0]; %do not flip any direction

pos = [0.5 0.5 0.5];

newdimorder = [4 1 3 2];

orthogonalresamplehelper(NO,P,flipv,pos,newdimorder);

SET(NO).ImageViewPlane = 'Sagittal';

%------------------------------
function intensitycrop_Callback
%------------------------------
%Crop intensities above and under HU threshold

global DATA SET NO

no = NO;
if ~isequal(SET(no).Modality,'CT')
  myfailed('Only CT are supported.');
  return;
end

%Ask for values
s = [];
s.MinHU = -1200;
s.MaxHU = 1200;
s = inputstruct(s,'Set HU limit');

if isa(SET(no).IM,'int16')
  SET(no).IM = max(int16(s.MinHU),SET(no).IM);
  SET(no).IM = min(int16(s.MaxHU),SET(no).IM);
  SET(no).minValue = single(min(SET(no).IM(:))); %s.MinHU;
  SET(no).maxValue = single(max(SET(no).IM(:))); %s.MaxHU;
else
  %Stored as normalized values

  %truedata is calculated as:
  %z = im*SET(no).IntensityScaling+SET(no).IntensityOffset;

  %maxim*scaling+offset = maxHU
  %minim*scaling+offset = minHU
  minim = (s.MinHU-SET(no).IntensityOffset)/SET(no).IntensityScaling;
  maxim = (s.MaxHU-SET(no).IntensityOffset)/SET(no).IntensityScaling;
  SET(no).IM = max(SET(no).IM,minim);
  SET(no).IM = min(SET(no).IM,maxim);

  %New normalization
  SET(no).IM = (SET(no).IM-minim)/(maxim-minim);
  SET(no).IntensityScaling = SET(no).IntensityScaling/(maxim-minim);
  SET(no).IntensityScaling = SET(no).IntensityScaling/(maxim-minim);
  SET(no).IntensityOffset = s.MinHU;
end

%Update graphically
sliderupdated = false;
calculatepreview = true;
drawfunctions('drawthumbnails',calculatepreview,sliderupdated);
DATA.ViewIM{DATA.CurrentPanel} = [];
drawfunctions('drawimages',DATA.CurrentPanel);

%----------------
function smartepi
%----------------
%Tool for smart epi segmentation from endo ballon tool
%activates by 'space' when you have current tool EndoBalloon

global SET NO

segmentation('epicopyfromlvendo_Callback'); %copy segmentation from endo to epi and expand 5mm
if ~isempty(SET(NO).StrainTagging) && isfield(SET(NO).StrainTagging, 'LVupdated')
  SET(NO).StrainTagging.LVupdated = 1;
end
if ~isempty(SET(NO).StrainMitt)
  updateredo(SET(NO).StrainMitt,'Epi');
end
calcfunctions('updatemarandscar',no);

%----------------
function closeallongoinginterpolations  %used in filemenu.m
%----------------
%go over all stacks and close all ongoing interpolations
global SET
silentmode = true;
for noloop = 1:length(SET)
  connectinterpolation(noloop,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'},silentmode)
end
%----------------
function numpoints = getnumpointsforno(no)
%----------------
% get nuber of points in the contour depending on the imageview
global SET DATA
if any(strcmpi(SET(no).ImageViewPlane,{'2CH','3CH','4CH'})) && ~(contains(SET(no).ImageType,'Flow','IgnoreCase',true))
  numpoints = DATA.NumPointsLax;
else
  numpoints = DATA.NumPoints;
end

%----------------
function numpoints = getnuminterppointsforno(no,type,slice,objectind)
%----------------
% get nuber of points for the interpolation depending on the imageview
arguments
  no
  type {mustBeMember(type,{'GeneralPenInterp','LAInterp','RAInterp',...
    'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'})}
  slice {mustBeInteger}
  objectind {mustBeInteger} = []; %used for general pen objects
end
global DATA SET

if strcmp(type,'GeneralPenInterp') && isempty(objectind)
  objectind = DATA.GeneralPenSettings.getcurrentobject;
end

isinterpongoing = helperfunctions('isinterpongoing',SET(no),type,objectind);
if isinterpongoing
  tf = SET(no).CurrentTimeFrame;
  typeInterpX = helperfunctions('parsesetfield',SET(no),type,'X',objectind,tf,slice);
  numinterppoints = length(typeInterpX);
  datanumpoints = getnumpointsforno(no);
  numpoints = floor(datanumpoints/15)*(numinterppoints-1);
  if numpoints > datanumpoints
    numpoints = datanumpoints;
  end
else
  if ismember(type,{'GeneralPenInterp','LAInterp','RAInterp'})
    % get number
    numpoints = SET(no).(type(1:end-6)).NumPoints;
  else
    numpoints = getnumpointsforno(no);
  end
end

%---------------------------------
function connectgeneralpeninterpolation(no)
%---------------------------------
% function to close interpolation line
global DATA SET

%Parameters
panels = find(DATA.ViewPanels == no);
timeframe = SET(no).CurrentTimeFrame;
slice = SET(no).CurrentSlice;
objectind = DATA.GeneralPenSettings.getcurrentobject;
type = 'GeneralPenInterp';

%Reset InterpOngoing flag
resetinterpolation(no,{type},objectind);

%Close contour
datax = SET(no).GeneralPenObjects(objectind).InterpX{timeframe,slice};
datay = SET(no).GeneralPenObjects(objectind).InterpY{timeframe,slice};
if length(datax) > 2

  numpoints = SET(no).GeneralPenObjects(objectind).getnumpoints;
  [x,y] = distributepointsalongcontour(datax,datay,numpoints);

  if DATA.ThisFrameOnly
    SET(no).GeneralPenObjects(objectind).X(:,timeframe,slice) = x;
    SET(no).GeneralPenObjects(objectind).Y(:,timeframe,slice) = y;
  else
    % duplicate new contours to all time frames
    SET(no).GeneralPenObjects(objectind).X(:,:,slice) = repmat(x,1,SET(no).TSize,1);
    SET(no).GeneralPenObjects(objectind).Y(:,:,slice) = repmat(y,1,SET(no).TSize,1);
  end

else
  SET(no).GeneralPenObjects(objectind).InterpX{timeframe,slice} = [];
  SET(no).GeneralPenObjects(objectind).InterpY{timeframe,slice} = [];
end

%Graphical update
try
  for numpanel = 1:length(panels)
    resetinterpolationcontours(panels(numpanel),'GeneralPen')
%     drawfunctions('drawcontoursgeneralpen',panels(numpanel))
    drawfunctions('drawno',no);
    if length(datax) < 3
      DATA.Handles.generalpencontour(panels(numpanel)).XData = nan;
      DATA.Handles.generalpencontour(panels(numpanel)).YData = nan;
    end
  end
catch me
  mydispexception(me)
end

%---------------------------------
function connectlarainterpolation(no,type)
%---------------------------------
% function to close interpolation line
global DATA SET

if (isempty(SET(no).(type(1:end-6))) || ...
    isempty(SET(no).(type(1:end-6)).InterpX)) || ...
    ~SET(no).(type(1:end-6)).InterpOngoing
  return
end

%Parameters
panels = find(DATA.ViewPanels == no);
timeframe = SET(no).CurrentTimeFrame;
slice = SET(no).CurrentSlice;

%Reset InterpOngoing flag
resetinterpolation(no,{type});

%Interpolation segmentation is done, save the contour
datax = SET(no).(type(1:end-6)).InterpX{timeframe,slice};
datay = SET(no).(type(1:end-6)).InterpY{timeframe,slice};
if length(datax) > 2

  opencontour = true;
  numpoints = SET(no).(type(1:end-6)).getnumpoints;
  [x,y] = distributepointsalongcontour(datax,datay,numpoints,opencontour);

  if DATA.ThisFrameOnly
    SET(no).(type(1:end-6)).X(:,timeframe,slice) = x;
    SET(no).(type(1:end-6)).Y(:,timeframe,slice) = y;
  else
    % duplicate new contours to all time frames
    if isrow(x)
      x = x'; % make sure column vector
    end
    if isrow(y)
      y = y'; % make sure column vector
    end
    SET(no).(type(1:end-6)).X(:,:,slice) = repmat(x,1,SET(no).TSize,1);
    SET(no).(type(1:end-6)).Y(:,:,slice) = repmat(y,1,SET(no).TSize,1);
  end

else
  SET(no).(type(1:end-6)).InterpX{timeframe,slice} = [];
  SET(no).(type(1:end-6)).InterpY{timeframe,slice} = [];
end

%Graphical update
try
  for numpanel = 1:length(panels)
    resetinterpolationcontours(panels(numpanel),type(1:end-6))
%     drawfunctions('drawcontourslara',panels(numpanel),type(1:end-6))
    drawfunctions('drawno',no);
    if length(datax) < 3
      DATA.Handles.([lower(type(1:end-6)),'contour'])(panels(numpanel)).XData = nan;
      DATA.Handles.([lower(type(1:end-6)),'contour'])(panels(numpanel)).YData = nan;
    end
    if (strcmp(type(1:end-6),'LA') || strcmp(type(1:end-6),'RA')) && ismember(timeframe,[SET(no).EST, SET(no).EDT])
      segment('updatemeasurement');
    end
  end
catch me
  mydispexception(me)
end

%----------------------------------
function connectallinterpolation(no)
%----------------------------------
%Function to connect all ongoing interpolation contour. Used when switching
%time frame or stack, for example.

if nargin < 1
  global DATA %#ok<TLEV> call DATA only is necessary
  panel = DATA.CurrentPanel;
  try
    no = DATA.ViewPanels(panel);
  catch
    global NO %#ok<TLEV> 
    no = NO;
  end
end

connectinterpolation(no,{'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'});
connectlarainterpolation(no,'LAInterp');
connectlarainterpolation(no,'RAInterp');

%---------------------------------
function connectinterpolation(no,typestoconnect,silentmode)
%---------------------------------
% function to close interpolation line

global DATA SET
if nargin < 3
  silentmode = false;
end
interptypes = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
if nargin == 1
  iconnames = {'interpendo','interpepi','interprvendo','interprvepi'};
  numicons = length(iconnames);
  iconstate = zeros(1,numicons);
  for loop = 1:numicons
    iconstate(loop) = DATA.Handles.configiconholder.findindented(iconnames{loop});
  end
  ind = find(iconstate);
else
  interptypes = interptypes(ismember(interptypes,typestoconnect));
  ind = [];
  for typeloop = 1:length(interptypes)
    type = interptypes{typeloop};
    if SET(no).([type(1:end-6),'InterpOngoing'])
      ind = [ind,typeloop]; %#ok<AGROW>
    end
  end
end

panels = find(DATA.ViewPanels == no);

if ~isempty(ind)
  for loop = 1:length(ind)
    type = interptypes{ind(loop)};
    resetinterpolation(no,{type})
    datax = SET(no).([type,'X']){SET(no).CurrentTimeFrame,SET(no).CurrentSlice};
    datay = SET(no).([type,'Y']){SET(no).CurrentTimeFrame,SET(no).CurrentSlice};
    if length(datax) > 2

      numpoints = getnumpointsforno(no);
      [x,y] = distributepointsalongcontour(datax,datay,numpoints);

      if DATA.ThisFrameOnly
        SET(no).([type(1:end-6),'X'])(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)= x;
        SET(no).([type(1:end-6),'Y'])(:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice)= y;
      else
        % duplicate new contours to all time frames
        SET(no).([type(1:end-6),'X'])(:,:,SET(no).CurrentSlice)= repmat(x,1,SET(no).TSize,1);
        SET(no).([type(1:end-6),'Y'])(:,:,SET(no).CurrentSlice)= repmat(y,1,SET(no).TSize,1);
      end

    else
      SET(no).([type,'X']){SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
      SET(no).([type,'Y']){SET(no).CurrentTimeFrame,SET(no).CurrentSlice} = [];
    end
    try
      for numpanel = 1:length(panels)
        resetinterpolationcontours(panels(numpanel),type(1:end-6))
%         drawfunctions('drawcontours',panels(numpanel),type(1:end-6))
        drawfunctions('drawno',no);
        if length(datax) < 3
          DATA.Handles.([lower(type(1:end-6)),'contour'])(panels(numpanel)).XData = nan;
          DATA.Handles.([lower(type(1:end-6)),'contour'])(panels(numpanel)).YData = nan;
        end
      end
    catch me
      mydispexception(me)
    end
  end
end
if ~silentmode
  segment('updatevolume');
end

%---------------------------------
function resetinterpolation(no,interptypes,objectind)
%---------------------------------
%reset interpolation points
global SET

if nargin < 2
  interptypes = {'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp'};
end

if nargin < 3
  objectind = []; % used for general pen objects
end

switch interptypes{1}
  case 'GeneralPenInterp'
    SET(no).GeneralPenObjects(objectind).InterpOngoing = false;
  case {'LAInterp','RAInterp'}
    SET(no).(interptypes{1}(1:end-6)).InterpOngoing = false;
  otherwise
    for loop = 1:length(interptypes)
      type = interptypes{loop};
      SET(no).([type(1:end-6),'InterpOngoing']) = false;
    end
end

%---------------------------------
function resetinterpolationcontours(panel,interptypes)
%---------------------------------
%reset interpolation contours
global DATA
if nargin == 1
  interptypes = {'Endo','Epi','RVEndo','RVEpi'};
else
  interptypes = {interptypes};
end

currentobjectind = 1;
if contains(interptypes,'GeneralPen')
  currentobjectind = DATA.GeneralPenSettings.getcurrentobject;
end

for loop = 1:length(interptypes)
  type = interptypes{loop};
  DATA.Handles.([lower(type),'contour'])(currentobjectind,panel).LineStyle = '-';
  DATA.Handles.([lower(type),'contour'])(currentobjectind,panel).Marker = 'none';
end

%---------------------------------
function cropupdate(nos)  %called from other functions
%---------------------------------
% graphical update after cropping was applied
global DATA
tools('disableundo')

DATA.ViewIM{DATA.CurrentPanel} = [];
for no = nos
  drawfunctions('drawno',no)
end
drawfunctions('drawselectedframe',DATA.CurrentPanel)

%also update thumbnail!
drawfunctions('drawthumbnails',1,0)

%---------------------------------
function [x,y] = distributepointsalongcontour(x,y,numpoints,opencontour,closethecurve)
%---------------------------------
% function to distribute points along whole contour using resampling
arguments
  x
  y
  numpoints = [];
  opencontour = false;
  closethecurve = true;
end

global DATA

if isempty(numpoints)
  numpoints = DATA.NumPoints;
end

%Remove NaN values
valueind = ~isnan(x);
x = x(valueind);
if isrow(x)
  x = x';
end
y = y(valueind);
if isrow(y)
  y = y';
end

%Remove duplicates
nonduplicates = unique([x,y],'rows','stable');
x = nonduplicates(:,1);
y = nonduplicates(:,2);

%Resample the points along the contour
if ~opencontour
  numpointstoresample = numpoints - 1;
else
  numpointstoresample = numpoints;
  closethecurve = false;
end
[x,y] = calcfunctions('resamplecurve',x,y,numpointstoresample,opencontour,closethecurve);

if ~opencontour
  x = [x,x(1)]';
  y = [y,y(1)]';
end