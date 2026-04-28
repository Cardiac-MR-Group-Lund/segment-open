function varargout = findfunctions(varargin)
% FINDFUNCTIONS
% Functions for finding image stacks or slices

% Moved out from segment_main by Nisse Lundahl

%#ok<*GVMIS>
%Invoke subfunction
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%------------------------------------------------------------
function slice = findslicewithlargestaera(heartpart,no,timeframe)
%------------------------------------------------------------
%Find slice with largest delineated endo area
arguments
  heartpart {mustBeMember(heartpart,{'LA','RA','LV','RV'})}
  no
  timeframe = []
end
global SET

if isempty(timeframe)
  timeframe = SET(no).CurentTimeFrame;
end

resy = SET(no).ResolutionY;
resx = SET(no).ResolutionX;

slice = NaN;
maxarea = 0;

for sliceloop = 1:SET(no).ZSize
  %extract X and Y coordinates
  switch heartpart
    case {'LA','RA'}
      if ~isempty(SET(no).(heartpart)) && ~isempty(SET(no).(heartpart).X)
        x = SET(no).(heartpart).X(:,timeframe,sliceloop);
        y = SET(no).(heartpart).Y(:,timeframe,sliceloop);
      else
        continue
      end
    case 'LV'
      if ~isempty(SET(no).EndoX)
        x = SET(no).EndoX(:,timeframe,sliceloop);
        y = SET(no).EndoY(:,timeframe,sliceloop);
      else
        continue
      end
    case 'RV'
      if ~isempty(SET(no).([heartpart,'EndoX']))
        x = SET(no).([heartpart,'EndoX'])(:,timeframe,sliceloop);
        y = SET(no).([heartpart,'EndoY'])(:,timeframe,sliceloop);
      else
        continue
      end
  end

  if all(isnan(x))
    continue
  end

  %calculate area
  area = stablepolyarea(resy * y, resx * x);
  if area > maxarea
    maxarea = area;
    slice = sliceloop;
  end
end

%------------------------------------------------------------
function [cineno,scarno,flowno,strainno,marno,varargout] = findno
%------------------------------------------------------------
%Find matching image stacks output

% normalno (normal short axis stack or closest equivalent)
% scarno (viability short axis stacks)
% flowno (flow image stack(s))
% marno (stacks with MaR)
% [stress] stress images
%
%Note that all scarno for instance is not garanteed to have non empty
%SET.Scar.

global SET

%Check what is possible to plot
cineno = [];
scarno = [];
flowno = [];
strainno = [];
marno = [];

if nargout>3
  varargout = cell(1,2);
end

if ~isfield(SET(1),'ImageType')
  SET(1).ImageType = '';
end
%Check with image types first!
for loop=1:length(SET)
  if isequal(lower(SET(loop).ImageType),'cine')
    cineno = [cineno loop]; %#ok<AGROW>
  end
  if isequal(lower(SET(loop).ImageType),'late enhancement') || ~isempty(SET(loop).Scar)
    scarno = [scarno loop]; %#ok<AGROW>
  end
  if isequal(lower(SET(loop).ImageType),'strain ffe') || isequal(lower(SET(loop).ImageType),'strain tfe') || ...
      ~isempty(SET(loop).Strain) || ~isempty(SET(loop).StrainTagging) || ~isempty(SET(loop).StrainMitt) || ...
      (~isempty(SET(loop).Flow) && isempty(SET(loop).Flow.PhaseNo) && ~isempty(SET(loop).Flow.PhaseX) && ~isempty(SET(loop).Flow.PhaseY) && isequal(SET(loop).Flow.MagnitudeNo,loop))
    strainno = [strainno loop]; %#ok<AGROW>
  end
  if isequal(SET(loop).ImageType,'Flow (magnitude)') && numel(SET(loop).Linked)> 1
    flowno = [flowno loop]; %#ok<AGROW>
  end
  if isequal(lower(SET(loop).ImageType),'perfusion rest') && ~isempty(SET(loop).MaR)
    marno = [marno loop]; %#ok<AGROW>
  end
  if isequal(lower(SET(loop).ImageType),'perfusion stress') && ~isempty(SET(loop).MaR)
    marno = [marno loop]; %#ok<AGROW>
  end

  %   if nargout>3
  %     if isequal(SET(loop).ImageType,'Stress baseline');
  %       varargout{1} = [varargout{1} loop];
  %     end;
  %   end;
  if nargout>4
    if isequal(lower(SET(loop).ImageType),'stress')
      varargout{2} = [varargout{2} loop];
    end
  end
end

%Continue with trying to find with automatic
for loop=1:length(SET)
  if not(isempty(SET(loop).Scar)) && isempty(find(scarno==loop,1)) && isempty(scarno)
    scarno = loop;
  else
    if not(isempty(SET(loop).Flow)) && isempty(find(flowno==loop,1))
      flowno = [flowno loop]; %#ok<AGROW>
    else
      if not(isempty(SET(loop).MaR)) && isempty(find(marno==loop,1))
        marno = [marno loop]; %#ok<AGROW>
      end %not mar

      %If neither scar,flow check if short axes?
      %Note it may be both mar and cine!
      if ~isempty(SET(loop).EndoX)||~isempty(SET(loop).EpiX)||~isempty(SET(loop).RVEndoX)||~isempty(SET(loop).RVEpiX)
        if (SET(loop).TSize>1)&&...
            (sum(sum(not(isnan([SET(loop).EndoX(:); SET(loop).EpiX(:); SET(loop).RVEndoX(:); SET(loop).RVEpiX(:)]))))>0)&&...
            isempty(find(cineno==loop,1))

          %             (sum(sum(not(isnan([SET(loop).EndoX(:) SET(loop).EpiX(:) SET(loop).RVEndoX(:) SET(loop).RVEpiX(:)]))))>0)&&...
          %             isempty(find(cineno==loop,1))
          cineno = [cineno loop]; %#ok<AGROW>
        end
      end

    end %not flow
  end %not scar
end

%-------------------------------------
function [saxno, success,cancelled] = findsaxno
%-------------------------------------
% find stack with short axis image

global SET NO DATA
cineno = findfunctions('findno');
saxno = find(strcmp({SET.ImageViewPlane},'Short-axis'));
zno = find([SET.ZSize] > 1);
saxno = intersect(cineno,union(saxno,zno));
success = 1;
cancelled = false;
if isempty(saxno)
  if ismember(NO,zno) && ~DATA.Silent 
    if yesno('Could not find image stack defined as Short-axis Cine. Use current stack?')
      saxno = NO;
      SET(saxno).ImageViewPlane = 'Short-axis';
      SET(saxno).ImageType = 'Cine';
      for p = find(DATA.ViewPanels==saxno)
        drawfunctions('drawpanel',p);
      end
    else
      success = 0;
      cancelled = 1;
      return
    end
  elseif ismember(NO,zno) && DATA.Silent
    saxno = NO;
    SET(saxno).ImageViewPlane = 'Short-axis';
    SET(saxno).ImageType = 'Cine';
  else
    if ~DATA.Silent
      myfailed('Could not find short-axis image stack');
    else
      disp('Could not find short-axis image stack');
    end
    success = 0;
    return
  end
end

if ismember(NO,saxno)
  saxno = NO;
else
  nos = cellfun(@(x)dprintf('Short-axis stack %d',x), ...
    num2cell(saxno),'UniformOutput',false);
  if ~DATA.Silent
    m = mymenu(dprintf('Select short-axis stack'), nos{:});
  else
    m = 1;
  end
  if m
    saxno = saxno(m);
  else
    success = 0;
    cancelled = 1;
    return;
  end
end

%---------------------------------------------
function cineshortaxisno = findcineshortaxisno(multiple, minscore)
%---------------------------------------------
%Find one [LV] or two [LV and RV] cine short axis stack
global SET

if nargin < 1
  multiple = false; %default, only return one no, true = return all cine sax
end
if nargin < 2
  minscore = 0; % default value to return any cine slice
end
%cineshortaxisno = [];
[cineno] = findno;

cineno = cineno(:)'; %ensure row vector

%If isempty then just return something.
if isempty(cineno)
  cineshortaxisno = [];
  return;
end

%Find maximal number of slices
maxslices = max(cat(1,SET(:).ZSize));

%Extract properties.
isshortaxis = false(size(cineno));
iscine = false(size(cineno));
haslvseg = false(size(cineno));
hasrvseg = false(size(cineno));
hasmaxslices = ones(size(cineno));
numberofsliceswithLV = zeros(size(cineno));
numberofsliceswithRV = zeros(size(cineno));
hasoverseven = zeros(size(cineno));

for loop=1:length(cineno)
  iscine(loop) = isequal(lower(SET(cineno(loop)).ImageType),'cine') || (SET(cineno(loop)).TSize>=20 && ~strncmp(SET(cineno(loop)).ImageType,'Perfusion',9));
  isshortaxis(loop) = isequal(lower(SET(cineno(loop)).ImageViewPlane),'short-axis');
  haslvseg(loop) = (~isempty(SET(cineno(loop)).EndoX)|~isempty(SET(cineno(loop)).EpiX));
  if haslvseg(loop)
    %check how many slices have lv segmentation
    if ~isempty(SET(cineno(loop)).EndoX)
      numberofsegmentedslices = findnumberofsegmentedslices(SET(cineno(loop)).EndoX);
    else
      numberofsegmentedslices = 0;
    end
    numberofsliceswithLV(loop) = numberofsliceswithLV(loop) + numberofsegmentedslices;
    if ~isempty(SET(cineno(loop)).EpiX)
      numberofsegmentedslices = findnumberofsegmentedslices(SET(cineno(loop)).EpiX);
    else
      numberofsegmentedslices = 0;
    end
    numberofsliceswithLV(loop) = numberofsliceswithLV(loop) + numberofsegmentedslices;
  end
  hasrvseg(loop) = (~isempty(SET(cineno(loop)).RVEndoX)|~isempty(SET(cineno(loop)).RVEpiX));
  if hasrvseg(loop)
    %check how many slices have lv segmentation
    if ~isempty(SET(cineno(loop)).RVEndoX)
      numberofsegmentedslices = findnumberofsegmentedslices(SET(cineno(loop)).RVEndoX);
    else
      numberofsegmentedslices = 0;
    end
    numberofsliceswithRV(loop) = numberofsliceswithRV(loop) + numberofsegmentedslices;
    if ~isempty(SET(cineno(loop)).RVEpiX)
      numberofsegmentedslices = findnumberofsegmentedslices(SET(cineno(loop)).RVEpiX);
    else
      numberofsegmentedslices = 0;
    end
    numberofsliceswithRV(loop) = numberofsliceswithRV(loop) + numberofsegmentedslices;
  end

  hasmaxslices(loop) = isequal(SET(cineno(loop)).ZSize,maxslices);
  hasoverseven(loop) = SET(cineno(loop)).ZSize > 7;
end

%Find if slices equals max
%Choose the best candidate in order of priority. Each row in A is one
%image stack, multiply with priority vector. Take largest sum.
A = [iscine' isshortaxis' haslvseg' hasoverseven' numberofsliceswithLV' hasrvseg' numberofsliceswithRV' hasmaxslices'];
score = A*[16;12;4;4;1;2;1;1]; %most important in priority cine, shortaxis, lvseg, numberofslicesLV, rvseg, numberofslicesLV, maxslices

if multiple
  %find lv and rv stacks that is cine and short-axis

  templvnos = find(haslvseg & isshortaxis);
  temprvnos = find(hasrvseg & isshortaxis);
  if ~isempty(templvnos)
    [~,lvind] = max(score(templvnos));
    lvno = cineno(templvnos(lvind));
  else
    highestvaluelv = max(score);
    lvnoinds = find(score==highestvaluelv);
    lvno = cineno(lvnoinds(end));
  end
  if ~isempty(temprvnos)
    [~,rvind] = max(score(temprvnos));
    rvno = cineno(temprvnos(rvind));
  else
    highestvaluerv = max(score);
    rvnoinds = find(score==highestvaluerv);
    rvno = cineno(rvnoinds(end));
  end
  cineshortaxisno = [lvno rvno];
else
  %   score = A*[16;8;4;2;1]; %most important in priority cine, shortaxis, lvseg, rvseg, maxslices
  [scorevalue,ind] = max(score);
  if scorevalue < minscore
    cineshortaxisno = [];
  else
    cineshortaxisno = cineno(ind);
  end
end

%------------------------------------------------
function numberofsegmentedslices = findnumberofsegmentedslices(arraywithsegmentation)
%------------------------------------------------
%Function to find number of slices in an 3D array
%(for example EndoX/EpiX/RVEndoX/RVEpiX)that contain segmentation

% permuting array so that the slice dimension is the first one
arraywithsegmentation   = permute(arraywithsegmentation,[3 1 2]);
% reshaping the array into 2D array with size [Number of Slices x rest]
arraywithsegmentation   = reshape(arraywithsegmentation, size(arraywithsegmentation,1),[]);
% calculating whether each slice contains segmentation points
numberofpointsperslice  = sum(~isnan(arraywithsegmentation),2);
% get the number of slices containing segmentations
numberofsegmentedslices = nnz(numberofpointsperslice);



%------------------------------------------------
function shortaxisno = findctshortaxisno(multiple)
%------------------------------------------------
%Find one [LV] or two [LV and RV] CT short axis stack
global SET

if nargin < 1
  multiple = false; %default, only return one no, true = return all cine sax
end

ctno = find(strcmp('CTheart',{SET.ImagingTechnique})); %ensure row vector

%If isempty then just return something.
if isempty(ctno)
  shortaxisno = [];
  return;
end

%Find maximal number of slices
maxslices = max(cat(1,SET(:).ZSize));

%Extract properties.
% iscine = false(size(ctno));
isshortaxis = false(size(ctno));
haslvseg = false(size(ctno));
hasrvseg = false(size(ctno));
hasmaxslices = ones(size(ctno));

for loop=1:length(ctno)
  %   iscine(loop) = isequal(lower(SET(ctno(loop)).ImageType),'cine') || (SET(ctno(loop)).TSize>=20 && ~strncmp(SET(ctno(loop)).ImageType,'Perfusion',9));
  isshortaxis(loop) = isequal(lower(SET(ctno(loop)).ImageViewPlane),'short-axis');
  haslvseg(loop) = (~isempty(SET(ctno(loop)).EndoX)|~isempty(SET(ctno(loop)).EpiX));
  hasrvseg(loop) = (~isempty(SET(ctno(loop)).RVEndoX)|~isempty(SET(ctno(loop)).RVEpiX));
  hasmaxslices(loop) = isequal(SET(ctno(loop)).ZSize,maxslices);
end

%Find if slices equals max
%Choose the best candidate in order of priority. Each row in A is one
%image stack, multiply with priority vector. Take largest sum.
A = [haslvseg' isshortaxis'  hasrvseg' hasmaxslices'];

if multiple
  %find lv and rv stacks that is cine and short-axis
  score = A*[8;4;2;1]; %most important in priority cine, lvseg, shortaxis, rvseg, maxslices
  templvnos = find(haslvseg);
  temprvnos = find(hasrvseg);
  if ~isempty(templvnos)
    [~,lvind] = max(score(templvnos));
    lvno = ctno(templvnos(lvind));
  else
    [~,lvnoind] = max(score);
    lvno = ctno(lvnoind);
  end
  if ~isempty(temprvnos)
    [~,rvind] = max(score(temprvnos));
    rvno = ctno(temprvnos(rvind));
  else
    [~,rvnoind] = max(score);
    rvno = ctno(rvnoind);
  end
  shortaxisno = [lvno rvno];
else
  score = A*[16;8;4;2;1]; %most important in priority cine, lvseg, shortaxis, rvseg, maxslices
  [~,ind] = max(score);
  shortaxisno = ctno(ind);
end

%----------------------------------------------------
function shortaxisno = findspectshortaxisno(multiple)
%----------------------------------------------------
%Find one [LV] or two [LV and RV] CT short axis stack
global SET

if nargin < 1
  multiple = false; %default, only return one no, true = return all cine sax
end

spectno = find(strcmp('NM',{SET.ImagingTechnique})); %ensure row vector

%If isempty then just return something.
if isempty(spectno)
  shortaxisno = [];
  return;
end

%Find maximal number of slices
maxslices = max(cat(1,SET(:).ZSize));

%Extract properties.
% iscine = false(size(spectno));
isshortaxis = false(size(spectno));
haslvseg = false(size(spectno));
hasrvseg = false(size(spectno));
hasmaxslices = ones(size(spectno));

for loop=1:length(spectno)
  %   iscine(loop) = isequal(lower(SET(spectno(loop)).ImageType),'cine') || (SET(spectno(loop)).TSize>=20 && ~strncmp(SET(spectno(loop)).ImageType,'Perfusion',9));
  isshortaxis(loop) = isequal(lower(SET(spectno(loop)).ImageViewPlane),'short-axis');
  haslvseg(loop) = (~isempty(SET(spectno(loop)).EndoX)|~isempty(SET(spectno(loop)).EpiX));
  hasrvseg(loop) = (~isempty(SET(spectno(loop)).RVEndoX)|~isempty(SET(spectno(loop)).RVEpiX));
  hasmaxslices(loop) = isequal(SET(spectno(loop)).ZSize,maxslices);
end

%Find if slices equals max
%Choose the best candidate in order of priority. Each row in A is one
%image stack, multiply with priority vector. Take largest sum.
A = [haslvseg' isshortaxis'  hasrvseg' hasmaxslices'];

if multiple
  %find lv and rv stacks that is cine and short-axis
  score = A*[8;4;2;1]; %most important in priority cine, lvseg, shortaxis, rvseg, maxslices
  templvnos = find(haslvseg);
  temprvnos = find(hasrvseg);
  if ~isempty(templvnos)
    [~,lvind] = max(score(templvnos));
    lvno = spectno(templvnos(lvind));
  else
    [~,lvnoind] = max(score);
    lvno = spectno(lvnoind);
  end
  if ~isempty(temprvnos)
    [~,rvind] = max(score(temprvnos));
    rvno = spectno(temprvnos(rvind));
  else
    [~,rvnoind] = max(score);
    rvno = spectno(rvnoind);
  end
  shortaxisno = [lvno rvno];
else
  score = A*[16;8;4;2;1]; %most important in priority cine, lvseg, shortaxis, rvseg, maxslices
  [~,ind] = max(score);
  shortaxisno = spectno(ind);
end

%------------------------------------------
function scarno = findscarshortaxisno
%------------------------------------------
%Find only one scar shortaxis stack
global SET DATA

[~,scarno] = findno;
preferMAG = true; %default value
if isfield(DATA.Pref, 'ChooseMAG')
  preferMAG = DATA.Pref.ChooseMAG;
end

%--- Check if multiple scardata.
if length(scarno)>1
  %Find best scar data to take. Take those with endo and scar data
  scar2use = zeros(size(scarno));
  for sloop=1:length(scarno)
    if existfunctions('existendo',scarno(sloop))
      scar2use(sloop) = scar2use(sloop)+0.5;
    end
    if existfunctions('existepi',scarno(sloop))
      scar2use(sloop) = scar2use(sloop)+0.5;
    end
    if not(isempty(SET(scarno(sloop)).Scar))
      scar2use(sloop) = scar2use(sloop)+4;
    end
    if isequal(lower(SET(scarno(sloop)).ImageViewPlane),'short-axis')
      scar2use(sloop) = scar2use(sloop)+2;
    end
    if SET(scarno(sloop)).ZSize > 5
      scar2use(sloop) = scar2use(sloop)+2;
    end
    if SET(scarno(sloop)).ZSize > 1
      scar2use(sloop) = scar2use(sloop)+1;
    end
    if contains(lower(SET(scarno(sloop)).SeriesDescription), 'psir')
      if preferMAG
        scar2use(sloop) = scar2use(sloop)-1;
      else
        scar2use(sloop) = scar2use(sloop)+1;
      end
    end
     if contains(lower(SET(scarno(sloop)).SeriesDescription), 'mag')
       if preferMAG
         scar2use(sloop) = scar2use(sloop)+1;
       else
         scar2use(sloop) = scar2use(sloop)-1;
       end
     end
     if contains(lower(SET(scarno(sloop)).DICOMImageType), 'm_ir')
       if preferMAG
         scar2use(sloop) = scar2use(sloop)+1;
       else
         scar2use(sloop) = scar2use(sloop)-1;
       end
     end
     if contains(lower(SET(scarno(sloop)).DICOMImageType), 'm_ffe')
       if preferMAG
         scar2use(sloop) = scar2use(sloop)-1;
       else
         scar2use(sloop) = scar2use(sloop)+1;
       end
     end
  end
  [maxval,~] = max(scar2use);
  allmax = find(scar2use == maxval);
  scarno = scarno(allmax); %#ok<FNDSB>

  oldscarno = scarno;
  take = true(size(scarno));
  for s = 1:length(scarno)
    if isempty(SET(scarno(s)).Scar)
      take(s) = false;
    end
  end
  scarno = scarno(take);
  if isempty(scarno)
    scarno = oldscarno;
  end

  if length(scarno)>1
    %mywarning('Detected multiple scar data. Taking data with maximal scar volume (arbitrary decision)');
    s = scarno(end); %if there is no scar data, take latest LGE stack
    maxml = 0;
    for sloop = 1:length(scarno)
      try
        ml = SET(scarno(sloop)).Scar.Percentage*SET(scarno(sloop)).LVM/100; %scar volume in ml
        if ml>maxml
          maxml = ml;
          s = scarno(sloop);
        end
      catch  %#ok<CTCH>
      end
    end
    scarno = s;
%   elseif length(scarno)>1
%     disp('Detected multiple scar data. Taking first image stack (arbitrary decision)');
%     scarno = scarno(1);
  end
end

%------------------------------------------
function marno = findmarshortaxisno
%------------------------------------------
%Find only one mar shortaxis stack
global SET

[~,~,~,~,marno] = findno;

%--- Check if multiple mardata.
if length(marno)>1
  %Find best mar data to take. Take those with endo and mar data
  mar2use = false(size(marno));
  for sloop=1:length(marno)
    mar2use(sloop) = existfunctions('existendo',marno(sloop)) && not(isempty(SET(marno(sloop)).MaR)) && max(SET(marno(sloop)).MaR.Percentage)>0;
  end
  marno = marno(mar2use);

  if length(marno)>1
    disp('Detected multiple MaR data. Taking data with maximal MaR volume (arbitrary decision)');
    s = marno(1);
    maxml = 0;
    for sloop = 1:length(marno)
      try
        ml = SET(marno(sloop)).MaR.Percentage*SET(marno(sloop)).LVM/100; %scar volume in ml
        if ml>maxml
          maxml = ml;
          s = marno(sloop);
        end
      catch  %#ok<CTCH>
      end
    end
    marno = s;
  end
end

%-----------------------------------------
function [flowno,flowroi] = findflowaxisno(flowno)
%-----------------------------------------
%Find only one flow image stack
global SET NO

if nargin < 1
  [~,~,flowno] = findno;
  magno = [];
  for noloop = 1:length(flowno)
    if ~isempty(SET(flowno(noloop)).Flow) && isfield(SET(flowno(noloop)).Flow,'MagnitudeNo')
      magno = [magno SET(flowno(noloop)).Flow.MagnitudeNo]; %#ok<AGROW> 
    end
  end
  flowno = unique(magno);
end

%--- Check if multiple flow image stacks.
if length(flowno)>1
  %Find best flow data to take according to ROI name
  flow2use = false(size(flowno));
  for sloop=1:length(flowno)
    flow2use(sloop) = (SET(flowno(sloop)).RoiN>0);
  end
  [maxval,~] = max(flow2use);
  allmax = find(flow2use == maxval);
  flowno = flowno(allmax); %#ok<FNDSB> 
  if length(flowno) > 1
    maxpoints = zeros(length(flowno),1);
    %find best image stack based on ROI names
    index = 1;
    for no = flowno
      points = zeros(SET(no).RoiN,1);
      for rloop = 1:SET(no).RoiN
        hasresult = (length(SET(no).Flow.Result)>=rloop && isempty(SET(no).Flow.Result(rloop)));
        if isequal(SET(no).Roi(rloop).Name,'Aortic ascending flow') && hasresult
          points(rloop) = 9;
        elseif isequal(SET(no).Roi(rloop).Name,'Pulmonary artery') && hasresult
          points(rloop) = 8;
        elseif not(isempty(strfind(SET(no).Roi(rloop).Name,'Static tissue'))) && hasresult
          points(rloop) = 7;
        elseif hasresult
          points(rloop) = 5;
        elseif isequal(SET(no).Roi(rloop).Name,'Aortic ascending flow')
          points(rloop) = 4;
        elseif isequal(SET(no).Roi(rloop).Name,'Pulmonary artery')
          points(rloop) = 3;
        elseif not(isempty(strfind(SET(no).Roi(rloop).Name,'Static tissue')))
          points(rloop) = 2;
        end
      end
      if ~isempty(points)
        maxpoints(index) = max(points);
      end
      index = index +1;
    end
    [maxval,~] = max(maxpoints);
    allmax = find(maxpoints == maxval);
    flowno = flowno(allmax);
  end
  if length(flowno)>1
    if ismember(NO,flowno)
      flowno = NO; %taking current image stack
    else
      %disp('Detected multiple flow image stacks.Taking first stack (arbitrary decision)');
      flowno = flowno(1);
    end
  end
end

%find flow ROI based on ROI name
if ~isempty(flowno)
  points = zeros(SET(flowno).RoiN,1);
  for rloop = 1:SET(flowno).RoiN
    hasresult = (length(SET(flowno).Flow.Result)>=rloop && isempty(SET(flowno).Flow.Result(rloop)));
    if isequal(SET(flowno).Roi(rloop).Name,'Aortic ascending flow') && hasresult
      points(rloop) = 9;
    elseif isequal(SET(flowno).Roi(rloop).Name,'Pulmonary artery') && hasresult
      points(rloop) = 8;
    elseif not(isempty(strfind(SET(flowno).Roi(rloop).Name,'Static tissue'))) && hasresult
      points(rloop) = 7;
    elseif hasresult
      points(rloop) = 5;
    elseif isequal(SET(flowno).Roi(rloop).Name,'Aortic ascending flow')
      points(rloop) = 4;
    elseif isequal(SET(flowno).Roi(rloop).Name,'Pulmonary artery')
      points(rloop) = 3;
    elseif not(isempty(strfind(SET(flowno).Roi(rloop).Name,'Static tissue')))
      points(rloop) = 2;
    end
  end
  if ~isempty(points)
    [maxval,~] = max(points);
  else
    maxval = 0;
  end
  flowroi = find(points == maxval);
  if length(flowroi)>1
    if ismember(SET(flowno).RoiCurrent,flowroi)
      flowroi = SET(flowno).RoiCurrent;
    else
      flowroi = flowroi(1);
    end
  end
else
  flowroi = [];
end


%------------------------------------------
function [txmapno, success] = findt1stackno
%------------------------------------------
% find stack applicable for T1 analysis
global SET NO DATA

t1stacks = findallt1stacks;
success = 1;

%Check if any T1 stack was found
if isempty(t1stacks)
  str = dprintf('Could not find T%s image stack','1');
  if ~DATA.Silent
    myfailed(str);
  else
    disp(str);
  end
  txmapno = [];
  success = 0;
  return
end

%if NO is a T1 stack, take that, if multiple select
if ismember(NO,t1stacks)
  txmapno = NO;
elseif length(t1stacks) > 1
  nos = cellfun(@(x)sprintf('%s %d  %s',dprintf('Stack'),x,SET(x).ImageType), ...
    num2cell(t1stacks),'UniformOutput',false);
  if ~DATA.Silent
    str = dprintf('Select T%s stack','1');
    m = mymenu(str, nos{:});
  else
    m = 1;
  end
  if m
    txmapno = t1stacks(m);
  else
    txmapno = [];
    success = 0;
    return;
  end
else
  txmapno = t1stacks(1);
end

%----------------------------------
function t1stacks = findallt1stacks
%----------------------------------
% find stack applicable for T1 analysis
global SET

%search for multiple, and unique, inversion times
points = zeros(1,length(SET));
for no = 1:length(SET)
  nbrinversiontime = size(SET(no).InversionTime,2);
  diffinversiontime = min(abs(diff(SET(no).InversionTime,1,2)),[],'all');
  if nbrinversiontime > 1
    points(no) = points(no)+2;
  end
  if not(isempty(diffinversiontime)) && diffinversiontime > 0
    points(no) = points(no)+1;
  end
end
maxpoints = max(points);

if maxpoints < 3
  t1stacks = [];
else
  t1stacks = find(points==maxpoints);
end

%------------------------------------------
function [txmapno, success] = findt2stackno
%------------------------------------------
% find stack applicable for T2 analysis
global SET NO DATA

txmapstacks = findallt2stacks;
success = 1;

%Check if any T2 stack was found
if isempty(txmapstacks)
  str = dprintf('Could not find T%s image stack','2');
  if ~DATA.Silent
    myfailed(str);
  else
    disp(str);
  end
  txmapno = [];
  success = 0;
  return
end

%if NO is a T2 stack, take that, if multiple select
if ismember(NO,txmapstacks)
  txmapno = NO;
elseif length(txmapstacks) > 1
  nos = cellfun(@(x)sprintf('%s %d  %s',dprintf('Stack'),x,SET(x).ImageType), ...
    num2cell(txmapstacks),'UniformOutput',false);
  if ~DATA.Silent
    str = dprintf('Select T%s stack','2');
    m = mymenu(str, nos{:});
  else
    m = 1;
  end
  if m
    txmapno = txmapstacks(m);
  else
    txmapno = [];
    success = 0;
    return;
  end
else
  txmapno = txmapstacks(1);
end

%------------------------------------------
function t2stacks = findallt2stacks
%------------------------------------------
% find all stacks applicable for T2 analysis
global SET

%search for multiple, and unique, t2prep or echotime
points = zeros(1,length(SET));
for no = 1:length(SET)
  pointsprep = 0;
  pointsecho = 0;
  if isfield(SET(no),'T2preptime')
    nbrpreptime = size(SET(no).T2preptime,2);
    diffpreptime = min(abs(diff(SET(no).T2preptime,1,2)),[],'all');
    if nbrpreptime > 1
      pointsprep = pointsprep+2;
    end
    if not(isempty(diffpreptime)) && diffpreptime > 0
      pointsprep = pointsprep+1;
    end
  end
  nbrechotime = size(SET(no).EchoTime,2);
  diffechotime = min(abs(diff(SET(no).EchoTime,1,2)),[],'all');
  if nbrechotime > 1
    pointsecho = pointsecho+2;
  end
  if not(isempty(diffechotime)) && diffechotime > 0
    pointsecho = pointsecho+1;
  end
  points(no) = max([pointsecho, pointsprep]);
end
maxpoints = max(points);
if maxpoints < 3
  t2stacks = [];
else
  t2stacks = find(points==maxpoints);
end

%----------------------------------------------
function [txmapno, success] = findt2starstackno
%----------------------------------------------
% find stack applicable for T2* analysis
global SET NO DATA

txmapstacks = findallt2starstacks;
success = 1;

%Check if any T2 stack was found
if isempty(txmapstacks)
  str = dprintf('Could not find T%s image stack','2*');
  if ~DATA.Silent
    myfailed(str);
  else
    disp(str);
  end
  txmapno = [];
  success = 0;
  return
end

%if NO is a T2* stack, take that, if multiple select
if ismember(NO,txmapstacks)
  txmapno = NO;
elseif length(txmapstacks) > 1
  nos = cellfun(@(x)sprintf('%s %d  %s',dprintf('Stack'),x,SET(x).ImageType), ...
    num2cell(txmapstacks),'UniformOutput',false);
  if ~DATA.Silent
    str = dprintf('Select T%s stack','2*');
    m = mymenu(str, nos{:});
  else
    m = 1;
  end
  if m
    txmapno = txmapstacks(m);
  else
    txmapno = [];
    success = 0;
    return;
  end
else
  txmapno = txmapstacks(1);
end

%----------------------------------------------
function txmapstacks = findallt2starstacks
%----------------------------------------------
% find all stacks applicable for T2* analysis
global SET

%search for multiple, and unique, echotime
points = zeros(1,length(SET));
for no = 1:length(SET)
  nbrechotime = size(SET(no).EchoTime,2);
  diffechotime = min(abs(diff(SET(no).EchoTime,1,2)),[],'all');
  if nbrechotime > 1
    points(no) = points(no)+2;
  end
  if not(isempty(diffechotime)) && diffechotime > 0
    points(no) = points(no)+1;
  end
end
maxpoints = max(points);
if maxpoints < 3
  txmapstacks = [];
else
  txmapstacks = find(points==maxpoints);
end

%--------------------------------------
function ind = findslicewithscarall(no)
%--------------------------------------
%Find and return slices with scar
global SET NO

if nargin < 1
  no = NO;
end

if isempty(SET(no).Scar)
  ind = false(SET(no).ZSize,1);
  return;
end
ind = squeeze(find(sum(sum(SET(no).Scar.Result(:,:,:)))~=0));

%--------------------------------------
function ind = findslicewithmarall(no)
%--------------------------------------
%Find and return slices with mar
global SET NO

if nargin < 1
  no = NO;
end

if isempty(SET(no).MaR)
  ind = false(SET(no).ZSize,1);
  return;
end

if SET(no).TSize>1
  ind = squeeze(find(sum(sum(sum(SET(no).MaR.Result(:,:,:,:))))~=0));
else
  ind = squeeze(find(sum(sum(SET(no).MaR.Result(:,:,:,:)))~=0));
end

%--------------------------------------
function ind = findslicewithendoall(no)
%--------------------------------------
%Find slices with endocard in all timeframes
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EndoX)
  ind = false(SET(no).ZSize,1);
  return;
end

if SET(no).TSize>1
  temp = not(isnan(squeeze(SET(no).EndoX(1,:,:))));
  if SET(no).ZSize==1
    ind = all(temp)';
  else
    ind = sum(temp,1)'==SET(no).TSize;
  end
else
  ind = squeeze(not(isnan(SET(no).EndoX(1,1,:))));
end

%-------------------------------------
function ind = findslicewithepiall(no)
%-------------------------------------
%Find slices with endocard in all timeframes
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).EpiX)
  ind = false(SET(no).ZSize,1);
  return;
end

if SET(no).TSize>1
  temp = not(isnan(squeeze(SET(no).EpiX(1,:,:))));
  if SET(no).ZSize==1
    ind = all(temp)';
  else
    ind = sum(temp,1)'==SET(no).TSize;
  end
else
  ind = squeeze(not(isnan(SET(no).EpiX(1,1,:))));
end

%---------------------------------------
function ind = findslicewithendo(no,tfs)
%---------------------------------------
%Find slices with endocard in any timeframe
global SET NO

if nargin<1
  no = NO;
end
if nargin < 2
  tfs = 1:SET(no).TSize;
end

if isempty(SET(no).EndoX)
  ind = false(SET(no).ZSize,1);
  return;
end

temp = not(isnan(squeeze(SET(no).EndoX(1,tfs,:))));
if length(tfs)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end
else
  ind = temp;
end

%-------------------------------------
function ind = findslicewithrvendo(no,tfs)
%-------------------------------------
%Find slices with RV endocard in any timeframe
global SET NO

if nargin<1
  no = NO;
end
if nargin < 2
  tfs = 1:SET(no).TSize;
end

if isempty(SET(no).RVEndoX)
  ind = false(SET(no).ZSize,1);
  return;
end

temp = not(isnan(squeeze(SET(no).RVEndoX(1,tfs,:))));
if length(tfs)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end
else
  ind = temp;
end

%-------------------------------------
function ind = findslicewithgeneralsegmentation(no,objectind)
%-------------------------------------
%Find slices with General Pen segmentation in any timeframe
global SET

numberoftimeframes = 1:SET(no).TSize;

if isempty(SET(no).GeneralPenObjects(objectind).X)
  ind = false(SET(no).ZSize,1);
  return
end

temp = not(isnan(squeeze(SET(no).GeneralPenObjects(objectind).X(1,numberoftimeframes,:))));
if length(numberoftimeframes)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end
else
  ind = temp;
end

%--------------------------------------
function ind = findslicewithrvendoall(no)
%--------------------------------------
%Find slices with RV endocard in all timeframes
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).RVEndoX)
  ind = false(SET(no).ZSize,1);
  return;
end

if SET(no).TSize>1
  temp = not(isnan(squeeze(SET(no).RVEndoX(1,:,:))));
  if SET(no).ZSize==1
    ind = all(temp);
  else
    ind = sum(temp,1)==SET(no).TSize;
  end
else
  ind = squeeze(not(isnan(SET(no).RVEndoX(1,1,:))));
end


%---------------------
function [measureind, pointind,measurex,measurey,measurez,mindist] = closestmeasure(panel,xclick,yclick)
%-----------------------
%finds the closest visible measure. returns the index of the point within
%the closest measure aswell as all corrdinates for the closest measure.
global DATA SET

no = DATA.ViewPanels(panel);
%No measures avaiable
if isempty(SET(no).Measure)
  measureind = [];
  pointind = [];
  measurex = [];
  measurey = [];
  measurez = [];
  mindist = inf;
  return
end

%retrieve coordinates of measures in panel
[measure,~] = viewfunctions('getmeasurecoords',panel);

%find which meaures to consider
slicestodo = viewfunctions('slicesinpanel',panel);
measurestocheck = find(cellfun(@(x,y) any(ismember(x,slicestodo)) && ...
  any(y==SET(no).CurrentTimeFrame) ,...
  {measure.Z},{measure.T}));

measurex = cell(1,length(measurestocheck));
measurey = measurex;
measurez = measurex;

for i = 1:length(measurestocheck)
  measurex{i} = measure(measurestocheck(i)).X;
  measurey{i} = measure(measurestocheck(i)).Y;
  measurez{i} = measure(measurestocheck(i)).Z;
end

if ~isempty(DATA.LastClickType) && strcmp(DATA.LastClickType ,'alt')
  [mindists,mininds] = cellfun(@(x,y) min(calcfunctions('calcdistpointoline',[xclick yclick 0],x,y)),measurey,measurex);
else
  [mindists,mininds] = cellfun(@(x,y) min(sqrt((xclick-x).^2+(yclick-y).^2)),measurey,measurex);
end
[mindist,ind] = min(mindists);

if mindist<DATA.Pref.ContourAdjustDistance
  measureind = measurestocheck(ind);
  pointind = mininds(ind);
  measurex = measurex{ind};
  measurey = measurey{ind};
  measurez = measurez{ind};
elseif ~isempty(mindist)
  measureind = measurestocheck(ind);
  pointind = [];
  measurex = [];
  measurey = [];
  measurez = [];
else
  measureind = [];
  pointind = [];
  measurex = [];
  measurey = [];
  measurez = [];
end

if isempty(mindist)
  mindist = inf;
end

%------------------------------
function [mindist, ind] = dist2contour(type,no,x,y,slice,tf)
%-----------------------------
%calcualte the minimal distance to contour
global SET

if ~isempty(SET(no).([type,'X']))
  contourx = SET(no).([type,'X'])(:,tf,slice);
  contoury = SET(no).([type,'Y'])(:,tf,slice);
  [mindist,ind] = min(sqrt((contourx-x).^2+(contoury-y).^2));
else
  mindist = inf;
  ind = 0;
end

%------------------------------
function [mindist, objind, pointind] = closestgeneralpen(no,x,y,slice,tf)
%-----------------------------
% Calcualte the minimal distance between the clicked position and a general
% pen contour.
% x, y are clicked coordinates in the given slice and time frame.
% mindist is the distance between clikced position and the closest general
% pen contour.
% objind is the object index in GeneralPenObjects
% pointind is the index of the closest point in the General Pen contour
global SET

if ~isempty(SET(no).GeneralPenObjects)
  objectslist = SET(no).GeneralPenObjects;
  n = numel(objectslist);
  contourx = cell(1,n);
  contoury = cell(1,n);
  for objind = 1:n
    contourx{objind} = objectslist(objind).X(:,tf,slice);
    contoury{objind} = objectslist(objind).Y(:,tf,slice);
  end
  [mindistarray,ind] = cellfun(@(contourx,contoury) min(sqrt((x-contourx).^2+(y-contoury).^2)),contourx,contoury);
  [mindist,objind] = min(mindistarray);
  pointind = ind(objind);
else
  mindist = inf;
  objind = 0;
  pointind = 0;
end

%------------------------------
function [mindist, objind, pointind] = closestgeneralpeninterp(no,x,y,slice,tf)
%-----------------------------
% Calcualte the minimal distance between the clicked position and a general
% pen interpolation point.
% x, y are clicked coordinates in the given slice and time frame.
% mindist is the distance between clikced position and the closest general
% pen contour.
% objind is the object index in GeneralPenObjects
% pointind is the index of the closest point in the General Pen contour
global SET

mindist = inf;
objind = 0;
pointind = 0;

if ~isempty(SET(no).GeneralPenObjects)
  objectslist = SET(no).GeneralPenObjects;
  n = numel(objectslist);
  contourx = cell(1,n);
  contoury = cell(1,n);
  for objind = 1:n
    if (~isempty(objectslist(objind).InterpX) && ~isempty(objectslist(objind).InterpX{tf,slice}))
      contourx{objind} = objectslist(objind).InterpX{tf,slice};
      contoury{objind} = objectslist(objind).InterpY{tf,slice};
    end
  end
  if any(cellfun(@isempty,contourx))
    return
  end
  [mindistarray,ind] = cellfun(@(contourx,contoury) min(sqrt((x-contourx).^2+(y-contoury).^2)),contourx,contoury);
  [mindist,objind] = min(mindistarray);
  pointind = ind(objind);
end

%-------------------------------------------------
function type = closestheartside(panel,x,y,slice,tf)
%--------------------------------------------------
%find the closest object
global DATA SET

no = DATA.ViewPanels(panel);

if nargin < 4
  slice = SET(no).CurrentSlice;
end

if nargin < 5
  tf = SET(no).CurrentTimeFrame;
end

%Find closest segmentation in current timeframe and slice
minendo = dist2contour('Endo',no,x,y,slice,tf);
minepi = dist2contour('Epi',no,x,y,slice,tf);
minrvendo = dist2contour('RVEndo',no,x,y,slice,tf);
minrvepi = dist2contour('RVEpi',no,x,y,slice,tf);
% distance to LA/RA contour
[mindistla,~] = closestlaracontour(panel,'LA',x,y,slice,tf);
[mindistra,~] = closestlaracontour(panel,'RA',x,y,slice,tf);

%Find closest interp in current timeframe and slice
[minendointerp,~] = closestinterp(panel,'EndoInterp',x,y,slice,tf);
[minepiinterp,~] = closestinterp(panel,'EpiInterp',x,y,slice,tf);
[minrvendointerp,~] = closestinterp(panel,'RVEndoInterp',x,y,slice,tf);
[minrvepiinterp,~] = closestinterp(panel,'RVEpiInterp',x,y,slice,tf);
[minlainterp,~] = closestlarainterp(panel,'LAInterp',x,y,slice,tf);
[minrainterp,~] = closestlarainterp(panel,'RAInterp',x,y,slice,tf);

%summarize which is closest
distright = min([minrvendo,minrvepi,...
  minrvendointerp,minrvepiinterp,...
  minrainterp,mindistra]);

distleft = min([minendo,minepi,...
  minendointerp,minepiinterp,...
  minlainterp,mindistla]);
if distleft > distright
  % right side is closer
  type = 'right';
else
  type = 'left';
end

%------------------------------
function [type,objind] = closestobject(panel,x,y,slice,tf)
%-----------------------
%find the closest object
global DATA SET

no = DATA.ViewPanels(panel);

if nargin < 4
  slice = SET(no).CurrentSlice;%clickedslice(panel);
end

if nargin < 5
  tf = SET(no).CurrentTimeFrame;
end

%Find closest segmentation in current timeframe and slice
minendo = dist2contour('Endo',no,x,y,slice,tf);
minepi = dist2contour('Epi',no,x,y,slice,tf);
minrvendo = dist2contour('RVEndo',no,x,y,slice,tf);
minrvepi = dist2contour('RVEpi',no,x,y,slice,tf);
if strcmp(DATA.ProgramName, 'Segment')
  mincenter = min(sqrt((SET(no).CenterX-x).^2+(SET(no).CenterY-y).^2));
else
  mincenter = Inf;
end

[mindistla,~] = closestlaracontour(panel,'LA',x,y,slice,tf);
[mindistra,~] = closestlaracontour(panel,'RA',x,y,slice,tf);

%Find closest measure in current timeframe and slice
[measureind,~,~,~,~,minmeasure] = closestmeasure(panel,y,x);

%Find closest roi in current timeframe and slice
[minroi,roiind] = closestroi(panel,x,y,slice,tf);

%Find closest point in current timeframe and slice
[minpoint,pointind] = closestpoint(panel,x,y,slice,tf);

%Find closest interp in current timeframe and slice
[minendointerp,endointerpind] = closestinterp(panel,'EndoInterp',x,y,slice,tf);
[minepiinterp,epiinterpind] = closestinterp(panel,'EpiInterp',x,y,slice,tf);
[minrvendointerp,rvendointerpind] = closestinterp(panel,'RVEndoInterp',x,y,slice,tf);
[minrvepiinterp,rvepiinterpind] = closestinterp(panel,'RVEpiInterp',x,y,slice,tf);
[minlainterp,lainterpind] = closestlarainterp(panel,'LAInterp',x,y,slice,tf);
[minrainterp,rainterpind] = closestlarainterp(panel,'RAInterp',x,y,slice,tf);


%Find closest general pen object in current timeframe and slice
[mingeneralpen, generalpenind] = closestgeneralpen(no,x,y,slice,tf);
[mingeneralpeninterp,generalpeninterpind] = closestgeneralpeninterp(no,x,y,slice,tf);

%summarize which is closest
dists = [minmeasure,minroi,minendo,minepi,minrvendo,minrvepi,...
  minendointerp,minepiinterp,minrvendointerp,minrvepiinterp,minpoint,... 
  mincenter,mingeneralpen,mingeneralpeninterp,minlainterp,minrainterp,mindistla,mindistra];

[mindist,typeind] = min(dists);

types = {'Measure','Roi','Endo','Epi','RVEndo','RVEpi','EndoInterp',...
  'EpiInterp','RVEndoInterp','RVEpiInterp','Point','Center',...
  'GeneralPen','GeneralPenInterp','LAInterp','RAInterp','LA','RA'};

contourswithinterpolation = {'Endo','Epi','RVEpi','RVEndo','LA','RA','GeneralPen'};

type = types{typeind};

%If selected type is contour then check if Interp is within intervall and
%set Interp version of contour as type
if any(matches(contourswithinterpolation,type))
  interptype = [type,'Interp'];
  if dists(matches(types,interptype)) <= DATA.Pref.ContourAdjustDistance
    type = interptype;
  end
end

%if smallest distance is more than contour adjust distance then return
%image type for context menu
if mindist > DATA.Pref.ContourAdjustDistance
  type = 'Image';
end

%if the closest clicked item is a measure roi or interp point return index to SET struct location
switch type
  case 'Measure'
    objind = measureind;
  case 'Roi'
    objind = roiind;
  case 'EndoInterp'
    objind = endointerpind;
  case'EpiInterp'
    objind = epiinterpind;
  case'RVEndoInterp'
    objind = rvendointerpind;
  case 'RVEpiInterp'
    objind = rvepiinterpind;
  case 'Point'
    objind = pointind;
  case 'GeneralPen'
    objind = generalpenind;
  case 'GeneralPenInterp'
    objind = generalpeninterpind;
  case 'LAInterp'
    objind = lainterpind;
  case 'RAInterp'
    objind = rainterpind;
  otherwise %objind is irrelevant
    objind = 0;
end

%------------------------------
function [dist,ind] = closestpoint(panel,x,y,slice,tf)
%-----------------------
%find the closest annotation point
global DATA SET

no = DATA.ViewPanels(panel);

if nargin < 4
  slice = SET(no).CurrentSlice;
end

if nargin < 5
  tf = SET(no).CurrentTimeFrame;
end

if ismember(DATA.ViewPanelsType{panel},{'trans3DP','viewport','sag3DP','cor3DP'})
  %If 3DP then all slices
  xi = SET(no).Point.X;
  yi = SET(no).Point.Y;
  zi = SET(no).Point.Z;
  [dist,ind] = min(sqrt((x-xi).^2+(y-yi).^2+(slice-zi).^2));

else
  %get coordinates of all relevant points
  if not(isempty(SET(no).Point.X))
    pointstodo = find(ismember(SET(no).Point.Z,slice)&...
      (ismember(SET(no).Point.T,tf)|isnan(SET(no).Point.T))); %nan means constant over time
    xi = SET(no).Point.X(pointstodo);
    yi = SET(no).Point.Y(pointstodo);
  else
    dist = inf;
    ind = [];
    return
  end

  %Compute distance
  [dist,ind] = min(sqrt((x-xi).^2+(y-yi).^2));
  ind = pointstodo(ind);

end


%---------------------------------------------------------------------
function [dist,ind,pointsinslice] = closestpoint3dp(panel,n,x,y,slice)
%---------------------------------------------------------------------
%find the closest annotation point in 3D for a the object n which is assumed to be a pointgroup

global DATA SET
no = DATA.ViewPanels(panel);

%Set defaults
dist = inf;
ind = [];
pointsinslice = [];

if isempty(n)
  return
end

%Get object
O = SET(no).LevelSet.Object;

if ~O.getvisible(n)
  return
end

%get coordinates of all relevant points
if not(isempty(SET(no).Point))
  %convert xyslice so that we can compare it to the point coordinates
  [r,g,b] = segment3dp.tools('xyz2rgb',SET(no).Point.X,SET(no).Point.Y,SET(no).Point.Z);

  %also get clicked position rgb coordinates
  [rc,gc,bc] = segment3dp.tools('xyz2rgb',x,y,slice);

  %Then find slice match and compare distance in slice plane for closest
  %point
  switch DATA.ViewPanelsType{panel}
    case 'trans3DP'
      pointstodo = find(round(r) == slice);
      [dist,ind] = min(sqrt((gc-g(pointstodo)).^2+(bc-b(pointstodo)).^2));
    case 'sag3DP'
      pointstodo = find(round(g) == slice);
      [dist,ind] = min(sqrt((rc-r(pointstodo)).^2+(bc-b(pointstodo)).^2));
    case 'cor3DP'
      pointstodo = find(round(b) == slice);
      [dist,ind] = min(sqrt((gc-g(pointstodo)).^2+(rc-r(pointstodo)).^2));
  end
else
  dist = inf;
  ind = [];
  return
end

ind = pointstodo(ind);

if nargout == 3
  pointsinslice = pointstodo;
end

%-----------------------------------------------------------
function [dist,ind] = closestinterp(panel,type,x,y,slice,tf)
%-----------------------------------------------------------
%find the closest interpolation point
global DATA SET

no = DATA.ViewPanels(panel);

if nargin<5
  slice = SET(no).CurrentSlice;%clickedslice(panel);
end

if nargin<6
  tf = SET(no).CurrentTimeFrame;
end

if isempty(SET(no).([type,'X'])) || isempty(SET(no).([type,'X']){tf,slice})
  dist = inf;
  ind = [];
  return
end

xi = SET(no).([type,'X']){tf,slice};
yi = SET(no).([type,'Y']){tf,slice};
[dist,ind] = min(sqrt((x-xi).^2+(y-yi).^2));

%------------------------------
function [dist,ind] = closestlaracontour(panel,type,x,y,slice,tf)
%-----------------------
%find the closest interpolation point
global DATA SET

no = DATA.ViewPanels(panel);

if nargin<5
  slice = SET(no).CurrentSlice;
end

if nargin<6
  tf = SET(no).CurrentTimeFrame;
end

if isempty(SET(no).(type)) || ...
    isempty(SET(no).(type).X)
  dist = inf;
  ind = [];
  return
end

xi = SET(no).(type).X(:,tf,slice);
yi = SET(no).(type).Y(:,tf,slice);
[dist,ind] = min(sqrt((x-xi).^2+(y-yi).^2));

%------------------------------
function [dist,ind] = closestlarainterp(panel,type,x,y,slice,tf)
%-----------------------
%find the closest interpolation point
global DATA SET

no = DATA.ViewPanels(panel);

type = type(1:end-6);

if nargin<5
  slice = SET(no).CurrentSlice;%clickedslice(panel);
end

if nargin<6
  tf = SET(no).CurrentTimeFrame;
end

if isempty(SET(no).(type)) || ...
    (isempty(SET(no).(type).InterpX) || ...
    isempty(SET(no).(type).InterpX{tf,slice}))
  dist = inf;
  ind = [];
  return
end

xi = SET(no).(type).InterpX{tf,slice};
yi = SET(no).(type).InterpY{tf,slice};
[dist,ind] = min(sqrt((x-xi).^2+(y-yi).^2));

%--------------------------------------------------------------------------
function [dist,roiind,pointind] = closestroi(panel,xclick,yclick,slice,tf,no)
%--------------------------------------------------------------------------
%The input coordinates are presumed to be in the same coordinate system as
%the stored rois i.e you need to scale and translate your input. pointind
%is optional and is the index of the closest point on the closest roi.

global DATA SET

if nargin<6
  no = DATA.ViewPanels(panel);
end

if not(isempty(SET(no).Flow)) && isfield(SET(no).Flow,'MagnitudeNo') && not(isempty(SET(no).Flow.MagnitudeNo))
  no = SET(no).Flow.MagnitudeNo;
end

if nargin<4 || isempty(slice)
  slice = SET(no).CurrentSlice;%clickedslice(panel);
end

if nargin<5 || isempty(tf)
  tf = SET(no).CurrentTimeFrame;
end

roistodo = find(cellfun(@(x,y) ~isempty(x) && ...
  any(x==slice) && any(y==tf)...
  ,{SET(no).Roi.Z},{SET(no).Roi.T}));

roix = cell(1,length(roistodo));
roiy = cell(1,length(roistodo));

for i = 1:length(roistodo)
  roix{i} = SET(no).Roi(roistodo(i)).X(:,tf);
  roiy{i} = SET(no).Roi(roistodo(i)).Y(:,tf);
end

%Find closest roi in current timeframe and slice
[dist,ind] = cellfun(@(x,y) min(sqrt((xclick-x).^2+(yclick-y).^2)),roix,roiy);
[~,ind_tmp] = min(dist);
dist = dist(ind_tmp);
roiind = roistodo(ind_tmp);

if isempty(dist) || isnan(dist)
  dist = inf;
end

if nargout == 3
  pointind = ind(ind_tmp);
end

%----------------------------------
function ind = findslicewithepi(no,tfs)
%----------------------------------
%Find slices with epicard in any timeframe
global SET NO

if nargin<1
  no = NO;
end
if nargin < 2
  tfs = 1:SET(no).TSize;
end

if isempty(SET(no).EpiX)
  ind = false(SET(no).ZSize,1);
  return;
end

temp = not(isnan(squeeze(SET(no).EpiX(1,tfs,:))));
if length(tfs)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end
else
  ind = temp;
end

%------------------------------------
function ind = findslicewithrvepi(no,tfs)
%------------------------------------
%Find slices with RV epicard in any timeframe
global SET NO

if nargin<1
  no = NO;
end

if nargin < 2
  tfs = 1:SET(no).TSize;
end

if isempty(SET(no).RVEpiX)
  ind = false(SET(no).ZSize,1);
  return;
end

temp = not(isnan(squeeze(SET(no).RVEpiX(1,tfs,:))));
if length(tfs)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end
else
  ind = temp;
end

%--------------------------------------
function ind = findslicewithrvepiall(no)
%--------------------------------------
%Find slices with RV endocard in all timeframes
global SET NO

if nargin<1
  no = NO;
end

if isempty(SET(no).RVEpiX)
  ind = false(SET(no).ZSize,1);
  return;
end

if SET(no).TSize>1
  temp = not(isnan(squeeze(SET(no).RVEpiX(1,:,:))));
  if SET(no).ZSize==1
    ind = all(temp);
  else
    ind = sum(temp,1)==SET(no).TSize;
  end
else
  ind = squeeze(not(isnan(SET(no).RVEpiX(1,1,:))));
end

%--------------------------------
function [x,y] = findlvcenter(no,slices)
%--------------------------------
%Finds the center of the LV, uses the autocrop function as a helper
%function.

%Einar Heiberg

global SET

debugplot = false;

%Call autocrop to get center of heart
[x,y] = autocrop(SET(no).IM);
mx = mean(x); %Calculate mean as this returns a vector
my = mean(y);

%Check if autocrop failed
autocropfailed=false;
if isempty(mx) || isnan(mx)
  mx = SET(no).XSize/2;
  my = SET(no).YSize/2;
  autocropfailed=true;
else
  maxdistfromcenter=140; %#ok<NASGU> %mm %estimate to (half) the radius of a large heart, half is enough but full to have a margin
end

%Takes some midventricular slices
len = 50/2; %50 is 50mm
slicethickness = (SET(no).SliceThickness+SET(no).SliceGap);

%This section could be improved by more actually knowing that
%it is midventricular slices
if nargin>1
  numslices=length(slices);%midslice=slices(round(length(slices)/4));
  midslicestart = slices(max(round(numslices/4),1));%max(round(midslice-len/slicethickness),1);
  midsliceend = slices(max(round(3*numslices/4),1));%min(round(midslice+len/slicethickness),SET(no).ZSize);
else
  midslicestart = max(round(SET(no).ZSize/2-len/slicethickness),1);
  midsliceend = min(round(SET(no).ZSize/2+len/slicethickness),SET(no).ZSize);
end

%Extract only midventricular slices
im = SET(no).IM(:,:,:,midslicestart:midsliceend);

%Average over time and slices
meanim = mean(mean(im,4),3);

%Find most representative timeframe and slice
numtimeframes=size(im,3);
numslices=size(im,4);
difftomeanim=zeros(numtimeframes,numslices);
for tloop=1:numtimeframes
  for zloop=1:numslices
    difftomeanim(tloop,zloop)=sum(sum((meanim-im(:,:,tloop,zloop)).^2));
  end
end

[~,sortindex]=sort(difftomeanim(:));
[tindex,zindex]=ind2sub([numtimeframes,numslices],sortindex);
numindex=min(round(0.1*numslices*numtimeframes),5);
tindex=tindex(1:numindex);
zindex=zindex(1:numindex);

%Take only center for temporary image to find percentile
xd = 30/SET(no).ResolutionX;
yd = 30/SET(no).ResolutionY;

cutim = im(...
  min(max(round(mx-xd):round(mx+xd),1),SET(no).XSize),...
  min(max(round(my-yd):round(my+yd),1),SET(no).YSize),tindex,zindex);

% cutim = smoothimage(cutim,3);
% cutim = cutim/176.8695;
% cutim = cutim(3:(end-2),3:(end-2));

if debugplot
  figure(12); %#ok<UNRCH> 
  imagesc(cutim(:,:,1,1));
  axis image off;
  colorbar
end

%Find percentile
sorted = sort(cutim(:));
initialclassification = 1+double(cutim(:)>(sorted(round(length(sorted)*0.5))));
tol = 0.0005; %default=0.0005
maxiter = 20;
[mu,sigma,~] = emalgorithm(cutim(:),initialclassification,tol,maxiter);

%Sort to ensure largest last
[mu,ind] = sort(mu);
sigma = sigma(ind);

%Mean minus std
thres = mean([mu(1)+2*sigma(1) mu(end)-2*sigma(end)]);

%Threshold over percentile
im=mean(mean(im(:,:,tindex,zindex),4),3);
thresim = im>thres;

%Find regions in thresholded image
[classim,numclass] = bwlabel(thresim);

%Find center of gravity and candidates
[x,y] = ndgrid(1:size(classim,1),1:size(classim,2));
cx = zeros(1,numclass);
cy = zeros(1,numclass);
log = nan(size(cx));

%define locations potential for lv center
potentialcenterlocation=false(SET(no).XSize,SET(no).YSize);
smallheartradius=100;
largeheartradius=140;
mindistfromxedge=round(min(max(SET(no).XSize-largeheartradius/SET(no).ResolutionX,0),smallheartradius/(SET(no).ResolutionX)/2));
mindistfromyedge=round(min(max(SET(no).YSize-largeheartradius/SET(no).ResolutionY,0),smallheartradius/(SET(no).ResolutionY)/2));
potentialcenterlocation(mindistfromxedge+1:end-mindistfromxedge,mindistfromyedge+1:end-mindistfromyedge)=true;

% %Move point 8mm left
% my = my-8/SET(no).ResolutionY; %Not very important...

for loop = 1:numclass

  regionim = (classim==loop); %Extract just that class
  if sum(regionim(:))>100 %Needs to be big enough
    meanregionim = mean(regionim(:).*im(:));
    cx(loop) = mean(regionim(:).*im(:).*x(:))/meanregionim;
    cy(loop) = mean(regionim(:).*im(:).*y(:))/meanregionim;

    %Calc dist to center of heart (center of image if autocrop failed)
    dist = sqrt((cy(loop)-my).^2+(cx(loop)-mx).^2);

    %If center is within the region of potential lv center and distance less than the radius of a large heart or the autocrop has
    %failed
    if potentialcenterlocation(round(cx(loop)),round(cy(loop))) && cy(loop)>my
      if dist<largeheartradius/mean([SET(no).ResolutionX SET(no).ResolutionY]) ||...
          autocropfailed
        log(loop) = dist; %mark as candidate by setting distance
      end
    end
  end
end

%Find closest candidate, i.e smallest number
[~,pos] = sort(log); %sort
pos = pos(1); %take smallest, nan's are found at the end...

%Here one could add some safeguards....

%Extract coordinates
x = cx(pos);
y = cy(pos);

%sum(classim(:)==pos)

SET(no).CenterX = x;
SET(no).CenterY = y;

if debugplot
  figure(6); %#ok<UNRCH> 
  imagesc(thresim.*im);
  hold on;
  plot(my,mx,'k*');
  plot(y,x,'rd');
  hold off;
  axis image off
end

%--------------------------------
function [x,y] = findrvcenter(no,slices)
%--------------------------------
%Finds the center of the RV, currenlty only using center cross definition

global SET NO
if nargin < 1
  no = NO;
end

x = SET(no).CenterX;
y = SET(no).CenterY;

%-------------------------------------------------
function tfs = findframeswithsegmentation(type,no,slices)
%-------------------------------------------------
%Find timeframes in no containing segmentation of type
global SET

if nargin ==2
  slices = 1:SET(no).ZSize;
end

contourx = SET(no).([type,'X']);

if ~isempty(contourx)
  tfs = squeeze(~all(isnan(contourx(1,:,slices)),3));
else
  tfs = zeros(1,SET(no).TSize);
  return
end

%-------------------------------------------------
function tfs = findframeswithinterpolationpoints(type,no,slices)
%-------------------------------------------------
%Find timeframes in no containing segmentation of type
global SET

if nargin ==2
  slices = 1:SET(no).ZSize;
end

typeInterpX = helperfunctions('parsesetfield',SET(no),[type,'Interp'],'X'); %get interp contour

if isempty(typeInterpX)
  tfs = zeros(1,SET(no).TSize);
else
  tfs = ~cellfun(@isempty,typeInterpX(:,slices));
end

if ~isrow(tfs)
  tfs = tfs';
end

if size(tfs,1)>1
  tfs = any(tfs);
end


%---------------------------------------------------------------
function [ind,varargout] = findoutflowtractslices(no,tfs,sectors)
%----------------------------------------------------------------
%Find slices with outflow tract (a wall thickness that is less than 2mm)
%tfs are the timeframes to search in-
%
%Optional output argument is computed wallthickness

global SET NO


if (SET(no).ResolutionX+SET(no).ResolutionY)/2<0.5
  warning('The threshold to find outflow tract (2mm) is species dependent, and might be incorrect for this images.');
end

tresh = 2; %2 mm threshold of wall thickness for slice to be considered outflow tract
ind = false(SET(no).ZSize,1);

if nargin<1
  no = NO;
end
if nargin < 2
  tfs = 1:SET(no).TSize;
end
if nargin < 3
  sectors = 24;
end

if isempty(SET(no).EndoX) && isempty(SET(no).EpiX)
  ind = false(SET(no).ZSize,1);
  return;
end

wallthickness = calcfunctions('calcwallthickness', sectors,no);

if nargout>1
  varargout = {wallthickness};
end

wallthickness = wallthickness(:,:,tfs);
for tf=1:size(wallthickness,3) %loop over timeframes
  s=squeeze(sum((wallthickness(:,:,tf)<=tresh),1)); %check in any slice has a sector with wallthickness less than 2 mm
  ind(find(s~=0))=true; %#ok<FNDSB>
end

%Ensure that no slices in the middle are removed
pos = find(ind,1,'last');
element = false; %start by putting false
for zloop=1:pos %loop to the last found one
  if ind(zloop)
    element = true; %As soon as found one, start putting ones
  end
  ind(zloop) = element;
end


%-------------------------------
function setstack_Callback(type)
%------------------------------
%open interface for user to manually define LV/RV/Flow image stack
global DATA SET
switch type
  case 'flow'
    stri = dprintf('Select image stack for flow report');
  case 'lv'
    stri = dprintf('Select image stack for LV report');
  case 'rv'
    stri = dprintf('Select image stack for RV report');
end
menuitems = cell(1,1);
nn = 1;
stacks = [];

for n = 1:length(SET)
  if isempty(SET(n).Parent) || (not(isempty(SET(n).Parent)) && SET(n).Parent == n)
    menuitems{nn} = sprintf('%d. %s',n,[SET(n).ImageType ' / ' SET(n).ImageViewPlane]);
    stacks = [stacks n]; %#ok<AGROW> 
    nn = nn +1;
  end
end

menuitems{nn} = dprintf('Unselect');
s = mymenu(stri,menuitems);

if ~isempty(s) && s~=0
  if s == length(menuitems)
    %unselect
    switch type
      case 'flow'
        DATA.FlowNO = [];
        DATA.FlowROI = [];
        if isfield(DATA.Handles,'flowstackpushbutton')
          set(DATA.Handles.flowstackpushbutton,'String',dprintf('Set Stack'));
        end
      case 'lv'
        DATA.LVNO = [];
        if isfield(DATA.Handles,'lvstackpushbutton')
          set(DATA.Handles.lvstackpushbutton,'String',dprintf('Set Stack'));
        end
      case 'rv'
        DATA.RVNO = [];
        if isfield(DATA.Handles,'rvstackpushbutton')
          set(DATA.Handles.rvstackpushbutton,'String',dprintf('Set Stack'));
        end
    end
  else
    no = stacks(s);
    switch type
      case 'flow'
        if SET(no).TSize == 1 || isempty(SET(no).Flow) %|| isempty(SET(no).Flow.Result)
          myfailed('Flow image stack need to be time resolved and contain flow analysis.',DATA.GUI.Segment);
          DATA.FlowNO = [];
          DATA.FlowROI = [];
          if isfield(DATA.Handles,'flowstackpushbutton')
            set(DATA.Handles.flowstackpushbutton,'String',dprintf('Set Stack'));
          end
        else
          DATA.FlowNO = no;
          [~,flowroi] = findflowaxisno(no); %identify flow ROI based on ROI names
          DATA.FlowROI = flowroi;
          if isfield(DATA.Handles,'flowstackpushbutton')
            set(DATA.Handles.flowstackpushbutton,'String',dprintf('Stack #%d',no));
          end
        end
      case 'lv'
        %If longaxis no then we need to find the appropriate set. will use
        %Karolinas code
        if any(strcmp(SET(no).ImageViewPlane,{'2CH','3CH','4CH'}))
          LAX_group = findlaxset;
          if all(LAX_group==0)
            DATA.LVNO = no;
            if isfield(DATA.Handles,'lvstackpushbutton')
              set(DATA.Handles.lvstackpushbutton,'String',dprintf('Stack #%d',no));
            end
          else
            LAX_group = LAX_group(LAX_group~=0);
            DATA.LVNO = LAX_group;
            str = [];
            for i = 1:length(LAX_group)
              str = [str, num2str(LAX_group(i)),',']; %#ok<AGROW> 
            end
            str(end) = [];
            if isfield(DATA.Handles,'lvstackpushbutton')
              set(DATA.Handles.lvstackpushbutton,'String',dprintf('Stack #%s',str));
            end
          end
        else
          DATA.LVNO = no;
          if isfield(DATA.Handles,'lvstackpushbutton')
            set(DATA.Handles.lvstackpushbutton,'String',dprintf('Stack #%d',no));
          end
        end
      case 'rv'
        DATA.RVNO = no;
        if isfield(DATA.Handles,'rvstackpushbutton')
          set(DATA.Handles.rvstackpushbutton,'String',dprintf('Stack #%d',no));
        end
    end
  end
end

switch type
  case 'flow'
    segment('updateflow');
  case {'lv','rv'}
    segment('updatevolume',true);
end


%------------------------------------------------------
function saxno = findctsaxwithsegmentation(type)
%-----------------------------------------------------
%Code for finding ct sax image stack with endosegmentation
global SET DATA

switch type
  case 'Endo'
    no = DATA.LVNO;
  case 'RVEndo'
    no = DATA.RVNO;
    %   case 'Epi'
    %      no = DATA.LVNO;
end
if ~isempty(no)
  segx = SET(no).([type,'X']);
else
  segx=[];
end
saxno=[];

if length(no)==1
  if strcmp(SET(no).ImageViewPlane,'Short-axis')&&  not(isempty(segx))
    saxno = no;
  end
end

if isempty(saxno)
  for loop = 1:length(SET)
    if (strcmp(SET(loop).ImageViewPlane,'Short-axis')&&  not(isempty(SET(loop).([type,'X']))))
      saxno= loop;
      break
    end

  end
end



%------------------------------------------------------
function transno = findcttransversalwithendo(type)
%-----------------------------------------------------
%Code for finding ct transversal image stack with endosegmentation

global SET DATA

switch type
  case 'Endo'
    no = DATA.LVNO;
  case 'RVEndo'
    no = DATA.RVNO;
end
transno=[];

if ~isempty(no) && (strcmp(SET(no).ImageViewPlane,'Transversal')) && not(isempty(SET(no).([type,'X'])))
  transno=no;
  return
else
  for loop = 1:length(SET)
    if (strcmp(SET(loop).ImageViewPlane,'Transversal') &&  not(isempty(SET(loop).([type,'X']))))
      transno = loop;
      break
    end
  end
end

%------------------------------------------------------
function transno = findmrtransversalwithendo(heartpart)
%-----------------------------------------------------
%Code for finding mr transversal image stack with endosegmentation
%heart part Lv, RV
global SET DATA

switch heartpart
  case 'LV'
    no = DATA.LVNO;
  case 'RV'
    no = DATA.RVNO;
end
transno=[];
type = 'Endo';

if ~isempty(no) && (strcmp(SET(no).ImageViewPlane,'Transversal')) &&  not(isempty(SET(no).([type,'X'])))
  transno=no;
  return
else
  for loop = 1:length(SET)
    if (strcmp(SET(loop).ImageViewPlane,'Transversal') &&   not(isempty(SET(loop).([type,'X']))))
      transno = loop;
      break
    end
  end
end

%-----------------------------------------------------
function LAX_group = findlaxset(segmentationatedes)
%-----------------------------------------------------
% Code for finding single slice LAX set with endo segmentation. Data is returned
% as [CH2, CH3, CH4], temporary data vector if no chamber can be found
% entry is zero
global SET

LAX_group=[0 0 0]; %% [CH2, CH3, CH4], temporary data vector

if nargin == 0
  segmentationatedes = 0;
end

for loop = 1:length(SET)
  if segmentationatedes
    validchecks = SET(loop).ZSize==1 && isempty(SET(loop).StrainTagging) &&  ...
      not(isempty(SET(loop).EndoX)) && ~all(isnan(SET(loop).EndoX(:,SET(loop).EDT))) && ~all(isnan(SET(loop).EndoX(:,SET(loop).EST)));
  else
    validchecks = SET(loop).ZSize==1 && isempty(SET(loop).StrainTagging) &&  ...
      not(isempty(SET(loop).EndoX));
  end

  if strcmp(SET(loop).ImageViewPlane,'2CH')&& validchecks%SET(loop).ZSize==1 && isempty(SET(loop).StrainTagging) &&  not(isempty(SET(loop).EndoX)))
    if (LAX_group(1)== 0)
      LAX_group(1)=loop;
    end
  end
  if strcmp(SET(loop).ImageViewPlane,'3CH')&& validchecks%SET(loop).ZSize==1 && isempty(SET(loop).StrainTagging) && not(isempty(SET(loop).EndoX)))
    if (LAX_group(2)== 0)
      LAX_group(2)=loop;
    end
  end
  if strcmp(SET(loop).ImageViewPlane,'4CH')&& validchecks%SET(loop).ZSize==1 && isempty(SET(loop).StrainTagging) && not(isempty(SET(loop).EndoX)))
    if (LAX_group(3)== 0)
      LAX_group(3)=loop;
    end
  end
end

%-----------------------------------------------------
function no = findstack(label)
%-----------------------------------------------------
%Finds the stacks with a specific ImageType or ImageViewPlane
%label: string with the ImageType or ImageViewPlane name sought
%no=the stacks with ImageType or ImageViewPlane to the corresponding label

global SET

ind=zeros(1,length(SET));
for loop=1:length(SET)
  if isequal(SET(loop).ImageType,label)
    ind(loop)=1;
  elseif isequal(SET(loop).ImageViewPlane,label)
    ind(loop)=1;
  end
end

no=find(ind==true);

% if isempty(no)
%   myfailed(['Stack named ', label, ' not found.'])
% end

%-----------------------------------------------------
function ind = findclosestannotationpoint(x,y,z)
%-----------------------------------------------------
%Finds the index of the annotation point closest to the coordinates (x,y,z)

global SET NO DATA
no=NO;

if ~contains(DATA.CurrentTheme,'3dp')
  T=SET(no).CurrentTimeFrame;

  points=find(SET(no).Point.T==T); %find annotation points in current timeframe

  if isempty(points)
    myfailed('No annotation points in current time frame.')
    return;
  end
else
  [r,g,b] = segment3dp.tools('xyz2rgb',SET(NO).Point.X,SET(NO).Point.Y, SET(NO).Point.Z);
  switch SET(NO).LevelSet.Pen.Color
    case 'r'
      points= find(round(r)==SET(NO).LevelSet.View.RSlice);
    case 'g'
      points= find(round(g)==SET(NO).LevelSet.View.GSlice);
    case 'b'
      points= find(round(b)==SET(NO).LevelSet.View.BSlice);
  end
end
dist=nan(size(points));
for loop=1:length(points) %loop over points and claculate distance from (x,y,z)
  dist(loop)=sqrt((x-SET(no).Point.X(points(loop)))^2+(y-SET(no).Point.Y(points(loop)))^2+(z-SET(no).Point.Z(points(loop)))^2);
end

ind=points(dist==min(dist));

%--------------------------------------------------------
function z = findcorrespondingslice(thisno,otherno,slice)
%--------------------------------------------------------
%Find slice in otherno stack based on slice in thisno
global SET

%Einar Heiberg

pos = calcfunctions('xyz2rlapfh',thisno,1,1,slice);
pos = calcfunctions('rlapfh2xyz',otherno,pos(1),pos(2),pos(3));
z = pos(3);
z = min(max(round(z),1),SET(otherno).ZSize);

%--------------------------------------------------------
function flowmagno = findflowmagno(flowno)
%--------------------------------------------------------
%Return only the stacks with magnitude information
global SET

if nargin < 1
  [~,~,flowno] = findno;
end

flowmagno = [];
for loop = flowno
  if loop == SET(loop).Flow.MagnitudeNo
    flowmagno = [flowmagno loop]; %#ok<AGROW>
  end
end

%-----------------------------------------------------------------------------------
function [laxgroup,laxnonedt] = findlaxgroup(type,withallexistingcontours,no,doauto)
%-----------------------------------------------------------------------------------
% Find the group of LAX images for segmentation long-axis images set
% type: 'tagging','cine','CT'
arguments
  type = 'cine'
  withallexistingcontours = false;
  no = [];
  doauto = false;
end

global SET DATA

laxnonedt = 0;
% Find image stacks with any long-axis view
laxstacks = findlaxall(type);

if any(laxstacks(:))
  % check if the segmentation exist
  segmentationexist = zeros(size(laxstacks));
  segmentationnonedt = segmentationexist;
  for chloop = 1:3
    for noloop = find(laxstacks(chloop,:))
      edt = SET(noloop).EDT;
      if withallexistingcontours
        laexists = existfunctions('existla',noloop);
        raexists = existfunctions('existra',noloop);
        for sl = 1: SET(noloop).ZSize
          if (~isempty(SET(noloop).EndoX) && ~all(isnan(SET(noloop).EndoX(:,edt,sl)))) ...
              || (~isempty(SET(noloop).RVEndoX) && ~all(isnan(SET(noloop).RVEndoX(:,edt,sl)))) ...
              || (raexists && ~all(isnan(SET(noloop).RA.X(:,edt,sl)))) ...
              || (laexists && ~all(isnan(SET(noloop).LA.X(:,edt,sl))))
            segmentationexist(chloop,noloop) = 1;
            break;
          end
        end
        for sl = 1: SET(noloop).ZSize
          if (~isempty(SET(noloop).EndoX) && all(isnan(SET(noloop).EndoX(:,edt,sl)))) ...
              && (~isempty(SET(noloop).RVEndoX) && all(isnan(SET(noloop).RVEndoX(:,edt,sl)))) ...
              && (raexists && all(isnan(SET(noloop).RA.X(:,edt,sl)))) ...
              && (laexists && all(isnan(SET(noloop).LA.X(:,edt,sl))))
            segmentationnonedt(chloop,noloop) = 1;
          else
            segmentationnonedt(chloop,noloop) = 0;
            break;
          end
        end
      else
        if (~isempty(SET(noloop).EndoX) && ~all(isnan(SET(noloop).EndoX(:,edt))))
          segmentationexist(chloop,noloop) = 1;
        end
      end
    end %end no loop
  end %end ch loop
else
  laxgroup = [];
  disp('No long-axis image found.');
  return
end

if withallexistingcontours 
  segmentationnonedt = laxstacks.*segmentationnonedt;
  laxstacks = laxstacks.*segmentationexist;

  if doauto
    laxgroup = find(sum(laxstacks));
  else
    if all(sum(laxstacks,2) <= 1)
      % exist max one image for each chamber view and there is segmentation
      laxgroup = find(sum(laxstacks));
    else
      laxgroup = chooseLAX('chooselax_helper',laxstacks,no);
    end
  end

  laxnonedt = find(sum(segmentationnonedt));

else
  if DATA.Autoloader
    % for Autoloader get best combination based on number of time frames
    laxgroup = findbestlaxcombination(laxstacks);
  else

    if doauto
      laxgroup = find(sum(laxstacks));
    else
      if all(sum(laxstacks,2) <= 1) && all(sum(segmentationexist,2) == 0)
        % exist max one image for each chamber view and there is no
        % segmentation
        laxgroup = find(sum(laxstacks));
      else
        laxgroup = chooseLAX('chooselax_helper',laxstacks,no);
      end
    end

  end
end

if isempty(laxgroup)
  return
end

if iscolumn(laxgroup)
  laxgroup = laxgroup';
end

% sort them 2CH,3CH,4CH
[~,ind] = sort({SET(laxgroup).ImageViewPlane});
laxgroup = laxgroup(ind);

%-------------------------------------------------------
function laxgroup = findbestlaxcombination(laxstacks)
%-------------------------------------------------------
% Function to find best LAX group combination based on available number of time
% frames

global SET

% Get indices of each type
laxidx{1} = find(laxstacks(1,:) == 1); % 2 chamber
laxidx{2} = find(laxstacks(2,:) == 1); % 3 chamber
laxidx{3} = find(laxstacks(3,:) == 1); % 4 chamber

bestcombo = [];
bestscore = -inf;
maxtfdiff = 4; % Time frame sizes can maximally differ by this value
maxtincr = 1;

% Check all 3-stack combinations
for ch2loop = 1:length(laxidx{1})
  for ch3loop = 1:length(laxidx{2})
    for ch4loop = 1:length(laxidx{3})
      combo = [laxidx{1}(ch2loop), laxidx{2}(ch3loop), laxidx{3}(ch4loop)];
      tsizes = [SET(combo(1)).TSize, SET(combo(2)).TSize, SET(combo(3)).TSize];
      tincr = [SET(combo(1)).TIncr, SET(combo(2)).TIncr, SET(combo(3)).TIncr];
      if max(tsizes) - min(tsizes) <= maxtfdiff && all(tincr < maxtincr)
        avgtsize = mean(tsizes);
        % '=' ensures the combination with the larger stack numbers is
        % taken, since probably that these stacks were re-acquired
        % because the earlier stacks were of lower quality
        if avgtsize >= bestscore 
          bestscore = avgtsize;
          bestcombo = combo;
        end
      end
    end
  end
end

% If no valid 3-stack combo, try 2-stack combos
if isempty(bestcombo)
  for type1 = 1:3 % Loop over all unique pairs of view types (2ch, 3ch, 4ch)
    for type2 = type1+1:3 % Each pair is only considered once
      for type1loop = 1:length(laxidx{type1}) % Loop over all available stacks for the first type (t1)
        for type2loop = 1:length(laxidx{type2}) % Loop over all available stacks for the second type (t2)
          combo = [laxidx{type1}(type1loop), laxidx{type2}(type2loop)];
          tsizes = [SET(combo(1)).TSize, SET(combo(2)).TSize];
          tincr = [SET(combo(1)).TIncr, SET(combo(2)).TIncr];
          if max(tsizes) - min(tsizes) <= maxtfdiff && all(tincr < maxtincr)
            avgtsize = mean(tsizes);
            if avgtsize >= bestscore
              bestscore = avgtsize;
              bestcombo = combo;
            end
          end
        end
      end
    end
  end
end

% If still empty, get stack with most time frames
if isempty(bestcombo)
  laxnos = cell2mat(laxidx);  
  tf = [SET(laxnos).TSize];
  ind = tf == max(tf);
  bestcombo = laxnos(ind);
end
laxgroup = bestcombo;

%-------------------------------------------------------
function nos = findsaxnowithheartpart(heartpart,inputnos)
%-------------------------------------------------------
%Return stacks number where the specified heart part exists in EDT. If
%inputnos is given as an input argument, then the function looks for hte
%heart part only in the provided stacks.
arguments
  heartpart
  inputnos = [];
end
nos = [];

functiontoexecute = ['exist' lower(heartpart) 'edt'];

if ~isempty(inputnos) %look only into the provided stacks
  numnos = numel(inputnos);
  %find stacks where heartpart exists
  nos = nan(1,numnos);
  for loop = 1:numnos
    if existfunctions(functiontoexecute,inputnos(loop))
      nos(loop) = inputnos(loop);
    end
  end
end
nos(isnan(nos)) = [];

%-------------------------------------------------------
function nos = findlaxnowithheartpart(heartpart,inputnos,timeframe)
%-------------------------------------------------------
%Return stacks number where the specified heart part exists in EDT or EST. If
%inputnos is given as an input argument, then the function looks for hte
%heart part only in the provided stacks.
arguments
  heartpart
  inputnos = [];
  timeframe {mustBeMember(timeframe,{'edt','est',''})} = 'edt'
end
global SET
nos = [];

functiontoexecute = ['exist' lower(heartpart) timeframe];

if ~isempty(inputnos) %look only into the provided stacks
  numnos = numel(inputnos);
  %find stacks where heartpart exists
  nos = nan(1,numnos);
  for loop = 1:numnos
    if existfunctions(functiontoexecute,inputnos(loop))
      nos(loop) = inputnos(loop);
    end
  end

else
  % Find image stacks with any long-axis view
  chamberexist = findlaxall;
  if any(chamberexist(:))
    %find stacks where heartpart exists
    nos = nan(1,3);
    for chloop = 1:3
      for no = find(chamberexist(chloop,:))
        if existfunctions(functiontoexecute,no)
          nos(chloop) = no;
        end
      end %end no loop
    end %end ch loop
  end
end

nos(isnan(nos)) = [];

% additionally exclude LA from 3CH
if strcmpi(heartpart,'LA')
  ind = strcmpi({SET(nos).ImageViewPlane},'3CH');
  nos = nos(~ind);
end
% additionally include RA from 4CH only
if strcmpi(heartpart,'RA')
  ind = strcmpi({SET(nos).ImageViewPlane},'4CH');
  nos = nos(ind);
end

%-------------------------------------------------------
function nos = findlaxnowithlv(inputnos)
%-------------------------------------------------------
%Return stacks number where LV exists in EDT
arguments
  inputnos = [];
end
nos = findlaxnowithheartpart('LV',inputnos);

%-------------------------------------------------------
function nos = findlaxnowithrv(inputnos)
%-------------------------------------------------------
%Return stacks number where RV exists in EDT
arguments
  inputnos = [];
end
nos = findlaxnowithheartpart('RV',inputnos);

%-------------------------------------------------------
function nos = findlaxnowithla(inputnos)
%-------------------------------------------------------
%Return stacks number where LA exists in EDT
arguments
  inputnos = [];
end
nos = findlaxnowithheartpart('LA',inputnos);

%-------------------------------------------------------
function nos = findlaxnowithlaest(inputnos)
%-------------------------------------------------------
%Return stacks number where LA exists in EDT
arguments
  inputnos = [];
end
nos = findlaxnowithheartpart('LA',inputnos,'est');

%-------------------------------------------------------
function nos = findlaxnowithra(inputnos)
%-------------------------------------------------------
%Return stacks number where RA exists in EDT
arguments
  inputnos = [];
end
nos = findlaxnowithheartpart('RA',inputnos);

%-------------------------------------------------------
function nos = findlaxnowithraest(inputnos)
%-------------------------------------------------------
%Return stacks number where RA exists in EST
arguments
  inputnos = [];
end
nos = findlaxnowithheartpart('RA',inputnos,'est');


%-------------------------------------------------------
function chamberexist = findlaxall(type)
%-------------------------------------------------------
% Find image stacks with any long-axis view
arguments
  type = 'cine'
end
global SET

typestring = getimagetypestring_helper(type);

chamberexist(1,:) = strcmp('2CH',{SET.ImageViewPlane}).* strcmp(typestring,{SET.ImageType});
chamberexist(2,:) = strcmp('3CH',{SET.ImageViewPlane}).* strcmp(typestring,{SET.ImageType});
chamberexist(3,:) = strcmp('4CH',{SET.ImageViewPlane}).* strcmp(typestring,{SET.ImageType});

%-------------------------------------------------------
function nos = findnoXch(viewplane,type)
%-------------------------------------------------------
% Find image stacks with specified image view plane
arguments
  viewplane {mustBeMember(viewplane,{'2CH','3CH','4CH'})}
  type {mustBeMember(type,{'tagging','cine','ct'})} = 'cine'
end

global SET

typestring = getimagetypestring_helper(type);

nos = find(strcmp(viewplane,{SET.ImageViewPlane}).* strcmp(typestring,{SET.ImageType}));

%-------------------------------------------------------
function typestring = getimagetypestring_helper(type)
%-------------------------------------------------------
%Helper function to get the full and correct string for image type
switch lower(type)
  case 'tagging'
    typestring = 'Strain from tagging';
  case 'cine'
    typestring = 'Cine';
  case 'ct'
    typestring = 'General';
end

%--------------------------------------------------------------------------
function [stacks,roinum] = findstackswithroiflow(roilabel,flowmagnitudeno)
%--------------------------------------------------------------------------
% find all magnitude stacks with a ROI labelled with the specified label 
% and the corresponding ROI numbers
global SET

%Init input arguments
if nargin < 1
  flowmagnitudeno = findflowmagno;
end

%Init output arguments
stacks = [];
roinum = [];

%Loop over all magnitude stacks and find ROIs
for no = flowmagnitudeno
  if SET(no).RoiN > 0
    rois = {};
    %Get all ROI's labels
    for roin = 1:SET(no).RoiN
      rois{roin} = SET(no).Roi(roin).Name; %#ok<AGROW>
    end
    %Find those labelled with the specified ROI
    paf = find(strcmp(roilabel,rois), 1);
    if ~isempty(paf)
      stacks = [stacks,no]; %#ok<AGROW> we don't know how many stacks contain the specified ROI
      roinum = [roinum,paf(1)]; %#ok<AGROW>
    end      
  end
end

%--------------------------------------------------------------------------
function [stack,roinum,success] = findonestackwithroiflow(roilabel,flowmagnitudeno)
%--------------------------------------------------------------------------
% find one stack with a ROI labelled with the specified label and the 
% corresponding ROI numbers. If multiple stacks are found then the user is
% prompted to choose one unique stack.

%Init input arguments
if nargin < 1
  flowmagnitudeno = findflowmagno;
end

%Init output arguments
success = false;
stack = [];
roinum = [];

%Find stacks
[stacks,roinums] = findstackswithroiflow(roilabel,flowmagnitudeno);

%Multiple stacks found
if numel(stacks) > 1
  str = dprintf('Found multiple stacks with a ROI labelled:');
  str = sprintf('%s "%s"',str,roilabel);

  stackcell = cell(1,numel(stacks));
  for no = 1:numel(stacks)
    stackcell{no} = sprintf('%s %d',dprintf('Stack'),stacks(no));
  end
  ind = mymenu(str,stackcell);
  if ind ~= 0
    stack = stacks(ind);
    roinum = roinums(ind);
    success = true;
  end
else
  %Unique stack found
  stack = stacks;
  roinum = roinums;
  success = true;
end

%-----------------------------------------------------------------
function [cinesaxlvno,cinesaxrvno] = findcinesaxreport
%-----------------------------------------------------------------
%Find shortaxis cine to use for report/export functionnalities
global DATA

%Find stacks with LV and RV segmentations
returnboth = true;
cinesaxnos = findcineshortaxisno(returnboth); %Two image stacks returned (LV and RV)

%If LVNO and RVNO exist, then use them. Otherwise use the stacks found before
if length(DATA.LVNO) == 1 && ~strcmp(DATA.ProgramName,'Segment')
  cinesaxlvno = DATA.LVNO;
else
  if ~isempty(cinesaxnos)
    cinesaxlvno = cinesaxnos(1);
  else
    cinesaxlvno = [];
  end
end

if length(DATA.RVNO) == 1 && ~strcmp(DATA.ProgramName,'Segment')
  cinesaxrvno = DATA.RVNO;
else
  if ~isempty(cinesaxnos)
    cinesaxrvno = cinesaxnos(2);
  else
    cinesaxrvno = [];
  end
end

%-----------------------------------------------------------------
function thisslice = findcorrespondingslicebasedonendocircumf(refslice,refno,thisno)
%-----------------------------------------------------------------
%Find corresponding slice based on circumference of endo
global SET
thisslice = [];
try
  timeframe = 1;
  x = squeeze([SET(refno).EndoX(:,timeframe, refslice); SET(refno).EndoX(1,timeframe, refslice)]);
  y = squeeze([SET(refno).EndoY(:,timeframe, refslice); SET(refno).EndoY(1,timeframe, refslice)]);
  circumferenceref = sum(sqrt((diff(x)*(SET(refno).ResolutionX)).^2+(diff(y)*SET(refno).ResolutionY).^2));

  slices = find(findfunctions('findslicewithendo',thisno,timeframe)); %find slices with endo segmentation
  for slice = slices %loop through all slices with endo
    x = [SET(thisno).EndoX(:,timeframe, slice); SET(thisno).EndoX(1,timeframe, slice)];
    y = [SET(thisno).EndoY(:,timeframe, slice); SET(thisno).EndoY(1,timeframe, slice)];
    circumferencethisno(slice) = sum(sqrt((diff(x)*(SET(thisno).ResolutionX)).^2+(diff(y)*SET(thisno).ResolutionY).^2)); %#ok<AGROW> 
  end
  diffcircumf = abs(circumferencethisno-circumferenceref);
  m = min(diffcircumf);
  thisslice = find(diffcircumf == m);
catch me
  mydispexception(me)
end

%-----------------------------------------------------------------
function nos = findall4dflowstacks
%-----------------------------------------------------------------
% find alls stacks with 4D flow
global SET

numstacks = length(SET);
nos = [];
for loop = 1:numstacks
  %if isequal(lower(SET(loop).ImageType),'cine') && ... %Removed this
  %requirement. EiH
  if isequal(length(SET(loop).Children),3) && ...
      ~isequal(length(SET(loop).VENC),0)

    nos = [nos,loop]; %#ok<AGROW>
  end
end

%----------------------
function nos = findecvno
%----------------------
% find ECV stacks
global SET
nos = find(cellfun(@(x) isfield(x,'reportecvim')&& not(isempty(x)),{SET.ECV}));

%--------------------------------
function nos = findrelevantstacks
%--------------------------------
global DATA SET

if strcmp(DATA.CurrentTheme, 'txmap')
  % Handle 'txmap' directly
  t1stacks = findallt1stacks;
  t2stacks = findallt2stacks;
  t2starstacks = findallt2starstacks;
  txmapstacks = findtxmapnosviaimagetype;
  % find stacks based on image type
  nos = unique([t1stacks, t2stacks, t2starstacks,txmapstacks]);
elseif strcmp(DATA.CurrentTheme, 'strain')
  nos = findallstrainstacks;
else
  [cineno,scarno,flowno,~,marno,~] = findno;
  switch DATA.CurrentTheme
    case 'function'
      nos = cineno;
    case 'roi'
      nos = flowno;
    case 'scar' %viability tab      
      nos = [scarno,marno];
      if ~isempty(DATA.LVNO)
        lvno = DATA.LVNO;
      else
        multiplenos = false;
        minscore = 28;
        lvno = findcineshortaxisno(multiplenos, minscore);
      end 
      nos = [nos,lvno];
    otherwise
      nos = 1:length(SET);
  end
end
if isempty(nos)
  nos = 1:length(SET);
  DATA.FoundRelevantStacks = false;
else
  % add all currently visible stacks
  visiblenos = nonzeros(DATA.ViewPanels);
  nos = union(visiblenos,nos);
  DATA.FoundRelevantStacks = true;
end
nos = sort(nos);

%---------------------------------------
function strainnos = findallstrainstacks
%---------------------------------------
% function to find all SAX stacks that are relevant for strain
global SET
persistent includetagging
numstacks = length(SET);
strainnos = zeros(1,numstacks);
imagetypes = {'Feature tracking','General','Cine','Strain FFE','Strain TFE'};
viewplanes = {'Short-axis','2CH','3CH','4CH'};
if isempty(includetagging)
  [~,includetagging] = segment('checkstrainlicense'); % check license for Strain 1st generation
end
if includetagging
  imagetypes{end+1} = 'Strain from tagging';
end

% check if strain mitt or stain tagging exist
for noloop = 1:numstacks
  if ~isempty(SET(noloop).StrainMitt) || ~isempty(SET(noloop).StrainTagging)
    strainnos(noloop) = noloop;
  elseif matches(SET(noloop).ImageViewPlane,viewplanes,"IgnoreCase",true) &&...
      matches(SET(noloop).ImageType,imagetypes,"IgnoreCase",true)...
      && SET(noloop).TSize > 8
    strainnos(noloop) = noloop;
  end
end
strainnos = nonzeros(strainnos);

%--------------------------------------
function nos = findtxmapnosviaimagetype
%--------------------------------------
% find tx-maps stacks using imge type only
global SET

imagetypes = {'T1 map','T2 map','T2* map','ECV'};

nos = find(contains({SET.ImageType},imagetypes,"IgnoreCase",true));

%---------------------------------------------------------------------------
function [foundstacks,noaorta,nopulmo,roinumaorta,roinumpulmo] = findqpqsno
%---------------------------------------------------------------------------
% Find stacks with Pulmonary artery and Aortic ascending flow
foundstacks = false(1,2);
noaorta = [];
nopulmo = [];
roinumaorta = [];
roinumpulmo = [];
%Find magnitude image on flow
[~,~,flowno] = findfunctions('findno');
flowmagnitudeno = findfunctions('findflowmagno',flowno);
[aostacks,aoroinum] = findfunctions('findstackswithroiflow','Aortic ascending flow',flowmagnitudeno);
[pulstacks,pulroinum] = findfunctions('findstackswithroiflow','Pulmonary artery',flowmagnitudeno);
if ~isempty(aostacks)
  % take the first ones
  idx = 1;
  noaorta = aostacks(idx);  
  roinumaorta = aoroinum(idx);
  foundstacks(1) = true;
end
if ~isempty(pulstacks)
  % take the first ones
  idx = 1;  
  nopulmo = pulstacks(idx);
  roinumpulmo = pulroinum(idx);
  foundstacks(2) = true;
end

%----------------------------------------------------------
function indroi = findallroiinstackexcept(no,excludednames)
%----------------------------------------------------------
% find index for all ROIs except the ROIs with provided names
global SET

roinames = {SET(no).Roi.Name};

indroi = find(~matches(roinames, excludednames) & strlength(roinames) > 0);


%-----------------------------
function nos = findgeneralscar
%-----------------------------
% Find stacks that have more than one slice, only one time frames and
% imageview plane is short-axis
global SET

tsize = [SET.TSize];
zsize = [SET.ZSize];
type = {SET.ImageViewPlane};

nos = find(tsize==1 & zsize > 1 & contains(type,'short-axis','IgnoreCase',true));


