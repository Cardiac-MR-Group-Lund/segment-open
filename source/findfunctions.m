function varargout = findfunctions(varargin)
% FINDFUNCTIONS
% Functions for finding image stacks or slices

% Moved out from segment_main by Nisse Lundahl

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
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
      ~isempty(SET(loop).Strain) || ...
      (~isempty(SET(loop).Flow) && isempty(SET(loop).Flow.PhaseNo) && ~isempty(SET(loop).Flow.PhaseX) && ~isempty(SET(loop).Flow.PhaseY) && isequal(SET(loop).Flow.MagnitudeNo,loop))
     strainno = [strainno loop]; %#ok<AGROW>
  end
  if isequal(SET(loop).ImageType,'Flow (magnitude)')
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
function [saxno, success] = findsaxno %#ok<DEFNU>
%-------------------------------------
% find current stack with short axis image
global SET NO DATA
cineno = findfunctions('findno');
saxno = find(strcmp({SET.ImageViewPlane},'Short-axis'));
zno = find([SET.ZSize] > 1);
saxno = intersect(cineno,union(saxno,zno));
success = 1;
if isempty(saxno)
  if ismember(NO,zno) && ~DATA.Silent &&...
      yesno('Could not find image stack defined as Short-axis Cine. Use current stack?')
    saxno = NO;
    SET(saxno).ImageViewPlane = 'Short-axis';
    SET(saxno).ImageType = 'Cine';
    for p = find(DATA.ViewPanels==saxno)
      drawfunctions('drawpanel',p);
    end
  elseif ismember(NO,zno) && DATA.Silent
    saxno = NO;
    SET(saxno).ImageViewPlane = 'Short-axis';
    SET(saxno).ImageType = 'Cine';
  else
    if ~DATA.Silent
      myfailed('Could not find short-axis image stack');
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
    m = mymenu('Do segmentation on short-axis stack?', nos{:});
  else
    m = 1;
  end
  if m
    saxno = saxno(m);
  else
    success = 0;
    return;
  end
end

%---------------------------------------------
function cineshortaxisno = findcineshortaxisno(multiple, minscore) %#ok<DEFNU>
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
end

%Find if slices equals max
%Choose the best candidate in order of priority. Each row in A is one
%image stack, multiply with priority vector. Take largest sum.
A = [iscine' isshortaxis' haslvseg' numberofsliceswithLV' hasrvseg' numberofsliceswithRV' hasmaxslices'];
score = A*[16;12;4;1;2;1;1]; %most important in priority cine, shortaxis, lvseg, numberofslicesLV, rvseg, numberofslicesLV, maxslices

if multiple
  %find lv and rv stacks that is cine and short-axis
  
  templvnos = find(haslvseg);
  temprvnos = find(hasrvseg);
  if ~isempty(templvnos)
    [~,lvind] = max(score(templvnos));
    lvno = cineno(templvnos(lvind));
  else
    [~,lvnoind] = max(score);
    lvno = cineno(lvnoind);
  end
  if ~isempty(temprvnos)
    [~,rvind] = max(score(temprvnos));
    rvno = cineno(temprvnos(rvind));
  else
    [~,rvnoind] = max(score);
    rvno = cineno(rvnoind);
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
function shortaxisno = findctshortaxisno(multiple) %#ok<DEFNU>
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
function shortaxisno = findspectshortaxisno(multiple) %#ok<DEFNU>
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
function scarno = findscarshortaxisno %#ok<DEFNU>
%------------------------------------------
%Find only one scar shortaxis stack
global SET

[~,scarno] = findno;

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
  end
  [maxval,maxind] = max(scar2use);
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
    s = scarno(1);
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
  elseif length(scarno)>1
    disp('Detected multiple scar data. Taking first image stack (arbitrary decision)');
    scarno = scarno(1);
  end
end

%------------------------------------------
function marno = findmarshortaxisno %#ok<DEFNU>
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
      magno = [magno SET(flowno(noloop)).Flow.MagnitudeNo];
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
  [maxval,maxind] = max(flow2use);
  allmax = find(flow2use == maxval);
  flowno = flowno(allmax);
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
    [maxval,maxind] = max(maxpoints);
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
    [maxval,maxind] = max(points);
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


%--------------------------------------
function ind = findslicewithscarall(no) %#ok<DEFNU>
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
function ind = findslicewithmarall(no) %#ok<DEFNU>
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
function ind = findslicewithendoall(no) %#ok<DEFNU>
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
function ind = findslicewithepiall(no) %#ok<DEFNU>
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
function ind = findslicewithendo(no,tfs) %#ok<DEFNU>
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
function ind = findslicewithrvendo(no,tfs) %#ok<DEFNU>
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

%--------------------------------------
function ind = findslicewithrvendoall(no) %#ok<DEFNU>
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
%-----------------------
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
function [type,objind] = closestobject(panel,x,y,slice,tf) %#ok<DEFNU>
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);

if nargin<4
  slice = SET(no).CurrentSlice;%clickedslice(panel);
end

if nargin<5
  tf = SET(no).CurrentTimeFrame;
end

minendo = dist2contour('Endo',no,x,y,slice,tf);
minepi = dist2contour('Epi',no,x,y,slice,tf);
minrvendo = dist2contour('RVEndo',no,x,y,slice,tf);
minrvepi = dist2contour('RVEpi',no,x,y,slice,tf);
if strcmp(DATA.ProgramName, 'Segment')
  mincenter = min(sqrt((SET(no).CenterX-x).^2+(SET(no).CenterY-y).^2));
else
  mincenter = Inf;
end

%this is the code to find closest measure
[measureind, ~,~,~,~,minmeasure] = closestmeasure(panel,y,x); 

%Find closest roi in current timeframe and slice
[minroi,roiind] = closestroi(panel,x,y,slice,tf);

%Find closest point in current timeframe and slice
[minpoint,pointind] = closestpoint(panel,x,y,slice,tf);

%Find closest interp
[minendointerp,endointerpind] = closestinterp(panel,'EndoInterp',x,y,slice,tf);
[minepiinterp,epiinterpind] = closestinterp(panel,'EpiInterp',x,y,slice,tf);
[minrvendointerp,rvendointerpind] = closestinterp(panel,'RVEndoInterp',x,y,slice,tf);
[minrvepiinterp,rvepiinterpind] = closestinterp(panel,'RVEpiInterp',x,y,slice,tf);

%summarize which is closest
dists = [minmeasure,minroi,minendo,minepi,minrvendo,minrvepi,...
    minendointerp,minepiinterp,minrvendointerp,minrvepiinterp,minpoint,mincenter];

[mindist,typeind] = min(dists);

types = {'Measure','Roi','Endo','Epi','RVEndo','RVEpi','EndoInterp',...
    'EpiInterp','RVEndoInterp','RVEpiInterp','Point', 'Center'};

type = types{typeind};

%If selected type is contour then check if Interp is within intervall and
%set Interp version of contour as type
if typeind>2 && typeind<7 && dists(typeind+4)<=DATA.Pref.ContourAdjustDistance %preference dist to object
  type = types{typeind+4};
end

%if smallest distance is more than contour adjust distance then return
%image type for context menu
if mindist >DATA.Pref.ContourAdjustDistance
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
  otherwise %objind is irrelevant
    objind = 0;
end


%------------------------------
function [dist,ind] = closestpoint(panel,x,y,slice,tf) 
%-----------------------
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


%------------------------------
function [dist,ind,pointsinslice] = closestpoint3dp(panel,x,y,slice) %#ok<DEFNU>
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);

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

%------------------------------
function [dist,ind] = closestinterp(panel,type,x,y,slice,tf) 
%-----------------------
global DATA SET

no = DATA.ViewPanels(panel);

if nargin<5
    slice = SET(no).CurrentSlice;%clickedslice(panel);
end

if nargin<6
    tf = SET(no).CurrentTimeFrame;
end

if isempty(SET(no).([type,'X'])) ||...
        isempty(SET(no).([type,'X']){tf,slice})
    dist = inf;
    ind = [];
    return
end

xi = SET(no).([type,'X']){tf,slice};
yi = SET(no).([type,'Y']){tf,slice};
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
function ind = findslicewithepi(no,tfs) %#ok<DEFNU>
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
function ind = findslicewithrvepi(no,tfs) %#ok<DEFNU>
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
function ind = findslicewithrvepiall(no) %#ok<DEFNU>
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
function [x,y] = findlvcenter(no,slices) %#ok<DEFNU>
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
	maxdistfromcenter=140;%mm %estimate to (half) the radius of a large heart, half is enough but full to have a margin
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
  figure(12);
  imagesc(cutim(:,:,1,1));
  axis image off;
  colorbar
end

%Find percentile
sorted = sort(cutim(:));
initialclassification = 1+double(cutim(:)>(sorted(round(length(sorted)*0.5))));
tol = 0.0005; %default=0.0005
maxiter = 20;
[mu,sigma,alpha] = emalgorithm(cutim(:),initialclassification,tol,maxiter);

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
  figure(6);
  imagesc(thresim.*im);
  hold on;
  plot(my,mx,'k*');
  plot(y,x,'rd');
  hold off;  
  axis image off
end

%--------------------------------
function [x,y] = findrvcenter(no,slices) %#ok<DEFNU>
%--------------------------------
%Finds the center of the RV, currenlty only using center cross definition

global SET NO
if nargin < 1
  no = NO;
end

x = SET(no).CenterX;
y = SET(no).CenterY;

%-------------------------------------------------
function tfs = findframeswithsegmentation(type,no,slices) %#ok<DEFNU>
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
function tfs = findframeswithinterpolationpoints(type,no,slices) %#ok<DEFNU>
%-------------------------------------------------
%Find timeframes in no containing segmentation of type
global SET

if nargin ==2
    slices = 1:SET(no).ZSize;
end

if isempty(SET(no).([type,'InterpX']))
    tfs = zeros(1,SET(no).TSize);
else
    tfs = ~cellfun(@isempty,SET(no).([type,'InterpX'])(:,slices));
end

if ~isrow(tfs)
    tfs = tfs';
end

if size(tfs,1)>1
    tfs = any(tfs);
end

 
%---------------------------------------
function [ind] = findoutflowtractslices(no,tfs) %#ok<DEFNU>
%---------------------------------------
%Find slices with outflow tract (a wall thickness that is less than 2mm)
%tfs are the timeframes to search in
global SET NO

if (SET(no).ResolutionX+SET(no).ResolutionY)/2<0.5
  warning('The threshold to find outflow tract (2mm) is species dependent, and might be incorrect for this images.');
end

sectors=24; %we use a 24 sector model
tresh=2; %2 mm threshold of wall thickness for slice to be considered outflow tract
ind=false(SET(no).ZSize,1);

if nargin<1
  no = NO;
end
if nargin < 2
  tfs = 1:SET(no).TSize;
end

if isempty(SET(no).EndoX) && isempty(SET(no).EpiX)
  ind = false(SET(no).ZSize,1);
  return;
end

wallthickness = calcfunctions('calcwallthickness', sectors,no);
wallthickness=wallthickness(:,:,tfs);
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
function setstack_Callback(type) %#ok<DEFNU>
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
    stacks = [stacks n];
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
                str = [str, num2str(LAX_group(i)),','];
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
global SET NO DATA

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

global SET NO DATA

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

if ~strcmp(DATA.CurrentTheme,'3dp')
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

ind=points(find(dist==min(dist)));

%--------------------------------------------------------
function z = findcorrespondingslice(thisno,otherno,slice) %#ok<DEFNU>
%--------------------------------------------------------
%Find slice in otherno stack based on slice in thisno
global SET

%Einar Heiberg

pos = calcfunctions('xyz2rlapfh',thisno,1,1,slice);
pos = calcfunctions('rlapfh2xyz',otherno,pos(1),pos(2),pos(3));
z = pos(3);
z = min(max(round(z),1),SET(otherno).ZSize);
