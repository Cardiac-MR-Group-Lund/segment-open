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
end;

if ~isfield(SET(1),'ImageType')
  SET(1).ImageType = '';
end;
%Check with image types first!
for loop=1:length(SET)
  if isequal(lower(SET(loop).ImageType),'cine')
    cineno = [cineno loop]; %#ok<AGROW>
  end;
  if isequal(lower(SET(loop).ImageType),'late enhancement') && ~isempty(SET(loop).Scar)
    scarno = [scarno loop]; %#ok<AGROW>
  elseif isempty(scarno) && (isequal(lower(SET(loop).ImageType),'late enhancement') || ~isempty(SET(loop).Scar) )
    scarno = loop;
  end;
  if isequal(lower(SET(loop).ImageType),'strain ffe') || isequal(lower(SET(loop).ImageType),'strain tfe') || ...
      ~isempty(SET(loop).Strain) || ...
      (~isempty(SET(loop).Flow) && isempty(SET(loop).Flow.PhaseNo) && ~isempty(SET(loop).Flow.PhaseX) && ~isempty(SET(loop).Flow.PhaseY) && isequal(SET(loop).Flow.MagnitudeNo,loop))
     strainno = [strainno loop]; %#ok<AGROW>
  end;
  if isequal(SET(loop).ImageType,'Flow (magnitude)')
    flowno = [flowno loop]; %#ok<AGROW>    
  end;
  if isequal(lower(SET(loop).ImageType),'perfusion rest') && ~isempty(SET(loop).MaR)
    marno = [marno loop]; %#ok<AGROW>
  end;
  if isequal(lower(SET(loop).ImageType),'perfusion stress') && ~isempty(SET(loop).MaR)
    marno = [marno loop]; %#ok<AGROW>
  end;
  
%   if nargout>3
%     if isequal(SET(loop).ImageType,'Stress baseline');
%       varargout{1} = [varargout{1} loop];
%     end;
%   end;
  if nargout>4
    if isequal(lower(SET(loop).ImageType),'stress')
      varargout{2} = [varargout{2} loop];
    end;    
  end;
end;
  
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
      end; %not mar
      
      %If neither scar,flow check if short axes?
      %Note it may be both mar and cine!
      if ~isempty(SET(loop).EndoX)||~isempty(SET(loop).EpiX)||~isempty(SET(loop).RVEndoX)||~isempty(SET(loop).RVEpiX)
        if (SET(loop).TSize>1)&&...
            (sum(sum(not(isnan([SET(loop).EndoX(:); SET(loop).EpiX(:); SET(loop).RVEndoX(:); SET(loop).RVEpiX(:)]))))>0)&&...
            isempty(find(cineno==loop,1))

%             (sum(sum(not(isnan([SET(loop).EndoX(:) SET(loop).EpiX(:) SET(loop).RVEndoX(:) SET(loop).RVEpiX(:)]))))>0)&&...
%             isempty(find(cineno==loop,1))
          cineno = [cineno loop]; %#ok<AGROW>
        end;
      end;
      
    end; %not flow
  end; %not scar
end;

%---------------------------------------------
function cineshortaxisno = findcineshortaxisno(multiple)
%---------------------------------------------
%Find one [LV] or two [LV and RV] cine short axis stack
global SET

if nargin < 1
  multiple = false; %default, only return one no, true = return all cine sax
end

%cineshortaxisno = [];
[cineno] = findno;

cineno = cineno(:)'; %ensure row vector

%If isempty then just return something.
if isempty(cineno)
  cineshortaxisno = [];
  return;
end;

%Find maximal number of slices
maxslices = max(cat(1,SET(:).ZSize));

%Extract properties.
isshortaxis = false(size(cineno));
iscine = false(size(cineno));
haslvseg = false(size(cineno));
hasrvseg = false(size(cineno));
hasmaxslices = ones(size(cineno));

for loop=1:length(cineno)
  iscine(loop) = isequal(lower(SET(cineno(loop)).ImageType),'cine') || (SET(cineno(loop)).TSize>=20 && ~strncmp(SET(cineno(loop)).ImageType,'Perfusion',9));
  isshortaxis(loop) = isequal(lower(SET(cineno(loop)).ImageViewPlane),'short-axis');  
  haslvseg(loop) = (~isempty(SET(cineno(loop)).EndoX)|~isempty(SET(cineno(loop)).EpiX));
  hasrvseg(loop) = (~isempty(SET(cineno(loop)).RVEndoX)|~isempty(SET(cineno(loop)).RVEpiX));
  hasmaxslices(loop) = isequal(SET(cineno(loop)).ZSize,maxslices);
end;

%Find if slices equals max
%Choose the best candidate in order of priority. Each row in A is one
%image stack, multiply with priority vector. Take largest sum.
A = [iscine' haslvseg' isshortaxis'  hasrvseg' hasmaxslices'];
 
if multiple
  %find lv and rv stacks that is cine and short-axis
  score = A*[16;8;4;2;1]; %most important in priority cine, lvseg, shortaxis, rvseg, maxslices
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
  score = A*[16;8;4;2;1]; %most important in priority cine, lvseg, shortaxis, rvseg, maxslices
  [~,ind] = max(score);
  cineshortaxisno = cineno(ind);
end

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
end;

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
end;

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
end;

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
end;

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
    if existfunctions('existendo',scarno(sloop));
      scar2use(sloop) = scar2use(sloop)+0.5;
    end
    if existfunctions('existepi',scarno(sloop));
      scar2use(sloop) = scar2use(sloop)+0.5;
    end
    if not(isempty(SET(scarno(sloop)).Scar));
      scar2use(sloop) = scar2use(sloop)+4;
    end
    if isequal(lower(SET(scarno(sloop)).ImageViewPlane),'short-axis');
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
    end;
  end;
  scarno = scarno(take);
  if isempty(scarno)
    scarno = oldscarno;
  end;
  
  if length(scarno)>1 
    mywarning('Detected multiple scar data. Taking data with maximal scar volume (arbitrary decision)');
    s = scarno(1);
    maxml = 0;
    for sloop = 1:length(scarno)
      try
        ml = SET(scarno(sloop)).Scar.Percentage*SET(scarno(sloop)).LVM/100; %scar volume in ml
        if ml>maxml
          maxml = ml;
          s = scarno(sloop);
        end;
      catch  %#ok<CTCH>
      end;
    end;
    scarno = s;
  elseif length(scarno)>1
    mywarning('Detected multiple scar data. Taking first image stack (arbitrary decision)');
    scarno = scarno(1);
  end
end;

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
  end;
  marno = marno(mar2use);
  
  if length(marno)>1
    mywarning('Detected multiple MaR data. Taking data with maximal MaR volume (arbitrary decision)');
    s = marno(1);
    maxml = 0;
    for sloop = 1:length(marno)
      try
        ml = SET(marno(sloop)).MaR.Percentage*SET(marno(sloop)).LVM/100; %scar volume in ml
        if ml>maxml
          maxml = ml;
          s = marno(sloop);
        end;
      catch  %#ok<CTCH>
      end;
    end;
    marno = s;
  end;
end;

%-----------------------------------------
function [flowno,flowroi] = findflowaxisno
%-----------------------------------------
%Find only one flow image stack
global SET

[~,~,flowno] = findno;
magno = [];
for noloop = 1:length(flowno)
  magno = [magno SET(flowno(noloop)).Flow.MagnitudeNo];
end
flowno = unique(magno);

%--- Check if multiple scardata.
if length(flowno)>1
  %Find best scar data to take. Take those with endo and scar data
  flow2use = false(size(flowno));
  for sloop=1:length(flowno)
    flow2use(sloop) = (SET(flowno(sloop)).RoiN>0);
  end;
  [maxval,maxind] = max(flow2use);
  allmax = find(flow2use == maxval);
  flowno = flowno(allmax);
  if length(flowno) > 1
    points = zeros(length(flowno),1);
    %find best image stack based on ROI names
    index = 1;
    for no = flowno
      for rloop = 1:SET(no).RoiN
        hasresult = (length(SET(no).Flow.Result)>=rloop && isempty(SET(no).Flow.Result(rloop)));
        if isequal(SET(no).Roi(rloop).Name,'Aortic ascending flow') && hasresult
          points(index) = 7;
        elseif isequal(SET(no).Roi(rloop).Name,'Pulmonary artery') && hasresult
          points(index) = 6;
        elseif isempty(strfind(SET(no).Roi(rloop).Name,'ROI-')) && hasresult
          points(index) = 5;
        elseif hasresult
          points(index) = 4;
        elseif isequal(SET(no).Roi(rloop).Name,'Aortic ascending flow')
          points(index) = 3;
        elseif isequal(SET(no).Roi(rloop).Name,'Pulmonary artery')
          points(index) = 2;
        elseif isempty(strfind(SET(no).Roi(rloop).Name,'ROI-'))
          points(index) = 1;
        end
      end
      index = index +1;
    end
    [maxval,maxind] = max(points);
    allmax = find(points == maxval);
    flowno = flowno(allmax);
  end  
  if length(flowno)>1
    disp('Detected multiple flow image stacks.Taking first stack (arbitrary decision)');
    flowno = flowno(1);
  end
end
  
%find flow ROI based on ROI name
if ~isempty(flowno)
  points = zeros(SET(flowno).RoiN,1);
  for rloop = 1:SET(flowno).RoiN
    hasresult = (length(SET(flowno).Flow.Result)>=rloop && isempty(SET(flowno).Flow.Result(rloop)));
    if isequal(SET(flowno).Roi(rloop).Name,'Aortic ascending flow') && hasresult
      points(rloop) = 7;
    elseif isequal(SET(flowno).Roi(rloop).Name,'Pulmonary artery') && hasresult
      points(rloop) = 6;
    elseif isempty(strfind(SET(flowno).Roi(rloop).Name,'ROI-')) && hasresult
      points(rloop) = 5;
    elseif hasresult
      points(rloop) = 4;
    elseif isequal(SET(flowno).Roi(rloop).Name,'Aortic ascending flow')
      points(rloop) = 3;
    elseif isequal(SET(flowno).Roi(rloop).Name,'Pulmonary artery')
      points(rloop) = 2;
    elseif isempty(strfind(SET(flowno).Roi(rloop).Name,'ROI-'))
      points(rloop) = 1;
    end
  end
  [maxval,maxind] = max(points);
  flowroi = find(points == maxval);
  if length(flowroi)>1
    disp('Detected multiple flow ROI. Taking first stack (arbitrary decision)');
    flowroi = flowroi(1);
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
end;
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
end;

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
end;

if isempty(SET(no).EndoX)
  ind = false(SET(no).ZSize,1);
  return;
end;

if SET(no).TSize>1
  temp = not(isnan(squeeze(SET(no).EndoX(1,:,:))));
  if SET(no).ZSize==1
    ind = all(temp);
  else
    ind = sum(temp,1)==SET(no).TSize;
  end;
else
  ind = squeeze(not(isnan(SET(no).EndoX(1,1,:))));
end;

%-------------------------------------
function ind = findslicewithepiall(no) %#ok<DEFNU>
%-------------------------------------
%Find slices with endocard in all timeframes
global SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EpiX)
  ind = false(SET(no).ZSize,1);
  return;
end;

if SET(no).TSize>1
  temp = not(isnan(squeeze(SET(no).EpiX(1,:,:))));
  if SET(no).ZSize==1
    ind = all(temp);
  else
    ind = sum(temp,1)==SET(no).TSize;
  end;
else
  ind = squeeze(not(isnan(SET(no).EpiX(1,1,:))));
end;

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
end;

temp = not(isnan(squeeze(SET(no).EndoX(1,tfs,:))));
if length(tfs)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end;
else
  ind = temp;
end;

%-------------------------------------
function ind = findslicewithrvendo(no,tfs) %#ok<DEFNU>
%-------------------------------------
%Find slices with RV endocard in any timeframe
global SET NO

if nargin<1
  no = NO;
end;
if nargin < 2
  tfs = 1:SET(no).TSize;
end

if isempty(SET(no).RVEndoX)
  ind = false(SET(no).ZSize,1);
  return;
end;

temp = not(isnan(squeeze(SET(no).RVEndoX(1,tfs,:))));
if length(tfs)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end;
else
  ind = temp;
end;

%----------------------------------
function ind = findslicewithepi(no,tfs) %#ok<DEFNU>
%----------------------------------
%Find slices with endocard in any timeframe
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
end;

temp = not(isnan(squeeze(SET(no).EpiX(1,tfs,:))));
if length(tfs)>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end;  
else
  ind = temp;
end;

%------------------------------------
function ind = findslicewithrvepi(no) %#ok<DEFNU>
%------------------------------------
%Find slices with RV epicard in any timeframe
global SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).RVEpiX)
  ind = false(SET(no).ZSize,1);
  return;
end;

temp = not(isnan(squeeze(SET(no).RVEpiX(1,:,:))));
if SET(no).TSize>1
  ind = (sum(temp,1)>0)';
  if SET(no).ZSize==1
    ind = max(ind(:));
  end;  
else
  ind = temp;
end;

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
end;

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
end;

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
		end;
	end;
end;

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
end;

%-------------------------------------------------
function tfs = findframeswithsegmentation(type,no) %#ok<DEFNU>
%-------------------------------------------------
%Find timeframes in no containing segmentation of type
global SET

switch type
  case 'endo'
    tfs = squeeze(~all(isnan(SET(no).EndoX(1,:,:)),3));
  case 'epi'
    tfs = squeeze(~all(isnan(SET(no).EpiX(1,:,:)),3));
end
 
%---------------------------------------
function [ind] = findoutflowtractslices(no,tfs) %#ok<DEFNU>
%---------------------------------------
%Find slices with outflow tract (a wall thickness that is less than 2mm)
%tfs are the timeframes to search in
global SET NO

if (SET(no).ResolutionX+SET(no).ResolutionY)/2<0.5
  warning('The threshold to find outflow tract (2mm) is species dependent, and might be incorrect for this images.');
end;

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
end;

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
  end;
  ind(zloop) = element;
end;


%-------------------------------
function setstack_Callback(type) %#ok<DEFNU>
%------------------------------
%open interface for user to manually define flow image stack
global DATA SET
switch type
  case 'flow'
    stri = 'Select image stack for flow report';
  case 'lv'
    stri = 'Salect image stack for LV report';
  case 'rv'
    stri = 'Select image stack for RV report';
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
s = mymenu(stri,menuitems);
if ~isempty(s) && s~=0
  no = stacks(s);
  switch type
    case 'flow'
      DATA.FlowNO = no;
      DATA.FlowROI = SET(no).RoiN;
      set(DATA.Handles.flowstackpushbutton,'String',sprintf('Stack #%d',no));
    case 'lv'
      DATA.LVNO = no;
      set(DATA.Handles.lvstackpushbutton,'String',sprintf('Stack #%d',no));
    case 'rv'
      DATA.RVNO = no;
      set(DATA.Handles.rvstackpushbutton,'String',sprintf('Stack #%d',no));
  end
end
switch type
  case 'flow'
    DATA.flowreportupdate;
    DATA.updateflowaxes;
  case {'lv','rv'}
    DATA.volumereportupdate;
    DATA.updatevolumeaxes;
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
