function [varargout] = ctautocontrast(varargin)
%----------------------------------------------------
%ctautocontrast: Automatic contrast adjustment for CT images

if nargin == 0
  varargin{1} = 'autocontrastcardiac';
end

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%-------------------------------
function autocontrastcardiac(no)
%------------------------------
%Automatic contrast adjustment for Cardiac CT images

%Map so the LV blood pool intensity is 0.9 and the LV myocardium intensity
%is 0.5. Use the HU information from the image if LV segmentation exist
%function based on segment_main.m/autocontrast_Callback

global SET NO DATA

set(DATA.Handles.autocontrasticon,'state','off');

if nargin == 0
  no = NO;
end

%use HU values from image if LV segmentation exist
if sum(~isnan(SET(no).EpiX(:)))
  t = find(any(~isnan(SET(no).EpiX(1,:,:)),3),1,'first'); %using first time frame with LV segmentation
  %find midslices with LV segmentation
  startslice = find(~isnan(SET(no).EpiX(1,t,:)),1,'first');
  endslice = find(~isnan(SET(no).EpiX(1,t,:)),1,'last');
  nbrofslices = endslice-startslice+1;
  midslices = round(startslice+0.25*nbrofslices):round(endslice-0.25*nbrofslices);
  %make myocardial mask
  myomask = zeros(size(SET(no).IM));
  bloodmask = myomask;
  for zloop = midslices
    endomask = segment('createmask',...
      [SET(no).XSize SET(no).YSize],...
      SET(no).EndoY(:,t,zloop),SET(no).EndoX(:,t,zloop));
    epimask = segment('createmask',...
      [SET(no).XSize SET(no).YSize],...
      SET(no).EpiY(:,t,zloop),SET(no).EpiX(:,t,zloop));
    myomask(:,:,t,zloop) = epimask-endomask;
    bloodmask(:,:,t,zloop) = endomask;
  end
  %estimate HU values for LV blood and LV myocardium
  originalintensity = calcfunctions('calctruedata',SET(no).IM,no);
  LVmyo = originalintensity(find(myomask));
  LVblood = originalintensity(find(bloodmask));
  medianLVmyo = median(LVmyo);
  %medianLVblood = median(LVblood);
  sortedLVmyo = sort(LVmyo);
  sortedLVblood = sort(LVblood);
  includedmyo = sortedLVmyo(abs(sortedLVmyo-medianLVmyo)<100);
  HULVmyo = mean(includedmyo);
  medianLVblood = median(sortedLVblood(round(length(sortedLVblood)/3):end));
  includedblood = sortedLVblood(abs(sortedLVblood-medianLVblood)<100);
  HULVblood = mean(includedblood);
else
  originalintensity = calcfunctions('calctruedata',SET(no).IM,no);
  HULVmyo = 150; %LV myocardium is about 150 HU
  HULVblood = 600;  %LV blood pool is about 600 HU
end

myworkon; %says that Segment is busy working

if isfield(SET(no),'CT') && isfield(SET(no).CT,'transverseNO')
  transverseintensity = calcfunctions('calctruedata',SET(SET(no).CT.transverseNO).IM,SET(no).CT.transverseNO);
  maxIntOriginal = max(transverseintensity(:)); %3071;
  minIntOriginal = min(transverseintensity(:)); %-3024;
else
  maxIntOriginal = max(originalintensity(:)); 
  minIntOriginal = min(originalintensity(:)); 
end
diffIntOriginal = maxIntOriginal-minIntOriginal;

x1 = (HULVmyo-minIntOriginal)/diffIntOriginal; 
x2 = (HULVblood-minIntOriginal)/diffIntOriginal;
y1 = 0.5;
y2 = 0.95;%0.9;
% if isa(SET(no).IM,'int16')
%   y1 = 0.5*(SET(no).maxValue-SET(no).minValue);
%   y2 = 0.9*(SET(no).maxValue-SET(no).minValue);
% end

k = (y2-y1)/(x2-x1);
m = y1-k*x1;

upperthresh = (1-m)/k;
lowerthresh = (0-m)/k;
if isa(SET(no).IM,'int16')
  upperthresh = (SET(no).maxValue-m)/k;
  lowerthresh = (SET(no).minValue-m)/k;
end

maxInt = single(max(SET(no).IM(:)));
minInt = single(min(SET(no).IM(:)));

if isempty(lowerthresh)
  lowerthresh=minInt;
end
if isempty(upperthresh)
  upperthresh=maxInt;
end

%calcualte contrast and brightness so that the intensities below
%lowerthresh and above upperthresh gets saturated
contrast=(maxInt-minInt)/(upperthresh-lowerthresh);
brightness=(maxInt+minInt)/2-(upperthresh+lowerthresh)/2*contrast+0.5;%(maxInt+minInt)/2-(upperthresh-lowerthresh)/2+0.5;%lowerthresh-minInt+0.5;

SET(no).IntensityMapping.Contrast=contrast;
SET(no).IntensityMapping.Brightness=brightness;

%update image
drawfunctions('drawcontrastimage',no);

myworkoff;


%------------------------------
function autocontrastcardiacall
%------------------------------
%Do auto contrast for all loaded image stacks

global SET DATA

set(DATA.Handles.autocontrastallicon,'state','off');

h = mywaitbarstart(length(SET)*2,'Auto contrast in process',0,DATA.GUI.Segment);
for no = 1:length(SET)
  h = mywaitbarupdate(h);
  autocontrastcardiac(no);
  h = mywaitbarupdate(h);
end
mywaitbarclose(h);

