function [varargout] = perfusion(varargin)
%Perfusion part of Segment.
%Einar Heiberg, Helen Soneson

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%-----------------------
function initdefault(no)
%-----------------------
global SET NO

if nargin==0
  no = NO;
end

SET(no).Perfusion.Auto = false(size(SET(no).IM));
SET(no).Perfusion.Manual = repmat(int8(0),size(SET(no).IM));
SET(no).Perfusion.Result = false(size(SET(no).IM));
SET(no).Perfusion.MyocardMask = false(size(SET(no).IM));
SET(no).Perfusion.Percentage = zeros(1,SET(no).TSize);
SET(no).Perfusion.Mode = 'manual';
SET(no).Perfusion.UpdateDirectly = 0;

SET(no).Perfusion.MR = [];

SET(no).Perfusion.MPS.ScoringIschManual = [];
SET(no).Perfusion.MPS.ScoringFixManual = [];
SET(no).Perfusion.MPS.TPD = [];
SET(no).Perfusion.MPS.TPDLAD = [];
SET(no).Perfusion.MPS.TPDLCx = [];
SET(no).Perfusion.MPS.TPDRCA = [];
SET(no).Perfusion.MPS.SeverityIndex = [];
SET(no).Perfusion.MPS.NadirIndex = [];
SET(no).Perfusion.MPS.TPDchange = [];
SET(no).Perfusion.MPS.TPDchangeLAD = [];
SET(no).Perfusion.MPS.TPDchangeLCx = [];
SET(no).Perfusion.MPS.TPDchangeRCA = [];


%-----------------------------
function createmyocardmask(no)
%-----------------------------
global SET NO

if nargin==0
  no=NO;
end

if isempty(SET(no).Perfusion)
  initdefault(no);
end

if isempty(SET(no).EndoX)
  myfailed('No LV endocardium available.');
  return;
end

%Create mask, allways update, segmentation could have changed.
SET(no).Perfusion.MyocardMask = false([SET(no).XSize SET(no).YSize SET(no).TSize SET(no).ZSize]);
slicetomask = find(findfunctions('findslicewithepi',no));%find(findfunctions('findslicewithendo',no)&findfunctions('findslicewithepi',no));
for zloop=1:length(slicetomask)
  for tloop=1:SET(no).TSize
    mask = false(SET(no).XSize,SET(no).YSize);
    if not(isnan(SET(no).EpiX(1,tloop,slicetomask(zloop))))
      mask = segment('createmask',...
        [SET(no).XSize SET(no).YSize],...
        SET(no).EpiY(:,tloop,slicetomask(zloop)),...
        SET(no).EpiX(:,tloop,slicetomask(zloop)));
    end
    if not(isnan(SET(no).EndoX(1,tloop,slicetomask(zloop))))
      endomask = segment('createmask',...
        [SET(no).XSize SET(no).YSize],...
        SET(no).EndoY(:,tloop,slicetomask(zloop)),...
        SET(no).EndoX(:,tloop,slicetomask(zloop)));
      mask = mask & not(endomask);
    end
    SET(no).Perfusion.MyocardMask(:,:,tloop,slicetomask(zloop)) = mask;

    %Update manual mask
    manualtemp=SET(no).Perfusion.Manual(:,:,tloop,slicetomask(zloop));
    manualtemp(~mask)=int8(0);
    SET(no).Perfusion.Manual(:,:,tloop,slicetomask(zloop)) = manualtemp;
  end
end

%-----------------------
function calcvolume(no) 
%-----------------------
global SET NO

if nargin==0
  no=NO;
end

%Calculate percentage
for tloop=1:SET(no).TSize
  SET(no).Perfusion.Percentage(tloop) = 100*sum(sum(sum(SET(no).Perfusion.Result(:,:,tloop,:))))/sum(sum(sum(SET(no).Perfusion.MyocardMask(:,:,tloop,:))));
end


%--------------------------
function showmanualinteraction_Callback(arg)
%--------------------------
%Toggle hide state 
global DATA

if nargin==0
  arg = 'toggle';
end

switch arg
  case 'toggle'
    c = get(DATA.Handles.Perfusionshowmanualinteractionmenu,'checked');
    if isequal(c,'on')
      set(DATA.Handles.Perfusionshowmanualinteractionmenu,'checked','off');
    else
      set(DATA.Handles.Perfusionshowmanualinteractionmenu,'checked','on');      
    end
  case 'on'
    set(DATA.Handles.Perfusionshowmanualinteractionmenu,'checked','off');
  case 'off'
    set(DATA.Handles.Perfusionshowmanualinteractionmenu,'checked','on');
  case 'update'
    %Do nothing
  otherwise
    error('Invalid option to showmanualinteraction_Callback, should not occur.');    
end

%Update
c = get(DATA.Handles.Perfusionshowmanualinteractionmenu,'checked');
if isequal(c,'on')
  %on, i.e view => no hide
  segment('unlighttool',DATA.Tools.Perfusionhide);  
else
  segment('highlighttool',DATA.Tools.Perfusionhide);
end

update;

%-------------------------
function clearall_Callback
%-------------------------
global DATA SET NO

tools('enableundo');
SET(NO).Perfusion = [];

update;
segment('updatevolume');
viewfunctions('setview');  %drawfunctions('drawimageno');

%--------------
function update(no,silent)
%--------------
%Update
global DATA NO SET

if nargin==0
  no=NO;
end

if nargin<2
  silent=false;
end

if not(isempty(SET(no).Perfusion))
  SET(no).Perfusion.Result = SET(no).Perfusion.Auto;
  SET(no).Perfusion.Result(SET(no).Perfusion.Manual==int8(1)) = true;
  SET(no).Perfusion.Result(SET(no).Perfusion.Manual==int8(-1)) = false;
  SET(no).Perfusion.Result = SET(no).Perfusion.Result&SET(no).Perfusion.MyocardMask;

  calcvolume(no);
  segment('updatevolume');

  if isequal(SET(no).Perfusion.Mode,'manual')&&not(silent)
    viewfunctions('setview');  %drawfunctions('drawimageno'); 
  end
end

%---------------------------------------------------
function mask = getPerfusionmask(no,timeframe,slice)
%---------------------------------------------------
%Returns the mask for Perfusion, if it does not exist, then it is created.
global DATA SET 

if isempty(SET(no).Perfusion)
  initdefault;
end

if nargin<1
  error('Expected one input argument to getPerfusionmask.');
end

if nargin<2
  timeframe = SET(no).CurrentTimeFrame;
end

if nargin<3
  slice = SET(no).CurrentSlice;
end

mask = SET(no).Perfusion.Result(:,:,timeframe,slice);


%----------------------------
function drawhelper(no,panel)
%----------------------------
%Function to draw Perfusion contour on screeen used from drawimageslice,
%drawimagemontage

global DATA SET

%Get mask
switch DATA.ViewPanelsType{panel}
  case {'montage','montagerow','montagefit','sax3'}
    result = segment('reshape2layout',squeeze(SET(no).Perfusion.Result(:,:,SET(no).CurrentTimeFrame,:)),no,panel);
  case 'one'
    result = SET(no).Perfusion.Result(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
  otherwise 
    myfailed('Unknown viewtype in Perfusion-drawhelper.');
end

if (sum(result(:))>0) 
  
  %Ensure that we delete stuff
  set(DATA.Handles.imageaxes(panel),'NextPlot','add');
    
  %Normal Perfusion outline
  [c,DATA.Handles.Perfusioncontour{panel}] = contour(DATA.Handles.imageaxes(panel),double(result),[0.5 0.5]);  
  set(DATA.Handles.Perfusioncontour{panel},'linecolor',[1 1 1]);
  if DATA.Pref.LineWidth>0
    set(DATA.Handles.Perfusioncontour{panel},'linewidth',DATA.Pref.LineWidth);
  else
    set(DATA.Handles.Perfusioncontour{panel},'visible','off');
  end
else
  DATA.Handles.Perfusioncontour{panel} = [];
end  

%---------------------
function showedits(no)
%---------------------
%Show viability edits on screen as a temporary overlay.
global DATA SET NO

if nargin==0
  no = NO;
end

if isempty(SET(no).Perfusion)
  return;
end

%Check if menu is enabled
if isequal(get(DATA.Handles.hideoverlayicon,'state'),'on')
  manualinteraction = false;
else
  manualinteraction = true;
end

if not(manualinteraction)
  return;
end

tempnos=no;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

if nargin==0
  panelstodo = find(DATA.ViewPanels==no);
else
  panelstodo = DATA.CurrentPanel;
end

for panel=panelstodo

  %Determine if it is a RGB image or not.
  temp = get(DATA.Handles.imagehandle(panel),'CData');
  if ndims(temp)>2
    isrgbimage = true;
  else
    isrgbimage = false;
  end
  
  clear temp;
  
  %Ok lets draw it
  switch DATA.ViewPanelsType{panel}
    case {'one','mmodespatial'}
      temp = SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      
      if isrgbimage  %RGB iamge
        tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp,no,SET(no).Colormap));
      else  %Default grayscale image
        tempuint8 = min(uint8(255),uint8(1)+calcfunctions('remapuint8',temp));
        tempuint8 = repmat(tempuint8,[1 1 3]); %EH:
      end
      clear temp;
      sz=size(tempuint8);
      Perfusionimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);

      tmp=SET(no).Perfusion.Manual(:,:,SET(no).CurrentTimeFrame,SET(no).CurrentSlice);
      tmp=tmp(:);
      Perfusionimrgb((tmp==int8( 1)),2) = uint8(255);
      Perfusionimrgb((tmp==int8(-1)),3) = uint8(255);
      Perfusionimrgb=reshape(Perfusionimrgb,[sz(1:2) 3]); %EH: added 3

    case {'montage','montagerow','montagefit','sax3'}
      % Convert to 2D and layout

      if isrgbimage  %RGB iamge
        tempuint8 = calcfunctions('remapuint8',...
          segment('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:))),...
          no,panel,SET(no).Colormap);
      else  %Default grayscale image
        tempuint8 = calcfunctions('remapuint8',...
          segment('reshape2layout',squeeze(SET(no).IM(:,:,SET(no).CurrentTimeFrame,:))),...
          no,panel,calcfunctions('returnmapping',no,true));
      end
      tempuint8 = min(uint8(255),uint8(1)+tempuint8);
      sz = size(tempuint8);
      Perfusionimrgb = reshape(tempuint8,[prod(sz(1:2)) 3]);
      
      tmp=segment('reshape2layout',SET(no).Perfusion.Manual(:,:,SET(no).CurrentTimeFrame,:),no,panel);
      tmp=tmp(:);
      Perfusionimrgb((tmp==int8( 1)),2) = uint8(255);
      Perfusionimrgb((tmp==int8(-1)),3) = uint8(255);
      Perfusionimrgb=reshape(Perfusionimrgb,sz);
  end

  set(DATA.Handles.imagehandle(panel),'CData',Perfusionimrgb);
end

%---------------------------------
function importfromscar_Callback
%---------------------------------
global SET NO DATA

no = NO;

if isempty(SET(no).Scar)
  myfailed('No Scar to import from.');
end

if isempty(SET(no).Perfusion)
  initdefault;
end

if SET(no).TSize>1
  myfailed('Can only import from Scar when Perfusion is non-timeresolved.');
end

createmyocardmask;

SET(no).Perfusion.Auto(:,:,1,:)=SET(no).Scar.Auto;
%Use NoReflow to add in manual
manual = SET(no).Scar.Manual;
manual(SET(no).Scar.NoReflow) = int8(1);
SET(no).Perfusion.Manual(:,:,1,:)=manual;

update;

%--------------------------------------------------------------
function importfromimagestackwithdifferentcrop_Callback(fromno)
%--------------------------------------------------------------
global SET NO DATA

no = NO;

if nargin==0
  %Find what imagestack
  s = myinputdlg({'Enter what image stack to take from '},'ImageStack',1,{sprintf('%d',1)});

  if isempty(s)
    myfailed('Invalid image stack.',DATA.GUI.Segment);
    return;
  else
    [fromno,ok] = str2num(s{1}); %#ok<ST2NM>
    if not(ok)
      myfailed('Invalid image stack.',DATA.GUI.Segment);
      return;
    end
  end
end

if fromno==no
  myfailed('Cannot import from same image stack.',DATA.GUI.Segment);
  return;
end
if (fromno>length(SET))||(fromno<1)
  myfailed('Invalid image stack selected.',DATA.GUI.Segment);
  return;
end

if isempty(SET(fromno).Perfusion)
  myfailed('Perfusion must be present in image stack to import form');
  return;
end

if isempty(SET(no).Perfusion)
  initdefault;
end

createmyocardmask;

maskindexinno = find(SET(no).Perfusion.MyocardMask(:));
maskindexinfromno = find(SET(fromno).Perfusion.MyocardMask(:));
SET(no).Perfusion.Auto(maskindexinno) = SET(fromno).Perfusion.Auto(maskindexinfromno);
SET(no).Perfusion.Manual(maskindexinno )= SET(fromno).Perfusion.Manual(maskindexinfromno);

update;

%---------------------------------
function exporttoscar_Callback
%---------------------------------
global SET NO DATA

no = NO;

if isempty(SET(no).Perfusion)
  myfailed('No Perfusion to export from.');
end

if SET(NO).TSize>1
  myfailed('Can only exprot to Scar when Perfusion is non-timeresolved.');
end

if isempty(SET(no).Scar)
  viability('viabilityreset_Callback','weighted');
end

SET(no).Scar.Mode='manual';
SET(no).Scar.Auto=squeeze(SET(no).Perfusion.Auto(:,:,1,:));
%Use NoReflow to add in manual
SET(no).Scar.Manual=squeeze(SET(no).Perfusion.Manual(:,:,1,:));

viability('viabilitycalc');

%----------------------------------------
function clearmanualinteraction_Callback
%----------------------------------------
global SET NO DATA

if isempty(SET(NO).Perfusion)
  initdefault;
end

tools('enableundo');
%Clear manual interaction
SET(NO).Perfusion.Manual = repmat(int8(0),size(SET(NO).Perfusion.Manual));
update;


%---------------------------------------------------------------
function [res,varargout] = calctransmuralityline(numsectors,no)
%---------------------------------------------------------------
%Calculates transmurality.
%
%Output in order
%- Mean transmurality in the number of sectors
%- Max transmurality in each sector
%- Mean transmurality calculated only over infarcted areas.
%- Total Extent 
%- "Start" transmurality
%- "End" transmurality

global DATA SET NO

if nargin<2
  no = NO;
end

pos = SET(no).StartSlice:SET(no).EndSlice;
numslices = length(pos);
res = zeros(numsectors,numslices);
startind = res;
endind = res;
if isempty(SET(no).Perfusion)
  myfailed('No Perfusion data available.',DATA.GUI.Segment);  
  return;
end

if isempty(SET(no).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end
if isempty(SET(no).EpiX)
  myfailed('No LV epicardium available.',DATA.GUI.Segment);
  return;
end
  
%Upsample model
if not(DATA.Pref.RadialProfiles==DATA.NumPoints)
  [endox,endoy] = calcfunctions('resamplemodel',SET(no).EndoX,SET(no).EndoY,DATA.Pref.RadialProfiles);
  [epix,epiy] = calcfunctions('resamplemodel',SET(no).EpiX,SET(no).EpiY,DATA.Pref.RadialProfiles);
else
  endox = SET(no).EndoX;
  endoy = SET(no).EndoY;
  epix = SET(no).EpiX;
  epiy = SET(no).EpiY;
end

maxtrans = zeros(numsectors,numslices);
transmurality = zeros(DATA.Pref.RadialProfiles,1,SET(no).ZSize);
sectorstartind = nan(size(transmurality));
sectorendind = sectorstartind;
tf = SET(no).CurrentTimeFrame;

for sloop=1:length(pos) %Loop over slices
  %Check if there is Perfusion else take next slide.
  
  %Find center
  mx = mean(endox(:,tf,pos(sloop)));
  my = mean(endoy(:,tf,pos(sloop)));
   
  %Follow along the epicardium
  angleendo = angle(complex(endox(:,tf,pos(sloop))-mx,endoy(:,tf,pos(sloop))-my))*180/pi;
  angles2find = angle(complex(epix(:,tf,pos(sloop))-mx,epiy(:,tf,pos(sloop))-my))*180/pi;  

  if max(max(SET(no).Perfusion.Result(:,:,tf,pos(sloop))))>0
    for loop=1:DATA.Pref.RadialProfiles %Loop over radial profiles

      %Find correspondance for endo      
      [trash,indendo] = min(abs(angleendo(:)-angles2find(loop)));

      %Find correspondance for epi
      indepi = loop; %Follow along the epicardium
      
      x = linspace(...
        endox(indendo,tf,pos(sloop)),...
        epix(indepi,tf,pos(sloop)),50);
      y = linspace(...
        endoy(indendo,tf,pos(sloop)),...
        epiy(indepi,tf,pos(sloop)),50);

      idx = sub2ind([SET(no).XSize SET(no).YSize],round(x),round(y));
      slicetempim = SET(no).Perfusion.Result(:,:,tf,pos(sloop));

      transmurality(indepi,1,pos(sloop)) = 100*sum(slicetempim(idx))/50;
      sectorind = find(slicetempim(idx));

      if length(sectorind)>2
        sectorendind(indepi,1,pos(sloop)) = (sectorind(end)-1)/49;
        sectorstartind(indepi,1,pos(sloop)) = (sectorind(1)-1)/49;
      end
    end %loop over datapoints
    tempres = calcfunctions('findmeaninsector','epi',transmurality,pos(sloop),numsectors,no);  % [res(:,sloop),tempa,tempa,tempa,maxtrans(:,sloop)] = 
    res(:,sloop) = tempres(:,tf);
    tempres = calcfunctions('findmeaninsector','epi',sectorendind,pos(sloop),numsectors,no);
    startind(:,sloop) = tempres(:,tf);
    tempres = calcfunctions('findmeaninsector','epi',sectorstartind,pos(sloop),numsectors,no);  
    endind(:,sloop) = tempres(:,tf);
  end %Perfusion exist
end %loop over slices

meantrans4infarct = sum(transmurality(:))./sum(transmurality(:)>0);
tempendorad = calcfunctions('calcendoradius',no);
endorad = tempendorad(:,tf,:);
totextent = 100*nansum(endorad(:).*(transmurality(:)>0))/nansum(endorad(:));

%Assign optional output
varargout =   cell(nargout-1,1);
if nargout>1
  varargout{1} = maxtrans;
end
if nargout>2
  varargout{2} = meantrans4infarct;
end
if nargout>3
  varargout{3} = totextent;
end
if nargout>4
  varargout{4} = startind;
end
if nargout>5
  varargout{5} = endind;
end



