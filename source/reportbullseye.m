function [varargout] = reportbullseye(varargin)
%GUI for reporting bullseye

%Adapted for use with mygui class by Nils Lundahl

%Commented for technical manual by Einar Heiberg
%#ok<*GVMIS>
if nargin == 0 || isempty(varargin{1})
  varargin{1} = 'init';
elseif strncmp(varargin{1},'get',3)
  [varargout{1:nargout}] = feval('getdata',varargin{:}); %#ok<FVAL>
  return
elseif strncmp(varargin{1},'set',3)
  [varargout{1:nargout}] = feval('setdata',varargin{:}); %#ok<FVAL>
  return
end
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init
%------------
%Initialize the GUI

global DATA SET NO

if isempty(SET(NO).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
elseif isempty(SET(NO).EpiX)
  myfailed('No LV epicardium available.',DATA.GUI.Segment);
  return;
end

tempnos = NO;
imissingle = classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end
startbullseye;

%---------------------
function startbullseye
%---------------------
%start up the bullseye analysis

global DATA NO SET
% 2025-10-01 check data and angle
no = NO;
segment('checkconsistency',1:SET(no).TSize,1:SET(no).ZSize);
% Update SectorRotation if RV insertion points exist
reportbullseye('sectorrotationhelper',NO);

% Ask for LV rotation if that is not defined
if SET(no).SectorRotation == 0
  calcfunctions('autosectorrotation_Callback');
  LVrotation;
end


%Find longaxis cine image. Order of preference is 4CH > 3CH > 2CH
cineno = findfunctions('findno');
laxno = [];
for cham = 2:4
  chno = find(strcmp({SET(cineno).ImageViewPlane},sprintf('%dCH',cham)),1);
  if ~isempty(chno) && chno ~= NO
    laxno = cineno(chno);
  end
end

%Initalize
DATA.initbullseyeslices; %Set start and end slice in GUI's with such fcnality
if isempty(SET(NO).StartSlice)
  myfailed('No slices selected.',DATA.GUI.Segment);
  %   close_Callback;
  return;
end

%ensure current time frame have LV seg
tf = SET(NO).CurrentTimeFrame;
if isnan(mynanmean(mynanmean(SET(NO).EndoX(:,tf,:))))
  %change to first time frame with LV seg
  tffirst = find(~isnan(mynanmean(SET(NO).EndoX(1,:,:),3)),1,'first');
  SET(NO).CurrentTimeFrame = tffirst;
  DATA.ViewIM{DATA.CurrentPanel} = [];
  drawfunctions('drawno',NO);
end

%KG, 2019-06-19
%Took away resampling below because function for removing duplicate points already exists in
%calcfunctions('resamplecurve')
%
%Resample model to ensure distributed points.
% for tloop=1:SET(NO).TSize
%   for zloop=1:SET(NO).ZSize
%
%     %Check endo
%     if not(isnan(SET(NO).EndoX(1,tloop,zloop)))
%       %Calculate distances between points
%       len = sqrt(...
%         conv2(SET(NO).EndoX(:,tloop,zloop),[1;-1],'valid').^2+...
%         conv2(SET(NO).EndoY(:,tloop,zloop),[1;-1],'valid').^2);
%       %Remove indices of duplicate points (len(i) == 0)
%       ind = [true; len > 0];
%       len = len(ind(2:end));
%       len = [0;len(:)]; %Add zero first
%       len = cumsum(len);
%       totallength = len(end);
%
%       if totallength>0
%         %Resample
%         SET(NO).EndoX(:,tloop,zloop) = interp1(len,SET(NO).EndoX(ind,tloop,zloop),linspace(0,totallength,DATA.NumPoints));
%         SET(NO).EndoY(:,tloop,zloop) = interp1(len,SET(NO).EndoY(ind,tloop,zloop),linspace(0,totallength,DATA.NumPoints));
%       end
%     end
%
%     %Check epi
%     if ~isempty(SET(NO).EpiX)
%       if not(isnan(SET(NO).EpiX(1,tloop,zloop)))
%         len = sqrt(...
%           conv2(SET(NO).EpiX(:,tloop,zloop),[1;-1],'valid').^2+...
%           conv2(SET(NO).EpiY(:,tloop,zloop),[1;-1],'valid').^2);
%          %Remove indices of duplicate points (len(i) == 0)
%          ind = [true; len > 0];
%          len = len(ind(2:end));
%          len = [0;len(:)]; %Add zero first
%          len = cumsum(len);
%          totallength = len(end);
%
%         if totallength>0
%           %Resample
%           SET(NO).EpiX(:,tloop,zloop) = interp1(len,SET(NO).EpiX(ind,tloop,zloop),linspace(0,totallength,DATA.NumPoints));
%           SET(NO).EpiY(:,tloop,zloop) = interp1(len,SET(NO).EpiY(ind,tloop,zloop),linspace(0,totallength,DATA.NumPoints));
%         end
%       end
%     end
%
%   end
% end

gui = mygui('bullseye.fig');
DATA.GUI.Bullseye = gui;
gui.no = NO;
set(gui.fig,'renderer','opengl');
set(gui.handles.textimagestack,'String',sprintf('%s %d %s',dprintf('Image stack'),gui.no,SET(gui.no).ImageType));
DATA.initbullseye;  %fix listbox for each software package
gui.laxno = laxno;
gui.Endocardium_percent = [];
gui.Epicardium_percent = [];

%Update
rotationangle = LVrotation('setrotationangle_helper');
set(gui.handles.rotationslider,'value',rotationangle);
addlistener(gui.handles.rotationslider, 'Value', 'PreSet', @slidermotion_Callback);
set(gui.handles.endocentercheckbox,'value',SET(NO).EndoCenter);
set(gui.handles.volumeconsistentcheckbox,'value',1);
gui.volumeconsistent = true;
temp = get(gui.handles.sectorslistbox,'String');
gui.slice = round(0.5*(SET(NO).StartSlice+SET(NO).EndSlice)); %for which image with spoke wheel is plotted.

numslices = SET(gui.no).StartSlice:SET(gui.no).EndSlice;
if length(numslices) >= 3
  set(gui.handles.aharadiobutton,'Value',0);
  set(gui.handles.aha16radiobutton,'Value',1);
  set(gui.handles.sectorradiobutton,'Value',0);
  set(gui.handles.smoothradiobutton,'Value',0);
  gui.numsectors = 24;
  set(gui.handles.sectorslistbox,'enable','off');
  set(gui.handles.volumeconsistentcheckbox,'enable','off');
else
  set(gui.handles.aha16radiobutton,'Value',0);
  set(gui.handles.aharadiobutton,'Value',0);
  set(gui.handles.sectorradiobutton,'Value',1);
  set(gui.handles.smoothradiobutton,'Value',0);
  gui.numsectors = str2num(temp{mygetlistbox(gui.handles.sectorslistbox)}); %#ok<ST2NM>
end

set(gui.handles.slicetext,'String',sprintf('%s %d',dprintf('Slice'),gui.slice));
if SET(NO).ZSize > 1
  set(gui.handles.sliceslider,'Min',1,'Max',...
    SET(NO).ZSize,'Value',SET(NO).ZSize-gui.slice+1,'SliderStep',...
    [1 3]/(SET(NO).ZSize));
else
  set(gui.handles.sliceslider,'Visible','off');
  set(gui.handles.slicetext,'Visible','off');
end
gui.outdata = []; %extracted outdata
gui.ahaoutdata = []; %in 17 segment format.

if isempty(gui.slice)
  gui.slice = 1;
end

if ~DATA.GUISettings.PointsEnabled
  set(gui.handles.rotationfromannotationcheckbox,'visible','off','value',0);
  set(gui.handles.endocentercheckbox,'visible','off','value',0);
end

%init colormap listbox
gui.colormaplist = {...
  'Jet',...
  'Hot',...
  'HSV',...
  'SPECT',...
  dprintf('Gray'),...
  'Gadgetron',...
  };
set(gui.handles.colormaplistbox,'string',gui.colormaplist);

%initialize lonaxis image, true long-axis or reconstructed from short-axis
initlaximage;

% %Initiate longaxis image, if available. Otherwise hide axis
% if isfield(gui,'laxno') && ~isempty(gui.laxno)
%   initlaximage;
% else
%   set(gui.handles.laximageaxes,'Visible','off');
% end

%Set selection change function for plot method panel
set(gui.handles.plotmethodpanel,'SelectionChangeFcn', ...
  'reportbullseye(''plotmethodpanel_SelectionChange'')');

updateall;

%--------------------
function slidermotion_Callback(~,~)
%--------------------
% update image when slider was dragged
updatesliceimage

%--------------------
function initlaximage
%--------------------
%Initiate long axis image in bullseye gui
global DATA SET
gui = DATA.GUI.Bullseye;

tf = SET(gui.no).CurrentTimeFrame;
if ~isempty(gui.laxno)
  %find closest tf in normalized heart beat this fixes if LAX image has
  %different number of timeframes.
  [~,tf] = min(abs(SET(gui.laxno).TimeVector/SET(gui.laxno).TimeVector(end) - SET(gui.no).TimeVector(tf)/SET(gui.no).TimeVector(end)));
  temp = SET(gui.laxno).IM(:,:,tf);
else
  if SET(gui.no).EndoCenter
    hslice = round(mynanmean(mynanmean(SET(gui.no).EndoX(:,tf,:))));
    if isnan(hslice)
      hslice = round(mynanmean(mynanmean(mynanmean(SET(gui.no).EndoX(:,:,:)))));
    end
  else
    hslice = round(mynanmean(mynanmean(SET(gui.no).EpiX(:,tf,:))));
    if isnan(hslice)
      hslice = round(mynanmean(mynanmean(mynanmean(SET(gui.no).EpiX(:,:,:)))));
    end
  end
  %hslice = round(SET(gui.no).XSize/2);
  temp = permute(SET(gui.no).IM(hslice,:,tf,:),[4 2 3 1]);
end
temp = min(max(temp,0),1);

%force true color
image(calcfunctions('remapuint8',temp,gui.no,calcfunctions('returnmapping',gui.no,true)),...
  'parent',gui.handles.laximageaxes);
if ~isempty(gui.laxno)
  set(gui.handles.laximageaxes,'PlotBoxAspectRatio', ...
    [SET(gui.laxno).YSize*SET(gui.laxno).ResolutionY ...
    SET(gui.laxno).XSize*SET(gui.laxno).ResolutionX 1]);
else
  set(gui.handles.laximageaxes,'PlotBoxAspectRatio', ...
    [SET(gui.no).YSize*SET(gui.no).ResolutionY ...
    SET(gui.no).ZSize*(SET(gui.no).SliceGap+SET(gui.no).SliceThickness) 1]);
end
axis(gui.handles.laximageaxes,'off');
hold(gui.handles.laximageaxes,'on');

%init vector of handles for intersection lines
zsz = SET(gui.no).ZSize;
gui.handles.intersectionlines = nan(zsz,1);

for loop = 1:zsz
  if ~isempty(gui.laxno)
    [x,y] = calcfunctions('calcplaneintersections',gui.laxno,gui.no,'one','one',loop);
  else
    [x,y] = calcfunctions('calcplaneintersections',gui.no,gui.no,'hla','one',loop,hslice);
  end
  if isempty(x) || isempty(y)
    x = nan;
    y = nan;
  end
  if loop >= SET(gui.no).StartSlice && loop <= SET(gui.no).EndSlice
    gui.handles.intersectionlines(loop) = plot(gui.handles.laximageaxes,y,x,'y');
  else
    gui.handles.intersectionlines(loop) = plot(gui.handles.laximageaxes,y,x,'w');
  end
end
% myset(gui.handles.intersectionlines,'Visible','off');

%-------------------------------
function defineslices(sign,part)
%-------------------------------
%manual define slices to include in the bullsye plot

global DATA SET
gui = DATA.GUI.Bullseye;

indendo = findfunctions('findslicewithendo',gui.no);
indepi = findfunctions('findslicewithepi',gui.no);
indlv = min(1,indendo+indepi);
switch part
  case 'basal'
    switch sign
      case 'plus'
        mostbasallv = find(indlv,1,'first');
        if ~isempty(mostbasallv)
          SET(gui.no).StartSlice = max(mostbasallv,SET(gui.no).StartSlice-1);
        end
      case 'minus'
        SET(gui.no).StartSlice = min(SET(gui.no).StartSlice+1,SET(gui.no).EndSlice);
    end
  case 'apical'
    switch sign
      case 'plus'
        mostapicallv = find(indlv,1,'last');
        if ~isempty(mostapicallv)
          SET(gui.no).EndSlice = min(mostapicallv,SET(gui.no).EndSlice+1);
        end
      case 'minus'
        SET(gui.no).EndSlice = max(SET(gui.no).StartSlice,SET(gui.no).EndSlice-1);
    end
end
%check so gui.slice is inside selected slices
gui.slice = min(max(gui.slice,SET(gui.no).StartSlice),SET(gui.no).EndSlice);
set(gui.handles.sliceslider,'Value',SET(gui.no).ZSize-gui.slice+1);
set(gui.handles.slicetext,'String',sprintf('%s %d',dprintf('Slice'),gui.slice));
updatesliceimage;

%update all views
updatelongaxisimage;
updateall;

%-----------------
function updateall
%-----------------
%Update entire GUI, calls rotationslide and listbox. Thats all.

updatesliceimage;
updateplot;

%---------------------------------------
function [stri,pos] = aha17nameandpos(i)
%---------------------------------------
%Get name of 17 segment

switch i
  case 1
    stri = 'Basal anteroseptal';
    pos = 2;
  case 2
    stri = 'Basal anterior';
    pos = 1;
  case 3
    stri = 'Basal anterolateral';
    pos = 6;
  case 4
    stri = 'Basal inferolateral';
    pos = 5;
  case 5
    stri = 'Basal inferior';
    pos = 4;
  case 6
    stri = 'Basal inferoseptal';
    pos = 3;
  case 7
    stri = 'Mid anteroseptal';
    pos = 8;
  case 8
    stri = 'Mid anterior';
    pos = 7;
  case 9
    stri = 'Mid anterolateral';
    pos = 12;
  case 10
    stri = 'Mid inferolateral';
    pos = 11;
  case 11
    stri = 'Mid inferior';
    pos = 10;
  case 12
    stri = 'Mid inferoseptal';
    pos = 9;
  case 13
    stri = 'Apical septal';
    pos = 14;
  case 14
    stri = 'Apical anterior';
    pos = 13;
  case 15
    stri = 'Apical lateral';
    pos = 16;
  case 16
    stri = 'Apical inferior';
    pos = 15;
  case 17
    stri = 'Apex';
    pos = 17;
end

%--------------------------
function varargout = export
%--------------------------
%Export data to clipboard.

global DATA SET
gui = DATA.GUI.Bullseye;

%Do graphical update
updateplot;

%--- Normal output
if mygetvalue(gui.handles.sectorradiobutton)
  %init cell
  nsectors = size(gui.outdata, 1);
  nslices = size(gui.outdata, 2);
  outcell = cell(nslices+1, nsectors+4);
  %print general header
  outcell{1, 1} = 'Parameter';
  outcell{1, 2} = 'StartSlice';
  outcell{1, 3} = 'EndSlice';
  outcell{1, 4} = 'Slice';
  %print general info
  outcell{2, 1} = gui.parameter;
  outcell{2, 2} = sprintf('%d',SET(gui.no).StartSlice);
  outcell{2, 3} = sprintf('%d',SET(gui.no).EndSlice);
  %print sector header
  for sectorloop = 1:nsectors
    outcell{1, nsectors+4-sectorloop+1} = sprintf('Sector%d',sectorloop);
  end
  %print sector data
  for sliceloop = 1:nslices
    if sliceloop == 1
      outcell{sliceloop+1, 4} = sprintf('%d (basal)',sliceloop);
    else
      if sliceloop == nslices
        outcell{sliceloop+1, 4} = sprintf('%d (apical)',sliceloop);
      else
        outcell{sliceloop+1, 4} = sprintf('%d',sliceloop);
      end
    end
    for sectorloop = 1:nsectors
      outcell{sliceloop+1, sectorloop+4} = ...
        sprintf('%0.5g',gui.outdata(sectorloop,sliceloop));
    end
  end
  %send to clipboard
  if nargout == 0
    segment('cell2clipboard',outcell);
  else
    varargout{1} = outcell;
  end

  %--- Smooth output to .bin file
elseif mygetvalue(gui.handles.smoothradiobutton)
  tempim = get(get(gui.handles.bullseyeaxes,'children'),'cdata');

  if nargout==0
    filename = myinputdlg({'Enter file name'},'File name',1,{sprintf('%s_%s',SET(gui.no).PatientInfo.ID,gui.parameter)});
    fid = fopen([filename{1},'.bin'],'w');
    fwrite(fid,tempim,'float32');
    status = fclose(fid);
    if status == 0
      mymsgbox(['Data saved in file ',filename{1},'.bin'],'Done!',DATA.GUI.Segment);
    else
      myfailed('Failed to save data to file');
    end
  else
    stri = sprintf('%s\t%s\t%d\t%d\n',gui.parameter,gui.outunit,SET(gui.no).StartSlice,SET(gui.no).EndSlice);
    stri = [stri sprintf('\t')];
    stri = [stri newline];
    for rowloop = 1:size(tempim,1)
      for columnloop = 1:size(tempim,2)
        stri = [stri sprintf('%0.5g\t',tempim(rowloop,columnloop))]; %#ok<AGROW>
      end
      stri = [stri newline]; %#ok<AGROW>
    end
    varargout{1} = stri;
  end

  %--- AHA output
elseif mygetvalue(gui.handles.aharadiobutton) || mygetvalue(gui.handles.aha16radiobutton)
  nsegments = 17;
  outcell = cell(2,nsegments+3);
  outcell{1, 1} = 'Parameter';
  outcell{1, 2} = 'StartSlice';
  outcell{1, 3} = 'EndSlice';
  outcell{2, 1} = gui.parameter;
  outcell{2, 2} = sprintf('%d',SET(gui.no).StartSlice);
  outcell{2, 3} = sprintf('%d',SET(gui.no).EndSlice);
  for loop = 1:nsegments
    [stri,pos] = aha17nameandpos(loop); %Get name and position of export
    outcell{1,pos+3} = sprintf('%s [%s]',stri,gui.outunit);
    outcell{2,pos+3} = gui.ahaoutdata(loop);
  end
  if nargout == 0
    segment('cell2clipboard',outcell);
  else
    varargout{1} = outcell;
  end
end

%---------------------------------------
function plotmethodpanel_SelectionChange
%---------------------------------------
%Callback for changing type of plot (sector, smooth or AHA model)
global DATA SET
gui = DATA.GUI.Bullseye;

%If AHA model then use 24 sectors 4*6
if mygetvalue(gui.handles.aharadiobutton) || mygetvalue(gui.handles.aha16radiobutton)
  tempslices = SET(gui.no).StartSlice:SET(gui.no).EndSlice;
  if length(tempslices)<3 %|| mygetvalue(gui.handles.thissliceonlycheckbox)
    myfailed('Expected at least three slices for AHA segment model.');
    set(gui.handles.sectorradiobutton,'Value',1);
  else
    gui.numsectors = 24;
    set(gui.handles.sectorslistbox,'enable','off');
    set(gui.handles.volumeconsistentcheckbox,'enable','off');
  end

else
  set(gui.handles.sectorslistbox,'enable','on');
  temp = get(gui.handles.sectorslistbox,'String');
  gui.numsectors = str2double(temp{mygetlistbox(gui.handles.sectorslistbox)});
  set(gui.handles.volumeconsistentcheckbox,'enable','on');
  gui.volumeconsistent = get(gui.handles.volumeconsistentcheckbox,'value');
end
updateall;

%---------------
function sectors
%---------------
%called when number of sectors is updated.
global DATA
gui = DATA.GUI.Bullseye;

temp = get(gui.handles.sectorslistbox,'String');
gui.numsectors = str2num(temp{mygetlistbox(gui.handles.sectorslistbox)}); %#ok<ST2NM>
updateall;

%------------------
function endocenter
%------------------
%Called when the endocenter checkbox is clicked. Updates SET.EndoCenter.
global DATA SET
gui = DATA.GUI.Bullseye;

SET(gui.no).EndoCenter = get(gui.handles.endocentercheckbox,'value');
updateall;

%------------------------
function volumeconsistent
%------------------------
%Called when volumeconsistent checkbox is called.
global DATA
gui = DATA.GUI.Bullseye;

gui.volumeconsistent = get(gui.handles.volumeconsistentcheckbox,'value');
updateall;

%--------------------------------------
function pos = sectorrotationhelper(no)
%--------------------------------------
%Find suitable sector rotation based on RV insertion points

global SET

%Find slices with RV insertion points
slices = false(1,SET(no).ZSize);
for loop = 1:length(SET(no).Point.Z)
  if isequal(SET(no).Point.Label{loop},'RV insertion') || isequal(SET(no).Point.Label{loop},'P1')|| isequal(SET(no).Point.Label{loop},'P2')
    slices(SET(no).Point.Z(loop)) = true;
  end
end

%Find slices
pos = find(slices);
if isempty(pos)
  logdisp(sprintf('No RV points found. Current sector rotation is %0.5g',SET(no).SectorRotation));
  return
end

%Check if endo or epi
doepi = true;
if isempty(SET(no).EpiX)
  doepi = false;
end

%Find sector rotation
sectorrot = 0; %mean sector rotation
for loop = 1:length(pos)

  %Find RV insertion points for this contour
  ind = [];
  for rloop = 1:length(SET(no).Point.Z)
    if (isequal(SET(no).Point.Label{rloop},'RV insertion') || isequal(SET(no).Point.Label{rloop},'P2') || isequal(SET(no).Point.Label{rloop},'P1')) && isequal(pos(loop),SET(no).Point.Z(rloop))
      ind = [ind rloop]; %#ok<AGROW>
    end
  end
  % check time frames where both points and contours exist
  if doepi
    lvtf = find(~isnan(mynanmean(SET(no).EpiX(1,:,pos(loop)),3)));
  else
    lvtf = find(~isnan(mynanmean(SET(no).EndoX(1,:,pos(loop)),3)));
  end

  pnttf = SET(no).Point.T(ind);
  if any(isnan(pnttf))
    pnttf = 1:SET(no).TSize;
    istimeresolved = true;
  else
    istimeresolved = false;
  end
  tf = intersect(lvtf,pnttf);

  numtf = length(tf);
  if numtf == 0
    fprintf('No RV insertion points and segmentation in one and the same time frame in slice %d\n',pos(loop));
    continue
  elseif numtf > 1
    edt = SET(no).EDT;
    if ismember(edt,tf)
      % take edt
      tf = edt;
    else
      [~,indmin] = min(abs(tf-edt));
      tf = tf(indmin(1));
    end
    fprintf('RV insertion points and segmentation are taken from time frame %d\n',tf);
    if ~istimeresolved
      ind = ind(pnttf == tf);
    end
  end

  if length(ind) < 2
    message = dprintf('Expected two insertion points, found %d.',length(ind));
    disp(message)
    myfailed(message);
    return;
  end

  if ~isequal(length(ind),2)
    fprintf('Expected two insertion points, found %d. Taking the two points with largest distance between.\n',length(ind));
    dist = nan(length(ind),length(ind));
    for indloop = 1:length(ind)
      dist(:,indloop) = sqrt((SET(no).Point.X(ind(indloop))-SET(no).Point.X(ind)).^2+ ...
        (SET(no).Point.Y(ind(indloop))-SET(no).Point.Y(ind)).^2);
    end
    [~,largestdist] = max(dist(:));
    [ind1,ind2] = ind2sub([length(ind) length(ind)],largestdist);
    ind = sort([ind(ind1) ind(ind2)]);
  end

  p1x = SET(no).Point.X(ind(1));
  p1y = SET(no).Point.Y(ind(1));
  p2x = SET(no).Point.X(ind(2));
  p2y = SET(no).Point.Y(ind(2));

  %Extract contour and find centre
  if doepi
    x = SET(no).EpiX(:,tf,pos(loop));
    y = SET(no).EpiY(:,tf,pos(loop));
  else
    x = SET(no).EndoX(:,tf,pos(loop));
    y = SET(no).EndoY(:,tf,pos(loop));
  end

  mx = mean(x);
  my = mean(y);

  %Find points on contour that are closest
  dist1 = sqrt((p1x-x).^2+(p1y-y).^2);
  dist2 = sqrt((p2x-x).^2+(p2y-y).^2);

  [~,indp1] = min(dist1(:));
  [~,indp2] = min(dist2(:));

  %Sort them in order
  temp = sort([indp1 indp2]);
  indp1 = temp(1);
  indp2 = temp(2);
  clear temp;

  %Extract contour
  if length(indp1:indp2)/length(x)<0.5
    xr = x(indp1:indp2); %this is then the closest
    yr = y(indp1:indp2);
  else
    xr = [x(indp2:end);x(1:indp1)];
    yr = [y(indp2:end);y(1:indp1)];
  end

  %Extract point
  x = xr(round(length(xr)/2));
  y = yr(round(length(yr)/2));

  %figure(23);
  %imagesc(SET(no).IM(:,:,SET(no).CurrentTimeFrame,SET(no).Point.Z(ind(1))));
  %hold on;
  %plot(yr,xr,'r-');
  %plot(y,x,'bo');
  %plot(p1y,p1x,'ko');
  %plot(p2y,p2x,'ko');
  %hold off;
  %colormap(gray);
  %axis image off;

  %Calculate angle to contour
  alpha = atan2(y-my,x-mx);
  alpha = alpha*180/pi;
  sectorrot = sectorrot-(90+alpha); %add all then divide
end

sectorrot = sectorrot/length(pos);
if isnan(sectorrot)
  sectorrot = 0;
  disp('Failed on rotation calculation');
else
  fprintf('LV rotated %d degrees\n',round(sectorrot));
end

SET(no).SectorRotation = sectorrot;

%---------------------------------------
function rotationfromannotation_Callback
%---------------------------------------
%Finds rotation by looking at RV insertion points.

global DATA

gui = DATA.GUI.Bullseye;

pos = LVrotation('rotationfromannotationhelper',gui.no);
if isempty(pos) %No points for rotation found
  return;
end
%update slider
rotationangle = LVrotation('setrotationangle_helper');
set(gui.handles.rotationslider,'value',rotationangle);

%call to graphically update rotation & bullseye
updatesliceimage;
updateplot;

%-------------------------------
function rotationslider_Callback
%-------------------------------
%Callback for rotation slider. Update slice image.
global DATA
gui = DATA.GUI.Bullseye;

updatesliceimage;

%Enable update button
set(gui.handles.updatepushbutton,...
  'enable','on',...
  'BackgroundColor',DATA.GUISettings.HighlightedButtonColor,...
  'ForegroundColor',DATA.GUISettings.HighlightedButtonText);

%---------------------------
function updatelongaxisimage
%---------------------------
%Update which slices are active in the longaxis image
global DATA SET
gui = DATA.GUI.Bullseye;

if true || ~isempty(gui.laxno)
  %   myset(gui.handles.intersectionlines,'visible','off','color','w');
  %   if 0 %mygetvalue(gui.handles.thissliceonlycheckbox)
  %     slicecolor = 'r';
  %     %If this slice only then set as current slice
  %     if gui.slice >= SET(gui.no).StartSlice && gui.slice <= SET(gui.no).EndSlice
  %       SET(gui.no).CurrentSlice = gui.slice;
  %     end
  %   else
  %     slicecolor = 'c';
  %   end
  myset(gui.handles.intersectionlines(:),'color','w');
  myset(gui.handles.intersectionlines(SET(gui.no).StartSlice:SET(gui.no).EndSlice),'color','y');
  myset(gui.handles.intersectionlines(gui.slice),'color','c');
  %   myset(gui.handles.intersectionlines(SET(gui.no).StartSlice:SET(gui.no).EndSlice),'Visible','on');
end

%------------------------
function updatesliceimage(plotseperate)
%------------------------
%Changed the rotationslider, new slice image
global DATA SET
gui = DATA.GUI.Bullseye;

if nargin == 0
  plotseperate = false;
end

tf = SET(gui.no).CurrentTimeFrame;
SET(gui.no).SectorRotation = mygetvalue(gui.handles.rotationslider);

%Get image
temp = SET(gui.no).IM(:,:,tf,gui.slice);
temp = min(max(temp,0),1);
imsize = size(temp);

if plotseperate == 0
  %clear objects before drawing new objects in the figure
  if isfield(gui.handles,'epicontour')
    try
      delete(gui.handles.epicontour);
    catch
    end
  end
  if isfield(gui.handles,'endocontour')
    try
      delete(gui.handles.endocontour);
    catch
    end
  end
  if isfield(gui.handles,'numtext')
    for loop = 1:length(gui.handles.numtext)
      try
        delete(gui.handles.numtext{loop});
      catch
      end
    end
  end
  if isfield(gui.handles,'line')
    for loop = 1:length(gui.handles.line)
      try
        delete(gui.handles.line{loop});
      catch
      end
    end
  end

  %force true color
  image(calcfunctions('remapuint8',temp,gui.no,calcfunctions('returnmapping',gui.no,true)),...
    'parent',gui.handles.imageaxes);
  axis(gui.handles.imageaxes,'image','off');

else  %plot in seperate window
  figure(23);
  set(23,'Name','LV rotation','numbertitle','off');
  ax = gca;
  image(calcfunctions('remapuint8',temp,gui.no,calcfunctions('returnmapping',gui.no,true)),...
    'parent',ax);
  axis(ax,'image','off');
  cmap = gray(256);
  colormap(ax,cmap);
  title(dprintf('LV rotation'));
end

%Update longaxisimage
updatelongaxisimage;

%Upsample model
if not(DATA.Pref.RadialProfiles == DATA.NumPoints)
  [endox,endoy] = calcfunctions('resamplemodel',SET(gui.no).EndoX(:,tf,gui.slice),SET(gui.no).EndoY(:,tf,gui.slice),DATA.Pref.RadialProfiles);
  if ~isempty(SET(gui.no).EpiX)
    [epix,epiy] = calcfunctions('resamplemodel',SET(gui.no).EpiX(:,tf,gui.slice), SET(gui.no).EpiY(:,tf,gui.slice),DATA.Pref.RadialProfiles);
  end
else
  endox = SET(gui.no).EndoX(:,tf,gui.slice);
  endoy = SET(gui.no).EndoY(:,tf,gui.slice);
  if ~isempty(SET(gui.no).EpiX)
    epix = SET(gui.no).EpiX(:,tf,gui.slice);
    epiy = SET(gui.no).EpiY(:,tf,gui.slice);
  end
end
if isempty(SET(gui.no).EpiX)
  epix = NaN;
  epiy = NaN;
end

%Plot contours
if ~plotseperate
  hold(gui.handles.imageaxes,'on');
  gui.handles.endocontour = plot(gui.handles.imageaxes,endoy,endox,'r-');
  if DATA.Pref.LineWidth==0
    set(gui.handles.endocontour,'linewidth','visible','off');
  else
    set(gui.handles.endocontour,'linewidth',DATA.Pref.LineWidth);
  end
  if ~isempty(SET(gui.no).EpiX)
    gui.handles.epicontour = plot(gui.handles.imageaxes,epiy,epix,'g-');
    if DATA.Pref.LineWidth==0
      set(gui.handles.epicontour,'linewidth','visible','off');
    else
      set(gui.handles.epicontour,'linewidth',DATA.Pref.LineWidth);
    end
  end

else %plot in seperate window
  hold(ax,'on');
  endocontour = plot(ax,endoy,endox,'r-');
  set(endocontour,'linewidth',DATA.Pref.LineWidth);
  if ~isempty(SET(gui.no).EpiX)
    epicontour = plot(ax,epiy,epix,'g-');
    set(epicontour,'linewidth',DATA.Pref.LineWidth);
  end
end


%Plot sectors
isahaplot = (get(gui.handles.aharadiobutton,'value') || mygetvalue(gui.handles.aha16radiobutton));
hold(gui.handles.imageaxes,'on');
if isnan(epix(1))
  %--- Only endo exist => draw only endo
  %get positions.
  [gui.meanx,gui.meany,gui.sectors] = calcfunctions('findmeaninsectorslice','endo',DATA.Pref.RadialProfiles,...
    tf,gui.slice,gui.numsectors,gui.no);
  for loop=1:gui.numsectors
    xpos = endoy(gui.sectors(loop));
    ypos = endox(gui.sectors(loop));
    xposm = ...
      0.5*endoy(gui.sectors(loop))+...
      0.5*endoy(gui.sectors(loop+1));
    yposm = ...
      0.5*endox(gui.sectors(loop))+...
      0.5*endox(gui.sectors(loop+1));

    if ~plotseperate
      if (loop==1)
        %Draw an extra long line
        h = plot(gui.handles.imageaxes,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,1.5*xpos-0.5*gui.meany))], ...
          [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,1.5*ypos-0.5*gui.meany))],'y-');
        if isahaplot
          % stop drawing sectors
          break;
        end
      else
        h = plot(gui.handles.imageaxes,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,xpos))], ...
          [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,ypos))],'w-');
      end
      if gui.numsectors>1 && gui.numsectors<=20
        gui.handles.numtext{loop} = text(xposm,yposm,sprintf('%d',gui.numsectors+1-loop),'parent',gui.handles.imageaxes);
        set(gui.handles.numtext{loop},'color','y','fontsize',12);
      end
      if DATA.Pref.LineWidth==0
        set(h,'visible','off');
      else
        set(h,'linewidth',DATA.Pref.LineWidth);
      end
    else %plot in seperate window
      if (loop==1)
        %Draw an extra long line
        h = plot(ax,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,1.5*xpos-0.5*gui.meany))], ...
          [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,1.5*ypos-0.5*gui.meany))],'y-');
        set(h,'linewidth',DATA.Pref.LineWidth);
      end
    end
  end

else
  %--- Epi exists => draw both
  %get positions.
  [gui.meanx,gui.meany,gui.sectors] = calcfunctions('findmeaninsectorslice','epi',DATA.Pref.RadialProfiles,...
    tf,gui.slice,gui.numsectors,gui.no);
  for loop=1:gui.numsectors
    xpos = epiy(gui.sectors(loop));

    ypos = epix(gui.sectors(loop));
    xposm = ...
      0.5*epiy(gui.sectors(loop))+...
      0.5*epiy(gui.sectors(loop+1));
    yposm = ...
      0.5*epix(gui.sectors(loop))+...
      0.5*epix(gui.sectors(loop+1));

    if ~plotseperate
      if (loop==1)
        %Draw an extra long line
        gui.handles.line{loop} = plot(gui.handles.imageaxes,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,1.5*xpos-0.5*gui.meany))], ...
          [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,1.5*ypos-0.5*gui.meanx))],'y-');
        if isahaplot
          % stop drawing sectors
          break;
        end
      else
        gui.handles.line{loop} = plot(gui.handles.imageaxes,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,xpos))], ...
          [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,ypos))],'w-');
      end
      if gui.numsectors>1 && gui.numsectors<=20
        gui.handles.numtext{loop} = text(xposm,yposm,sprintf('%d',gui.numsectors+1-loop),'parent',gui.handles.imageaxes);
        set(gui.handles.numtext{loop},'color','y','fontsize',12);
      end
      if DATA.Pref.LineWidth==0
        set(gui.handles.line{loop},'visible','off');
      else
        set(gui.handles.line{loop},'linewidth',DATA.Pref.LineWidth);
      end

    else %plot in seperate window
      if (loop==1)
        %Draw an extra long line
        line = plot(ax,[min(imsize(2),max(1,gui.meany)) min(imsize(2),max(1,1.5*xpos-0.5*gui.meany))], ...
          [min(imsize(1),max(1,gui.meanx)) min(imsize(1),max(1,1.5*ypos-0.5*gui.meanx))],'y-');
        set(line,'linewidth',DATA.Pref.LineWidth);
      end
    end

  end
end
if ~plotseperate
  hold(gui.handles.imageaxes,'off');
else
  hold(ax,'off');
end

%----------------------
function close_Callback
%----------------------
%Properly close the GUI.

global DATA
try
  close(DATA.GUI.Bullseye);  %close the flow gui
catch
  close(gcbf);
end
DATA.GUI.Bullseye= [];

%-----------------------
function drawsectorimage
%-----------------------
%Draw image of sector division in new figure.

global DATA SET
gui = DATA.GUI.Bullseye;

figure(325);
%Changed the rotationslider, new image
tf = SET(gui.no).CurrentTimeFrame;
SET(gui.no).SectorRotation = mygetvalue(gui.handles.rotationslider);

%Plot image
temp = SET(gui.no).IM(:,:,tf,gui.slice);
temp = min(max(temp,0),1);

%force true color
image(calcfunctions('remapuint8',temp,gui.no,calcfunctions('returnmapping',gui.no,true)));
axis('image','off');

%Upsample model
if not(DATA.Pref.RadialProfiles==DATA.NumPoints)
  [endox,endoy] = calcfunctions('resamplemodel',SET(gui.no).EndoX(:,tf,gui.slice),SET(gui.no).EndoY(:,tf,gui.slice),DATA.Pref.RadialProfiles);
  if ~isempty(SET(gui.no).EpiX)
    [epix,epiy] = calcfunctions('resamplemodel',SET(gui.no).EpiX(:,tf,gui.slice), SET(gui.no).EpiY(:,tf,gui.slice),DATA.Pref.RadialProfiles);
  end
else
  endox = SET(gui.no).EndoX(:,tf,gui.slice);
  endoy = SET(gui.no).EndoY(:,tf,gui.slice);
  if ~isempty(SET(gui.no).EpiX)
    epix = SET(gui.no).EpiX(:,tf,gui.slice);
    epiy = SET(gui.no).EpiY(:,tf,gui.slice);
  end
end
if isempty(SET(gui.no).EpiX)
  epix = NaN;
  epiy = NaN;
end

%Plot contours
hold on;
h = plot(endoy,endox,'r-');
if DATA.Pref.LineWidth==0
  set(h,'linewidth','visible','off');
else
  set(h,'linewidth',DATA.Pref.LineWidth);
end

if ~isempty(SET(gui.no).EpiX)
  h = plot(epiy,epix,'g-');
  if DATA.Pref.LineWidth==0
    set(h,'linewidth','visible','off');
  else
    set(h,'linewidth',DATA.Pref.LineWidth);
  end
end

%Plot sectors
hold on;
if isnan(epix(1))
  %--- Only endo exist => draw only endo
  %get positions.
  [gui.meanx,gui.meany,gui.sectors] = calcfunctions('findmeaninsectorslice','endo',DATA.Pref.RadialProfiles,...
    tf,gui.slice,gui.numsectors,gui.no);
  for loop=1:gui.numsectors
    xpos = endoy(gui.sectors(loop));
    ypos = endox(gui.sectors(loop));
    xposm = ...
      0.5*endoy(gui.sectors(loop))+...
      0.5*endoy(gui.sectors(loop+1));
    yposm = ...
      0.5*endox(gui.sectors(loop))+...
      0.5*endox(gui.sectors(loop+1));
    if (loop==1)
      h = plot([gui.meany 1.5*xpos-0.5*gui.meany],[gui.meanx 1.5*ypos-0.5*gui.meany],'y-');
    else
      h = plot([gui.meany xpos],[gui.meanx ypos],'y-');
    end
    if DATA.Pref.LineWidth==0
      set(h,'linewidth','visible','off');
    else
      set(h,'linewidth',DATA.Pref.LineWidth);
    end
    if gui.numsectors<=20
      h = text(xposm,yposm,sprintf('%d',gui.numsectors+1-loop));
      set(h,'color','y','fontsize',12);
    end
  end

else
  %--- Epi exists => draw both
  %get positions.
  [gui.meanx,gui.meany,gui.sectors] = calcfunctions('findmeaninsectorslice','epi',DATA.Pref.RadialProfiles,...
    tf,gui.slice,gui.numsectors,gui.no);
  for loop=1:gui.numsectors
    xpos = epiy(gui.sectors(loop));
    ypos = epix(gui.sectors(loop));
    xposm = ...
      0.5*epiy(gui.sectors(loop))+...
      0.5*epiy(gui.sectors(loop+1));
    yposm = ...
      0.5*epix(gui.sectors(loop))+...
      0.5*epix(gui.sectors(loop+1));
    if (loop==1)
      h = plot([gui.meany 1.5*xpos-0.5*gui.meany],[gui.meanx 1.5*ypos-0.5*gui.meanx],'y-');
    else
      h = plot([gui.meany xpos],[gui.meanx ypos],'y-');
    end
    if DATA.Pref.LineWidth==0
      set(h,'linewidth','visible','off');
    else
      set(h,'linewidth',DATA.Pref.LineWidth);
    end
    if gui.numsectors<=20
      h = text(xposm,yposm,sprintf('%d',gui.numsectors+1-loop));
      set(h,'color','y','fontsize',12);
    end
  end
end
hold off;

%------------------------------
function thissliceonly_Callback
%------------------------------
%Callback for this slice only checkbox
updatelongaxisimage;
updateplot;

%--------------------------------
function bullseyelistbox_Callback
%--------------------------------
%Callback for bullseye listbox. Updates plot
fromlistbox = true;
updateplot(fromlistbox);

%--------------------------------
function colormaplistbox_Callback
%--------------------------------
%Callback for colormap listbox. Updates plot
updateplot;

%-----------------------------
function invertcolors_Callback
%-----------------------------
%Callback for invert colors checkbox. Updates plot
updateplot;

%----------------------
function nedit_Callback
%----------------------
%Callback for edit to change value of n. Updates plot
updateplot;

%------------------------
function minedit_Callback
%------------------------
%Callback for edit to change min value. Updates plot
updateplot;

%------------------------
function maxedit_Callback
%------------------------
%Callback for edit to change max value. Updates plot
updateplot;

%-------------------------------
function separatewindow_Callback
%-------------------------------
%Callback for separate window checkbox. Updates plot
updateplot;

%---------------------------------
function updatepushbutton_Callback
%---------------------------------
%Callback for update pushbutton. Updates plot
updateplot;

%------------------
function updateplot(fromlistbox)
%------------------
%Plot different type of data. Calculate / retrieve
%data, and perform graphical update. This is a main workhorse.

global DATA SET
gui = DATA.GUI.Bullseye;

if nargin == 0
  fromlistbox = false;
end

%Grey out update button
set(gui.handles.updatepushbutton,'enable','off',...
  'BackgroundColor',DATA.GUISettings.BackgroundColor,'ForegroundColor',DATA.GUISettings.ForegroundColor);

%Get String
temp = mygetlistbox(gui.handles.bullseyelistbox);
tempstri = get(gui.handles.bullseyelistbox,'String');
gui.parameter = tempstri{temp};

ind = false(1,SET(gui.no).ZSize);
if 1 %~mygetvalue(gui.handles.thissliceonlycheckbox)
  ind(SET(gui.no).StartSlice:SET(gui.no).EndSlice) = true;
else
  ind(gui.slice) = true;
end
if sum(ind) < 3
  set([gui.handles.aharadiobutton, gui.handles.aha16radiobutton],'enable','off');
  if mygetvalue(gui.handles.aharadiobutton) || mygetvalue(gui.handles.aha16radiobutton)
    set(gui.handles.sectorradiobutton,'value',1);
    plotmethodpanel_SelectionChange
    return % returning since the method above includes the whole update
  end
else
  set([gui.handles.aharadiobutton, gui.handles.aha16radiobutton],'enable','on');
end

flipx = false;
valuetype = 'mean'; %Could be mean value or a max value. Used for subsequent
%handling such as calculating AHA bullseye.

m = zeros(gui.numsectors,sum(ind));%1);
outdata = m;
outunit = 'Nothing';
tf = 1;
outtf = [];

engparameter = translation.dictionary(tempstri{temp},'English',DATA.Pref.Language); %ensure that it is English
switch lower(engparameter)
  case 'maximal expansion velocity (temporal max)'
    radvel = calcfunctions('calcradialvelocity',gui.no);
    meanradvel = calcfunctions('findmeaninsector','endo',radvel,find(ind),gui.numsectors,gui.no);
    outdata = max(meanradvel,[],3);
    tf = SET(gui.no).EDT;
    outunit = 'cm/s';
  case 'maximal contraction velocity (temporal max)'
    radvel = calcfunctions('calcradialvelocity',gui.no);
    meanradvel = calcfunctions('findmeaninsector','endo',radvel,find(ind),gui.numsectors,gui.no);
    outdata = -min(meanradvel,[],3);
    tf = SET(gui.no).EDT;
    outunit = 'cm/s';
  case 'expansion velocity at pfr'
    radvel = calcfunctions('calcendoradius',gui.no);
    meanradvel = calcfunctions('findmeaninsector','endo',radvel,find(ind),gui.numsectors,gui.no);
    outdata = -meanradvel(:,:,SET(gui.no).PFRT);
    tf = SET(gui.no).PFRT;
    outunit = 'cm/s';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'contraction velocity at per'
    radvel = calcfunctions('calcendoradius',gui.no);
    meanradvel = calcfunctions('findmeaninsector','endo',radvel,find(ind),gui.numsectors,gui.no);
    outdata = meanradvel(:,:,SET(gui.no).PERT);
    tf = SET(gui.no).PERT;
    outunit = 'cm/s';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'maximal wallthickness (temporal max)'
    wallthickness = calcfunctions('calcwallthickness',gui.numsectors,gui.no);
    outdata = squeeze(max(wallthickness(:,ind,:),[],3));
    tf = SET(gui.no).EDT;
    outunit = 'mm';
  case 'mean ed wallthickness'
    wallthickness = calcfunctions('calcwallthickness',gui.numsectors,gui.no);
    outdata = wallthickness(:,ind,SET(gui.no).EDT);
    tf = SET(gui.no).EDT;
    outunit = 'mm';
  case 'mean es wallthickness'
    wallthickness = calcfunctions('calcwallthickness',gui.numsectors,gui.no);
    outdata = wallthickness(:,ind,SET(gui.no).EST);
    tf = SET(gui.no).EST;
    outunit = 'mm';
  case 'fractional wallthickening'
    if isequal(SET(gui.no).EDT,SET(gui.no).EST)
      mywarning('Warning, end-diastole occurs at the same time as end-systole. Use autodetect under edit menu.',DATA.GUI.Segment);
    end
    wallthickness = calcfunctions('calcwallthickness',gui.numsectors,gui.no);
    outdata = (wallthickness(:,ind,SET(gui.no).EST)-wallthickness(:,ind,SET(gui.no).EDT))./wallthickness(:,ind,SET(gui.no).EDT);
    outdata = outdata*100;
    tf = [SET(gui.no).EDT SET(gui.no).EST];
    outunit = '%';
  case 'wallthickening'
    wallthickness = calcfunctions('calcwallthickness',gui.numsectors,gui.no);
    outdata = wallthickness(:,ind,SET(gui.no).EST)-wallthickness(:,ind,SET(gui.no).EDT);
    tf = [SET(gui.no).EDT SET(gui.no).EST];
    outunit = 'mm';
  case 'myocard volume'
    valuetype = 'sum';
    tf = SET(gui.no).CurrentTimeFrame;
    outdata = calcfunctions('calcmyocardvolume',gui.numsectors,gui.no);
    outdata = outdata(:,ind);
    outunit = 'ml';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'myocard intensity'
    tf = SET(gui.no).CurrentTimeFrame;
    outdata = calcfunctions('calcintensityanddefect',SET(gui.no).IM,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no, ...
      SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],[SET(gui.no).ResolutionX SET(gui.no).ResolutionY]);

    outunit = '';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'myocard intensity (unnormalized values)'
    tf = SET(gui.no).CurrentTimeFrame;
    outdata = calcfunctions('calcintensityanddefect',SET(gui.no).IM,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no, ...
      SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],[SET(gui.no).ResolutionX SET(gui.no).ResolutionY]);

    outdata = calcfunctions('calctruedata',outdata,gui.no);
    outunit = '';
    if SET(gui.no).TSize > 1, outtf = tf; end

  case 'mapping values'
    %- numsectors:                   number of sectors
    %- nprofiles:                         number of profiles
    %- pos:                                  slices to use
    %- endox,endoy,epix,epiy: myocardial borders
    %sz:                                        size of image stack
    %resolution:                         resolution of image stack
    %defect:                                (EH: I do not know what it is, works if empty)
    %numwidth:                        specify if it should calculate several sectors across wallthickness.
    %timeresolved:                   true if timeresolvde

    if isempty(gui.Endocardium_percent) || fromlistbox
      s = [];
      s.Endocardium_percent = 20;
      s.Epicardium_percent = 20;
      [s,ok] = inputstruct(s,'Selected coverage.');
      if ~ok
        %       myfailed('Aborted.');
        return;
      end
      gui.Endocardium_percent = s.Endocardium_percent;
      gui.Epicardium_percent = s.Epicardium_percent;
    else
      s.Endocardium_percent = gui.Endocardium_percent;
      s.Epicardium_percent = gui.Epicardium_percent;
    end


    endo = s.Endocardium_percent/100;
    epi = s.Epicardium_percent/100;

    tf = SET(gui.no).CurrentTimeFrame;
    outdata = calcfunctions('calcintensityanddefect',...
      SET(gui.no).IM,...
      tf,...
      DATA.Pref.RadialProfiles,...
      gui.numsectors,...
      DATA.Pref.RadialProfiles, ...
      find(ind),...
      gui.no, ...
      SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],...
      [SET(gui.no).ResolutionX SET(gui.no).ResolutionY],...
      [],...
      [endo 1-epi],...
      false);

    outdata = calcfunctions('calctruedata',outdata,gui.no);
    outunit = '';
    if SET(gui.no).TSize > 1, outtf = tf; end

  case 'scar transmurality area based'
    if isempty(SET(gui.no).Scar)
      myfailed('No infarct data available, use report myocardial intensity before.',DATA.GUI.Segment);
      return;
    end
    tf = SET(gui.no).CurrentTimeFrame;
    [~,outdata] = calcfunctions('calcintensityanddefect',SET(gui.no).IM,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no, ...
      SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],[SET(gui.no).ResolutionX SET(gui.no).ResolutionY],SET(gui.no).Scar.Result,1,false);
    outunit = '%';
  case 'scar transmurality line based'
    if isempty(SET(gui.no).Scar)
      myfailed('No infarct data available, use report myocardial intensity before.',DATA.GUI.Segment);
      return;
    end
    tf = SET(gui.no).CurrentTimeFrame;
    outdata = viability('calctransmuralityline',gui.numsectors,gui.no);
    outunit = '%';
  case 'max scar transmurality (spatial max)'
    if isempty(SET(gui.no).Scar)
      myfailed('No infarct data available, use report myocardial intensity before.',DATA.GUI.Segment);
      return;
    end
    valuetype = 'max';
    tf = SET(gui.no).CurrentTimeFrame;
    [~,outdata] = viability('calctransmuralityline',gui.numsectors,gui.no);
    outunit = '%';
  case 'scar transmurality area based and myocardial intensity'
    if isempty(SET(gui.no).Scar)
      myfailed('No Scar data available.',DATA.GUI.Segment);
      return;
    end
    tf = SET(gui.no).CurrentTimeFrame;
    outdata = calcfunctions('calcintensityanddefect',SET(gui.no).IM,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no,SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],[SET(gui.no).ResolutionX SET(gui.no).ResolutionY],SET(gui.no).Scar.Result);
    [~,defect] = calcfunctions('calcintensityanddefect',SET(gui.no).IM,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no,SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],[SET(gui.no).ResolutionX SET(gui.no).ResolutionY],SET(gui.no).Scar.Result,1,false);
    %     %fill out the nan values in the outflow tract region by smoothing
    %     [filloutrow,filloutcol] = find(isnan(outdata));
    %     for filloop = 1:length(filloutrow)
    %       [row,col] = find(~isnan(outdata(filloutrow(filloop),filloutcol(filloop):end)),1,'first');
    %       if isempty(row)
    %         [~,col] = find(~isnan(outdata(filloutrow(filloop),1:filloutcol(filloop))),1,'last');
    %         if ~isempty(col)
    %           outdata(filloutrow(filloop),filloutcol(filloop)) = outdata(filloutrow(filloop),col);
    %         end
    %       else
    %         outdata(filloutrow(filloop),filloutcol(filloop)) = outdata(filloutrow(filloop),col+filloutcol(filloop)-1);
    %       end
    %     end
    defect(defect < 20) = 1;
    defect(defect > 50) = 0;
    defect(defect > 1) = 0.5;
    outdata = outdata.*defect;
    outunit = '';
  case 'weighted infarct transmurality'
    if isempty(SET(gui.no).Scar)
      myfailed('No Scar data available.',DATA.GUI.Segment);
      return;
    end
    %Call viabilityweight to get weighting of infarct
    tf = SET(gui.no).CurrentTimeFrame;
    [~,infarctweightmap] = viability('viabilityweight',gui.no);
    if isempty(infarctweightmap)
      infarctweightmap = ones(size(SET(gui.no).Scar.Result));
    end
    meaninfarctweight = calcfunctions('calcintensityanddefect', ...
      infarctweightmap.*SET(gui.no).Scar.Result,tf,DATA.Pref.RadialProfiles, ...
      gui.numsectors,DATA.Pref.RadialProfiles,find(ind),gui.no, ...
      SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize], ...
      [SET(gui.no).ResolutionX SET(gui.no).ResolutionY],[],1,false);
    outdata = 100*meaninfarctweight;
    outunit = '%';
  case 'grayzone transmurality area based'
    if isempty(SET(gui.no).Scar)
      myfailed('No infarct data available, use report myocardial intensity before.',DATA.GUI.Segment);
      return;
    end
    mapsz = size(SET(gui.no).Scar.GreyZone.map);
    rescale = @(x)(0.5+mapsz(1)/SET(gui.no).XSize*(x-0.5));
    im = zeros(mapsz); %imresize(SET(gui.no).IM,mapsz,'bilinear');
    tf = SET(gui.no).CurrentTimeFrame;
    [~,outdata] = calcfunctions('calcintensityanddefect',im,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no, ...
      rescale(SET(gui.no).EndoX),rescale(SET(gui.no).EndoY), ...
      rescale(SET(gui.no).EpiX),rescale(SET(gui.no).EpiY), ...
      mapsz(1:2),[SET(gui.no).ResolutionX SET(gui.no).ResolutionY],SET(gui.no).Scar.GreyZone.map==1,1,false);
    outunit = '%';
  case 'core transmurality area based'
    if isempty(SET(gui.no).Scar)
      myfailed('No infarct data available, use report myocardial intensity before.',DATA.GUI.Segment);
      return;
    end
    mapsz = size(SET(gui.no).Scar.GreyZone.map);
    rescale = @(x)(0.5+mapsz(1)/SET(gui.no).XSize*(x-0.5));
    im = zeros(mapsz); %imresize(SET(gui.no).IM,mapsz,'bilinear');
    tf = SET(gui.no).CurrentTimeFrame;
    [~,outdata] = calcfunctions('calcintensityanddefect',im,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no, ...
      rescale(SET(gui.no).EndoX),rescale(SET(gui.no).EndoY), ...
      rescale(SET(gui.no).EpiX),rescale(SET(gui.no).EpiY), ...
      mapsz(1:2),[SET(gui.no).ResolutionX SET(gui.no).ResolutionY],SET(gui.no).Scar.GreyZone.map==2,1,false);
    outunit = '%';
  case 'mar transmurality area based'
    if isempty(SET(gui.no).MaR)
      myfailed('No MaR data available.',DATA.GUI.Segment);
      return;
    end
    tf = SET(gui.no).CurrentTimeFrame;
    [~,outdata] = calcfunctions('calcintensityanddefect',SET(gui.no).IM,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no, ...
      SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],[SET(gui.no).ResolutionX SET(gui.no).ResolutionY],SET(gui.no).MaR.Result);
    outunit = '%';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'mar transmurality line based'
    if isempty(SET(gui.no).MaR)
      myfailed('No MaR data available.',DATA.GUI.Segment);
      return;
    end
    tf = SET(gui.no).CurrentTimeFrame;
    outdata = mar('calctransmuralityline',gui.numsectors,gui.no);
    outunit = '%';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'max mar transmurality (spatial max)'
    if isempty(SET(gui.no).MaR)
      myfailed('No MaR data available.',DATA.GUI.Segment);
      return;
    end
    valuetype = 'max';
    tf = SET(gui.no).CurrentTimeFrame;
    [~,outdata] = mar('calctransmuralityline',gui.numsectors,gui.no);
    outunit = '%';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'mar transmurality area based and myocardial intensity'
    if isempty(SET(gui.no).MaR)
      myfailed('No MaR data available.',DATA.GUI.Segment);
      return;
    end
    tf = SET(gui.no).CurrentTimeFrame;
    [outdata,defect] = calcfunctions('calcintensityanddefect',SET(gui.no).IM,tf,DATA.Pref.RadialProfiles,gui.numsectors,DATA.Pref.RadialProfiles, ...
      find(ind),gui.no, ...
      SET(gui.no).EndoX,SET(gui.no).EndoY,SET(gui.no).EpiX,SET(gui.no).EpiY, ...
      [SET(gui.no).XSize SET(gui.no).YSize],[SET(gui.no).ResolutionX SET(gui.no).ResolutionY],SET(gui.no).MaR.Result);
    %     %fill out the nan values in the outflow tract region by smoothing
    %     [filloutrow,filloutcol] = find(isnan(outdata));
    %     for filloop = 1:length(filloutrow)
    %       [row,col] = find(~isnan(outdata(filloutrow(filloop),filloutcol(filloop):end)),1,'first');
    %       if isempty(row)
    %         [~,col] = find(~isnan(outdata(filloutrow(filloop),1:filloutcol(filloop))),1,'last');
    %         if ~isempty(col)
    %           outdata(filloutrow(filloop),filloutcol(filloop)) = outdata(filloutrow(filloop),col);
    %         end
    %       else
    %         outdata(filloutrow(filloop),filloutcol(filloop)) = outdata(filloutrow(filloop),col+filloutcol(filloop)-1);
    %       end
    %     end
    defect(defect < 20) = 1;
    defect(defect > 50) = 0;
    defect(defect > 1) = 0.5;
    outdata = outdata.*defect;
    outunit = '';
    if SET(gui.no).TSize > 1, outtf = tf; end
  case 'clipboard data'
    outdata = clipboard('paste');
    [outdata,ok] = str2num(outdata);
    outdata = outdata';
    flipx = true;
    if not(ok)
      myfailed('Could not import clipboard data.',DATA.GUI.Segment);
      return;
    end
    tf = 1;
    outunit = '';
  otherwise
    gui.parameter = '';
    outunit = '';
end

%Determine max and min value
[maxv,ok] = str2num(mygetedit(gui.handles.maxedit));
if not(ok) || isempty(maxv)
  maxv = max(outdata(:));
end
[minv,ok] = str2num(mygetedit(gui.handles.minedit));
if not(ok) || isempty(minv)
  minv = min(outdata(:));
end

%Store
gui.outdata = outdata;
gui.outunit = outunit;
gui.outtf = outtf;

%Find axes to plot in
isseparetewindow = false;
if get(gui.handles.separatewindowcheckbox,'value')
  separatewindowfig = figure(22);
  setupicon(separatewindowfig);
  set(separatewindowfig,'Name','Bullseye plot','numbertitle','off');
  set(separatewindowfig,'Color',DATA.GUISettings.BackgroundColor);
  tempax = gca;
  isseparetewindow = true;
else
  tempax = gui.handles.bullseyeaxes;
end

%Colormap
n = 256;
selectedcolormap = gui.colormaplist{mygetlistbox(gui.handles.colormaplistbox)};
switch selectedcolormap
  case 'Jet'
    cmap = jet(n);
  case 'Hot'
    cmap = hot(n);
  case 'HSV'
    cmap = hsv(n);
  case 'SPECT'
    cmap = spect(n);
  case dprintf('Gray')
    cmap = gray(n);
  case 'Gadgetron'
    temp = load('colormapgadgetron');
    cmapgt = temp.colormapgadgetron;
    cmap = zeros(n,3);
    cmap(:,1) = interp1(1:length(cmapgt),cmapgt(:,1),linspace(1,length(cmapgt),n))';
    cmap(:,2) = interp1(1:length(cmapgt),cmapgt(:,2),linspace(1,length(cmapgt),n))';
    cmap(:,3) = interp1(1:length(cmapgt),cmapgt(:,3),linspace(1,length(cmapgt),n))';
end

if get(gui.handles.invertcolorscheckbox,'value')
  cmap = flipud(cmap);
end
colormap(tempax,cmap);

%get number of points in image
[res,ok] = str2num(mygetedit(gui.handles.nedit));
if not(ok)
  mywarning('Not a valid number for n.',DATA.GUI.Segment);
  res = 200;
end

%Determine if smooth bulleye version or not.
gui.ahaoutdata = [];
if mygetvalue(gui.handles.sectorradiobutton)
  bullseye(gui.outdata,tempax,res,gui.volumeconsistent,gui.no,tf);
elseif mygetvalue(gui.handles.smoothradiobutton)
  % check if data is double, scatteredInterpolant expects double
  if ~isa(gui.outdata,'double')
    gui.outdata = double(gui.outdata);
  end
  bullseye2(gui.outdata,tempax,res,flipx,gui.volumeconsistent,gui.no);
elseif mygetvalue(gui.handles.aharadiobutton) || mygetvalue(gui.handles.aha16radiobutton)
  gui.ahaoutdata = bullseyeaha(gui.outdata,tempax,res,valuetype);
  %Determine max and min value
  [maxv,ok] = str2num(mygetedit(gui.handles.maxedit));
  if not(ok) || max(gui.ahaoutdata(:)) > maxv
    maxv = round(max(gui.ahaoutdata(:)),4);
  end
  [minv,ok] = str2num(mygetedit(gui.handles.minedit));
  if not(ok) || min(gui.ahaoutdata(:)) < minv
    minv = round(min(gui.ahaoutdata(:)),4);
  end
end
colormap(tempax,cmap);

%Adjust max/min
if ~isempty(minv) && ~isempty(maxv) && not(isnan(minv)) && not(isnan(maxv)) && (maxv>minv)
  set(tempax,'clim',[minv maxv]);
%   set(gui.handles.maxedit,'String',maxv);
%   set(gui.handles.minedit,'String',minv);
end

%Add title
if isempty(gui.outunit) && isempty(gui.outtf)
  titlestr = sprintf('%s',gui.parameter);
elseif isempty(gui.outunit) && ~isempty(gui.outtf)
  titlestr = sprintf('%s %s',gui.parameter,dprintf('in time frame %d',gui.outtf));
elseif ~isempty(gui.outunit) && isempty(gui.outtf)
  titlestr = sprintf('%s [%s]',gui.parameter,gui.outunit);
else
  titlestr = sprintf('%s [%s] %s',gui.parameter,gui.outunit,dprintf('in time frame %d',gui.outtf));
end
title(tempax,titlestr,'Color',DATA.GUISettings.ForegroundColor);
if isseparetewindow
  tempax.Title.FontSize = 14;
end

%Add colorbar
gui.handles.colorbar = colorbar('peer',tempax);
if isseparetewindow
  fontweight = 'normal';
  fontsize = 12;
  linewidth = 1;
else
  fontweight = 'bold';
  fontsize = 10;
  linewidth = 1.5;
end
set(gui.handles.colorbar,'Color',DATA.GUISettings.ForegroundColor, ...
  'LineWidth',linewidth,'FontSize',fontsize,'FontWeight',fontweight);

%----------------------------
function sliceslider_Callback
%----------------------------
%Callback for slider to toggle slice
global DATA SET
gui = DATA.GUI.Bullseye;

gui.slice = SET(gui.no).ZSize-round(mygetvalue(gui.handles.sliceslider))+1;
set(gui.handles.sliceslider,'Value',round(mygetvalue(gui.handles.sliceslider)));
set(gui.handles.slicetext,'String',sprintf('%s %d',dprintf('Slice'),gui.slice));
updatesliceimage;
updatelongaxisimage
if 0%mygetvalue(gui.handles.thissliceonlycheckbox)
  updateplot;
end

%---------------------------------------------------
function [varargout] = bullseyeaha(m,ax,n,valuetype,linecolor,includingvalues)
%---------------------------------------------------
%Calculate and/or plot AHA 17 segment model.
%- m is a matrix in polar coordinates. It could also be a vector of 17
%  segments. Please note in such cases then the order is not the same as
%  the standard AHA numbering. The names are given by the function
%  aha17nameandpos.
%
%- ax is axis to plot it in.
%- n is the number of pixels
%- valuetype is how to treat the merging of data into sectors; mean, sum,
%  max.

global DATA
gui = DATA.GUI.Bullseye;

if nargin<4
  valuetype = 'mean'; %default;
end

if nargin<5
  linecolor = DATA.GUISettings.ForegroundColor;
end

if nargin<6
  includingvalues=true;
end

varargout = cell(1,nargout);
if ~isempty(gui)
  is16segments = mygetvalue(gui.handles.aha16radiobutton);
else
  is16segments = false; %17-segment model by default, for report for example
end

if numel(m)==17
  %This must be just the 17 segments. Added this EH:
  v = m;
else
  %Calculate v from the polar matrix m

  if size(m,1)~=24
    myfailed('Expected 24 sectors as input size.',DATA.GUI.Segment);
    return;
  end

  numslices = size(m,2);
  m = fliplr(m); %Basal slice is first column
  %Reshape the transmurality by replicating each slice into three
  newind = 1:numslices;
  newind = repmat(newind,3,1);
  newind = newind(:)';
  m = m(:,newind);

  if numslices == 3 || is16segments
    firstslice = 1;
    ratios = [1/3, 2/3];
    mlength=size(m,2);
  else
    firstslice = 4;
    ratios = [1/3, 2/3];
    mlength=size(m,2)-firstslice+1;
  end

  %Calculate mean
  switch valuetype
    case 'mean'
      tempm = zeros(size(m,1),4);
      if numslices == 3 || is16segments
        tempm(:,1) = NaN;  %apex
      else
        tempm(:,1) = mean(m(:,1:3),2);  %apex %for 17-segment model, apex is always the most apical slice
      end
      tempm(:,2) = mean(m(:,firstslice:round(mlength*ratios(1))+firstslice-1),2);  %apical part
      tempm(:,3) = mean(m(:,firstslice-1+round(mlength*ratios(1))+1:firstslice-1+round(mlength*ratios(2))),2);  %mid-ventricular part
      tempm(:,4) = mean(m(:,firstslice-1+round(mlength*ratios(2))+1:end),2);  %basal part
      m = tempm;

      %Make AHA
      v = [];
      temp = m(:,4);
      temp = sum(reshape(temp(:),[4 6]))/4;  %last on reshape is how many sectors to get
      v = [v;temp(:)]; %Add basal slices
      temp = m(:,3);
      temp=sum(reshape(temp(:),[4 6]))/4;
      v = [v;temp(:)]; %Add mid slices
      temp = m(:,2);
      %shift it
      temp = [temp(23:24);temp(1:22)];
      temp=sum(reshape(temp(:),[6 4]))/6;
      v = [v;temp(:)]; %Add apical slices
      v = [v;mean(m(:,1))];  %Add apex
    case 'sum'
      %Calculate sum
      tempm = zeros(size(m,1),4);
      if numslices == 3 || is16segments
        tempm(:,1) = NaN;  %apex
      else
        tempm(:,1) = sum(m(:,1:3),2);  %apex %for 17-segment model, apex is always the most apical slice
      end
      tempm(:,2) = sum(m(:,firstslice:round(size(m,2)*ratios(1))),2);  %apical part
      tempm(:,3) = sum(m(:,round(size(m,2)*ratios(1))+1:round(size(m,2)*ratios(2))),2);  %mid-ventricular part
      tempm(:,4) = sum(m(:,round(size(m,2)*ratios(2))+1:end),2);  %basal part
      m = tempm;
      m = m/3; %since we copied each slice 3 times
      %Make AHA
      v = [];
      temp = m(:,4);
      temp=sum(reshape(temp(:),[4 6]));
      v = [v;temp(:)]; %Add basal slices
      temp = m(:,3);
      temp=sum(reshape(temp(:),[4 6]));
      v = [v;temp(:)]; %Add mid slices
      temp = m(:,2);
      %shift it
      temp = [temp(23:24);temp(1:22)];
      temp=sum(reshape(temp(:),[6 4]));
      v = [v;temp(:)]; %Add apical slices
      v = [v;sum(m(:,1))];
    case 'max'
      tempm = zeros(size(m,1),4);
      if numslices == 3 || is16segments
        tempm(:,1) = NaN;  %apex
      else
        tempm(:,1) = squeeze(max(m(:,1:round(size(m,2)*0.15)),[],2));  %apex
      end
      tempm(:,2) = squeeze(max(m(:,firstslice:round(size(m,2)*ratios(1))),[],2));  %apical part
      tempm(:,3) = squeeze(max(m(:,round(size(m,2)*ratios(1))+1:round(size(m,2)*ratios(2))),[],2));  %mid-ventricular part
      tempm(:,4) = squeeze(max(m(:,round(size(m,2)*ratios(2))+1:end),[],2));  %basal part
      m = tempm;
      %Make AHA
      v = [];
      temp = m(:,4);
      temp=max(reshape(temp(:),[4 6]),[],1);
      v = [v;temp(:)]; %Add basal slices
      temp = m(:,3);
      temp=max(reshape(temp(:),[4 6]),[],1);
      v = [v;temp(:)]; %Add mid slices
      temp = m(:,2);
      %shift it
      temp = [temp(23:24);temp(1:22)];
      temp=max(reshape(temp(:),[6 4]),[],1);
      v = [v;temp(:)]; %Add apical slices
      v = [v;max(m(:,1))];  %apex
  end
end %End of v calculation

doplot = true;

if (nargin<2)||isempty(ax)
  doplot = false;
end

%Unpack vector to matrix
m = zeros(12,4);
temp = v(1:6);
temp=repmat(temp(:),[1 2])';
m(:,4) = temp(:); %Basal slice
temp = v(7:12);
temp=repmat(temp(:),[1 2])';
m(:,3) = temp(:); %mid slice
temp = v(13:16);
temp=repmat(temp(:),[1 3])';
m(:,2) = temp(:); %apical slice
temp = v(17);
temp=repmat(temp(:),[1 12])';
m(:,1) = temp(:); %apical slice

numslices = size(m,2)-1;
numsectors = 12;

if doplot || (nargout>1)

  [x,y] = ndgrid(...
    linspace(-numslices-1,numslices+1,2*n+1),...
    linspace(-numslices-1,numslices+1,2*n+1));
  rad = sqrt(x.*x+y.*y);

  %Createidx outer
  ang = angle(complex(y,x))+pi;
  ang = numsectors*ang/(2*pi);
  idxouter = 1+min(floor(ang),(numsectors-1))+(numsectors)*min(floor(rad),numslices);

  %Createidx inner
  ang = mod(angle(complex(y,x))+pi+pi/4,2*pi);
  ang = numsectors*ang/(2*pi);
  idxinner = 1+min(floor(ang),(numsectors-1))+(numsectors)*min(floor(rad),numslices);

  idx = idxouter;
  idx(rad<2) = idxinner(rad<2);

  im = m(idx);
  im(rad>(numslices+1)) = NaN;
  %im = rad;
end

if doplot
  scale = n/(numslices+1);

  %View data
  alpha = double(not(isnan(im)));
  im(isnan(im)) = 0;
  cmap = ax.Colormap;
  h = imagesc(im,'parent',ax);
  set(h,'alphadata',alpha,'AlphaDataMapping','scaled');
  axis(ax,'image','off');

  %Draw circles
  om = linspace(0,2*pi,100);
  xc = sin(om);
  yc = cos(om);
  hold(ax,'on');
  for loop=1:(numslices+1)
    h = plot(ax,n+1+scale*loop*xc,n+1+scale*loop*yc,'color',linecolor,'linewidth',2);
  end
  hold(ax,'off');

  %Draw lines
  hold(ax,'on');
  b = sqrt(0.75);
  a = 0.5;
  c = 1/sqrt(2);
  h = plot(ax,scale*[0 2],scale*[4 4],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[6 8],scale*[4 4],'color',linecolor,'linewidth',2);
  %h = plot(scale*[4 4],scale*[5 6],'w-'); set(h,'linewidth',2);
  %h = plot(scale*[4 4],scale*[2 3],'w-'); set(h,'linewidth',2);
  h = plot(ax,scale*[4-c 4-2*c],scale*[4-c 4-2*c],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[4+c 4+2*c],scale*[4+c 4+2*c],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[4-c 4-2*c],scale*[4+c 4+2*c],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[4+c 4+2*c],scale*[4-c 4-2*c],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[4-4*a 4-2*a],scale*[4-4*b 4-2*b],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[4-4*a 4-2*a],scale*[4+4*b 4+2*b],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[4+4*a 4+2*a],scale*[4+4*b 4+2*b],'color',linecolor,'linewidth',2);
  h = plot(ax,scale*[4+4*a 4+2*a],scale*[4-4*b 4-2*b],'color',linecolor,'linewidth',2);

  if includingvalues
    %define segment positions in im and corresponding indices in matrice
    segmentpositions = [...
      4,   4,   1, 1; ... % Apex
      2.5, 4,   1, 2; ... % Apical 1
      4,   2.5, 4, 2; ... % Apical 2
      5.5, 4,   7, 2; ... % Apical 3
      4,   5.5, 10,2; ... % Apical 4
      1.75,3,   1, 3; ... % Mid 1
      4,   1.5, 3, 3; ... % Mid 2
      6.25,3,   5, 3; ... % Mid 3
      6.25,5.25,7, 3; ... % Mid 4
      4,   6.5, 9, 3; ... % Mid 5
      1.75,5.25,11,3; ... % Mid 6
      1,   2,   1, 4; ... % Basal 1
      4,   0.5, 3, 4; ... % Basal 2
      7,   2,   5, 4; ... % Basal 3
      7,   6,   7, 4; ... % Basal 4
      4,   7.5, 9, 4; ... % Basal 5
      1,   6,   11,4  ... % Basal 6
      ];

    %loop over segments
    for loop = 1:length(segmentpositions)
      row = segmentpositions(loop,3);
      col = segmentpositions(loop,4);
      valuestr = getplotstring(m(row,col));
      if ~strcmp(valuestr,'NaN')
        x = scale*segmentpositions(loop,1);
        y = scale*segmentpositions(loop,2);
        textcolor = gettextcolor_helper(str2double(valuestr),cmap,v);
        %display segment value in bullseye plot
        text(x,y,valuestr,'parent',ax,'HorizontalAlignment','center','Color',textcolor,'FontWeight','bold');
      end
    end
  end

  hold(ax,'off');
end

%returns value if necessary.
if nargout>0
  varargout{1} = v;
end

if nargout>1
  varargout{2} = im;
end

%---------------------------------------------------
function textcolor = gettextcolor_helper(value,cmap,data)
%---------------------------------------------------
arguments
  value     %value of the segment to display in bullseye plot
  cmap      %bullseye's colormap
  data = [] %bullseyes's data
end
global DATA
gui = DATA.GUI.Bullseye;

if isempty(data)
  data = gui.outdata(:);
end

% get max/min values in bullseye plot
if ~isempty(gui)
  maxvalue = str2double(mygetedit(gui.handles.maxedit));
  minvalue = str2double(mygetedit(gui.handles.minedit));
else
  maxvalue = NaN;
  minvalue = NaN;
end
if isnan(maxvalue) || value > maxvalue
  maxvalue = ceil(max(data));
end
if isnan(minvalue) || value < minvalue
  minvalue = floor(min(data));
end
if maxvalue == 0 %avoid division by zero
  maxvalue = 1;
end

%normalize the segment value
segmentvalue = (value - minvalue) / (maxvalue - minvalue);

%get segment color and corresponding text color
segmentcolorRGB = getsegmentcolor(cmap,segmentvalue);
textcolor = gettextcolor(segmentcolorRGB);

%---------------------------------------------------
function segmentcolorRGB = getsegmentcolor(cmap,segmentvalue)
%---------------------------------------------------
%Return segment color based on segment's value and plot's colormap, used in bullseye plot.
n = size(cmap,1) - 1;
%scale the value to the index in the colormap
valueind = round(segmentvalue*n) + 1;
%retrieve RGB color from colormap
segmentcolorRGB = cmap(valueind,:);

%---------------------------------------------------
function textcolor = gettextcolor(backgroundcolorRGB)
%---------------------------------------------------
%Return text color based on brightness of background color, used in bullseye plot.

%calculate brightness of background color
brightness = colorfunctions.getbrightness(backgroundcolorRGB);

%return text color based on brightness
if brightness < 0.5
  textcolor = [1 1 1]; %white for dark background
else
  textcolor = [0 0 0]; %black for light background
end

%---------------------------------------------------
function plotstring = getplotstring(value)
%---------------------------------------------------
% get plot string with correct formatting depending on the value
if value >= 1000
  formatstring = '%0.4g';
else
  formatstring = '%0.3g';
end
plotstring = sprintf(formatstring,value);

%---------------------------------------------------
function [varargout] = bullseye2(m,ax,n,flipx,vc,no)
%---------------------------------------------------
%Generate bullesye data. Optional output argument is im
%m:     values in the sectors
%ax:    figure axes
%n:     number of pixels in the resulting image
%flipx: flip in x-direction (true or false)
%vc:    volume consistent (true or false)

global SET NO;

%Set constants
if nargin < 6
  no = NO;
  if nargin < 5
    vc = false;
    if nargin < 4
      flipx = false;
      if nargin < 3
        n = 200;
      end
    end
  end
end

%Fix orientation of m
m = m';
m = flipud(m);

%Make sure enough sectors. Pad array before upsampling to avoid sharp edge
minnbrsectors = 40;
if size(m,2) < minnbrsectors
  m = imresize(repmat(m,1,3),minnbrsectors,'bilinear');
  m = m(1:minnbrsectors:end,size(m,2)/3+1:2*size(m,2)/3);
end

%Determine size
sectors = size(m,2);
slices = size(m,1);

%Sampled positions in polar coordinates
omega = 2*pi*(0:(sectors-1))/(sectors-1)-pi;
%was omega = 2*pi*(0:(sectors-1))/(sectors)-pi;
omega = repmat(omega,[slices 1]);

if vc
  %calculate the myocardial volume for each sector
  ind = false(1,SET(no).ZSize);
  ind(SET(no).StartSlice:SET(no).EndSlice) = true;
  myocardvolume = calcfunctions('calcmyocardvolume',sectors,no);
  myocardvolume = myocardvolume(:,ind);
  myocardvolume = fliplr(myocardvolume);  %flip to have apex first
  myocardvolumeslice = sum(myocardvolume,'omitnan');
  myocardvolumesliceper = myocardvolumeslice./(sum(myocardvolumeslice(:)));
  %calculate the radius in the bullseye for each slice
  rtotal = slices;
  r = ones(1,slices);
  Atotal = rtotal^2*pi;
  A = Atotal*myocardvolumesliceper;
  for rloop = 1:slices
    r(rloop) = sqrt(sum(A(1:rloop))/pi)-sum(r(1:rloop-1));
  end
  rcumsum = cumsum(r./rtotal);
  rcumsum(end) = 1;  %correction for roundoff errors
  rcumsum = [0 rcumsum];
  rmid = zeros(1,length(rcumsum)-1);
  for loop = 1:length(rcumsum)-1
    rmid(loop) = mean(rcumsum(loop:loop+1));
  end
  rmid(end) = 1;
  r = repmat(rmid',[1 sectors]);

  %Convert to cartesian
  x = r.*sin(omega);
  y = r.*cos(omega);

  %Interpolation positions
  [xi,yi] = ndgrid(...
    linspace(-1,1,2*n+1),...
    linspace(-1,1,2*n+1));
else
  r = 1:slices;
  r = repmat(r',[1 sectors]);

  %Convert to cartesian
  x = r.*sin(omega);
  y = r.*cos(omega);

  %Interpolation positions
  [xi,yi] = ndgrid(...
    linspace(-slices,slices,2*n+1),...
    linspace(-slices,slices,2*n+1));
end

%Interpolate
%im = griddata(x,y,m,xi,yi,'linear');
%Avoid sharp edge caused by duplicate values very near angle 0
im = griddata([zeros(slices,1) x(:,2:end-1)], ...
  y(:,1:end-1),[(m(:,1)+m(:,end))/2 m(:,2:end-1)],xi,yi,'linear');

alpha = double(not(isnan(im)));

if flipx
  im = fliplr(im);
end

if nargout==0
  h = imagesc(im,'parent',ax);
  set(h,'alphadata',alpha,'AlphaDataMapping','scaled');
  axis(ax,'image','off');
else
  varargout{1} = im;
end

%-----------------------------------
function im = bullseye(m,ax,n,vc,no,tf)
%-----------------------------------
%Calculate bullseye from matrix m. Ax is optional axis where the output
%image im should be displayed.
%m:   values in the sectors
%ax:  figure axes
%n:   number of pixels in the resulting image
%vc:  volume consistent (true or false)
%im:  resulting bullseye image

global SET NO;

%Set constants
if nargin < 6
  tf = SET(NO).CurrentTimeFrame;
  if nargin < 5
    no = NO;
    if nargin < 4
      vc = false;
      if nargin < 3
        n = 200;
      end
    end
  end
end

m = fliplr(m); %Basal slice is first column

%Determine size
sectors = size(m,1);
slices = size(m,2);

%Add extra column (slice) for apex
m = [nan(sectors,1) m];

%calculate the myocardial volume for each sector
ind = false(1,SET(no).ZSize);
ind(SET(no).StartSlice:SET(no).EndSlice) = true;
if length(tf) > 1
  for tfloop = 1:length(tf)
    myocardvolumetf(:,:,tfloop) = calcfunctions('calcmyocardvolume',sectors,no,tf(tfloop));
  end
  myocardvolume = mynanmean(myocardvolumetf,3);
else
  myocardvolume = calcfunctions('calcmyocardvolume',sectors,no,tf);
end
myocardvolume = myocardvolume(:,ind);
myocardvolume = fliplr(myocardvolume);  %flip to have apex first
myocardvolumeslice = nansum(myocardvolume,1);
ismyocardium = repmat([false logical(myocardvolumeslice)],[sectors 1]); %add extra slice for apex
m(~ismyocardium) = NaN;
totalmyocardvolume = nansum(myocardvolumeslice(:));
myocardvolumeslice(myocardvolumeslice==0) = totalmyocardvolume/length(myocardvolumeslice);
myocardvolumesliceper = myocardvolumeslice./(nansum(myocardvolumeslice(:))+1^2*pi);
myocardvolumesliceper = myocardvolumesliceper/sum(myocardvolumesliceper);

if vc
  %calculate the radius in the bullseye for each slice
  rtotal = slices+1;
  r = ones(1,slices+1);
  r(1) = 1;
  Atotal = rtotal^2*pi;
  A = [r(1)^2*pi Atotal*myocardvolumesliceper];
  for rloop = 1:slices
    r(rloop+1) = sqrt(sum(A(1:rloop+1))/pi)-sum(r(1:rloop));
  end
  [x,y] = ndgrid(...
    linspace(-1,1,2*n+1),...
    linspace(-1,1,2*n+1));
  ang = angle(complex(y,x))+pi;
  ang = sectors*ang/(2*pi);
  rad = sqrt(x.*x+y.*y);
  rcumsum = cumsum(r./rtotal);
  rcumsum(end) = 1;  %correction for roundoff errors
  rad(abs(rad)<rcumsum(1)) = NaN;
  rad(abs(rad)>rcumsum(end)) = NaN;
  for radloop = 2:slices+1
    rad(abs(rad)<rcumsum(radloop)) = radloop-1;
  end
  %Create idx
  idx = 1+min(floor(ang),(sectors-1))+(sectors)*min(floor(rad),slices);
  im = m(idx);
  im(isnan(rad)) = NaN;
  scale = n/(slices+1);

else

  scale = n/(slices+1);
  [x,y] = ndgrid(...
    linspace(-slices-1,slices+1,2*n+1),...
    linspace(-slices-1,slices+1,2*n+1));
  ang = angle(complex(y,x))+pi;
  ang = sectors*ang/(2*pi);
  rad = sqrt(x.*x+y.*y);
  %Create idx
  idx = 1+min(floor(ang),(sectors-1))+(sectors)*min(floor(rad),slices);
  im = m(idx);
  im(rad>(slices+1)) = NaN;
end

doplot = true;
if (nargin<2)||isempty(ax)
  doplot = false;
end

if doplot
  global DATA %#ok<TLEV> call DATA only if doplot because we need to know which freground color to use
  %View data
  alpha = double(not(isnan(im)));
  im(isnan(im)) = 0;
  cmap = ax.Colormap;
  h = imagesc(im,'parent',ax);
  set(h,'alphadata',alpha,'AlphaDataMapping','scaled');
  axis(ax,'image','off');

  %Draw circles
  om = linspace(0,2*pi,100);
  xc = sin(om);
  yc = cos(om);
  hold(ax,'on');
  linecolor = DATA.GUISettings.ForegroundColor;
  for loop=1:(slices+1)
    if vc
      h = plot(ax,n+1+n*rcumsum(loop)*xc,n+1+n*rcumsum(loop)*yc,'color',linecolor);
    else
      h = plot(ax,n+1+scale*loop*xc,n+1+scale*loop*yc,'color',linecolor);
    end
    set(h,'linewidth',2);
  end
  hold(ax,'off');

  %Label sectors
  hold(ax,'on');
  omega = 2*pi*(0:(sectors-1))/sectors+pi;
  omega = omega+0.5*(2*pi/sectors); %omega(2)-omega(1)

  if sectors < 21
    for sliceloop=1:(slices)
      if vc
        slicescale = n*mean([rcumsum(sliceloop) rcumsum(sliceloop+1)]);
      else
        slicescale = scale*mean([sliceloop sliceloop+1]);
      end
      for loop=1:sectors
        valuestr = getplotstring(m(loop,1+sliceloop));
        x = n+1+slicescale*1*cos(omega(loop));
        y = n+1+slicescale*1*sin(omega(loop));
        textcolor = gettextcolor_helper(str2double(valuestr),cmap);
        text(x,y,valuestr, ...
          'HorizontalAlignment','center', ...
          'parent',ax,'color',textcolor,'FontWeight','bold');
      end
    end
  end
  hold(ax,'off');
end

%-----------------------------
function value = getdata(type)
%-----------------------------
%Helper function to extract data from the module. Used for instance
%from reportsheet generator.

global DATA
gui = DATA.GUI.Bullseye;
switch type
  case 'getahadata'
    value = gui.ahaoutdata;
  case 'getmainframe'
    im = frame2im(mygetframe(gui.fig));
    value = im;
  case 'getbullseyeframe'
    im = frame2im(mygetframe(gui.handles.bullseyeaxes));
    value = im;
  case 'getimageframe'
    im = frame2im(mygetframe(gui.handles.imageaxes));
    value = im;
  case 'getcolorbar'
    im = frame2im(mygetframe(gui.handles.colorbar));
    value = im;
end

%----------------------------
function setdata(type,datain)
%----------------------------
%interface function to gui to set data from code, for instance
%reportsheet generator.

global DATA
gui = DATA.GUI.Bullseye;

switch type
  case 'setlistbox'
    set(gui.handles.bullseyelistbox,'value',datain);
  case 'setseparate'
    set(gui.handles.separatewindowcheckbox,'value',datain);
  case 'setnormal'
    set(gui.handles.sectorradiobutton,'value',1);
  case 'setsmooth'
    set(gui.handles.smoothradiobutton,'value',1);
  case 'setaha'
    set(gui.handles.aharadiobutton,'value',1);
  case 'setsectors'
    set(gui.handles.sectorslistbox,'value',datain);
    temp = get(gui.handles.sectorslistbox,'String');
    gui.numsectors = str2num(temp{mygetlistbox(gui.handles.sectorslistbox)});     %#ok<ST2NM>
end