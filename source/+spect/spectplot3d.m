function spectplot3d(varargin)
%----------------------------
%Plot a 3D view of the left ventricle in the spect image
%written by Helen Soneson 2008-03-03

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:});

%---------------------
function init_Callback %#ok<DEFNU>
%---------------------
%open the 3D view

global SET NO DATA

if isempty(SET(NO).EndoX) || isempty(SET(NO).EpiX)
  myfailed('Need LV segmentation');
  return;
end

if SET(NO).TSize > 1
  myfailed('Not yet implemented for time resolved images.');
  return;
end

tempnos = NO;
imissingle = classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

scar = [];
%finding slices with segmentation
startslices = zeros(1,SET(NO).TSize);
endslices = zeros(1,SET(NO).TSize);
for tf = 1:SET(NO).TSize
  startslices(tf) = find(~isnan(SET(NO).EndoX(1,tf,:)),1,'first');  %the most basal slice
  endslices(tf) = find(~isnan(SET(NO).EndoX(1,tf,:)),1,'last');  %the most apical slice
end
nbrofslices = endslices-startslices+1;
nbroftimeframes = SET(NO).TSize;
if isfield(SET(NO),'Scar') && isfield(SET(NO).Scar,'Auto') && ~isempty(SET(NO).Scar.Auto) && sum(SET(NO).Scar.Auto(:)) > 0
  isscar = 1;
else 
  isscar = 0;
end
meanresolution = mean([SET(NO).ResolutionX SET(NO).ResolutionY]);
% resampledIM = spect.specttools('upsamplevolume',samplefactor,samplefactor,SET(NO).IM);  %resample the image stack to a resolution of 1x1 mm
nbrofpoints = 100;  %number of points around the myocardium in each slice
nbrofminpoints = 10;  %number of points in the search area for the min index
midmuralcount = NaN(nbrofpoints,SET(NO).TSize,max(nbrofslices));
midmuralx = midmuralcount;
midmuraly = midmuralcount;
myocardium = midmuralcount;
if isscar
  scar = midmuralcount;
end
%calculate the mean radial count
for tf = 1:SET(NO).TSize
  for slice = startslices(tf):endslices(tf)  %loop over all slices with segmentation
    [endox,endoy,epix,epiy] = spect.specttools('centerlinemethod',slice,NO,tf,nbrofpoints,nbrofminpoints);
    dist = sqrt((endox-epix).^2+(endoy-epiy).^2);
    for point = 1:nbrofpoints
      if dist(point) > 1/meanresolution
        [cx,cy,count] = improfile(SET(NO).IM(:,:,1,slice),[endoy(point); epiy(point)],[endox(point); epix(point)]);
        midmuralcount(point,tf,slice-startslices(tf)+1) = mean(count);
        midmuralx(point,tf,slice-startslices(tf)+1) = mean(cx);
        midmuraly(point,tf,slice-startslices(tf)+1) = mean(cy);
        myocardium(point,tf,slice-startslices(tf)+1) = 1;
        if isscar
          scartemp = improfile(SET(NO).Scar.Result(:,:,slice),[endoy(point); epiy(point)],[endox(point); epix(point)]);
          scar(point,tf,slice-startslices(tf)+1) = round(mean(scartemp)+0.3);
        end
      else
        midmuralx(point,tf,slice-startslices(tf)+1) = mean([endoy(point) epiy(point)]);
        midmuraly(point,tf,slice-startslices(tf)+1) = mean([endox(point) epix(point)]);
      end
    end
  end
end
midmuralz(:,1:nbroftimeframes,:) = repmat(nbrofslices:-1:1,[nbrofpoints nbroftimeframes 1]);
maxcount = max(midmuralcount(:));
midmuralcount = 100*midmuralcount/maxcount;

%fill in holes in the myocardium mask
for tf = 1:nbroftimeframes
  myocardium(:,tf,:) = removeholes(squeeze(myocardium(:,tf,:)));
  if isscar
    scar(:,tf,:) = removeholes(squeeze(scar(:,tf,:)));
  end
end

[outflowtractrow outflowtractcol] = find(isnan(squeeze(myocardium(:,1,:))));
nbrofbasal = max(outflowtractcol);
if isempty(nbrofbasal)
  nbrofbasal = 0;
end
%create a smooth outflow tract
if isscar
  [myocardium,midmuralcount,scar] = smoothoutflow(myocardium,nbroftimeframes,nbrofbasal,midmuralcount,scar);
else
  [myocardium,midmuralcount] = smoothoutflow(myocardium,nbroftimeframes,nbrofbasal,midmuralcount);
end

%mask with the myocardium mask
midmuralx = midmuralx.*myocardium;
midmuraly = midmuraly.*myocardium;
midmuralz = midmuralz.*myocardium;
midmuralcount = midmuralcount.*myocardium;

%rotate the LV so the image "starts" in the outflow tract in order to have
%a hole there
if isempty(outflowtractrow)
  shift = round(nbrofpoints/2);
else
  shift = round(mean(outflowtractrow));  %middle of the outflow tract
end
midmuralx = circshift(midmuralx,[-shift 0 0]);
midmuraly = circshift(midmuraly,[-shift 0 0]);
midmuralz = circshift(midmuralz,[-shift 0 0]);
midmuralcount = circshift(midmuralcount,[-shift 0 0]);
if isscar
  scar = circshift(scar,[-shift 0 0]);
end

%connect apex
for tf = 1:nbroftimeframes
  apicalcenterx = mean(midmuralx(:,tf,end));
  apicalcentery = mean(midmuraly(:,tf,end));
  apicalcount = mean(midmuralcount(:,tf,end));
%   apicalradius = sqrt((apicalcenterx-midmuralx(:,tf,end-2:end)).^2+(apicalcentery-midmuraly(:,tf,end-2:end)).^2);
%   apicalradiusdeformed = interp2([3 2 1],apicalradius,[3 2.5 2 1.5 1 0.5])*exp(linspace(-0.05,-1.0,6));
  midmuralx(:,tf,end+1) = apicalcenterx;
  midmuraly(:,tf,end+1) = apicalcentery;
  midmuralz(:,tf,end+1) = 0;
  midmuralcount(:,tf,end+1) = apicalcount;
  if isscar
    scar(:,tf,end+1) = median(scar(:,tf,end));
  end
end

%connect the surface below the outflow tract
midmuralx(nbrofpoints+1,:,nbrofbasal+1:end) = midmuralx(1,:,nbrofbasal+1:end);
midmuraly(nbrofpoints+1,:,nbrofbasal+1:end) = midmuraly(1,:,nbrofbasal+1:end);
midmuralz(nbrofpoints+1,:,nbrofbasal+1:end) = midmuralz(1,:,nbrofbasal+1:end);
midmuralcount(nbrofpoints+1,:,nbrofbasal+1:end) = midmuralcount(1,:,nbrofbasal+1:end);
if isscar
  scar(nbrofpoints+1,:,nbrofbasal+1:end) = scar(1,:,nbrofbasal+1:end);
end
midmuralx(nbrofpoints+1,:,1:nbrofbasal) = NaN;
midmuraly(nbrofpoints+1,:,1:nbrofbasal) = NaN;
midmuralz(nbrofpoints+1,:,1:nbrofbasal) = NaN;
midmuralcount(nbrofpoints+1,:,1:nbrofbasal) = NaN;
if isscar
  scar(nbrofpoints+1,:,1:nbrofbasal) = NaN;
end

%fix correct proportions
midmuralx = midmuralx*SET(NO).ResolutionX;
midmuraly = midmuraly*SET(NO).ResolutionY;
midmuralz = midmuralz*(SET(NO).SliceThickness+SET(NO).SliceGap);

%move the surface to be in the middle of the image
minx = min(midmuralx(:));
miny = min(midmuraly(:));
minz = min(midmuralz(:));
midmuralx = midmuralx-minx;
midmuraly = midmuraly-miny;
midmuralz = midmuralz-minz;

if isscar
  %invers the scar image
  scar(scar == 0) = 2;
  scar(scar == 1) = 0;
  scar(scar == 2) = 1;
end

%erase the outflow tract
for tf = 1:nbroftimeframes
  for slice = 1:nbrofbasal
    startmyo = find(~isnan(midmuralx(:,tf,slice)),1,'first');
    endmyo = find(~isnan(midmuralx(:,tf,slice)),1,'last');
    lengthmyo = endmyo-startmyo+1;
    xx = linspace(0,1,lengthmyo);
    xi = linspace(0,1,nbrofpoints+1);
    tempmyo = interp1(xx,midmuralx(startmyo:endmyo,tf,slice),xi);
    midmuralx(:,tf,slice) = tempmyo;
    tempmyo = interp1(xx,midmuraly(startmyo:endmyo,tf,slice),xi);
    midmuraly(:,tf,slice) = tempmyo;
    tempmyo = interp1(xx,midmuralz(startmyo:endmyo,tf,slice),xi);
    midmuralz(:,tf,slice) = tempmyo;
    tempmyo = interp1(xx,midmuralcount(startmyo:endmyo,tf,slice),xi);
    midmuralcount(:,tf,slice) = tempmyo;
    if isscar
      tempmyo = interp1(xx,scar(startmyo:endmyo,tf,slice),xi);
      scar(:,tf,slice) = tempmyo;
    end
  end
end

%plot the 3d image
if isopengui(['+spect' filesep 'spectplot3d.fig'])
  gui = DATA.GUI.SpectPlot3d;
  figure(gui.fig);
else
	DATA.GUI.SpectPlot3d = mygui(['+spect' filesep 'spectplot3d.fig']);
  gui = DATA.GUI.SpectPlot3d;
  myadjust(gui.fig,DATA.GUI.Segment);
end
gui.NO = NO;
gui.handles.tf = 1;
gui.handles.midmuralcount = midmuralcount;
gui.handles.midmuralx = midmuralx;
gui.handles.midmuraly = midmuraly;
gui.handles.midmuralz = midmuralz;
gui.handles.scar = scar;
if isscar
  gui.handles.isscar = 1;
else
  gui.handles.isscar = 0;
end
gui.handles.axes1image = surf(squeeze(gui.handles.midmuraly(:,gui.handles.tf,:)),squeeze(gui.handles.midmuralx(:,gui.handles.tf,:)), ...
  squeeze(gui.handles.midmuralz(:,gui.handles.tf,:)),squeeze(gui.handles.midmuralcount(:,gui.handles.tf,:)), ...
  'FaceColor','interp','parent',gui.handles.axes1);
set(gui.handles.axes1,'clim',[0 100]);
colormap('spect');
colorbar('peer',gui.handles.axes1);
shading(gui.handles.axes1,'interp');
axis(gui.handles.axes1,'image','off');
h = rotate3d(gui.handles.axes1);
set(h,'Enable','on');

if isscar
  %initial not show scar
  set(gui.handles.radiobuttonshow,'value',0);
  set(gui.handles.radiobuttonhide,'value',1);
else
  set(gui.handles.radiobuttonshow,'enable','off');
  set(gui.handles.radiobuttonhide,'enable','off');
end

if SET(NO).TSize == 1
  %not possibly to play movie
  set(gui.handles.pushbuttonPrev,'enable','off');
  set(gui.handles.pushbuttonNext,'enable','off');
  set(gui.handles.pushbuttonPlay,'enable','off');
  set(gui.handles.pushbuttonStop,'enable','off');
end


%--------------------------------
function mask = removeholes(mask)
%--------------------------------
%remove holes
%
%written by Helen Soneson 2010-05-05

smallarea = 10;  %the size of the small area in pixels
n = 4;
sz = size(mask);
%define the neighborhood
tempmask = mask;
tempmask(isnan(tempmask)) = 0;
[index,numindex] = bwlabeln(tempmask,n);  %make clusters of the holes

for labelloop = 1:numindex  %loop over all holes
  holearea= length(find(index==labelloop));
  if holearea <= smallarea  %the hole is smaller than the threshold
    [xx,yy] = find(index==labelloop);
    for p = 1:length(xx)
      medianvalue = median([mask(max([1 xx(p)-1]),yy(p)) mask(min([sz(1) xx(p)+1]),yy(p)) ...
        mask(xx(p),max([1 yy(p)-1])) mask(xx(p),min([sz(2) yy(p)+1]))]);
      mask(xx(p),yy(p)) = medianvalue;  %#ok<FNDSB> %fill in the small holes
    end
  end
end

tempmask = mask;
tempmask(isnan(tempmask)) = 1;
[index,numindex] = bwlabeln(tempmask,n);  %make clusters of the holes

for labelloop = 1:numindex  %loop over all holes
  holearea= length(find(index==labelloop));
  if holearea <= smallarea  %the hole is smaller than the threshold
    [xx,yy] = find(index==labelloop);
    for p = 1:length(xx)
      medianvalue = median([mask(max([1 xx(p)-1]),yy(p)) mask(min([sz(1) xx(p)+1]),yy(p)) ...
        mask(xx(p),max([1 yy(p)-1])) mask(xx(p),min([sz(2) yy(p)+1]))]);
      mask(xx(p),yy(p)) = medianvalue;  %#ok<FNDSB> %fill in the small holes
    end
  end
end

tempmask = mask;
tempmask(tempmask == 0) = 2;
tempmask(tempmask == 1) = 0;
tempmask(tempmask == 2) = 1;
tempmask(isnan(tempmask)) = 1;
[index,numindex] = bwlabeln(tempmask,n);  %make clusters of the holes

for labelloop = 1:numindex  %loop over all holes
  holearea= length(find(index==labelloop));
  if holearea <= smallarea  %the hole is smaller than the threshold
    [xx,yy] = find(index==labelloop);
    for p = 1:length(xx)
      medianvalue = median([mask(max([1 xx(p)-1]),yy(p)) mask(min([sz(1) xx(p)+1]),yy(p)) ...
        mask(xx(p),max([1 yy(p)-1])) mask(xx(p),min([sz(2) yy(p)+1]))]);
      mask(xx(p),yy(p)) = medianvalue;  %#ok<FNDSB> %fill in the small holes
    end
  end
end

tempmask = mask;
tempmask(tempmask == 0) = 2;
tempmask(tempmask == 1) = 0;
tempmask(tempmask == 2) = 1;
tempmask(isnan(tempmask)) = 0;
[index,numindex] = bwlabeln(tempmask,n);  %make clusters of the holes

for labelloop = 1:numindex  %loop over all holes
  holearea= length(find(index==labelloop));
  if holearea <= smallarea  %the hole is smaller than the threshold
    [xx,yy] = find(index==labelloop);
    for p = 1:length(xx)
      medianvalue = median([mask(max([1 xx(p)-1]),yy(p)) mask(min([sz(1) xx(p)+1]),yy(p)) ...
        mask(xx(p),max([1 yy(p)-1])) mask(xx(p),min([sz(2) yy(p)+1]))]);
      mask(xx(p),yy(p)) = medianvalue;  %#ok<FNDSB> %fill in the small holes
    end
  end
end


%--------------------------------------------------------------------------
function [mask,midmuralcount,scar] = smoothoutflow(mask,nbroftimeframes,nbrofbasal,midmuralcount,scar)
%--------------------------------------------------------------------------
%smooth the outflow tract
%
%written by Helen Soneson 2010-05-05

if nargin == 5
  isscar = 1;
else 
  isscar = 0;
  scar = [];
end

outflowtract = zeros(nbrofbasal*2,nbroftimeframes);
for tf = 1:nbroftimeframes
  for slice = 1:nbrofbasal
    outflowtracttemp = find(isnan(squeeze(mask(:,tf,slice))));
    outflowtract(2*slice-1,tf) = min(outflowtracttemp);
    outflowtract(2*slice,tf) = max(outflowtracttemp);
  end
  middleoutflowtract = mean(outflowtract);
  move = outflowtract(2:2:end)-middleoutflowtract;
  outflowtract(2:2:end) = middleoutflowtract-move;
  %fit a 3-dimensional curve to the outflow tract border
  x = floor(1:0.5:nbrofbasal+0.5);
  y = outflowtract(:,tf)';
  p = polyfit(x,y,2);
  borders = round(p(1)*x.^2+p(2)*x+p(3)); 
  outflowtractstart = borders(1:2:end);
  move = middleoutflowtract-borders(1:2:end);
  outflowtractend = round(borders(1:2:end)+2*move);
  for slice = 1:nbrofbasal
    mask(outflowtractstart(slice):outflowtractend(slice),tf,slice) = NaN;
    midmuralcount(outflowtractstart(slice):outflowtractend(slice),tf,slice) = NaN;
    if isscar
      scar(outflowtractstart(slice):outflowtractend(slice),tf,slice) = NaN;
    end
    fillout = find(isnan(mask(1:outflowtractstart(slice),tf,slice)));
    if ~isempty(fillout)
      mask(fillout,tf,slice) = 1;
      midmuralcount(fillout,tf,slice) = midmuralcount(min(fillout)-1,tf,slice);
      if isscar
        fillout = find(isnan(scar(1:outflowtractstart(slice),tf,slice)));
        scar(fillout,tf,slice) = scar(min(fillout)-1,tf,slice);
      end
    end
    fillout = find(isnan(mask(outflowtractend(slice):end,tf,slice)));    
    if ~isempty(fillout)
      fillout = fillout+outflowtractend(slice)-1;
      mask(fillout,tf,slice) = 1;
      midmuralcount(fillout,tf,slice) = midmuralcount(max(fillout)+1,tf,slice);
      if isscar
        fillout = find(isnan(scar(outflowtractend(slice):end,tf,slice)));
        fillout = fillout+outflowtractend(slice)-1;
        scar(fillout,tf,slice) = scar(max(fillout)+1,tf,slice);
      end
    end
  end
end


%--------------------
function prev_Callback %#ok<DEFNU>
%--------------------
%one time frame backwards

global DATA SET

gui = DATA.GUI.SpectPlot3d;

if SET(gui.NO).TSize == 1
  myfailed('There is only one time frame',DATA.GUI.SpectPlot3d);
  return;
end

gui.handles.timeframe = gui.handles.timeframe-1;
SET(gui.NO).CurrentTimeFrame = gui.handles.timeframe;
if gui.handles.timeframe < 1
  gui.handles.timeframe = SET(gui.NO).TSize;
end;
spect.spectplot3d('trackupdate');


%--------------------
function next_Callback %#ok<DEFNU>
%--------------------
%one time frame forwards

global DATA SET

gui = DATA.GUI.SpectPlot3d;

if SET(gui.NO).TSize == 1
  myfailed('There is only one time frame',DATA.GUI.SpectPlot3d);
  return;
end

gui.handles.timeframe = gui.handles.timeframe+1;
SET(gui.NO).CurrentTimeFrame = gui.handles.timeframe;
if gui.handles.timeframe > SET(gui.NO).TSize
  gui.handles.timeframe = 1;
end;
spect.spectplot3d('trackupdate');


%--------------------
function play_Callback %#ok<DEFNU>
%--------------------
%play a movie

global DATA SET

gui = DATA.GUI.SpectPlot3d;

if SET(gui.NO).TSize == 1
  myfailed('There is only one time frame',DATA.GUI.SpectPlot3d);
  return;
end

DATA.StartFrame = gui.handles.timeframe;
DATA.StartTime = now;
try
  while get(gui.handles.pushbuttonPlay,'value')
    gui.handles.timeframe = segment('getframenumber');
    SET(gui.NO).CurrentTimeFrame = gui.handles.timeframe;
    spect.spectplot3d('trackupdate');
    pause(0.5*SET(gui.NO).BeatTime/SET(gui.NO).TSize);
  end;
catch %#ok<CTCH>
end


%-----------------------
function stop_Callback
%-----------------------
set(gui.handles.pushbuttonPlay,'value',0);


%-------------------
function trackupdate
%-------------------
%update the tracking

global DATA

gui = DATA.GUI.SpectPlot3d;

tempnos = gui.NO;
imissingle = classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

set(gui.handles.axes1image,'xdata',[]);
set(gui.handles.axes1image,'ydata',[]);
set(gui.handles.axes1image,'zdata',[]);
set(gui.handles.axes1image,'cdata',[]);
set(gui.handles.axes1image,'xdata',gui.handles.midmuralx{gui.handles.timeframe});
set(gui.handles.axes1image,'ydata',gui.handles.midmuraly{gui.handles.timeframe});
set(gui.handles.axes1image,'zdata',gui.handles.midmuralz{gui.handles.timeframe});
set(gui.handles.axes1image,'cdata',gui.handles.midmuralcolor{gui.handles.timeframe});
shading interp;
axis image off;
try
  if DATA.Record
    drawnow;
    DATA.MovieFrame = mygetframe(get(get(gui.handles.image,'parent'),'parent'));
    export('exportmovierecorder_Callback','newframe');
  end
catch %#ok<CTCH>
end


%--------------------------
function close_Callback %#ok<DEFNU>
%--------------------------
%close the gui

global DATA

try
  DATA.GUI.SpectPlot3d=close(DATA.GUI.SpectPlot3d);
catch %#ok<CTCH>
  close(gcf)
end


%-------------------------
function showscar_Callback %#ok<DEFNU>
%-------------------------
%show the scar segmentation

global DATA

gui = DATA.GUI.SpectPlot3d;
tf = gui.handles.tf;

if gui.handles.isscar
  set(gui.handles.axes1image,'cdata',squeeze(gui.handles.midmuralcount(:,tf,:).*gui.handles.scar(:,tf,:)));
  set(gui.handles.radiobuttonshow,'value',1);
  set(gui.handles.radiobuttonhide,'value',0);
else
  %no scar
  myfailed('There is no scar to show',DATA.GUI.SpectPlot3d);
end


%-------------------------
function hidescar_Callback %#ok<DEFNU>
%-------------------------
%hide the scar segmentation

global DATA

gui = DATA.GUI.SpectPlot3d;
tf = gui.handles.tf;

if gui.handles.isscar
  set(gui.handles.axes1image,'cdata',squeeze(gui.handles.midmuralcount(:,tf,:)));
  set(gui.handles.radiobuttonshow,'value',0);
  set(gui.handles.radiobuttonhide,'value',1);
else
  %no scar
  myfailed('There is no scar to show',DATA.GUI.SpectPlot3d);
end




% % text('String','apex','Position', ...
% %   [mean([min(min(epix))/SET(NO).ResolutionX max(max(epix))/SET(NO).ResolutionX]) ...
% %      mean([min(min(epix))/SET(NO).ResolutionX max(max(epix))/SET(NO).ResolutionX]) min(min(z))-0.05*abs(min(min(z)))], ....
% %   'FontSize',13)
% % text('String','septum','Position', ...
% %   [min(min(epix))/SET(NO).ResolutionX mean([min(min(epix))/SET(NO).ResolutionX max(max(epix))/SET(NO).ResolutionX]) max(max(z))], ....
% %   'FontSize',13)
% % text('String','lateral','Position', ...
% %   [max(max(epix))/SET(NO).ResolutionX mean([min(min(epix))/SET(NO).ResolutionX max(max(epix))/SET(NO).ResolutionX]) max(max(z))], ....
% %   'FontSize',13)
% % text('String','inferior','Position', ...
% %   [mean([min(min(epix))/SET(NO).ResolutionX max(max(epix))/SET(NO).ResolutionX]) min(min(epix))/SET(NO).ResolutionX max(max(z))], ....
% %   'FontSize',13)
% % text('String','anterior','Position', ...
% %   [mean([min(min(epix))/SET(NO).ResolutionX max(max(epix))/SET(NO).ResolutionX]) max(max(epix))/SET(NO).ResolutionX max(max(z))], ....
% %   'FontSize',13)
% %plot3(cornerx,cornery,cornerz,'k');
% %mesh(cornerx,cornery,cornerz);
% % hold off;