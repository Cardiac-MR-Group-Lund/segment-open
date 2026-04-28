function varargout = flow(varargin)
%2D flow tracking and analysis in MR images

%Einar Heiberg

%#ok<*GVMIS>

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%--------------------------------
function floweddycurrent_Callback 
%--------------------------------
%This function have been moved to floweddycurrent.m. Keeping stub for
%backwards compability for scripts and testing software.

floweddycurrent('init');

%-------------------------------------
function flowsegmentroi_Callback(name) 
%-------------------------------------
%Do segmentation of a flow roi

global SET NO DATA  

if isempty(SET(NO).Flow)
  myfailed('No flow data.',DATA.GUI.Segment);
  return;
end

if isempty(SET(NO).Flow.PhaseNo)
  myfailed('No through plane velocity data found.',DATA.GUI.Segment);
  return;
end

tempnos=[SET(NO).Flow.PhaseNo,SET(NO).Flow.MagnitudeNo];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

tools('enableundo');

%Find index to phase image
no = SET(NO).Flow.PhaseNo;

if isequal(no,NO)
  NO = SET(no).Flow.MagnitudeNo;
  viewfunctions('switchimagestack',NO);
end

p = [];

%Find what slices
if SET(NO).ZSize==1
  nslices = 1;
  ind = 1;
  p.zofs = 0;
else
  nslices = SET(NO).EndSlice-SET(NO).StartSlice+1;
  ind = false(1,SET(NO).ZSize);
  ind(SET(NO).StartSlice:SET(NO).EndSlice)=true;
  ind = find(ind);
  if isempty(ind)
    myfailed('No slices selected.',DATA.GUI.Segment);
    return;
  end
  p.zofs = ind(1)-1;
end

%Find peak velocity timeframe
d = 3;
x = round(SET(NO).CenterX);
y = round(SET(NO).CenterY);
meanvel = zeros(1,SET(NO).TSize);
for tloop=1:SET(NO).TSize
  temp = SET(no).IM(...
    (x-d):(x+d),(y-d):(y+d),tloop,ind);
  meanvel(tloop) = abs(mean(temp(:))-0.5);
end

[temp,maxtf] = max(meanvel);
magim = SET(NO).IM(:,:,maxtf,ind);
phaseim = abs(SET(no).IM(:,:,maxtf,ind)-0.5);

if temp>0.1
  speedim = magim.*phaseim*(1/temp)-0.2;
else
  speedim = magim-0.3;
end

%Find magnitude edges
[e0,e1,e2,e3] = newfindedges(...
    speedim,...
    [-1;2;-1],...
    [1;2;1]/4,...
    [1;2;1]/4,...
    0,0,0);

%Set forces
p.M_f = 0; %was 0
p.D_f = 0; %0;
%Get parameters
p.maxiter = 20;
if isfield(SET(NO).Flow,'parameter')
  p.curv_f = SET(NO).Flow.parameter.edgeforce;
  p.balloon_f = SET(NO).Flow.parameter.balloonforce;
  p.edge_f = SET(NO).Flow.parameter.edgeforce;
else
  p.curv_f = 2;
  p.balloon_f = 1;
  p.edge_f = 20;
end
p.edgesmootht = 0;
p.maxiter = 300;

%--- Initialize model
n = 80; 
omega = ((2*pi/n)*(1:n))'; %linspace(0,2*pi,n)';
r0 = 2;
p.x = repmat(r0*sin(omega)+(SET(NO).CenterX),[1 1 nslices]);
p.y = repmat(r0*cos(omega)+(SET(NO).CenterY),[1 1 nslices]);

[x,y] = deformit_newlv(p,speedim,e0,e1,e2,e3);

x = repmat(x,[1 SET(NO).TSize 1]);
y = repmat(y,[1 SET(NO).TSize 1]);

%Store the ROI
for rloop=1:length(ind)
  SET(NO).RoiN = SET(NO).RoiN+1;
  SET(NO).Roi(SET(NO).RoiN).X = x(:,:,rloop);
  SET(NO).Roi(SET(NO).RoiN).Y = y(:,:,rloop);
  SET(NO).Roi(SET(NO).RoiN).T = find(not(isnan(x(1,:,rloop))));
  SET(NO).Roi(SET(NO).RoiN).Z = ind(rloop);
  SET(NO).Roi(SET(NO).RoiN).Sign = 1;
  if nargin==0
    SET(NO).Roi(SET(NO).RoiN).Name = 'Flow ROI';
  else
    SET(NO).Roi(SET(NO).RoiN).Name = name;
  end
  
  if isempty(SET(NO).RoiCurrent)
    SET(NO).Roi(SET(NO).RoiN).LineSpec = 'b-';
  else
    SET(NO).Roi(SET(NO).RoiN).LineSpec = SET(NO).Roi(SET(NO).RoiCurrent(end)).LineSpec;
  end
end
SET(NO).RoiCurrent = SET(NO).RoiN;

flowtrackroi_Callback;

if nargin==0
  roi('roisetlabel_Callback',SET(NO).RoiCurrent);
end

%--------------------------------
function flowtrackroi_Callback(m)
%--------------------------------
%Track a flow ROI through all time frames
global DATA SET NO

if isempty(SET(NO).Flow)
  myfailed('No flow data.',DATA.GUI.Segment);
  return;
end

if isempty(SET(NO).Flow.PhaseNo)
  myfailed('Could not find through plane velocity data.',DATA.GUI.Segment);
  return;
end

tempnos=[SET(NO).Flow.PhaseNo,SET(NO).Flow.MagnitudeNo];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

nom = SET(NO).Flow.MagnitudeNo;

if nargin<1
  m = SET(nom).RoiCurrent;
end

if isempty(m)
  myfailed('No ROIs selected.',DATA.GUI.Segment);
  return
end

if any((m>SET(nom).RoiN)|(m==0))
  myfailed('No ROIs available.',DATA.GUI.Segment);
  return;
end

if SET(nom).TSize<8
  myfailed('Too few timeframes.',DATA.GUI.Segment);
  return;
end

tools('enableundo');

mm = m;
for m = mm
if length(SET(nom).Roi(m).T)==1
  roi('roicopyalltimeframes_Callback',m);
end

%--- Take ROI from current time frame or from frame with roi closest to roi
tfrep=[SET(nom).Roi(m).T-SET(NO).TSize SET(nom).Roi(m).T SET(nom).Roi(m).T+SET(nom).TSize];
[tdiff,tindex] = min(abs(tfrep -SET(nom).CurrentTimeFrame));
tf = rem(tfrep(tindex),length(SET(nom).Roi(m).T));
if tf==0
  tf=SET(NO).TSize;
end
xroi = SET(nom).Roi(m).X(:,tf);
yroi = SET(nom).Roi(m).Y(:,tf);

%Calculate image to take
xmin = round(min(xroi(:)));
ymin = round(min(yroi(:)));
xmax = round(max(xroi(:)));
ymax = round(max(yroi(:)));

%Calculate how much extra to take
extra = max(round(0.3*mean([ymax-ymin,xmax-xmin])),4);
xmin = max(xmin-extra,1);
ymin = max(ymin-extra,1);
xmax = min(xmax+extra,SET(nom).XSize);
ymax = min(ymax+extra,SET(nom).YSize);

%--- Check sizes
xsize = xmax-xmin+1;
ysize = ymax-ymin+1;

if isequal(xsize/2,round(xsize/2))
  %even number of rows
  %disp('even rows');
  rowofs = xsize/2+1;
else
  rowofs = (xsize+1)/2;
end

if isequal(ysize/2,round(ysize/2))
  %even number of cols
  %disp('even cols');
  colofs = ysize/2+1;  
else
  colofs = (ysize+1)/2;
end

newxsize = xsize*2;
newysize = ysize*2;
rowofs = xsize+1;
colofs = ysize+1;

if (newxsize<1)||(newysize<1)
  myfailed('Too small ROI to track.',DATA.GUI.Segment);
  return;
end

%--- Loop over timeframes
dx = zeros(1,SET(nom).TSize);
dy = dx;
h = mywaitbarstart(SET(nom).TSize,'Tracking ROI.',5);
for tloop=1:SET(nom).TSize

  %Find next timeframe
  tfp1 = tf+1;
  if tfp1>SET(nom).TSize
    tfp1 = 1;
  end
  
  %--- Extract image
  im = SET(nom).IM(xmin:xmax,ymin:ymax,tf,SET(nom).CurrentSlice);
  imp1 = SET(nom).IM(xmin:xmax,ymin:ymax,tfp1,SET(nom).CurrentSlice);
  
  %--- upsample
  im = imresize(im,[newxsize newysize],'bilinear');
  imp1 = imresize(imp1,[newxsize newysize],'bilinear');

  %--- Calc edge image
  im = sqrt(conv2(im,[1;0;-1],'same').^2+conv2(im,[1 0 -1],'same').^2);
  imp1 = sqrt(conv2(imp1,[1;0;-1],'same').^2+conv2(imp1,[1 0 -1],'same').^2);
  im(1,:)=0; im(end,:)=0; im(:,1)=0; im(:,end)=0;
  imp1(1,:)=0; imp1(end,:)=0; imp1(:,1)=0; imp1(:,end)=0;

  %--- Upsample
  %a = imresize(im,2,'bilinear');
  %b = imresize(imp1,2,'bilinear');
  
  %--- Perform correlation
  c = fftshift(abs(ifft2(fft2(im).*conj(fft2(imp1)))));
  
  %--- Find peak
  [cmax,ind] = max(c(:));
  [row,col] = ind2sub(size(c),ind);

  %Calculate position of peak
  drow=row-rowofs;
  dcol=col-colofs;
  drow = drow/2;
  dcol = dcol/2;
  dx(tfp1) = drow;
  dy(tfp1) = dcol;      
  
  %Update position
  drow = round(drow);
  dcol = round(dcol);
  xmin = min(max(xmin-drow,1),SET(nom).XSize);
  xmax = min(max(xmax-drow,1),SET(nom).XSize);
  ymin = min(max(ymin-dcol,1),SET(nom).YSize);
  ymax = min(max(ymax-dcol,1),SET(nom).YSize);
    
  %Take next timeframe
  tf = tf+1;
  if tf>SET(nom).TSize
    tf=1;
  end
  
  h = mywaitbarupdate(h);
end %Loop over timeframes
mywaitbarclose(h);

x = cumsum(dx);
y = cumsum(dy);
x = x-linspace(0,1,length(x))*(x(end)-x(1));
y = y-linspace(0,1,length(y))*(y(end)-y(1));

x = x-x(SET(nom).CurrentTimeFrame);
y = y-y(SET(nom).CurrentTimeFrame);

%Copy the ROI
for tloop=1:SET(nom).TSize
  SET(nom).Roi(m).X(:,tloop) = xroi-x(tloop);
  SET(nom).Roi(m).Y(:,tloop) = yroi-y(tloop);  
end

%flowshrink_Callback(m,0.8);
old = DATA.ThisFrameOnly;
  
% if SET(nom).TSize<60
%Normal case
tf = SET(nom).CurrentTimeFrame;

DATA.ThisFrameOnly = true;
h = mywaitbarstart(SET(nom).TSize,'Refining ROI',5);
for tloop=1:SET(nom).TSize
  SET(nom).CurrentTimeFrame = tloop;
  flowrefine_Callback(m,true); %true to make it silent
  flowrefine_Callback(m,true); %true to make it silent
  if SET(NO).TSize<=60
    flowrefine_Callback(m,true); %true to make it silent
    flowrefine_Callback(m,true); %true to make it silent
  end
  h = mywaitbarupdate(h);
end

%Expand outward
if isfield(SET(nom).Flow,'parameter')
  expandflowroi(nom,m,SET(nom).Flow.parameter.expandoutward);
end
mywaitbarclose(h);
end
% else
%  %%% Special case...
% 
%   if false
%     cert = zeros(1,SET(nom).TSize);
%     rx = round(mean(SET(nom).RoiX(:,SET(nom).CurrentTimeFrame,m)));
%     ry = round(mean(SET(nom).RoiY(:,SET(nom).CurrentTimeFrame,m)));
%     startintensity = SET(nom).IM((rx-2):(rx+2),(ry-2):(ry+2),SET(nom).CurrentTimeFrame,SET(nom).StartSlice:SET(nom).EndSlice);
%     startintensity = mean(startintensity(:));
% 
%     DATA.ThisFrameOnly = true;
%     h = mywaitbarstart(SET(nom).TSize,'Refining ROI first pass',5);
%     for tloop=1:SET(nom).TSize;
%       SET(nom).CurrentTimeFrame = tloop;
%       rx = round(mean(SET(nom).RoiX(:,SET(nom).CurrentTimeFrame,m)));
%       ry = round(mean(SET(nom).RoiY(:,SET(nom).CurrentTimeFrame,m)));
%       tempintensity = SET(nom).IM((rx-2):(rx+2),(ry-2):(ry+2),SET(nom).CurrentTimeFrame,SET(nom).StartSlice:SET(nom).EndSlice);
%       tempintensity = mean(tempintensity(:));
%       cert(tloop) = max(tempintensity/startintensity,0.7);
%       if tempintensity>0.7*startintensity
%         flowrefine_Callback(m,true); %true to make it silent
%         flowrefine_Callback(m,true); %true to make it silent
%         flowrefine_Callback(m,true); %true to make it silent
%         flowrefine_Callback(m,true); %true to make it silent
%       end;
% 
%       h = mywaitbarupdate(h);
%     end;
%     mywaitbarclose(h);
%     cert = cert-0.7;
% 
%     %--- Move back
%     for tloop=1:SET(nom).TSize
%       SET(nom).RoiX(:,tloop,m) = SET(nom).RoiX(:,tloop,m)+x(tloop);
%       SET(nom).RoiY(:,tloop,m) = SET(nom).RoiY(:,tloop,m)+y(tloop);
%     end;
% 
%     %--- Smooth
%     fx = linspace(-4,4,31); %EH: 2*round(SET(nom).TSize*0.15)+1);
%     f = exp(-fx.^2);
%     f = f./sum(f(:));
%     %figure(88);
%     %plot(f);
% 
%     %figure(72);
%     %plot(cert);
% 
%     cert = [cert cert cert];
%     cert = repmat(cert,size(SET(nom).RoiX,1),1);
%     tempx = [SET(nom).RoiX(:,:,m) SET(nom).RoiX(:,:,m) SET(nom).RoiX(:,:,m)];
%     tempy = [SET(nom).RoiY(:,:,m) SET(nom).RoiY(:,:,m) SET(nom).RoiY(:,:,m)];
%     temp = eps+conv2(cert,f,'same');
%     tempx = conv2(tempx.*cert,f,'same')./temp;
%     tempy = conv2(tempy.*cert,f,'same')./temp;
%     SET(nom).RoiX(:,:,m) = tempx(:,(SET(nom).TSize+1):(SET(nom).TSize*2));
%     SET(nom).RoiY(:,:,m) = tempy(:,(SET(nom).TSize+1):(SET(nom).TSize*2));
% 
%     %--- Move again
%     for tloop=1:SET(nom).TSize
%       SET(nom).RoiX(:,tloop,m) = SET(nom).RoiX(:,tloop,m)-x(tloop);
%       SET(nom).RoiY(:,tloop,m) = SET(nom).RoiY(:,tloop,m)-y(tloop);
%     end;
% 
%     DATA.ThisFrameOnly = true;
%     h = mywaitbarstart(SET(nom).TSize,'Refining ROI second pass',5);
%     for tloop=1:SET(nom).TSize;
%       SET(nom).CurrentTimeFrame = tloop;
%       flowrefine_Callback(m,true); %true to make it silent
%       flowrefine_Callback(m,true); %true to make it silent
%       h = mywaitbarupdate(h);
%     end;
%     mywaitbarclose(h);
%   end;
% end; %TSize<60!? Special hack for Mikael

SET(nom).CurrentTimeFrame = tf;
DATA.ThisFrameOnly = old;
viewfunctions('setview');  %drawfunctions('drawimageno');

%---------------------------------
function flowpropagate_Callback(m) 
%---------------------------------
%Propagate a flow ROI to next time frame and refine
global DATA SET NO
persistent startframe

if isempty(SET(NO).Flow)
  myfailed('No flow data.',DATA.GUI.Segment);
  return;
end

if isempty(SET(NO).Flow.PhaseNo)
  myfailed('Could not find through plane velocity data.',DATA.GUI.Segment);
  return;
end

tempnos=[SET(NO).Flow.PhaseNo,SET(NO).Flow.MagnitudeNo];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

nom=SET(NO).Flow.MagnitudeNo;

if nargin<1
  m = SET(nom).RoiCurrent;
end

if isempty(m)
  myfailed('No ROIs selected.',DATA.GUI.Segment);
  return
end

if any((m>SET(nom).RoiN)|(m==0))
  myfailed('No ROIs available.',DATA.GUI.Segment);
  return;
end

if isempty(startframe)||isnan(startframe)||DATA.Testing
  startframe = SET(nom).CurrentTimeFrame;
else
  if isequal(SET(nom).CurrentTimeFrame,startframe)
    startframe = [];
    mywarning('Heart cycle completed',DATA.GUI.Segment);
    segment('updateflow');
    return;
  end
end
 
%Find list of affected image stacks
nos = [];
if ~isempty(SET(NO).Flow)
  nos = [nos ...
    SET(NO).Flow.MagnitudeNo ...    
    SET(NO).Flow.PhaseNo ...
    SET(NO).Flow.PhaseX ...
    SET(NO).Flow.PhaseY ...
    SET(NO).Flow.Angio ...
    SET(NO).Flow.VelMag];
end

tools('enableundo');

for roin = m
  
  %find starting timeframe, which must contain a roi
  tfrep=[SET(nom).Roi(roin).T-SET(NO).TSize SET(nom).Roi(roin).T SET(nom).Roi(roin).T+SET(nom).TSize];
  [tdiff,tindex] = min(abs(tfrep -SET(nom).CurrentTimeFrame));
  tm1 = rem(tfrep(tindex),SET(NO).TSize);
  t=tm1+1;
  if tm1==0
    tm1=SET(NO).TSize;
  end
  % t = SET(nom).CurrentTimeFrame+1;
  % tm1 = SET(nom).CurrentTimeFrame;
  % if t>SET(nom).TSize
  %   t = 1;
  %   tm1 = SET(nom).TSize;
  % end;
  
  
  %Get parameters
  p.maxiter = 20;
  if isfield(SET(nom).Flow,'parameter')
    p.curv_f = SET(nom).Flow.parameter.curveforce;
    p.balloon_f = SET(nom).Flow.parameter.balloonforce;
    p.edge_f = SET(nom).Flow.parameter.edgeforce;
  else
    p.curv_f = 2;
    p.balloon_f = 1;
    p.edge_f = 20;
  end
  
  p.M_f = 0;
  p.D_f = 0;
  p.Z_f = 0;
  p.along = 1;
  p.across = 1;
  p.unsignededge = 0;
  if ~DATA.EndoEdgeDetected
    lv('findedges',1,p,nom); %Use segmentation parameters 1==endo
  end
  
  %Get model
  z = SET(nom).Roi(roin).Z;
  p.x = SET(nom).Roi(roin).X(1:(DATA.NumPoints-1),tm1);
  p.y = SET(nom).Roi(roin).Y(1:(DATA.NumPoints-1),tm1);
  
  %Calculate balloon image, ignore p.centerintensity since endocardium
  rx = round(mean(p.x(:)));
  ry = round(mean(p.y(:)));
  temp = SET(nom).IM((rx-2):(rx+2),(ry-2):(ry+2),t,z);
  meanint = mean(temp(:));
  speedim = SET(nom).IM(:,:,t,z)-single(meanint-0.2);
  
  [p.x,p.y] = deformit_newlv(p,...
    speedim,...
    DATA.EndoEDGE0(:,:,t,z),...
    DATA.EndoEDGE1(:,:,t,z),...
    DATA.EndoEDGE2(:,:,t,z),...
    DATA.EndoEDGE3(:,:,t,z));
  
  %Store the result
  SET(nom).Roi(roin).X(1:(DATA.NumPoints-1),t) = p.x;
  SET(nom).Roi(roin).Y(1:(DATA.NumPoints-1),t) = p.y;
  SET(nom).Roi(roin).X(end,t) = p.x(1,:);
  SET(nom).Roi(roin).Y(end,t) = p.y(1,:);
  SET(nom).Roi(roin).T=find(not(isnan(SET(nom).Roi(roin).X(1,:))));
  
  SET(nom).CurrentTimeFrame = t;
  for loop = 1:length(nos)
    SET(nos(loop)).CurrentTimeFrame = t;
  end
  [~,SET(nom).Roi(roin).Area] = ...
    calcfunctions('calcroiarea',nom,roin);
  [mn,sd]=calcfunctions('calcroiintensity',nom,roin);
  SET(nom).Roi(roin).Mean = mn;
  SET(nom).Roi(roin).StD = sd;
end

%update flow result panel
if ~isempty(DATA.FlowNO) && ismember(DATA.FlowNO,SET(nom).Linked) && ~isempty(DATA.FlowROI)
  if SET(nom).RoiN <1
    DATA.FlowROI = [];
  else
    calcfunctions('calcflow',nom);
  end
end
viewfunctions('switchtimeframe',0, DATA.Synchronize); 
% % DATA.updateaxestables('area',nom);
% panels = find(ismember(DATA.ViewPanels,SET(nom).Linked));
% for p = panels
%   drawfunctions('drawroi',p);
% end
% segment('updateflow');
% viewfunctions('updatetimebars');

%-------------------------------------
function flowrefine_Callback(m,silent)
%-------------------------------------
%Refine segmentation of a flow ROI
global DATA SET NO

if isempty(SET(NO).Flow)
  myfailed('No flow data.',DATA.GUI.Segment);
  return;
end

nom = SET(NO).Flow.MagnitudeNo;
nop = SET(NO).Flow.PhaseNo;

if isempty(SET(nom).Flow.PhaseNo)
  myfailed('Could not find through plane velocity data.',DATA.GUI.Segment);
  return;
end

tempnos=[SET(NO).Flow.PhaseNo,SET(NO).Flow.MagnitudeNo];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

if nargin<1
  m = SET(nom).RoiCurrent;
end

if nargin<2
  silent = false;
end

if any((m>SET(nom).RoiN)|(m==0))
  myfailed('No ROIs available.',DATA.GUI.Segment);
  return;
end

% if strcmp(DATA.ProgramName,'Segment CMR')
%   DATA.ThisFrameOnly = false;
% end

if ~silent
  tools('enableundo');
end

for roin = m
  
  %Get parameters
  p.maxiter = 20;
  if isfield(SET(nom).Flow,'parameter')
    p.curv_f = SET(nom).Flow.parameter.curveforce;
    p.balloon_f = SET(nom).Flow.parameter.balloonforce;
    p.edge_f = SET(nom).Flow.parameter.edgeforce;
  else
    p.curv_f = 2;
    p.balloon_f = 1;
    p.edge_f = 20;
  end
  p.M_f = 0;
  p.D_f = 0;
  p.Z_f = 0;
  p.along = 1;
  p.across = 1;
  p.unsignededge = 0;
  if ~DATA.EndoEdgeDetected
    lv('findedges',1,p,nom); %Use segmentation parameters 1==endo
  end
  
  if DATA.ThisFrameOnly
    timeframes = SET(nom).CurrentTimeFrame;
    p.balloon_f = p.balloon_f;
    if isempty(find(SET(nom).Roi(roin).T==timeframes,1))
      myfailed('To be able to refine flow ROI, ROI must be present in the current timeframe');
      return;
    end
  else
    timeframes = 1:SET(nom).TSize;
    if not(isequal(SET(nom).Roi(roin).T,timeframes))
      myfailed('To be able to refine flow ROI, ROI must be present in all timeframes or single frame mode should be used');
      return;
    end
  end
  
  %Get model
  p.x = SET(nom).Roi(roin).X(1:(DATA.NumPoints-1),timeframes);
  p.y = SET(nom).Roi(roin).Y(1:(DATA.NumPoints-1),timeframes);
  
  if length(timeframes)>1
    p.D_f = 0.5;
    p.M_f = 1;
    p.maxiter = 200;
  end
  
  % Set range so that calcbint does not return NaN, which happens if they
  % are not set!
  SET(nom).StartSlice = SET(nom).CurrentSlice;
  SET(nom).EndSlice = SET(nom).CurrentSlice;
  
  %Calculate balloon image, ignore p.centerintensity since "endocardium"
  if DATA.ThisFrameOnly
    rx = round(mean(p.x(:)));
    ry = round(mean(p.y(:)));
    temp = SET(nom).IM((rx-2):(rx+2),(ry-2):(ry+2),timeframes,SET(nom).StartSlice:SET(nom).EndSlice);
    meanint = mean(temp(:));
  else
    meanint=rv('calcbint',mean(p.x(:)),mean(p.y(:)));
  end
  %calcballoon(DATA.BpInt,p);
  %DATA.BpInt
  DATA.BALLOON = SET(nom).IM-single(meanint-0.2);
  
  %Start deform (mex'ed)
  if not(DATA.ThisFrameOnly)
    if isequal(size(DATA.BALLOON),size(DATA.EndoEDGE0))
      [p.x,p.y] = deformit_newlv(p,...
        DATA.BALLOON,...
        DATA.EndoEDGE0,...
        DATA.EndoEDGE1,...
        DATA.EndoEDGE2,...
        DATA.EndoEDGE3);
    else
      DATA.EndoEDGE0 = repmat(DATA.EndoEDGE0,[1 1 1 size(DATA.BALLOON,4)]);
      DATA.EndoEDGE1 = repmat(DATA.EndoEDGE1,[1 1 1 size(DATA.BALLOON,4)]);
      DATA.EndoEDGE2 = repmat(DATA.EndoEDGE2,[1 1 1 size(DATA.BALLOON,4)]);
      DATA.EndoEDGE3 = repmat(DATA.EndoEDGE3,[1 1 1 size(DATA.BALLOON,4)]);
      if isequal(size(DATA.BALLOON),size(DATA.EndoEDGE0))
        [p.x,p.y] = deformit_newlv(p,...
          DATA.BALLOON,...
          DATA.EndoEDGE0,...
          DATA.EndoEDGE1,...
          DATA.EndoEDGE2,...
          DATA.EndoEDGE3);
      end
    end
  else
    if isequal(size(DATA.BALLOON),size(DATA.EndoEDGE0))
      [p.x,p.y] = deformit_newlv(p,...
        DATA.BALLOON(:,:,timeframes,:),...
        DATA.EndoEDGE0(:,:,timeframes,:),...
        DATA.EndoEDGE1(:,:,timeframes,:),...
        DATA.EndoEDGE2(:,:,timeframes,:),...
        DATA.EndoEDGE3(:,:,timeframes,:));
    else
      DATA.EndoEDGE0 = repmat(DATA.EndoEDGE0,[1 1 1 size(DATA.BALLOON,4)]);
      DATA.EndoEDGE1 = repmat(DATA.EndoEDGE1,[1 1 1 size(DATA.BALLOON,4)]);
      DATA.EndoEDGE2 = repmat(DATA.EndoEDGE2,[1 1 1 size(DATA.BALLOON,4)]);
      DATA.EndoEDGE3 = repmat(DATA.EndoEDGE3,[1 1 1 size(DATA.BALLOON,4)]);
      if isequal(size(DATA.BALLOON),size(DATA.EndoEDGE0))
        [p.x,p.y] = deformit_newlv(p,...
          DATA.BALLOON(:,:,timeframes,:),...
          DATA.EndoEDGE0(:,:,timeframes,:),...
          DATA.EndoEDGE1(:,:,timeframes,:),...
          DATA.EndoEDGE2(:,:,timeframes,:),...
          DATA.EndoEDGE3(:,:,timeframes,:));
      end
    end
  end
  
  %Store the result
  SET(nom).Roi(roin).X(1:(DATA.NumPoints-1),timeframes) = p.x;
  SET(nom).Roi(roin).Y(1:(DATA.NumPoints-1),timeframes) = p.y;
  SET(nom).Roi(roin).X(end,timeframes) = p.x(1,:,:);
  SET(nom).Roi(roin).Y(end,timeframes) = p.y(1,:,:);
  
  if ~silent
    [~,SET(nom).Roi(roin).Area] = ...
      calcfunctions('calcroiarea',nom,roin);
    [mn,sd]=calcfunctions('calcroiintensity',nom,roin);
    SET(nom).Roi(roin).Mean = mn;
    SET(nom).Roi(roin).StD = sd;
  end
end
calcfunctions('calcflow',nom);
if ~silent
  segment('updateflow');
%     DATA.updateaxestables('area',nom);
  %Get all linked stacks and draw the rois in them as well.
  panels = find(ismember(DATA.ViewPanels,SET(nom).Linked));
  for p = panels
    drawfunctions('drawroi',p);
  end
end

%-------------------------------
function flowplotquiver_Callback 
%-------------------------------
%Plots a quiver plot of current image stack
global SET NO DATA

if isempty(SET(NO).Flow)
  myfailed('Current image stack does not contain flow data.',DATA.GUI.Segment);
  return;
end

%Find image stacks.
nom = SET(NO).Flow.MagnitudeNo;
nox = SET(NO).Flow.PhaseX;
noy = SET(NO).Flow.PhaseY;

if isempty(nom)||isempty(nox)||isempty(noy)
  myfailed('Flow is not bidirectional.',DATA.GUI.Segment);
  return;
end

tempnos=[nom nox noy];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

cs = SET(nom).CurrentSlice;
ct = SET(nom).CurrentTimeFrame;

vx = SET(nox).IM(:,:,ct,cs)-0.5;
vy = SET(noy).IM(:,:,ct,cs)-0.5;
mag = SET(nom).IM(:,:,ct,cs);

if ~isnan(SET(nom).EndoX(1,ct,cs))
  mask = segment('createmask',...
    [SET(nom).XSize SET(NO).YSize],...
    SET(nom).EndoY(:,ct,cs),...
    SET(nom).EndoX(:,ct,cs));
  vx(~mask) = NaN;
  vy(~mask) = NaN;
end

figure(12);
imagesc(mag);
colormap(gray);

hold on;
quiver(vx,vy,3);
plot(SET(nom).EndoY(:,ct,cs),SET(nom).EndoX(:,ct,cs),'r-');
hold off;
axis image off;

%--------------------------------
function flowflipdirs_Callback
%-------------------------------
%Flips LR and UP direction
global SET NO DATA

if isempty(SET(NO).Flow)
  myfailed('Current image stack does not contain flow data.',DATA.GUI.Segment);
  return;
end

%Find image stacks.
nom = SET(NO).Flow.MagnitudeNo;
nox = SET(NO).Flow.PhaseX;
noy = SET(NO).Flow.PhaseY;

if isempty(nom)||isempty(nox)||isempty(noy)
  myfailed('Flow is not bidirectional.',DATA.GUI.Segment);
  return;
end

tempnos=[nom nox noy];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

temp = SET(nox).IM;
SET(nox).IM = SET(noy).IM;
SET(noy).IM = temp;

%segment('viewrefresh_Callback');
viewfunctions('setview');

%----------------------------------
function flowflipupdown_Callback
%--------------------------------
%Flips up/down direction'
global SET NO DATA

if isempty(SET(NO).Flow)
  myfailed('Current image stack does not contain flow data.',DATA.GUI.Segment);
  return;
end

%Find image stacks.
nom = SET(NO).Flow.MagnitudeNo;
nox = SET(NO).Flow.PhaseX;
noy = SET(NO).Flow.PhaseY;

if isempty(nom)||isempty(nox)||isempty(noy)
  myfailed('Flow is not bidirectional.',DATA.GUI.Segment);
  return;
end

tempnos=[nom nox noy];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

SET(noy).IM = 1-SET(noy).IM;

%segment('viewrefresh_Callback');
viewfunctions('setview')

%--------------------------------------
function flowflipleftright_Callback
%------------------------------------
%Flips up/down direction'
global SET NO DATA

if isempty(SET(NO).Flow)
  myfailed('Current image stack does not contain flow data.',DATA.GUI.Segment);
  return;
end

%Find image stacks.
nom = SET(NO).Flow.MagnitudeNo;
nox = SET(NO).Flow.PhaseX;
noy = SET(NO).Flow.PhaseY;

if isempty(nom)||isempty(nox)||isempty(noy)
  myfailed('Flow is not bidirectional.',DATA.GUI.Segment);
  return;
end

tempnos=[nom nox noy];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

SET(nox).IM = 1-SET(nox).IM;

%segment('viewrefresh_Callback');
viewfunctions('setview')

%------------------
function flowfix(n)
%------------------
%temp function to flip 2d flow directions
global SET NO DATA

if isempty(SET(NO).Flow)
  myfailed('Current image stack does not contain flow data.',DATA.GUI.Segment);
  return;
end

%Find image stacks.
nom = SET(NO).Flow.MagnitudeNo;
nox = SET(NO).Flow.PhaseX;
noy = SET(NO).Flow.PhaseY;

if isempty(nom)||isempty(nox)||isempty(noy)
  myfailed('Flow is not bidirectional.',DATA.GUI.Segment);
  return;
end

tempnos=[nom nox noy];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

if n>4
  n = n-4;
  flowflipdirs_Callback;
end

switch n
  case 1
    %nothing
  case 2
    flowflipupdown_Callback;
  case 3
    flowflipleftright_Callback;    
  case 4
    flowflipupdown_Callback;    
    flowflipleftright_Callback;
  otherwise
    myfailed('Invalid argument.',DATA.GUI.Segment);
    return;
end

%-----------------------------------
function setflowpreferences_Callback
%-----------------------------------
%open a box where it is possibly to set the ROI propagation parameters:
%balloon force, edge force and curve force

%written by Helen Soneson 2008-08-01

global SET NO DATA

s = [];
if ~isfield(SET(NO).Flow,'parameter')
  s.BalloonForce = 1;
  s.EdgeForce = 20;
  s.CurveForce = 10;  
  s.ExpandOutward = 2.3;
else
  s.BalloonForce = SET(NO).Flow.parameter.balloonforce;
  s.EdgeForce = SET(NO).Flow.parameter.edgeforce;
  s.CurveForce = SET(NO).Flow.parameter.curveforce;
  s.ExpandOutward = SET(NO).Flow.parameter.expandoutward;
end

[s,ok] = inputstruct(s,'ROI propagation');
if not(ok)
  myfailed('Failed to set the parameters',DATA.GUI.Segment);
  return;
else
  SET(NO).Flow.parameter.balloonforce = max(0,min(100,s.BalloonForce));
  SET(NO).Flow.parameter.edgeforce = max(0,min(100,s.EdgeForce));
  SET(NO).Flow.parameter.curveforce = max(0,min(100,s.CurveForce));
  SET(NO).Flow.parameter.expandoutward = max(-20,min(20,s.ExpandOutward));
end

%---------------------
function [ok,pulmno,pulmroi,aortno,aortroi] = flowfindqpqs  
%---------------------
%Find data sets with pulmonary and aortic flow, to use for Qp/Qs analysis
global SET NO DATA

ok = false;
%Find data set with pulmonary flow
pulmno = NaN;
for loop=1:length(NO)
  for rloop=1:length(SET(loop).RoiN)
    if isequal(lower(SET(loop).Roi(rloop).Name),'pulmonary flow')
      pulmno = loop;
      pulmroi = rloop;
    end
  end
end

%Find data set with aortic flow
aortno = NaN;
for loop=1:length(NO)
  for rloop=1:length(SET(loop).RoiN)
    if isequal(lower(SET(loop).Roi(rloop).Name),'aortic ascending flow')
      aortno = loop;
      aortroi = rloop;
    end
  end
end

if (not(isnan(pulmno)))&&(not(isnan(aortno)))
else
  %They do not both exist
  myfailed('Could not find both pulmonary flow and ascending aortic flow.',DATA.GUI.Segment);
  return;
end

%-----------------------------
function expandflowroi(no,m,d)
%-----------------------------
%Expand the roi M outwards by D mm's

global SET

if nargin<3
  myfailed('Expected three input arguments.');
  return;
end

nom = SET(no).Flow.MagnitudeNo;

dsdx = conv2(cat(1,SET(nom).Roi(m).X(end,:),SET(nom).Roi(m).X(:,:),SET(nom).Roi(m).X(1,:)),...
  [1;0;-1],'valid');
dsdy = conv2(cat(1,SET(nom).Roi(m).Y(end,:),SET(nom).Roi(m).Y(:,:),SET(nom).Roi(m).Y(1,:)),...
  [1;0;-1],'valid');

%Make rotation by a simple mirror operation
xprim = -dsdy;
yprim = dsdx;

%Normalize
primnorm = sqrt(xprim.^2+yprim.^2);
xprim = xprim./primnorm;
yprim = yprim./primnorm;

%Add distance
SET(nom).Roi(m).X = SET(nom).Roi(m).X+(d/SET(nom).ResolutionX)*xprim;
SET(nom).Roi(m).Y = SET(nom).Roi(m).Y+(d/SET(nom).ResolutionY)*yprim;

%Force last & first point to be equal
mx = 0.5*(SET(nom).Roi(m).X(1,:)+SET(nom).Roi(m).X(end,:));
my = 0.5*(SET(nom).Roi(m).Y(1,:)+SET(nom).Roi(m).Y(end,:));
SET(nom).Roi(m).X(1,:) = mx;
SET(nom).Roi(m).Y(1,:) = my;

%---------------------------------------------------------------
function flow3dmovie_Callback(arg,scaling,no,up,down,left,right)
%---------------------------------------------------------------
%Displays 3D movie of velocity in a vessel.
global DATA SET NO
persistent handles

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

if not(ischar(arg))

  if isfield(handles,'fig') && ishandle(handles.fig) 
    figure(handles.fig);
  else
    %--- Init
    fig = mygui('flow3movie.fig');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = fig.handles;
    handles.fig = handles.figure1;
    handles.tf = 1;
    handles.timeframes = SET(NO).StartAnalysis:SET(NO).EndAnalysis;
    handles.mask = arg;
    handles.scaling = scaling;
    handles.no = no;
    handles.up = up;
    handles.down = down;
    handles.left = left;
    handles.right = right;
    handles.x = (handles.up:handles.down)-handles.up+1;
    handles.y = (handles.left:handles.right)-left+1;  
    handles.surf = [];
  end
  flow3dmovie_Callback('update');
else
  switch arg
    case 'next'
      if handles.tf<length(handles.timeframes)
        handles.tf = handles.tf+1;
      else
        handles.tf = 1;
      end
      flow3dmovie_Callback('update');
    case 'prev'
      if handles.tf>1
        handles.tf = handles.tf-1;
      else
        handles.tf = length(handles.timeframes);
      end
      flow3dmovie_Callback('update');      
    case 'update'
      if isempty(SET(handles.no).Flow.PhaseCorr)
        %No phase correction
        temp = SET(NO).Roi(1).Sign*handles.mask(:,:,handles.tf).*...
          (SET(handles.no).IM(:,:,handles.timeframes(handles.tf),SET(NO).Roi(1).Z)-0.5)*2*SET(SET(NO).Flow.PhaseNo).VENC;        
      else
        slice = SET(handles.no).CurrentSlice;
        if SET(handles.no).Flow.PhaseCorrTimeResolved
          %Time resolved phase correction          
          temp = SET(NO).Roi(1).Sign*handles.mask(:,:,handles.tf).*...
            (SET(handles.no).IM(:,:,handles.timeframes(handles.tf),SET(NO).Roi(1).Z)-0.5-...
            SET(handles.no).Flow.PhaseCorr(:,:,handles.tf,slice))*2*SET(SET(NO).Flow.PhaseNo).VENC;          
        else
          %Stationary phase correction
          temp = SET(NO).Roi(1).Sign*handles.mask(:,:,handles.tf).*...
            (SET(handles.no).IM(:,:,handles.timeframes(handles.tf),SET(NO).Roi(1).Z)-0.5-...
            SET(handles.no).Flow.PhaseCorr)*2*SET(SET(NO).Flow.PhaseNo).VENC;
        end
      end
      temp = double(temp);
      temp = temp(handles.up:handles.down,handles.left:handles.right);
      if isempty(handles.surf)
        handles.surf = surf(handles.plotaxes,handles.y,handles.x,0.25*temp*handles.scaling);
        axis image off;
        cameratoolbar(handles.fig,'show');   
        set(handles.plotaxes,'clim',...
          0.25*handles.scaling*[-SET(SET(NO).Flow.PhaseNo).VENC SET(SET(NO).Flow.PhaseNo).VENC]);
        set(handles.surf,'facelighting','phong');
      else
        set(handles.surf,...
          'cdata',0.25*temp*handles.scaling,...
          'zdata',0.25*temp*handles.scaling);
      end
      if DATA.Record
        drawnow;
        DATA.MovieFrame = mygetframe(handles.fig);
        export('exportmovierecorder_Callback','newframe');
      end

  end
end


%----------------------------
function flowsetvenc_Callback 
%----------------------------
%Set VENC for current image stack.
global SET NO DATA

if isempty(SET(NO).Flow)
  myfailed('No flow data.',DATA.GUI.Segment);
  return;
end

nos = [SET(NO).Flow.MagnitudeNo ...
  SET(NO).Flow.PhaseNo ...
  SET(NO).Flow.PhaseX ...
  SET(NO).Flow.PhaseY ...
  SET(NO).Flow.Angio ...
  SET(NO).Flow.VelMag ...
  ];

%nop = SET(NO).Flow.PhaseNo;
%nox = SET(NO).Flow.PhaseX;
%noy = SET(NO).Flow.PhaseY;
%nom = SET(NO).Flow.MagnitudeNo;

s = myinputdlg({'Enter VENC (velocity encoding range) [cm/s]'},'VENC',1,{sprintf('%0.5g',SET(nos(1)).VENC)});
if isempty(s)
  return;
else
  [venc,ok] = str2num(s{1}); 
  if not(ok)
    myfailed('Invalid VENC.',DATA.GUI.Segment);
    return;
  end
  if (venc<0)
    rangestr = '> 0';
    errmsg = dprintf('Needs to be %s',rangestr);
    myfailed(errmsg,DATA.GUI.Segment);
    return;
  end

  for noloop=nos
    SET(noloop).VENC = venc;
    % update flow report
    segment('updateflow');
  end
end

%-----------------------------------------------
function [varargout] = flowvesselenhancement(no) 
%-----------------------------------------------
%Beta function to show vessels with time depending flow.
global DATA SET NO

if nargin==0
  if isempty(SET(NO).Flow)
    myfailed('No flow data.',DATA.GUI.Segment);
    return;
  end

  if isempty(SET(NO).Flow.PhaseNo)
    myfailed('No through plane velocity data found.',DATA.GUI.Segment);
    return;
  end
  
  %Find index to phase image
  no = SET(NO).Flow.PhaseNo;
  
end

tempnos=[no NO];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

%Reserve memory
speedim = repmat(single(0),[SET(NO).XSize SET(NO).YSize SET(NO).ZSize]);

%Set up time to do Fourier analysis
t = linspace(0,2*pi,SET(no).TSize);
t = reshape(t,[1 1 SET(no).TSize]);
sint = sin(t);
cost = cos(t);
sint = repmat(sint,[SET(no).XSize SET(no).YSize 1]);
cost = repmat(cost,[SET(no).XSize SET(no).YSize 1]);

%--- Loop over the slices
for zloop=1:SET(NO).ZSize
  
  %Extract image
  im = SET(no).IM(:,:,:,zloop)-0.5; %Subtract with 0.5 to get zero velocity to zero.
  imm = mean(SET(NO).IM(:,:,:,zloop),3); %Take mean over timeframes
  
  %- Do fourier analysis
  %V = v0+a*sin(t)+b*cos(t)+d*sin(2*t)+e*cos(2*t);
  ima = sum(im.*sint,3);
  imb = sum(im.*cost,3);
  imfas = sum(sqrt(ima.^2+imb.^2),3);
  imfas(imm<0.05) = 0;
  speedim(:,:,zloop) = imfas.*imm; %Store
end

if nargout>0
  varargout = cell(1,1);
  varargout{1} = speedim;
else
  figure(15);
  imagesc(speedim(:,:,SET(NO).CurrentSlice));
  axis image off;
  set(15,'Numbertitle','off','name','Vessel enhancement.');
end

%--------------------------------
function flowcreateangio_Callback 
%--------------------------------
%Create angio image
flowcreate_helper(true); %angio

%---------------------------------
function flowcreatevelmag_Callback
%---------------------------------
%Create velocity magnitude image.
flowcreate_helper(false); %normal mag no angio

%------------------------------
function flowcreate_helper(angio)
%------------------------------
%Helper function to create flow derived image stacks.
global DATA SET NO

if isempty(SET(NO).Flow)
  myfailed('No velocity data set.',DATA.GUI.Segment);
  return;
end

if  angio && ~isempty(SET(NO).Flow.Angio)
  myfailed('Angio already exists.',DATA.GUI.Segment);
  return;
elseif  ~angio && ~isempty(SET(NO).Flow.VelMag)
  myfailed('Velocity Magnitude already exists.',DATA.GUI.Segment);
  return;
end

tempnos=[SET(NO).Flow.PhaseNo SET(NO).Flow.PhaseX SET(NO).Flow.PhaseY];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

k = length(SET)+1;
SET(k) = SET(NO);
temp = 0*SET(NO).IM;

%Add the different velocity coding directions
if ~isempty(SET(NO).Flow.PhaseNo)
  temp = (SET(SET(NO).Flow.PhaseNo).IM-0.5).^2;
end
if ~isempty(SET(NO).Flow.PhaseX)
  temp = (SET(SET(NO).Flow.PhaseX).IM-0.5).^2;
end
if ~isempty(SET(NO).Flow.PhaseY)
  temp = (SET(SET(NO).Flow.PhaseY).IM-0.5).^2;
end
temp = sqrt(temp);

%Normalize
tempmin = min(temp(:));
temp = temp-tempmin;
tempmax = max(temp(:));
temp = temp./tempmax;

if angio
  temp = temp.*SET(SET(NO).Flow.MagnitudeNo).IM;
end

%Create references
nos = [...
  SET(NO).Flow.MagnitudeNo ...
  SET(NO).Flow.PhaseNo ...
  SET(NO).Flow.PhaseX ...
  SET(NO).Flow.PhaseY ...
  SET(NO).Flow.VelMag ...
  SET(NO).Flow.Angio];
if angio
  for loop=1:length(nos)
    SET(nos(loop)).Flow.Angio = k;
  end  
else
  for loop=1:length(nos)
    SET(nos(loop)).Flow.VelMag = k;
  end
end

%Normalize and store
minv = min(temp(:));
temp = temp-minv;
maxv = max(temp(:));

SET(k).IntensityOffset = 0;
SET(k).IntensityScaling = maxv;
SET(k).IM = temp/maxv;
SET(k).Flow=SET(SET(k).Flow.MagnitudeNo).Flow;

if angio
  SET(k).ImageType = 'Angio image';
else
  SET(k).ImageType = 'Velocity magnitude';
end

viewfunctions('setview');  %drawfunctions('drawimageno',NO);

%--------------------------------
function flowdeleteangio_Callback 
%--------------------------------
%Delete angio image
flowdelete_helper(true); %angio

%---------------------------------
function flowdeletevelmag_Callback 
%---------------------------------
%Delete velocity magnitude image
flowdelete_helper(false); %normal mag no angio

%------------------------------
function flowdelete_helper(angio)
%------------------------------
%Helper function to delete angio images.
global DATA SET NO

if isempty(SET(NO).Flow)
  myfailed('No velocity data set.',DATA.GUI.Segment);
  return;
end

if  angio && isempty(SET(NO).Flow.Angio)
  myfailed('No Angio image exists.',DATA.GUI.Segment);
  return;
elseif  ~angio && isempty(SET(NO).Flow.VelMag)
  myfailed('No Velocity Magnitude image exists.',DATA.GUI.Segment);
  return;
end

flow = SET(NO).Flow;

nos = [...
  SET(NO).Flow.MagnitudeNo ...
  SET(NO).Flow.PhaseNo ...
  SET(NO).Flow.PhaseX ...
  SET(NO).Flow.PhaseY ...
  SET(NO).Flow.VelMag ...
  SET(NO).Flow.Angio];

% shuffle down.
if angio
  kill=SET(NO).Flow.Angio;
  flow.Angio=[];
  if SET(NO).Flow.VelMag>kill
    flow.VelMag=flow.VelMag-1;
  end
else
  kill=SET(NO).Flow.VelMag;
  flow.VelMag=[];
  if SET(NO).Flow.Angio>kill
    flow.Angio=flow.Angio-1;
  end
end

for loop=1:length(nos)
  SET(nos(loop)).Flow = flow;
end 

ind=true(size(SET));
ind(kill)=false;
SET = SET(ind);

NO=kill;
%Take away from viewpanels
zer = zeros(1,length(DATA.ViewPanels));
tempind = not(DATA.ViewPanels==NO);
DATA.ViewPanels = DATA.ViewPanels(tempind);
DATA.ViewPanelsType = DATA.ViewPanelsType(tempind);
DATA.ViewPanelsMatrix = DATA.ViewPanelsMatrix(tempind);
DATA.ViewIM = DATA.ViewIM(tempind);
DATA.ViewMatrix = [];

%DATA.ViewPanels have internal references to old image stack ordering.
for loop=1:length(DATA.ViewPanels)
  temp = DATA.ViewPanels(loop);
  if temp>0
    newpos = zer; %zer frome above;
    newpos(temp) = 1;
    newpos = newpos(tempind); %tempind from above
    newpos = find(newpos);
    if ~isempty(newpos)
      DATA.ViewPanels(loop) = newpos;
    else
      DATA.ViewPanels(loop) = 0;
    end
  end
end

%Update datasetpreview
ind = repmat(ind,DATA.GUISettings.ThumbnailSize,1);
% vertical
DATA.DATASETPREVIEW = DATA.DATASETPREVIEW(ind(:)',:);

NO = NO-1;
if NO<1
  NO = 1;
end

DATA.CurrentPanel = DATA.CurrentPanel-1;
if DATA.CurrentPanel<1
  DATA.CurrentPanel = 1;
end

%Switch data does some backup before...
SET(NO).StartSlice = SET(NO).StartSlice;
SET(NO).EndSlice = SET(NO).EndSlice;
SET(NO).CurrentTimeFrame = SET(NO).CurrentTimeFrame;

drawfunctions('drawno',NO)
% drawfunctions('drawall'); %update
% DATA.switchtoimagestack(NO,true); %force

%------------------------
function coupleflowstacks
%------------------------
%Couple matching magnitude and phase stacks of each loaded flow study.
global SET

magnos = [];
phasenos = [];
flownos = [];
[s,i,j] = unique(round([SET.AcquisitionTime]));
for t = 1:numel(s)
  nos = find(j == t);
  nos = nos(cellfun(@isempty,{SET(nos).Flow}))';
  flownos = [flownos nos];
  if numel(nos) == 2 && ~all(ismember(nos,SET(nos(1)).Linked))
    %Found a pair, do some checks
    if isequal(size(SET(nos(1)).IM),size(SET(nos(2)).IM)) && ...
        isequal(SET(nos(1)).ImagePosition,SET(nos(2)).ImagePosition) && ...
        isequal(SET(nos(1)).ImageOrientation,SET(nos(2)).ImageOrientation)
    
    %See if we can find which one is magnitude
    ix = regexp(SET(nos(2)).SeriesDescription,SET(nos(1)).SeriesDescription,'end');
    if isempty(ix)
      ix = regexp(SET(nos(1)).SeriesDescription,SET(nos(2)).SeriesDescription,'end');
      if ~isempty(ix) && ...
          ~isempty(regexp(SET(nos(1)).SeriesDescription(ix+1:end),'P', 'once'))
        phasenos = [phasenos nos(1)];
        magnos = [magnos nos(2)];
      end
    else
      if ~isempty(regexp(SET(nos(2)).SeriesDescription(ix+1:end),'P', 'once'))
        phasenos = [phasenos nos(2)];
        magnos = [magnos nos(1)];
      end
    end
    if isempty(ix) 
      %No indication found of which stack is phase
      %Guess that mean of IM closest to 1/2 is phase
      no1mean = mean(SET(nos(1)).IM(:));
      no2mean = mean(SET(nos(2)).IM(:));
      [~,ind] = sort(abs([no1mean no2mean]-0.5));
      phasenos = [phasenos nos(ind(1))];
      magnos = [magnos nos(ind(2))];
    end
    end
  end
end

if isempty(magnos)
  if length(unique(flownos)) > 1
    %ask user to manually define magnitude and phase image stack
    stri = dprintf('Define stacks out of  [');
    for loop = 1:length(flownos)
      stri = [stri sprintf(' %d ',flownos(loop))];
    end
    stri =[stri sprintf(']\n\n')];
    stri1 = dprintf('Magnitude stack');
    stri2 = dprintf('Phase stack');
    stri3 = dprintf('Set stacks');
    nosman = myinputdlg({[stri stri1], stri2}, stri3, 1);
    if isempty(nosman)
      return;
    end
    magnos = str2num(nosman{1});
    phasenos = str2num(nosman{2});
    if isempty(magnos) || isempty(phasenos) || isempty(find(flownos==magnos)) || isempty(find(flownos==phasenos))
      myfailed('Incorrect entered image stack numbers');
      return;
    end
  else
    myfailed('Found no stacks to couple');
    return
  end
end

magnostri = '';
phasenostri = '';
for i = 1:numel(magnos)
  magno = magnos(i);
  phaseno = phasenos(i);
  
  %Update
  flowstruct = struct('MagnitudeNo',magno,'PhaseNo',phaseno, ...
    'PhaseX',[],'PhaseY',[],'PhaseCorr',[],'Angio',[],'VelMag',[]);
  linkies = sort([magno phaseno]);
  SET(magno).Flow = flowstruct;
  SET(magno).ImageType = 'Flow (magnitude)';
  SET(magno).Children = phaseno;
  SET(magno).Linked = linkies;
  SET(phaseno).Flow = flowstruct;
  SET(phaseno).ImageType = 'Flow (through plane)';
  SET(phaseno).Parent = magno;
  SET(phaseno).Linked = linkies;
  
  no = phaseno;
  im = calcfunctions('calctruedata',SET(no).IM,no);
  
  %Normalize phase image
  
  if isequal(SET(no).Scanner,'Siemens')
    if isempty(SET(magno).VENC) || SET(magno).VENC == 0
      SET(magno).VENC = SET(no).VENC;
    end
    %venc = SET(magno).VENC;
    %SET(no).IntensityScaling = 2*venc;
    %SET(no).IntensityOffset = -venc;
    if max(im(:))>3000
      if min(im(:))<-3000
        %Newer Siemens files!?
        im = im/single(8192)+single(0.5);
        donenorm = true;
      else
        %Old Siemens files does not specify correct rescale and intercept in
        %DICOM files.
        im = im/4096;
        donenorm = true;
      end
    end
  end
  
  SET(phaseno).IM = im;
  
  magnostri = [magnostri sprintf('%d,',magno)];
  phasenostri = [phasenostri sprintf('%d,',phaseno)];
end

drawfunctions('drawthumbnailframes');
mymsgbox(dprintf('Coupled magnitude stacks %s with phase stacks %s', ...
  magnostri(1:end-1),phasenostri(1:end-1)));

%------------------------
function couplemagtophasestack
%------------------------
%Couple matching magnitude and phase stacks of each loaded flow study.
global SET NO

prompt = {sprintf('%s: ',dprintf('Magnitude stack')),sprintf('%s: ',dprintf('Phase stack'))};
dlg_title = dprintf('Enter magnitude and phase stacks to couple (2D-PC only):');
num_lines = 1;
def = {num2str(NO),num2str(NO+1)};
answer = myinputdlg(prompt,dlg_title,num_lines,def);
nopair = cellfun(@str2num, answer);

  magno = nopair(1);
  phaseno = nopair(2);
  
  %Update
  flowstruct = struct('MagnitudeNo',magno,'PhaseNo',phaseno, ...
    'PhaseX',[],'PhaseY',[],'PhaseCorr',[],'Angio',[],'VelMag',[]);
  linkies = sort([magno phaseno]);
  SET(magno).Flow = flowstruct;
  SET(magno).ImageType = 'Flow (magnitude)';
  SET(magno).Children = phaseno;
  SET(magno).Linked = linkies;
  SET(phaseno).Flow = flowstruct;
  SET(phaseno).ImageType = 'Flow (through plane)';
  SET(phaseno).Parent = magno;
  SET(phaseno).Linked = linkies;
  
  no = phaseno;
  im = calcfunctions('calctruedata',SET(no).IM,no);
  
  %Normalize phase image
  
  if isequal(SET(no).Scanner,'Siemens')
    if isempty(SET(magno).VENC) || SET(magno).VENC == 0
      SET(magno).VENC = SET(no).VENC;
    end
    %venc = SET(magno).VENC;
    %SET(no).IntensityScaling = 2*venc;
    %SET(no).IntensityOffset = -venc;
    if max(im(:))>3000
      if min(im(:))<-3000
        %Newer Siemens files!?
        im = im/single(8192)+single(0.5);
        donenorm = true;
      else
        %Old Siemens files does not specify correct rescale and intercept in
        %DICOM files.
        im = im/4096;
        donenorm = true;
      end
    end
  end
  
  SET(phaseno).IM = im;


drawfunctions('drawthumbnailframes');
mymsgbox(dprintf('Coupled magnitude stacks %d with phase stacks %d', ...
  magno,phaseno));