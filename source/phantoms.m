function varargout = phantoms(varargin)
%PHANTOMS helper function to Segment to generate digital test phantoms.

%Einar Heiberg

%The mighty main clause goes...
macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
  
%--------------------------
function filephantomhelper
%--------------------------
%Helper fcn used when creating image phantoms. Set fake patient data.
global DATA SET NO

DATA.Preview.PatientInfo.Name = 'Image phantom.';
DATA.Preview.PatientInfo.ID = '';
DATA.Preview.PatientInfo.BirthDate = '';
DATA.Preview.PatientInfo.Sex = '';
DATA.Preview.PatientInfo.Age = '';
DATA.Preview.HeartRate = 60;
DATA.Preview.PatientInfo.AcquisitionDate = datestr(now,'yyyy-mm-dd HH:MM:SS');
DATA.Preview.AcquisitionTime = '00.00';
DATA.Preview.PatientInfo.BSA = 0;
DATA.Preview.PatientInfo.Weight = 0;
DATA.Preview.PatientInfo.Length = 0;
DATA.Preview.ResolutionX = 1;
DATA.Preview.ResolutionY = 1;
DATA.Preview.SliceThickness = 5;
DATA.Preview.SliceGap = 0;
DATA.Preview.XSize = size(SET(NO).IM,1);
DATA.Preview.YSize = size(SET(NO).IM,2);
DATA.Preview.TSize = size(SET(NO).IM,3);
DATA.Preview.ZSize = size(SET(NO).IM,4);
DATA.Preview.TIncr = 1;
DATA.Preview.TimeVector = 0:1:size(SET(NO).IM,3)-1;
DATA.Preview.Bitstored = 12;
DATA.Preview.TriggerTime = DATA.Preview.TimeVector;

%----------------------------
function filephantom_Callback %#ok<DEFNU>
%----------------------------
%Displays menu of different built in computer models.
global DATA

c = mymenu('Choose What Phantom to Load/Create.',{...
  'Phantom1',...
  'Cylinder',...
  'Grayscale Test',...
  'Size Test',...
  'Plain Cylinder',...
  'Torus Test (Rotated Image Stack)',...
  '3D Flow Test',...
  'QFlow Test',...
  'PWV Test',...
  'Strain Test 1', ...
  'Strain Test 2', ...
  'T2* phantom',...
  'Patchinessphantom',...
  'Tagging phantom',...
  'PET phantom',...
  'Abort Operation'},DATA.GUI.Segment);

switch c
  case 1
    phantom1;
  case 2
    cylinder;
  case 3
    phantom3;
  case 4
    phantom4;
  case 5
    plaincylinder;
  case 6
    torus;
  case 7
    flow3test;
  case 8
    qflowtest;
  case 9
    pwvtest;
  case 10
    straintest1;
  case 11
    straintest2;
  case 12
    t2starphantom;
  case 13
    patchinessphantom;
  case 14
    taggingphantom;
  case 15
    petphantom;
  otherwise
    mywarning('Aborted by user.',DATA.GUI.Segment);
end;

%----------------
function phantom1
%----------------
global SET NO

%phantom1
nx = 64;
ny = 64;
nz = 12;
nt = 16;
[x,y] = ndgrid(1:nx,1:nx);
NO = length(SET)+1;
SET(NO).IM = repmat(single(0),[nx ny nt nz]);
for loop=1:size(SET(NO).IM,3)
  SET(NO).IM(:,:,loop,:) = repmat(single(sqrt((x-nx/2).^2+(y-ny/2).^2)<4*sqrt(loop)),[1 1 1 nz]);
end;

SET(NO).IM(:,34,:,:) = 1-SET(NO).IM(:,34,:,:);
SET(NO).CenterX = nx/2;
SET(NO).CenterY = nx/2;

filephantomhelper;

%Intialize all variables
openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%-----------------------------
function cylinder(inrad,outrad)
%-----------------------------
%Cylinder phantom.

global DATA SET NO

if nargin==0
  %ask inner radius
  s = inputdlg({'Enter inner radius (endocard)'},'InnerRadius',1,{sprintf('%d',30)});
  if isempty(s)
    myfailed('Invalid inner radius.',DATA.GUI.Segment);
    return;
  else
    [inrad,ok] = str2num(s{1}); %#ok<ST2NM>
    if not(ok)
      myfailed('Invalid inner radius.',DATA.GUI.Segment);
      return;
    end;
  end;
  %ask outer radius
  s = inputdlg({'Enter outer radius (epicard)'},'OuterRadius',1,{sprintf('%d',50)});
  if isempty(s)
    myfailed('Invalid outer radius.',DATA.GUI.Segment);
    return;
  else
    [outrad,ok] = str2num(s{1}); %#ok<ST2NM>
    if not(ok)
      myfailed('Invalid outer radius.',DATA.GUI.Segment);
      return;
    end;
  end;
end;

nx = 128;
ny = 128;
nz = 20;
nt = 1;
[x,y] = ndgrid(1:nx,1:ny);
NO = length(SET)+1;
temp = 0.3*single(sqrt((x-nx*0.5).^2+(y-ny*0.5).^2)<=outrad);
temp = temp+0.7*single(sqrt((x-nx*0.5).^2+(y-ny*0.4).^2)<=inrad);
SET(NO).IM = repmat(single(temp),[1 1 nt nz]);

filephantomhelper;

%Initalize all variables

openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%----------------
function phantom3
%----------------
global SET NO

nx = 256;
x = ndgrid(1:nx,1:nx);
NO = length(SET)+1;
SET(NO).IM = single(x)/nx;
filephantomhelper;

openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%----------------
function phantom4
%----------------
global DATA SET NO

nx = 25;
ny = 100;
NO = length(SET)+1;
im = single(zeros(nx,ny));
im = repmat(im,[1 1 1 10]);
for loop=1:10;
  im((loop):(loop+10),40:60,1,loop) = 1;
end;
SET(NO).IM = single(im);
filephantomhelper;
DATA.Preview.ResolutionX = 1;
DATA.Preview.ResolutionY = 0.25;

openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%----------------------------
function im = helpfilter(im,n)
%----------------------------
xi = linspace(-2,2,n);
f = exp(-xi.^2);
f = f./sum(f(:));

im = conv2(im,f,'same');
im = conv2(im,f','same');

%---------------------------------
function im = helppattern(d,nx,ny)
%---------------------------------

%Create patch image
[x,y] = ndgrid(linspace(-d/2,d/2,d),linspace(-d/2,d/2,d));
rad = sqrt(x.^2+y.^2);
pim = (rad>=(d/2));

%Create im by replicating it
im = repmat(pim,[ceil(nx/d)+1 ceil(ny/d)+1]);

%Cut it
im = single(im(1:nx,1:ny));

%smooth
im = helpfilter(im,5);

%-------------------------
function patchinessphantom
%-------------------------
%Phantom to check patchiness index

global DATA SET NO

nx = 101;
ny = 101;
nz = 16;

endorad = 35;
epirad = 45;

NO = length(SET)+1;

%initalize image
im = zeros(nx,ny,1,nz);

%Create mask
[xi,yi] = ndgrid((1:nx),1:ny);
xi = xi-nx/2;
yi = yi-ny/2;
rad = sqrt(xi.^2+yi.^2);
mask = (rad<=epirad) & (rad>endorad);

%--- Fix the image
tempim = ones(nx,ny);
tempim(1:round(nx/2),1:round(ny/2)) = 1;

im(:,:,1) = helpfilter(tempim.*single(rad<37),5);
im(:,:,2) = helpfilter(tempim.*single(rad<39),5);
im(:,:,3) = helpfilter(tempim.*single(rad<41),5);
im(:,:,4) = helpfilter(tempim.*single(rad<44),5);

%quadrant
im(1:round(nx/2),1:round(ny/2),5) = 1; 
im(:,:,5) = im(:,:,5).*helppattern(8,nx,ny);

%semi-circle
im(1:round(nx/2),:,6) = 1; 
im(:,:,6) = im(:,:,6).*helppattern(8,nx,ny);

%three quadrants
im(:,:,7) = 1;
im(1:round(nx/2),1:round(ny/2),7) = 0; 
im(:,:,7) = im(:,:,7).*helppattern(8,nx,ny);

%Full circle
im(:,:,8) = 1;
im(:,:,8) = im(:,:,8).*helppattern(8,nx,ny);

%pattern
im(:,:,9) = helppattern(8,nx,ny);
im(:,:,10) = helppattern(7,nx,ny);
im(:,:,11) = helppattern(6,nx,ny);
im(:,:,12) = helppattern(5,nx,ny);

im(:,:,13) = helpfilter(rand(nx,ny),10);
im(:,:,14) = helpfilter(rand(nx,ny),5);
im(:,:,15) = helpfilter(rand(nx,ny),2);
im(:,:,16) = rand(nx,ny);

%Apply mask
for loop = 1:size(im,4)
  im(:,:,loop) = im(:,:,loop).*single(mask);
end;

%Store image
SET(NO).IM = single(im);

%Set it up
filephantomhelper;
DATA.Preview.ResolutionX = 1;
DATA.Preview.ResolutionY = 1;

openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);

%Add endo & epicardium
omega = linspace(0,2*pi,DATA.NumPoints);
omega = omega(:);
x = sin(omega);
y = cos(omega);
mx = nx/2;
my = ny/2;

SET(NO).EndoX = repmat(mx+x*endorad,[1 1 nz]);
SET(NO).EndoY = repmat(my+y*endorad,[1 1 nz]);
SET(NO).EpiX = repmat(mx+x*epirad,[1 1 nz]);
SET(NO).EpiY = repmat(my+y*epirad,[1 1 nz]);

viability('viabilityautoweighted');

%Remove automated no reflow detection
SET(NO).Scar.NoReflow = false(size(SET(NO).Scar.NoReflow));

[~,fullindex] = viability('calcpatchinessindex',NO);

SET(NO).Point = [];
SET(NO).Point.X = repmat(nx/2,[1 nz]);
SET(NO).Point.Y = repmat(ny/2,[1 nz])-10;
SET(NO).Point.T = nan(1,nz);
SET(NO).Point.Z = 1:nz;

l = cell(1,nz);
for loop = 1:nz
  l{loop} = sprintf('%0.3g',fullindex(loop));
end;

SET(NO).Point.Label = l;

segment('makeviewim');
drawfunctions('drawall');



%-----------------------------------
function plaincylinder(inrad,outrad)
%-----------------------------------
global DATA SET NO

%plain cylinder

if nargin==0
  %ask inner radius
  s = inputdlg({'Enter inner radius (endocard)'},'InnerRadius',1,{sprintf('%d',20)});
  if isempty(s)
    myfailed('Invalid inner radius.',DATA.GUI.Segment);
    return;
  else
    [inrad,ok] = str2num(s{1}); %#ok<ST2NM>
    if not(ok)
      myfailed('Invalid inner radius.',DATA.GUI.Segment);
      return;
    end;
  end;
  %ask outer radius
  s = inputdlg({'Enter outer radius (epicard)'},'OuterRadius',1,{sprintf('%d',30)});
  if isempty(s)
    myfailed('Invalid outer radius.',DATA.GUI.Segment);
    return;
  else
    [outrad,ok] = str2num(s{1}); %#ok<ST2NM>
    if not(ok)
      myfailed('Invalid outer radius.',DATA.GUI.Segment);
      return;
    end;
  end;
end;

nx = 100;
ny = 100;
nz = 3;
nt = 3;
[x,y] = ndgrid(1:nx,1:ny);
NO = length(SET)+1;
temp = 0.3*single(sqrt((x-nx*0.5).^2+(y-ny*0.5).^2)<=outrad);
temp = temp+0.7*single(sqrt((x-nx*0.5).^2+(y-ny*0.5).^2)<=inrad);
SET(NO).IM = repmat(single(temp),[1 1 nt nz]);
filephantomhelper;

openfile('setupstacksfromdicom',NO);

SET(NO).TIncr = 1;
segment('renderstacksfromdicom',NO);

%-------------
function torus
%-------------
global DATA SET NO

%torus test
inrad = 10;
outrad = 15;
nx = 100;
ny = 100;
nz = 3;
nt = 3;
[x,y] = ndgrid(1:nx,1:ny);
NO = length(SET)+1;
temp = 0.3*single(sqrt((x-nx*0.5).^2+(y-ny*0.8).^2)<=outrad);
temp = temp+0.7*single(sqrt((x-nx*0.5).^2+(y-ny*0.8).^2)<=inrad);
SET(NO).IM = repmat(single(temp),[1 1 nt nz]);
temp = 0.3*single(sqrt((x-nx*0.5).^2+(y-ny*0.2).^2)<=outrad);
temp = temp+0.7*single(sqrt((x-nx*0.5).^2+(y-ny*0.2).^2)<=inrad);
SET(NO).IM = SET(NO).IM+repmat(single(temp),[1 1 nt nz]);
filephantomhelper;

openfile('setupstacksfromdicom',NO);

SET(NO).TIncr = 0.3;
segment('renderstacksfromdicom',NO);

c = 30/10; %cm radius to tube center center
alpha = 10/10; %cm radius of tube
volume = 2*(pi^2)*(alpha^2)*c;
mymsgbox(...
  sprintf('Innerradius=10, Outerradius=15, Radius of torso (center)=30, volume=%0.5g ml',...
  volume),'',DATA.GUI.Segment);

%-----------------
function qflowtest
%-----------------
global SET NO

%Qflow test
nx = 64;
ny = 64;
nz = 1;
nt = 16;
rad = 20;
statsz = 10;
[x,y] = ndgrid(1:nx,1:ny);

%Magnitude image
mag = single(sqrt((x-nx*0.5).^2+(y-ny*0.5).^2)<=rad);
mag = repmat(mag,[1 1 nt nz]);
mag(1:statsz,1:statsz,:,:) = 1;
mag(1:statsz,(ny-statsz):end,:,:) = 1;
mag((nx-statsz):end,1:statsz,:,:) = 1;
mag((nx-statsz):end,(ny-statsz):end,:,:) = 1;
NO = length(SET)+1; magno = NO;
SET(NO).IM = mag+0.05*randn(size(mag));
SET(NO).VENC = 200;
filephantomhelper;

openfile('setupstacksfromdicom',NO);

%Phase image
t = 0.5*cos(linspace(0,2*pi,nt));
%t = repmat(0.5,[1 nt]);
phase = mag;
for tloop=1:nt
  phase(:,:,tloop,:) = (mag(:,:,tloop,:)>0.3)*t(tloop)+0.5;
end;
log = (mag<0.3);
phase(log) = 1;

%phase(log) = phase(log).*rand(size(log));
phase(1:statsz,1:statsz,:,:) = 0.5;
phase(1:statsz,(ny-statsz):end,:,:) = 0.5;
phase((nx-statsz):end,1:statsz,:,:) = 0.5;
phase((nx-statsz):end,(ny-statsz):end,:,:) = 0.5;
NO = length(SET)+1; phaseno = NO;
SET(NO).IM = phase; %+0.05*randn(size(phase));
filephantomhelper;

openfile('setupstacksfromdicom',NO);

%Set phase and magnitude info
SET(magno).Flow.PhaseNo = phaseno;
SET(magno).Flow.PhaseX = [];
SET(magno).Flow.PhaseY = [];
SET(magno).Flow.Angio = [];
SET(magno).Flow.VelMag = [];

SET(magno).Flow.MagnitudeNo = magno;
SET(phaseno).Flow.PhaseNo = phaseno;
SET(phaseno).Flow.MagnitudeNo = magno;
SET(phaseno).Flow.PhaseX = [];
SET(phaseno).Flow.PhaseY = [];
SET(phaseno).Flow.Angio = [];
SET(phaseno).Flow.VelMag = [];
SET(magno).VENC = 200;
SET(magno).ImageType = 'Flow (Magnitude)';

SET(phaseno).VENC = 200;
SET(phaseno).ImageType = 'Flow (Through plane)';

segment('renderstacksfromdicom',NO);

%------------------------------
function yi = upsamplehelper(y,n)
%------------------------------
%Helper functionk to upsample signal to n points

x = linspace(1,length(y),length(y));
xi = linspace(1,length(y),n);
yi = interp1(x,y,xi,'pchip');

%---------------------------------
function [s,deltat] = impulse(l,t)
%---------------------------------
%Reponse function with length l samples and delay t in ms fraction of signal
%(0.1 = 10%))

n = 10000;
%Generate high resolution signal
s = zeros(1,n);

indend = round(t*n/1000);
s(indend) = 1;

%Add capacitance
x = linspace(0,1000,n);
f = exp(-x/(t/5)); %down to almost e-1 after t/5

%Convolve capacitance
s = conv(s,f);
s = s(1:n);

%Resample to desired length
s = upsamplehelper(s,l);

%normalize
s = 0.6*s./sum(s(:));

%Center of gravity
deltat = sum(x.*s)/sum(s);

%----------------------------------
function pwhelper(rad,flow,roiname)
%----------------------------------
%Helper function to add image stacks

global SET NO

nx = 256;
ny = 256;
nz = 1;
nt = length(flow);
venc = 150;

[x,y] = ndgrid(1:nx,1:ny);

%Magnitude image
mag = single(sqrt((x-nx*0.5).^2+(y-ny*0.5).^2)<=rad);
mag = repmat(mag,[1 1 nt nz]);
NO = length(SET)+1; magno = NO;
SET(NO).IM = mag; %+0.05*randn(size(mag));
SET(NO).VENC = venc;
filephantomhelper;

openfile('setupstacksfromdicom',NO);

%Convert flow to pixel values
numpixels = sum(sum(mag(:,:,1,1)));
area = 0.1*0.1*numpixels; %cm^2

%vol = area*vel*deltat => vel = vol/(area*deltat)
velocity = flow/(area*1); %cm/s

%Convert to phasevalue 0..1 (0.5=0 cm/s);
phasevalue = 0.5+0.5*velocity/venc;

if max(phasevalue)>1
  error('Wrapping.');
end;

if max(phasevalue)<0
  error('Wrapping.');
end;

%Phase image
phase = mag;
for tloop=1:nt
  phase(:,:,tloop,:) = phasevalue(tloop)*mag(:,:,tloop,:);
end;

phase(phase==0) = 0.5;
phase(1) = 0; %Ensure normalization
phase(2) = 1;

%Store
NO = length(SET)+1; phaseno = NO;
SET(NO).IM = phase; %+0.05*randn(size(phase));
filephantomhelper;

openfile('setupstacksfromdicom',NO);

%Set phase and magnitude info
SET(magno).Flow.MagnitudeNo = magno;
SET(magno).Flow.PhaseNo = phaseno;
SET(magno).Flow.PhaseX = [];
SET(magno).Flow.PhaseY = [];
SET(magno).Flow.Angio = [];
SET(magno).Flow.VelMag = [];

SET(phaseno).Flow.PhaseNo = phaseno;
SET(phaseno).Flow.MagnitudeNo = magno;
SET(phaseno).Flow.PhaseX = [];
SET(phaseno).Flow.PhaseY = [];
SET(phaseno).Flow.Angio = [];
SET(phaseno).Flow.VelMag = [];
SET(magno).VENC = venc;
SET(magno).ImageType = 'Flow (Magnitude)';

SET(phaseno).VENC = venc;
SET(phaseno).ImageType = 'Flow (Through plane)';

%Add ROI
SET(magno).Roi = [];
omega = linspace(0,2*pi,80)';
SET(magno).Roi.X = repmat(rad*sin(omega),[1 nt])+nx/2;
SET(magno).Roi.Y = repmat(rad*cos(omega),[1 nt])+ny/2;
SET(magno).Roi.T = 1:nt;
SET(magno).Roi.Z = 1;
SET(magno).Roi.Sign = 1;
SET(magno).Roi.Name = roiname;
SET(magno).Roi.LineSpec = 'b-';
SET(magno).Roi.InterpX = [];
SET(magno).Roi.InterpY = [];
SET(magno).RoiN = 1;

segment('renderstacksfromdicom',NO);

%Store again to avoid normalization
SET(phaseno).IM = phase;

%Change timevector
SET(magno).TimeVector = linspace(0,1,nt);
SET(magno).TIncr = SET(magno).TimeVector(2)-SET(magno).TimeVector(1);

SET(phaseno).TimeVector = linspace(0,1,nt);
SET(phaseno).TIncr = SET(phaseno).TimeVector(2)-SET(phaseno).TimeVector(1);


%---------------
function pwvtest
%---------------
%Create phantoms for different resolutions

global SET

aorticlength = 240;

points = [20:40]; % 25 30 35 40];

medianerrorpercent = nan(1,length(points));
medianerrorvel = nan(1,length(points));
errorpercent20 = nan(1,length(points));
errorpercent14 = nan(1,length(points));
errorpercent10 = nan(1,length(points));
errorpercent4 = nan(1,length(points));
maxpwv = 20;

for loop = 1:length(points);
  
  ti = aorticlength./linspace(2,maxpwv,10);
  
  ttrue = zeros(1,length(ti));
  tmeas = zeros(1,length(ti));  
  
  for tloop = 1:length(ti)
    
    aorticflow = [2.598908 35.102913 143.620758 269.466339 345.688629 377.106323 ....
      370.276184 345.908478 315.445404 281.697540 238.217163 191.375366 138.121368 ...
      71.464691 -11.375866 -46.683861 -28.849838 -15.191250 2.068810 19.120804 ...
      26.928499 21.993345 13.702498 11.135820 7.779805 6.833653 6.753801 3.557485 ...
      2.488269 0.022708 -1.034425 0.541035 -0.553828 1.434078 4.693777 4.872540 ...
      5.345801 2.741043 -1.921593 -5.615703];
    
    aorticflow2 = [-12.269211 1.010513 104.714584 232.957077 334.462646 390.486237 405.852081 ...
      393.985596 361.134277 311.153259 252.534485 191.633606 117.670441 27.972794 -7.453537 ...
      3.009796 8.271790 21.235657 27.738190 19.450378 8.251190 1.799011 -1.717758 -2.268219 ...
      -0.423431 -0.128174 0.008011 0.706100 0.537872 -0.161362 0.516129 1.579285 0.289536 ...
      0.018311 1.115799 2.776337 4.616547 5.564117 1.755524 -5.348969];    
    
    %Sample up
    n = 10*1000;
    aorticflow = upsamplehelper(aorticflow,n);
    
    %Generate impulse response
    %[imp,deltat] = impulse(n,ti(tloop));   
    imp = zeros(size(aorticflow));
    imp(round(ti(tloop)*10)) = 1; %*10 as
    
    output = conv(aorticflow,imp);
    output = output(1:length(aorticflow));
    
    %figure(23+loop);
    %plot(imp);
    
    ttrue(tloop) = ti(tloop); %was deltat
    
    %Sample down
    n = points(loop);
    
    aorticflow = upsamplehelper(aorticflow,n);
    abdominalflow = upsamplehelper(output,n);
    
    %Sample up again
    n = 200;
    aorticflow = upsamplehelper(aorticflow,n);
    abdominalflow = upsamplehelper(abdominalflow,n);
    
    if 1==0
      figure(78);
      clf;
      subplot(2,1,1);
      plot(aorticflow);
      hold on;
      plot(abdominalflow);
      hold off;
      subplot(2,1,2);
      plot(imp);
      pause;
    end;
    
    tmeas(tloop) = pwtestcalchelper(aorticflow,abdominalflow);
  end;
  
  truevel = aorticlength./ttrue;
  measvel = aorticlength./tmeas;
  
  medianerrorpercent(loop) = median(100*abs(truevel-measvel)./truevel);  
  medianerrorvel(loop) = median(truevel-measvel);
  errorpercent20(loop) = 100*(truevel(end)-measvel(end))./truevel(end);  
  errorpercent14(loop) = 100*(truevel(7)-measvel(7))./truevel(7);
  errorpercent10(loop) = 100*(truevel(5)-measvel(5))./truevel(5);
  errorpercent4(loop) = 100*(truevel(2)-measvel(2))./truevel(2);
  
  fs = 15;
    
  if loop==16
    1
  end;
  
  disp(sprintf('%d timeframes: %0.5g+-%0.5g',points(loop),mean(truevel-measvel),std(truevel-measvel)));
  
  if 1==0
  figure(69+loop);
  %subplot(4,3,loop);
  h = plot(truevel,measvel,'ko');  
  set(h,'linewidth',2);
  set(h,'markersize',8);
  hold on;
  h = plot([0 maxpwv],[0 maxpwv],'k-'); set(h,'linewidth',2);
  hold off;
  axis image 
  set(gca,'xlim',[0 maxpwv],'ylim',[0 maxpwv],'fontsize',fs);
  h = xlabel('True pulse wave velocity [m/s]'); set(h,'fontsize',fs);
  h = ylabel('Measured pulse wave velocity [m/s]'); set(h,'fontsize',fs);
  h = title(dprintf('Derived from %d timeframes',points(loop))); set(h,'fontsize',fs);  
  set(69+loop,'color',[1 1 1]);
  end;
  
end;

figure(190);
h = plot(points,errorpercent20,'r-');set(h,'linewidth',2); set(h,'markersize',6);
hold on;
h = plot(points,errorpercent14,'g-');set(h,'linewidth',2); set(h,'markersize',6);
h = plot(points,errorpercent10,'m-.');set(h,'linewidth',2); set(h,'markersize',6);
h = plot(points,errorpercent4,'b-.');set(h,'linewidth',2); set(h,'markersize',6);
plot([20 40],[-6 -6],'k:');
plot([20 40],[6 6],'k:');
plot([35 35],[-60 60],'k-');
hold off;
legend('PWV 20 m/s','PWV 14 m/s','PWV 10 m/s','PWV 4 m/s');
h = xlabel('Time points'); set(h,'fontsize',fs);
h = ylabel('PWV error [%]'); set(h,'fontsize',fs);
set(gca,'ylim',[-60 60],'fontsize',fs);
set(190,'color',[1 1 1]);

outdata = cell(21,4);
for loop=20:40
  outdata{loop-19,1} = loop;
  outdata{loop-19,2} = 1000/loop;
  outdata{loop-19,3} = medianerrorvel(loop-19);
  outdata{loop-19,4} = medianerrorpercent(loop-19);
  outdata{loop-19,6} = errorpercent20(loop-19);
  outdata{loop-19,7} = errorpercent14(loop-19);
  outdata{loop-19,8} = errorpercent10(loop-19);
  outdata{loop-19,9} = errorpercent4(loop-19);
end;

segment('cell2clipboard',outdata);

return;

aortaarea = 4.6; %cm2

%Calculate radius
rad = sqrt(aortaarea/pi); %A=pi*r^2
rad = rad*10; %convert from cm to mm

%Set up image stacks
pwhelper(rad,aorticflow,'Aortic ascending flow')

%Set up image stacks
pwhelper(rad,abdominalflow,'Abdominal Aorta')

SET(1).Measure.X = [1 200];
SET(1).Measure.Y = [50 50 ];
SET(1).Measure.Z = 1;
SET(1).Measure.Length = 200;
SET(1).Measure.Name = 'AL';
SET(1).Measure.LongName = 'Aortic Length';
SET(1).Measure.T = 1;
     
%----------------------------------------------
function deltat = pwtestcalchelper(flow1,flow2)
%----------------------------------------------
%Excerpt from code in PWV

tmax = 0;
ymin = inf;
ymax = -inf;
slopet = [0 0];

for i = 1:2
  if i==1
    flowcurve = flow1;
  else
    flowcurve = flow2;
  end;
  
  timevec = linspace(0,1000,length(flow1));
  tmax = max(tmax,timevec(end));
  ymin = min(ymin,min(flowcurve));
  ymax = max(ymax,max(flowcurve));
  
  sigma = 0.025;
  
  x = linspace(-1,1,length(timevec));
  f = exp(-x.^2/sigma.^2);
  smoothcurve = conv(flowcurve,f,'same')./conv(ones(size(flowcurve)),f,'same');
  %smoothcurve = flowcurve;
  
  [maxslope,maxix] = max(diff(smoothcurve));
  k = maxslope/mean(diff(timevec));
  maxtime = timevec(maxix);
  m = (flowcurve(maxix)-k*maxtime);
  slopet(i) = -m/k;
  
  
end

deltat = diff(slopet);

%-----------------
function flow3test
%-----------------
global SET NO

%3D flow test
nx = 64;
ny = 64;
nz = 4;
nt = 16;
rad = 20;
statsz = 10;
[x,y] = ndgrid(1:nx,1:ny);

%Magnitude image
mag = single(sqrt((x-nx*0.5).^2+(y-ny*0.5).^2)<=rad);
mag = repmat(mag,[1 1 nt nz]);
mag(1:statsz,1:statsz,:,:) = 1;
mag(1:statsz,(ny-statsz):end,:,:) = 1;
mag((nx-statsz):end,1:statsz,:,:) = 1;
mag((nx-statsz):end,(ny-statsz):end,:,:) = 1;
NO = length(SET)+1; magno = NO;
SET(NO).IM = mag+0.05*randn(size(mag));
SET(NO).VENC = 200;
filephantomhelper;

openfile('setupstacksfromdicom',NO);

%Phase image
t = 0.5*cos(linspace(0,2*pi,nt));
phase = mag;
for tloop=1:nt
  phase(:,:,tloop,:) = (mag(:,:,tloop,:)>0.3)*t(tloop)+0.5;
end;
log = find((mag<0.3));
phase(log) = 1;
phase(log) = phase(log).*rand(size(log));
phase(1:statsz,1:statsz,:,:) = 0.5;
phase(1:statsz,(ny-statsz):end,:,:) = 0.5;
phase((nx-statsz):end,1:statsz,:,:) = 0.5;
phase((nx-statsz):end,(ny-statsz):end,:,:) = 0.5;
NO = length(SET)+1; phaseno = NO;
SET(NO).IM = phase+0.05*randn(size(phase));
filephantomhelper;

openfile('setupstacksfromdicom',NO);

%Set phase and magnitude info
SET(magno).Flow.PhaseNo = phaseno;
SET(magno).Flow.PhaseX = [];
SET(magno).Flow.PhaseY = [];
SET(magno).Flow.Angio = [];
SET(magno).Flow.VelMag = [];
SET(magno).Flow.Result = [];
SET(magno).Flow.MagnitudeNo = magno;
SET(phaseno).Flow.PhaseNo = phaseno;
SET(phaseno).Flow.MagnitudeNo = magno;
SET(phaseno).Flow.PhaseX = [];
SET(phaseno).Flow.PhaseY = [];
SET(phaseno).Flow.Angio = [];
SET(phaseno).Flow.VelMag = [];
SET(phaseno).Flow.Result = [];

SET(magno).VENC = 200;
SET(magno).ImageType = 'Flow (Magnitude)';

SET(phaseno).VENC = 200;
SET(phaseno).ImageType = 'Flow (Through plane)';

segment('renderstacksfromdicom',NO);

%-------------------
function straintest1
%-------------------
global DATA SET NO

%Strain test
SET = [];
nx = 64;
ny = 64;
nz = 1;
nt = 20;
venc = 0.1; %cm/s

%Setup preview structure
DATA.Preview.ResolutionX = 1; %mm
DATA.Preview.ResolutionY = 1;
DATA.Preview.TIncr = 50/1000; %50ms

%Magnitude image
SET(1).IM = repmat(single(0.5),[nx ny nt nz]);
NO = 1;
magno = NO;
filephantomhelper;
openfile('setupstacksfromdicom',NO);

%Flow image x
SET(2).IM = SET(1).IM;
v = 0.5-0.5*single(sin(2*pi*(0:(nt-1))/nt));
for loop=1:nt
  SET(2).IM(:,:,loop,1) = v(loop);
end;
phasenox = 2;
NO = 2;     
filephantomhelper;
openfile('setupstacksfromdicom',NO);

%Flow image y
SET(3).IM = single(0.5+0.1*rand(size(SET(1).IM)));
SET(3).IM(1) = 0;
SET(3).IM(2) = 1;
phasenoy = 3;
NO = 3;
filephantomhelper;
openfile('setupstacksfromdicom',NO);

%---  Set flow data structure
SET(magno).Flow.PhaseNo = [];
SET(magno).Flow.PhaseX = phasenox;
SET(magno).Flow.PhaseY = phasenoy;
SET(magno).Flow.Angio = [];
SET(magno).Flow.VelMag = [];
SET(magno).Flow.MagnitudeNo = magno;
SET(magno).Flow.Result = [];

SET(magno).VENC = venc; %cm/s
SET(phasenox).Flow.PhaseNo = [];
SET(phasenox).Flow.MagnitudeNo = magno;
SET(phasenox).Flow.PhaseX = phasenox;
SET(phasenox).Flow.PhaseY = phasenoy;
SET(phasenox).Flow.Angio = [];
SET(phasenox).Flow.VelMag = [];

SET(phasenox).Flow.Result = [];
SET(phasenox).VENC = venc; %cm/s
SET(phasenoy).Flow.PhaseNo = [];
SET(phasenoy).Flow.MagnitudeNo = magno;
SET(phasenoy).Flow.PhaseX = phasenox;
SET(phasenoy).Flow.PhaseY = phasenoy;
SET(phasenoy).Flow.Angio = [];
SET(phasenoy).Flow.VelMag = [];
SET(phasenoy).Flow.Result = [];
SET(phasenoy).VENC = venc; %cm/s
SET(magno).ImageType = 'Flow (Magnitude)';
SET(phasenox).ImageType = 'Flow (up/down)';
SET(phasenoy).ImageType = 'Flow (left/right)';

SET(1).ImageType = 'Strain FFE';
SET(1).ImageViewPlane = '2CH';

segment('renderstacksfromdicom',NO);

%-------------------
function straintest2
%-------------------
global DATA SET NO

%Strain test
SET = [];
nx = 64;
ny = 64;
nz = 1;
nt = 20;
venc = 20; %cm/s
vscale = 2; %Amplitude of the velocity
noisescale = 1/100;

%Setup preview structure
DATA.Preview.ResolutionX = 1; %Pixel resolution mm
DATA.Preview.ResolutionY = 1; %mm
DATA.Preview.TIncr = 1/nt;

%Generate timevector
timevec = linspace(0,2*pi,nt); %Generate "time"...

%Generate velocity vector
vel = vscale*sin(timevec); %in cm/s
pos = cumsum(vel)*DATA.Preview.TIncr*10; %position of "piston" in mm, *10 due to get it into mm.

%----------------------
% for veloop=1:nt
%   if (veloop == 1) 
%     pos(veloop) = 0; % cm
%   else 
%     accel = (vel(veloop) - vel(veloop-1))/DATA.Preview.TIncr;
%     disp(veloop) = vel(veloop-1)*DATA.Preview.TIncr + accel*(DATA.Preview.TIncr^2)/2;
%     pos(veloop) = pos(veloop-1)+ disp(veloop); %cm
%   end  
% end
% pos = pos * 10; % mm
%----------------------------

%figure(12);
%plot(1:nt,vel);
%title('velocity [cm/s]');

%figure(13);
%plot(1:nt,pos);
%title('position mm');

%Magnitude image, a bit uggly code since it assumes resolution to be 1mm.
SET(1).IM = repmat(single(0),[nx ny nt nz]);
for tloop=1:nt
  SET(1).IM(16:round(48-pos(tloop)),16:48,tloop) = 1;
end;

NO = 1;
magno = NO;
filephantomhelper;

openfile('setupstacksfromdicom',NO);
%----------------------------------------------
%Flow image x
SET(2).IM = 0*SET(1).IM; %This is just to create space for it.
SET(2).IM = SET(2).IM+0.5; %The normal "zero" value is 0.5
SET(2).IM(1) = 1;
SET(2).IM(2) = 0;

%Loop over time.
for tloop=1:nt
  sweep = 0.5+zeros(nx,1); %nx is size of image
  templength = length(16:round(48-pos(tloop))); %Take out only the portion of where the box is...
  
  %vel(t) is velocity in cm/s over time.
  
  %the next line is the critical line. -minus to get correction direction
  sweep(16:round(48-pos(tloop))) = vel2segment(linspace(0,-vel(tloop),templength),venc); %Generate velocity field....

  
  sweep = sweep(:); %Reshape it.
  SET(2).IM(:,16:48,tloop) = repmat(sweep,[1 33]); %Store into the matrix.
end;

max(SET(2).IM(:))
min(SET(2).IM(:))

SET(2).IM = SET(2).IM+noisescale*randn(size(SET(2).IM));
phasenox = 2;
NO = 2;     
filephantomhelper;
openfile('setupstacksfromdicom',NO);

%-------------------------------------------------------
%Flow image y
%Simple since it should be constant.
SET(3).IM = 0.5+repmat(single(0),size(SET(1).IM));
SET(3).IM(1) = 0; %only add this for normalization issue, up in the left 
SET(3).IM(2) = 1; %corner.
phasenoy = 3;
NO = 3;
filephantomhelper;
openfile('setupstacksfromdicom',NO);

%-----------------------------------------------------
%---  Set flow data structure
SET(magno).Flow.PhaseNo = [];
SET(magno).Flow.PhaseX = phasenox;
SET(magno).Flow.PhaseY = phasenoy;
SET(magno).Flow.Angio = [];
SET(magno).Flow.VelMag = [];
SET(magno).Flow.MagnitudeNo = magno;
SET(magno).Flow.Result = [];


SET(magno).VENC = venc; %cm/s
SET(phasenox).Flow.PhaseNo = [];
SET(phasenox).Flow.MagnitudeNo = magno;
SET(phasenox).Flow.PhaseX = phasenox;
SET(phasenox).Flow.PhaseY = phasenoy;
SET(phasenox).Flow.Angio = [];
SET(phasenox).Flow.VelMag = [];
SET(phasenox).Flow.Result = [];

SET(phasenox).VENC = venc; %cm/s
SET(phasenoy).Flow.PhaseNo = [];
SET(phasenoy).Flow.MagnitudeNo = magno;
SET(phasenoy).Flow.PhaseX = phasenox;
SET(phasenoy).Flow.PhaseY = phasenoy;
SET(phasenoy).Flow.Angio = [];
SET(phasenoy).Flow.VelMag = [];
SET(phasenoy).Flow.Result = [];

SET(phasenoy).VENC = venc; %cm/s
SET(magno).ImageType = 'Flow (Magnitude)';
SET(phasenox).ImageType = 'Flow (up/down)';
SET(phasenoy).ImageType = 'Flow (left/right)';

SET(1).ImageType = 'Strain FFE';
SET(1).ImageViewPlane = '2CH';

tincr = 50/1000; %50ms Time increment between frames, expressed in seconds.
SET(magno).TIncr = tincr;
SET(phasenox).TIncr = tincr;
SET(phasenoy).TIncr = tincr;

drawfunctions('drawall');
segment('renderstacksfromdicom',NO);

%-----------------------------
function t2starphantom(noisef)
%-----------------------------
global DATA SET NO

if nargin==0
  noisef = 0.01;
end;

nx = 100;
ny = nx;
nt = 12;
NO = length(SET)+1;
SET(NO).IM = single(zeros(nx,nx,nt)); %reserve memory

SET(NO).ResolutionX = 1;
SET(NO).ResolutionY = 1;

%Set standard fields
filephantomhelper;

%Set special fields
DATA.Preview.TIncr = 1;

%True t2 values
%t2map = ndgrid(linspace(1,nx,nx),ones(1,nx));
t2map = 40*ones(nx,nx); %ms
t2 = linspace(1,50,100);
pos = 1;
for xloop=1:10
  for yloop=1:10
    t2map((1+(xloop-1)*10):(xloop*10),(1+(yloop-1)*10):(yloop*10)) = t2(pos);
    pos = pos+1;    
  end;
end;

r2map = 1./t2map;

%echo times
echotimes = [1 4 8 12 16 20 24 28 32 36 40 44];

%Loop over timeframes
for tloop=1:nt
  for xloop=1:nx
    for yloop=1:ny
      SET(NO).IM(xloop,yloop,tloop) = single(exp(-r2map(xloop,yloop)*echotimes(tloop)));
    end;
  end;
end;

%Add noise
SET(NO).IM = SET(NO).IM+noisef*randn(size(SET(NO).IM));
SET(NO).IM = abs(SET(NO).IM);

openfile('setupstacksfromdicom',NO);
SET(NO).EchoTime = echotimes;

segment('renderstacksfromdicom',NO);

%----------------------
function taggingphantom
%----------------------
%Create phantom for tagging images
global SET NO

%Define limits
x = -1:0.01:1;
y = x;
t = 0:0.05:1;

%Define function for pumping motion
%d = inputdlg({'Radial power (0-1)','Circumferential power (0-1)'});
%pwrs = str2double(d);
pwrs = [0.1 0]; % parameter to change
if pwrs(1) >= 0 && pwrs(1) <= 1 && pwrs(2) >= 0 && pwrs(2) <= 1
  radialpower = pwrs(1);
  circpower = pwrs(2);
else
  myfailed('Invalid input of strain power');
  return
end

%R = diastolic diameter of center of myocardium
basicradpumpfcn = @(R,r,t)(R-r-sin(pi*t)*radialpower);
radpumpfcn = @(r,t)basicradpumpfcn(0.8,r,t);
circpumpfcn = @(r,t)(r+sin(pi*t)*circpower);
spreadpumpfcn = @(t)(10*circpumpfcn(-5,t));
%radialsin = @(r,t)sin(pumpfcn(r,t));
radialexp = @(r,t)exp(-50*radpumpfcn(r,t).^2);
radialexp1 = @(r,t)exp(-300*basicradpumpfcn(0.9,r,t).^2);
radialexp2 = @(r,t)exp(-300*radpumpfcn(r,t).^2);
radialexp4 = @(r,t)exp(-300*basicradpumpfcn(0.7,r,t).^2);
radialexp5 = @(r,t)exp(spreadpumpfcn(t).*radpumpfcn(r,t).^2); %thickening

ang = @(x,y,t)circpumpfcn(atan(y./x),t);
angleexp = @(ang)(exp(-50*(mod(ang,pi/4)-pi/8).^2));
angleexp1 = @(ang)(exp(-400*(mod(ang+pi/32,pi/16)-pi/32).^2)); %32 points, 300 decide size of points
angleexp2 = @(ang)(exp(-400*(mod(ang,pi/32)-pi/64).^2)); %64 points
rr = @(x,y)sqrt(x.^2+y.^2);

%Create phantom image
[X,Y,T] = meshgrid(x,y,t);
Z = radialexp1(rr(X,Y),T).*angleexp1(ang(X,Y,T)) + ...
  radialexp2(rr(X,Y),T).*angleexp2(ang(X,Y,T));
Z = radialexp2(rr(X,Y),0).*angleexp1(ang(X,Y,T)) + ...
  radialexp4(rr(X,Y),T).*angleexp2(ang(X,Y,T));
% Z = radialexp5(rr(X,Y),T).*angleexp1(ang(X,Y,T));

NO = length(SET)+1;
SET(NO).IM = single(Z);

SET(NO).ResolutionX = 1;
SET(NO).ResolutionY = 1;

%Set standard fields
filephantomhelper;

%Setup and render stacks
openfile('setupstacksfromdicom',NO);
SET(NO).TimeVector = t;
SET(NO).TIncr = t(2)-t(1);
SET(NO).IntensityOffset = 0;
SET(NO).IntensityScaling = 255;
SET(NO).ImageType = 'Strain from tagging';
SET(NO).ImageViewPlane = 'Short-axis';
segment('renderstacksfromdicom',NO);

%-------------------------------------
function [res] = vel2segment(vel,venc)
%-------------------------------------
%Converts from velocity in cm/s to venc

res = 0.5+vel/(2*venc);

%------------------
function petphantom
%------------------
% PET phantom
% John Heerfordt Sjöqvist & Mattias Albinsson
global SET NO

%Set image size
nx = 64;
ny = 64;
nt = 15;
nz = 1;

NO = length(SET) + 1;
no = NO;

%Selection of number of compartments
choice = questdlg('Choose model:', ...
	'Compartment model selection', ...
	'1C','2C','1C');
switch choice
    case '1C'
        nbrOfCompartments = 1;
    case '2C'
        nbrOfCompartments = 2;
end

%Declaration of time interval
t = 1:1:240; %[s]

%Start and end times of each frame
frame_timeint = [1 10; 11 20; 21 30; 31 40; 41 50; 51 60; 61 70;...
    71 80; 81 90; 91 100; 101 110; 111 120; 121 150; 151 180; 181 240];

%Midpoints of frames
t_frames = [5:10:115 135 165 210]; %[s]
nbrOfImages = length(t_frames);

%Set parameter values (at rest from litterature, min -> seconds)
K_1 = 0.8/60; %[ml/(g*s)]
k_2 = 0.2/60; %[s^-1]
TBV = 0.4; %[dimensionless]
if (nbrOfCompartments == 2) 
    k_3 = 0.11/60; %[s^-1]
end

%Calculate arterial input function. Associates the midpoint 
%of each frame with the mean of the activity.
C_a = input_function(t);
C_a_mean = zeros(1,nbrOfImages);
for n = 1:nbrOfImages
    C_a_mean(n) = mean(C_a(frame_timeint(n,1):frame_timeint(n,2)));
end

%Addition of Gaussian noise to input
%C_a_mean = C_a_mean + 20000*randn(1,nbrOfImages);

%Integration of o.d.e.
if (nbrOfCompartments == 1)
    dCt = @(t,Ct) K_1*C_a(round(t)) - k_2*Ct;
    [~,C_T] = ode45(dCt,t,0);
    C_T = C_T';
else
    dCe = @(t,Ce) K_1*C_a(round(t)) - (k_2 + k_3)*Ce;
    [~,C_E] = ode45(dCe,t,0);
    C_E = C_E';
    dCg = @(t,Cg) k_3*C_E(round(t));
    [~,C_G] = ode45(dCg,t,0);
    C_G = C_G';
    C_T = C_E + C_G;
end

%Calculate the measured activity
C_m = TBV*C_a + (1-TBV)*C_T;
C_m_mean = zeros(1,nbrOfImages);
for n = 1:nbrOfImages
    C_m_mean(n) = mean(C_m(frame_timeint(n,1):frame_timeint(n,2)));
end

%Set radius of blood pool and annular myocardium
rad_outer = 25;
rad_inner = 15;

%Image formation
[x,y] = ndgrid(1:nx,1:ny);

im = sqrt((x-nx*0.5).^2+(y-ny*0.5).^2) <= rad_outer;
im2 = sqrt((x-nx*0.5).^2+(y-ny*0.5).^2) <= rad_inner;
im = im + im2;

im = repmat(im,[1 1 nt nz]);

for k = 1:nt
    for m = 1:ny
        for n = 1:nx
            if (im(n,m,k,nz) == 2)
                im(n,m,k,nz) = C_a_mean(k);
            else
                if (im(n,m,k,nz) == 1)
                im(n,m,k,nz) = C_m_mean(k);
                end
            end
        end
    end 
end

% %Optional addition of noise
% choice = questdlg('Would you like to add noise?', ...
% 	'Noise addition', ...
% 	'Yes','No','No');
% switch choice
%     case 'Yes'
%         im = im + 10000*randn(nx,ny,nt,nz); %With gaussian noise
%     case 'No'
% end

SET(NO).IM = single(im);
filephantomhelper;
openfile('setupstacksfromdicom',NO);
segment('renderstacksfromdicom',NO);
SET(no).TimeVector = t_frames;
SET(no).TIncr = [];

%-------------------------------
function C_a = input_function(t)
%-------------------------------
%Helper function to PET phantom.
%Approximates the arterial concentration

if t < 15
    C_a = 0;
else 
    C_a = 10^5*log(0.1*t).*exp(-0.03*(t-15));
end
indices = C_a < 0; %find the elements which are negative
C_a(indices) = 0;  %set all the elements which are negative to zero.


