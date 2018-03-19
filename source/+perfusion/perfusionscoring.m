function varargout = perfusionscoring(varargin)
%GUI for fast MR perfusion analysis
macro_helper(varargin{:});
if nargin == 0 || isempty(varargin{1})
  varargin{1} = 'init';
end
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------
function init %#ok<DEFNU>
%------------
%Initiate GUI
global DATA SET


% Add no and scoringedit array to mygui properties
allnos = 1:numel(SET);
stressnos = find(strcmp({SET.ImageType},'Perfusion Stress Aligned'));
if isempty(stressnos)
  stressnos = find(strcmp({SET.ImageType},'Perfusion Stress'));
  hasalignedstress = false;
else
  hasalignedstress = true;
end

restnos = find(strcmp({SET.ImageType},'Perfusion Rest Aligned'));
if isempty(restnos)
  restnos = find(strcmp({SET.ImageType},'Perfusion Rest'));
  hasalignedrest = false;
else
  hasalignedrest = true;
end

restonly = false;
stressonly = false;
if hasalignedstress && hasalignedrest 
elseif hasalignedstress
  if yesno(['Found only Perfusion Stress Aligned image stack. Proceed ' ...
      'using only one stack for perfusion analysis?'])
    stressonly = true;
  else
    return
  end

  
elseif hasalignedrest
  if yesno(['Found only Perfusion Rest Aligned image stack. Proceed ' ...
      'using only one stack for perfusion analysis?'])
    %analysed = true;
    restonly = true;
    %stressnos = restnos;
  else
    return
  end
end

scarnos = findfunctions('findscarshortaxisno');

if (numel(restnos) ~= 1 || numel(stressnos) ~= 1)
  if isempty(restnos) && isempty(stressnos)
    myfailed(['Could not find either of stress/rest image stacks. Please ' ...
      'set image description of two stacks to Perfusion Stress and ' ...
      'Perfusion Rest before launching Perfusion module']);
    return
  elseif isempty(restnos) && ~hasalignedstress
    if yesno(['Found only Perfusion Stress image stack. Proceed ' ...
        'using only one stack for perfusion analysis?'])
      %restnos = stressnos;
      stressonly = true;
    else
      return
    end
  elseif isempty(stressnos) && ~hasalignedrest
    if yesno(['Found only Perfusion Rest image stack. Proceed ' ...
        'using only one stack for perfusion analysis?'])
      %stressnos = restnos;
      restonly = true;
    else
      return
    end
  else %none is empty
    mywarning(['Found more than one rest or stress image stack. ' ...
      'Taking first (arbitrary decision).']);
  end
end

DATA.GUI.PerfusionScoring = mygui('+perfusion/perfusionscoring.fig');
gui=DATA.GUI.PerfusionScoring;
gui.currenttag='stress';
gui.stressonly = stressonly;
gui.restonly = restonly;
gui.hidelv=0;


if stressonly
  set(gui.handles.stresstimebaraxes,'visible','on')
elseif restonly
  set(gui.handles.resttimebaraxes,'visible','on')
else
  set(gui.handles.stresstimebaraxes,'visible','on')
  set(gui.handles.resttimebaraxes,'visible','on')
end

if ~isempty(stressnos)
  gui.stressno = stressnos(1);
else
  gui.stressno = [];
end
if ~isempty(restnos)
  gui.restno = restnos(1);
else
  gui.restno = [];
end

if ~isempty(scarnos)
  gui.scarno = scarnos(1);
else
  gui.scarno = [];
  set(gui.handles.scaraxes, 'Visible', 'off')
end

%Determine which slices to use
inds=[gui.stressno, gui.restno];
[~,nslicestack]=min([SET(inds).ZSize]);
nslicestack = inds(nslicestack);
gui.stressslices=[];
gui.restslices=[];
gui.scarslices=[];

if gui.stressonly
  gui.stressslices=1:SET(gui.stressno).ZSize;
  
  gui.mrest = zeros(24,4);
  if isfield(SET(gui.stressno).PerfusionScoring,'mstress') && ~isempty(SET(gui.stressno).PerfusionScoring.mstress)
    gui.mstress = SET(gui.stressno).PerfusionScoring.mstress;
  else
    gui.mstress = zeros(24,4);
  end
  
elseif gui.restonly;
  gui.restslices=1:SET(gui.restno).ZSize;

  gui.mstress = zeros(24,4);
  
  if isfield(SET(gui.restno).PerfusionScoring,'mrest') && ~isempty(SET(gui.restno).PerfusionScoring.mrest)
    gui.mrest = SET(gui.restno).PerfusionScoring.mrest;
  else
    gui.mrest = zeros(24,4);
  end
  
else
  gui.stressslices=1:SET(nslicestack).ZSize;
  gui.restslices=1:SET(nslicestack).ZSize;
  
  [gui.stressslices,~,~,~,~,~] = segmentation('findmatchingslices', ...
    nslicestack,gui.stressno,0,0,0,0,0,0);
  [gui.restslices,~,~,~,~,~] = segmentation('findmatchingslices', ...
    nslicestack,gui.restno,0,0,0,0,0,0);

  if isfield(SET(gui.stressno).PerfusionScoring,'mstress') && ~isempty(SET(gui.stressno).PerfusionScoring.mstress)
    gui.mstress = SET(gui.stressno).PerfusionScoring.mstress;
  else
    gui.mstress = zeros(24,4);
  end
  
  if isfield(SET(gui.restno).PerfusionScoring,'mrest') && ~isempty(SET(gui.restno).PerfusionScoring.mrest)
    gui.mrest = SET(gui.restno).PerfusionScoring.mrest;
  else
    gui.mrest = zeros(24,4);
  end
  
end

if ~isempty(gui.scarno)
  [gui.scarslices,~,~,~,~,~] = segmentation('findmatchingslices', ...
    nslicestack,gui.scarno,0,0,0,0,0,0);
end

if ~isrow(gui.stressslices)
  gui.stressslices=gui.stressslices';
end
if ~isrow(gui.restslices)
  gui.restslices=gui.restslices';
end
if ~isrow(gui.scarslices)
  gui.scarslices=gui.scarslices';
end

if ~isempty(gui.scarno)
  if isfield(SET(gui.scarno).PerfusionScoring,'mscar') && ~isempty(SET(gui.scarno).PerfusionScoring.mscar)
    gui.mscar=SET(gui.scarno).PerfusionScoring.mscar;
  else
    gui.mscar = zeros(24,4);
  end
else
    gui.mscar = zeros(24,4);
end

gui.mdiff = gui.mstress-gui.mrest;

fcn = @(hObject,eventdata)perfusion.perfusionscoring('timebaraxes_ButtonDownFcn',hObject,eventdata);
set(gui.fig,'WindowButtonDownFcn',fcn);

segment('recursekeypressfcn',gui.fig,@(hObject,eventdata)perfusion.perfusionscoring('keypress_Callback',eventdata))
load('newicons.mat','newicons')

gui.Icons=newicons;
gui.iconholder = myiconplaceholder(gui.handles.iconaxes,0,1,gui.fig);
initiconholder;

%indent autozoom
indent(gui.iconholder,'autozoom',0);

% %estimate rotation;
% impos1 = SET(gui.stressno).ImagePosition;
% impos2 = SET(gui.restno).ImagePosition;
% imo1 = SET(gui.stressno).ImageOrientation(1:3);
% imo2 = SET(gui.restno).ImageOrientation(1:3);
% 
% theta = acos(imo1/norm(imo1)*imo2'/norm(imo2))/pi*180;
% 
% im1 = SET(gui.stressno).IM(:,:,1,1);
% im2 = SET(gui.restno).IM(:,:,1,1);
% im3 = SET(gui.scarno).IM(:,:,1,1);

% fig = figure;
% subplot(121)
% imagesc(im1);
% subplot(122)
% imagesc(imrotate(im2,rad))
 gui.autozoom=1;

set(gui.fig,'WindowButtonMotionFcn','perfusion.perfusionscoring(''motionfunc'')');
gui.resttf=1;
gui.stresstf=1;
gui.scartf=1;
gui.restn=[];
gui.stressn=[];
gui.scarn=[];
%gui.cinetf=1;

if ~isempty(gui.restno)
  gui.restn=SET(gui.restno).TSize;
end
if ~isempty(gui.stressno)
  gui.stressn=SET(gui.stressno).TSize;
end
if ~isempty(gui.scarno)
  gui.scarn=SET(gui.scarno).TSize;
end
% if ~isempty(gui.cineno)
%   gui.cinen=SET(gui.cineno).TSize;
% end

%gui.cineplay=0;
gui.restplay=0;
gui.stressplay=0;
gui.play=0;

% for type={'stress','rest','scar'}
%   generateimages(type{1});
% end


%Initiate imageaxis and timebar
for imname = {'stress','rest'}
  autozoom(imname{1});
  generateimages(imname{1});
  initimageaxis(imname{:});
  inittimebar(imname{1})
  settimeframe('currenttime',1,imname{1})
  %inittimebar(imname{:},analysed);
end

autozoom('scar');
generateimages('scar');
initimageaxis('scar');
generatebullseye('stress')
generatebullseye('rest')
%generatebullseye('scar')
generatebullseye('diff')

% for type={'Stress','Rest'}
% title(gui.handles.(['bullseyeaxes',lower(type{1})]),['\fontsize{14}',type{1}])
% end
% title(gui.handles.diffaxes,['\fontsize{14}','Diff'])

save2set


%---------------------------
function generateimages(type) 
%---------------------------

global SET DATA
gui = DATA.GUI.PerfusionScoring;
no = gui.([type,'no']);
if isempty(no)
  return;
end
      
    
xsz=SET(no).XSize;
ysz=SET(no).YSize;
T = gui.([type,'n']);
slices = gui.([type,'slices']);
%this was taken from a really nicely orientated cine SAX
%  0.7598    0.6030   -0.2431   -0.4615    0.2369   -0.8549
imoref = [0.7598,    0.6030,   -0.2431];
imorefn= [-0.4579    0.7617    0.4583];
imo = SET(no).ImageOrientation(1:3)/norm(SET(no).ImageOrientation(1:3));

    
%     figure;
%     plot3([0 imo(1)],[0 imo(2)],[0 imo(3)],'r');
imo=imo-(imorefn*imo')*imorefn;
%       hold on;
%     plot3([0 imo(1)],[0 imo(2)],[0 imo(3)],'c');
%
%     plot3([0 imorefn(1)],[0 imorefn(2)],[0 imorefn(3)],'g')


theta= acos(imoref/norm(imoref)*imo'/norm(imo))/pi*180;
xlim = gui.([type, 'xlim']);
ylim = gui.([type, 'ylim']);
% xscale=gui.([type,'xscale']);
% yscale=gui.([type,'yscale']);
scale =gui.([type,'scale']);
im=[];
for tf=1:T
  tmp=[];
  for i = slices
    tmp= cat(2,tmp, imresize(imrotate(squeeze(SET(no).IM(xlim,ylim,tf,i)),theta,'bilinear','crop'),scale*[xsz ysz],'bilinear'));
  end
  im=cat(3,im,tmp);
end

%     for slice = slices
%       im=imresize(imrotate(SET(no).IM(xlim,ylim,:,slice),theta),[]);
%     end
    gui.([type 'im'])=im;
%     figure;subplot(121)
%     imagesc(im(:,:,1))
%     subplot(122)
%     imagesc(SET(1).IM(:,:,8,1))
    
%---------------------
function resize_fcn
%---------------------

try
  global DATA
  persistent chk
  
  if ~isempty(chk)
    return
  else
    chk=1;
  end
  
  if isfield(DATA.GUI,'PerfusionScoring')
    gui = DATA.GUI.PerfusionScoring;
    try
      render(gui.iconholder);
      drawnow
    catch
      %it wasnt there
    end
  end
  
  chk=[];
catch
%no data loaded  
end

%---------------------
function motionfunc
%-----------------------
global DATA

gui = DATA.GUI.PerfusionScoring;
motion(gui.iconholder);

% %-------------------------
%  function hidecontour
% %--------------------------
% global DATA
% gui = DATA.GUI.PerfusionScoring;
% gui.hidelv=~gui.hidelv;
% for type = {'stress','rest'};
%   drawimages(type{1})
% end

%---------------------
function initiconholder
%-----------------------
global DATA
gui = DATA.GUI.PerfusionScoring;

iconcell={};
iconcell{1,end+1}=myicon('playstress',gui.iconholder,gui.Icons.playstress,'Play stress stack',@() perfusion.perfusionscoring('play','stress'),2,1);
iconcell{1,end+1}=myicon('playrest',gui.iconholder,gui.Icons.playrest,'Play rest stack',@() perfusion.perfusionscoring('play','rest'),2,1);
iconcell{1,end+1}=myicon('playall',gui.iconholder,gui.Icons.play,'Play all stacks',@()  perfusion.perfusionscoring('playall'),2,1);
%iconcell{1,end+1}=myicon('hidesegmentation',gui.iconholder,gui.Icons.hidelv,'Hide icons',@() perfusion.perfusionscoring('hidecontour'),2);
iconcell{1,end+1}=myicon('autozoom',gui.iconholder,gui.Icons.autozoom,'Auto zoom',@() perfusion.perfusionscoring('autozoomtoggle'),2);
add(gui.iconholder,iconcell)%gui.iconholder.add(iconcell);
%render(gui.iconholder)

%----------------------------
function autozoomtoggle
%---------------------- --
 global DATA SET
 gui = DATA.GUI.PerfusionScoring;
 
 if gui.autozoom
   gui.autozoom=0;
   if ~isempty(gui.scarno)
     gui.scarxlim=1:SET(gui.scarno).XSize;
     gui.scarylim=1:SET(gui.scarno).YSize;
   else
     gui.scarxlim=[];
     gui.scarylim=[];
   end
   
   if ~isempty(gui.stressno)
     gui.stressxlim=1:SET(gui.stressno).XSize;
     gui.stressylim=1:SET(gui.stressno).YSize;
   else
     gui.stressxlim=[];
     gui.stressylim=[];
   end
   
   if ~isempty(gui.restno)
     gui.restxlim=1:SET(gui.restno).XSize;
     gui.restylim=1:SET(gui.restno).YSize;
   else
     gui.restxlim=[];
     gui.restylim=[];
   end
   
   for type = {'scar', 'rest', 'stress'}
 generateimages(type{1});
     initimageaxis(type{1})
   end
   
 else
   
   gui.autozoom=1;
   for type = {'scar', 'rest', 'stress'}
     autozoom(type{1});
     generateimages(type{1});
      initimageaxis(type{1});
   end
 end
 colormap([153/255 204/255 1;1 0.95 0;1 0 0])
%----------------------------
 function autozoom(type)
%-----------------------
%function which finds larges epicardial contour and crops image a fixed
%distance from it.
global DATA %SET

gui = DATA.GUI.PerfusionScoring;

no=gui.([type,'no']);

if isempty(no)
  return
end
cinesaxno = findfunctions('findcineshortaxisno');
if ~isempty(cinesaxno)
  %[~,xlim,ylim] = segment('getbox',cinesaxno,no);
  [xlim,ylim] = segment('getbox',cinesaxno,no,1);
else
  [xlim,ylim] = segment('getbox',no,no,1);
end
gui.([type,'xlim']) = xlim(1):xlim(end);
gui.([type,'ylim']) = ylim(1):ylim(end);
if ylim(end)-ylim(1)>=xlim(end)-xlim(1)
   gui.([type,'scale']) = 180/(ylim(end)-ylim(1)+1);
 else
   gui.([type,'scale']) = 180/(xlim(end)-xlim(1)+1);
  end

% 
% slices=gui.([type,'slices']);
% xsz=SET(no).XSize;
% ysz=SET(no).YSize;
% 
% if ~isempty(SET(no).EpiX) && ~all(all(all(isnan(SET(no).EpiX(:,:,slices))))) 
% 
% epix = SET(no).EpiX(:,:,slices);
% epiy = SET(no).EpiY(:,:,slices);
% 
% xmin = floor(min(epix(:)));
% ymin = floor(min(epiy(:)));
% xmax = ceil(max(epix(:)));
% ymax = ceil(max(epiy(:)));
% 
% dia = round(norm([xmax,ymax]-[xmin,ymin])/2*0.75);
% xmin=xmin-dia;
% ymin=ymin-dia;
% xmax=xmax+dia;
% ymax=ymax+dia;
% 
% %make square cut?
% %min(abs([xmin-xmax,ymin-ymax]))
% 
% if xmin<1
%   xmin=1;
% end
% 
% if ymin<1
%   ymin=1;
% end
% 
% if xmax>xsz
%   xmax=xsz;
% end
% 
% if ymax>ysz
%   ymax=ysz;
% end
% 
% %this catches segmentation blow ups
% if ymax-ymin<20
%   ymin=1;
%   ymax=ysz;
% end
% 
% if xmax-xmin<20
%   xmin=1;
%   xmax=xsz;
% end
% 
% if ymax-ymin>=xmax-xmin
%   gui.([type,'yscale']) = 180/(ymax-ymin+1);
%   gui.([type,'xscale']) = gui.([type,'yscale']);%*(ymax-ymin+1)/(xmax-xmin+1);%120/(ymax-ymin+1);
% else
%   gui.([type,'xscale']) = 180/(xmax-xmin+1);
%   gui.([type,'yscale']) =  gui.([type,'xscale']);%*(xmax-xmin+1)/(ymax-ymin+1);
% end
%   gui.([type,'xlim']) = xmin:xmax;
%   gui.([type,'ylim']) = ymin:ymax;
% 
% else
%   im= SET(no).IM(:,:,:,slices);
%   [xlim,ylim] = autocrop(im);
%   
%   if isempty(xlim)
%     xmin=1;
%     xmax=xsz;
%   else
%     xmin=xlim(1);
%     xmax=xlim(end);
%   end
%   
%   if isempty(ylim)
%     ymin=1;
%     ymax=ysz;
%   else
%     ymin=ylim(1);
%     ymax=ylim(end);
%   end
% 
% if ymax-ymin>=xmax-xmin
%   gui.([type,'yscale']) = 180/(ymax-ymin+1);
%   gui.([type,'xscale']) = gui.([type,'yscale']);%*(ymax-ymin+1)/(xmax-xmin+1);%120/(ymax-ymin+1);
% else
%   gui.([type,'xscale']) = 180/(xmax-xmin+1);
%   gui.([type,'yscale']) =  gui.([type,'xscale']);%*(xmax-xmin+1)/(ymax-ymin+1);
% end
% 
%   gui.([type,'xlim']) = xmin:xmax;
%   gui.([type,'ylim']) = ymin:ymax;
% end  

%---------------------------------------
function togglemode(type)
%---------------------------------------
global DATA
gui = DATA.GUI.PerfusionScoring;
gui.currenttag=type;
     


%----------------------------
function initimageaxis(field)
%----------------------------
%Initiate image axis with images from current stack
global DATA SET
gui = DATA.GUI.PerfusionScoring;
handles = gui.handles;
h = handles.([field 'axes']);
no = gui.([field 'no']);
tf = gui.([field 'tf']);
slices = gui.([field,'slices']);

cla(h)

if isempty(no)
  return
end

for hloop = h'
  hold(hloop,'on');
  %axis(hloop,[0 1 0 1]);
  axis(hloop,'ij');%'equal'
  if ~isempty(SET(no).Colormap)
    colormap(hloop,SET(no).Colormap)
  else
    colormap(hloop,'gray')
  end
end

% xlim = gui.([field, 'xlim']);
% ylim = gui.([field, 'ylim']);
% xsz=xlim(end)-xlim(1);
% ysz=ylim(end)-ylim(1);
% 
% xscale=gui.([field,'xscale']);
% yscale=gui.([field,'yscale']);
% %gui.([field,'im'])=[];
% im=[];
%for i = slices
 %im = [im, imresize(squeeze(SET(no).IM(xlim,ylim,tf,i)),[xscale*xsz yscale*ysz],'bilinear')];
%end
im = gui.([field, 'im']);
im=im(:,:,tf);
%plot images
cmap = gray(256);
c = SET(no).IntensityMapping.Contrast;
b = SET(no).IntensityMapping.Brightness;
rim = segment('remap',im,cmap(:,1),c,b);
gim = segment('remap',im,cmap(:,2),c,b);
bim = segment('remap',im,cmap(:,3),c,b);
im = cat(3,rim, gim, bim);
 gui.([field,'imhandle']) = image(im,'parent',h);
 set(gui.([field,'imhandle']),'buttondownfcn',sprintf('perfusion.perfusionscoring(''togglemode'',''%s'')',field))
 axis(h,'image') 
% xlim(h,gui.([field,'xlim']))
% ylim(h,gui.([field,'ylim']))
 
if ~strcmp('scar',field)
pos = plotboxpos(gui.handles.([field,'axes'])); 
set(gui.handles.([field,'timebaraxes']),'position',[pos(1)+pos(3)/5,pos(2)-0.027,3*pos(3)/5,0.02])
end

%initcontours
  %zsz = SET(gui.([field,'no'])).ZSize;
%   xsz = SET(no).XSize;
%  ysz = SET(gui.([field,'no'])).YSize;
% xres = SET(gui.([field,'no'])).ResolutionX;
% yres = SET(gui.([field,'no'])).ResolutionY;
% set(gui.handles.([field,'axes']),'plotboxaspectratio',[zsz*ysz*yres xsz*xres 1]);
%     xlim = gui.([field,'xlim']);
%     ylim = gui.([field,'ylim']);
%     xmin = xlim(1);
%     ymin = ylim(1);
%     ymax = ylim(end);
%     
%     xscale = gui.([field,'xscale']);
%     yscale = gui.([field,'yscale']);
% if ~isempty(SET(no).EndoX)
%   x=[];
%   y=[];
%   for i = 1:length(slices)
%     x = [x;nan;xscale*(SET(no).EndoX(:,tf,slices(i))-xmin)+1];
%     y = [y;nan;yscale*(SET(no).EndoY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
%   end
%  gui.([field,'endoseghandle']) =  plot(h,y,x,'r');
% else
%   gui.([field,'endoseghandle']) =  plot(h,nan,nan,'r');
% end
% 
% if ~isempty(SET(no).EpiX)
%   x=[];
%   y=[];
%   for i = 1:length(slices)
%     x = [x;nan;xscale*(SET(no).EpiX(:,tf,slices(i))-xmin)+1];
%     y = [y;nan;yscale*(SET(no).EpiY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
%   end
%   gui.([field,'episeghandle']) =  plot(h,y,x,'g');
% else
%   gui.([field,'episeghandle']) =  plot(h,nan,nan,'g');
% end
% 
% if ~isempty(SET(no).RVEndoX)
%   x=[];
%   y=[];
%   for i = 1:length(slices)
%     x = [x;nan;xscale*(SET(no).RVEndoX(:,tf,slices(i))-xmin)+1];
%     y = [y;nan;yscale*(SET(no).RVEndoY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
%   end
%   gui.([field,'rvendoseghandle']) =  plot(h,y,x,'m');
% else
%   gui.([field,'rvendoseghandle']) =  plot(h,nan,nan,'m');
% end
% 
% if ~isempty(SET(no).RVEpiX)
%   x=[];
%   y=[];
%   for i = 1:length(slices)
%     x = [x;nan;xscale*(SET(no).RVEpiX(:,tf,slices(i))-xmin)+1];
%     y = [y;nan;yscale*(SET(no).RVEpiY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
%   end
%   gui.([field,'rvepiseghandle']) =  plot(h,y,x,'c');
% else
%   gui.([field,'rvepiseghandle']) =  plot(h,nan,nan,'c');
% end


 %drawimages(field);
%-----------------------------------------
function autodetectstartstop
%------------------------------------
global DATA SET

gui = DATA.GUI.PerfusionScoring;

activation = zeros(1,gui.stressn);

for t=1:gui.stressn
  activation(t) = sum(sum(sum(squeeze(SET(gui.stressno).IM(:,:,t,gui.stressslices)))));
end

%-----------------------------------
function inittimebar(field)
%-----------------------------------
%Initiate timebar axis for image specified by input parameter 'field'.
global DATA SET
gui = DATA.GUI.PerfusionScoring;
h = gui.handles.([field 'timebaraxes']);
no = gui.([field 'no']);

if isempty(no)
  return
end

if gui.stressonly && strcmp(field,'rest')
return
end

if gui.restonly && strcmp(field,'stress')
return
end

delete(get(h,'Children'));

tvec = SET(no).TimeVector;

hold(h,'on');
fcn = @(hObject,eventdata)perfusion.perfusionscoring('timebar_ButtonDownFcn',hObject,eventdata,field);

%Draw timebar (red) and set its buttondown fcn
timebar = plot(h,tvec(SET(no).CurrentTimeFrame)*[1 1],[0 1],'r','Tag','currenttime','linewidth',2);
set(timebar,'ButtonDownFcn',fcn);
set(h,'fontsize',14)
%Draw start and end bars (blue) and set buttondown fcns

%   endtime = tvec(SET(no).EndAnalysis);
%   plot(h,endtime*[1 1],[0 1],'b','Tag','endtime','ButtonDownFcn',fcn,'linewidth',2);
%   text(endtime,0.9,'End','Parent',h,'Tag','endtext','ButtonDownFcn',fcn);
% 
%   starttime = tvec(SET(no).StartAnalysis);
%   plot(h,starttime*[1 1],[0 1],'b','Tag','starttime','ButtonDownFcn',fcn,'linewidth',2);
%   text(starttime,0.9,'Start','Parent',h,'Tag','starttext',...
%     'ButtonDownFcn',fcn,'HorizontalAlignment','right');


updatetimebar(field);

%----------------------------------
function settimeframe(tag,tf,field)
%----------------------------------
%Sets current/start/end timeframe for image stack specified by input 
%parameter 'field'
global DATA SET
gui = DATA.GUI.PerfusionScoring;
no = gui.([field 'no']);

if isempty(no)
  return
end

switch tag
  case 'currenttime'
    SET(no).CurrentTimeFrame = max(min(tf,SET(no).TSize),1);
    gui.([field,'tf']) = SET(no).CurrentTimeFrame;
  case {'starttime','starttext'}
    SET(no).StartAnalysis = max(min(tf,SET(no).EndAnalysis),1);
  case {'endtime','endtext'}
    SET(no).EndAnalysis = max(min(tf,SET(no).TSize),SET(no).StartAnalysis);
end
drawimages(field);
updatetimebar(field);


%----------------------------------------------------
  function keypress_Callback(evnt)
    %---------------------------------------------------
    global DATA
    gui = DATA.GUI.PerfusionScoring;
    currenttag = get(get(gui.fig,'currentaxes'),'tag');
    
    if isempty(currenttag) || ~strcmp(currenttag,'iconaxes')
      currenttag=gui.currenttag;
    end
    
    key=getkey(evnt);
    
    switch key
      case 'rightarrow'
        
        if regexp(currenttag,'stress')
          
          tf = gui.stresstf+1;
          if tf>gui.stressn
            tf=1;
          end
          settimeframe('currenttime',tf,'stress')
        end
        
        if regexp(currenttag,'rest')
          tf=gui.resttf+1;
          if tf>gui.restn
            tf=1;
          end
          settimeframe('currenttime',tf,'rest')
        end
        
      case 'leftarrow'
        
        if regexp(currenttag,'stress')
          tf = gui.stresstf-1;
          if tf<1
            tf=gui.stressn;
          end
          settimeframe('currenttime',tf,'stress')
        end
        
        if regexp(currenttag,'rest')
          tf=gui.resttf-1;
          if tf<1
            tf=gui.restn;
          end
          settimeframe('currenttime',tf,'rest')
        end
        
      case  'shift-rightarrow'
        if gui.restonly
          ref_tf=(gui.resttf+1)/gui.restn;
          
          if ref_tf>1;
            ref_tf=1/gui.restn;
          end
          
        else
          ref_tf=(gui.stresstf+1)/gui.stressn;
          
          if ref_tf>1;
            ref_tf=1/gui.stressn;
          end
        end
        for type = {'stress','rest'}%,'cine'}
          % pick tf closest to current stress or rest if restonly frame
          T=gui.([type{1}, 'n']);
          [~,tf] = min(abs(ref_tf-(1:T)/T));
          settimeframe('currenttime',tf,type{1})
        end
      case  'shift-leftarrow'
        if gui.restonly
          ref_tf=(gui.resttf-1)/gui.restn;
          
          if ref_tf<1/gui.restn;
            ref_tf=1;
          end
          
        else
          ref_tf=(gui.stresstf-1)/gui.stressn;
          
          if ref_tf<1/gui.stressn;
            ref_tf=1;
          end
        end
        for type = {'stress','rest'}%,'cine'}
          % pick tf closest to current stress or rest if restonly frame
          T=gui.([type{1}, 'n']);
          [~,tf] = min(abs(ref_tf-(1:T)/T));
          settimeframe('currenttime',tf,type{1})
        end
      case 'p'
        if ~gui.play
          indent(gui.iconholder,'playall',1)
        else
          undent(gui.iconholder,'playall',1)
        end
    end

%----------------------------------------------------
function timebaraxes_ButtonDownFcn(hObject, ~)
%----------------------------------------------------
%Buttondown function for timebar axes of image specified by input 
%parameter 'field'. Changes current timeframe to the one closest to
%position of clicked point
global DATA SET
gui = DATA.GUI.PerfusionScoring;
handleAddress=hittest(hObject);
try
switch get(handleAddress,'tag')
  case 'stresstimebaraxes'
    field='stress';
  case 'resttimebaraxes'
    field='rest';
  otherwise
    return
end
catch
return;  
end
gui.currenttag=field;
no = gui.([field 'no']);
[x,y] = mygetcurrentpoint(handleAddress);
[~,tvix] = min(abs(SET(no).TimeVector-x));
settimeframe('currenttime',tvix,field);
kids=get(handleAddress,'children');
obj=kids(strcmp('currenttime',get(kids,'tag')));

if isempty(obj)
  return 
end

motionfcn = @(hObject,eventdata)perfusion.perfusionscoring('timebaraxes_MotionFcn',hObject,eventdata,obj,no,field,1);
set(gui.fig,'WindowButtonMotionFcn',motionfcn);
buttonupfcn = @(hObject,eventdata)perfusion.perfusionscoring('timebaraxes_ButtonUpFcn',hObject,eventdata);
set(gui.fig,'WindowButtonUpFcn', buttonupfcn);

%------------------------------------------------
function timebar_ButtonDownFcn(hObject, ~, field)
%------------------------------------------------
%Buttondown function for graphical timebar object of image specified 
%by input parameter 'field'. Activates dragging of timebars.
global DATA
gui = DATA.GUI.PerfusionScoring;
no = gui.([field 'no']);
obj = hObject;
motionfcn = @(hObject,eventdata)perfusion.perfusionscoring('timebaraxes_MotionFcn',hObject,eventdata,obj,no,field);
set(gui.fig,'WindowButtonMotionFcn',motionfcn);
buttonupfcn = @(hObject,eventdata)perfusion.perfusionscoring('timebaraxes_ButtonUpFcn',hObject,eventdata);
set(gui.fig,'WindowButtonUpFcn', buttonupfcn);


%-----------------------------------------------------------
function timebaraxes_MotionFcn(hObject, ~, tbobj, no, field,axesclick)
%-----------------------------------------------------------
%Mouse motion function for timebar axes of image specified by input 
%parameter 'field'. Used for dragging timebars to change current
%timeframe or start/end points of timeframes in which to align images.
global SET

if nargin<6
  axesclick=0;
end

x = mygetcurrentpoint(get(tbobj,'Parent'));
[~,tvix] = min(abs(SET(no).TimeVector-x));
settimeframe(get(tbobj,'Tag'),tvix,field);
  
%-------------------------------------------
function timebaraxes_ButtonUpFcn(hObject, ~)
%-------------------------------------------
%Buttonup function for timebar axes of image specified by input 
%parameter 'field'. Deactivates dragging of timebar.
set(hObject,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[]);



%---------------------------
function updatetimebar(field)
%---------------------------
%Update timebar axis specified by handle 'h' from input arguments
global DATA SET
gui = DATA.GUI.PerfusionScoring;
h = gui.handles.([field 'timebaraxes']);
no = gui.([field,'no']);

fcn = get(h,'ButtonDownFcn');
tvec = SET(no).TimeVector;
tmin = tvec(1);
tmax = tvec(end);

starttime = tvec(SET(no).StartAnalysis); 
endtime = tvec(SET(no).EndAnalysis); 
currenttime = tvec(SET(no).CurrentTimeFrame); %#ok<NASGU>

%Update timebars
for kid = get(h,'Children')'
  tag = get(kid,'Tag');
  switch tag
    case {'starttime','endtime','currenttime'}
      eval(sprintf('set(kid,''XData'',%s*[1 1])',tag));
    case 'starttext'
      set(kid,'Position',[starttime 0.9 0]);
    case 'endtext'
      set(kid,'Position',[endtime 0.9 0]);
  end
end

%Set axes options
marg = (tmax-tmin)/100;
axis(h,[tmin-marg tmax+marg 0 1]);
tstep = 5*ceil((tmax-tmin)/40);
tickvec = [ceil(tmin/tstep)*tstep:tstep:floor(tmax/tstep)*tstep];
set(h, 'Xtick', tickvec,'YTick',[],'XMinorTick','on');
set(h,'ButtonDownFcn',fcn);

%-------------------------
function drawimages(field)
%-------------------------
%Do an update of all image axes
global DATA SET
gui = DATA.GUI.PerfusionScoring;
no = gui.([field 'no']);
tf=gui.([field,'tf']);

if isempty(no)
  return
end

if strcmp(field,'stress') && gui.restonly || ...
    strcmp(field,'rest') && gui.stressonly
  return
end

h = gui.handles.([field 'axes']);

%gui.([field,'im'])=[];

xlim = gui.([field, 'xlim']);
ylim = gui.([field, 'ylim']);
% xsz=xlim(end)-xlim(1);
% ysz=ylim(end)-ylim(1);

% xscale=gui.([field,'xscale']);
% yscale=gui.([field,'yscale']);

% for i = gui.([field,'slices'])
%  gui.([field,'im']) = [gui.([field,'im']),  imresize(squeeze(SET(no).IM(gui.([field,'xlim']),gui.([field,'ylim']),tf,i)),[xscale*xsz yscale*ysz],'bilinear')];%squeeze(SET(gui.([field,'no'])).IM(xlim,ylim,tf,i))];
% end
im=gui.([field,'im']);
im=im(:,:,tf);

%plot images
cmap = gray(256);
c = SET(gui.([field,'no'])).IntensityMapping.Contrast;
b = SET(gui.([field,'no'])).IntensityMapping.Brightness;
rim = segment('remap',im,cmap(:,1),c,b);
gim = segment('remap',im,cmap(:,2),c,b);
bim = segment('remap',im,cmap(:,3),c,b);
im = cat(3,rim, gim, bim);
 set(gui.([field,'imhandle']),'Cdata', im)
%  xlim(h,gui.([field,'xlim']))
%  ylim(h,gui.([field,'ylim']))
 %hidelv button indented

   %drawcontours(field)
%    gui.endoseghandle =  plot(SET(gui.([field,'no'])).EndoX,SET(gui.([field,'no'])).EndoY);
%    gui.episeghandle =  plot(SET(gui.([field,'no'])).EpiX,SET(gui.([field,'no'])).EpiY);
%    gui.rvendoseghandle =  plot(SET(gui.([field,'no'])).RVEndoX,SET(gui.([field,'no'])).RVEndoY);
%    gui.rvepiseghandle =  plot(SET(gui.([field,'no'])).RVEpiX,SET(gui.([field,'no'])).RVEpiY);
%  end
 %--------------------------------------
  function drawcontours(field)
    %-----------------------------------
    global DATA SET
    
    gui = DATA.GUI.PerfusionScoring;
    no = gui.([field 'no']);
    tf = gui.([field 'tf']);
    slices = gui.([field, 'slices']);
    xlim = gui.([field,'xlim']);
    ylim = gui.([field,'ylim']);
    xmin = xlim(1);
    ymin = ylim(1);
    ymax = ylim(end);
    
    scale = gui.([field,'scale']);
    
    if ~isempty(SET(no).EndoX) && ~gui.hidelv
      x=[];
      y=[];
      for i = 1:length(slices)
        x = [x;nan;scale*(SET(no).EndoX(:,tf,slices(i))-xmin)+1];
        y = [y;nan;scale*(SET(no).EndoY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
      end
        set(gui.([field,'endoseghandle']),'XData',y,'YData',x)%plot(h,x,y);
    else
        set(gui.([field,'endoseghandle']), 'XData',nan,'YData',nan)
    end
    
    if ~isempty(SET(no).EpiX) && ~gui.hidelv
      x=[];
      y=[];
      for i = 1:length(slices)
        x = [x;nan;scale*(SET(no).EpiX(:,tf,slices(i))-xmin)+1];
        y = [y;nan;scale*(SET(no).EpiY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
      end
        set(gui.([field,'episeghandle']),'XData',y,'YData',x)%plot(h,x,y);
    else
      set(gui.([field,'episeghandle']), 'XData',nan,'YData',nan)
    end
    
    if ~isempty(SET(no).RVEndoX) && ~gui.hidelv
      x=[];
      y=[];
      for i = 1:length(slices)
        x = [x;nan;scale*(SET(no).RVEndoX(:,tf,slices(i))-xmin)+1];
        y = [y;nan;scale*(SET(no).RVEndoY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
      end
      set(gui.([field,'rvendoseghandle']),'XData',y,'YData',x)%plot(h,x,y);
    else
      set(gui.([field,'rvendoseghandle']), 'XData',nan,'YData',nan)
    end
    
    if ~isempty(SET(no).RVEpiX) && ~gui.hidelv
      x=[];
      y=[];
      for i = 1:length(slices)
        x = [x;nan;scale*(SET(no).RVEpiX(:,tf,slices(i))-xmin)+1];
        y = [y;nan;scale*(SET(no).RVEpiY(:,tf,slices(i))+(ymax-ymin)*(i-1)-ymin)+1+(i-1)];
      end
      set(gui.([field,'rvepiseghandle']),'XData',y,'YData',x)%plot(h,x,y);
    else
      set(gui.([field,'rvepiseghandle']), 'XData',nan,'YData',nan)
    end
    
    %-------------------
    function playcheck(type)
    %-----------------
    global DATA 
      gui = DATA.GUI.PerfusionScoring;

    if ~strcmp(type,'rest') && gui.restplay
        %undent and stop rest play
        gui.restplay=0;
        undent(gui.iconholder.iconCell{2})
        render(gui.iconholder)
      end
      
      if ~strcmp(type,'stress') && gui.stressplay
        gui.stressplay=0;
        undent(gui.iconholder.iconCell{1})
        render(gui.iconholder)  
      end
      
      if ~strcmp(type,'play') && gui.play
        gui.play=0;
        undent(gui.iconholder.iconCell{3})
        render(gui.iconholder)
      end
      
      %---------------------------------
      function play(type)
        %--------------------------------
        global DATA SET
        
        gui = DATA.GUI.PerfusionScoring;
        no=gui.([type, 'no']);
        
        gui.([type,'play'])=~gui.([type,'play']);
        
        gui.currenttag=type;
        playcheck(type)
        
        if isempty(no)
          return
        end
        
        try
        while gui.([type,'play'])
          gui.([type, 'tf']) = gui.([type, 'tf'])+1;
          
          if gui.([type, 'tf'])>gui.([type, 'n'])
            gui.([type, 'tf'])=1;
          end
          
          SET(no).CurrentTimeFrame=gui.([type, 'tf']);
          
          drawimages(type);
          updatetimebar(type);
          drawnow
          %pause(0.001)
        end
        catch
          %gui has been closed.
          
           mydisp('closed while playing')
          return
        end

%---------------------------------
function playall
%--------------------------------
global DATA SET

gui = DATA.GUI.PerfusionScoring;

% if gui.stressonly
%   play('stress')
%   return
% end
% 
% if gui.restonly
%   play('rest')
%   return
% end

gui.play=~gui.play;

playcheck('play')

try
while gui.play
  if gui.restonly
    ref_tf=(gui.resttf+1)/gui.restn;
  
  if ref_tf>1;
    ref_tf=1/gui.restn;
  end
  
  else
    ref_tf=(gui.stresstf+1)/gui.stressn;
  
  if ref_tf>1;
    ref_tf=1/gui.stressn;
  end
  end
  
  if gui.stressonly
    types = {'stress'};
  elseif gui.restonly
    types={'rest'};
  else
    types={'stress','rest'};
  end
  
  for type = types%,'cine'}
    % pick tf closest to current stress or rest if restonly frame
    T=gui.([type{1}, 'n']);
    [~,gui.([type{1}, 'tf'])] = min(abs(ref_tf-(1:T)/T));%SET(gui.([type{1}, 'no'])).TSize));
    SET(gui.([type{1}, 'no'])).CurrentTimeFrame=gui.([type{1}, 'tf']);
    drawimages(type{1});
    updatetimebar(type{1});
  end
  drawnow
  %pause(0.001)
end
catch
  %closed while playing
  mydisp('closed while playing')
end

%--------------------
function score(type)
%--------------
global DATA

gui=DATA.GUI.PerfusionScoring;

switch type
  case 'stress'
    gui.currenttag=type;
    bullseyehandle=gui.handles.bullseyeaxesstress;
  case 'rest'
    gui.currenttag=type;
    bullseyehandle=gui.handles.bullseyeaxesrest;
  case 'scar'
    return
    %bullseyehandle=gui.handles.bullseyeaxesscar;
end
%first identify which region was clicked then update the correct m values 
[x,y] = mygetcurrentpoint(bullseyehandle);

xc=x-gui.center;
yc=y-gui.center;

[theta,rho] = cart2pol(xc,yc);

%determine slice

%apex case completely handled here
if rho < 50
 %apex adjust m(1:4) then return
 region=1;
 sector='apex';
end

if rho >50 && rho <100
 region=2;
end

if rho >100 && rho <150
region=3;
end

if rho >150 && rho <200
  region=4;
end

if rho>200
  return
end

%apical
if region == 2
  if  abs(theta)<pi/4
    %lateral
    sector = 'lateral';
  end
  
  if theta >pi/4 && theta <3*pi/4
    %inferior
    sector = 'inferior';
  end
  
  if abs(theta) > 3*pi/4% || (theta >-pi/4 && theta < 0)
    %septal
    sector = 'septal';
  end
  
  if theta <-pi/4 && theta >-3*pi/4
    %anterior
    sector ='anterior';
  end
end


%basal mid
if region >2
  if theta<pi/3 && theta >0
    %inferolateral
    sector = 'inferolateral';
  end
  
  if theta >pi/3 && theta <2*pi/3
    %inferior
    sector = 'inferior';
  end
  
  if theta >2*pi/3 && theta <pi
    %inferoseptal
    sector = 'inferoseptal';
  end
  
  if theta < 0 && theta > -pi/3
    %anteroseptal
    sector = 'anterolateral';
  end
  
  if theta < -pi/3 && theta >-2*pi/3
    %anterior
    sector = 'anterior';
  end
  
  if theta < -2*pi/3 && theta >-pi
    %anterolateral
    sector = 'anteroseptal';
  end
end

clicktype = get(get(bullseyehandle,'parent'),'selectiontype');

switch clicktype
   case 'open'
      setsectionscore(type,region,sector,2)
  case 'normal' %ordinary mouseclick add
    setsectionscore(type,region,sector,1)
  case {'alt','extend'}
    setsectionscore(type,region,sector,0)
end

generatebullseye(type)
gui.mdiff = gui.mstress - gui.mrest; 
generatebullseye('diff')
save2set

%-------------------------------
function setsectionscore(type,region,sector, incr)
%-----------------------------
global DATA
    %Unpack vector to matrix
    
    gui=DATA.GUI.PerfusionScoring;

    switch sector
      case 'anterior'
        if region==2
          inds=[1, 2, 21, 22, 23, 24];
        else
          inds=1:4;
        end
        
      case 'anteroseptal'
        inds=5:8;
      
      case 'inferoseptal'
        inds=9:12;
        
      case 'inferior'
        if region==2
          inds = 9:14;
        else
          inds = 13:16;
        end
        
      case 'inferolateral'
        inds=17:20;
      
      case 'anterolateral'
        inds=21:24;

      case 'septal'
        inds=3:8;
      
      case 'lateral'
        inds=15:20;
      
      case 'apex'
        inds=1:24;
    end
    
    m = gui.(['m',type]);
    if m(inds(1),region)==1 && incr==1
      incr=2;
    elseif m(inds(1),region)==2 && incr==1
      incr=0;
    end
    
    switch type
      case 'rest'
        gui.mrest(inds,region)=incr;%gui.mrest(inds(1),region)+incr;
      case 'stress'
        gui.mstress(inds,region)=incr;%gui.mstress(inds(1),region)+incr;
      case 'scar'
        gui.mscar(inds,region)=incr;%gui.mscar(inds(1),region)+incr;
    end
    
 %---------------------
  function save2set
  %---------------------
  global DATA SET
  
  gui=DATA.GUI.PerfusionScoring;

  outdata=cell(3,17);
   for loop=1:17
    [stri,pos] = reportbullseye('aha17nameandpos',loop); %Get name and position of export
    outdata{1,loop+1}=stri;
   end
   
   counter=1;
   for type={'stress','rest','diff'}
     m=gui.(['m',type{1}]);
     counter=counter+1;
     outdata{counter,1}=type{1};
   outdata{counter,2} = mynanmean(m(1:4,4));
      outdata{counter,3} = mynanmean(m(5:8,4));
      outdata{counter,4} = mynanmean(m(9:12,4));
      outdata{counter,5} = mynanmean(m(13:16,4));
      outdata{counter,6} = mynanmean(m(17:20,4));
      outdata{counter,8} = mynanmean(m(21:24,4));
      %mid
      outdata{counter,8} = mynanmean(m(1:4,3));
      outdata{counter,9} = mynanmean(m(5:8,3));
      outdata{counter,10} = mynanmean(m(9:12,3));
      outdata{counter,11} = mynanmean(m(13:16,3));
      outdata{counter,12} = mynanmean(m(17:20,3));
      outdata{counter,13} = mynanmean(m(21:24,3));
      %apical
      outdata{counter,14} = mynanmean([m(1:2,2) ; m(21:24,2)]);
      outdata{counter,15} = mynanmean(m(3:8,2));
      outdata{counter,16} = mynanmean(m(9:14,2));
      outdata{counter,17} = mynanmean(m(15:20,2));
      %apex
      outdata{counter,18} = mynanmean(m(:,1));
   end
   for type = {'stress','rest','scar'}
     no=gui.([type{1},'no']);
     if ~isempty(no)
       SET(no).PerfusionScoring.export=outdata;
       SET(no).PerfusionScoring.(['m',type{1}])=gui.(['m',type{1}]);
       if any(strcmp(type{1},{'stress','rest'}))
         SET(no).PerfusionScoring.mdiff=gui.mdiff;
       end
       SET(no).PerfusionScoring.stressno=gui.stressno;
       SET(no).PerfusionScoring.restno=gui.restno;
       SET(no).PerfusionScoring.scarno=gui.scarno;
     end
   end
   
   
    %---------------------------
  function reset(type)
    %------------------------
global DATA 
gui=DATA.GUI.PerfusionScoring;
gui.(['m',type])= zeros(24,4);
generatebullseye(type);

gui.mdiff = gui.mstress - gui.mrest; 
generatebullseye('diff');
save2set
 
 %-----------
function generatebullseye(type,handle,no)
%-----------
global DATA 


%Set up
ahanumslices = 3;
n = 200;
scale = n/(ahanumslices+1);
numsectors = 24;

if nargin == 1;
gui=DATA.GUI.PerfusionScoring;
switch type
  case 'rest'
    m = gui.mrest;
    bullseyehandle=gui.handles.bullseyeaxesrest;
  case 'stress'
    m=gui.mstress;
    bullseyehandle=gui.handles.bullseyeaxesstress;
  case 'diff'
    m=gui.mdiff;
    bullseyehandle=gui.handles.diffaxes;
     case 'scar'
    m=gui.mscar;
    bullseyehandle=gui.handles.bullseyeaxesscar;   
end
gui.center=n+1;
else
  global SET
  try
  m=SET(no).PerfusionScoring.(['m',type]);
  bullseyehandle=handle;
  catch
    mydisp('Perfusion Scoring not available')
    return;
  end
end

set(bullseyehandle,'Color',[0 0 0],'Visible','on');
[x,y] = ndgrid(...
  linspace(-ahanumslices-1,ahanumslices+1,2*n+1),...
  linspace(-ahanumslices-1,ahanumslices+1,2*n+1));
rad = sqrt(x.*x+y.*y);

%Createidx outer
ang = angle(complex(y,x))+pi;
ang = numsectors*ang/(2*pi);
ang = mod(-(ang-numsectors/3),numsectors); %orient it iaccording to sectors in AHA 17-segment model
idxouter = 1+min(floor(ang),(numsectors-1))+(numsectors)*min(floor(rad),ahanumslices);

%Createidx inner
ang = mod(angle(complex(y,x))+pi+pi/4,2*pi);
ang = numsectors*ang/(2*pi);
ang = mod(-(ang-numsectors/3),numsectors); %orient it iaccording to sectors in AHA 17-segment model
idxinner = 1+min(floor(ang),(numsectors-1))+(numsectors)*min(floor(rad),ahanumslices);

idx = idxouter;
idx(rad<2) = idxinner(rad<2);

im = m(idx);
im(rad>(ahanumslices+1)) = NaN;
%im = rad;

%View data
alpha = double(not(isnan(im)));
im(isnan(im)) = 0;
h = imagesc(im,'parent',bullseyehandle);

set(h,'alphadata',alpha,'AlphaDataMapping','scaled');
axis(bullseyehandle,'image','off');
%cmap=[0.4259    0.2759    0.2759;   0.6045    0.3935    0.3935;    0.7412    0.4833    0.4833];
%a=pink;
%cmap = [a(90,:);a(140,:);a(190,:)];

%pink one
%cmap=[0.7412    0.4833    0.4833; 0.8356    0.7230    0.6040;0.9107    0.9107    0.7043];
cmap=[153/255 204/255 1;1 0.95 0;1 0 0];
%cmap=[0.0598    0.6841    0.7247; 0.9696    0.7300    0.2706; 1          0         0];
colormap(bullseyehandle,cmap);
caxis(bullseyehandle,[0 2])

if nargin == 1
  set(h,'ButtonDownFcn',sprintf('perfusion.perfusionscoring(''score'',''%s'')',type));
end

%Draw circles
om = linspace(0,2*pi,100);
xc = sin(om);
yc = cos(om);
hold(bullseyehandle,'on');
for loop=1:(ahanumslices+1)
  h = plot(bullseyehandle,n+1+scale*loop*xc,n+1+scale*loop*yc,'w-');
  set(h,'linewidth',2);
end;
hold(bullseyehandle,'off');

%Draw lines
hold(bullseyehandle,'on');
b = sqrt(0.75);
a = 0.5;
c = 1/sqrt(2);
h = plot(bullseyehandle,scale*[0 2],scale*[4 4],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[6 8],scale*[4 4],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4-c 4-2*c],scale*[4-c 4-2*c],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4+c 4+2*c],scale*[4+c 4+2*c],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4-c 4-2*c],scale*[4+c 4+2*c],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4+c 4+2*c],scale*[4-c 4-2*c],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4-4*a 4-2*a],scale*[4-4*b 4-2*b],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4-4*a 4-2*a],scale*[4+4*b 4+2*b],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4+4*a 4+2*a],scale*[4+4*b 4+2*b],'w-'); set(h,'linewidth',2);
h = plot(bullseyehandle,scale*[4+4*a 4+2*a],scale*[4-4*b 4-2*b],'w-'); set(h,'linewidth',2);
hold(bullseyehandle,'off');

%print out the strain values in the bullseye
rloops = [1 4 6 6];
for cloop=1:(ahanumslices+1)
  for rloop = 1:rloops(cloop)
    if cloop > 2 %basal and mid
      r = max(1,mod(round(rloop*length(xc)/rloops(cloop)+2/6*length(xc)),length(xc)));
    else
      r = max(1,mod(round(rloop*length(xc)/rloops(cloop)+1/2*length(xc)),length(xc)));
    end
    if rloops(cloop) == 6 %basal and mid
      mx = 1+(rloop-1)*4:4+(rloop-1)*4;
    elseif rloops(cloop) == 4 %apical
      if rloop == 4
        mx = [1:2 21:24];
      else
        mx = 3+(rloop-1)*6:8+(rloop-1)*6;
      end
    else %apex
      mx = 1:24;
    end
    value = (mynanmean(m(mx,cloop)));
    if ~isnan(value)
      text_handle=text(n+1+(scale/2*min(1,max(0,cloop-1))+scale*(cloop-1))*xc(r), ...
        n+1+(scale/2*min(1,max(0,cloop-1))+scale*(cloop-1))*yc(r), ...
        sprintf('%0.3g',value),'Parent',bullseyehandle,'HorizontalAlignment','center','fontweight','bold');
      
      if nargin == 1
        set(text_handle,'ButtonDownFcn',sprintf('perfusion.perfusionscoring(''score'',''%s'')',type));
      end
      
    end
  end
end;

switch type
  case 'stress'
    title(bullseyehandle,dprintf('Stress'),'FontSize',14)
  case 'rest'
    title(bullseyehandle,dprintf('Rest/LGE'),'FontSize',14)
  case 'diff'
    title(bullseyehandle,dprintf('Diff'),'FontSize',14)
end
% colorbar('peer',gui.handles.bullseyeaxes);
% maxvalue = max(max(m(:)),abs(min(m(:))));
% if isnan(maxvalue) || maxvalue==0
%   maxvalue = 1;
% end;
% set(gui.handles.bullseyeaxes,'clim',[-maxvalue maxvalue]);


%--------------------------------------------
function close_Callback
%--------------------------------------------
global DATA SET
try
gui=DATA.GUI.PerfusionScoring;
gui.play=0;
gui.stressplay=0;
gui.restplay=0;

for type = {'stress' 'rest'}
  SET(gui.([type{1},'no'])).PerfusionScoring.diffbullseye=frame2im(mygetframe(gui.handles.diffaxes));
end

for type = {'stress' 'rest' 'scar'}
  SET(gui.([type{1},'no'])).PerfusionScoring.([type{1}, 'bullseye'])=frame2im(mygetframe(gui.handles.(['bullseyeaxes',type{1}])));
end 

save2set

close(gui)
catch
close(gcf)  
end