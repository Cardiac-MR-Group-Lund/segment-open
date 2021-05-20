function [varargout] = longaxistools(varargin)
%Utility for functionality requested in Doppler-CIP project

if nargin == 0
  calcbiplanevolume
else
  macro_helper(varargin{:});
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
end

%--------------------------------------------
function [didcalc,outnos] = calcbiplanevolume
%--------------------------------------------
global SET

didcalc = false;
outnos = [];
fields = {'Endo','Epi','RVEndo','RVEpi'};
for segloop = 1:length(fields)
  segfield = fields{segloop};
  xfield = [segfield 'X'];
  volfield = '';
  switch segfield
    case 'Endo'
      volfield = 'LVV';
    case 'Epi'
      volfield = 'EPV';
    case 'RVEndo'
      volfield = 'RVV';
    case 'RVEpi'
      volfield = 'RVEPV';
  end
  
  ch4no = find(strcmp({SET.ImageViewPlane},'4CH'),1);
  if ~isempty(ch4no) && ~isempty(SET(ch4no).(xfield))
    seg4tfs = find(~isnan(SET(ch4no).(xfield)(1,:)));
    if ~isempty(seg4tfs)
      seg4tfs(seg4tfs == SET(ch4no).EST) = 0;
      seg4tfs(seg4tfs == SET(ch4no).EDT) = -1;
      seg4tfs(seg4tfs>0) = [];
      if numel(seg4tfs)<2
          seg4tfs = [];
      end
    end
  else
    seg4tfs = [];
    ch4no = [];
  end
  
  ch3no = find(strcmp({SET.ImageViewPlane},'3CH'),1);
  if ~isempty(ch3no) && ~isempty(SET(ch3no).(xfield))
    seg3tfs = find(~isnan(SET(ch3no).(xfield)(1,:)));
    if ~isempty(seg3tfs)
      seg3tfs(seg3tfs == SET(ch3no).EST) = 0;
      seg3tfs(seg3tfs == SET(ch3no).EDT) = -1;
      seg3tfs(seg3tfs>0) = [];
      if numel(seg3tfs)<2
          seg3tfs = [];
      end
    end
  else
    seg3tfs = [];
    ch3no = [];
  end
  
  ch2no = find(strcmp({SET.ImageViewPlane},'2CH'),1);
  if ~isempty(ch2no) && ~isempty(SET(ch2no).(xfield))
    seg2tfs = find(~isnan(SET(ch2no).(xfield)(1,:)));
    if ~isempty(seg2tfs)
      seg2tfs(seg2tfs == SET(ch2no).EST) = 0;
      seg2tfs(seg2tfs == SET(ch2no).EDT) = -1;
      seg2tfs(seg2tfs>0) = [];
      if numel(seg2tfs)<2
          seg2tfs = [];
      end
    end
  else
    seg2tfs = [];
    ch2no = [];
  end
  
  nos = [ch4no ch3no ch2no];
  multisliceinds = [SET(nos).ZSize] > 1;
  for no = nos(multisliceinds)
    fprintf(['%s image contains multiple slices. Remove slices to ' ...
      'include in calculation\n'],SET(no).ImageViewPlane);
  end
  nos = nos(~multisliceinds);
  outnos = union(outnos,nos);
  
  if ~isempty(volfield)
    for i = 1:length(nos)
      no = nos(i);
      SET(no).(volfield) = nan(1,SET(no).TSize);
    end
  end
  
  intsec = intersect(seg4tfs,intersect(seg3tfs,seg2tfs));
  if ~isempty(intsec)
    for tf = intsec
      [~,didcalc] = calcbiplanevolume_helper(nos,tf,segfield);
    end
  return  
  end
    
  for tf = setdiff(intersect(seg2tfs,seg3tfs),seg4tfs)
    [~,didcalc] = calcbiplanevolume_helper([ch2no ch3no],tf,segfield);
  end
  for tf = setdiff(intersect(seg3tfs,seg4tfs),seg2tfs)
    [~,didcalc] = calcbiplanevolume_helper([ch3no ch4no],tf,segfield);
  end
  for tf = setdiff(intersect(seg4tfs,seg2tfs),seg3tfs)
    [~,didcalc] = calcbiplanevolume_helper([ch4no ch2no],tf,segfield);
  end
  
end

%------------------------------------------------------------------------
function [volume,stus] = calcbiplanevolume_helper(nos,timeframe,segfield)
%------------------------------------------------------------------------
%Calculates area and biplane volume

global SET
stus = 2; %nbr of image stacks used to calculate, default 2

if numel(timeframe) > 1
  tfs = timeframe;
elseif isempty(timeframe)
  volume = nan;
  return
else
  if timeframe == 0
    tfs = [SET(nos).EST];
  elseif timeframe == -1
    tfs = [SET(nos).EDT];
  else
    tfs = timeframe*ones(size(nos));
  end
end

xfield = [segfield 'X'];
yfield = [segfield 'Y'];
volfield = '';
switch segfield
  case 'Endo'
    volfield = 'LVV';
  case 'Epi'
    volfield = 'EPV';
  case 'RVEndo'
    volfield = 'RVV';
  case 'RVEpi'
    volfield = 'RVEPV';
end

volume=zeros(length(nos),1);
if length(nos) == 2
  for noloop = 1:length(nos)
    no = nos(noloop);
    timeframe = tfs(noloop);
    [ix,iy] = calcfunctions('calcplaneintersections',no,setdiff(nos,no));
    if isempty(ix) || isempty(iy)
      disp('Image planes do not intersect, cannot calculate.')
      volume = nan; %fix for #1333 was volume(no) = nan; 
      stus = false;
      return
    end
    X = [SET(no).ResolutionX*(SET(no).(xfield)(:,timeframe)-ix(1)) ...
      SET(no).ResolutionY*(SET(no).(yfield)(:,timeframe)-iy(1))]';
  
    %do extra check that a horse shoe is drawn
    %extract convhull
    inds = convhull(X');
    a1 = calcfunctions('stablepolyarea',X(1,:),X(2,:));
    a2 = calcfunctions('stablepolyarea',X(1,inds),X(2,inds));
    if a2>a1*(1.5)
      % jelena: skip calculation for now
      volume = nan; 
      stus = false;
      return
      %this is likely a horseshoe segmentation
      [endoindexes, epiindexes] = splithorseshowsegmentation(X, inds); 
       X = X(:,endoindexes);
    end
  
    ivec = [ix(2)-ix(1);iy(2)-iy(1)];
    ivec = ivec/norm(ivec);
    S = [ivec ([0 1;-1 0]*ivec)];
    Xhat = S*X;
    x = Xhat(1,:);
    r = Xhat(2,:);
    slicearea = pi/4*sign(r).*r.^2;
    volume(noloop) = (1/1000)*stablepolyarea(x,slicearea);
  end
  volume = sum(volume);
elseif length(nos) > 2
  for noloop = 1:length(nos)
    no = nos(noloop);
    noix = nos ~= no;
    newnos = nos(noix);
    newtfs = tfs(noix);
    volume(noloop) = calcbiplanevolume_helper(newnos,newtfs,segfield);
  end
  if ~any(isnan(volume))
    stus = 3;
  end
  volume = mynanmean(volume);
else
  myfailed('Need to draw contours in at least two images');
  return
end

%Store it
if ~isempty(volfield)
  for noloop = 1:length(nos)
    no = nos(noloop);
    timeframe = tfs(noloop);
    SET(no).(volfield)(timeframe) = volume;
  end
end

%--------------------------------------
function [endoindexes, epiindexes] = splithorseshowsegmentation(originalcontour, epiindexes)
% find jumps in indices larger than 10% of the contour length

minjump = floor(length(originalcontour)/10);
maxjump = length(originalcontour) - 10;
jumps = find (abs(diff(epiindexes)) >= minjump & abs(diff(epiindexes)) <= maxjump);
% loop over all jumps
contourlength = zeros(length(jumps),1);
for loop = 1: length(jumps)
  % calculate contour length between points
  startindex = epiindexes(jumps(loop));
  endindex   = epiindexes(jumps(loop) + 1);
  if startindex > endindex
    tmp = startindex;
    startindex = endindex;
    endindex = tmp;
  end
  x = originalcontour(1,startindex : endindex); 
  y = originalcontour(2,startindex : endindex);
  contourlength(loop) = sum(sqrt(diff(x).^2 + diff(y).^2)); 
end
[epilength,index] = max(contourlength);
startindex = epiindexes(jumps(index));
endindex   = epiindexes(jumps(index) + 1);
if startindex > endindex
  tmp = startindex;
  startindex = endindex;
  endindex = tmp;
end
epiindexes = startindex:endindex; 
endoindexes = setdiff(1:length(originalcontour),epiindexes);
x = originalcontour(1,endoindexes); 
y = originalcontour(2,endoindexes);
endolength = sum(sqrt(diff(x).^2 + diff(y).^2));

if endolength > epilength
  %switch
  tmp = epiindexes;
  epiindexes = endoindexes;
  endoindexes = tmp;
end



%--------------------------------------
function showpointinallviews(clickedim) %#ok<DEFNU>
%--------------------------------------
global DATA SET NO
segment_main('normal_Buttondown',clickedim);
switch DATA.ViewPanelsType{clickedim}
  case 'hla'
    [y,z,x] = segment('getclickedcoords');
  case 'vla'
    [x,z,y] = segment('getclickedcoords');
  case 'gla'
    [ytemp,xtemp] = segment('getclickedcoords');
    [x,y,z] = calcfunctions('gla2sax',xtemp,ytemp,NO);
  otherwise
    [y,x,z] = segment('getclickedcoords');
end
h = DATA.Handles.imageaxes(DATA.ViewPanels==NO);
drawmode = 'yo';
doline = false;
if SET(NO).ZSize == 1
  doline = true;
  linx = [x x];
  liny = [y y];
  linz = [-300 300];
  lindrawmode = 'y:';
  linpos = calcfunctions('xyz2rlapfh',NO,linx,liny,linz);
end
for i = 1:numel(h)
  hloop = h(i);
  hold(hloop,'on');
  switch DATA.ViewPanelsType{i}
    case 'hla'
      plot(hloop,y,z,'yo');
    case 'vla'
      plot(hloop,x,z,'yo');
    case 'gla'
      [xo,yo] = calcfunctions('sax2gla',x,y,z,NO);
      plot(hloop,yo,xo,'yo');
    case {'one','orth'}
      plot(hloop,y,x,'yo');
    otherwise
      [xofs,yofs] = calcfunctions('calcoffset',z,[],NO);
      plot(hloop,y+yofs,x+xofs,'yo');
  end
end
pos = calcfunctions('xyz2rlapfh',NO,x,y,z);

for imloop = find(DATA.ViewPanels ~= NO & DATA.ViewPanels > 0)
  no = DATA.ViewPanels(imloop);
  h = DATA.Handles.imageaxes(imloop);
  hold(h,'on');
  xyz = calcfunctions('rlapfh2xyz',no,pos(:,1),pos(:,2),pos(:,3));
  plot(h,xyz(2,:),xyz(1,:),drawmode);
  if doline
    linxyz = calcfunctions('rlapfh2xyz',no,linpos(:,1),linpos(:,2),linpos(:,3));
    plot(h,linxyz(2,:),linxyz(1,:),lindrawmode);
  end
end

%Reset to normal buttondown function
% for loop = 1:length(DATA.Handles.imagehandle)
%   set(DATA.Handles.imagehandle(loop),'ButtonDownFcn',...
%     sprintf('segment_main(''normal_Buttondown'',%d)',loop));
% end



%-----------------------------------
function cros = dointersect(no1,no2)
%-------------------------------
%Checks if image stacks numbered nos intersect
[ix,iy] = calcfunctions('calcplaneintersections',no1,no2);
cros = ~isempty(ix) && ~isempty(iy);