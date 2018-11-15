%--------------------------------
function [varargout] = plugin_exporttofourflow(fcn,varargin)
%-------------------------------------------------
%Plugin for export flow to fourflow
if nargin==0
  myfailed('Need at least one input argument.');
end;

switch fcn
  case 'getname'
    varargout = cell(1,1);
    varargout{1} = 'Export to FourFlow';

    %Set the main "summarize" menu to not perform a callback.
    set(varargin{1},'Callback','');

    %Check if valid license, then display options
    %if isequal(3,getmodule(2,'3',[],true)) %'' => make license check silent.

    %Register submenus
    uimenu(varargin{1},'Label','Background correction...','Callback','plugin_exporttofourflow(''background'')');
    uimenu(varargin{1},'Label','Unwrap flow...','Callback','flowunwrap');
    uimenu(varargin{1},'Label','Remove outer "zero-border" (if needed)','Callback', 'plugin_exporttofourflow(''fixisoborders'')');
    uimenu(varargin{1},'Label','Duplicate heartcycle for 4D stack','Callback', 'plugin_exporttofourflow(''repeat4dcycles'')');
    uimenu(varargin{1},'Label','Export to FourFlow','Callback','fourFlowExportGUI');
  case 'getdependencies'
    %Here: List all depending files.
    varargout = cell(1,4);

    %M-files
    varargout{1} = {...
      'matlab2ensightgolddeformingmesh.m' ...
      'matlab2ensightgold_multiple_parts.m' ...
      'fourFlowExportGUI.m' ...
      'getangles.m' ...
      'matlab2vtk.m'};

    %Fig-files
    varargout{2} = {'fourFlowExportGUI.fig'};

    %Mat-files
    varargout{3} = {};
%     varargout{3} = {'4dstacknbrs.mat'};

    %Mex-files
    varargout{4} = {};
      
  otherwise
    macro_helper(fcn,varargin{:});
    [varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard
end;

%--------------------------------
function fixisoborders %#ok<DEFNU>
%--------------------------------
%
%Code which removes the high-gradient along the outer borders of some 4D
%dataset which otherwise result in poor visibility in FourFlow
%iso-surfaces
%
%
global SET

%Get 4D image stacks from user:
nom = find4dflowno();
nops = [SET(nom).Flow.PhaseX SET(nom).Flow.PhaseY SET(nom).Flow.PhaseNo];

if nnz(isnan([nom(1), nops]))>0
	mymsgbox('Invalid input: Exiting.')
	return
end

%set all zero values in magnitude image to 0.5 in velocity images:
nbrofvels = length(nops);
if max(nops)>length(SET)
	mymsgbox('Not enough stacks loaded: Exiting')
	return
end

for v = 1:nbrofvels
	SET(nops(v)).IM(SET(nom(1)).IM==0) = 0.5;
end

%--------------------------------
function repeat4dcycles %#ok<DEFNU>
%--------------------------------
%
%Implement this function first in this file. However, if it works fine,
%implement it in fourflowexportGUI instead so that the stacks does not get
%replicated each time. 
%
%NB: This code does not support ROI:s or ENDO/EPI contours in the 4D
%dataset. 
%

global SET

%Get 4d stack numbers: 
nom = find4dflowno();
nops = [SET(nom).Flow.PhaseX SET(nom).Flow.PhaseY SET(nom).Flow.PhaseNo];

nbrofvels = length(nops);
if max([nom, nops])>length(SET)
	mymsgbox('Invalid input: Maximum stack number exceeded. Exiting')
	return
end

%Get number of cycles from user: 
question = {'Select Nbr of heart cycles: '};
dlgtitle = 'Please select desired number of heartbeats: ';
default = {'2'};
userinput = inputdlg(question, dlgtitle, 1, default);
reps = round(str2double(userinput{1}));

if reps>1
	%Duplicate timeframes "reps" times in the entire 4D stack. 
	%magnitude stack:
	SET(nom).IM = repmat(SET(nom).IM, [1, 1, reps, 1]);
	SET(nom).TSize = reps*SET(nom).TSize;
	SET(nom).TimeVector = 0:SET(nom).TIncr:(SET(nom).TSize-1)*SET(nom).TIncr;
	
	%Velocity stacks:
	for v = 1:nbrofvels
		SET(nops(v)).IM = repmat(SET(nops(v)).IM, [1, 1, reps, 1]);
		SET(nops(v)).TSize = reps*SET(nops(v)).TSize;
		SET(nops(v)).TimeVector = 0:SET(nops(v)).TIncr:(SET(nops(v)).TSize-1)*SET(nops(v)).TIncr;
	end
	mymsgbox('Replications completed.');
else
	mymsgbox('1 heartbeat selected: No repetitions were performed. ');
	return
end


%----------------------
function background(varargin) %#ok<DEFNU>
%----------------------
global SET DATA

if nargin == 1
  no = varargin{1};
elseif nargin == 0
  stacks = cell(1);
  for noloop = 1:length(SET)
    stacks{noloop} = ...
      sprintf('%i: %s', noloop, strtrim(SET(noloop).SeriesDescription));
  end
  no = mymenu('Choose flow (2D or 4D)', stacks);
else
  no = 0;
end

if no < 1 || no > length(SET) || isempty(SET(no).Flow)
  myfailed('Error need to be flow data to do automatic background correction.');
  return;
end;

%Ensure that no is magnitude image
no = SET(no).Flow.MagnitudeNo;

%Extract
nx = SET(no).Flow.PhaseX;
ny = SET(no).Flow.PhaseY;
nz = SET(no).Flow.PhaseNo;

if isempty([nx ny nz])
  myfailed('Error need to be flow data to do automatic background correction.');
  return;
end;

%--- User input
dostore = yesno('Do you want to store background correction?', [], DATA.GUI.Segment);
%dostaticphantom = yesno('Are you using a phantom that is static in time?', [], DATA.GUI.Segment);

datatypes = {'in vivo', 'static phantom', 'vortex tank'};
typeind = mymenu('Type of data', datatypes);
datatype = datatypes{typeind};

dim = mymenu('Select order of background correction',...
  '0 constant',...
  '1 linear',....
  '2nd order',...
  '3rd order',....
  '4th order');

if isequal(dim,0)
  myfailed('Aborted.');
  return;
end;

dim = dim-1; %since first choice should return 0

%Ask settings
instruct = [];
instruct.ThresholdPercentage = 50;
instruct.VoxelEdgeCut = 2;
instruct.MagnitudeThresholdPercentage = 10;
[instruct,ok] = inputstruct(instruct,'Enter settings');

if ~ok
  myfailed('Aborted.');
  return;
end;

if SET(no).ZSize == 1
  if instruct.VoxelEdgeCut > 0.15*min([SET(no).XSize SET(no).YSize])
    myfailed('Too large edge cut. Aborted.');
    return;
  end
else
  if instruct.VoxelEdgeCut > 0.15*min([SET(no).XSize SET(no).YSize SET(no).ZSize])
    myfailed('Too large edge cut. Aborted.');
    return;
  end
end;
if instruct.MagnitudeThresholdPercentage>70
  myfailed('Too large magnitude threshold (max 70). Aborted.');
  return;
end;

%If dostore then 4 new image stacks are formed with mag,and phasecorrection.
%If dostaticphantom is used to do background correction for flow phantoms
%that are static in time.

%Prepare to loop over the three velocity components.
nop = {nx, ny, nz};

%--- Setup
sz = size(SET(no).IM);
if SET(no).ZSize == 1
  sz = [sz 1];
end

%Create a grid
h = waitbar(0,'Creating grid');
myworkon;
[x,y,z] = ndgrid(...
  sz(1)/2*linspace(-1,1,sz(1)),...
  sz(2)/2*linspace(-1,1,sz(2)),...
  sz(4)/2*linspace(-1,1,sz(4))); % NB Z==4
x = x(:);
y = y(:);
z = z(:);

switch dim
  case 0
    XV = ones(size(x));
  case 1
    XV = [ones(size(x)) ...
      x y z];
  case 2
    XV = [ones(size(x)) ...
      x y z ...
      x.^2 y.^2 z.^2 x.*y x.*z y.*z];
  case 3
    XV = [ones(size(x)) ...
      x y z ...
      x.^2 y.^2 z.^2 x.*y x.*z y.*z ...
      x.^3 y.^3 z.^3 x.^2.*y x.^2.*z y.^2.*x y.^2.*z z.^2.*x z.^2.*y x.*y.*z];
  case 4
    XV  = [ones(size(x)) ...
      x y z ...
      x.^2 y.^2 z.^2 x.*y x.*z y.*z ...
      x.^3 y.^3 z.^3 x.^2.*y x.^2.*z y.^2.*x y.^2.*z z.^2.*x z.^2.*y x.*y.*z ...
      x.^4 y.^4 z.^4 x.^3.*y x.^3.*z y.^3.*x y.^3.*z z.^3.*x z.^3.*y x.^2.*y.^2 x.^2.*z.^2 x.^2.*y.*z y.^2.*z.^2 y.^2.*x.*z z.^2.*x.*y ...
      ];
end;

switch datatype
  case {'in vivo', 'static phantom'}
    %--- Calculate weight matrix
    var_t = single(zeros([SET(no).XSize,SET(no).YSize 1 SET(no).ZSize]));
    for nloop=1:length(nop)
      if ~isempty(nop{nloop})
        var_t = var_t + var(SET(nop{nloop}).IM, [], 3); % variance along temporal dimension
      end
    end
    sd = sqrt(var_t);
    
    % Convert to physical velocities, normalize to 1.0 m/s.
    sd = sd*2*SET(no).VENC/100;

    %Cut away pixels with too low magnitude
    magthreshold = instruct.MagnitudeThresholdPercentage/100;
    sd(squeeze(mean(SET(no).IM,3))<magthreshold) = NaN;
    
    % Find threshold for standard deviation. Cut away with too large std
    threspercentage = instruct.ThresholdPercentage/100;
    sd(sd==0) = NaN; %Cut away pixels with exact zero.
    sdsort = sort(sd(:)); %NaN's are stored last upon sorting
    sdsort = sdsort(~isnan(sdsort));
    thres = sdsort(round(numel(sdsort)*threspercentage));
    
    sdw = 1./(sd.^2);
    sdw(sd>thres) = NaN;

  case 'vortex tank'
    % Use everything that's >5cm from centerline
    R0 = 75; % stationary fluid begins D mm from centerline
    xx = (1:SET(no).XSize)*SET(no).ResolutionX;
    yy = (1:SET(no).YSize)*SET(no).ResolutionY;
    zz = (1:SET(no).ZSize)*(SET(no).SliceThickness + SET(no).SliceGap);
    cx = SET(no).CenterX*SET(no).ResolutionX;
    cy = SET(no).CenterY*SET(no).ResolutionY;
    
    [XX, YY, ~, ~] = ndgrid(xx, yy, 1, zz);
    R = sqrt((XX-cx).^2 + (YY-cy).^2);
    
    sdw = ones([SET(no).XSize SET(no).YSize 1 SET(no).ZSize]);
    sdw(R < R0) = 0;
end

% Make sure single pixels are not weighted too high -
% cap at 95th percentile
sdw(isnan(sdw)) = 0; % Make NaN's as 0
sdwsort = sort(sdw(:));
sdwmax = sdwsort(round(numel(sdwsort)*0.95));
sdw = min(sdw,sdwmax);

sdw(isinf(sdw)) = NaN;

% Cut away edges in the volume - problems seem to arise with the weight map
% there.
if instruct.VoxelEdgeCut >=1
  e = instruct.VoxelEdgeCut;
  mask = zeros(size(sdw));

  if SET(no).ZSize == 1
    mask(...
      (1+e):(SET(no).XSize-e), ...
      (1+e):(SET(no).YSize-e), ...
      :) = 1;
  else
    mask(...
      (1+e):(SET(no).XSize-e), ...
      (1+e):(SET(no).YSize-e), ...
      (1+e):(SET(no).ZSize-e)) = 1;
  end
  
  sdw(~mask) = 0;
end

%If static over time then cut away what is stationary tissue.
if isequal(datatype, 'static phantom')
  velmag = repmat(single(0),[SET(no).XSize,SET(no).YSize 1 SET(no).ZSize]);
  for nloop=1:length(nop)
    velmag = velmag+max(abs(SET(nop(nloop)).IM-single(0.5)),[],3).^2;
  end;
  velmag = sqrt(velmag);
  sdw(velmag>0.04) = 0; %0.05 = 10% of VENC
  clear velmag;
end;

%Extract which ones to take
logind = (sdw~=0);
mydisp(dprintf('%d percent usage',round(100*sum(logind(:))/numel(sdw))));

close(h);
myworkoff;


% Display diagnostics
plotsl = round(size(sdw,4)/2);
tfr = round(size(sdw,3)/2);
im = SET(no).IM(:,:,tfr,plotsl);
wim = sdw(:,:,1,plotsl).*im;

figure(22);

imagesc(wim)
title('weight image')
axis equal
axis tight
colorbar;

if ~yesno('Do you want to proceed?', [], DATA.GUI.Segment);
  myworkoff;
  return;
end;

if dostore
  %If dostore create new image stacks
  k = length(SET);
  nw = k+1;

  SET(nw) = SET(no);
  
  nextno = nw+1;
  for ii = 1:length(nop)
    if isempty(nop{ii})
      continue
    end
    
    newnop{ii} = nextno; %#ok<AGROW>
    SET(newnop{ii}) = SET(nop{ii});
   
    nextno = nextno+1;
  end
end;

%Loop over the three velocity components.
h = waitbar(0,'Please wait background correction.');
savewarn = warning;
warning off; %#ok<WNOFF>

XVcut = XV(logind,:);

%Loop over time
condM = zeros(1, SET(no).TSize);
ndirections = ~isempty(nx) + ~isempty(ny) + ~isempty(nz);
steps = SET(no).TSize*ndirections;
stepsdone = 0;
for tloop=1:SET(no).TSize

  w = sdw.*SET(no).IM(:,:,tloop,:); % (sd weight)*(image mag)
 
  wcut = w(logind);
  wcut = wcut(:);

  M = XVcut.*repmat(wcut,1,size(XVcut,2));
  condM(tloop) = cond(M);

  if dostore
    SET(nw).IM(:,:,tloop,:) = w;
  end;

  % Loop over directions
  for nloop=1:3
    if isempty(nop{nloop})
      continue
    end
    
    v = SET(nop{nloop}).IM(:,:,tloop,:) - single(0.5);
    vcut = v(logind);

    % Compute coefficients A for fitting polynomial
    A = M \ (vcut(:).*wcut);

    % Compute correction and subtract
    C = XV*A;
    C = reshape(C,[SET(nop{nloop}).XSize SET(nop{nloop}).YSize 1 SET(nop{nloop}).ZSize]);
    SET(nop{nloop}).IM(:,:,tloop,:) = SET(nop{nloop}).IM(:,:,tloop,:) - C;

    if dostore
      SET(newnop{nloop}).IM(:,:,tloop,:) = C;
    end
    
    stepsdone = stepsdone+1;
  end
  
  waitbar(stepsdone/steps,h);
end;

close(h);
warning(savewarn);
myworkoff;

disp(sprintf('Fitting condition number: %g +- %2.2f%%', ...
  mean(condM), 100*std(condM)/mean(condM))); %#ok<DSPS>

%Additional storing and data manipulation
if dostore
  % For all images
  for nn = [nw newnop{1} newnop{2} newnop{3}]
    SET(nn).Flow.MagnitudeNo = nw;
    SET(nn).Flow.PhaseX = newnop{1};
    SET(nn).Flow.PhaseY = newnop{2};
    SET(nn).Flow.PhaseNo = newnop{3};
    SET(nn).Linked = [nw newnop{1} newnop{2} newnop{3}];
    
    SET(nn).IntensityMapping.Brightness = 0.5;
    SET(nn).IntensityMapping.Contrast = 1;
  end
  
  % Weight image
  SET(nw).Children = [newnop{1} newnop{2} newnop{3}];
  SET(nw).Parent = [];
  scalefactor = max(SET(nw).IM(:));
  SET(nw).IM = SET(nw).IM/scalefactor;
  SET(nw).IntensityScaling = scalefactor;
  SET(nw).IntensityOffset = 0;
  SET(nw).ImageType = sprintf('bgc order %i weight', dim);

  % Phase images
  maxv = 0;
  for ii = 1:length(newnop)
    if isempty(newnop{ii})
      continue
    end
    
    maxv_here = max(abs(SET(newnop{ii}).IM(:))); % 0.5 not yet added
    maxv = max(maxv, maxv_here);
    SET(newnop{ii}).Children = [];
    SET(newnop{ii}).Parent = nw;
  end
  
  % Rescale VENC for display
  R = 5;
  newVENC = R*ceil(maxv*2*SET(no).VENC/R);
  imtype = { ...
    sprintf('bgc order %i (Up/Down), scale %i cm/s', dim, newVENC), ...
    sprintf('bgc order %i (Left/right), scale %i cm/s', dim, newVENC), ...
    sprintf('bgc order %i (Through plane), scale %i cm/s', dim, newVENC)};
  
  for ii = 1:length(newnop)
    if ~isempty(newnop{ii})
      SET(newnop{ii}).IM = SET(newnop{ii}).IM*SET(no).VENC/newVENC + single(0.5);
      SET(newnop{ii}).VENC = newVENC;
      SET(newnop{ii}).ImageType = imtype{ii};
    end
  end
  SET(nw).VENC = newVENC;
end % if dostore

% Mark flow stacks as corrected
for noloop = [no nop{1} nop{2} nop{3}]
    SET(noloop).SeriesDescription = ...
    sprintf('%s, bgc order %i', SET(noloop).SeriesDescription, dim);
end

drawfunctions('drawthumbnails');


%--------------------------------
function result = numfromstr(str)
%--------------------------------
%Extract number from string

%replace period (.) with -
str = regexprep(str,'\.','-');
pat = '\d+-\d+|\d+';
result = regexp(str,pat,'match');
result = regexprep(result,'-','\.');


%--------------------------------
function flowno = find4dflowno
%--------------------------------
% Find 4D flow stack, ask the user if there are several.
%

global SET DATA

flowno = find4dflowstacks;

if isempty(flowno)
  myfailed('No 3D flow image stacks found',DATA.GUI.Segment);
  return;
end;

%Ask in a menu if more than one 4D flow
if length(flowno)>1
  names = cell(1,length(flowno));
  for loop=1:length(flowno)
    names{loop} = sprintf('Stack %d: %s',flowno(loop),SET(flowno(loop)).ImageType);
  end;
  m = mymenu('Choose 3D image stack to take from',names{:});
  if isequal(m,0)
    myfailed('Aborted.');
    return;
  end;
  flowno = flowno(m);
end;


%--------------------------------
function flownos = find4dflowstacks
%--------------------------------
% Find 4D flow stacks
%

global SET
flownos = [];

for nloop = 1:length(SET)
  if not(isempty(SET(nloop).Flow)) && ...
      not(isempty(SET(nloop).Flow.PhaseNo)) && ...
      not(isempty(SET(nloop).Flow.PhaseX)) && ...
      not(isempty(SET(nloop).Flow.PhaseY)) && ...
      isequal(SET(nloop).Flow.MagnitudeNo,nloop)
    flownos = [flownos nloop]; %#ok<AGROW>
  end;
end


%--------------------------------------------
function outim = resampletimeframes(z,frames)
%--------------------------------------------
%Simple so far...

%--- Resample the newim to correct number of timeframes
if not(isequal(size(z,3),frames))
  if isequal(size(z,3),1)
    outim = repmat(z,[1 1 frames]);
  else
    %--- Time resolved clause

    %Resample image to correct number of frames,.
    intpos = linspace(1,size(z,3)-1,frames);
    intleft = floor(intpos);
    intright = floor(intpos)+1;

    %Set up space for output data
    if isa(z,'double')
      outim = zeros([size(z,1) size(z,2) frames size(z,4)]);
    else
      outim = repmat(single(0),[size(z,1) size(z,2) frames size(z,4)]);
    end;
    switch ndims(z)
      case 3
        for tloop=1:length(intpos)
          outim(:,:,tloop) = ...
            (1-intpos(tloop)+intleft(tloop))*z(:,:,intleft(tloop))+...
            (1-intright(tloop)+intpos(tloop))*z(:,:,intright(tloop));
        end;
      case 4
        for tloop=1:length(intpos)
          outim(:,:,tloop,:) = ...
            (1-intpos(tloop)+intleft(tloop))*z(:,:,intleft(tloop),:)+...
            (1-intright(tloop)+intpos(tloop))*z(:,:,intright(tloop),:);
        end;
      otherwise
        myfailed('Upsampling in time not defined for this dimensionality.');
    end;
  end; %Timeresolved clause
else
  %Copy
  outim = z;
end;