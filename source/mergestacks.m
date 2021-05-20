function [varargout] = mergestacks(varargin)
%MERGESTACKS
%Function to merge several image stacks into one, provided they have the
%same image shape. Avoids slice duplication and stores the merged stack in
%a new SET element

%Nils Lundahl
%Minor bugfixing by Einar Heiberg


global SET DATA

if nargin < 1
  %Launch GUI
  gui = mygui('mergestacks.fig');
  DATA.GUI.MergeStacks = gui;
  handles = gui.handles;
  stackc = cell(1,numel(SET));
  for no = 1:numel(SET)
    stackc{no} = sprintf('%d. %s, %s',no,SET(no).ImageType,SET(no).ImageViewPlane);
  end
  set(handles.imagestackslistbox,'String',stackc);
elseif isnumeric(varargin{1})
  domerge(varargin{:})
else
  arg = varargin{1};
  handles = guihandles(gcbf);
  %to enable multiple selection, set max = 2.0
  nos = mygetvalue(handles.imagestackslistbox);
  if strcmp(arg,'setfactor_Callback')
    setfactor_Callback;
  elseif strcmp(arg,'update')
    kids = [SET(nos).Children];
    set(handles.imagestackslistbox,'Value',union(nos,kids));
  elseif strcmp(arg,'domerge')
    if numel(nos) < 2
      myfailed('Need to select at least two image stacks');
      return
    end
    
    %do timeframe check
    Ts=zeros(1,length(nos));
    
    counter=1;
    for no=nos
      Ts(counter)=SET(no).TSize;
      counter=counter+1;
    end
    
    
    if ~all(Ts(1)==Ts)
      answer=yesno('Number of timeframes does not match. Do you wish to resample stacks?');
      if answer
        initwhichstack_Callback(nos)
%       if ~DATA.GUI.whichstack.proceed
%         return
%       end
%       
      else
        return
      end
      
    end
    
    parents = unique([SET(nos).Parent]); %Bugfix, EH: Was without (nos) here.
    if numel(parents) > 1
      domerge(parents);
      parno = numel(SET);
      SET(parno).Parent = [];
      SET(parno).Children = parno+(1:numel(SET(parents(1)).Linked)-1);
      SET(parno).Linked = [parno SET(parno).Children];
      for lix = 1:numel(SET(parents(1)).Linked)-1
        mergenos = [];
        for pno = parents
          linkies = setdiff(SET(pno).Linked,pno);
          mergenos = [mergenos linkies(lix)]; %#ok<AGROW>
        end
        domerge(mergenos);
        tempno = numel(SET);
        SET(tempno).Parent = parno;
        SET(tempno).Children = [];
        SET(tempno).Linked = SET(parno).Linked;
      end
    else
      domerge(nos);
    end
    viewfunctions('setview',1,1,numel(SET),{'one'})%segment('switchtoimagestack',numel(SET),true);
    flushlog;
    close(handles.figure1);
  end
end

%---------------------------------
function domerge(nostomerge,force)
%---------------------------------
global SET DATA

%Check number of nos before, make sure to only generate one new stack
nosbefore = numel(SET);

if nargin < 2
  force = false;
end

sliceloc = zeros(size(nostomerge));
for i = 1:numel(nostomerge)
  no = nostomerge(i);
  zdir = cross(...
    SET(no).ImageOrientation(1:3),...
    SET(no).ImageOrientation(4:6));
  sliceloc(i) = dot(zdir,SET(no).ImagePosition);
end
[~,ord] = sort(sliceloc,'descend');
nostomerge = nostomerge(ord);

slicegap=zeros(1,length(nostomerge)-1);
counter=1;

while numel(nostomerge) > 1
  no1 = nostomerge(1);
  no2 = nostomerge(2);
  %Check ResolutionX and ResolutionY
  xres = abs(SET(no1).ResolutionX-SET(no2).ResolutionX)/...
    (0.5*(SET(no1).ResolutionX+SET(no2).ResolutionX));
  yres = abs(SET(no1).ResolutionY-SET(no2).ResolutionY)/...
    (0.5*(SET(no1).ResolutionY+SET(no2).ResolutionY));
  if (xres>1e-6) || (yres>1e-6)
    area1 = SET(no1).ResolutionX*SET(no1).ResolutionY;
    area2 = SET(no2).ResolutionX*SET(no2).ResolutionY;
    if ~yesno(dprintf('Not the same ResolutionX or Resolution Y, difference in percentage of area in pixels is: %0.5g %%. This will introduce and error in the quantification, proceed?',...
        100*min([area1 area2])/max([area1 area2])))
      return;
    end
  end
  
  thesenos = [no1 no2];
  if all([SET(thesenos).ZSize] > 1)
    slicedist = [SET(thesenos).SliceGap] + [SET(thesenos).SliceThickness];
    if std(slicedist) > 1e-4 %Was e-5 changed to e-4 (EH)
      myfailed('Slice distance does not match.');
      return
    end
  elseif SET(no1).ZSize == 1 && SET(no2).ZSize > 1
    SET(no1).SliceThickness = SET(no2).SliceThickness;
  end
  
  %Image orientation
  zdir1 = cross(...
    SET(no1).ImageOrientation(1:3),...
    SET(no1).ImageOrientation(4:6));
  zdir2 = cross(...
    SET(no2).ImageOrientation(1:3),...
    SET(no2).ImageOrientation(4:6));
  angl = acos(zdir1*zdir2');
  if angl > 1e-2 && ~yesno(dprintf('Image orientations do not match. Angle offset is %0.2f degrees. Proceed?',angl*180/pi))
    return
  end
  
  % Extract and check image sizes, calculate new ZSize value
  orsz1 = size(SET(no1).IM);
  orsz1 = [orsz1 1 1]; %add trailing ones if t=1 or z=1;
  orsz1 = orsz1(1:4); %Take only first four
  xm1 = ceil(orsz1(1)/2);
  ym1 = ceil(orsz1(2)/2);
  zsz1 = orsz1(4);
  
  orsz2 = size(SET(no2).IM);
  orsz2 = [orsz2 1 1]; %add trailing ones if t=1 or z=1;
  orsz2 = orsz2(1:4); %Take only first four
  xm2 = ceil(orsz2(1)/2);
  ym2 = ceil(orsz2(2)/2);
  zsz2 = orsz2(4);
  newzsz = zsz1+zsz2;
  
  if ~isequal(orsz1(3),orsz2(3))
        myfailed('Number of timeframes does not match.')
        return
  end
  
  %Check slice order by looking at fh coordinates of end slices
  if ~force
    overlap = true;
    start1 = 1;
    start2 = 1;
    
    mostbasal=1e6;
    mostapical=-1e6;
    while overlap
      sortord = '';
      p1_1 = calcfunctions('xyz2rlapfh',no1,xm1,ym1,start1);
      if zsz1 > 1
        p1_2 = calcfunctions('xyz2rlapfh',no1,xm1,ym1,zsz1);
      else
        p1_2 = [];
        sortord = 'x';
      end
      p2_1 = calcfunctions('xyz2rlapfh',no2,xm2,ym2,start2);
      if zsz2 > 1
        p2_2 = calcfunctions('xyz2rlapfh',no2,xm2,ym2,zsz2);
      else
        p2_2 = [];
        sortord = [sortord 'y'];
      end
      
      pmat = [p1_1;p1_2;p2_1;p2_2]';

      zvec = zdir1*pmat;
      mostbasal=min([zvec,mostbasal]);
      mostapical=max([zvec,mostapical]);
      
      [~,sortind] = sort(zvec);
      sortord = [sortord char(sortind+'a'-1)];
      switch sortord
        case {'abcd','dcba','xabc','yabc','xcba','ycba','xyab','xyba'} %everything is fine
          noswitch = false;
          overlap = false;
        case {'acbd','dbca','xbac','xcab','yacb','ybca'} %overlap
          noswitch = false;
          overlap = true;
          start2 = start2+1;
        case {'cdab','badc','xacb','xbca','ycab','ybac'} %no overlap, but stacks in wrong order
          noswitch = true;
          overlap = false;
        case {'bdac','cadb'} %wrong order and overlap
          noswitch = true;
          overlap = true;
          start1 = start1+1;
        otherwise
          myfailed('Orientations of slice order do not match');
          return
      end
    end
    
    iters = max(start1,start2)-1;
    if iters > 0
      mywarning(dprintf(['%d slices appear to overlap. Consider ' ...
        'removing slices and reperforming this merge operation.'],iters));
    end
           
    if noswitch
      tmp = no2;
      no2 = no1;
      no1 = tmp;
    end
  end
  
  %Create new SET struct and get indices to fill
  setstruct = SET(no1);
  set2 = SET(no2);
  zpre = setstruct.ZSize;
  xext = 1:size(set2.IM,1);
  yext = 1:size(set2.IM,2);
  zext = zpre+1:newzsz;
  
  %estimation of slice gap    
  %If all merged stacks are single slice then we need to estimate
  %slicegap since zero in original images

  if newzsz>1
    slicegap(counter)=abs(mostapical-mostbasal)/(newzsz-1)-setstruct.SliceThickness;
    counter=counter+1;
    %setstruct.SliceGap=abs(mostapical-mostbasal)/(newzsz-1)-setstruct.SliceThickness;
    
  end
  
  %--- Update variables
  
  %Convert to true data and then back again to get correct scaling
  if isa(setstruct.IM,'single')
    setstruct.IM = calcfunctions('calctruedata',setstruct.IM,no1);
    setstruct.IM(xext,yext,:,zext) = calcfunctions('calctruedata',set2.IM,no2);
    setstruct.IM = calcfunctions('truedata2im',setstruct.IM,no1);
  else
    %If not single then just add, assume CT or true count values.
    setstruct.IM(xext,yext,:,zext) = set2.IM;
  end
  setstruct.ZSize = newzsz;
  setstruct.XSize = max(SET(no1).XSize,set2.XSize);
  setstruct.YSize = max(SET(no1).YSize,set2.YSize);
  setstruct.StartSlice = 1;
  setstruct.EndSlice = newzsz;

  %Scar
  if isempty(setstruct.Scar) && ~isempty(set2.Scar)
    %no1 is empty but no2 is not
    viability('viabilityreset_Callback','weighted',no1);
    setstruct.Scar = SET(no1).Scar;
  end
  if ~isempty(setstruct.Scar) && isempty(set2.Scar)
    %no1 have scar but no2 not
    viability('viabilityreset_Callback','weighted',no2);
    setstruct.Scar = SET(no2).Scar;
  end
  
  if not(isempty(set2.Scar))
    %Assign if exist in set2
    setstruct.Scar.IM(xext,yext,zext) = set2.Scar.IM;
    setstruct.Scar.Auto(xext,yext,zext) = set2.Scar.Auto;
    setstruct.Scar.Result(xext,yext,zext) = set2.Scar.Result;
    setstruct.Scar.Manual(xext,yext,zext) = set2.Scar.Manual;
    %SET(no).Scar.Undo = SET(no).Scar.Undo(:,:,ind);
    setstruct.Scar.MyocardMask(xext,yext,zext) = set2.Scar.MyocardMask;
    setstruct.Scar.NoReflow(xext,yext,zext) = set2.Scar.NoReflow;
  end
  
  %MaR
  if isempty(setstruct.MaR) && ~isempty(set2.MaR)
    %no1 is empty but no2 is not
    mar('initdefault',no1);
    setstruct.MaR = SET(no1).MaR;
  end
  if ~isempty(setstruct.MaR) && isempty(set2.MaR)
    %no1 have mar but no2 not
    mar('initdefault',no2);
    setstruct.MaR = SET(no2).MaR;
  end
  
  if not(isempty(setstruct.MaR))
    setstruct.MaR.Auto(xext,yext,:,zext) = set2.MaR.Auto;
    setstruct.MaR.Result(xext,yext,:,zext) = set2.MaR.Result;
    setstruct.MaR.Manual(xext,yext,:,zext) = set2.MaR.Manual;
    setstruct.MaR.MyocardMask(xext,yext,:,zext) = set2.MaR.MyocardMask;
  end
  
  %Merge contours
  contfields = {'Endo','Epi','RVEndo','RVEpi',...
    'EndoPin','EpiPin','RVEndoPin','RVEpiPin',...
    'EndoInterp','EpiInterp','RVEndoInterp','RVEpiInterp',};
  for floop = 1:numel(contfields)
    xprop = [contfields{floop} 'X'];
    yprop = [contfields{floop} 'Y'];
    cond = (~isempty(setstruct.(xprop))) + 2*(~isempty(set2.(xprop)));
    if ismember(contfields{floop},{'Endo','Epi','RVEndo','RVEpi'})
      %Field contains a matrix
      switch cond
        case 0 %property is empty in both stacks
          %do nothing
        case 1 %empty only in set2
          setstruct.(xprop)(:,:,zext) = nan;
          setstruct.(yprop)(:,:,zext) = nan;
        case 2 %empty only in setstruct
          setstruct.(xprop)(:,:,zext) = set2.(xprop);
          setstruct.(yprop)(:,:,zext) = set2.(yprop);
          setstruct.(xprop)(:,:,1:zpre) = nan;
          setstruct.(yprop)(:,:,1:zpre) = nan;
        case 3 %non empty
          setstruct.(xprop)(:,:,zext) = set2.(xprop);
          setstruct.(yprop)(:,:,zext) = set2.(yprop);
      end
    else
      %Field contains a cell
      switch cond
        case 0 %property is empty in both stacks
          %do nothing
        case 1 %empty only in set2
          setstruct.(xprop)(:,zext) = {[]};
          setstruct.(yprop)(:,zext) = {[]};
        case 2 %empty only in setstruct
          setstruct.(xprop) = cell(set2.TSize,newzsz);
          setstruct.(yprop) = cell(set2.TSize,newzsz);
          setstruct.(xprop)(:,zext) = set2.(xprop);
          setstruct.(yprop)(:,zext) = set2.(yprop);
        case 3 %non empty
          setstruct.(xprop)(:,zext) = set2.(xprop);
          setstruct.(yprop)(:,zext) = set2.(yprop);
      end
    end
  end
  
  setstruct.EndoDraged(:,zext) = set2.EndoDraged;
  setstruct.EpiDraged(:,zext) = set2.EpiDraged;
  
  % ROIs
  for loop=1:set2.RoiN
    set2.Roi(loop).Z = zext(set2.Roi(loop).Z);
    setstruct.Roi = [setstruct.Roi set2.Roi(loop)];
  end
  setstruct.RoiN = setstruct.RoiN + set2.RoiN;
  if setstruct.RoiN>0
    setstruct.RoiCurrent = 1;
  else
    setstruct.RoiCurrent = [];
  end
  
  % Measure
  for loop=1:numel(set2.Measure)
    set2.Measure(loop).Z = zext(set2.Measure(loop).Z);
  end
  setstruct.Measure = [setstruct.Measure set2.Measure];
  
  %Set no of generated image stack and remove links
  newno = nosbefore+1;
  setstruct.Parent = [];
  setstruct.Children = [];
  SET(newno) = setstruct;
  SET(newno).Linked = newno;
  nostomerge = [newno nostomerge(3:end)];
end

%use last slicegap calculation as slicegap.
SET(newno).SliceGap=slicegap(end);

if max(diff(slicegap))>1
  ok = yesno('More than one millimeter difference in spacing of slices. Proceed any way?','Not Equidistant slices',DATA.GUI.Segment);
  if ~ok
    SET(newno) = [];
    return;  
  end
end

%remove Strain tagging analysis
if not(isempty(SET(newno).StrainTagging))
  disp('Merging image stacks clears the Strain tagging analysis');
  if ~isempty(DATA.GUI.StrainTagging)
    straintagging.straintagging('close_Callback');
  end
  SET(newno).StrainTagging = [];
end

drawfunctions('drawthumbnails');

close_Callback; 


%--------------------------------------------------------------------
function initwhichstack_Callback(nostomerge)
%--------------------------------------------------------------------

global DATA SET
%Launch GUI
gui = mygui('whichstack.fig');
DATA.GUI.WhichStack = gui;
handles = gui.handles;
%DATA.GUI.whichstack.handles=handles;
handles.nostomerge=nostomerge;
stackc = cell(1,numel(nostomerge));
counter=1;
for no = nostomerge
  stackc{counter} = sprintf('%d. Number of timeframes: %d.',no,SET(no).TSize);
  counter=counter+1;
end
set(handles.imagestackslistbox,'String',stackc);
%alters if we should
%DATA.GUI.whichstack.proceed=0; 
uiwait(gcf)

%--------------------------------------------------------------------
function setfactor_Callback
%--------------------------------------------------------------------
global SET DATA NO

handles=guihandles(gcbf);
gui=DATA.GUI.whichstack;
ind=get(handles.imagestackslistbox,'value');
nom=SET(gui.nostomerge(ind)).TSize;
oldNO=NO;
for no=gui.nostomerge
  NO=no;
  tools('upsampletemporal_Callback',(nom/SET(no).TSize));
end
NO=oldNO;
close(gcbf)

%----------------------
function close_Callback 
%----------------------
global DATA

try
  DATA.GUI.MergeStacks = close(DATA.GUI.MergeStacks);
catch me
  mydispexception(me)%#ok<CTCH>
  DATA.GUI.MergeStacks =[];
  delete(gcbf);
end
