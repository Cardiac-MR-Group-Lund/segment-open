function [varargout] = export(varargin)
%Export functionality in Segment

%Einar Heiberg

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%----------------------------------
function exportslicevolume_Callback %#ok<DEFNU>
%----------------------------------
%Export slicebased volume
global SET NO DATA

if SET(NO).Rotated
  myfailed('Rotated image stacks currently not supported for this operation.',DATA.GUI.Segment);
  return;
end;

endo = true;
epi = true;
rvendo = true;
rvepi = true;

if isempty(SET(NO).EndoX)
  endo = false;
end;
if isempty(SET(NO).EpiX)
  epi = false;
end;
if isempty(SET(NO).RVEndoX)
  rvendo = false;
end;
if isempty(SET(NO).RVEpiX)
  rvepi = false;
end;

numcontours = endo+epi+rvendo+rvepi;

if (SET(NO).TSize>1)&&(SET(NO).Longaxis>1)
  mywarning('Automated long axis compensation not included in the export.',DATA.GUI.Segment);
end;
  
%Reserve memory
outdata = cell(SET(NO).ZSize+8,1+SET(NO).TSize*numcontours);

%Write header
outdata{1,1} = 'FileName';
outdata{1,2} = SET(NO).FileName;
outdata{1,4} = 'Time';
outdata{1,5} = datestr(now);

%Export total volume
outdata{2,1} = 'TimeFrame';
outdata{3,1} = 'TotalEndo';
outdata{4,1} = 'TotalEpi';
outdata{5,1} = 'TotalRVEndo';
outdata{6,1} = 'TotalRVEpi';

for tloop=1:SET(NO).TSize
  outdata{2,1+tloop} = sprintf('%d',tloop);
  if endo
    outdata{3,1+tloop} = SET(NO).LVV(tloop);
  end;
  if epi
    outdata{4,1+tloop} = SET(NO).EPV(tloop);
  end;
  if rvendo
    outdata{5,1+tloop} = SET(NO).RVV(tloop);    
  end;
  if rvepi
    outdata{6,1+tloop} = SET(NO).RVEPV(tloop);    
  end;  
end;

%Find number of ROI's
temproinames=cell(1,SET(NO).RoiN);
for rloop=1:SET(NO).RoiN
  temproinames{rloop}=SET(NO).Roi(rloop).Name;
end
roinames = union(temproinames,{});

rowoffset = 8+length(roinames);

for zloop=1:SET(NO).ZSize
  outdata{rowoffset+zloop,1} = sprintf('Slice:%d',zloop);
end;

coloffset = 1;
if endo
  outdata = exportslicehelper(outdata,rowoffset,coloffset,'Endo',...
    SET(NO).EndoX,SET(NO).EndoY);
  coloffset = coloffset+SET(NO).TSize;
end;
if epi
  outdata = exportslicehelper(outdata,rowoffset,coloffset,'Epi',...
    SET(NO).EpiX,SET(NO).EpiY);
  coloffset = coloffset+SET(NO).TSize;  
end;
if rvendo
  outdata = exportslicehelper(outdata,rowoffset,coloffset,'RVEndo',...
    SET(NO).RVEndoX,SET(NO).RVEndoY);
  coloffset = coloffset+SET(NO).TSize;
end;
if rvepi
  outdata = exportslicehelper(outdata,rowoffset,coloffset,'RVEpi',...
    SET(NO).RVEpiX,SET(NO).RVEpiY);
  coloffset = coloffset+SET(NO).TSize;  
end;

%Find ROI's
for loop=1:length(roinames)
  
  %Write ROI name
  outdata{6+loop,1} = roinames{loop};
  
  %Calculate volume for that ROI
  roivolume = nan(SET(NO).ZSize,SET(NO).TSize);
  for rloop=1:SET(NO).RoiN
    if isequal(SET(NO).Roi(rloop).Name,roinames{loop})
      for zloop=1:SET(NO).ZSize
        if ~isempty(find(SET(NO).Roi(rloop).Z==zloop,1))
          for tloop=1:SET(NO).TSize
            if ~isempty(find(SET(NO).Roi(rloop).T==tloop,1))
              temp = (SET(NO).SliceThickness+SET(NO).SliceGap)*stablepolyarea(...
                SET(NO).ResolutionY*SET(NO).Roi(rloop).Y(:,tloop),...
                SET(NO).ResolutionX*SET(NO).Roi(rloop).X(:,tloop))/1000;
              if ~isnan(temp)
                if isnan(roivolume(zloop,tloop))
                  roivolume(zloop,tloop) = 0;
                end;
                roivolume(zloop,tloop) = roivolume(zloop,tloop)+temp;
              end;
            end
          end;
        end;
      end;
    end;
  end;
  
  %Export ROI volume
  temp = roivolume;
  temp(isnan(temp)) = 0;
  temp = sum(temp,1);
  
  for tloop=1:SET(NO).TSize
    outdata{rowoffset,tloop+coloffset} = sprintf('%s:TF%d',roinames{loop},tloop);
    outdata{6+loop,1+tloop} = temp(tloop);    
  end;
  for zloop=1:SET(NO).ZSize
    for tloop=1:SET(NO).TSize
      outdata{rowoffset+zloop,tloop+coloffset} = roivolume(zloop,tloop);
    end;
  end;
  
  %Increae column offset
  coloffset = coloffset+SET(NO).TSize;    
end;

segment('cell2clipboard',outdata);

%-------------------------------------------------------------------------
function outdata = exportslicehelper(outdata,rowoffset,coloffset,type,x,y)
%-------------------------------------------------------------------------
%Helper function to export slice based data.
global SET NO

for tloop=1:SET(NO).TSize
  outdata{rowoffset,tloop+coloffset} = sprintf('%s:TF%d',type,tloop);
end;

for zloop=1:SET(NO).ZSize
  for tloop=1:SET(NO).TSize
    outdata{rowoffset+zloop,tloop+coloffset} = (SET(NO).SliceThickness+SET(NO).SliceGap)*stablepolyarea(...
      SET(NO).ResolutionY*y(1:end-1,tloop,zloop),...
      SET(NO).ResolutionX*x(1:end-1,tloop,zloop))/1000;
  end;
end;

%------------------------------
function exportcontour_Callback %#ok<DEFNU>
%------------------------------
%Export contour. Ask what contour to take.
global SET NO DATA

[x,y,name] = segment('askcontour','Choose which contour to export.');
      
if isempty(x)
  myfailed('User aborted or no contour available.',DATA.GUI.Segment);
  return;
end;

%Find slices with contour
totinslice = size(x,1)*size(x,2);
slices = find(squeeze(sum(sum(isnan(x),1),2))~=totinslice);

%l1 <filename>: <contour>  <resolutionx> <resolutiony>
%l2 filename contourname xres yres
%l3 empty line
%l4 <slice> <xtf1> <ytf1> <xtf2> <ytf3>
%l5 
outdata = cell(5+length(slices)*size(x,1),1+2*size(x,2));

%write header
outdata{1,1} = 'Filename:';
outdata{2,1} = SET(NO).FileName;
outdata{1,2} = 'Contour:';
outdata{2,2} = name;
outdata{1,3} = 'ResolutionX';
outdata{2,3} = SET(NO).ResolutionX;
outdata{1,4} = 'ResolutionY';
outdata{2,4} = SET(NO).ResolutionY;

outdata{4,1} = 'Slice';
for tloop=1:size(x,2)
  outdata{4,2+2*(tloop-1)} = sprintf('X_tf%02d',tloop);
  outdata{4,3+2*(tloop-1)} = sprintf('Y_tf%02d',tloop);  
end;
for zloop=1:length(slices)
  %write slice number
  for nloop=1:size(x,1)
    outdata{4+(zloop-1)*size(x,1)+nloop,1} = slices(zloop);
  end;
  for tloop=1:size(x,2)
    for nloop=1:size(x,1)
      outdata{4+(zloop-1)*size(x,1)+nloop,2+2*(tloop-1)} = SET(NO).ResolutionX*x(nloop,tloop,slices(zloop));
      outdata{4+(zloop-1)*size(x,1)+nloop,3+2*(tloop-1)} = SET(NO).ResolutionY*y(nloop,tloop,slices(zloop));      
    end;
  end;
end;

segment('cell2clipboard',outdata);

%-----------------------------------------
function exportimage_Callback(image2store) %#ok<DEFNU>
%-----------------------------------------
%Export an image to an image file. If no input image exist, export current image with current slice and timeframe 
global DATA SET NO

no = NO;

multipleslices = false;
templatename = 'image';
if nargin==0
  if (SET(NO).EndSlice-SET(NO).StartSlice)>0
    if yesno('Several images selected. Export all?')
      multipleslices = true;
      templatename = 'image_slice000';
    end;
  end;
end;

try
  temppwd = pwd;
  if exist(DATA.Pref.exportpath,'dir')
    cd(DATA.Pref.exportpath);
  else
    mydisp('Warning: Export path does not exist, please check preferences.');
  end;
  
  [filename, pathname,filterindex] = myuiputfile(...
    { '*.png','PNG image (*.png)';...
    '*.jpg','JPEG image (*.jpg)';...
    '*.bmp','BMP image (*.bmp)';...
    '*.tif','TIFF image (*.tif)'},...
    'Save file as',templatename);
  cd(temppwd);
  
  if isequal(filename,0)
    myfailed('Operation cancelled.',DATA.GUI.Segment);
    return;
  end;
  
  if multipleslices
    filename = removenumerics(filename);    
    slices = SET(no).StartSlice:SET(no).EndSlice;
  else
    slices = SET(no).CurrentSlice;
  end;  
  
  if multipleslices
    h = waitbar(0,'Please wait.');
  end;
  for loop = 1:length(slices);
  
    %Get the image
    if nargin < 1
      if multipleslices
        image2store = SET(no).IM(:,:,SET(no).CurrentTimeFrame,slices(loop));
        
        if DATA.Pref.ViewInterpolated
          scale = 2;
          imxsz = SET(no).XSize*scale;
          imysz = SET(no).YSize*scale;
          image2store = imresize(image2store,[imxsz imysz],'bilinear');      
        end;
        
        image2store = calcfunctions('remapuint8viewim',image2store,no);
      else
        image2store = squeeze(DATA.ViewIM{DATA.CurrentPanel}(:,:,SET(no).CurrentTimeFrame,:));
      end;
    end
  
    %Add extension if necessary
    [~,filename,~] = fileparts(filename);
    
    if multipleslices
      f = fullfile(pathname,sprintf('%s%03d',filename,slices(loop)));
    else
      f = fullfile(pathname,filename);
    end;
  
    switch filterindex
      case 1
        f = [f '.png']; %#ok<AGROW>
      case 2
        f = [f '.jpg']; %#ok<AGROW>
      case 3
        f = [f '.bmp']; %#ok<AGROW>
      case 4
        f = [f '.tif']; %#ok<AGROW>
    end;    
    
    switch filterindex
      case 1
        imwrite(image2store,f,'png','bitdepth',8,'software','Segment',...
          'creationtime',datestr(now));
      case 2
        imwrite(image2store,f,'jpg','quality',100);
      case 3
        imwrite(image2store,f,'bmp');
      case 4
        imwrite(image2store,f,'tif');
    end;
    
    if multipleslices
      waitbar(loop/length(slices),h);
    end;
  
  end; %Loop over slices (either only one or many)
  
  if multipleslices
    close(h);
  end;
    
catch me
  mydispexception(me);
  myfailed('Export of image failed.');
  cd(temppwd);
  return;
end
  
cd(temppwd); %Just to be sure
mymsgbox('Export successful.','Done!');

%-----------------------------------
function stri = removenumerics(stri)
%-----------------------------------
%Remove numbers from a filename
stri = stri( (stri<47) | (stri>57) );

%---------------------------
function screenshot_Callback %#ok<DEFNU>
%---------------------------
%Create an image file of a screenshot of the main axis
global DATA SET

temp = pwd;
if exist(DATA.Pref.exportpath,'dir')
  cd(DATA.Pref.exportpath);
else
  mydisp('Warning: Export path does not exist, please check preferences.');
end;

m = mymenu('Save to file or to PACS?','Save to file','Save to PACS',...
  'Save for inclusion in report');

switch m
  case 1
    [filename, pathname,filterindex] = myuiputfile(...
      { '*.png','PNG image (*.png)';...
      '*.jpg','JPEG image (*.jpg)';...
      '*.bmp','BMP image (*.bmp)';...
      '*.tif','TIFF image (*.tif)'},...
      'Save file as','screenshot');
  case 2 
    filename = 'screenshot.dcm';
    pathname = getpreferencespath;
    filterindex = 5;
  case 3
    name = removeforbiddenchars(SET(1).PatientInfo.Name);
    if isempty(name)
      name = 'Hidden';
    end;
    pathname = fullfile(DATA.Pref.Pacs.ReportsheetPath, ...
      name,'Screenshots');
    if ~exist(pathname,'dir')
      sucess = mkdir(pathname);
      if ~sucess
        myfailed('Could not create screenshot directory. Aborted');
        return
      end
    end
    reportdir = dir(fullfile(pathname,'screenshot*.png'));
    nbrs = zeros(1,numel(reportdir));
    for i = 1:numel(reportdir)
      nbrs(i) = sscanf(reportdir(i).name(12:15),'%f');
    end
    newnbr = min(setdiff(1:9999,nbrs));
    filename = sprintf('screenshot_%04.0f.png',newnbr);
    filterindex = 1;
  case 0
    filename = 0;
end

cd(temp);
set(DATA.Handles.screenshoticon,'state','off');
if isequal(filename,0)
  myfailed('Operation cancelled.',DATA.GUI.Segment);
  return;
end;

h = DATA.Handles.imageaxes(DATA.CurrentPanel);
frame = mygetframe(h);
im = frame2im(frame);
f = fullfile(pathname,filename);

%Add extension if necessary
[~,~,ext] = fileparts(f);
if isempty(ext)
  switch filterindex
    case 1
      f = [f '.png'];
    case 2
      f = [f '.jpg'];
    case 3
      f = [f '.bmp'];
    case 4
      f = [f '.tif'];
    case 5
      f = [f '.dcm'];
  end;
end;

switch filterindex
  case 1
    imwrite(im,f,'png','bitdepth',8,'software','Segment',...
      'creationtime',datestr(now));
  case 2
    imwrite(im,f,'jpg','quality',100);
  case 3
    imwrite(im,f,'bmp');
  case 4
    imwrite(im,f,'tif');
  case 5
    makeimagedicom(im,f,DATA.ViewPanels(DATA.CurrentPanel));
    try
      pacs('savetopacs_helper',{f},getpreferencespath);
    catch me
      mydispexception(me);
      myfailed(me.message);
    end
end;
%catch
%  myfailed('Export of image failed.');
%  return;
%end

mymsgbox('Export successful.','Done!');

%--------------------------------------------------------
function ok = exportsavemovie(mov,left,right,up,down,fps)
%--------------------------------------------------------
%Exports a move as an avi file or a set of png-files.
%- mov is a movie struct.
%- left,right,up,down are crop coordinates.
%- fps is frame rate.
%
%Allows user to select different codecs.

global DATA 

if exist(DATA.Pref.exportpath,'dir')
  epath = DATA.Pref.exportpath;
else
  mydisp('Warning: Export path does not exist, please check preferences.');
  epath = pwd;
end;

if nargin<2
  left = 1;
end;

if nargin<3
  right = size(mov(1).cdata,2);
end;

if nargin<4
  up = 1;
end;

if nargin<5
  down = size(mov(1).cdata,2);
end;

if nargin<6
  fps = 15;
end;

ok = false; %Pessimistic

m = mymenu('Save as',{'avi','png-files','animated gif'},DATA.GUI.Segment);
switch m
  case 1
    %avi
    m = mymenu('Select video compressor',...
      {'None' ...
      'Motion JPEG' ...
      'Motion JPEG 2000' ...
      'Motion JPEG 2000 (lossless)' ...
      },'Motion JPEG',DATA.GUI.Segment);
    
    if isequal(m,0)
      myfailed('User pressed cancel',DATA.GUI.Segment);
      return;
    end;

    [filename, pathname] = myuiputfile(...
      '*.avi',...
      'Save avi-file as',...
      fullfile(epath,'movie.avi'));

    if isequal(filename,0)
      myfailed('Operation cancelled.',DATA.GUI.Segment);
      return;
    end;
    
    %Set the compressor
    switch m
      case 1
        profile = 'Uncompressed AVI';
      case 2
        profile = 'Motion JPEG AVI';
      case 3
        profile = 'Motion JPEG 2000';
      case 4
        profile = 'Archival';
      otherwise
        profile = 'None';
    end;
    
    %Open the file
    fid = VideoWriter(fullfile(pathname,filename), profile);
    fid.FrameRate = fps;
    fid.open;
    
    %Add frames
    myworkon;
    h = mywaitbarstart(length(mov),'Please wait, storing movie.');
    myadjust(h.h,DATA.GUI.Segment);
    for loop=1:length(mov)
      temp = mov(loop).cdata;
      if nargin>1
        temp = temp(up:down,left:right,:); %EH: 2009-08-31
      end;
      fid.writeVideo(temp);
      h = mywaitbarupdate(h);
    end;
    mywaitbarclose(h);
    myworkoff;
    
    %Close the file
    fid.close; 

  case 2
    %png-files
    [filename, pathname] = myuiputfile(...
      '*.*',...
       'Save files as',...
      [epath filesep 'im']);
    handles.filename = fullfile(pathname,filename);    
    [pathname,filename] = fileparts(handles.filename);
    handles.filename = fullfile(pathname,filename);
    
    myworkon;
    h = mywaitbarstart(length(mov),'Please wait, storing movie.');
    myadjust(h,DATA.GUI.Segment)
    for loop=1:length(mov)
      [temp] = frame2im(mov(loop));
      if nargin>1
        temp = temp(up:down,left:right,:);
      end;
      imwrite(temp,[handles.filename sprintf('%05d.png',loop)],'png');
      h = mywaitbarupdate(h);
    end;
    mywaitbarclose(h);
    myworkoff;
  case 3
    %animated gif
    
    %Select filename
    [filename, pathname] = myuiputfile(...
      '*.gif',...
      'Save animated gif-file as',...
      fullfile(epath,'movie.gif'));

    if isequal(filename,0)
      myfailed('Operation cancelled.',DATA.GUI.Segment);
      return;
    end;
    
    dither = yesno('Do you want to use dithering for better colormap?');
    if dither
      ditherstring = 'dither';
    else
      ditherstring = 'nodither';
    end;
    
    myworkon;
    h = mywaitbarstart(length(mov),'Please wait, storing movie.');
    
    for loop = 1:length(mov)
      [imrgb] = frame2im(mov(loop));
      imrgb = imrgb(up:down,left:right,:);

      if isequal(loop,1)
        [im,map] = rgb2ind(imrgb,256,ditherstring);
        im(1,1,1,length(mov)) = 0; %Reserve memory
      else
        im(:,:,1,loop) = rgb2ind(imrgb,map,ditherstring);
      end;
      
      h = mywaitbarupdate(h);
    end;

    imwrite(im,map,[pathname filesep filename],'DelayTime',0,'LoopCount',inf) %g443800

    mywaitbarclose(h);
    myworkoff;
    
  otherwise
    myfailed('Saving aborted.',DATA.GUI.Segment);
    return;
end;

ok = true;

%-----------------------------------------
function exportmovierecorder_Callback(arg)
%-----------------------------------------
%Movie recorder GUI

global DATA SET NO

persistent handles

if nargin==0
  
  if SET(NO).TSize==1
    myfailed('Need more than one timeframe to record movie.',DATA.GUI.Segment);
    return;
  end;
  
  set(DATA.Handles.movierecordericon,'state','off');
  
  fig = openfig('movierecorder.fig','reuse');
  myadjust(fig,DATA.GUI.Segment);
  blockfig(fig);
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

  % Generate a structure of handles to pass to callbacks, and store it.
  handles = guihandles(fig);
  handles.fig = fig;
  handles.recording = false;
  handles.frames = SET(NO).TSize;
  handles.fps = 15;
  handles.storedframes = 0;
  handles.crop = false;
  handles.size = [];
  DATA.Record = true;
  axis(handles.previewaxes,'off');
  set(handles.framesedit,'string',sprintf('%d',handles.frames));
  mymsgbox([...
    'Force an image update by pressing ''next'' button in the GUI that ' ...
    'you want to record. This will then to display a movie recording ' ...
    'preview. If required set cropping limits, and adjust the number of '...
    'timeframes recorded'],'');
else
  switch arg
    case 'fps'
      s = mygetedit(handles.fpsedit);
      temp = str2double(s);
      if isnan(temp)
        set(handles.fpsedit,'String',sprintf('%d',handles.fps));
        myfailed('Could not interpret number.',handles.fig);
        return;
      end;
      handles.fps = temp;
      set(handles.fpsedit,'String',sprintf('%d',handles.fps));      
    case 'frames'
      s = mygetedit(handles.framesedit);
      temp = str2double(s);
      if isnan(temp)
        set(handles.framesedit,'String',sprintf('%d',handles.frames));
        myfailed('Could not interpret number.',handles.fig);
        return;
      end;
      handles.frames = temp;
      set(handles.framesedit,'String',sprintf('%d',handles.frames));
    case 'newframe'
      if not(handles.recording)
        %Not recording, show preview.
        temp = frame2im(DATA.MovieFrame);
        handles.size = size(temp);
        image(temp,'parent',handles.previewaxes);
        axis(handles.previewaxes,'image','off');
      else
        %Recording
        if handles.storedframes<handles.frames
          if handles.storedframes==0
            handles.movie = DATA.MovieFrame;
            handles.storedframes = handles.storedframes+1;
          else
            handles.storedframes = handles.storedframes+1;
            handles.movie(handles.storedframes) = DATA.MovieFrame;
          end;
          mydisp(dprintf('Frame %d stored.',handles.storedframes));
        end;
        if handles.storedframes==handles.frames
          exportmovierecorder_Callback('store');
          return;
        end;
      end;
    case 'crop'
      if isempty(handles.size)
        myfailed('No preview available.',handles.fig);
        return;
      end;
      mymsgbox('Click in the upper-left corner and the bottom-right corner to set crop limits.','');
      [x,y] = ginput(2);
      hold(handles.previewaxes,'on');
      plot(handles.previewaxes,...
        [x(1) x(2) x(2) x(1) x(1)],...
        [y(1) y(1) y(2) y(2) y(1)],'r-');
      hold(handles.previewaxes,'off');
      x = round(x);
      y = round(y);
      handles.left = max(min(x),1);
      handles.right = min(max(x),handles.size(2));
      handles.up = max(min(y),1);
      handles.down = min(max(y),handles.size(1));
      handles.crop = true;
    case 'nocrop'
      handles.left = 1;
      handles.right = handles.size(2);
      handles.up = 1;
      handles.down = handles.size(1);
      handles.crop = false;
    case 'start'
      DATA.Record = true;
      handles.recording = true;
      mymsgbox('Starting to record. Each image update will cause a new frame in the image.','');
    case 'stop'
      handles.recording = false;
      DATA.Record = false;
      if handles.storedframes>1
        exportmovierecorder_Callback('store');
      end;
    case 'store'
      mydisp('Storing movie to disc.');
      if not(handles.crop)
        handles.left = 1;
        handles.right = handles.size(2);
        handles.up = 1;
        handles.down = handles.size(1);
      else
        handles.left = max(handles.left,1);
        handles.right = min(handles.right,handles.size(2));
        handles.up = max(handles.up,1);
        handles.down = min(handles.down,handles.size(1));
        handles.crop = true;
      end;
      ok = exportsavemovie(handles.movie,handles.left,handles.right,handles.up,handles.down,handles.fps);
      if ~ok
        myfailed('Export of movie failed.',handles.fig);
        return;
      end;
      handles.movie = [];
      handles.storedframes = 0;
      DATA.Record = false;
    case 'close'
      DATA.Record = false;
      close(handles.fig);
      handles = [];
    otherwise
      myfailed(dprintf('Unknown option %s for exportmovierecorder',arg),handles.fig);
      return;
  end;
end;

%---------------------------
function exportmovie_Callback %#ok<DEFNU>
%---------------------------
%Function to export a movie without contours. This is a quick method to
%generate a movie of the current image stack.

global DATA SET NO

no = DATA.ViewPanels(DATA.CurrentPanel);
myworkon;
h = mywaitbarstart(SET(no).TSize,'Please wait, generating movie.');
for tloop=1:SET(no).TSize
  %Extract image and convert to index image
  temp = squeeze(DATA.ViewIM{DATA.CurrentPanel}(:,:,tloop,:,:));
  if isempty(SET(NO).Colormap)
    temp = repmat(temp,[1 1 3]);
    mov(tloop) = im2frame(temp); %#ok<AGROW>
  else
    mov(tloop) = im2frame(temp,SET(NO).Colormap); %#ok<AGROW>    
  end;
  h = mywaitbarupdate(h);
end;
mywaitbarclose(h);

if ~DATA.Testing
  s = inputdlg({'Number of frames per second'},'fps',1,{sprintf('%d',15)});
else
  s{1} = popfrombuffer('String');
end
if isempty(s)
  myfailed('Invalid fps.',DATA.GUI.Segment);
  return;
else
  [fps,ok] = str2num(s{1}); %#ok<ST2NM>
  if not(ok)
    myfailed('Invalid fps.',DATA.GUI.Segment);
    return;
  end;
end;

%Call to save to disk
ok = exportsavemovie(mov,1,size(temp,2),1,size(temp,1),fps);
myworkoff;

if ok
  mymsgbox('Export successful.','Done!');
end;

%----------------------------------
function exportvolumecurve_Callback %#ok<DEFNU>
%----------------------------------
%Export volume curve to clipboard.
global SET NO

stri = [];
stri = [stri ...
  sprintf('Name:\t%s\n',SET(NO).PatientInfo.Name) ...
  sprintf('\n\n') ...
  sprintf('Time[ms]\tEndovol[ml]\tEpivol[ml]\tPapvol\tRVEndovol[ml]\tRVEpivol[ml]\n')];
for loop=1:SET(NO).TSize
  stri = [stri ...
    sprintf('%f\t%f\t%f\t%f\t%f\t%f\n',...
    SET(NO).TimeVector(loop)*1000,...
    SET(NO).LVV(loop),...
    SET(NO).EPV(loop),...
    SET(NO).PV(loop),...
    SET(NO).RVV(loop),...
    SET(NO).RVEPV(loop))...
    ]; %#ok<AGROW>
end;

clipboard('copy',stri);
mymsgbox('Results copied to clipboard','Done!');

%----------------------------------
function mmodemeasurements_Callback %#ok<DEFNU>
%----------------------------------
%Export mmode measurements
global SET

stri = [];
stri = [stri ...
  sprintf('Name:\t%s\n',SET(1).PatientInfo.Name) ...
  sprintf('\n\n') ...
  sprintf('Image\tDistance[mm]\tTime[ms]\n')];
for no=1:numel(SET)
  [spacedist,timedist] = calcfunctions('calcmmodedists',no);
  stri = [stri ...
    sprintf('%s\t%f\t%f\n',...
    SET(no).ImageViewPlane, spacedist, timedist) ...
    ]; %#ok<AGROW>
end;

clipboard('copy',stri);
mymsgbox('Results copied to clipboard','Done!');

%---------------------------------
function exportendoepitoensight_Callback %#ok<DEFNU>
%---------------------------------
%Exports endocardium and epicardium in Ensight format.

global DATA SET NO

if isempty(SET(NO).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end;

% Ask user for directory and filename
[outfile, outdir] = myuiputfile;

timesteps = SET(NO).TimeVector;

if not(isempty(SET(NO).EndoX))
	x = SET(NO).EndoX;
	y = SET(NO).EndoY;
	z = (0:size(x,3)-1);
	matlab2ensightgolddeformingmesh(outdir, outfile, 'Endo', x, y, z, timesteps);
end;

if not(isempty(SET(NO).EpiX))
	x = SET(NO).EpiX;
	y = SET(NO).EpiY;
	z = (0:size(x,3)-1);
	matlab2ensightgolddeformingmesh(outdir, outfile, 'Epi', x, y, z, timesteps);
end;

if not(isempty(SET(NO).RVEndoX))
	x = SET(NO).RVEndoX;
	y = SET(NO).RVEndoY;
	z = (0:size(x,3)-1);
	matlab2ensightgolddeformingmesh(outdir, outfile, 'RVEndo', x, y, z, timesteps);
end;

if not(isempty(SET(NO).RVEpiX))
	x = SET(NO).RVEpiX;
	y = SET(NO).RVEpiY;
	z = (0:size(x,3)-1);
	matlab2ensightgolddeformingmesh(outdir, outfile, 'RVEpi', x, y, z, timesteps);
end;

%---------------------------------
function exportlv2ensight_Callback %#ok<DEFNU>
%---------------------------------
%Export left ventricle to Ensight.

global DATA SET NO

if isempty(SET(NO).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end;

% Started by Einar Heiberg (?)
% Seriously by Johannes Toger Jan 20, 2009

% Extract coordinates from Segment
% Result: size(x) == size(y) == [<points along contour> <z> <time>]
x = permute(SET(NO).EndoX,[1 3 2]);
y = permute(SET(NO).EndoY,[1 3 2]);

Npoints = size(x,1);
Nz = size(x,2);
Ntimesteps = size(x,3);

% Get spacing info
dz = (SET(NO).SliceThickness + SET(NO).SliceGap);
dz = dz*1e-3; % mm -> m
dx = SET(NO).ResolutionX*1e-3;
dy = SET(NO).ResolutionY*1e-3;

% Rescale x and y
x = x*dx;
y = y*dy;

% Shift -(1/2,1/2) in (x,y).
x = x - dx/2;
y = y - dy/2;

% Construct Z array
zcoords = (0:Nz-1)*dz;
zcoords = fliplr(zcoords); 
z = repmat(zcoords,[Npoints,1,Ntimesteps]);

% Segment directions in LPS
% LPS = Left, Posterior, Superior: DICOM patient coordinate system, which
% we use in Ensight and all 3D visualization.
xdir = SET(NO).ImageOrientation(4:6);
ydir = SET(NO).ImageOrientation(1:3);
zdir = cross(ydir, xdir);

% Get origin
%origin = SET(NO).ImagePosition*1e-3;

% Z direction is flipped, fix it.
origin = SET(NO).ImagePosition/1e3 - zdir*((Nz-1)*dz); 

% Convert to LPS
L = x*xdir(1) + y*ydir(1) + z*zdir(1) + origin(1);
P = x*xdir(2) + y*ydir(2) + z*zdir(2) + origin(2);
S = x*xdir(3) + y*ydir(3) + z*zdir(3) + origin(3);

timesteps = SET(NO).TimeVector;

% Since all slices aren't used typically, x and y have NaN:s for those
% slices (just what Segment does, it seems). These are still in L,P and S.
% Find them and remove them.
xrow = x(1,:,1);
usedslices = ~isnan(xrow); 
usedslices = find(usedslices); % list of used slice indices

L = L(:,usedslices,:);
P = P(:,usedslices,:);
S = S(:,usedslices,:);

% Two periods (Quick hack)
writetwoperiods = yesno('Extend to two periods?'); % ask user

if writetwoperiods
    % Extend everything
    L = repmat(L,[1 1 2]);
    P = repmat(P,[1 1 2]);
    S = repmat(S,[1 1 2]);
    timesteps = SET(NO).TIncr*(0:(2*SET(NO).TSize-1));
end

% Ask user for directory and filename
[outfile, outdir] = myuiputfile;

% Output to Ensight 5 format
% MATLAB2ENSIGHT(X, Y, Z, U, FILEPATH, FILENAME, PART_DESC, VAR_DESC, 
%                TIME_VECTOR, VARIABLE_GEOMETRY);
matlab2ensight(...
  L, ... %x
  P, ... %y
  S, ... %z
  zeros(size(L)), ... %data
  outdir, ... %pathname
  outfile, ... %filename
  'LV Endocardium', ... %part desc
  'nothing', ... %var desc
  timesteps, ...
  1);

%-----------------------------------------------------
function exporttoclipboardthisstack_Callback(doheader) %#ok<DEFNU>
%-----------------------------------------------------
%Exports data for current image stack to clipboard.
global NO

if nargin<1
  doheader = false;
end;

exporttoclipboard_Callback(doheader,NO);

%-----------------------------------------------
function exporttoclipboard_Callback(doheader,no) 
%-----------------------------------------------
%Export data to clipboard. Calls another function to do the export.
if nargin==0
  doheader = true;
end;

includenormalized = yesno('Do you want to include BSA normalized values?');

if nargin<2
  exportdata(doheader,includenormalized);
else
  exportdata(doheader,includenormalized,no);
end;  

%---------------------------------------
function [outdata,ind] = header(onlyone)
%---------------------------------------
%Helper function to write header when exporting data.
if ~onlyone
  outdata{1, 1} = 'FileName';
else
  outdata{1, 1} = 'ImageStack';
end;
outdata{1, 2} = 'PatientName';
outdata{1, 3} = 'PatientID';
outdata{1, 4} = 'AcquisitionDate';
outdata{1, 5} = 'Age';
outdata{1, 6} = 'Height[cm]';
outdata{1, 7} = 'Weight[kg]';
outdata{1, 8} = 'Sex';
outdata{1, 9} = 'BSA[m2]';
outdata{1,10} = 'HeartRate[bpm]';
outdata{1,11} = 'R-R[ms]';
outdata{1,12} = 'LVM[ml]';
outdata{1,13} = 'LVM[g]';
outdata{1,14} = 'LVMI[g/m2]';
outdata{1,15} = 'EDV[ml]';
outdata{1,16} = 'EDVI[ml/m2]';
outdata{1,17} = 'ESV[ml]';
outdata{1,18} = 'ESVI[ml/m2]';
outdata{1,19} = 'SV[ml]';
outdata{1,20} = 'SVI[ml/m2]';
outdata{1,21} = 'EF[%]';
outdata{1,22} = 'CO[l/min]';
outdata{1,23} = 'CI[l/min]';
outdata{1,24} = 'PFR[ml/s]';
outdata{1,25} = 'PER[ml/s]';
outdata{1,26} = 'RVM[ml]';
outdata{1,27} = 'RVM[g]';
outdata{1,28} = 'RVMI[g/m2]';
outdata{1,29} = 'RVEDV[ml]';
outdata{1,30} = 'RVEDVI[ml/m2]';
outdata{1,31} = 'RVESV[ml]';
outdata{1,32} = 'RVESVI[ml/m2]';
outdata{1,33} = 'RVSV[ml]';
outdata{1,34} = 'RVSVI[ml/m2]';
outdata{1,35} = 'RVEF[%]';
outdata{1,36} = 'LVM_DE[ml]';
outdata{1,37} = 'Scar[%]';
outdata{1,38} = 'Scar[ml](fromLVMDE)';
outdata{1,39} = 'TotExt[%]';
outdata{1,40} = 'MeanTransmurality[%]';
outdata{1,41} = 'MaxTransmurality[%]';
outdata{1,42} = 'TotExt[%](Weighted)';
outdata{1,43} = 'MeanTransmurality[%](Weighted)';
outdata{1,44} = 'MaxTransmurality[%](Weighted)';
outdata{1,45} = 'MO region[%]';
outdata{1,46} = 'MO[%]';
outdata{1,47} = 'MaR[%] ED';
outdata{1,48} = 'MaR[%] ES';
outdata{1,49} = 'MaR[ml] ED';
outdata{1,50} = 'MaR[ml] ES';
outdata{1,51} = 'MaR TPD';
outdata{1,52} = 'MaR TPD in LAD';
outdata{1,53} = 'MaR TPD in LCx';
outdata{1,54} = 'MaR TPD in RCA';
outdata{1,55} = 'Roi Name';
outdata{1,56} = 'Roi Total Flow-Volume [ml]';
outdata{1,57} = 'Roi Total Flow [ml/min]'; %SBT20160620
outdata{1,58} = 'Flow Heart-Rate [bpm]'; %SBT20160620

outdata{1,59} = 'Measurement Name';
outdata{1,60} = 'Measurement Distance [mm]';
outdata{1,61} = 'RVPFR[ml/s]';
outdata{1,62} = 'RVPER[ml/s]';
outdata{1,63} = 'LVM measured at ED';
outdata{1,64} = 'LVM measured at ES';


ind = [9 14 16 18 20 23 28 30 32 34];

%-------------------------------------------
function exportmultiple_Callback(dosegdicom) %#ok<DEFNU>
%-------------------------------------------
%Creaty summary of multiple matfiles in one folder
%This function is very useful for research. The user 
%performs all delineations and then exports all data to
%one spread sheet.

global DATA SET NO

if nargin < 1
  dosegdicom = false;
end
if dosegdicom
  suffix = 'segdicom';
else
  suffix = 'mat';
end

%Select path
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,sprintf('Select a folder with .%s files',suffix));
if isequal(pathname,0)
  myfailed('Aborted.',DATA.GUI.Segment);
  return;
end;

%Find files to process
files2load = dir([pathname filesep sprintf('*.%s',suffix)]);
numfiles = length(files2load);

if numfiles==0
  myfailed('Found no files to summarize.',DATA.GUI.Segment);
  return;
end;

includenormalized = yesno('Do you want to include BSA normalized values?');

%Create output matrix
outdata = cell(numfiles+1,60); %+1 since header, 58 since header size %SBT20160620: Changed colsize from 58 to 60
[temp,indforbsa] = header(false); %False means not only one image stack

%--- Write header
for loop=1:length(temp)
  outdata{1,loop} = temp{loop};
end;

%Loop over all files
numrows=1;
h = mywaitbarstart(numfiles,'Please wait, loading and summarizing files.',1);
for fileloop=1:numfiles

  %--- Load file
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
  
  disp(dprintf('Loading %s.',files2load(fileloop).name)); %#ok<DSPS>
  
  %Set filename
  outdata{numrows+1,1} = files2load(fileloop).name;
  
  SET = []; % %Make sure a fresh start
  
  %---- try
  % This try-catch clause used to be commented out, which caused errors
  % in one file to crash the whole export process. Reintroducing it. /NL
  try
    %Load
    SET=[]; 
    if ~dosegdicom
      load([pathname filesep files2load(fileloop).name],'-mat');
    else
      loader = segloader();
      loader.adddicomfiles({fullfile(pathname,files2load(fileloop).name)});
      [~, setstruct] = loader.render(pathname, []);      
      clear loader;
    end
    
    %Assign
    SET = setstruct;
    clear setstruct;
    
    NO = 1;
    
    %Call to intialize all variables correcly after loaded data.
    openfile('setupstacksfrommat',1);
    segment('renderstacksfrommat');
    
    %Call one to get info for current file
    temp = exportdata(false,true); %false means no header; true normalized will be removed later
  
    %Copy the data
		for col=2:size(temp,2);
			for row=1:size(temp,1);
				outdata{numrows+row,col} = temp{row,col};
			end;
		end
		
		numrows=numrows+size(temp,1);
    
  catch me
    %--- Some thing went wrong
    mydispexception(me);
    outdata{fileloop+1,2} = 'FAILED.';
  end
  h = mywaitbarupdate(h);
end; %loop over files
mywaitbarclose(h);

%If not wanted remove non normalized values.
ind = true(1,size(outdata,2));
if ~includenormalized
  ind(indforbsa) = false;
end;
outdata = outdata(:,ind);

%--- Output to a string
segment('cell2clipboard',outdata);

%Make sure starting with something fresh.
segment('filecloseall_Callback',true);

%Stop the silent mode.
DATA.Silent = false;

%-------------------------------------------
function exportmultipleROI_Callback(dosegdicom) %#ok<DEFNU>
%-------------------------------------------
%Creaty summary of multiple matfiles in one folder
%This function is very useful for research. The user 
%performs all delineations and then exports all data to
%one spread sheet.

global DATA SET NO

if nargin < 1
  dosegdicom = false;
end
if dosegdicom
  suffix = 'segdicom';
else
  suffix = 'mat';
end

%Select path
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,sprintf('Select a folder with .%s files',suffix));
if isequal(pathname,0)
  myfailed('Aborted.',DATA.GUI.Segment);
  return;
end;

%Find files to process
files2load = dir([pathname filesep sprintf('*.%s',suffix)]);
numfiles = length(files2load);

if numfiles==0
  myfailed('Found no files to summarize.',DATA.GUI.Segment);
  return;
end;
% 
% includenormalized = yesno('Do you want to include BSA normalized values?');

%Create output matrix
%outdata = cell(numfiles+1,60); %+1 since header, 58 since header size %SBT20160620: Changed colsize from 58 to 60
%[temp,indforbsa] = header(false); %False means not only one image stack

%--- Write header
%for loop=1:length(temp)
%  outdata{1,loop} = temp{loop};
%end;

% if ~onlyone
%   outdata{1, 1} = 'FileName';
% else
%   outdata{1, 1} = 'ImageStack';
% end;
% outdata{1, 2} = 'PatientName';
% outdata{1, 3} = 'PatientID';
% outdata{1, 4} = 'AcquisitionDate';
% outdata{1, 5} = 'Age';
% outdata{1, 6} = 'Height[cm]';
% outdata{1, 7} = 'Weight[kg]';
% outdata{1, 8} = 'Sex';
% outdata{1, 9} = 'BSA[m2]';
% outdata{1,10} = 'HeartRate[bpm]';
outdata{1, 1}= 'FileName';
outdata{1, 2}= 'Stack number';
outdata{1, 3}= 'ImageViewPlane';
outdata{1, 4}= 'RoiName';
outdata{1, 5}= 'Area';
outdata{1, 6}= 'Mean';
outdata{1, 7}= 'StD';
outdata{1, 8}= 'Flow';


%Loop over all files
numrows=1;
h = mywaitbarstart(numfiles,'Please wait, loading and summarizing files.',1);
for fileloop=1:numfiles

  %--- Load file
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
  
  disp(dprintf('Loading %s.',files2load(fileloop).name)); %#ok<DSPS>
  
  %Set filename
  outdata{numrows+1,1} = files2load(fileloop).name;
  
  SET = []; % %Make sure a fresh start
  
  %---- try
  % This try-catch clause used to be commented out, which caused errors
  % in one file to crash the whole export process. Reintroducing it. /NL
  try
    %Load
    SET=[]; 
    if ~dosegdicom
      load([pathname filesep files2load(fileloop).name],'-mat');
    else
      loader = segloader();
      loader.adddicomfiles({fullfile(pathname,files2load(fileloop).name)});
      [~, setstruct] = loader.render(pathname, []);      
      clear loader;
    end
    
    %Assign
    SET = setstruct;
    clear setstruct;
    
    %traverse SET to get all rois
    for no=1:length(SET)
      Rois=SET(no).Roi; 
      %We want them in alphabetical order
      for i=1:length(Rois)
        names{i,1}=Rois.Name;
      end
      [~,sortind]=sortrows(names);
      names={};
      Rois=Rois(sortind);
      if ~(length(Rois)==1 && isempty(Rois(1).X)) 
        outdata{numrows+1,2}=no;
        outdata{numrows+1,3}=SET(no).ImageViewPlane;
        for roi=Rois
          outdata{numrows+1,4}=roi.Name;
          outdata{numrows+1,5}=nanmean(roi.Area);
          outdata{numrows+1,6}=nanmean(roi.Mean);
          outdata{numrows+1,7}=nanmean(roi.StD);
          if ~isempty(roi.Flow)
            outdata{numrows+1,8}=nanmean(roi.Flow.meanflow);
          end
          numrows=numrows+1;
        end
      end
      %numrows=numrows+3;
    end
    
    %NO = 1;
%     
%     %Call to intialize all variables correcly after loaded data.
%     openfile('setupstacksfrommat',1);
%     segment('renderstacksfrommat');
%     
% 		
%		numrows=numrows+size(temp,1);
    
  catch me
    %--- Some thing went wrong
    mydispexception(me);
    outdata{numrows+1,2} = 'FAILED.';
  end
  h = mywaitbarupdate(h);
end; %loop over files
mywaitbarclose(h);

%If not wanted remove non normalized values.
%ind = true(1,size(outdata,2));
%if ~includenormalized
%  ind(indforbsa) = false;
%end;
%outdata = outdata(:,ind);

%--- Output to a string
segment('cell2clipboard',outdata);

%Make sure starting with something fresh.
segment('filecloseall_Callback',true);

%Stop the silent mode.
DATA.Silent = false;

%-----------------------------------
function exportmultipleinfo_Callback %#ok<DEFNU>
%-----------------------------------
%Exports information of image stacks for a folder of mat files.
%This function is useful for debugging and checking purposes of
%the integrity of multiple .mat files.

global DATA SET

%Select path
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder with .mat files');
if isequal(pathname,0)
  myfailed('Aborted.',DATA.GUI.Segment);
  return;
end;

%Find files to process
files2load = dir([pathname filesep '*.mat']);
numfiles = length(files2load);

if numfiles==0
  myfailed('Found no files to export information about.',DATA.GUI.Segment);
  return;
end;

%Create output matrix
nbrdataperstack = 18;
outdata = cell(numfiles+1,nbrdataperstack*6); %6 is just a guess, 18 is data per image stack
maxno = 1;

%Loop over all files
h = mywaitbarstart(numfiles,'Please wait, loading and summarizing files.',1);
for fileloop=1:numfiles

  %--- Load file
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
  
  disp(dprintf('Loading %s.',files2load(fileloop).name)); %#ok<DSPS>
  
  %Set filename
  outdata{fileloop+1,1} = files2load(fileloop).name;
  
  SET = []; %#ok<NASGU> %Make sure a fresh start
  
  %---- try 
  %try
    %Load
    SET=[]; %#ok<NASGU>
    load([pathname filesep files2load(fileloop).name],'-mat');
    
    %Assign
    SET = setstruct;
    clear setstruct;
    
    %Call to intialize all variables correcly after loaded data.
    openfile('setupstacksfrommat',1);
    segment('renderstacksfrommat');
    
    maxno = max(length(SET),maxno); %Find maximum number of image stacks
    
    for sloop=1:length(SET)
      if isempty(SET(sloop).Scar)
        isscar = 'NO';
      else
        isscar = 'YES';
      end;
      if isempty(SET(sloop).Flow)
        isflow = 'NO';
      else
        isflow = 'YES';
      end;
      if isempty(SET(sloop).EndoX)&&isempty(SET(sloop).EpiX)
        isseg = 'NO';
      else
        isseg = 'YES';
      end;      
      coloffset = 1+(sloop-1)*nbrdataperstack;
      row = fileloop+1;
      outdata{row, 1+coloffset} = SET(sloop).ImageType;
      outdata{row, 2+coloffset} = SET(sloop).XSize;
      outdata{row, 3+coloffset} = SET(sloop).YSize;
      outdata{row, 4+coloffset} = SET(sloop).ZSize;
      outdata{row, 5+coloffset} = SET(sloop).TSize;
      outdata{row, 6+coloffset} = SET(sloop).ResolutionX;
      outdata{row, 7+coloffset} = SET(sloop).ResolutionY;
      outdata{row, 8+coloffset} = SET(sloop).SliceThickness;
      outdata{row, 9+coloffset} = SET(sloop).SliceGap;
      outdata{row,10+coloffset} = SET(sloop).TIncr;
      outdata{row,11+coloffset} = isscar;
      outdata{row,12+coloffset} = isflow;      
      outdata{row,13+coloffset} = isseg;
      outdata{row,14+coloffset} = SET(sloop).FlipAngle;
      outdata{row,15+coloffset} = SET(sloop).EchoTime;
      outdata{row,16+coloffset} = SET(sloop).RepetitionTime;
      outdata{row,17+coloffset} = SET(sloop).InversionTime;
      outdata{row,18+coloffset} = SET(sloop).HeartRate;
    end;
    
  %catch
  %  %--- Some thing went wrong
  %  outdata{fileloop+1,2} = 'FAILED.';
  %end
  h = mywaitbarupdate(h);
end; %loop over files
mywaitbarclose(h);

%--- Write header
for sloop=1:maxno
  coloffset = 1+(sloop-1)*nbrdataperstack;
  row = 1; %First line is header line  
  outdata{row, 1+coloffset} = 'ImageType';
  outdata{row, 2+coloffset} = 'XSize';
  outdata{row, 3+coloffset} = 'YSize';
  outdata{row, 4+coloffset} = 'ZSize';
  outdata{row, 5+coloffset} = 'TSize';
  outdata{row, 6+coloffset} = 'ResX';
  outdata{row, 7+coloffset} = 'ResY';
  outdata{row, 8+coloffset} = 'SliceThickness';
  outdata{row, 9+coloffset} = 'SliceGap';
  outdata{row,10+coloffset} = 'TIncr';
  outdata{row,11+coloffset} = 'ScarData';
  outdata{row,12+coloffset} = 'FlowData';
  outdata{row,13+coloffset} = 'Segmentation';
  outdata{row,14+coloffset} = 'FlipAngle';
  outdata{row,15+coloffset} = 'EchoTime';
  outdata{row,16+coloffset} = 'RepetitionTime';
  outdata{row,17+coloffset} = 'InversionTime';
  outdata{row,18+coloffset} = 'HeartRate';
end;

%--- Output to a string
segment('cell2clipboard',outdata);

%Make sure starting with something fresh.
segment('filecloseall_Callback',true);

%Stop the silent mode.
DATA.Silent = false;

%-------------------------------------------
function exportmultiplestraingraphs_Callback(type) %#ok<DEFNU>
%-------------------------------------------
%Creaty summary of strain result (tagging or cine strain analysis)
%from multiple matfiles in one folder.
%This function is very useful for research. The user 
%performs all strain analysis and then exports all data to
%one spreadsheet.

global DATA SET NO

suffix = 'mat';

%Ask if wnat to save before closing current image stack
if ~isempty(SET)
  if yesno('Would you like to store current open file before closing it?')
    %store file
    segment('filesaveallas_Callback');
  end
  %close current file
  segment('filecloseall_Callback');
end

%Select path
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,sprintf('Select a folder with .%s files',suffix));
if isequal(pathname,0)
  myfailed('Aborted.',DATA.GUI.Segment);
  return;
end;

%Find files to process
files2load = dir([pathname filesep sprintf('*.%s',suffix)]);
numfiles = length(files2load);

if numfiles==0
  myfailed('Found no files to summarize.',DATA.GUI.Segment);
  return;
end;

% do strain?
output=questdlg('Do you wish to redo strain analysis?');
switch output
  case 'Yes'
    doStrain=1;
  case 'No'
    doStrain=0;
  case 'Cancel'
    return;
end

currentsilent = DATA.Silent;
% %Create output matrix
outdata = cell(1,1); %+1 since header, 58 since header size

%Loop over all files
h = mywaitbarstart(numfiles,'Please wait, loading and summarizing files.',1);

% we need to do some roboclicking!
import java.awt.*;
import java.awt.event.*;
%Create a Robot-object to do the key-pressing
rob=Robot;

line=2;
  
for fileloop=1:numfiles
  %--- Load file
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
  disp(dprintf('Loading %s.',files2load(fileloop).name));
  
  SET = []; % %Make sure a fresh start
  load([pathname filesep files2load(fileloop).name],'-mat');
  %Assign
  SET = setstruct;
  clear setstruct;
  
  openfile('setupstacksfrommat',1);
  segment('renderstacksfrommat');
  
  %For long axis its sufficient to consider one of the images to obtain
  %all information. therefore checked exists
  checked=[];
    
  bullseyeradial = cell(1,1);
  bullseyecirc = cell(1,1);
  for no=1:length(SET)
    if doStrain
      if ~isfield(SET(no),'StrainTagging') || isempty(SET(no).StrainTagging) ||~isfield(SET(no).StrainTagging,'globalrad') %|| ismember(no,checked)
        %Skip this no
      else
        try
          %Do Strain calculations
          switch SET(no).ImageViewPlane
            case 'Short-axis'
              imageviewplane = 'shortaxis';
            case {'2CH','3CH','4CH'}
              imageviewplane = 'longaxis';              
            otherwise
              disp('Unknown image view plane');
          end
          
          switch SET(no).ImageType;
            case {'Cine', 'Feature tracking'}
              imagetype='cine';
            case 'Strain from tagging'
              imagetype='tagging';
          end
          
          if isequal(type,imagetype)
            NO=no;
            SET(no).StrainTagging.LVupdated=1;
            straintagging.straintagging('init',imagetype,imageviewplane);
            disp(['Performed strain analysis on no=',num2str(no)]);
            %store bullseye strain values
            gui = DATA.GUI.StrainTagging;
            set(gui.handles.bullseyepopupmenu,'value',1);
            straintagging.straintagging('updatebullseye');
            bullseyecirc{no} = gui.bullseyeplot;
            set(gui.handles.bullseyepopupmenu,'value',2);
            straintagging.straintagging('updatebullseye');
            bullseyeradial{no} = gui.bullseyeplot;
            straintagging.straintagging('close_Callback');
          end
        catch
          disp(['Failed to redo strain analysis, skipping no=',num2str(no)]);
        end
        if isfield(SET(NO).StrainTagging,'taggroup')
          checked=[checked,SET(NO).StrainTagging.taggroup];
        end
      end
    end %end of do calc strain
  end %end of loop over image stacks
    
  include2ch=0;
  include3ch=0;
  include4ch=0;
  Tlax=1;
  for no=1:length(SET)
    if ~isfield(SET(no),'StrainTagging')|| isempty(SET(no).StrainTagging) || ~isfield(SET(no).StrainTagging,'globalrad') %|| ismember(no,checked)
      %Skip no
    else
      
      switch SET(no).ImageType;
        case {'Cine', 'Feature tracking'}
          imagetype='cine';
        case 'Strain from tagging'
          imagetype='tagging';
      end
      
      if isequal(type,imagetype)
        %store global strain values
         switch SET(no).ImageViewPlane
          case '2CH'
            tf2ch = SET(no).StrainTagging.peaktf;
            globalcirc2ch = mynanmean(SET(no).StrainTagging.globalcirc(tf2ch,:));
            globalrad2ch = mynanmean(SET(no).StrainTagging.globalrad(tf2ch,:));
            segmentalcirc2ch = SET(no).StrainTagging.segmentcirc;
            segmentalrad2ch = SET(no).StrainTagging.segmentrad;
            Tlax=SET(no).TSize;
            include2ch=1;
          case '3CH'
            tf3ch = SET(no).StrainTagging.peaktf;
            globalcirc3ch = mynanmean(SET(no).StrainTagging.globalcirc(tf3ch,:));
            globalrad3ch = mynanmean(SET(no).StrainTagging.globalrad(tf3ch,:));
            segmentalcirc3ch = SET(no).StrainTagging.segmentcirc;
            segmentalrad3ch = SET(no).StrainTagging.segmentrad;
            Tlax=SET(no).TSize;
            include3ch=1;
          case '4CH'
            tf4ch = SET(no).StrainTagging.peaktf;
            globalcirc4ch = mynanmean(SET(no).StrainTagging.globalcirc(tf4ch,:));
            globalrad4ch = mynanmean(SET(no).StrainTagging.globalrad(tf4ch,:));
            segmentalcirc4ch = SET(no).StrainTagging.segmentcirc;
            segmentalrad4ch = SET(no).StrainTagging.segmentrad;
            Tlax=SET(no).TSize;
            include4ch=1;
        end
      end
    end
  end
  %calculate global mean strain
  if ~include2ch
    segmentalcirc2ch = NaN*ones(Tlax,7);
    segmentalrad2ch = NaN*ones(Tlax,7);
  end
  
  if ~include3ch
    segmentalcirc3ch = NaN*ones(Tlax,7);
    segmentalrad3ch = NaN*ones(Tlax,7);
  end
  
  if ~include4ch
    segmentalcirc4ch = NaN*ones(Tlax,7);
    segmentalrad4ch = NaN*ones(Tlax,7);
  end
  
  bullseyecirclax=cell(1,Tlax);
  bullseyeradiallax=cell(1,Tlax);
  for t=1:Tlax
    bullseyecirclax{t} = straintagging.straintagging('getbullseyevalueslax',segmentalcirc2ch,segmentalcirc3ch,segmentalcirc4ch,t,t,t);
    bullseyeradiallax{t} = straintagging.straintagging('getbullseyevalueslax',segmentalrad2ch,segmentalrad3ch,segmentalrad4ch,t,t,t);
  end
  %   line=2;
%   %output
checked= [];
  for no=1:length(SET)    
    if ~isfield(SET(no),'StrainTagging')|| isempty(SET(no).StrainTagging) || ~isfield(SET(no).StrainTagging,'globalrad') || ismember(no,checked) 
      %Skip no
    else
          %--- output ---
      %patient info
      outdata{line,1} = 'Patient name';
      outdata{line+1,1} = 'Patient ID';
      outdata{line+2,1} = 'Heart Rate';
      outdata{line+3,1} = 'Image Type';
      
      switch SET(no).ImageType;
        case {'Cine', 'Feature tracking'}
          imagetype='cine';
        case 'Strain from tagging'
          imagetype='tagging';
      end
      
      
      radcol = 1;
      circcol = SET(no).TSize+3;
      %aha sections
      for loop=1:17
        [stri,pos] = reportbullseye('aha17nameandpos',loop); %Get name and position of export
        outdata{line+8+loop,radcol} = stri;
        outdata{line+8+loop,circcol} = stri;
      end
       outdata{line+7,1} = 'Radial strain';
     if strcmp(SET(no).ImageViewPlane,'Short-axis')
        outdata{line+7,circcol} = 'Circumferential strain';
     else
        outdata{line+7,circcol} = 'Longitudinal strain';
     end
    
 %     if isequal(type,imagetype)
        outdata{line,2} = SET(no).PatientInfo.Name;
        outdata{line+1,2} = SET(no).PatientInfo.ID;
        outdata{line+2,2} = SET(no).HeartRate;
        if strcmp(SET(no).ImageViewPlane,'Short-axis')
          outdata{line+3,2} = sprintf('%s %s',SET(no).ImageType,SET(no).ImageViewPlane);
        else
          %currently used slices
          cus=[];
          for tagno=SET(no).StrainTagging.taggroup;
            cus=[cus, ' ', SET(tagno).ImageViewPlane];
          end
          outdata{line+3,2} = sprintf('%s %s',SET(no).ImageType,cus);
        end
        outdata{line+5,radcol}='Time [s]';
        outdata{line+5,circcol}='Time [s]';
        for t=1:SET(no).TSize
          outdata{line+5,radcol+t}=SET(no).TimeVector(t);
          outdata{line+5,circcol+t}=SET(no).TimeVector(t);
        end
        line=line+9;
        switch SET(no).ImageViewPlane
          case 'Short-axis'
            %get output values for bullseye
            for t=1:SET(no).TSize
              bullseyeplot = straintagging.straintagging('getbullseyevaluessax',SET(no).StrainTagging.segmentrad,t,SET(no).StrainTagging.saslices);
            %basal
            outdata{line,radcol+t} = mynanmean(bullseyeplot(1:4,4));
            outdata{line+1,radcol+t} = mynanmean(bullseyeplot(5:8,4));
            outdata{line+2,radcol+t} = mynanmean(bullseyeplot(9:12,4));
            outdata{line+3,radcol+t} = mynanmean(bullseyeplot(13:16,4));
            outdata{line+4,radcol+t} = mynanmean(bullseyeplot(17:20,4));
            outdata{line+5,radcol+t} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line+6,radcol+t} = mynanmean(bullseyeplot(1:4,3));
            outdata{line+7,radcol+t} = mynanmean(bullseyeplot(5:8,3));
            outdata{line+8,radcol+t} = mynanmean(bullseyeplot(9:12,3));
            outdata{line+9,radcol+t} = mynanmean(bullseyeplot(13:16,3));
            outdata{line+10,radcol+t} = mynanmean(bullseyeplot(17:20,3));
            outdata{line+11,radcol+t} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line+12,radcol+t} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line+13,radcol+t} = mynanmean(bullseyeplot(3:8,2));
            outdata{line+14,radcol+t} = mynanmean(bullseyeplot(9:14,2));
            outdata{line+15,radcol+t} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line+16,radcol+t} = mynanmean(bullseyeplot(:,1));
            
            bullseyeplot = straintagging.straintagging('getbullseyevaluessax',SET(no).StrainTagging.segmentcirc,t,SET(no).StrainTagging.saslices);
            %basal
            outdata{line,circcol+t} = mynanmean(bullseyeplot(1:4,4));
            outdata{line+1,circcol+t} = mynanmean(bullseyeplot(5:8,4));
            outdata{line+2,circcol+t} = mynanmean(bullseyeplot(9:12,4));
            outdata{line+3,circcol+t} = mynanmean(bullseyeplot(13:16,4));
            outdata{line+4,circcol+t} = mynanmean(bullseyeplot(17:20,4));
            outdata{line+5,circcol+t} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line+6,circcol+t} = mynanmean(bullseyeplot(1:4,3));
            outdata{line+7,circcol+t} = mynanmean(bullseyeplot(5:8,3));
            outdata{line+8,circcol+t} = mynanmean(bullseyeplot(9:12,3));
            outdata{line+9,circcol+t} = mynanmean(bullseyeplot(13:16,3));
            outdata{line+10,circcol+t} = mynanmean(bullseyeplot(17:20,3));
            outdata{line+11,circcol+t} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line+12,circcol+t} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line+13,circcol+t} = mynanmean(bullseyeplot(3:8,2));
            outdata{line+14,circcol+t} = mynanmean(bullseyeplot(9:14,2));
            outdata{line+15,circcol+t} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line+16,circcol+t} = mynanmean(bullseyeplot(:,1));
            end
            
            line=line+18;

          case {'2CH','3CH','4CH'}
            %get output values for bullseye
            %counterclockwise from 12 o'clock if you are looking at the plot
            for t=1:SET(no).TSize
            bullseyeplot = bullseyeradiallax{t};
            %basal
            outdata{line,radcol+t} = mynanmean(bullseyeplot(1:4,4));
            outdata{line+1,radcol+t} = mynanmean(bullseyeplot(5:8,4));
            outdata{line+2,radcol+t} = mynanmean(bullseyeplot(9:12,4));
            outdata{line+3,radcol+t} = mynanmean(bullseyeplot(13:16,4));
            outdata{line+4,radcol+t} = mynanmean(bullseyeplot(17:20,4));
            outdata{line+5,radcol+t} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line+6,radcol+t} = mynanmean(bullseyeplot(1:4,3));
            outdata{line+7,radcol+t} = mynanmean(bullseyeplot(5:8,3));
            outdata{line+8,radcol+t} = mynanmean(bullseyeplot(9:12,3));
            outdata{line+9,radcol+t} = mynanmean(bullseyeplot(13:16,3));
            outdata{line+10,radcol+t} = mynanmean(bullseyeplot(17:20,3));
            outdata{line+11,radcol+t} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line+12,radcol+t} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line+13,radcol+t} = mynanmean(bullseyeplot(3:8,2));
            outdata{line+14,radcol+t} = mynanmean(bullseyeplot(9:14,2));
            outdata{line+15,radcol+t} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line+16,radcol+t} = mynanmean(bullseyeplot(:,1));
            
            bullseyeplot = bullseyecirclax{t};
            %basal
            outdata{line,circcol+t} = mynanmean(bullseyeplot(1:4,4));
            outdata{line+1,circcol+t} = mynanmean(bullseyeplot(5:8,4));
            outdata{line+2,circcol+t} = mynanmean(bullseyeplot(9:12,4));
            outdata{line+3,circcol+t} = mynanmean(bullseyeplot(13:16,4));
            outdata{line+4,circcol+t} = mynanmean(bullseyeplot(17:20,4));
            outdata{line+5,circcol+t} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line+6,circcol+t} = mynanmean(bullseyeplot(1:4,3));
            outdata{line+7,circcol+t} = mynanmean(bullseyeplot(5:8,3));
            outdata{line+8,circcol+t} = mynanmean(bullseyeplot(9:12,3));
            outdata{line+9,circcol+t} = mynanmean(bullseyeplot(13:16,3));
            outdata{line+10,circcol+t} = mynanmean(bullseyeplot(17:20,3));
            outdata{line+11,circcol+t} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line+12,circcol+t} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line+13,circcol+t} = mynanmean(bullseyeplot(3:8,2));
            outdata{line+14,circcol+t} = mynanmean(bullseyeplot(9:14,2));
            outdata{line+15,circcol+t} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line+16,circcol+t} = mynanmean(bullseyeplot(:,1));
            
            end
           line=line+18;
           
            checked=[checked,SET(no).StrainTagging.taggroup];  
        end
        
        %line = line+1;
%       end
    end
  end; %loop over image stack
  if doStrain
    %store file
    
    %Create thumbnails before storing.
    calcfunctions('calcdatasetpreview');
    
    % Set view settings
    DATA.ViewPanels = 1;
    DATA.ViewPanelsType = {'one'};
    DATA.ViewMatrix = [1 1];
    DATA.ThisFrameOnly = 0;
    DATA.CurrentPanel = 1;
    DATA.CurrentTheme = 'lv';
    DATA.CurrentTool = 'select';

    %Save the file.
    %segment('filesaveall_Callback');
    corrupted=segment('checkcorrupteddataforautomaticsave');
    if corrupted
      corruptedfiles=sprintf('%s, %s',corruptedfiles,filename);
      mywarning(dprintf('Image file %s seems to be corrupted from last save. Please load and manually re-analyse strain to ensure that the image is not corrupted before saving',filename));
    else
      filemenu('saveallas_helper',pathname,files2load(fileloop).name);
      disp(sprintf('Saving %s',[pathname filesep files2load(fileloop).name]));
    end    
    
    rob.keyPress(KeyEvent.VK_SPACE);
    pause(0.1);
    rob.keyRelease(KeyEvent.VK_SPACE);
  end
  h = mywaitbarupdate(h);  
  %close current file
  segment('filecloseall_Callback');
end %loop over files
mywaitbarclose(h);

%--- Output to a string
segment('cell2clipboard',outdata);
%Stop the silent mode.
DATA.Silent = currentsilent;
  

%-------------------------------------------
function exportmultiplestrain_Callback(type) %#ok<DEFNU>
%-------------------------------------------
%Creaty summary of strain result (tagging or cine strain analysis)
%from multiple matfiles in one folder.
%This function is very useful for research. The user 
%performs all strain analysis and then exports all data to
%one spreadsheet.

global DATA SET NO

suffix = 'mat';

%Ask if wnat to save before closing current image stack
if ~isempty(SET)
  if yesno('Would you like to store current open file before closing it?')
    %store file
    segment('filesaveallas_Callback');
  end
  %close current file
  segment('filecloseall_Callback');
end

%Select path
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,sprintf('Select a folder with .%s files',suffix));
if isequal(pathname,0)
  myfailed('Aborted.',DATA.GUI.Segment);
  return;
end;

%Find files to process
files2load = dir([pathname filesep sprintf('*.%s',suffix)]);
numfiles = length(files2load);

if numfiles==0
  myfailed('Found no files to summarize.',DATA.GUI.Segment);
  return;
end;

% do strain?
output=questdlg('Do you wish to redo strain analysis?');
switch output
  case 'Yes'
    doStrain=1;
  case 'No'
    doStrain=0;
  case 'Cancel'
    return;
end


currentsilent = DATA.Silent;
% %Create output matrix
outdata = cell(1,1); %+1 since header, 58 since header size

%Loop over all files
h = mywaitbarstart(numfiles,'Please wait, loading and summarizing files.',1);

% we need to do some roboclicking!
import java.awt.*;
import java.awt.event.*;
%Create a Robot-object to do the key-pressing
rob=Robot;

% %--- output ---
% %patient info
% outdata{2,1} = 'Patient name';
% outdata{2,2} = 'Patient ID';
% outdata{2,3} = 'Heart Rate';
% outdata{2,4} = 'Image Type';
% %global strain values
% outdata{2,5} = 'Peak mean radial strain [%]';
% outdata{2,6} = 'Peak mean circ./longit. strain [%]';
% outdata{2,7} = 'Peak Time [ms]';
% outdata{2,8} = 'Peak mean radial strain [%] (entire LV)';
% outdata{2,9} = 'Peak mean circ./longit. strain [%] (entire LV)';
% %segmental strain values from bullseye
% outdata{1,10} = 'Segmental peak radial strain [%]';
% outdata{1,27} = 'Segmental peak circ./longit. strain [%]';
% %aha sections
% for loop=1:17
%   [stri,pos] = reportbullseye('aha17nameandpos',loop); %Get name and position of export
%   outdata{2,9+pos} = stri;
%   outdata{2,26+pos} = stri;
% end
% line = 3;

line=2;
  
for fileloop=1:numfiles
  %--- Load file
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
  disp(dprintf('Loading %s.',files2load(fileloop).name));
  
  SET = []; % %Make sure a fresh start
  load([pathname filesep files2load(fileloop).name],'-mat');
  %Assign
  SET = setstruct;
  clear setstruct;
  
  openfile('setupstacksfrommat',1);
  segment('renderstacksfrommat');
  
  %For long axis its sufficient to consider one of the images to obtain
  %all information. therefore checked exists
  checked=[];
    
  bullseyeradial = cell(1,1);
  bullseyecirc = cell(1,1);
  for no=1:length(SET)
    if doStrain
      if ~isfield(SET(no),'StrainTagging') || isempty(SET(no).StrainTagging) ||~isfield(SET(no).StrainTagging,'globalrad') %|| ismember(no,checked)
        %Skip this no
      else
        try
          %Do Strain calculations
          switch SET(no).ImageViewPlane
            case 'Short-axis'
              imageviewplane = 'shortaxis';
            case {'2CH','3CH','4CH'}
              imageviewplane = 'longaxis';              
            otherwise
              disp('Unknown image view plane');
          end
          
          switch SET(no).ImageType;
            case {'Cine', 'Feature tracking'}
              imagetype='cine';
            case 'Strain from tagging'
              imagetype='tagging';
          end
          
          if isequal(type,imagetype)
            NO=no;
            SET(no).StrainTagging.LVupdated=1;
            straintagging.straintagging('init',imagetype,imageviewplane);
            disp(['Performed strain analysis on no=',num2str(no)]);
            %store bullseye strain values
            gui = DATA.GUI.StrainTagging;
            set(gui.handles.bullseyepopupmenu,'value',1);
            straintagging.straintagging('updatebullseye');
            bullseyecirc{no} = gui.bullseyeplot;
            set(gui.handles.bullseyepopupmenu,'value',2);
            straintagging.straintagging('updatebullseye');
            bullseyeradial{no} = gui.bullseyeplot;
            straintagging.straintagging('close_Callback');
          end
        catch
          disp(['Failed to redo strain analysis, skipping no=',num2str(no)]);
        end
        if isfield(SET(NO).StrainTagging,'taggroup')
          checked=[checked,SET(NO).StrainTagging.taggroup];
        end
      end
    end %end of do calc strain
  end %end of loop over image stacks
    
  %calculate global mean strain
  globalcirc2ch = NaN;
  globalrad2ch = NaN;
  globalcirc3ch = NaN;
  globalrad3ch = NaN;
  globalcirc4ch = NaN;
  globalrad4ch = NaN;
  segmentalcirc2ch = NaN*ones(1,7);
  segmentalrad2ch = NaN*ones(1,7);
  segmentalcirc3ch = NaN*ones(1,7);
  segmentalrad3ch = NaN*ones(1,7);
  segmentalcirc4ch = NaN*ones(1,7);
  segmentalrad4ch = NaN*ones(1,7);
  tf2ch = 1;
  tf3ch = 1;
  tf4ch = 1;
%   no2ch=[];
%   no3ch=[];
%   no4ch=[];
  
  for no=1:length(SET)
    if ~isfield(SET(no),'StrainTagging')|| isempty(SET(no).StrainTagging) || ~isfield(SET(no).StrainTagging,'globalrad') %|| ismember(no,checked)
      %Skip no
    else
      
      switch SET(no).ImageType;
        case {'Cine', 'Feature tracking'}
          imagetype='cine';
        case 'Strain from tagging'
          imagetype='tagging';
      end
      
      if isequal(type,imagetype)
        %store global strain values
         switch SET(no).ImageViewPlane
          case '2CH'
            tf2ch = SET(no).StrainTagging.peaktf;
            globalcirc2ch = mynanmean(SET(no).StrainTagging.globalcirc(tf2ch,:));
            globalrad2ch = mynanmean(SET(no).StrainTagging.globalrad(tf2ch,:));
            segmentalcirc2ch = SET(no).StrainTagging.segmentcirc;
            segmentalrad2ch = SET(no).StrainTagging.segmentrad;
%             if isfield(SET(no).StrainTagging,'strainratecircum') && ~isempty(SET(no).StrainTagging.strainratecircum)
%               segmentalcircSR2CH = SET(no).StrainTagging.strainratecircum;
%               [~, tfcirc2chdown]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.downslopecircum(:,2)/1000)));
%               [~, tfcirc2chup]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.upslopecircum(:,2)/1000)));
%               
%               segmentalradSR2CH = SET(no).StrainTagging.strainraterad;
%               [~, tfrad2chdown]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.downsloperad(:,2)/1000)));
%               [~, tfrad2chup]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.upsloperad(:,2)/1000)));
%               no2ch=no;
%             else
%               segmentalcircSR2CH = nan(SET(no).ZSize,1);
%               tfcirc2chup=1;
%               tfcirc2chdown=1;
%               
%                segmentalradSR2CH = nan(SET(no).ZSize,1);
%                tfrad2chup=1;
%                tfrad2chdown=1;
%             end
           case '3CH'
             tf3ch = SET(no).StrainTagging.peaktf;
             globalcirc3ch = mynanmean(SET(no).StrainTagging.globalcirc(tf3ch,:));
             globalrad3ch = mynanmean(SET(no).StrainTagging.globalrad(tf3ch,:));
             segmentalcirc3ch = SET(no).StrainTagging.segmentcirc;
             segmentalrad3ch = SET(no).StrainTagging.segmentrad;
%              if isfield(SET(no).StrainTagging,'strainratecircum') && ~isempty(SET(no).StrainTagging.strainratecircum)
%                segmentalcircSR3CH = SET(no).StrainTagging.strainratecircum;
%                [~, tfcirc3chdown]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.downslopecircum(:,2)/1000)));
%                [~, tfcirc3chup]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.upslopecircum(:,2)/1000)));
%               
%                segmentalradSR3CH = SET(no).StrainTagging.strainraterad;
%               [~, tfrad3chdown]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.downsloperad(:,2)/1000)));
%               [~, tfrad3chup]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.upsloperad(:,2)/1000)));
%               
%                no3ch=no;
%              else
%                segmentalcircSR3CH = nan(SET(no).ZSize,1);
%                tfcirc3chup=1;
%                tfcirc3chdown=1;
%                segmentalradSR3CH = nan(SET(no).ZSize,1);
%                tfrad3chup=1;
%                tfrad3chdown=1;
%              
%              end
           case '4CH'
             tf4ch = SET(no).StrainTagging.peaktf;
             globalcirc4ch = mynanmean(SET(no).StrainTagging.globalcirc(tf4ch,:));
             globalrad4ch = mynanmean(SET(no).StrainTagging.globalrad(tf4ch,:));
             segmentalcirc4ch = SET(no).StrainTagging.segmentcirc;
             segmentalrad4ch = SET(no).StrainTagging.segmentrad;
%              if isfield(SET(no).StrainTagging,'strainratecircum') && ~isempty(SET(no).StrainTagging.strainratecircum)
%                segmentalcircSR4CH = SET(no).StrainTagging.strainratecircum;
%                [~, tfcirc4chdown]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.downslopecircum(:,2)/1000)));
%                [~, tfcirc4chup]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.upslopecircum(:,2)/1000)));
%                
%                segmentalradSR4CH = SET(no).StrainTagging.strainraterad;
%                [~, tfrad4chdown]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.downsloperad(:,2)/1000)));
%                [~, tfrad4chup]=min(abs(SET(no).StrainTagging.strainrateTvect-mean(SET(no).StrainTagging.upsloperad(:,2)/1000)));
%                no4ch=no;
%              else
%                segmentalcircSR4CH = nan(SET(no).ZSize,1);
%                tfcirc4chup=1;
%                tfcirc4chdown=1;
%                segmentalradSR4CH = nan(SET(no).ZSize,1);
%                tfrad4chup=1;
%                tfrad4chdown=1;
%              end
         end
      end
    end
  end
  bullseyecirclax = straintagging.straintagging('getbullseyevalueslax',segmentalcirc2ch,segmentalcirc3ch,segmentalcirc4ch,tf2ch,tf3ch,tf4ch);
  bullseyeradiallax = straintagging.straintagging('getbullseyevalueslax',segmentalrad2ch,segmentalrad3ch,segmentalrad4ch,tf2ch,tf3ch,tf4ch);
  
%   bullseyeSRcirclaxdown=nan(1,2);
%   bullseyeSRcirclaxup=nan(1,2);
%   bullseyeSRradlaxdown=nan(1,2);
%   bullseyeSRradlaxup=nan(1,2);
%   
%   bullseyeSRcirclaxdown = straintagging.straintagging('getbullseyevalueslax',segmentalcircSR2CH,segmentalcircSR3CH,segmentalcircSR4CH,...
%     tfcirc2chdown,tfcirc3chdown,tfcirc4chdown);
%   bullseyeSRcirclaxup = straintagging.straintagging('getbullseyevalueslax',segmentalcircSR2CH,segmentalcircSR3CH,segmentalcircSR4CH,...
%     tfcirc2chup,tfcirc3chup,tfcirc4chup);
%   
%   bullseyeSRradlaxdown = straintagging.straintagging('getbullseyevalueslax',segmentalradSR2CH,segmentalradSR3CH,segmentalradSR4CH,...
%     tfrad2chdown,tfrad3chdown,tfrad4chdown);
%   bullseyeSRradlaxup = straintagging.straintagging('getbullseyevalueslax',segmentalradSR2CH,segmentalradSR3CH,segmentalradSR4CH,...
%     tfrad2chup,tfrad3chup,tfrad4chup);
%   line=2;
%   %output
  for no=1:length(SET)    
    if ~isfield(SET(no),'StrainTagging')|| isempty(SET(no).StrainTagging) || ~isfield(SET(no).StrainTagging,'globalrad') %|| ismember(no,checked) 
      %Skip no
    else
          %--- output ---
      %patient info
      outdata{line,1} = 'Patient name';
      outdata{line,2} = 'Patient ID';
      outdata{line,3} = 'Heart Rate';
      outdata{line,4} = 'Image Type';
      %global strain values
      outdata{line,5} = 'Peak mean radial strain [%]';
      outdata{line,6} = 'Peak mean circ./longit. strain [%]';
      outdata{line,7} = 'Peak Time [ms]';
      outdata{line,8} = 'Peak mean radial strain [%] (entire LV)';
      outdata{line,9} = 'Peak mean circ./longit. strain [%] (entire LV)';
      if isfield(SET(no).StrainTagging,'globalRVstrain') && ~isempty(SET(no).StrainTagging.globalRVstrain)
        outdata{line,10}='Peak mean circ./longit. strain [%] (entire RV)';
        %segmental strain values from bullseye
        outdata{line-1,11} = 'Segmental peak radial strain [%]';
        outdata{line-1,28} = 'Segmental peak circ./longit. strain [%]';
        %aha sections
        for loop=1:17
          [stri,pos] = reportbullseye('aha17nameandpos',loop); %Get name and position of export
          outdata{line,10+pos} = stri;
          outdata{line,27+pos} = stri;
        end
      else
        %segmental strain values from bullseye
        outdata{line-1,10} = 'Segmental peak radial strain [%]';
        outdata{line-1,27} = 'Segmental peak circ./longit. strain [%]';
        %aha sections
        for loop=1:17
          [stri,pos] = reportbullseye('aha17nameandpos',loop); %Get name and position of export
          outdata{line,9+pos} = stri;
          outdata{line,26+pos} = stri;
        end
      end
      line = line+1;

      switch SET(no).ImageType;
        case {'Cine', 'Feature tracking'}
          imagetype='cine';
        case 'Strain from tagging'
          imagetype='tagging';
      end
      
      %if isequal(type,imagetype)
        outdata{line,1} = SET(no).PatientInfo.Name;
        outdata{line,2} = SET(no).PatientInfo.ID;
        outdata{line,3} = SET(no).HeartRate;
        outdata{line,4} = sprintf('%s %s',SET(no).ImageType,SET(no).ImageViewPlane);
        
        %Global peak strain
        peaktf = SET(no).StrainTagging.peaktf;
        outdata{line,5} = mynanmean(SET(no).StrainTagging.globalrad(peaktf,:));
        outdata{line,6} = mynanmean(SET(no).StrainTagging.globalcirc(peaktf,:));
        outdata{line,7} = SET(no).TimeVector(peaktf)*1000;
        switch SET(no).ImageViewPlane
          case 'Short-axis'
            outdata{line,8} = mynanmean(SET(no).StrainTagging.globalrad(peaktf,:));
            outdata{line,9} = mynanmean(SET(no).StrainTagging.globalcirc(peaktf,:));
            
            if isfield(SET(no).StrainTagging,'globalRVstrain') && ~isempty(SET(no).StrainTagging.globalRVstrain)
              outdata{line,10} = mynanmean(SET(no).StrainTagging.globalRVstrain(1,peaktf,:));
              radcol=11;
              circcol=28;
              col = 45;
            else
              radcol = 10;
              circcol = 27;
              col = 44;
            end
              
            %get output values for bullseye
            bullseyeplot = straintagging.straintagging('getbullseyevaluessax',SET(no).StrainTagging.segmentrad,peaktf,SET(no).StrainTagging.saslices);
            
            %basal
            %radcol = 10;
            outdata{line,radcol} = mynanmean(bullseyeplot(1:4,4));
            outdata{line,radcol+1} = mynanmean(bullseyeplot(5:8,4));
            outdata{line,radcol+2} = mynanmean(bullseyeplot(9:12,4));
            outdata{line,radcol+3} = mynanmean(bullseyeplot(13:16,4));
            outdata{line,radcol+4} = mynanmean(bullseyeplot(17:20,4));
            outdata{line,radcol+5} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line,radcol+6} = mynanmean(bullseyeplot(1:4,3));
            outdata{line,radcol+7} = mynanmean(bullseyeplot(5:8,3));
            outdata{line,radcol+8} = mynanmean(bullseyeplot(9:12,3));
            outdata{line,radcol+9} = mynanmean(bullseyeplot(13:16,3));
            outdata{line,radcol+10} = mynanmean(bullseyeplot(17:20,3));
            outdata{line,radcol+11} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line,radcol+12} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line,radcol+13} = mynanmean(bullseyeplot(3:8,2));
            outdata{line,radcol+14} = mynanmean(bullseyeplot(9:14,2));
            outdata{line,radcol+15} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line,radcol+16} = mynanmean(bullseyeplot(:,1));
            
            bullseyeplot = straintagging.straintagging('getbullseyevaluessax',SET(no).StrainTagging.segmentcirc,peaktf,SET(no).StrainTagging.saslices);
            %basal
            %circcol = 27;
            outdata{line,circcol} = mynanmean(bullseyeplot(1:4,4));
            outdata{line,circcol+1} = mynanmean(bullseyeplot(5:8,4));
            outdata{line,circcol+2} = mynanmean(bullseyeplot(9:12,4));
            outdata{line,circcol+3} = mynanmean(bullseyeplot(13:16,4));
            outdata{line,circcol+4} = mynanmean(bullseyeplot(17:20,4));
            outdata{line,circcol+5} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line,circcol+6} = mynanmean(bullseyeplot(1:4,3));
            outdata{line,circcol+7} = mynanmean(bullseyeplot(5:8,3));
            outdata{line,circcol+8} = mynanmean(bullseyeplot(9:12,3));
            outdata{line,circcol+9} = mynanmean(bullseyeplot(13:16,3));
            outdata{line,circcol+10} = mynanmean(bullseyeplot(17:20,3));
            outdata{line,circcol+11} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line,circcol+12} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line,circcol+13} = mynanmean(bullseyeplot(3:8,2));
            outdata{line,circcol+14} = mynanmean(bullseyeplot(9:14,2));
            outdata{line,circcol+15} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line,circcol+16} = mynanmean(bullseyeplot(:,1));
            
            %slice based export of strain and strain rate
            %col = 44;
            nbrslices = length(SET(no).StrainTagging.saslices);
            outdata{line-2,col} = 'Peak radial strain [%]';
            for sliceloop = 1:nbrslices
              outdata{line-1,col} = sprintf('Slice %d',sliceloop);
              outdata{line,col} = SET(no).StrainTagging.globalrad(SET(no).StrainTagging.peaktf,sliceloop);
              col = col+1;
            end
            outdata{line-2,col} = 'Peak circ. strain [%]';
            for sliceloop = 1:nbrslices
              outdata{line-1,col} = sprintf('Slice %d',sliceloop);
              outdata{line,col} = SET(no).StrainTagging.globalcirc(SET(no).StrainTagging.peaktf,sliceloop);
              col = col+1;
            end
            
            if isfield(SET(no).StrainTagging,'globalRVstrain') && ~isempty(SET(no).StrainTagging.globalRVstrain)
              nbrslices = length(SET(no).StrainTagging.saslices);
              outdata{line-2,col} = 'Peak circ. strain RV [%]';
              for sliceloop = 1:nbrslices
                outdata{line-1,col} = sprintf('Slice %d',sliceloop);
                outdata{line,col} = SET(no).StrainTagging.globalRVstrain(1,SET(no).StrainTagging.peaktf,sliceloop);
                col = col+1;
              end
            end
            
            outdata{line-2,col} = 'Radial strain rate [%/s]';
            
            outdata{line-1,col} ='Mean Upslope';
            outdata{line,col}=nanmean(SET(no).StrainTagging.upsloperad(:,1));
            col=col+1;
            outdata{line-1,col} ='Mean Upslope time';
            outdata{line,col}=nanmean(SET(no).StrainTagging.upsloperad(:,2));
            col=col+1;
            outdata{line-1,col} ='Mean Downslope ';
            outdata{line,col}=nanmean(SET(no).StrainTagging.downsloperad(:,1));
            col=col+1;
            outdata{line-1,col} ='Mean Downslope time';
            outdata{line,col}=nanmean(SET(no).StrainTagging.downsloperad(:,2));
            col=col+1;
            for sliceloop = 1:nbrslices
              outdata{line-1,col} = sprintf('Slice %d Upslope',sliceloop);
              outdata{line,col} = SET(no).StrainTagging.upsloperad(sliceloop,1);
              outdata{line-1,col+1} = sprintf('Slice %d Upslope time',sliceloop);
              outdata{line,col+1} = SET(no).StrainTagging.upsloperad(sliceloop,2);
              outdata{line-1,col+2} = sprintf('Slice %d Downslope',sliceloop);
              outdata{line,col+2} = SET(no).StrainTagging.downsloperad(sliceloop,1);
              outdata{line-1,col+3} = sprintf('Slice %d Downslope time',sliceloop);
              outdata{line,col+3} = SET(no).StrainTagging.downsloperad(sliceloop,2);
              col = col+4;
            end
            outdata{line-2,col} = 'Circ. strain rate [%/s]';
            outdata{line-1,col} ='Mean Upslope';
            outdata{line,col}=nanmean(SET(no).StrainTagging.upslopecircum(:,1));
            col=col+1;
            outdata{line-1,col} ='Mean Upslope time';
            outdata{line,col}=nanmean(SET(no).StrainTagging.upslopecircum(:,2));
            col=col+1;
            outdata{line-1,col} ='Mean Downslope ';
            outdata{line,col}=nanmean(SET(no).StrainTagging.downslopecircum(:,1));
            col=col+1;
            outdata{line-1,col} ='Mean Downslope time';
            outdata{line,col}=nanmean(SET(no).StrainTagging.downslopecircum(:,2));
            col=col+1;
            
            for sliceloop = 1:nbrslices
              outdata{line-1,col} = sprintf('Slice %d Upslope',sliceloop);
              outdata{line,col} = SET(no).StrainTagging.upslopecircum(sliceloop,1);
              outdata{line-1,col+1} = sprintf('Slice %d Upslope time',sliceloop);  
              outdata{line,col+1} = SET(no).StrainTagging.upslopecircum(sliceloop,2);            
              outdata{line-1,col+2} = sprintf('Slice %d Downslope',sliceloop);
              outdata{line,col+2} = SET(no).StrainTagging.downslopecircum(sliceloop,1);
              outdata{line-1,col+3} = sprintf('Slice %d Downslope time',sliceloop);
              outdata{line,col+3} = SET(no).StrainTagging.downslopecircum(sliceloop,2);
              col = col+4;
            end
            
            line=line+3;

          case {'2CH','3CH','4CH'}
            outdata{line,8} = mynanmean([globalrad2ch globalrad3ch globalrad4ch]);
            outdata{line,9} = mynanmean([globalcirc2ch globalcirc3ch globalcirc4ch]);
            
             if isfield(SET(no).StrainTagging,'globalRVstrain') && ~isempty(SET(no).StrainTagging.globalRVstrain)
              outdata{line,10} = mynanmean(SET(no).StrainTagging.globalRVstrain(1,peaktf,:));
              radcol=11;
              circcol=28;
              col = 45;
            else
              radcol = 10;
              circcol = 27;
              col = 44;
            end
             
            
            %get output values for bullseye
            %counterclockwise from 12 o'clock if you are looking at the plot
            bullseyeplot = bullseyeradiallax;
            %basal
            %radcol = 10;
            outdata{line,radcol} = mynanmean(bullseyeplot(1:4,4));
            outdata{line,radcol+1} = mynanmean(bullseyeplot(5:8,4));
            outdata{line,radcol+2} = mynanmean(bullseyeplot(9:12,4));
            outdata{line,radcol+3} = mynanmean(bullseyeplot(13:16,4));
            outdata{line,radcol+4} = mynanmean(bullseyeplot(17:20,4));
            outdata{line,radcol+5} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line,radcol+6} = mynanmean(bullseyeplot(1:4,3));
            outdata{line,radcol+7} = mynanmean(bullseyeplot(5:8,3));
            outdata{line,radcol+8} = mynanmean(bullseyeplot(9:12,3));
            outdata{line,radcol+9} = mynanmean(bullseyeplot(13:16,3));
            outdata{line,radcol+10} = mynanmean(bullseyeplot(17:20,3));
            outdata{line,radcol+11} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line,radcol+12} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line,radcol+13} = mynanmean(bullseyeplot(3:8,2));
            outdata{line,radcol+14} = mynanmean(bullseyeplot(9:14,2));
            outdata{line,radcol+15} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line,radcol+16} = mynanmean(bullseyeplot(:,1));
            
            bullseyeplot = bullseyecirclax;
            %basal
            %circcol = 27;
            outdata{line,circcol} = mynanmean(bullseyeplot(1:4,4));
            outdata{line,circcol+1} = mynanmean(bullseyeplot(5:8,4));
            outdata{line,circcol+2} = mynanmean(bullseyeplot(9:12,4));
            outdata{line,circcol+3} = mynanmean(bullseyeplot(13:16,4));
            outdata{line,circcol+4} = mynanmean(bullseyeplot(17:20,4));
            outdata{line,circcol+5} = mynanmean(bullseyeplot(21:24,4));
            %mid
            outdata{line,circcol+6} = mynanmean(bullseyeplot(1:4,3));
            outdata{line,circcol+7} = mynanmean(bullseyeplot(5:8,3));
            outdata{line,circcol+8} = mynanmean(bullseyeplot(9:12,3));
            outdata{line,circcol+9} = mynanmean(bullseyeplot(13:16,3));
            outdata{line,circcol+10} = mynanmean(bullseyeplot(17:20,3));
            outdata{line,circcol+11} = mynanmean(bullseyeplot(21:24,3));
            %apical
            outdata{line,circcol+12} = mynanmean([bullseyeplot(1:2,2) ; bullseyeplot(21:24,2)]);
            outdata{line,circcol+13} = mynanmean(bullseyeplot(3:8,2));
            outdata{line,circcol+14} = mynanmean(bullseyeplot(9:14,2));
            outdata{line,circcol+15} = mynanmean(bullseyeplot(15:20,2));
            %apex
            outdata{line,circcol+16} = mynanmean(bullseyeplot(:,1));
            
            
            %slice based export of strain and strain rate
            col = 44;
           % nbrslices = length(SET(no).StrainTagging.saslices);
            outdata{line-2,col} = 'Peak radial strain [%]';
            %for sliceloop = 1:nbrslices
              %outdata{line-1,col} = sprintf('Slice %d',sliceloop);
              outdata{line,col} = SET(no).StrainTagging.globalrad(SET(no).StrainTagging.peaktf,1);
              col = col+1;
           % end
            outdata{line-2,col} = 'Peak circ. strain [%]';
            %for sliceloop = 1:nbrslices
              %outdata{line-1,col} = sprintf('Slice %d',sliceloop);
              outdata{line,col} = SET(no).StrainTagging.globalcirc(SET(no).StrainTagging.peaktf,1);
              col = col+1;
            %end
            outdata{line-2,col} = 'Radial strain rate [%/s]';
              outdata{line-1,col} = 'Upslope';
              outdata{line,col} = SET(no).StrainTagging.upsloperad(1,1);
              outdata{line-1,col+1} = 'Upslope time';
              outdata{line,col+1} = SET(no).StrainTagging.upsloperad(1,2);
              outdata{line-1,col+2} = 'Downslope';
              outdata{line,col+2} = SET(no).StrainTagging.downsloperad(1,1);
              outdata{line-1,col+3} = 'Downslope time';
              outdata{line,col+3} = SET(no).StrainTagging.downsloperad(1,2);
              col = col+4;
            outdata{line-2,col} = 'Longit. strain rate [%/s]';
              outdata{line-1,col} = 'Upslope';
              outdata{line,col} = SET(no).StrainTagging.upslopecircum(1,1);
              outdata{line-1,col+1} = 'Upslope time';  
              outdata{line,col+1} = SET(no).StrainTagging.upslopecircum(1,2);            
              outdata{line-1,col+2} = 'Downslope';
              outdata{line,col+2} = SET(no).StrainTagging.downslopecircum(1,1);
              outdata{line-1,col+3} = 'Downslope time';
              outdata{line,col+3} = SET(no).StrainTagging.downslopecircum(1,2);
            col = col+4; 
            line=line+3;
            
        end
        %line = line+1;
      %end
    end
  end; %loop over image stack
  if doStrain
    %store file
    
    %Create thumbnails before storing.
    calcfunctions('calcdatasetpreview');
    
    % Set view settings
    DATA.ViewPanels = 1;
    DATA.ViewPanelsType = {'one'};
    DATA.ViewMatrix = [1 1];
    DATA.ThisFrameOnly = 0;
    DATA.CurrentPanel = 1;
    DATA.CurrentTheme = 'lv';
    DATA.CurrentTool = 'select';

    %Save the file.
    %segment('filesaveall_Callback');
    corrupted=segment('checkcorrupteddataforautomaticsave');
    if corrupted
      corruptedfiles=sprintf('%s, %s',corruptedfiles,filename);
      mywarning(dprintf('Image file %s seems to be corrupted from last save. Please load and manually re-analyse strain to ensure that the image is not corrupted before saving',filename));
    else
      filemenu('saveallas_helper',pathname,files2load(fileloop).name);
      disp(sprintf('Saving %s',[pathname filesep files2load(fileloop).name]));
    end    
    
    rob.keyPress(KeyEvent.VK_SPACE);
    pause(0.1);
    rob.keyRelease(KeyEvent.VK_SPACE);
  end
  h = mywaitbarupdate(h);  
  %close current file
  segment('filecloseall_Callback');
end %loop over files
mywaitbarclose(h);

%--- Output to a string
segment('cell2clipboard',outdata);
%Stop the silent mode.
DATA.Silent = currentsilent;
      
  


%   outdata{1,8} = 'Radial strain rate';
%   outdata{2,8} = 'Systole Max inclination [%/s]';
%   outdata{2,9} = 'Time Systole Max inclination [s]';
%   outdata{2,8} = 'Diastole Max inclination [%/s]';
%   outdata{2,9} = 'Time Diastole Max inclination [s]';

        
%       %Get the different outputs depending on imageviewplane and
%       %imagetype
%       
%       imageviewplane = SET(no).ImageViewPlane;
%       if strcmp(imageviewplane,'Short-axis')
%         outdata{end+1,1} =imageviewplane;
%         hline=size(outdata,1);
%         outdata{end+2,1} = 'Strain';
%         saslices=['b','m','a'];
%         
%         %In order to get info on which segment are used
%         for ahaloop=1:17
%           [ahastri{ahaloop},pos(ahaloop)] = reportbullseye('aha17nameandpos',ahaloop); %Get name and position of export
%         end
%         
%         hline=size(outdata,1);
%         
%         %Segmental radial strain peak
%         segmentradptf=SET(no).StrainTagging.segmentrad(peaktf,:);
%         outdata{hline+2,2} = 'Segmental Radial Strain';
%         outdata{hline+3,1} = nan;
%         
%         for slice=1:SET(no).ZSize;
%           outdata{hline+3,3+2*(slice-1)} = ['Slice', ' ', num2str(slice)];
%         end
%         
%         line=size(outdata,1);
%         for slice=1:SET(no).ZSize;
%           if ismember(saslices(slice),'b')
%             sectors = 1:6;
%           elseif ismember(saslices(slice),'m') %mid or basal slices
%             sectors = 7:12;
%           else
%             sectors = 13:16;
%           end
%           for i = 1:length(sectors)
%             outdata{line+i,2+2*(slice-1)} = ahastri{sectors(i)};
%             outdata{line+i,3+2*(slice-1)} = segmentradptf(i);
%           end
%         end
%         
%         
%         vline=2*(length(sectors));
%         %Segmental circumferential strain peak
%         segmentcircptf=SET(no).StrainTagging.segmentcirc(peaktf,:);
%         outdata{hline+2,vline+1} = 'Segmental Circumferential Strain';
%         outdata{hline+3,vline} = nan;
%         
%         for slice=1:SET(no).ZSize;
%           outdata{hline+3,vline+2+2*(slice-1)} = ['Slice', ' ', num2str(slice)];
%         end
%         
%         for slice=1:SET(no).ZSize;
%           if ismember(saslices(slice),'b')
%             sectors = 1:6;
%           elseif ismember(saslices(slice),'m') %mid or basal slices
%             sectors = 7:12;
%           else
%             sectors = 13:16;
%           end
%           for i = 1:length(sectors)
%             outdata{hline+3+i,vline+1+2*(slice-1)} = ahastri{sectors(i)};
%             outdata{hline+3+i,vline+2+2*(slice-1)} = segmentcircptf(i);
%           end
%         end
%         
%         %outdataSR is radial outdataSRcirc circumferential
%         %Strain rate
%         vline=vline+2*(length(sectors)-1);
%         outdata{hline+2,vline+2} = 'Radial Strain Rate';
%         outdata{hline+3,vline+2} = 'Systole';
%         outdata{end+1,1} = nan;
%         for slice =1:SET(no).ZSize
%           outdata{hline+3,vline+3+2*(slice-1)} = ['Slice', ' ', num2str(slice)];
%         end
%         outdata{hline+4,vline+2} = 'Max inclination [%/s]';
%         for slice =1:SET(no).ZSize
%           outdata{hline+4,vline+3+2*(slice-1)} = SET(no).StrainTagging.upsloperad(slice,1);
%         end
%         outdata{hline+5,vline+2} = 'Time [ms]';
%         for slice =1:SET(no).ZSize
%           outdata{hline+5,vline+3+2*(slice-1)} = SET(no).StrainTagging.upsloperad(slice,2);
%         end
%         outdata{hline+7,vline+2} = 'Diastole';        
%         outdata{hline+8,vline+2} = 'Max inclination [%/s]';
%         for slice =1:SET(no).ZSize
%           outdata{hline+8,vline+3+2*(slice-1)} = SET(no).StrainTagging.downsloperad(slice,1);
%         end
%         
%         outdata{hline+9,vline+2} = 'Time [ms]';
%         
%         for slice =1:SET(no).ZSize
%           outdata{hline+9,vline+3+2*(slice-1)} = SET(no).StrainTagging.downsloperad(slice,2);
%         end
%         
%         vline=vline+SET(no).ZSize*2+1;
%         
%         %Circumferential strainrate
%         outdata{hline+2,vline+2} = 'Circumferential Strain Rate';
%         outdata{hline+3,vline+2} = 'Systole';
%         for slice =1:SET(no).ZSize
%           outdata{hline+3,vline+3+2*(slice-1)} = ['Slice', ' ', num2str(slice)];
%         end
%         outdata{hline+4,vline+2} = 'Max inclination [%/s]';
%         for slice =1:SET(no).ZSize
%           outdata{hline+4,vline+3+2*(slice-1)} = SET(no).StrainTagging.downslopecircum(slice,1);
%         end
%         outdata{hline+5,vline+2} = 'Time [ms]';
%         for slice =1:SET(no).ZSize
%           outdata{hline+5,vline+3+2*(slice-1)} = SET(no).StrainTagging.downslopecircum(slice,2);
%         end
%         
%         outdata{hline+7,vline+2} = 'Diastole';        
%         outdata{hline+8,vline+2} = 'Max inclination [%/s]';
%         
%         for slice =1:SET(no).ZSize
%           outdata{hline+8,vline+3+2*(slice-1)} = SET(no).StrainTagging.upslopecircum(slice,1);
%         end
%         
%         outdata{hline+9,vline+2} = 'Time [ms]';
%         
%         for slice =1:SET(no).ZSize
%           outdata{hline+9,vline+3+2*(slice-1)} = SET(no).StrainTagging.upslopecircum(slice,2);
%         end
%         
%         vline=vline+2*SET(no).ZSize+2;
%         
% %         if strcmp(SET(no).ImageType,'Strain from tagging')
% %           %torsion and rotation max values
% %           [mRot,mRot_i] =  max(SET(no).StrainTagging.slice_rotation{end}-SET(no).StrainTagging.slice_rotation{1});
% %           outdata{hline+2,vline+1} = 'Rotational Difference';
% %           outdata{hline+3,vline+1} = 'Peak [R_a-R_b]';
% %           outdata{hline+3,vline+2} = mRot;
% %           %                       outdata{hline+4,vline+1} = 'Frame';
% %           %                       outdata{hline+4,vline+2} = mRot_i;
% %           outdata{hline+4,vline+1} = 'Time [ms]';
% %           outdata{hline+4,vline+2} = SET(no).TimeVector(mRot_i)*1000;
% %           
% %           vline=vline+3;
% %           [mTor,mTor_i] = max(SET(no).StrainTagging.globaltorsion);
% %           outdata{hline+2,vline+1} = 'Mean Torsion';
% %           outdata{hline+3,vline+1} = 'd = distance between slices in long axis direction, r = radius';
% %           outdata{hline+4,vline+1} = 'Peak [(R_a-R_b)(r_a+r_b)/(2d)]';
% %           outdata{hline+4,vline+2} = mTor;
% %           %                       outdata{hline+5,vline+1} = 'Frame';
% %           %                       outdata{hline+5,vline+2} = mTor_i;
% %           outdata{hline+5,vline+1} = 'Time [ms]';
% %           outdata{hline+5,vline+2} = SET(no).TimeVector(mTor_i)*1000;
% %         end
%         
%         
%         
%         
%         
%         
%       else
%         outdata{end+1,1} = 'Long axis';
%         outdata{end+2,1} = 'Strain';
%         counter=1;
%         
%         ordered_taggroup = zeros(1,3);
%         ch_i_rad=[];
%         ch_i_circ=[];
%         for tagno= SET(no).StrainTagging.taggroup
%           peaktf=SET(tagno).StrainTagging.peaktf;
%           ch_i_rad(counter) = SET(tagno).StrainTagging.globalrad(peaktf);%(1:floor(length(SET(tagno).StrainTagging.globalrad)*0.7));
%           ch_i_circ(counter) = SET(tagno).StrainTagging.globalcirc(peaktf);%(1:floor(length(SET(tagno).StrainTagging.globalcirc)*0.7));
%           
%           %taggroup ordered as 2ch 3ch 4ch
%           switch SET(tagno).ImageViewPlane
%             case '2CH'
%               ordered_taggroup(1) = tagno;
%             case '3CH'
%               ordered_taggroup(2) = tagno;
%             case '4CH'
%               ordered_taggroup(3) = tagno;
%           end
%           counter=counter+1;
%         end
%         
%         outdata{end+1,2} ='Peak mean radial strain [%]';
%         chmeanrad =mynanmean(ch_i_rad);
%         outdata{end,3} = chmeanrad;
%         %outdata{end+1,2} = 'Time [ms]';
%         %outdata{end,3} = SET(tagno).TimeVector(peaktfrad)*1000;
%         
%         %Global circumferentiell strain peak
%         outdata{end+1,2} = 'Peak mean circumferential strain [%]';
%         %[chmeancirc, peaktfcirc] = min(mynanmean(ch_i_circ,2));
%         chmeancirc = mynanmean(ch_i_circ);
%         outdata{end,3} = chmeancirc;
%         %outdata{end+1,2} = 'Time [ms]';
%         %outdata{end,3} = SET(no).TimeVector(peaktfcirc)*1000;
%         hline=7;
%         vline=9;
%         
%         %Strain rate
%         outdata{end+2,1} = 'Strain Rate';
%         outdata{end+1,2} = 'Radial';
%         outdata{end+1,3} = 'Systole';
%         
%         %                   for tagno =  SET(no).StrainTagging.taggroup
%         %                         tmp_ind=find(ordered_taggroup==tagno)-1;
%         %                       outdata{end,4+tmp_ind*2} = SET(tagno).ImageViewPlane;
%         %                   end
%         outdata{end,4} = '2CH';
%         outdata{end,6} = '3CH';
%         outdata{end,8} = '4CH';
%         
%         outdata{end+1,3} = 'Max inclination [%/s]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{end,4+tmp_ind*2} = SET(tagno).StrainTagging.upsloperad(1);
%         end
%         
%         outdata{end+1,3} = 'Time [ms]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{end,4+tmp_ind*2} = SET(tagno).StrainTagging.upsloperad(2);
%         end
%         
%         outdata{end+2,3} = 'Diastole';
%         outdata{end+1,3} = 'Max inclination [%/s]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{end,4+tmp_ind*2} = SET(tagno).StrainTagging.downsloperad(1);
%         end
%         
%         outdata{end+1,3} = 'Time [ms]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{end,4+2*tmp_ind} = SET(tagno).StrainTagging.downsloperad(2);
%           counter = counter+1;
%         end
%         
%         %horisontal and vertical line
%         %hline=11;
%         %hline=8;
%         hline=11;
%         vline=9;
%         
%         outdata{hline+2,vline+4} = '2CH';
%         outdata{hline+2,vline+6} = '3CH';
%         outdata{hline+2,vline+8} = '4CH';
%         
%         %Circumferential strainrate
%         outdata{hline+1,vline+2} = 'Circumferential';
%         outdata{hline+2,vline+3} = 'Systole';
%         outdata{hline+3,vline+3} = 'Max inclination [%/s]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{hline+3,vline+4+2*tmp_ind} = SET(tagno).StrainTagging.downslopecircum(1);
%         end
%         
%         outdata{hline+4,vline+3} = 'Time [ms]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{hline+4,vline+4+2*tmp_ind} = SET(tagno).StrainTagging.downslopecircum(2);
%         end
%         hline=16;
%         %vline=vline+9;
%         
%         outdata{hline+1,vline+3} = 'Diastole';
%         outdata{hline+2,vline+3} = 'Max inclination [%/s]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{hline+2,vline+4+2*tmp_ind} = SET(tagno).StrainTagging.upslopecircum(1);
%         end
%         
%         outdata{hline+3,vline+3} = 'Time [ms]';
%         
%         for tagno =  SET(no).StrainTagging.taggroup
%           tmp_ind=find(ordered_taggroup==tagno)-1;
%           outdata{hline+3,vline+4+2*tmp_ind} = SET(tagno).StrainTagging.upslopecircum(2);
%         end
%         
%         checked=[checked,SET(no).StrainTagging.taggroup];
%       end
%       outdata{end+1,1}=nan;
%     end
%   end
% %   catch me
% %     %--- Some thing went wrong
% %     mydispexception(me);
% %     outdata{fileloop+1,2} = 'FAILED.';
% % end
% % h = mywaitbarupdate(h);
% % end; %loop over files
% % mywaitbarclose(h);
% % 
% % %--- Output to a string
% % segment('cell2clipboard',outdata);
% % %Stop the silent mode.
% % DATA.Silent = currentsilen;
% % SET=oldSet;
% % NO=oldNo;
% %Make sure starting with something fresh.
% %segment('filecloseall_Callback',true);



%---------------------------------
function exportthisstack(doheader) %#ok<DEFNU>
%---------------------------------
%Export data from current image stack to clipboard.
global NO

includenormalized = yesno('Do you want to include BSA normalized values?');

exportdata(doheader,includenormalized,NO);

%---------------------------------------------------------------
function varargout = exportdata(doheader,includenormalized,no)
%---------------------------------------------------------------
%This is the workhorse of export functions. 
%- doheader tells whether to include a header.
%- includenormalized tells whether to include BSA normalized data. 
%- no tells whether to export ONLY image stack no. 

global SET NO

if nargin<3
  [~,~,flowno,~,marno] = findfunctions('findno');
  no=findfunctions('findcineshortaxisno');
  scarno=findfunctions('findscarshortaxisno');
  if length(marno) > 1, marno = marno(1); end
  onlyone = false;
else
  onlyone = true;  
  if isempty(SET(no).Scar)
    scarno = [];
  else
    scarno = no;
  end;
  if isempty(SET(no).Flow)
    flowno = [];
  else
    flowno = no;
  end;
  if isempty(SET(no).MaR)
    marno = [];
  else
    marno = no;
  end;
end;

if isnan(no)
  no = NO;
end;

if isempty(no)
  no = NO;
end;

if nargin<2
  includenormalized = yesno('Do you want to include BSA normalized values?');
end;

outdata = cell(1,20);

if nargin==0
  doheader = true;
end;

if doheader
  %Write Header
  outdata = header(onlyone);
  row = 2;
else
  row = 1;
end;

if nargout>0
  varargout = cell(1,nargout);
end;

%Write data
if ~onlyone
  outdata{row, 1} = SET(no).FileName;
else
  outdata{row, 1} = sprintf('%d-%s %s',no,SET(no).ImageType,SET(no).ImageViewPlane);
end;

outdata{row, 2} = SET(no).PatientInfo.Name;
outdata{row, 3} = SET(no).PatientInfo.ID;
outdata{row, 4} = SET(no).PatientInfo.AcquisitionDate;
outdata{row, 5} = SET(no).PatientInfo.Age;
outdata{row, 6} = SET(no).PatientInfo.Length;
outdata{row, 7} = SET(no).PatientInfo.Weight;
outdata{row, 8} = SET(no).PatientInfo.Sex;
outdata{row, 9} = SET(no).PatientInfo.BSA;
outdata{row,10} = SET(no).HeartRate;
outdata{row,11} = 1000*SET(no).TSize*SET(no).TIncr;

lvmed = SET(no).LVM(SET(no).EDT);
lvmes = SET(no).LVM(SET(no).EST);
lvm = 0.5*(lvmed+lvmes);
try
  rvm = 0.5*(SET(no).RVM(SET(no).EDT)+SET(no).RVM(SET(no).EST));
catch %#ok<CTCH>
  rvm = NaN;
end;

if numel(SET(no).PatientInfo.BSA) ~= 1 || SET(no).PatientInfo.BSA==0
  bsa_1 = NaN;
else
  bsa_1 = 1/SET(no).PatientInfo.BSA;
end;

outdata{row,12} = lvm;
outdata{row,13} = lvm*1.05;
outdata{row,14} = bsa_1*lvm*1.05;
outdata{row,15} = SET(no).EDV;
outdata{row,16} = bsa_1*SET(no).EDV;
outdata{row,17} = SET(no).ESV;
outdata{row,18} = bsa_1*SET(no).ESV;
if isequal(SET(no).SV,0)
  outdata{row,19} = SET(no).EDV-SET(no).ESV;
  outdata{row,20} = bsa_1*(SET(no).EDV-SET(no).ESV);
else
  outdata{row,19} = SET(no).SV;
  outdata{row,20} = bsa_1*SET(no).SV;
end;

if isequal(SET(no).EF,0) && (SET(no).EDV>0)
  outdata{row,21} = 100*(SET(no).EDV-SET(no).ESV)/SET(no).EDV;
else
  outdata{row,21} = 100*SET(no).EF;
end;
outdata{row,22}= SET(no).HeartRate*SET(no).SV/1000;
outdata{row,23}= bsa_1*SET(no).HeartRate*SET(no).SV/1000;

%PEF PFR
outdata{row,24} = SET(no).PFR;
outdata{row,25} = SET(no).PER;

%rmv...
outdata{row,26} = rvm;
outdata{row,27} = rvm*1.05;
outdata{row,28} = bsa_1*rvm*1.05;
outdata{row,29} = SET(no).RVEDV;
outdata{row,30} = bsa_1*SET(no).RVEDV;
outdata{row,31} = SET(no).RVESV;
outdata{row,32} = bsa_1*SET(no).RVESV;
outdata{row,33} = SET(no).RVSV;
outdata{row,34} = bsa_1*SET(no).RVSV;
outdata{row,35} = 100*SET(no).RVEF;

if ~isempty(scarno)
  lvmde = SET(scarno).LVM(SET(scarno).EDT);
  outdata{row,36} = lvmde;
  if ~isempty(SET(scarno).Scar) %Moved this check here to get LVM out.
    outdata{row,37} = SET(scarno).Scar.Percentage;
    outdata{row,38} = lvmde*SET(scarno).Scar.Percentage/100;

    %Total extent
    tempstart = SET(scarno).StartSlice;
    tempend = SET(scarno).EndSlice;
    SET(scarno).StartSlice = 1;
    SET(scarno).EndSlice = SET(scarno).ZSize;
    [~,maxtransmurality,meantrans4infarct,totextent] = ...
      viability('calctransmuralityline',24,scarno);
    SET(scarno).StartSlice = tempstart;
    SET(scarno).EndSlice = tempend;
    outdata{row,39} = totextent;
    outdata{row,40} = meantrans4infarct;
    outdata{row,41} = max(maxtransmurality(:));
    
    if SET(scarno).Scar.UseWeighting
      %Total extent for weighted
      tempstart = SET(scarno).StartSlice;
      tempend = SET(scarno).EndSlice;
      SET(scarno).StartSlice = 1;
      SET(scarno).EndSlice = SET(scarno).ZSize;
      [~,~,maxtransmurality,meantrans4infarct,totextent] = ...
        viability('calctransmuralityweighted',24,scarno);
      SET(scarno).StartSlice = tempstart;
      SET(scarno).EndSlice = tempend;
      outdata{row,42} = totextent;
      outdata{row,43} = meantrans4infarct;
      outdata{row,44} = max(maxtransmurality(:));
    else
      outdata{row,42} = NaN;
      outdata{row,43} = NaN;
      outdata{row,44} = NaN;
    end

    %MO
    volscale = SET(scarno).ResolutionX*SET(scarno).ResolutionY*(SET(scarno).SliceThickness+SET(scarno).SliceGap)/1e3;
    outdata{row,45} = sum(SET(scarno).Scar.NoReflow(:))*volscale; %Region that is marked as where MO might reside.
    
    %For backwards compability.
    if isnan(SET(scarno).Scar.MOPercentage)
      viability('viabilitycalcvolume',scarno);
    end;
    
    outdata{row,46} = SET(scarno).Scar.MOPercentage;
  end;
else
  outdata{row,36} = NaN;
  outdata{row,37} = NaN;
  outdata{row,38} = NaN;
  outdata{row,39} = NaN;
  outdata{row,40} = NaN;
  outdata{row,41} = NaN;
  outdata{row,42} = NaN;
  outdata{row,43} = NaN;
  outdata{row,44} = NaN;
  outdata{row,45} = NaN;
  outdata{row,46} = NaN;
end;

if ~isempty(marno) && ~isempty(SET(marno).MaR)
  lvmmaredt = SET(marno).LVM(SET(marno).EDT);
  lvmmarest = SET(marno).LVM(SET(marno).EST);
  marpctedt = SET(marno).MaR.Percentage(SET(marno).EDT);
  if SET(marno).TSize > 1
    marpctest = SET(marno).MaR.Percentage(SET(marno).EST);
  else
    marpctest = NaN;
  end
  outdata{row,47} = marpctedt;
  outdata{row,48} = marpctest;
  outdata{row,49} = lvmmaredt*marpctedt/100;
  outdata{row,50} = lvmmarest*marpctest/100;
  
  if strcmp(SET(marno).ImagingTechnique,'NM')  %Spect data
    if ~isempty(SET(marno).MaR.MPS)
      outdata{row,51} = SET(marno).MaR.MPS.TPD;
      outdata{row,52} = SET(marno).MaR.MPS.TPDLAD;
      outdata{row,53} = SET(marno).MaR.MPS.TPDLCx;
      outdata{row,54} = SET(marno).MaR.MPS.TPDRCA;
    else
      outdata{row,51} = NaN;
      outdata{row,52} = NaN;
      outdata{row,53} = NaN;
      outdata{row,54} = NaN;
    end
  end;
else
  outdata{row,47} = NaN;
  outdata{row,48} = NaN;
  outdata{row,49} = NaN;
  outdata{row,50} = NaN;
end;

%Flow
if ~isempty(flowno)
  rowoffset = 0;
  for loop=1:length(flowno)
    NO = flowno(loop);

    if SET(flowno(loop)).RoiN>0
      %Calculate flow, call to get total flow.
      flowgo = reportflow;
      if flowgo
        flowgui = reportflow('getalldata');%SBT20160620
        tots = flowgui.nettotvol;%SBT20160620
        hr = flowgui.hr;
        co = flowgui.nettotvol.*flowgui.hr;%SBT20160620
        
        
        for rloop=1:length(tots)
          outdata{row+rowoffset,55} = SET(NO).Roi(rloop).Name;
          outdata{row+rowoffset,56} = tots(rloop);
          outdata{row+rowoffset,57} = co(rloop); %SBT20160620
          outdata{row+rowoffset,58} = hr;%SBT20160620
          rowoffset = rowoffset+1;
        end;
        reportflow('close_Callback');%SBT20160620
      end
    end;

  end;
end;

%Measurements
if ~onlyone
  loopover = 1:length(SET);
else
  loopover = no;
end;

msdone = 0;
for sloop=loopover
  for mloop=1:length(SET(sloop).Measure)
		outdata{row+msdone,59} = SET(sloop).Measure(mloop).Name;
    outdata{row+msdone,60} = SET(sloop).Measure(mloop).Length;
    msdone = msdone + 1;
  end;
end;

outdata{row,61} = SET(no).RVPFR;
outdata{row,62} = SET(no).RVPER;
outdata{row,63} = lvmed;
outdata{row,64} = lvmes;

%If not wanted remove non normalized values.
ind = true(1,size(outdata,2));
if ~includenormalized
  [~,indforbsa] = header(onlyone);
  ind(indforbsa) = false;
end;
outdata = outdata(:,ind);

%If called with no output arguments then copy to clipboard.
if nargout==0
  segment('cell2clipboard',outdata);
else
  varargout{1} = outdata;
end;

%------------------------------
function [x,y] = meshfixer(x,y)
%------------------------------
%Fix a mesh if the "layers" are twisted.
%The mesh is assumed to be rows of coordinates.

%Einar Heiberg

%Extract first row
xrow = x(1,:);
yrow = y(1,:);
for loop = 2:size(x,1)
  
  %initiate first solution
  mindist = 1e10;
  xmin = x(loop,:);
  ymin = y(loop,:);
  xtest = xmin;
  ytest = ymin;
  
  %Loop over all possible twistst to find best match
  for nloop = 1:size(x,2)
    dist = sum(sqrt((xtest-xrow).^2+(ytest-yrow).^2));
    if dist<mindist
      mindist = dist;
      xmin = xtest;
      ymin = ytest;
    end;
    
    %rotate for next iteration
    xtest = [xtest(2:end) xtest(1)];
    ytest = [ytest(2:end) ytest(1)];
  end;
  
  %Store the best match
  x(loop,:) = xmin;
  y(loop,:) = ymin;
  
  %Take the new row
  xrow = x(loop,:);
  yrow = y(loop,:);
end;
  
%---------------------------------------------------------------------------------
function export2stl_helper(no,x,y,resolution,tf,closeapex,fignr,pathname,filename)
%---------------------------------------------------------------------------------
%Helper function that fixes with coordinates etc before exporting

global SET 

%Extract timeframe
x = x(:,tf,:);
y = y(:,tf,:);
x = squeeze(x)'; %Now each row is indeed a row in the surface
y = squeeze(y)'; %Now each row is indeed a row in the surface

%Create slice data
z = 1:SET(no).ZSize;
z = repmat(z',[1 size(x,2)]);

%Cut vector
logind = ~isnan(x(:,1));
first = find(logind,1,'first');
last = find(logind,1,'last');
if ~isequal(sum(logind),last-first+1)
  myfailed('None consecutive segmentation detected (slices)');
  return;
end;
x = x(first:last,:);
y = y(first:last,:);
z = z(first:last,:);

%Fix the mesh nice
[x,y] = meshfixer(x,y);
    
%Check if only one slice

%remember the size
sz = size(x);

%Convert to MR coordinates
pos = calcfunctions('xyz2rlapfh',no,x(:),y(:),z(:));
rl = pos(:,1);
ap = pos(:,2);
fh = pos(:,3);

%Convert back to mesh size
rl = reshape(rl,sz);
ap = reshape(ap,sz);
fh = reshape(fh,sz);

%Run the export
if nargin==8
  fid = pathname;
  mesh2stl(rl,ap,fh,resolution,closeapex,fignr,fid); %Called with fid rather than pathname and filename  
else
  mesh2stl(rl,ap,fh,resolution,closeapex,fignr,pathname,filename);
end;

%---------------------------
function export2stl_Callback %#ok<DEFNU>
%---------------------------
%Exports mesh in current image stack as STL file. Asks for surface to
%export and then takes current timeframe to export. 

%Einar Heiberg

global SET NO

no = NO;

%Ask user to select surface
l = {'LV endocardium','LV epicardium','RV endocardium','RV epicardium','Abort'};
c = mymenu('Please select surface to export',l{:});

switch c
  case 1 %LV endo
    x = SET(no).EndoX;
    y = SET(no).EndoY;
  case 2 %LV epi
    x = SET(no).EpiX;
    y = SET(no).EpiY;
  case 3 %RV endo
    x = SET(no).RVEndoX;
    y = SET(no).RVEndoY;
  case 4 %RV epi
    x = SET(no).RVEpiX;
    y = SET(no).RVEpiY;
  otherwise
    myfailed('Aborted.');
    return;
end;

if isempty(x) || isempty(y)
  myfailed(dprintf('No data in %s',l{c}));
  return;
end;

%--- Calculate suggested resolution
resolution = (SET(no).ResolutionX+SET(no).ResolutionY+SET(no).SliceThickness+SET(no).SliceGap)/3;
s.Resolution_mm = resolution;
s.CloseApex = false;
[s,ok] = inputstruct(s,'Set resolution of STL file [mm]');
if ~ok
  mywarning('Invalid input, using default.');
end;
resolution = s.Resolution_mm;
closeapex = s.CloseApex;

[filename,pathname,~,ok] = myuiputfile('*.stl','Select filename for STL file export');
if ~ok
  myfailed('Aborted or illegal file selected.');
  return;
end;

figure(99);
clf;
set(99,'numbertitle','off','name','STL output display');

export2stl_helper(no,x,y,resolution,SET(no).CurrentTimeFrame,closeapex,99,pathname,filename);

mymsgbox('File exported.');

%-----------------------------------------------------------------------------------------
function exportclosesurfaces(no1,x1,y1,no2,x2,y2,closebase,closeapex,resolution,fignr,fid)
%-----------------------------------------------------------------------------------------
%Helper function to close to surfaces (basal) and export to STL file.
%Please not that this takes a mesh in form of points * slices (i.e
%timeframe already selected. This is different from elsewere in STL export.
%no1 and no2 are required to be able to calculate 3D coordinates (patient
%system)

%Einar Heiberg

if nargin<11
  error('Expected 11 input arguments.');
end;

if isempty(fignr)
  doplot = false;
else
  doplot = true;
  figure(fignr);
  hold on;
end;

%Extract basal rows from surface 1
basex1 = x1;
basey1 = y1;
baseslice1 = find(~isnan(basex1(1,:)),1,'first');
basex1 = basex1(:,baseslice1);
basey1 = basey1(:,baseslice1);

%Extract basal rows from surface 2
basex2 = x2;
basey2 = y2;
baseslice2 = find(~isnan(basex2(1,:)),1,'first');
basex2 = basex2(:,baseslice2);
basey2 = basey2(:,baseslice2);

%Extract apex rows from surface 1
apexx1 = x1;
apexy1 = y1;
apexslice1 = find(~isnan(apexx1(1,:)),1,'last');
apexx1 = apexx1(:,apexslice1);
apexy1 = apexy1(:,apexslice1);

%Extract basal rows from surface 2
apexx2 = x2;
apexy2 = y2;
apexslice2 = find(~isnan(apexx2(1,:)),1,'last');
apexx2 = apexx2(:,apexslice2);
apexy2 = apexy2(:,apexslice2);

for typeloop = 1:1 %base and apex
  
  switch typeloop
    case 1
      x1 = basex1;
      y1 = basey1;
      slice1 = baseslice1;
      x2 = basex2;
      y2 = basey2;      
      slice2 = baseslice2;
      closeit = closebase;
    case 2
      x1 = apexx1;
      y1 = apexy1;      
      slice1 = apexslice1;      
      x2 = apexx2;
      y2 = apexy2;            
      slice2 = apexslice2;
      closeit = closeapex;
  end;
  
  if closeit
    %Convert to 3D coordinates for surface 1
    %remember the size
    
    %Convert to MR coordinates
    pos1 = calcfunctions('xyz2rlapfh',no1,x1(:),y1(:),repmat(slice1,1,length(x1)));
    pos2 = calcfunctions('xyz2rlapfh',no2,x2(:),y2(:),repmat(slice2,1,length(x2)));
    rl1 = pos1(:,1);
    ap1 = pos1(:,2);
    fh1 = pos1(:,3);
    rl2 = pos2(:,1);
    ap2 = pos2(:,2);
    fh2 = pos2(:,3);
    
    %Calculate lengths of rows
    d1 = sum(sqrt((rl1(2:end)-rl1(1:(end-1))).^2 + (ap1(2:end)-ap1(1:(end-1))).^2 + (fh1(2:end)-fh1(1:(end-1))).^2));
    d2 = sum(sqrt((rl2(2:end)-rl2(1:(end-1))).^2 + (ap2(2:end)-ap2(1:(end-1))).^2 + (fh2(2:end)-fh2(1:(end-1))).^2));
    
    %Make sure d1>d2
    if d1<d2
      %--- Need to switch
      
      %Backup
      orl1 = rl1;
      oap1 = ap1;
      ofh1 = fh1;
      od1 = d1;
      
      %One direction
      rl1 = rl2;
      ap1 = ap2;
      fh1 = fh2;
      d1 = d2;
      
      %Use backup
      rl2 = orl1;
      ap2 = oap1;
      fh2 = ofh1;
      d2 = od1;
      
    end;
    
    %Resample
    newn1 = round(d1/resolution);
    newn2 = round(d2/resolution);
    n1 = 1:length(rl1);
    n2 = 1:length(rl2);
    if newn1<1
      newn1 = 1;
    end;
    if newn2<1
      newn2 = 1;
    end;
    if newn1>length(rl1)
      newn1 = length(rl1);
    end;
    if newn2>length(rl2)
      newn2 = length(rl2);
    end;
    ni1 = linspace(1,length(rl1)-(length(rl1)-1)/newn1,newn1);
    ni2 = linspace(1,length(rl2)-(length(rl2)-1)/newn2,newn2);
    rl1 = interp1(n1,rl1,ni1);
    ap1 = interp1(n1,ap1,ni1);
    fh1 = interp1(n1,fh1,ni1);
    rl2 = interp1(n2,rl2,ni2);
    ap2 = interp1(n2,ap2,ni2);
    fh2 = interp1(n2,fh2,ni2);
    
    %Try to match the surfaces, keep surface 1 fix
    mindist = 1e10;
    rl2best = rl2;
    ap2best = ap2;
    fh2best = fh2;
    rltest = rl2;
    aptest = ap2;
    fhtest = fh2;
    
    %Possible later check rotation
    
    %We know that d1>d2 => length(rl1)>=length(rl2)
    ind = round(linspace(1,length(rl2),length(rl1)));
    for loop = 1:length(rl2) %Loop over the number of points in r2
      
      %Calculate distance between surfaces of different shifts
      dist = sum(sqrt( (rltest(ind)-rl1).^2 + (aptest(ind)-ap1).^2 + (fhtest(ind)-fh1).^2 ));
      
      if dist<mindist
        %Store it
        rl2best = rltest;
        ap2best = aptest;
        fh2best = fhtest;
        mindist = dist;
      end;
      
      %Shift it
      rltest = [rltest(2:end) rltest(1)];
      aptest = [aptest(2:end) aptest(1)];
      fhtest = [fhtest(2:end) fhtest(1)];
      
    end;
    
    %Loop over r1 to create triangles (r1 got more or equal points as r2)
    nextind = [1:length(rl1) 1];
    ind = round(linspace(1,length(rl2),length(rl1)))';
    for loop = 1:length(rl1)
      
      %Extract point 1
      x1 = rl1(nextind(loop));
      y1 = ap1(nextind(loop));
      z1 = fh1(nextind(loop));
      
      %Extract point 2 (one row up)
      x2 = rl2best(ind(loop));
      y2 = ap2best(ind(loop));
      z2 = fh2best(ind(loop));
      
      %Extract point 3 (to the right)
      x3 = rl1(nextind(loop+1));
      y3 = ap1(nextind(loop+1));
      z3 = fh1(nextind(loop+1));
      
      %Write triangle 1->2->3
      normal = -cross([x2-x1,y2-y1,z2-z1],[x3-x1,y3-y1,z3-z1]);
      fprintf(fid, 'facet normal %f %f %f \n', normal );
      fprintf(fid, 'outer loop \n');
      fprintf(fid, 'vertex %f %f %f \n', x1,y1,z1);
      fprintf(fid, 'vertex %f %f %f \n', x2,y2,z2);
      fprintf(fid, 'vertex %f %f %f \n', x3,y3,z3);
      fprintf(fid, 'endloop \n');
      fprintf(fid, 'endfacet \n');
      
      if doplot
        %Graphical update
        patch([x1;x2;x3],[y1;y2;y3],[z1;z2;z3],zeros(3,1));
      end;
        
    end;
    
  end; %Closeit
end; %Loop over base,apex (1:2)

if doplot
  hold off;
end;

%-----------------------------------------------
function exportlv2stl(no,tf,fignr,filename,resolution)
%-----------------------------------------------
%Exports LV surface to STL (with inner and outer contour, closed in apex
%and surfaces merged in the base.

%Einar Heiberg

global SET

fid = fopen(filename,'w');
if isequal(fid,-1)
  myfailed(dprintf('Could not open the file %s for writing.',[pathname filename]));
  return;
end;
  
%Extract endocardium
endox = SET(no).EndoX;
endoy = SET(no).EndoY;

%Extract epicardium
epix = SET(no).EpiX;
epiy = SET(no).EpiY;

%Start the file
fprintf(fid,'solid segment stl \n');

%Endocardial export. Flip it to get other direction
closeapex = true;
export2stl_helper(no,flip(endox,1),flip(endoy,1),resolution,tf,closeapex,fignr,fid);

%Ensure that the epicardial surface is closed
export2stl_helper(no,epix,epiy,resolution,tf,closeapex,fignr,fid);

%Close the two surfaces together.
closebase = true;
closeapex = false;
exportclosesurfaces(...
  no,...
  squeeze(endox(:,tf,:)),...
  squeeze(endoy(:,tf,:)),...
  no,...
  squeeze(epix(:,tf,:)),...
  squeeze(epiy(:,tf,:)),...
  closebase,...
  closeapex,...
  resolution,...
  fignr,...
  fid);

%Close the file
fprintf(fid,'endsolid stl\n');
fclose(fid);

%-----------------------------
function exportlv2stl_Callback %#ok<DEFNU>
%-----------------------------
%Exports LV (endo and epi) as one surface and takes care of closing the
%surface suitable to read into CAD software.

%Einar Heiberg

global SET NO DATA

no = NO;

%Check if data exists
if isempty(SET(no).EndoX) || isempty(SET(no).EpiX)
  myfailed('Either endocardium or epicardium is missing.');
  return;
end;

%Fix resolution
resolution = 5; %Later ask

%Open a file
pathname = DATA.Pref.exportpath; %Later asks  
pathname = myuigetdir(pathname,'Select a folder with .mat files');
if isequal(pathname,0)
  myfailed('Aborted.');
  return;
end;
filename = 'lv_export.stl';

%Call the function to do it.
tf = SET(no).CurrentTimeFrame;
fignr = 98;
exportlv2stl(no,tf,fignr,[pathname filesep filename],resolution);

%Message
mymsgbox('LV exported.');

%-----------------------------------------------
function [xnew,ynew] = innersurface(x,y,no,dist)
%-----------------------------------------------
%Takes a surface x,y and returns a new surface inside of this surface with
%a distance of dist mm. x are assumed to be size NumDataPoints*Z

%Einar Heiberg

global SET

xnew = x;
ynew = y;

%Loop over slices
for zloop = 1:size(x,3)
  for tloop = 1:size(x,2)
  
    if ~isnan(x(1,tloop,zloop))
    
      xp = conv([x(:,tloop,zloop);x(1,tloop,zloop)],[-1;1],'valid');
      yp = conv([y(:,tloop,zloop);y(1,tloop,zloop)],[-1;1],'valid');
      
      %Normalize
      n = sqrt(xp.^2+yp.^2);
      xp = xp./n;
      yp = yp./n;
  
      %Move points
      xnew(:,tloop,zloop) = x(:,tloop,zloop) - yp*dist/SET(no).ResolutionX;
      ynew(:,tloop,zloop) = y(:,tloop,zloop) + xp*dist/SET(no).ResolutionY;
      
      %Make sure last point and first point are the same
      xnew(end,tloop,zloop) = xnew(1,tloop,zloop);
      ynew(end,tloop,zloop) = ynew(1,tloop,zloop);
      
      %figure(98);
      %plot(y(:,tloop,zloop),x(:,tloop,zloop));
      %hold on;
      %plot(ynew(:,tloop,zloop),xnew(:,tloop,zloop));
      %hold off;
      %pause
      
    end;
    
  end;
  
end;

%-----------------------------------------------------
function exportrv2stl(no,tf,fignr,filename,resolution)
%-----------------------------------------------------
%Exports RV surface to STL (with inner and outer contour, closed in apex
%and surfaces merged in the base. OBS this function ignores RV epicardium
%and fakes in an epicardium. Wall thickness is set to 1mm.

%Einar Heiberg

global SET

fid = fopen(filename,'w');
if isequal(fid,-1)
  myfailed(dprintf('Could not open the file %s for writing.',[pathname filename]));
  return;
end;
  
%Extract endocardium
endox = SET(no).RVEndoX;
endoy = SET(no).RVEndoY;

%Extract surface
endox = endox(:,tf,:);
endoy = endoy(:,tf,:);

%Calc epicardium
dist = -1; %mm
[epix,epiy] = innersurface(endox,endoy,no,dist);

%Start the file
fprintf(fid,'solid segment stl \n');

%Endocardial export. Flip it to get other direction
closeapex = true;
tf = 1; %as we cropped above
export2stl_helper(no,flip(endox,1),flip(endoy,1),resolution,tf,closeapex,fignr,fid);

%Ensure that the epicardial surface is closed
export2stl_helper(no,epix,epiy,resolution,tf,closeapex,fignr,fid);

%Close the two surfaces together.
closebase = true;
closeapex = false;
exportclosesurfaces(...
  no,...
  squeeze(endox(:,tf,:)),...
  squeeze(endoy(:,tf,:)),...
  no,...
  squeeze(epix(:,tf,:)),...
  squeeze(epiy(:,tf,:)),...
  closebase,...
  closeapex,...
  resolution,...
  fignr,...
  fid);

%Close the file
fprintf(fid,'endsolid stl\n');
fclose(fid);

%-----------------------------
function exportrv2stl_Callback %#ok<DEFNU>
%-----------------------------
%Exports RV to STL file. This file ignores epicardium and fakes in a
%epicardium.

%Einar Heiberg

global SET NO DATA

no = NO;

%Check if data exists
if isempty(SET(no).RVEndoX) 
  myfailed('RV endocardium is missing.');
  return;
end;

%Fix resolution
resolution = 5; %Later ask

%Open a file
pathname = DATA.Pref.exportpath; 
pathname = myuigetdir(pathname,'Select a folder with .mat files');
if isequal(pathname,0)
  myfailed('Aborted.');
  return;
end;
filename = 'rv_export.stl';

%Call the function to do it.
tf = SET(no).CurrentTimeFrame;
fignr = 98;
exportrv2stl(no,tf,fignr,[pathname filesep filename],resolution);

%Message
mymsgbox('RV exported.');

%-------------------------------------------------------------------
function exportall2stl(no,pathname,filetemplate,fignr,resolution,tf)
%-------------------------------------------------------------------
%Helper function to exportall2stl

%Einar Heiberg

global SET DATA 

%Export LV
disp('Exporting LV');
if (~isempty(SET(no).EndoX)) && (~isempty(SET(no).EpiX))
  exportlv2stl(no,tf,fignr,[pathname filesep filetemplate '-lv.stl'],resolution)
end;

%Export RV
disp('Exporting RV');
if (~isempty(SET(no).RVEndoX))
  exportrv2stl(no,tf,fignr,[pathname filesep filetemplate '-rv.stl'],resolution)
end;

%Export LV atria
disp('Exporting left atria (RV Epi)');
if (~isempty(SET(no).RVEpiX))
  exportleftatria2stl(no,tf,fignr,[pathname filesep filetemplate '-leftatria.stl'],resolution); 
end;

%**** Export ROI's ****
 
%--- Check ROI's for unique names
names = cell(1,SET(no).RoiN);
for loop = 1:length(names)
  names{loop} = SET(no).Roi(loop).Name;
end;
names = union(names,{});
  
%--- Loop over roi uniquenames and export
for namesloop = 1:length(names)
  templatename = names{namesloop};
  
  x = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
  y = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
    
  %Loop over ROI's
  for loop = 1:SET(no).RoiN
    if isequal(SET(no).Roi(loop).Name,templatename);
      if isnan(x(1,1,SET(no).Roi(loop).Z))
        x(:,:,SET(no).Roi(loop).Z) = SET(no).Roi(loop).X;
        y(:,:,SET(no).Roi(loop).Z) = SET(no).Roi(loop).Y;
      else
        mywarning('Two ROIs in the same slice. Ambigous result.',DATA.GUI.Segment);
      end;
    end;
  end;
    
  %Loop over ED and ES
  for tloop = 1:2
    
    if (tloop==1)
      %Export
      roifilename = sprintf('%s_%s_%s.stl',filetemplate,'ED',templatename);
    else
      roifilename = sprintf('%s_%s_%s.stl',filetemplate,'ES',templatename);      
    end;
    

    %Open file
    fid = fopen([pathname filesep roifilename],'w');
    if isequal(fid,-1)
      myfailed(dprintf('Could not open the file %s for writing.',[pathname filesep roifilename]));
      return;
    end;

    %Start the file
    fprintf(fid,'solid segment stl \n');
   
    %"endocardium",
    if (tloop==1)
      tf = SET(no).EDT;
    else
      tf = SET(no).EST;      
    end;
    
    endox = x(:,tf,:);
    endoy = y(:,tf,:);

  disp(dprintf('Exporting ROI %s',templatename)); %#ok<DSPS>
  
    closeapex = false;
    %flip it to get other direction (normal inwards)
    export2stl_helper(no,flip(endox,1),flip(endoy,1),resolution,1,closeapex,fignr,fid); %1 = timeframe as already cropped
     
    %"epicardium"
    [epix,epiy] = innersurface(endox,endoy,no,-1); %-1 is 1mm thickness
  
    %Export "epicardium"
    export2stl_helper(no,epix,epiy,resolution,1,closeapex,fignr,fid); %tf = 1 as we already have croppend

    %Close the two surfaces together.
    closebase = true;
    closeapex = true;
    exportclosesurfaces(...
      no,...
      endox,...
      endoy,...
      no,...
      endox,...
      endoy,...
      closebase,...
      closeapex,...
      resolution,...
      fignr,...
      fid);
  
    %Close the file
    fprintf(fid,'endsolid stl\n');
    fclose(fid);
  end; %tloop, i.e ED and ES
end; %names (roi's)

%------------------------------
function exportall2stl_Callback %#ok<DEFNU>
%------------------------------
%Exports all contours and ROI's to stl file(s)
%
%LV endo & epi => left ventricle
%RV endo => right ventricle with 1mm thickness
%RV epi => left atrium closed in base
%
%ROI's => tubes with 1 mm thickness

%Einar Heiberg

global DATA SET NO

no = NO;

%--- Calculate suggested resolution
resolution = (SET(no).ResolutionX+SET(no).ResolutionY+SET(no).SliceThickness+SET(no).SliceGap)/3;
s.Resolution_mm = resolution;
s.File_Name = SET(1).PatientInfo.Name;
[s,ok] = inputstruct(s,'Set resolution and name of STL file');
if ~ok
  mywarning('Invalid input, using default.');
end;
resolution = s.Resolution_mm;
filetemplate = s.File_Name;

%Ask for pathname
pathname = DATA.Pref.exportpath; 
pathname = myuigetdir(pathname,'Select a folder where to the output files.');
if isequal(pathname,0)
  myfailed('Aborted.');
  return;
end;

tf = SET(no).CurrentTimeFrame;

%Prepare graphical output
fignr = 99;
figure(fignr);
clf;
set(fignr,'numbertitle','off','name','STL export display');

exportall2stl(no,pathname,filetemplate,fignr,resolution,tf)

mymsgbox('Files exported.');

figure(fignr);
axis image off;
cameratoolbar(fignr);

%------------------------------------------------------------
function exportleftatria2stl(no,tf,fignr,filename,resolution)
%------------------------------------------------------------
%Exports RV epicardium surface to STL but names it LV atria. Closed in
%base. Wall thickness is set to 1mm.

%Einar Heiberg

global SET

fid = fopen(filename,'w');
if isequal(fid,-1)
  myfailed(dprintf('Could not open the file %s for writing.',[pathname filename]));
  return;
end;
  
%Extract endocardium
endox = SET(no).RVEpiX;
endoy = SET(no).RVEpiY;

%Extract surface
endox = endox(:,tf,:);
endoy = endoy(:,tf,:);

%Calc epicardium
dist = -1; %mm
[epix,epiy] = innersurface(endox,endoy,no,dist);

%Start the file
fprintf(fid,'solid segment stl \n');

%Endocardial export. Flip it to get other direction
closeapex = false;
tf = 1; %as we cropped above
export2stl_helper(no,flip(endox,1),flip(endoy,1),resolution,tf,closeapex,fignr,fid);

%Ensure that the epicardial surface is closed
export2stl_helper(no,epix,epiy,resolution,tf,closeapex,fignr,fid);

%Close the two surfaces together.
closebase = true;
closeapex = false;
exportclosesurfaces(...
  no,...
  squeeze(endox(:,tf,:)),...
  squeeze(endoy(:,tf,:)),...
  no,...
  squeeze(epix(:,tf,:)),...
  squeeze(epiy(:,tf,:)),...
  closebase,...
  closeapex,...
  resolution,...
  fignr,...
  fid);

%Close the file
fprintf(fid,'endsolid stl\n');
fclose(fid);

%--------------------------------
function exportbatch2stl_Callback %#ok<DEFNU>
%--------------------------------
%Export all2stl for multiple .mat files.

%Einar Heiberg

global DATA SET NO

pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder with .mat files');
if isequal(pathname,0)
  myfailed('Aborted.');
  return;
end;

files2load = dir([pathname filesep '*.mat']);
numfiles = length(files2load);

if numfiles == 0
  myfailed('No files for STL export.');
  return;
end;

%Ask for resolution
s.Resolution_mm = 3;
[s,ok] = inputstruct(s,'Set resolution and name of STL file');
if ~ok
  mywarning('Invalid input, using default.');
end;
resolution = s.Resolution_mm;

fignr = []; %no plotting.

%--- Loop over files
h = mywaitbarstart(numfiles,'Please wait, processing.');
for fileloop=1:numfiles

  %--- Load file
  DATA.Silent = true;
  disp(dprintf('Loading %s.',files2load(fileloop).name)); %#ok<DSPS>
  SET = []; %#ok<NASGU>
  load([pathname filesep files2load(fileloop).name],'-mat');
  SET = setstruct;
  clear setstruct;
  
  openfile('setupstacksfrommat',1);
  segment('renderstacksfrommat');
  
  %Find cine stack 
  no = findfunctions('findcineshortaxisno');
  
  NO = no; %for sure.
  
  %---Export
  filetemplate = SET(no).PatientInfo.Name;
  
  %ED
  tf = SET(no).EDT;
  exportall2stl(no,pathname,[filetemplate '_ED'],fignr,resolution,tf);   
  
  %ES
  tf = SET(no).EST;
  exportall2stl(no,pathname,[filetemplate '_ES'],fignr,resolution,tf);  
  
  h = mywaitbarupdate(h);
  
end;
mywaitbarclose(h);

DATA.Silent = false;