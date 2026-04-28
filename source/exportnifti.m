function [varargout] = exportnifti(varargin) %skiptranslation
%Function that batch export images and segmentation to nifti files.
%User is asked what to export and select a folder of .mat-files
%This function is slightly a misnormer in the filename as it also can
%export annotation points as a .csv-file

%Einar Heiberg

%skiptranslation

%#ok<*GVMIS>

if nargin==0
  exportinit;
else
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
end

%------------------
function exportinit
%------------------
%This is the main function that asks user and loops over files. Calls
%helper function exportthis to do the exporting.

global DATA SET

if ~yesno('Experimental code for research only, continue?')
  return
end

%Check if Segment 3DPrint
if isequal(DATA.ProgramName,'Segment 3DPrint')
  issegment3dp = true;
else
  issegment3dp = false;
end

%issegment3dp = false; %for debugging

%Check what to export
d = askwhattoexport(issegment3dp);

%Check if aborted
if isempty(d)
  return
end

%Store also issegment3dp to struct
d.IsSegment3DP = issegment3dp;

%Ask where to load files from
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder with .mat files');
if isequal(pathname,0)
  myfailed('Aborted.');
  return;
end

%Find files to process
files2load = dir([pathname filesep '*.mat']);
numfiles = length(files2load);

if numfiles==0
  myfailed('Found no files to summarize.');
  return;
end

%Loop over files to load

%Loop over all files
h = maingui.mywaitbarstart(numfiles,'Please wait, loading and summarizing files.',1);
for fileloop=1:numfiles

  %--- Load file
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.

  logdisp(sprintf('Loading %s.',files2load(fileloop).name));


  %---- try
  try
    SET = []; %Make sure a fresh start
    %Load
    load([pathname filesep files2load(fileloop).name],'-mat'); %#ok<LOAD>

    %Assign
    SET = setstruct;
    clear setstruct;

    %Call to intialize all variables correcly after loaded data.
    openfile('setupstacksfrommat',1);
    segment('renderstacksfrommat');

    %--- Process the loaded file

    %Get filename
    filename = files2load(fileloop).name;
    [~,name,~] = fileparts(filename);
    if ~isempty(pathname)
      filename = [pathname filesep name];
    else
      filename = name;
    end

    %Export for this file
    exportthis(filename,d); %d is a struct with what to do

  catch %#ok<CTCH>
    %--- Some thing went wrong
    logdisp(sprintf('Something went wrong with file %s.',files2load(fileloop).name));
  end

  h = maingui.mywaitbarupdate(h);
end %loop over files
maingui.mywaitbarclose(h);

%Make sure starting with something fresh.
segment('filecloseall_Callback',true);

%Stop the silent mode.
DATA.Silent = false;

mymsgbox('Export finished');

%------------------------------
function exportthis(filename,d)
%------------------------------
%Export this image stack with the data specified by the struct d that
%contains what to export.

global SET

%Find image stack
if d.IsSegment3DP
  no = 1; %For now only take from the first SET if Segment3DP.
else
  no = findfunctions('findcineshortaxisno'); %Get which is the shortaxis stack
end

if isempty(no)
  logdisp(['No cine SAX to export from file ',filename]);
  return
end

if d.doEDESTime
  %Check if ED&ES was selected; update text in the filename
  doEDES = true;
else
  doEDES = false;
end

%--- Export image volume
if d.doImages
  %Write image
  %donifiti_helper('Image',no,filename, doEDES)
  imagefilename = [filename '_IMAGE'];
  vol = squeeze(SET(no).IM); %Remove singleton dimensions
  if ~isempty(SET(no).IntensityScaling)
    if not(isa(vol,'int16'))
      vol = vol*SET(no).IntensityScaling+SET(no).IntensityOffset;
    else
      vol = single(vol);
    end
  end

  savenifitii_helper(vol,no, imagefilename, doEDES)
end

%--- Export LV endocardium
if d.doLVEndocardium
  vol = getbinaryvolumehelper('Endo',no);
  vol = uint8(255)*uint8(vol); %Convert to uint8
  imagefilename = [filename '_LVENDOCARDIUM'];

  savenifitii_helper(vol,no, imagefilename, doEDES)
end

%--- Export LV epicardium, abit of duplicated code here, but ok since short
%and simple
if d.doLVEpicardium
  vol = getbinaryvolumehelper('Epi',no);
  vol = uint8(255)*uint8(vol); %Convert to uint8
  imagefilename = [filename '_LVEPICARDIUM'];

  savenifitii_helper(vol,no, imagefilename, doEDES)
end

%--- Export RV endocardium
if d.doRVEndocardium
  %donifiti_helper('RVEndo',no, filename, doEDES)
  vol = getbinaryvolumehelper('RVEndo',no);
  vol = uint8(255)*uint8(vol); %Convert to uint8
  imagefilename = [filename '_RVENDOCARDIUM'];

  savenifitii_helper(vol,no, imagefilename, doEDES)
end

%--- Export RV epicardium
if d.doRVEpicardium
  %donifiti_helper('RVEpi',no, filename, doEDES)
  vol = getbinaryvolumehelper('RVEpi',no);
  vol = uint8(255)*uint8(vol); %Convert to uint8
  imagefilename = [filename '_RVEPICARDIUM'];

  savenifitii_helper(vol,no, imagefilename, doEDES)
end

%--- Export myocardium and blood pool (LV and RV)
if d.doMyocardium
  %0 is background
  %85 is LV blood pool
  %170 is RV blood pool
  %255 is LV myocardium

  %Get binary volumes
  lvendovol = getbinaryvolumehelper('Endo',no);
  lvepivol = getbinaryvolumehelper('Epi',no);
  rvendovol = getbinaryvolumehelper('RVEndo',no);

  %create empty volume, with uint8(0)
  vol = uint8(0)*uint8(lvendovol);

  %Assign
  vol(lvendovol) = uint8(85);
  vol(rvendovol) = uint8(170);
  vol(and(lvepivol,~lvendovol)) = uint8(255); %this is LV myocardium

  imagefilename = [filename '_MYOCARDIUM'];
  savenifitii_helper(vol,no, imagefilename, doEDES)
end

%--- Export Flow ROIs
if d.doFlowROIs
  [flowno,flowroi] = findfunctions('findflowaxisno'); %Returns only one

  if isempty(flowno) || isempty(flowroi)
    %No flow found
    logdisp(sprintf('No flow found in %s',filename));
  else

    %Mag image
    imagefilename = [filename '_FLOWIMAGE'];
    %savenifitii_helper(vol,flowno, imagefilename, doEDES)
    vol = (SET(flowno).IM); %squeeze(SET(flowno).IM);
    %savenifitii_helper(vol,flowno, imagefilename, doEDES)
    info = getinfostruct(vol,flowno);
    vol = volumefixer(vol,flowno);
    niftiwrite(vol,imagefilename,info);
    logdisp(sprintf('Wrote file %s',imagefilename));

    %Flow segmentation
    volmask = false(size(vol)); %Create emptymask
    for tloop = 1:SET(flowno).TSize
      for zloop = 1:SET(flowno).ZSize %it is extremly rare to have Z>1, but why not prepare for it
        mask = segment('createmask',...
          [SET(flowno).XSize SET(flowno).YSize],...
          SET(flowno).Roi(flowroi).Y(:,tloop,zloop),... %y first since matlab order of x & y
          SET(flowno).Roi(flowroi).X(:,tloop,zloop));
        volmask(:,:,zloop,tloop) = mask; %assign it to volume
      end
    end

    vol = uint8(255)*uint8(volmask); %Convert to uint8

    imagefilename = [filename '_FlOWMASK'];
    info = getinfostruct(vol,flowno);
    %vol = volumefixer(vol,flowno);
    niftiwrite(vol,imagefilename,info);
    logdisp(sprintf('Wrote file %s',imagefilename));
  end %end of there is flow

end

%--- Export levelset objects
if d.doObjects

  %Check if there are any objects
  hasobjects = false;
  if ~isempty(SET(no).LevelSet)
    if segment3dp.tools('getnumobjects')>0
      hasobjects = true;
    end
  end

  %If there are objects, lets export them
  if hasobjects
    O = SET(no).LevelSet.Object; %Extract the object for easier accesss.
    bwindx = find(cellfun(@(x)isequal(x,'BW'),{O.Objs.Type})); %index of BW objects
    for oloop = bwindx  %length(O.Objs)

      %Set filename
      imagefilename = [filename '_' removeforbiddenchars(O.Objs(oloop).Name)]; %new filename

      %Extract volume and convert to 0 for background and 255 for object
      vol = segment3dp.tools('getobjectbw',oloop,no) > uint8(127);

      vol = uint8(255)*uint8(vol);
      vol = volumefixer(vol,no); %convert time and z

      %Write image
      info = getinfostruct(vol,no);
      niftiwrite(vol,imagefilename,info);
      logdisp(sprintf('Wrote file %s',imagefilename));

    end %end of oloop

  else
    %There are no objects, report it.
    logdisp(sprintf('No objects in file %s',filename));
  end

end

%--- Export points (as .csv-file)
if d.doPoints

  %Check if points in this file, loop over all stacks
  containpoints = [];
  for loop = 1:length(SET)
    if ~isempty(SET(loop).Point.X)
      containpoints = [containpoints loop]; %#ok<*AGROW>
    end
  end

  if ~isempty(containpoints)

    %Open file
    fid = fopen([filename '.csv'],"w");
    if fid<0
      myfailed(sprinf('Could not open file %s for writing',[filename '.csv']));
      return
    end

    %Write header
    fprintf(fid,'Stack,Label,X,Y,Z,T\n');

    %Loop over stacks
    for sloop = 1:length(containpoints)
      sno = containpoints(sloop);

      %Loop over points in the stack
      for loop = 1:length(SET(sno).Point.X)
        p = SET(sno).Point;
        fprintf(fid,'%d,"%s",%f,%f,%d,%d\n',...
          sno,...
          removeinvalid(p.Label{loop}),... %internal that removes , and "
          p.Y(loop),...
          p.X(loop),...
          SET(sno).ZSize-p.Z(loop)+1,... %flip in Z to get consistent with exported nifti images that are flipped.
          p.T(loop));
      end
    end %Loop over all stacks with points

    %Close file
    fclose(fid);
    logdisp(sprintf('Wrote %s.csv',filename));
  else
    logdisp(sprintf('No points in filename %s.',filename));
  end

end

%----------------------------------
function stri = removeinvalid(stri)
%----------------------------------
%Removes , and " characters

logind = (stri==',') | (stri=='"');
stri = stri(~logind);

%------------------------------------
function info = getinfostruct(vol,no, edesonly)
%------------------------------------
%Computes struct for nifti export

global DATA SET

if nargin < 3
  edesonly = false;
end

%Define info struct
info = [];
info.Description = sprintf('Created by %s %s.',DATA.ProgramName, DATA.ProgramVersion);

%Compute imagesize, note this is in the order Nifti expects, need to call
%the helper function volumefix later to get the volumes in order with this.
imagesize = [SET(no).XSize SET(no).YSize];
if SET(no).ZSize >= 1
  imagesize = [imagesize SET(no).ZSize];
end
if SET(no).TSize >= 1 && ~edesonly
  imagesize = [imagesize SET(no).TSize];
else
  imagesize = [imagesize 1];
end

info.ImageSize =  imagesize; % SET(no).TSize];
info.Version = 'NIfTI1';

%Set pixel dimensions
info.PixelDimensions = [SET(no).ResolutionY SET(no).ResolutionX ];
if SET(no).ZSize >= 1
  info.PixelDimensions = [info.PixelDimensions (SET(no).SliceThickness+SET(no).SliceGap)]; %add slicethickness and gap is 3D
end

if SET(no).TSize > 1 && ~edesonly
  info.PixelDimensions = [info.PixelDimensions SET(no).TIncr];
  info.TimeUnits = 'Second';
else
  info.PixelDimensions = [info.PixelDimensions 1];
  info.TimeUnits = 'None';
end

info.Datatype = class(vol);
info.SpaceUnits = 'Millimeter';
info.SliceCode = 'Unknown';
info.FrequencyDimension = 0;
info.PhaseDimension = 0;
info.Qfactor = 1; %Scaling of transform
if isequal(info.Datatype,'single')
  info.IntentDescription = 'Image';
else
  info.IntentDescription = 'Segmentation';
end

%Currently we do not add image rotation information.
%This can later be done in Transform which is an affine3d struct

if SET(no).ZSize>1
  info.SpatialDimension = 3; %3D data
else
  info.SpatialDimension = 2; %2D data
end

%---------------------------------------------
function bvol = getbinaryvolumehelper(type,no)
%---------------------------------------------
%Returns binary volume with segmentation for the type. Type is either
%Endo,Epi,RVEndo,RVEpi

global SET

%Always return something
bvol = false(size(SET(no).IM));

%Check if there are any segmentation, if not return (which will return
%false volume).
if isempty(SET(no).([type 'X'])) %add X to check for the variable for instance EndoX when type is Endo
  logdisp(sprintf('%s is empty.',type));
  return
end


%---Ok now we know there is segmentation

%Loop over time and z. Although the other code is not prepared for time
%resolved images yet, we might as well do it here.

for tloop = 1:SET(no).TSize
  for zloop = 1:SET(no).ZSize

    %Create mask for one slice
    mask = segment('createmask',...
      [SET(no).XSize SET(no).YSize],...
      SET(no).([type 'Y'])(:,tloop,zloop),... %y first since matlab order of x & y
      SET(no).([type 'X'])(:,tloop,zloop));

    %Store mask
    bvol(:,:,tloop,zloop) = mask;
  end
end

%--------------------------------------
function volout = volumefixer(volin,no)
%--------------------------------------
%Reorder so that time is last and that z is flipped (to get right hand
%coordinate system which is dictacted by nifti format.

global SET

%Extract image
if SET(no).TSize > 1
  volout = permute(volin,[2 1 4 3]); %Segment stores as Y * X * T * Z whereas Nifti is X * Y * Z * T
else
  volin = permute(volin,[2 1 3]);
  volout = squeeze(volin);
end

%Extract image
if SET(no).ZSize >= 1
  volout = flip(volout,3);
end

%Ensure squeezed
%volout = squeeze(volout);

%-----------------------------------------
function d = askwhattoexport(issegment3dp)
%-----------------------------------------
%Select what to export

d = []; %Always return something

s = [];
n = 1;
s(n).Field = 'Images';
s(n).Label = 'Images';
s(n).Default = true;
n = n+1;

if ~issegment3dp
  s(n).Field = 'LVEndocardium';
  s(n).Label = 'LV Endocardium';
  s(n).Default = false;
  n = n+1;

  s(n).Field = 'LVEpicardium';
  s(n).Label = 'LV Epicardium';
  s(n).Default = false;
  n = n+1;

  s(n).Field = 'RVEndocardium';
  s(n).Label = 'RV Endocardium';
  s(n).Default = false;
  n = n+1;

  s(n).Field = 'RVEpicardium';
  s(n).Label = 'RV Epicardium';
  s(n).Default = false;
  n = n+1;

  s(n).Field = 'Myocardium';
  s(n).Label = 'Myocardium and blood pool';
  s(n).Default = false;
  n = n+1;

  s(n).Field = 'FlowROIs';
  s(n).Label = 'Flow ROIs';
  s(n).Default = false;
  n = n+1;

  s(n).Field = 'Points';
  s(n).Label = 'Points';
  s(n).Default = false;
  n = n+1;

  s(n).Field = 'EDESTime';
  s(n).Label = 'ED/ES Time Only (for ALL expect Points and Flow)';
  s(n).Default = false;
  %n = n+1;

else
  %Segment 3DPrint

  s(n).Field = 'Objects';
  s(n).Label = 'Objects';
  s(n).Default = true;
  n = n+1;

  s(n).Field = 'Points';
  s(n).Label = 'Points';
  s(n).Default = false;
  %n = n+1;

end

%Call to ask
[outs,ok] = myinputstruct(s,'Select what to export',10); %10 is 10 characters minimum

if not(ok)
  return
end

d.doImages = getfieldhelper(outs,'Images');
d.doLVEndocardium = getfieldhelper(outs,'LVEndocardium');
d.doLVEpicardium = getfieldhelper(outs,'LVEpicardium');
d.doRVEndocardium = getfieldhelper(outs,'RVEndocardium');
d.doRVEpicardium = getfieldhelper(outs,'RVEpicardium');
d.doMyocardium = getfieldhelper(outs,'Myocardium');
d.doFlowROIs = getfieldhelper(outs,'FlowROIs');
d.doPoints = getfieldhelper(outs,'Points');
d.doObjects = getfieldhelper(outs,'Objects');
d.doEDESTime = getfieldhelper(outs,'EDESTime');

%-----------------------------------
function t = getfieldhelper(s,fname)
%-----------------------------------
%Helper function to extract

if isfield(s,fname)
  t = s.(fname);
else
  t = false;
end


%--------------------------------
function savenifitii_helper(vol,no,filename, doEDES)
%----------------------------------
%Funtion that wraps all needed functions to create NIFTI
%It creates volumes and saves nifit according to given variables
%Type - type of contours, data,
%no - number of stack
%doEDES - only images with ED and ES markers
global SET

if nargin < 4
  doEDES = false;
end


if doEDES
  edestime= [SET(no).EDT, SET(no).EST];
  edestext{SET(no).EDT} = 'ED'; %it will be added at the end of the file's name
  edestext{SET(no).EST} = 'ES';

  if SET(no).TSize == 1
    edestext='notcine';%it is no cine, ED and ES markers are not existing in that case
  end

  for tloop = edestime
    volcurrent=(vol(:,:,tloop,:));
    imagefilename = [filename '_' edestext{tloop}];
    volcurrent = volumefixer(volcurrent,no);
    info = getinfostruct(volcurrent,no,doEDES);
    niftiwrite(volcurrent,imagefilename,info);
    logdisp(sprintf('Wrote file %s',imagefilename));
  end
else
  %imagefilename = [filename];
  vol = volumefixer(vol,no);
  info = getinfostruct(vol,no);

  niftiwrite(vol,filename,info);
  logdisp(sprintf('Wrote file %s',filename));
end


