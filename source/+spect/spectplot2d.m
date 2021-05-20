function [varargout] = spectplot2d(varargin)
%-------------------------------------------
%Plot a 2D view of the left ventricle for a spect image
%written by Helen Soneson 2008-03-03

if nargin == 0
  varargin{1} = 'init_Callback';
end

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:});

%----------------------------------------
function init_Callback(mode,makeviewmode)
%----------------------------------------
%open the 2D visualization view for either
%mode=visualization: stress and rest ungated visualization
%mode=scoring: stress and rest for manual scoring of ischemia and infarct

global SET DATA

if nargin < 2
  makeviewmode = false;
end

%find stress and rest image stack number
[restno,stressno,restgatedno,stressgatedno] = spect.findimagestacks;

switch mode
  case 'visualization'
    if isempty(restno) && isempty(stressno) && isempty(restgatedno) && isempty(stressgatedno)
      myfailed('Could not find Perfusion Rest or Perfusion Stress image stacks');
      return;
    end
  case 'scoring'
    if isempty(restno) && isempty(stressno)
      myfailed('Missing one or more of Perfusion Rest and Perfusion Stress image stacks');
      return;
    end
end

%find unique image stacks
restno = spect.specttools('uniqueimagestack',restno,'Rest');
stressno = spect.specttools('uniqueimagestack',stressno,'Stress');
restgatedno = spect.specttools('uniqueimagestack',restgatedno,'Rest gated');
stressgatedno = spect.specttools('uniqueimagestack',stressgatedno,'Stress gated');
nos = [restno stressno restgatedno stressgatedno];

%check patientID
if nos > 1 %#ok<BDSCI>
  for a = 1:length(nos)
    for b = a+1:length(nos)
      if ~strcmp(SET(nos(a)).PatientInfo.ID,SET(nos(b)).PatientInfo.ID)
        myfailed('Different Patient ID for the image stacks',DATA.GUI.Segment);
        return;
      end
    end
  end
end

nos = [];
nos.rest = restno;
nos.stress = stressno;
nos.restgated = restgatedno;
nos.stressgated = stressgatedno;
fnames = fieldnames(nos);

h = mywaitbarstart(5,'Calculating',1,DATA.GUI.Segment);

IM = [];
LVsegmentation = [];
EndoX = [];
EndoY = [];
EpiX = [];
EpiY = [];
%crop the images
for i = 1:4
  field = fnames{i};
  no = nos.(field);
  if ~isempty(no)
    if ~makeviewmode
      [IM.(field),LVsegmentation.(field),EndoX.(field),EndoY.(field),EpiX.(field),EpiY.(field)] = calccrop(no);
    else
      IM.(field) = SET(no).IM;
      LVsegmentation.(field) = ~isempty(SET(no).EndoX);
      EndoX.(field) = SET(no).EndoX;
      EndoY.(field) = SET(no).EndoY;
      EpiX.(field) = SET(no).EpiX;
      EpiY.(field) = SET(no).EpiY;
    end
  end
end

%check so LV segmentation exist for all image stacks
for i = 1:4
  field = fnames{i};
  no = nos.(field);
  if ~isempty(no) && isempty(SET(no).EndoX)
    myfailed(sprintf('Need LV segmentation in %s image stack',field));
  return;
  end
end

%define slices to include in the visualization
[isLVseg,startslices,endslices,SAslices,HLAslice,VLAslice] = deal([]);
for i = 1:4
  field = fnames{i};
  no = nos.(field);
  if ~isempty(no)
    if LVsegmentation.(field)
      %find time frames and slices with segmentation
      isLVseg.(field) = find(sum(~isnan(squeeze(EpiX.(field)(1,:,:)))') > 0);
      [startslices.(field),endslices.(field),SAslices.(field), ...
        HLAslice.(field),VLAslice.(field)] = ...
        calcslicesfromsegmentation( ...
        EpiX.(field),EpiY.(field),isLVseg.(field),mode);
    else
      %finding marked slices
      [startslices.(field),endslices.(field),SAslices.(field), ...
        HLAslice.(field),VLAslice.(field)] = ...
        calcslicesfromstack(IM.(field),no);
      isLVsegrest = [];
    end
  end
end

%calculate number of time frames in the upsampled image stacks
if ~isempty(nos.restgated) || ~isempty(nos.stressgated)
  if ~isempty(nos.restgated) && ~isempty(nos.stressgated)
    if SET(nos.restgated).TSize == SET(nos.stressgated).TSize
      bothtf = SET(nos.restgated).TSize;
      if bothtf < 11
        nbroftimeframes = floor(15/SET(nos.restgated).TSize)*SET(nos.restgated).TSize;%-1;
      else
        nbroftimeframes = bothtf;
      end
    else
      nbroftimeframes = (SET(nos.restgated).TSize-1)*(SET(nos.stressgated).TSize-1)+1;
      bothtf = nbroftimeframes;
    end
  elseif ~isempty(nos.restgated)
    nbroftimeframes = SET(nos.restgated).TSize;
    bothtf = nbroftimeframes;
  elseif ~isempty(nos.stressgated)
    nbroftimeframes = SET(nos.stressgated).TSize;
    bothtf = nbroftimeframes;
  end
else
  nbroftimeframes = 1;
  bothtf = 1;
end

h = mywaitbarupdate(h);
  
[resamplexy,resamplez,SAapicalimage,SAbasalimage,SAmidimage, ...
  HLAimage,VLAimage,HLAmyocardiumx,HLAmyocardiumy,VLAmyocardiumx, ...
  VLAmyocardiumy,maxcount,mincount] = deal([]);
for i = 1:4
  field = fnames{i};
  no = nos.(field);
  if ~isempty(no)
    tnbr = 1;
    if strfind(field,'gated')
      tnbr = nbroftimeframes;
    end
    if ~makeviewmode
      [IM.(field), ...
        EndoX.(field),EndoY.(field),EpiX.(field),EpiY.(field), ...
        isLVseg.(field), startslices.(field),endslices.(field), ...
        SAslices.(field),HLAslice.(field),VLAslice.(field), ...
        resamplexy.(field),resamplez.(field)] = ...
        resamplecrop(no,IM.(field),SAslices.(field),HLAslice.(field), ...
        VLAslice.(field), LVsegmentation.(field),EndoX.(field), ...
        EndoY.(field),EpiX.(field),EpiY.(field), ...
        startslices.(field),endslices.(field),isLVseg.(field),tnbr,mode);
    else
      resamplexy.(field) = mean([SET(no).ResolutionX SET(no).ResolutionY]);
    end
    
    [IM.(field),SAapicalimage.(field),SAbasalimage.(field), ...
      SAmidimage.(field),HLAimage.(field),VLAimage.(field), ...
      EndoX.(field),EndoY.(field),EpiX.(field),EpiY.(field), ...
      isLVseg.(field), ...
      startslices.(field),endslices.(field),SAslices.(field), ...
      HLAslice.(field),VLAslice.(field),HLAmyocardiumx.(field), ...
      HLAmyocardiumy.(field),VLAmyocardiumx.(field), ...
      VLAmyocardiumy.(field), maxcount.(field),mincount.(field)] = ...
      initimages(IM.(field),SAslices.(field),HLAslice.(field), ...
      VLAslice.(field),LVsegmentation.(field),EndoX.(field), ...
      EndoY.(field),EpiX.(field),EpiY.(field), startslices.(field), ...
      endslices.(field),isLVseg.(field),tnbr,mode,resamplexy.(field));
  end 
  h = mywaitbarupdate(h);
end

mywaitbarclose(h);

if makeviewmode
  for i = 1:4
    field = fnames{i};
    no = nos.(field);
    if ~isempty(no)
      sax3 = [];
      sax3.IM = cat(4, ...
        SAbasalimage.(field),SAmidimage.(field),SAapicalimage.(field));
      sax3.slices = flipud(cell2mat(SAslices.(field)));
      sax3.EndoX = EndoX.(field)(:,:,sax3.slices);
      sax3.EndoY = EndoY.(field)(:,:,sax3.slices);
      sax3.EpiX = EpiX.(field)(:,:,sax3.slices);
      sax3.EpiY = EpiY.(field)(:,:,sax3.slices);
      SET(no).SAX3 = sax3;
      
      hla = [];
      vla = [];
      hla.IM = HLAimage.(field);
      hla.slice = HLAslice.(field);
      hla.EndoX = HLAmyocardiumx.(field);
      hla.EndoY = HLAmyocardiumy.(field);
      vla.IM = VLAimage.(field);
      vla.slice = VLAslice.(field);
      vla.EndoX = VLAmyocardiumx.(field);
      vla.EndoY = VLAmyocardiumy.(field);
      SET(no).HLA = hla;
      SET(no).VLA = vla;
    end
  end
  return
end

%%%%%%%%%% OPEN THE GUI %%%%%%%%%

if isempty(restno) && isempty(stressno) && isempty(restgatedno) && isempty(stressgatedno)
  myfailed('No rest or stress perfusion image exist.');
  return
end
if isopengui('spectplot2d.fig');%['+spect' filesep 'spectplot2d.fig'])
  gui = DATA.GUI.SpectPlot2d;
  figure(gui.fig);
else
  DATA.GUI.SpectPlot2d = mygui(['+spect' filesep 'spectplot2d.fig']);
  gui = DATA.GUI.SpectPlot2d;
  myadjust(gui.fig,DATA.GUI.Segment);
end
gui.handles.tf = 1;
gui.mode = mode;
gui.bothtf = bothtf;
gui.nbroftimeframes = nbroftimeframes;
if ~isempty(nos.rest)
  gui.restno = nos.rest;
  gui.restIM = IM.rest;
  gui.LVsegmentationrest = LVsegmentation.rest;
else
  gui.restno = [];
  gui.restIM = [];
  gui.LVsegmentationrest = [];
end
if ~isempty(stressno)
  gui.stressno = nos.stress;
  gui.stressIM = IM.stress;
  gui.LVsegmentationstress = LVsegmentation.stress;
else
  gui.stressno = [];
  gui.stressIM = [];
  gui.LVsegmentationstress = [];
end
if ~isempty(restgatedno)
  gui.restgatedno = nos.restgated;
  gui.LVsegmentationrestgated = LVsegmentation.restgated;
else
  gui.restgatedno = [];
  gui.LVsegmentationrestgated = [];
end
if ~isempty(stressgatedno)
  gui.stressgatedno = nos.stressgated;
  gui.LVsegmentationstressgated = LVsegmentation.stressgated;
else
  gui.stressgatedno = [];
  gui.LVsegmentationstressgated = [];
end

if ~isempty(nos.rest)
  gui.handles.StringRest = sprintf(['LVM: ',num2str(round(SET(nos.rest).LVM*1.05)),' g\nEDV: ',num2str(round(SET(nos.rest).EDV)),' ml']);
else
  gui.handles.StringRest = [];
end
if ~isempty(nos.stress)
  gui.handles.StringStress = sprintf(['LVM: ',num2str(round(SET(nos.stress).LVM*1.05)),' g\nEDV: ',num2str(round(SET(nos.stress).EDV)),' ml']);
else
  gui.handles.StringStress = [];
end
if ~isempty(nos.restgated)
  gui.handles.StringRestgated = sprintf(['LVM: ',num2str(round(mean(SET(nos.restgated).LVM*1.05))),' g\nEDV: ',num2str(round(SET(nos.restgated).EDV)),' ml\nEF: ',num2str(round(100*SET(nos.restgated).EF)),' %']);
else
  gui.handles.StringRestgated = [];
end
if ~isempty(nos.stressgated)
  gui.handles.StringStressgated = sprintf(['LVM: ',num2str(round(mean(SET(nos.stressgated).LVM*1.05))),' g\nEDV: ',num2str(round(SET(nos.stressgated).EDV)),' ml\nEF: ',num2str(round(100*SET(nos.stressgated).EF)),' %']);
else
  gui.handles.StringStressgated = [];
end

switch mode
  case 'visualization'
    if ~isempty(restno) || ~isempty(stressno)
      if ~isempty(restno) && ~isempty(stressno)
        gui.showno = [stressno restno];
        gui.showimagetype = {'stress', 'rest'};
      elseif ~isempty(stressno)
        gui.showno = [stressno 0];
        gui.showimagetype = {'stress', ''};
      else
        gui.showno = [0 restno];
        gui.showimagetype = {'', 'rest'};
      end
      if ~isempty(restgatedno) || ~isempty(stressgatedno)
        set(gui.handles.radiobuttonungated,'value',1);
        set(gui.handles.radiobuttongated,'value',0);
      else
        set(gui.handles.radiobuttonungated,'enable','off');
        set(gui.handles.radiobuttongated,'enable','off');
      end
    else
      set(gui.handles.radiobuttonungated,'enable','off');
      set(gui.handles.radiobuttongated,'enable','off');      
      if ~isempty(restgatedno) && ~isempty(stressgatedno)
        gui.showno = [stressgatedno restgatedno];
        gui.showimagetype = {'stressgated','restgated'};
      elseif ~isempty(stressgatedno)
        gui.showno = [stressgatedno 0];
        gui.showimagetype = {'stressgated', ''};
      else
        gui.showno = [0 restgatedno];
        gui.showimagetype = {'', 'restgated'};
      end
    end
  case 'scoring'
    if ~isempty(stressno) && ~isempty(restno), gui.showno = [stressno restno]; gui.showimagetype = {'stress', 'rest'}; end
    if isempty(stressno) && ~isempty(restno), gui.showno = [0 restno]; gui.showimagetype = {'', 'rest'}; end
    if ~isempty(stressno) && isempty(restno), gui.showno = [stressno 0]; gui.showimagetype = {'stress', ''}; end
    if ~isempty(restgatedno)
      set(gui.handles.radiobuttonungated,'enable','on');
      set(gui.handles.radiobuttongated,'enable','on');
    else
      set(gui.handles.radiobuttonungated,'enable','off');
      set(gui.handles.radiobuttongated,'enable','off');
    end
    set(gui.handles.radiobuttonungated,'value',1);
    set(gui.handles.radiobuttongated,'value',0);
end  

gui.framerate = 0.5;
set(gui.handles.sliderframerate,'value',0.5);
set(gui.handles.radiobuttonnormeachtf,'value',1);  %0
set(gui.handles.radiobuttonnormall,'value',0);  %1
if ~get(gui.handles.radiobuttongated,'value') 
  set(gui.handles.pushbuttonnext,'enable','off');
  set(gui.handles.pushbuttonprev,'enable','off');
  set(gui.handles.pushbuttonplay,'enable','off');
  set(gui.handles.pushbuttonstop,'enable','off');
  set(gui.handles.radiobuttonnormeachtf,'enable','off');
  set(gui.handles.radiobuttonnormall,'enable','off');
  set(gui.handles.sliderframerate,'enable','off');
end
if isequal(mode,'scoring')
  set(gui.handles.pushbuttonexit,'String','Save and Exit');
else
  set(gui.handles.pushbuttonexit,'String','Save and Exit');
end
%keypressed
set(gui.fig,'keypressfcn',@keypressed);
allchildren = get(gui.fig,'children');
for loopa = 1:length(allchildren)
  try
    set(allchildren(loopa),'keypressfcn',@keypressed);
  catch
  end
  allgrandchildren = get(allchildren(loopa),'children');
  for loopb = 1:length(allgrandchildren)
    try
      set(allgrandchildren(loopb),'keypressfcn',@keypressed);
    catch
    end
  end
end
%default scoring values
if isequal(mode,'scoring')
  gui.scoringvalues = cell(3,17);
  for k = 1:size(gui.scoringvalues,2)
    if ~isempty(stressno) && isfield(SET(stressno).Perfusion,'MPS') && isfield(SET(stressno).Perfusion.MPS,'ScoringIschManual') && ~isempty(SET(stressno).Perfusion.MPS.ScoringIschManual)
      gui.scoringvalues{1,k} = SET(stressno).Perfusion.MPS.ScoringIschManual(k);
    else
      gui.scoringvalues{1,k} = 0;
    end
    if ~isempty(restno) && isfield(SET(restno).Perfusion,'MPS') && isfield(SET(restno).Perfusion.MPS,'ScoringIschManual') && ~isempty(SET(restno).Perfusion.MPS.ScoringIschManual)
      gui.scoringvalues{2,k} = SET(restno).Perfusion.MPS.ScoringIschManual(k);
    else
      gui.scoringvalues{2,k} = 0;
    end
    if ~isempty(restno) && isfield(SET(restno).Perfusion,'MPS') && isfield(SET(restno).Perfusion.MPS,'ScoringFixManual') && ~isempty(SET(restno).Perfusion.MPS.ScoringFixManual)
      gui.scoringvalues{3,k} = SET(restno).Perfusion.MPS.ScoringFixManual(k);
    else
      gui.scoringvalues{3,k} = 0;
    end
  end
end

%%%%%%%%% BOTH %%%%%%%%%%

if (~isempty(gui.LVsegmentationrest) && gui.LVsegmentationrest) || ...
    (~isempty(gui.LVsegmentationstress) && gui.LVsegmentationstress) || ...
    (~isempty(gui.LVsegmentationrestgated) && gui.LVsegmentationrestgated) || ...
    (~isempty(gui.LVsegmentationstressgated) && gui.LVsegmentationstressgated)
  set(gui.handles.radiobuttonshowsegmentation,'value',1);
  set(gui.handles.radiobuttonhidesegmentation,'value',0);
else
  set(gui.handles.radiobuttonshowsegmentation,'enable','off');
  set(gui.handles.radiobuttonhidesegmentation,'enable','off');
end


%%%%%%%%% REST %%%%%%%%%%

if ~isempty(restno)
  %images
  gui.SAbasalimagerest = SAbasalimage.rest;
  gui.SAmidimagerest = SAmidimage.rest;
  gui.SAapicalimagerest = SAapicalimage.rest;
  gui.HLAimagerest = HLAimage.rest;
  gui.VLAimagerest = VLAimage.rest;
  %segmentations
  gui.EndoXrest = EndoX.rest;
  gui.EndoYrest = EndoY.rest;
  gui.EpiXrest = EpiX.rest;
  gui.EpiYrest = EpiY.rest;
  gui.HLAmyocardiumxrest = HLAmyocardiumx.rest;
  gui.HLAmyocardiumyrest = HLAmyocardiumy.rest;
  gui.VLAmyocardiumxrest = VLAmyocardiumx.rest;
  gui.VLAmyocardiumyrest = VLAmyocardiumy.rest;
  %slices
  gui.startslicesrest = startslices.rest;
  gui.endslicesrest = endslices.rest; 
  gui.SAslicesrest = SAslices.rest;
  gui.HLAslicerest = HLAslice.rest;
  gui.VLAslicerest = VLAslice.rest;
  gui.isLVsegrest = isLVseg.rest;
  %min- and maxcounts
  gui.maxcountrest = maxcount.rest;
  gui.mincountrest = mincount.rest; 
  %resample factors
  gui.resamplexyrest = resamplexy.rest;
  gui.resamplezrest = resamplez.rest;  
end

if ~isempty(restgatedno)
  %images
  gui.SAbasalimagerestgated = SAbasalimage.restgated;
  gui.SAmidimagerestgated = SAmidimage.restgated;
  gui.SAapicalimagerestgated = SAapicalimage.restgated;
  gui.HLAimagerestgated = HLAimage.restgated;
  gui.VLAimagerestgated = VLAimage.restgated;
  %segmentations
  gui.EndoXrestgated = EndoX.restgated;
  gui.EndoYrestgated = EndoY.restgated;
  gui.EpiXrestgated = EpiX.restgated;
  gui.EpiYrestgated = EpiY.restgated;
  gui.HLAmyocardiumxrestgated = HLAmyocardiumx.restgated;
  gui.HLAmyocardiumyrestgated = HLAmyocardiumy.restgated;
  gui.VLAmyocardiumxrestgated = VLAmyocardiumx.restgated;
  gui.VLAmyocardiumyrestgated = VLAmyocardiumy.restgated;
  %slices
  gui.startslicesrestgated = startslices.restgated;
  gui.endslicesrestgated = endslices.restgated;
  gui.SAslicesrestgated = SAslices.restgated;
  gui.HLAslicerestgated = HLAslice.restgated;
  gui.VLAslicerestgated = VLAslice.restgated;
  gui.isLVsegrestgated = isLVseg.restgated;
  %min- and maxcounts
  gui.maxcountrestgated = maxcount.restgated;
  gui.mincountrestgated = mincount.restgated; 
  %time frames
  originaltf = linspace(1,gui.nbroftimeframes,SET(gui.restgatedno).TSize);
  gui.EDTrestgated = originaltf(SET(restgatedno).EDT);
  gui.ESTrestgated = originaltf(SET(restgatedno).EST);
  %resample factors
  gui.resamplexyrestgated = resamplexy.restgated;
  gui.resamplezrestgated = resamplez.restgated;  
end

%%%%%%%%% STRESS %%%%%%%%%%

if ~isempty(stressno)
  %images
  gui.SAbasalimagestress = SAbasalimage.stress;
  gui.SAmidimagestress = SAmidimage.stress;
  gui.SAapicalimagestress = SAapicalimage.stress;
  gui.HLAimagestress = HLAimage.stress;
  gui.VLAimagestress = VLAimage.stress;
  %segmentations
  gui.EndoXstress = EndoX.stress;
  gui.EndoYstress = EndoY.stress;
  gui.EpiXstress = EpiX.stress;
  gui.EpiYstress = EpiY.stress;
  gui.HLAmyocardiumxstress = HLAmyocardiumx.stress;
  gui.HLAmyocardiumystress = HLAmyocardiumy.stress;
  gui.VLAmyocardiumxstress = VLAmyocardiumx.stress;
  gui.VLAmyocardiumystress = VLAmyocardiumy.stress;
  %slices
  gui.startslicesstress = startslices.stress;
  gui.endslicesstress = endslices.stress; 
  gui.SAslicesstress = SAslices.stress;
  gui.HLAslicestress = HLAslice.stress;
  gui.VLAslicestress = VLAslice.stress;
  gui.isLVsegstress = isLVseg.stress;
  %min- and maxcounts
  gui.maxcountstress = maxcount.stress;
  gui.mincountstress = mincount.stress; 
  %resample factors
  gui.resamplexystress = resamplexy.stress;
  gui.resamplezstress = resamplez.stress;  
end

if ~isempty(stressgatedno)
  %images
  gui.SAbasalimagestressgated = SAbasalimage.stressgated;
  gui.SAmidimagestressgated = SAmidimage.stressgated;
  gui.SAapicalimagestressgated = SAapicalimage.stressgated;
  gui.HLAimagestressgated = HLAimage.stressgated;
  gui.VLAimagestressgated = VLAimage.stressgated;
  %segmentations
  gui.EndoXstressgated = EndoX.stressgated;
  gui.EndoYstressgated = EndoY.stressgated;
  gui.EpiXstressgated = EpiX.stressgated;
  gui.EpiYstressgated = EpiY.stressgated;
  gui.HLAmyocardiumxstressgated = HLAmyocardiumx.stressgated;
  gui.HLAmyocardiumystressgated = HLAmyocardiumy.stressgated;
  gui.VLAmyocardiumxstressgated = VLAmyocardiumx.stressgated;
  gui.VLAmyocardiumystressgated = VLAmyocardiumy.stressgated;
  %slices
  gui.startslicesstressgated = startslices.stressgated;
  gui.endslicesstressgated = endslices.stressgated;
  gui.SAslicesstressgated = SAslices.stressgated;
  gui.HLAslicestressgated = HLAslice.stressgated;
  gui.VLAslicestressgated = VLAslice.stressgated;
  gui.isLVsegstressgated = isLVseg.stressgated;
  %min- and maxcounts
  gui.maxcountstressgated = maxcount.stressgated;
  gui.mincountstressgated = mincount.stressgated; 
  %time frames
  originaltf = linspace(1,gui.nbroftimeframes,SET(gui.stressgatedno).TSize);
  gui.EDTstressgated = originaltf(SET(stressgatedno).EDT);
  gui.ESTstressgated = originaltf(SET(stressgatedno).EST);
  %resample factors
  gui.resamplexystressgated = resamplexy.stressgated;
  gui.resamplezstressgated = resamplez.stressgated;  
end

%%%%%%%%%% BOTH %%%%%%%%%%%%

plotimages;
plotLVsegmentation;

switch gui.mode
  case 'visualization'
    switch gui.showimagetype{1}
      case 'stress'
        set(gui.handles.textvolumesleft,'String',gui.handles.StringStress);
      case 'stressgated'
        set(gui.handles.textvolumesleft,'String',gui.handles.StringStressgated);
      case 'restgated'
        set(gui.handles.textvolumesleft,'String',gui.handles.StringRestgated);
    end
    switch gui.showimagetype{2}
      case 'rest'
        set(gui.handles.textvolumesright,'String',gui.handles.StringRest);
      case 'restgated'
        set(gui.handles.textvolumesright,'String',gui.handles.StringRestgated);
    end
    plotSAintersections('both');
    plotahasections('none');
    hideahavalues('both');
  case 'scoring'
    if ~get(gui.handles.radiobuttongated,'value')
      plotSAintersections('none');
      if ~isempty(gui.showimagetype{1}) && ~isempty(gui.showimagetype{2})
        plotahavalues('both');
        if get(gui.handles.radiobuttonshowsegmentation,'value')
          plotahasections('both');
        else
          plotahasections('none');
        end
      elseif ~isempty(gui.showimagetype{1})
        plotahavalues('left')
        hideahavalues('right');
        if get(gui.handles.radiobuttonshowsegmentation,'value')
          plotahasections('left');
        else
          plotahasections('none');
        end
      elseif ~isempty(gui.showimagetype{2})
        plotahavalues('right');
        hideahavalues('left');
        if get(gui.handles.radiobuttonshowsegmentation,'value')
          plotahasections('right');
        else
          plotahasections('none');
        end
      end
    else
      plotSAintersections('left');
      plotahavalues('right'); 
      hideahavalues('left');     
      if get(gui.handles.radiobuttonshowsegmentation,'value')
        plotahasections('right');
      else
        plotahasections('none');
      end
    end
end

if isequal(gui.mode,'scoring')
  set(get(gui.handles.axesHLAright,'children'),'ButtonDownFcn','spect.spectplot2d(''moveSAline'',''down'',''HLAright'')');
  set(get(gui.handles.axesVLAright,'children'),'ButtonDownFcn','spect.spectplot2d(''moveSAline'',''down'',''VLAright'')');
  set(get(gui.handles.axesHLAleft,'children'),'ButtonDownFcn','spect.spectplot2d(''moveSAline'',''down'',''HLAleft'')');
  set(get(gui.handles.axesVLAleft,'children'),'ButtonDownFcn','spect.spectplot2d(''moveSAline'',''down'',''VLAleft'')');
  gui.LAline.xstart = [];
  gui.LAline.ystart = [];
  gui.LAline.xend = [];
  gui.LAline.yend = [];
  gui.LAline.imagestart = [];
end
plotvolumecurve;


%----------------------------------------------------------------
function [IM,LVsegmentation,EndoX,EndoY,EpiX,EpiY] = calccrop(no)
%----------------------------------------------------------------
%calculate crop size
global SET

if isempty(no)
  LVsegmentation = 0;
  cropsizemm = [];
else
  %crop images and segmentations
  cropframesizemm = 20;
  %set cropped size
  if isempty(SET(no).EndoX) && isempty(SET(no).EpiX)
    LVsegmentation = 0;
  else
    LVsegmentation = 1;
  end
  if LVsegmentation
    xcenter = mynanmean(SET(no).EndoX(:));
    ycenter = mynanmean(SET(no).EndoY(:));
    xminborder = max([1 min(SET(no).EpiX(:))-cropframesizemm/SET(no).ResolutionX]);
    xmaxborder = min([SET(no).XSize max(SET(no).EpiX(:))+cropframesizemm/SET(no).ResolutionX]);
    yminborder = max([1 min(SET(no).EpiY(:))-cropframesizemm/SET(no).ResolutionY]);
    ymaxborder = min([SET(no).YSize max(SET(no).EpiY(:))+cropframesizemm/SET(no).ResolutionY]);
    cropboxlengthmm = min([(xmaxborder-xminborder)*SET(no).ResolutionX (ymaxborder-yminborder)*SET(no).ResolutionY]);
  else
    xcenter = SET(no).XSize/2;
    ycenter = SET(no).YSize/2;
    cropboxlengthmm = 140;
  end
  cropsizemm = min([cropboxlengthmm SET(no).XSize*SET(no).ResolutionX SET(no).YSize*SET(no).ResolutionY]);
end
%crop image
cropsizex = min([cropsizemm/SET(no).ResolutionX SET(no).XSize]);
cropsizey = min([cropsizemm/SET(no).ResolutionY SET(no).YSize]);
xmin = max([1 round(xcenter-cropsizex/2)-1]);
xmax = min([SET(no).XSize round(xcenter+cropsizex/2)+1]);
ymin = max([1 round(ycenter-cropsizey/2)-1]);
ymax = min([SET(no).YSize round(ycenter+cropsizey/2)+1]);
IM = SET(no).IM(xmin:xmax,ymin:ymax,:,:);
%crop segmentation
xofs = -xmin+1;
yofs = -ymin+1;
EndoX = SET(no).EndoX+xofs;
EndoY = SET(no).EndoY+yofs;
EpiX = SET(no).EpiX+xofs;
EpiY = SET(no).EpiY+yofs;

%--------------------------------------------------------------------------
function [IM, ...
  EndoX,EndoY,EpiX,EpiY,isLVseg, ...
  startslices,endslices,SAslices,HLAslice,VLAslice, ...
  resamplexy,resamplez] = ...
  resamplecrop(no,IM,SAslices,HLAslice,VLAslice, ...
  LVsegmentation,EndoX,EndoY,EpiX,EpiY,startslices,endslices,isLVseg,nbroftimeframes,mode)
%--------------------------------------------------------------------------
global SET

%resolution in the resampled image stack
resampledresolution = 1; %1.5;
%smoothing parameters
sigma = [1 1 1 1];%[2 2 1.3 2]; %
cropframesizemm = 20;

if ~isempty(no)
  %set upsampling factor
  if mean([SET(no).ResolutionX SET(no).ResolutionY]) < resampledresolution
    resamplexy = 1;
    fs = mean([SET(no).ResolutionX SET(no).ResolutionY]);
  else
    resamplexy = resampledresolution;
    fs = mean([SET(no).ResolutionX SET(no).ResolutionY])/resampledresolution;
  end
  if SET(no).SliceThickness+SET(no).SliceGap < resampledresolution
    resamplez = 1;
    fsz = SET(no).SliceThickness+SET(no).SliceGap;
  else
    resamplez = resampledresolution;
    fsz = (SET(no).SliceThickness+SET(no).SliceGap)/resampledresolution;
  end
  if SET(no).TSize > 1
    ft = nbroftimeframes/SET(no).TSize;
  else
    ft = 1;
  end
  %upsample images spatially
  IM = tools('upsampleslices',fsz,IM,1);
  IM = tools('upsamplevolume2',[fs fs 1 1],IM,1);
  %upsample images temporally
  if ft ~= 1
    IM = tools('upsampletemporal',ft,IM);
  end

  %adjust segmentation for upsampling
  if LVsegmentation
    EndoX = adjustsegmentation(fs,fsz,ft,EndoX,startslices,endslices,isLVseg);
    EndoY = adjustsegmentation(fs,fsz,ft,EndoY,startslices,endslices,isLVseg);
    EpiX = adjustsegmentation(fs,fsz,ft,EpiX,startslices,endslices,isLVseg);
    [EpiY,isLVseg] = adjustsegmentation(fs,fsz,ft,EpiY,startslices,endslices,isLVseg);
  end
    
  %update slices after upsampling
  cropLAslices = round(cropframesizemm*fsz/(SET(no).SliceThickness+SET(no).SliceGap));
  if LVsegmentation
    [startslices,endslices,SAslices,HLAslice,VLAslice,cropstartslice,cropendslice] ...
      = calcslicesfromsegmentation(EpiX,EpiY,isLVseg,mode,cropLAslices,size(IM,4)); %,SAslices,fsz);
  else
    [startslices,endslices,SAslices,HLAslice,VLAslice,cropstartslice,cropendslice] ...
      = calcslicesfromstack(IM,no);
  end
  
  %crop and smooth image
  IM = IM(:,:,:,cropstartslice:cropendslice);
  [IM,cropx,cropy,cropz] = smoothimage(IM,sigma);
  cropx = min([floor(min(EpiX(:)))-1 cropx size(IM,1)-ceil(max(EpiX(:)))-1]);
  cropy = min([floor(min(EpiY(:)))-1 cropy size(IM,1)-ceil(max(EpiY(:)))-1]);
  cropz = min([min(startslices)-1 cropz cropendslice-max(endslices)]);  
  IM = IM(cropx+1:end-cropx,cropy+1:end-cropy,:,cropz+1:end-cropz);
  cropstartslice = cropstartslice+cropz;
  cropendslice = cropendslice-cropz;
  HLAslice = HLAslice-cropx;
  VLAslice = VLAslice-cropy;
  for k = 1:size(SAslices,1)
    for l = 1:size(SAslices,2)
      SAslices{k,l} = SAslices{k,l}-cropstartslice+1;
    end
  end
  startslices = startslices-cropstartslice+1;
  endslices = endslices-cropstartslice+1;
  EndoX = EndoX(:,:,cropstartslice:cropendslice)-cropx;
  EndoY = EndoY(:,:,cropstartslice:cropendslice)-cropy;
  EpiX = EpiX(:,:,cropstartslice:cropendslice)-cropx;
  EpiY = EpiY(:,:,cropstartslice:cropendslice)-cropy;
end
%--------------------------------------------------------------------------
function [IM,SAapicalimage,SAbasalimage,SAmidimage,HLAimage,VLAimage, ...
  EndoX,EndoY,EpiX,EpiY,isLVseg, ...
  startslices,endslices,SAslices,HLAslice,VLAslice, ...
  HLAmyocardiumx,HLAmyocardiumy,VLAmyocardiumx,VLAmyocardiumy, ...
  maxcount,mincount] = ...
  initimages(IM,SAslices,HLAslice,VLAslice, ...
  LVsegmentation,EndoX,EndoY,EpiX,EpiY,startslices,endslices,isLVseg,nbroftimeframes,mode,resamplexy)
%--------------------------------------------------------------------------

%initialise images and segmentations
switch mode
  case 'visualization'
    SAapicalimage = squeeze(IM(:,:,:,SAslices{1}));
    SAmidimage = squeeze(IM(:,:,:,SAslices{2}));
    SAbasalimage = squeeze(IM(:,:,:,SAslices{3}));
  case 'scoring'
    SAapicalimage = squeeze(mean(IM(:,:,:,SAslices{1}),4));
    SAmidimage = squeeze(mean(IM(:,:,:,SAslices{2}),4));
    SAbasalimage = squeeze(mean(IM(:,:,:,SAslices{3}),4));
end
HLAimage = flipdim((squeeze(permute(IM(HLAslice,:,:,:),[4 2 3 1]))),1);
VLAimage = squeeze(permute(IM(:,VLAslice,:,:),[1 4 3 2]));

%find intersection points (LV segmentation in LA)
if LVsegmentation
  nbrofslices = size(HLAimage,1);
  %find segmentation for long axis images
  [HLAmyocardiumx,HLAmyocardiumy,VLAmyocardiumx,VLAmyocardiumy] = ...
    findintersectionpoints(EndoX,EndoY,EpiX,EpiY, ...
    startslices,endslices,HLAslice,VLAslice,nbrofslices,nbroftimeframes,isLVseg,resamplexy);
else
  HLAmyocardiumx = [];
  HLAmyocardiumy = [];
  VLAmyocardiumx = [];
  VLAmyocardiumy = [];
end

%calculate max- and mincounts
if LVsegmentation
  %find mincounts in the image stack and maxcounts within LV
  SAbasalmask = zeros(size(SAbasalimage));
  SAmidmask = zeros(size(SAmidimage));
  SAapicalmask = zeros(size(SAapicalimage));
  HLAmask = zeros(size(HLAimage));
  VLAmask = zeros(size(VLAimage));
  for tloop = isLVseg
    if ~isnan(EndoY(1,tloop,round(mean(SAslices{3,tloop}))))
      endomask = poly2mask(EndoY(:,tloop,round(mean(SAslices{3,tloop}))),EndoX(:,tloop,round(mean(SAslices{3,tloop}))), ...
        size(SAbasalimage,1),size(SAbasalimage,2));
    else
      endomask = zeros(size(SAbasalimage,1),size(SAbasalimage,2));
    end
    epimask = poly2mask(EpiY(:,tloop,round(mean(SAslices{3,tloop}))),EpiX(:,tloop,round(mean(SAslices{3,tloop}))), ...
      size(SAbasalimage,1),size(SAbasalimage,2));
    SAbasalmask(:,:,tloop) = epimask-endomask;
    if ~isnan(EndoY(1,tloop,round(mean(SAslices{2,tloop}))))
      endomask = poly2mask(EndoY(:,tloop,round(mean(SAslices{2,tloop}))),EndoX(:,tloop,round(mean(SAslices{2,tloop}))), ...
        size(SAmidimage,1),size(SAmidimage,2));
    else
      endomask = zeros(size(SAmidimage,1),size(SAmidimage,2));
    end
    epimask = poly2mask(EpiY(:,tloop,round(mean(SAslices{2,tloop}))),EpiX(:,tloop,round(mean(SAslices{2,tloop}))), ...
      size(SAmidimage,1),size(SAmidimage,2));
    SAmidmask(:,:,tloop) = epimask-endomask;
    if ~isnan(EndoY(1,tloop,round(mean(SAslices{1,tloop}))))
      endomask = poly2mask(EndoY(:,tloop,round(mean(SAslices{1,tloop}))),EndoX(:,tloop,round(mean(SAslices{1,tloop}))), ...
        size(SAapicalimage,1),size(SAapicalimage,2));
    else
      endomask = zeros(size(SAapicalimage,1),size(SAapicalimage,2));
    end
    epimask = poly2mask(EpiY(:,tloop,round(mean(SAslices{1,tloop}))),EpiX(:,tloop,round(mean(SAslices{1,tloop}))), ...
      size(SAapicalimage,1),size(SAapicalimage,2));
    SAapicalmask(:,:,tloop) = epimask-endomask;
    %HLA
    myocardiummask = poly2mask(HLAmyocardiumy{tloop},HLAmyocardiumx{tloop}, ...
      size(HLAimage,1),size(HLAimage,2));
    HLAmask(:,:,tloop) = myocardiummask;
    %VLA
    myocardiummask = poly2mask(VLAmyocardiumy{tloop},VLAmyocardiumx{tloop}, ...
      size(VLAimage,1),size(VLAimage,2));
    VLAmask(:,:,tloop) = myocardiummask;
  end
  SAbasalmyocardium = SAbasalimage.*SAbasalmask;
  SAmidmyocardium = SAmidimage.*SAmidmask;
  SAapicalmyocardium = SAapicalimage.*SAapicalmask;
  HLAmyocardium = HLAimage.*HLAmask;
  VLAmyocardium = VLAimage.*VLAmask;
  maxcount = max([reshape(SAbasalmyocardium,size(SAbasalmyocardium,1)*size(SAbasalmyocardium,2),size(SAbasalmyocardium,3)); ...
    reshape(SAmidmyocardium,size(SAmidmyocardium,1)*size(SAmidmyocardium,2),size(SAmidmyocardium,3)); ...
    reshape(SAapicalmyocardium,size(SAapicalmyocardium,1)*size(SAapicalmyocardium,2),size(SAapicalmyocardium,3)); ...
    reshape(HLAmyocardium,size(HLAmyocardium,1)*size(HLAmyocardium,2),size(HLAmyocardium,3)); ...
    reshape(VLAmyocardium,size(VLAmyocardium,1)*size(VLAmyocardium,2),size(VLAmyocardium,3))]);
  mincount = min([reshape(SAbasalimage,size(SAbasalimage,1)*size(SAbasalimage,2),size(SAbasalimage,3)); ...
    reshape(SAmidimage,size(SAmidimage,1)*size(SAmidimage,2),size(SAmidimage,3)); ...
    reshape(SAapicalimage,size(SAapicalimage,1)*size(SAapicalimage,2),size(SAapicalimage,3)); ...
    reshape(HLAimage,size(HLAimage,1)*size(HLAimage,2),size(HLAimage,3)); ...
    reshape(VLAimage,size(VLAimage,1)*size(VLAimage,2),size(VLAimage,3))]);
  if ~isempty(find(not(maxcount),1))
    maxcount = interp1(isLVseg,maxcount(isLVseg),1:length(maxcount));
    mincount = interp1(isLVseg,mincount(isLVseg),1:length(mincount));
  end
else
  maxcount = max([reshape(SAbasalimage,size(SAbasalimage,1)*size(SAbasalimage,2),size(SAbasalimage,3)); ...
    reshape(SAmidimage,size(SAmidimage,1)*size(SAmidimage,2),size(SAmidimage,3)); ...
    reshape(SAapicalimage,size(SAapicalimage,1)*size(SAapicalimage,2),size(SAapicalimage,3)); ...
    reshape(HLAimage,size(HLAimage,1)*size(HLAimage,2),size(HLAimage,3)); ...
    reshape(VLAimage,size(VLAimage,1)*size(VLAimage,2),size(VLAimage,3))]);
  mincount = min([reshape(SAbasalimage,size(SAbasalimage,1)*size(SAbasalimage,2),size(SAbasalimage,3)); ...
    reshape(SAmidimage,size(SAmidimage,1)*size(SAmidimage,2),size(SAmidimage,3)); ...
    reshape(SAapicalimage,size(SAapicalimage,1)*size(SAapicalimage,2),size(SAapicalimage,3)); ...
    reshape(HLAimage,size(HLAimage,1)*size(HLAimage,2),size(HLAimage,3)); ...
    reshape(VLAimage,size(VLAimage,1)*size(VLAimage,2),size(VLAimage,3))]);
  if ~isempty(find(not(maxcount),1))
    maxcount = interp1(isLVseg,maxcount(isLVseg),1:length(maxcount));
    mincount = interp1(isLVseg,mincount(isLVseg),1:length(mincount));
  end
end
%   mincount = min([min(SAbasalimage(:)) min(SAmidimage(:)) ...
%     min(SAapicalimage(:)) min(HLAimage(:)) min(VLAimage(:))]);



%--------------------------------------------------------------
function [smoothedim,cropx,cropy,cropz] = smoothimage(im,sigma)
%--------------------------------------------------------------
%smothing of the image im

kernelsize = ceil(3*sigma/2)*2+1;
f = fspecial('gaussian',kernelsize(1),sigma(1));
fx = f(:,ceil(length(f)/2));
fx = fx/sum(fx(:));
f = fspecial('gaussian',kernelsize(2),sigma(2));
fy = f(ceil(length(f)/2),:);
fy = fy/sum(fy(:));
f = fspecial('gaussian',kernelsize(3),sigma(3));
tempft = f(:,ceil(length(f)/2));
ft = reshape(tempft,[1 1 length(f) 1]);
ft = ft/sum(ft(:));
f = fspecial('gaussian',kernelsize(4),sigma(4));
tempfz = f(:,ceil(length(f)/2));
fz = reshape(tempfz,[1 1 1 length(f)]);
fz = fz/sum(fz(:));
temp = convn(single(im),single(fx),'same');
temp = convn(temp,single(fy),'same');
temp = convn(temp,single(ft),'same');
temp = convn(temp,single(fz),'same');
smoothedim = temp;
cropx = (kernelsize(1)-1)/2;
cropy = (kernelsize(2)-1)/2;
cropz = (kernelsize(4)-1)/2;


%--------------------------------------------------------------------------
function [HLAmyocardiumx,HLAmyocardiumy,VLAmyocardiumx,VLAmyocardiumy] = ...
  findintersectionpoints(endox,endoy,epix,epiy,startslices,endslices, ...
  HLAslice,VLAslice,nbrofslices,nbroftimeframes,isLVseg,resamplexy)
%--------------------------------------------------------------------------
%find the segmentation in the LA image projections
%input
%endox /endoy / epix / epiy:  segmentation in the SA image stack
%startslices:                 most apical slice in each time frame
%endslices:                   most basal slice in each time frame
%HLAslice:                    position for the HLA projection
%VLAslice:                    position for the VLA projection
%nbrofslices:                 number of SA slices in the whole image stack
%fsz:                         resampling parameter in the z-direction
%output
%HLAmyocardiumx / y:          myocardial border for the HLA projection
%VLAmyocardiumx / y:          myocardial border for the VLA projection
%LAbasalslice:                most basal slice in the croped image stack
%
%written by Helen Soneson 2010-06-02

m = max(endslices)-min(startslices)+1;
n = size(endox,1);
nerase = 3;
HLAendox = nan(m,nbroftimeframes,2);
HLAendoy = HLAendox;
HLAepix = HLAendox;
HLAepiy = HLAendox;
VLAendox = HLAendox;
VLAendoy = HLAendox;
VLAepix = HLAendox;
VLAepiy = HLAendox;
nbroftf = length(startslices);  % number of time frames in current image stack
for t = isLVseg
  tHLAepix = nan(endslices(t)-startslices(t)+1,2);
  tHLAepiy = tHLAepix;
  tVLAepix = tHLAepix;
  tVLAepiy = tHLAepix;
  tHLAendox = [];
  tHLAendoy = tHLAendox;
  tVLAendox = tHLAendox;
  tVLAendoy = tHLAendox;
  for z = startslices(t):endslices(t)    
    p = z-startslices(t)+1;
    HLAdistendo = HLAslice-endox(:,t,z);
    HLAdistepi = HLAslice-epix(:,t,z);
    VLAdistendo = VLAslice-endoy(:,t,z);
    VLAdistepi = VLAslice-epiy(:,t,z);
    %HLA endo
    [~,HLAdistendoindex] = min(abs(HLAdistendo));
    y1 = endoy(HLAdistendoindex,t,z);
    tempindex = [1:n 1:n 1:n];
    eraseindex = tempindex(HLAdistendoindex+n-nerase:HLAdistendoindex+n+nerase);
    HLAdistendo(eraseindex) = NaN;
    [~,HLAdistendoindex] = min(abs(HLAdistendo));
    y2 = endoy(HLAdistendoindex,t,z);
    if ~isnan(y1) && ~isnan(y2)
      if y1< y2
        tHLAendox(p,1) = nbrofslices-z;
        tHLAendoy(p,1) = y1;
        tHLAendox(p,2) = nbrofslices-z;
        tHLAendoy(p,2) = y2;
      else
        tHLAendox(p,1) = nbrofslices-z;
        tHLAendoy(p,1) = y2;
        tHLAendox(p,2) = nbrofslices-z;
        tHLAendoy(p,2) = y1;
      end
    end
    %HLA epi
    [~,HLAdistepiindex] = min(abs(HLAdistepi));
    y1 = epiy(HLAdistepiindex,t,z);    
    tempindex = [1:n 1:n 1:n];
    eraseindex = tempindex(HLAdistepiindex+n-nerase:HLAdistepiindex+n+nerase);
    HLAdistepi(eraseindex) = NaN;
    [~,HLAdistepiindex] = min(abs(HLAdistepi));
    y2 = epiy(HLAdistepiindex,t,z);
    if y1< y2
      tHLAepix(p,1) = nbrofslices-z;
      tHLAepiy(p,1) = y1;
      tHLAepix(p,2) = nbrofslices-z;
      tHLAepiy(p,2) = y2;
    else
      tHLAepix(p,1) = nbrofslices-z;
      tHLAepiy(p,1) = y2;
      tHLAepix(p,2) = nbrofslices-z;
      tHLAepiy(p,2) = y1;
    end
    %VLA endo
    [~,VLAdistendoindex] = min(abs(VLAdistendo));
    y1 = endox(VLAdistendoindex,t,z);
    tempindex = [1:n 1:n 1:n];
    eraseindex = tempindex(VLAdistendoindex+n-nerase:VLAdistendoindex+n+nerase);
    VLAdistendo(eraseindex) = NaN;
    [~,VLAdistendoindex] = min(abs(VLAdistendo));
    y2 = endox(VLAdistendoindex,t,z);
    if ~isnan(y1) && ~isnan(y2)
      if y1< y2
        tVLAendox(p,1) = y1;
        tVLAendoy(p,1) = z;
        tVLAendox(p,2) = y2;
        tVLAendoy(p,2) = z;
      else
        tVLAendox(p,1) = y2;
        tVLAendoy(p,1) = z;
        tVLAendox(p,2) = y1;
        tVLAendoy(p,2) = z;
      end
    end
    %VLA epi
    [~,VLAdistepiindex] = min(abs(VLAdistepi));
    y1 = epix(VLAdistepiindex,t,z);
    tempindex = [1:n 1:n 1:n];
    eraseindex = tempindex(VLAdistepiindex+n-nerase:VLAdistepiindex+n+nerase);
    VLAdistepi(eraseindex) = NaN;
    [~,VLAdistepiindex] = min(abs(VLAdistepi));
    y2 = epix(VLAdistepiindex,t,z);
    if y1 < y2
      tVLAepix(p,1) = y1;
      tVLAepiy(p,1) = z;
      tVLAepix(p,2) = y2;
      tVLAepiy(p,2) = z;
    else
      tVLAepix(p,1) = y2;
      tVLAepiy(p,1) = z;
      tVLAepix(p,2) = y1;
      tVLAepiy(p,2) = z;
    end
  end
  %erase outliers
  threshold = 3;
  xiendo = linspace(1,size(tHLAendox,1),m);
  xiepi = linspace(1,size(tHLAepix,1),m);
  %HLA
  distoutliers = max( ...
    [0 ; sqrt((tHLAendox(2:end-2,1)-tHLAendox(3:end-1,1)).^2+(tHLAendoy(2:end-2,1)-tHLAendoy(3:end-1,1)).^2) ; 0 ; 0], ...
    [0 ; 0 ; sqrt((tHLAendox(3:end-1,1)-tHLAendox(2:end-2,1)).^2+(tHLAendoy(3:end-1,1)-tHLAendoy(2:end-2,1)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  HLAendox(:,t,1) = interp1(inliersindex,tHLAendox(inliersindex,1),xiendo,'spline');
  HLAendoy(:,t,1) = interp1(inliersindex,tHLAendoy(inliersindex,1),xiendo,'spline');  
  distoutliers = max( ...
    [0 ; sqrt((tHLAendox(2:end-2,2)-tHLAendox(3:end-1,2)).^2+(tHLAendoy(2:end-2,2)-tHLAendoy(3:end-1,2)).^2) ; 0 ; 0], ...
    [0 ; 0 ; sqrt((tHLAendox(3:end-1,2)-tHLAendox(2:end-2,2)).^2+(tHLAendoy(3:end-1,2)-tHLAendoy(2:end-2,2)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  HLAendox(:,t,2) = interp1(inliersindex,tHLAendox(inliersindex,2),xiendo,'spline');
  HLAendoy(:,t,2) = interp1(inliersindex,tHLAendoy(inliersindex,2),xiendo,'spline');  
  distoutliers = max( ...
    [0 ; sqrt((tHLAepix(2:end-2,1)-tHLAepix(3:end-1,1)).^2+(tHLAepiy(2:end-2,1)-tHLAepiy(3:end-1,1)).^2) ; 0 ; 0], ...
    [0 ; 0 ; sqrt((tHLAepix(3:end-1,1)-tHLAepix(2:end-2,1)).^2+(tHLAepiy(3:end-1,1)-tHLAepiy(2:end-2,1)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  HLAepix(:,t,1) = interp1(inliersindex,tHLAepix(inliersindex,1),xiepi,'spline');
  HLAepiy(:,t,1) = interp1(inliersindex,tHLAepiy(inliersindex,1),xiepi,'spline');
  distoutliers = max( ...
    [0 ; sqrt((tHLAepix(2:end-2,2)-tHLAepix(3:end-1,2)).^2+(tHLAepiy(2:end-2,2)-tHLAepiy(3:end-1,2)).^2) ; 0 ; 0], ...
    [0 ; 0 ; sqrt((tHLAepix(3:end-1,2)-tHLAepix(2:end-2,2)).^2+(tHLAepiy(3:end-1,2)-tHLAepiy(2:end-2,2)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  HLAepix(:,t,2) = interp1(inliersindex,tHLAepix(inliersindex,2),xiepi,'spline');
  HLAepiy(:,t,2) = interp1(inliersindex,tHLAepiy(inliersindex,2),xiepi,'spline');
  
  %VLA
  distoutliers = max( ...
    [0 ; sqrt((tVLAendox(2:end-2,1)-tVLAendox(3:end-1,1)).^2+(tVLAendoy(2:end-2,1)-tVLAendoy(3:end-1,1)).^2) ; 0 ; 0], ...
    [0 ; 0 ; sqrt((tVLAendox(3:end-1,1)-tVLAendox(2:end-2,1)).^2+(tVLAendoy(3:end-1,1)-tVLAendoy(2:end-2,1)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  VLAendox(:,t,1) = interp1(inliersindex,tVLAendox(inliersindex,1),xiendo,'spline');
  VLAendoy(:,t,1) = interp1(inliersindex,tVLAendoy(inliersindex,1),xiendo,'spline');
  distoutliers = max( ...
    [0 ; sqrt((tVLAendox(2:end-2,2)-tVLAendox(3:end-1,2)).^2+(tVLAendoy(2:end-2,2)-tVLAendoy(3:end-1,2)).^2) ; 0 ; 0], ...
    [0 ; 0 ; sqrt((tVLAendox(3:end-1,2)-tVLAendox(2:end-2,2)).^2+(tVLAendoy(3:end-1,2)-tVLAendoy(2:end-2,2)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  VLAendox(:,t,2) = interp1(inliersindex,tVLAendox(inliersindex,2),xiendo,'spline');
  VLAendoy(:,t,2) = interp1(inliersindex,tVLAendoy(inliersindex,2),xiendo,'spline');  
  distoutliers = max( ...
    [0 ; sqrt((tVLAepix(2:end-2,1)-tVLAepix(3:end-1,1)).^2+(tVLAepiy(2:end-2,1)-tVLAepiy(3:end-1,1)).^2) ; 0 ; 0], ...
    [0 ; 0 ;sqrt((tVLAepix(3:end-1,1)-tVLAepix(2:end-2,1)).^2+(tVLAepiy(3:end-1,1)-tVLAepiy(2:end-2,1)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  VLAepix(:,t,1) = interp1(inliersindex,tVLAepix(inliersindex,1),xiepi,'spline');
  VLAepiy(:,t,1) = interp1(inliersindex,tVLAepiy(inliersindex,1),xiepi,'spline');
  distoutliers = max( ...
    [0 ; sqrt((tVLAepix(2:end-2,2)-tVLAepix(3:end-1,2)).^2+(tVLAepiy(2:end-2,2)-tVLAepiy(3:end-1,2)).^2) ; 0 ; 0], ...
    [0 ; 0 ; sqrt((tVLAepix(3:end-1,2)-tVLAepix(2:end-2,2)).^2+(tVLAepiy(3:end-1,2)-tVLAepiy(2:end-2,2)).^2) ; 0]);
  inliersindex = find(distoutliers < threshold);
  VLAepix(:,t,2) = interp1(inliersindex,tVLAepix(inliersindex,2),xiepi,'spline');
  VLAepiy(:,t,2) = interp1(inliersindex,tVLAepiy(inliersindex,2),xiepi,'spline');
end
% 
% HLAendox = 2*HLAendox;
% HLAepix = 2*HLAepix;
% % 
% VLAendoy = 2*VLAendoy;
% VLAepiy = 2*VLAepiy;


%%%%%%%%%% HLA %%%%%%%%%%

%erase endo in apex where there is no lumen and
%erase endo and epi in the outflow tract
lumendiameter = abs(HLAendoy(:,:,1)-HLAendoy(:,:,2));
wallthickness = abs(HLAendoy-HLAepiy);
lastendotemp = nan(nbroftf,1);
firstslice = nan(nbroftf,2);
lastepi = nan(nbroftf,1);
for t = isLVseg
  lastendotemp(t) = min(round(0.15*length(lumendiameter)),find(flipud(lumendiameter(:,t)) > 8,1,'first'));  % -> depend on LV lumen volume
  firstslice(t,1) = find(wallthickness(:,t,1) > 3/resamplexy,1,'first');  % -> depend on wallthickness for each time frame
  firstslice(t,2) = find(wallthickness(:,t,2) > 3/resamplexy,1,'first');
  lastepi(t) = find(~isnan(HLAepix(:,t,1)),1,'last');
end
lastendo = m-lastendotemp+1;
%put together HLA myocardium segmentation
HLAmyocardiumx = cell(1,nbroftf);
HLAmyocardiumy = HLAmyocardiumx;
for t = isLVseg
  apexsmoothingx = [HLAepix(lastepi(t),t,2)-1 HLAepix(lastepi(t),t,1)-1];
  myocardiumx = [HLAendox(firstslice(t,1):lastendo(t),t,1)' HLAendox(lastendo(t):-1:firstslice(t,2),t,2)' ...
    HLAepix(firstslice(t,2):lastepi(t),t,2)' apexsmoothingx HLAepix(lastepi(t):-1:firstslice(t,1),t,1)' HLAendox(firstslice(t,1),t,1)'];
  myocardiumx = interp1(linspace(0,1,length(myocardiumx)),myocardiumx,linspace(0,1,2*length(myocardiumx)));
  myocardiumx = [myocardiumx myocardiumx(1)];
  HLAmyocardiumx{t} = myocardiumx;
  apexsmoothingydist = abs(HLAepiy(lastepi(t),t,2)-HLAepiy(lastepi(t),t,1));
  apexsmoothingy = [HLAepiy(lastepi(t),t,2)-0.25*apexsmoothingydist HLAepiy(lastepi(t),t,1)+0.25*apexsmoothingydist];
  myocardiumy = [HLAendoy(firstslice(t,1):lastendo(t),t,1)' HLAendoy(lastendo(t):-1:firstslice(t,2),t,2)' ...
    HLAepiy(firstslice(t,2):lastepi(t),t,2)' apexsmoothingy HLAepiy(lastepi(t):-1:firstslice(t,1),t,1)' HLAendoy(firstslice(t,1),t,1)'];
  myocardiumy = interp1(linspace(0,1,length(myocardiumy)),myocardiumy,linspace(0,1,2*length(myocardiumy)));
  for k = 1:4  %smoothing
    myocardiumy = myocardiumy-0.2*(6*myocardiumy-circshift(myocardiumy,[0 1])-circshift(myocardiumy,[0 -1]) ...
      -circshift(myocardiumy,[0 2])-circshift(myocardiumy,[0 -2])-circshift(myocardiumy,[0 3])-circshift(myocardiumy,[0 -3]));
  end
  myocardiumy = [myocardiumy myocardiumy(1)];
  HLAmyocardiumy{t} = myocardiumy;
end

%%%%%%%%%% VLA %%%%%%%%%%

%erase endo in apex where there is no lumen and
%erase endo and epi in the outflow tract
lumendiameter = abs(VLAendox(:,:,1)-VLAendox(:,:,2));
wallthickness = abs(VLAendox-VLAepix);
lastendotemp = nan(nbroftf,1);
firstslice = nan(nbroftf,2);
lastepi = nan(nbroftf,1);
for t = isLVseg
  lastendotemp(t) = min(round(0.15*length(lumendiameter)),find(flipud(lumendiameter(:,t)) > 8,1,'first'));
  firstslice(t,1) = find(wallthickness(:,t,1) > 3/resamplexy,1,'first');
  firstslice(t,2) = find(wallthickness(:,t,2) > 3/resamplexy,1,'first');
  lastepi(t) = find(~isnan(VLAepix(:,t,1)),1,'last');
end
lastendo = m-lastendotemp+1;
%put together VLA myocardium segmentation
VLAmyocardiumx = cell(1,nbroftf);
VLAmyocardiumy = VLAmyocardiumx;
for t = isLVseg
  apexsmoothingxdist = abs(VLAepix(lastepi(t),t,2)-VLAepix(lastepi(t),t,1));
  apexsmoothingx = [VLAepix(lastepi(t),t,2)-0.25*apexsmoothingydist VLAepix(lastepi(t),t,1)+0.25*apexsmoothingydist];
  myocardiumx = [VLAendox(firstslice(t,1):lastendo(t),t,1)' VLAendox(lastendo(t):-1:firstslice(t,2),t,2)' ...
    VLAepix(firstslice(t,2):lastepi(t),t,2)' apexsmoothingx VLAepix(lastepi(t):-1:firstslice(t,1),t,1)' VLAendox(firstslice(t,1),t,1)'];
  myocardiumx = interp1(linspace(0,1,length(myocardiumx)),myocardiumx,linspace(0,1,2*length(myocardiumx)));
  for k = 1:4  %smoothing
    myocardiumx = myocardiumx-0.2*(6*myocardiumx-circshift(myocardiumx,[0 1])-circshift(myocardiumx,[0 -1]) ...
      -circshift(myocardiumx,[0 2])-circshift(myocardiumx,[0 -2])-circshift(myocardiumx,[0 3])-circshift(myocardiumx,[0 -3]));
  end
  myocardiumx = [myocardiumx myocardiumx(1)];
  VLAmyocardiumx{t} = myocardiumx;
  apexsmoothingy = [VLAepiy(lastepi(t),t,2)+1 VLAepiy(lastepi(t),t,1)+1];
  myocardiumy = [VLAendoy(firstslice(t,1):lastendo(t),t,1)' VLAendoy(lastendo(t):-1:firstslice(t,2),t,2)' ...
    VLAepiy(firstslice(t,2):lastepi(t),t,2)' apexsmoothingy VLAepiy(lastepi(t):-1:firstslice(t,1),t,1)' VLAendoy(firstslice(t,1),t,1)'];
  myocardiumy = interp1(linspace(0,1,length(myocardiumy)),myocardiumy,linspace(0,1,2*length(myocardiumy)));
  myocardiumy = [myocardiumy myocardiumy(1)];
  VLAmyocardiumy{t} = myocardiumy;
end
  
%------------------
function plotimages
%------------------
%plot the images
global DATA

gui = DATA.GUI.SpectPlot2d;

if gui.showno(1) ~= 0
  switch gui.showimagetype{1}
    case 'stress'
      tf = 1;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount = gui.maxcountstress(tf);
      else
        normmaxcount = max(gui.maxcountstress);
      end
      normmincount = gui.mincountstress;
      gui.handles.axesSAbasalimageleft = imagesc(gui.SAbasalimagestress(:,:,tf),'parent',gui.handles.axesSAbasalleft);
      gui.handles.axesSAmidimageleft= imagesc(gui.SAmidimagestress(:,:,tf),'parent',gui.handles.axesSAmidleft);
      gui.handles.axesSAapicalimageleft = imagesc(gui.SAapicalimagestress(:,:,tf),'parent',gui.handles.axesSAapicalleft);
      gui.handles.axesHLAimageleft = imagesc(gui.HLAimagestress(:,:,tf),'parent',gui.handles.axesHLAleft);
      gui.handles.axesVLAimageleft = imagesc(gui.VLAimagestress(:,:,tf),'parent',gui.handles.axesVLAleft);
      textleft = 'Stress';
    case 'stressgated'
      tf = gui.handles.tf;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount = gui.maxcountstressgated(tf);
        normmincount = gui.mincountstressgated(tf);
      else
        normmaxcount = max(gui.maxcountstressgated);
        normmincount = min(gui.mincountstressgated);
      end
      gui.handles.axesSAbasalimageleft = imagesc(gui.SAbasalimagestressgated(:,:,tf),'parent',gui.handles.axesSAbasalleft);
      gui.handles.axesSAmidimageleft= imagesc(gui.SAmidimagestressgated(:,:,tf),'parent',gui.handles.axesSAmidleft);
      gui.handles.axesSAapicalimageleft = imagesc(gui.SAapicalimagestressgated(:,:,tf),'parent',gui.handles.axesSAapicalleft);
      gui.handles.axesHLAimageleft = imagesc(gui.HLAimagestressgated(:,:,tf),'parent',gui.handles.axesHLAleft);
      gui.handles.axesVLAimageleft = imagesc(gui.VLAimagestressgated(:,:,tf),'parent',gui.handles.axesVLAleft);
      textleft = 'Stress gated';
    case 'restgated'
      tf = gui.handles.tf;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount = gui.maxcountrestgated(tf);
        normmincount = gui.mincountrestgated(tf);
      else
        normmaxcount = max(gui.maxcountrestgated);
        normmincount = min(gui.mincountrestgated);
      end
      gui.handles.axesSAbasalimageleft = imagesc(gui.SAbasalimagerestgated(:,:,tf),'parent',gui.handles.axesSAbasalleft);
      gui.handles.axesSAmidimageleft= imagesc(gui.SAmidimagerestgated(:,:,tf),'parent',gui.handles.axesSAmidleft);
      gui.handles.axesSAapicalimageleft = imagesc(gui.SAapicalimagerestgated(:,:,tf),'parent',gui.handles.axesSAapicalleft);
      gui.handles.axesHLAimageleft = imagesc(gui.HLAimagerestgated(:,:,tf),'parent',gui.handles.axesHLAleft);
      gui.handles.axesVLAimageleft = imagesc(gui.VLAimagerestgated(:,:,tf),'parent',gui.handles.axesVLAleft);
      textleft = 'Rest gated';
  end
  
  switch gui.showimagetype{1}    
    case {'stress','stressgated','restgated'}
      colormap(gui.handles.axesSAapicalleft,'spect');
      colormap(gui.handles.axesSAmidleft,'spect');
      colormap(gui.handles.axesSAbasalleft,'spect');      
      colormap(gui.handles.axesHLAleft,'spect');
      colormap(gui.handles.axesVLAleft,'spect');
      
      axis(gui.handles.axesSAbasalleft,'image','off');
      set(gui.handles.axesSAbasalleft,'clim',[normmincount normmaxcount]);
      axis(gui.handles.axesSAmidleft,'image','off');
      set(gui.handles.axesSAmidleft,'clim',[normmincount normmaxcount]);
      axis(gui.handles.axesSAapicalleft,'image','off');
      set(gui.handles.axesSAapicalleft,'clim',[normmincount normmaxcount]);
      %plot HLA
      axis(gui.handles.axesHLAleft,'image','off');
      set(gui.handles.axesHLAleft,'clim',[normmincount normmaxcount]);
      %plot VLA
      axis(gui.handles.axesVLAleft,'image','off');
      set(gui.handles.axesVLAleft,'clim',[normmincount normmaxcount]);
      set(gui.handles.textleft,'String',textleft);

    otherwise 
      %no images to plot      
      gui.handles.axesSAbasalimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAbasalleft);
      gui.handles.axesSAmidimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAmidleft);
      gui.handles.axesSAapicalimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAapicalleft);
      set(gui.handles.axesSAbasalimageleft,'cdata',[]);
      set(gui.handles.axesSAmidimageleft,'cdata',[]);
      set(gui.handles.axesSAapicalimageleft,'cdata',[]);
      axis(gui.handles.axesSAbasalleft,'off');
      axis(gui.handles.axesSAmidleft,'off');
      axis(gui.handles.axesSAapicalleft,'off');
      set(gui.handles.textSAbasalleft,'String',[]);
      set(gui.handles.textSAmidleft,'String',[]);
      set(gui.handles.textSAapicalleft,'String',[]);
      gui.handles.axesHLAimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesHLAleft);
      gui.handles.axesVLAimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesVLAleft);
      set(gui.handles.axesHLAimageleft,'cdata',[]);
      set(gui.handles.axesVLAimageleft,'cdata',[]);
      axis(gui.handles.axesHLAleft,'off');
      axis(gui.handles.axesVLAleft,'off');
      set(gui.handles.textHLAleft,'String',[]);
      set(gui.handles.textVLAleft,'String',[]);
      set(gui.handles.textleft,'String',[]);
  end
else
  %no images to plot
  gui.handles.axesSAbasalimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAbasalleft);
  gui.handles.axesSAmidimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAmidleft);
  gui.handles.axesSAapicalimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAapicalleft);
  set(gui.handles.axesSAbasalimageleft,'cdata',[]);
  set(gui.handles.axesSAmidimageleft,'cdata',[]);
  set(gui.handles.axesSAapicalimageleft,'cdata',[]);
  axis(gui.handles.axesSAbasalleft,'off');
  axis(gui.handles.axesSAmidleft,'off');
  axis(gui.handles.axesSAapicalleft,'off');
  set(gui.handles.textSAbasalleft,'String',[]);
  set(gui.handles.textSAmidleft,'String',[]);
  set(gui.handles.textSAapicalleft,'String',[]);
  gui.handles.axesHLAimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesHLAleft);
  gui.handles.axesVLAimageleft = imagesc([1 1 ; 1 1],'parent',gui.handles.axesVLAleft);
  set(gui.handles.axesHLAimageleft,'cdata',[]);
  set(gui.handles.axesVLAimageleft,'cdata',[]);
  axis(gui.handles.axesHLAleft,'off');
  axis(gui.handles.axesVLAleft,'off');
  set(gui.handles.textHLAleft,'String',[]);
  set(gui.handles.textVLAleft,'String',[]);
  set(gui.handles.textleft,'String',[]);
  set(gui.handles.uipanelcorrectLVleft,'visible','off');
end

if gui.showno(2) ~= 0
  switch gui.showimagetype{2}
    case 'rest'
      tf = 1;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount =  gui.maxcountrest(tf);
      else
        normmaxcount = max(gui.maxcountrest);
      end
      normmincount = gui.mincountrest;
      gui.handles.axesSAbasalimageright = imagesc(gui.SAbasalimagerest(:,:,tf),'parent',gui.handles.axesSAbasalright);
      gui.handles.axesSAmidimageright = imagesc(gui.SAmidimagerest(:,:,tf),'parent',gui.handles.axesSAmidright);
      gui.handles.axesSAapicalimageright = imagesc(gui.SAapicalimagerest(:,:,tf),'parent',gui.handles.axesSAapicalright);
      gui.handles.axesHLAimageright = imagesc(gui.HLAimagerest(:,:,tf),'parent',gui.handles.axesHLAright);
      gui.handles.axesVLAimageright = imagesc(gui.VLAimagerest(:,:,tf),'parent',gui.handles.axesVLAright);
      textright = 'Rest';
    case 'restgated'
      tf = gui.handles.tf;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount =  gui.maxcountrestgated(tf);
      else
        normmaxcount = max(gui.maxcountrestgated);
      end
      normmincount = gui.mincountrestgated(tf);
      gui.handles.axesSAbasalimageright = imagesc(gui.SAbasalimagerestgated(:,:,tf),'parent',gui.handles.axesSAbasalright);
      gui.handles.axesSAmidimageright = imagesc(gui.SAmidimagerestgated(:,:,tf),'parent',gui.handles.axesSAmidright);
      gui.handles.axesSAapicalimageright = imagesc(gui.SAapicalimagerestgated(:,:,tf),'parent',gui.handles.axesSAapicalright);
      gui.handles.axesHLAimageright = imagesc(gui.HLAimagerestgated(:,:,tf),'parent',gui.handles.axesHLAright);
      gui.handles.axesVLAimageright = imagesc(gui.VLAimagerestgated(:,:,tf),'parent',gui.handles.axesVLAright);
      textright = 'Rest gated';
  end  
  switch gui.showimagetype{2}
    case {'rest','restgated'}
      colormap(gui.handles.axesSAapicalright,'spect');
      colormap(gui.handles.axesSAmidright,'spect');
      colormap(gui.handles.axesSAbasalright,'spect');      
      colormap(gui.handles.axesHLAright,'spect');
      colormap(gui.handles.axesVLAright,'spect');
      axis(gui.handles.axesSAbasalright,'image','off');
      set(gui.handles.axesSAbasalright,'clim',[normmincount normmaxcount]);
      axis(gui.handles.axesSAmidright,'image','off');
      set(gui.handles.axesSAmidright,'clim',[normmincount normmaxcount]);
      axis(gui.handles.axesSAapicalright,'image','off');
      set(gui.handles.axesSAapicalright,'clim',[normmincount normmaxcount]);      
      %plot HLA
      axis(gui.handles.axesHLAright,'image','off');
      set(gui.handles.axesHLAright,'clim',[normmincount normmaxcount]);
      %VLA
      axis(gui.handles.axesVLAright,'image','off');
      set(gui.handles.axesVLAright,'clim',[normmincount normmaxcount]);
      set(gui.handles.textright,'String',textright);
      
    otherwise
      %no images to plot
      gui.handles.axesSAbasalimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAbasalright);
      gui.handles.axesSAmidimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAmidright);
      gui.handles.axesSAapicalimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAapicalright);
      set(gui.handles.axesSAbasalimageright,'cdata',[]);
      set(gui.handles.axesSAmidimageright,'cdata',[]);
      set(gui.handles.axesSAapicalimageright,'cdata',[]);
      axis(gui.handles.axesSAbasalright,'off');
      axis(gui.handles.axesSAmidright,'off');
      axis(gui.handles.axesSAapicalright,'off');
      set(gui.handles.textSAbasalright,'String',[]);
      set(gui.handles.textSAmidright,'String',[]);
      set(gui.handles.textSAapicalright,'String',[]);
      gui.handles.axesHLAimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesHLAright);
      gui.handles.axesVLAimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesVLAright);
      set(gui.handles.axesHLAimageright,'cdata',[]);
      set(gui.handles.axesVLAimageright,'cdata',[]);
      axis(gui.handles.axesHLAright,'off');
      axis(gui.handles.axesVLAright,'off');
      set(gui.handles.textHLAright,'String',[]);
      set(gui.handles.textVLAright,'String',[]);
      set(gui.handles.textright,'String',[]);
  end
else
  %no images to plot
  gui.handles.axesSAbasalimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAbasalright);
  gui.handles.axesSAmidimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAmidright);
  gui.handles.axesSAapicalimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesSAapicalright);
  set(gui.handles.axesSAbasalimageright,'cdata',[]);
  set(gui.handles.axesSAmidimageright,'cdata',[]);
  set(gui.handles.axesSAapicalimageright,'cdata',[]);
  axis(gui.handles.axesSAbasalright,'off');
  axis(gui.handles.axesSAmidright,'off');
  axis(gui.handles.axesSAapicalright,'off');
  set(gui.handles.textSAbasalright,'String',[]);
  set(gui.handles.textSAmidright,'String',[]);
  set(gui.handles.textSAapicalright,'String',[]);
  gui.handles.axesHLAimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesHLAright);
  gui.handles.axesVLAimageright = imagesc([1 1 ; 1 1],'parent',gui.handles.axesVLAright);
  set(gui.handles.axesHLAimageright,'cdata',[]);
  set(gui.handles.axesVLAimageright,'cdata',[]);
  axis(gui.handles.axesHLAright,'off');
  axis(gui.handles.axesVLAright,'off');
  set(gui.handles.textHLAright,'String',[]);
  set(gui.handles.textVLAright,'String',[]);
  set(gui.handles.textright,'String',[]);
  set(gui.handles.uipanelcorrectLVright,'visible','off');
end


%---------------------------
function plotLVsegmentation
%---------------------------
%draw the endo- and epicardial segmentation

global DATA

gui = DATA.GUI.SpectPlot2d;

if gui.showno(1) ~= 0
  LVplotting = 0;
  if get(gui.handles.radiobuttonshowsegmentation,'value')
    switch gui.showimagetype{1}
      case 'stress'
        if gui.LVsegmentationstress
          endoX = gui.EndoXstress;
          endoY = gui.EndoYstress;
          epiX = gui.EpiXstress;
          epiY = gui.EpiYstress;
          HLAmyocardiumx = gui.HLAmyocardiumxstress;
          HLAmyocardiumy = gui.HLAmyocardiumystress;
          VLAmyocardiumx = gui.VLAmyocardiumxstress;
          VLAmyocardiumy = gui.VLAmyocardiumystress;
          SAslices = gui.SAslicesstress;
          LVplotting = 1;
          tf = 1;
        end
      case 'stressgated'
        if gui.LVsegmentationstressgated
          endoX = gui.EndoXstressgated;
          endoY = gui.EndoYstressgated;
          epiX = gui.EpiXstressgated;
          epiY = gui.EpiYstressgated;
          HLAmyocardiumx = gui.HLAmyocardiumxstressgated;
          HLAmyocardiumy = gui.HLAmyocardiumystressgated;
          VLAmyocardiumx = gui.VLAmyocardiumxstressgated;
          VLAmyocardiumy = gui.VLAmyocardiumystressgated;
          SAslices = gui.SAslicesstressgated;
          LVplotting = 1;
          tf = gui.handles.tf;
        end
      case 'restgated'
        if gui.LVsegmentationrestgated
          endoX = gui.EndoXrestgated;
          endoY = gui.EndoYrestgated;
          epiX = gui.EpiXrestgated;
          epiY = gui.EpiYrestgated;
          HLAmyocardiumx = gui.HLAmyocardiumxrestgated;
          HLAmyocardiumy = gui.HLAmyocardiumyrestgated;
          VLAmyocardiumx = gui.VLAmyocardiumxrestgated;
          VLAmyocardiumy = gui.VLAmyocardiumyrestgated;
          SAslices = gui.SAslicesrestgated;
          LVplotting = 1;
          tf = gui.handles.tf;
        end
    end
  end
  if LVplotting
    basalendox = endoX(:,tf,round(mean(SAslices{3})));
    basalendoy = endoY(:,tf,round(mean(SAslices{3})));
    basalepix = epiX(:,tf,round(mean(SAslices{3})));
    basalepiy = epiY(:,tf,round(mean(SAslices{3})));    
    midendox = endoX(:,tf,round(mean(SAslices{2})));
    midendoy = endoY(:,tf,round(mean(SAslices{2})));
    midepix = epiX(:,tf,round(mean(SAslices{2})));
    midepiy = epiY(:,tf,round(mean(SAslices{2})));    
    apicalendox = endoX(:,tf,round(mean(SAslices{1})));
    apicalendoy = endoY(:,tf,round(mean(SAslices{1})));
    apicalepix = epiX(:,tf,round(mean(SAslices{1})));
    apicalepiy = epiY(:,tf,round(mean(SAslices{1})));
    HLAx = HLAmyocardiumx{tf};
    HLAy = HLAmyocardiumy{tf};
    VLAx = VLAmyocardiumx{tf};
    VLAy = VLAmyocardiumy{tf};
  else
    basalendox = [];
    basalendoy = [];
    basalepix = [];
    basalepiy = [];    
    midendox = [];
    midendoy = [];
    midepix = [];
    midepiy = [];    
    apicalendox = [];
    apicalendoy = [];
    apicalepix = [];
    apicalepiy = [];
    HLAx = [];
    HLAy = [];
    VLAx = [];
    VLAy = [];
  end
else
  basalendox = [];
  basalendoy = [];
  basalepix = [];
  basalepiy = [];
  midendox = [];
  midendoy = [];
  midepix = [];
  midepiy = [];
  apicalendox = [];
  apicalendoy = [];
  apicalepix = [];
  apicalepiy = [];
  HLAx = [];
  HLAy = [];
  VLAx = [];
  VLAy = [];
end
%basal
hold(gui.handles.axesSAbasalleft,'on');
gui.handles.axesSAbasalendoleft = plot(gui.handles.axesSAbasalleft,basalendoy,basalendox,'w-');
set(gui.handles.axesSAbasalendoleft,'linewidth',2);
gui.handles.axesSAbasalepileft = plot(gui.handles.axesSAbasalleft,basalepiy,basalepix,'w-');
set(gui.handles.axesSAbasalepileft,'linewidth',2);
hold(gui.handles.axesSAbasalleft,'off');
%mid
hold(gui.handles.axesSAmidleft,'on');
gui.handles.axesSAmidendoleft = plot(gui.handles.axesSAmidleft,midendoy,midendox,'w-');
set(gui.handles.axesSAmidendoleft,'linewidth',2);
gui.handles.axesSAmidepileft = plot(gui.handles.axesSAmidleft,midepiy,midepix,'w-');
set(gui.handles.axesSAmidepileft,'linewidth',2);
hold(gui.handles.axesSAmidleft,'off');
%apical
hold(gui.handles.axesSAapicalleft,'on');
gui.handles.axesSAapicalendoleft = plot(gui.handles.axesSAapicalleft,apicalendoy,apicalendox,'w-');
set(gui.handles.axesSAapicalendoleft,'linewidth',2);
gui.handles.axesSAapicalepileft = plot(gui.handles.axesSAapicalleft,apicalepiy,apicalepix,'w-');
set(gui.handles.axesSAapicalepileft,'linewidth',2);
hold(gui.handles.axesSAapicalleft,'off');
%HLA
hold(gui.handles.axesHLAleft,'on');
gui.handles.axesHLAmyocardiumleft = plot(gui.handles.axesHLAleft,HLAy,HLAx,'w-');
set(gui.handles.axesHLAmyocardiumleft,'linewidth',2);
hold(gui.handles.axesHLAleft,'off');
%VLA
hold(gui.handles.axesVLAleft,'on');
gui.handles.axesVLAmyocardiumleft = plot(gui.handles.axesVLAleft,VLAy,VLAx,'w-');
set(gui.handles.axesVLAmyocardiumleft,'linewidth',2);
hold(gui.handles.axesVLAleft,'off');

if gui.showno(2) ~= 0
  LVplotting = 0;
  if get(gui.handles.radiobuttonshowsegmentation,'value')
    switch gui.showimagetype{2}
      case 'rest'
        if gui.LVsegmentationrest
          endoX = gui.EndoXrest;
          endoY = gui.EndoYrest;
          epiX = gui.EpiXrest;
          epiY = gui.EpiYrest;
          HLAmyocardiumx = gui.HLAmyocardiumxrest;
          HLAmyocardiumy = gui.HLAmyocardiumyrest;
          VLAmyocardiumx = gui.VLAmyocardiumxrest;
          VLAmyocardiumy = gui.VLAmyocardiumyrest;
          SAslices = gui.SAslicesrest;
          LVplotting = 1;
          tf = 1;
        end
      case 'restgated'
        if gui.LVsegmentationrestgated
          endoX = gui.EndoXrestgated;
          endoY = gui.EndoYrestgated;
          epiX = gui.EpiXrestgated;
          epiY = gui.EpiYrestgated;
          HLAmyocardiumx = gui.HLAmyocardiumxrestgated;
          HLAmyocardiumy = gui.HLAmyocardiumyrestgated;
          VLAmyocardiumx = gui.VLAmyocardiumxrestgated;
          VLAmyocardiumy = gui.VLAmyocardiumyrestgated;
          SAslices = gui.SAslicesrestgated;
          LVplotting = 1;
          tf = gui.handles.tf;
        end
    end
  end
  if LVplotting
    basalendox = endoX(:,tf,round(mean(SAslices{3})));
    basalendoy = endoY(:,tf,round(mean(SAslices{3})));
    basalepix = epiX(:,tf,round(mean(SAslices{3})));
    basalepiy = epiY(:,tf,round(mean(SAslices{3})));    
    midendox = endoX(:,tf,round(mean(SAslices{2})));
    midendoy = endoY(:,tf,round(mean(SAslices{2})));
    midepix = epiX(:,tf,round(mean(SAslices{2})));
    midepiy = epiY(:,tf,round(mean(SAslices{2})));    
    apicalendox = endoX(:,tf,round(mean(SAslices{1})));
    apicalendoy = endoY(:,tf,round(mean(SAslices{1})));
    apicalepix = epiX(:,tf,round(mean(SAslices{1})));
    apicalepiy = epiY(:,tf,round(mean(SAslices{1})));
    HLAx = HLAmyocardiumx{tf};
    HLAy = HLAmyocardiumy{tf};
    VLAx = VLAmyocardiumx{tf};
    VLAy = VLAmyocardiumy{tf};
  else
    basalendox = [];
    basalendoy = [];
    basalepix = [];
    basalepiy = [];    
    midendox = [];
    midendoy = [];
    midepix = [];
    midepiy = [];    
    apicalendox = [];
    apicalendoy = [];
    apicalepix = [];
    apicalepiy = [];
    HLAx = [];
    HLAy = [];
    VLAx = [];
    VLAy = [];
  end
else
  basalendox = [];
  basalendoy = [];
  basalepix = [];
  basalepiy = [];
  midendox = [];
  midendoy = [];
  midepix = [];
  midepiy = [];
  apicalendox = [];
  apicalendoy = [];
  apicalepix = [];
  apicalepiy = [];
  HLAx = [];
  HLAy = [];
  VLAx = [];
  VLAy = [];
end
%basal
hold(gui.handles.axesSAbasalright,'on');
gui.handles.axesSAbasalendoright = plot(gui.handles.axesSAbasalright,basalendoy,basalendox,'w-');
set(gui.handles.axesSAbasalendoright,'linewidth',2);
gui.handles.axesSAbasalepiright = plot(gui.handles.axesSAbasalright,basalepiy,basalepix,'w-');
set(gui.handles.axesSAbasalepiright,'linewidth',2);
hold(gui.handles.axesSAbasalright,'off');
%mid
hold(gui.handles.axesSAmidright,'on');
gui.handles.axesSAmidendoright = plot(gui.handles.axesSAmidright,midendoy,midendox,'w-');
set(gui.handles.axesSAmidendoright,'linewidth',2);
gui.handles.axesSAmidepiright = plot(gui.handles.axesSAmidright,midepiy,midepix,'w-');
set(gui.handles.axesSAmidepiright,'linewidth',2);
hold(gui.handles.axesSAmidright,'off');
%apical
hold(gui.handles.axesSAapicalright,'on');
gui.handles.axesSAapicalendoright = plot(gui.handles.axesSAapicalright,apicalendoy,apicalendox,'w-');
set(gui.handles.axesSAapicalendoright,'linewidth',2);
gui.handles.axesSAapicalepiright = plot(gui.handles.axesSAapicalright,apicalepiy,apicalepix,'w-');
set(gui.handles.axesSAapicalepiright,'linewidth',2);
hold(gui.handles.axesSAapicalright,'off');
%HLA
hold(gui.handles.axesHLAright,'on');
gui.handles.axesHLAmyocardiumright = plot(gui.handles.axesHLAright,HLAy,HLAx,'w-');
set(gui.handles.axesHLAmyocardiumright,'linewidth',2);
hold(gui.handles.axesHLAright,'off');
%VLA
hold(gui.handles.axesVLAright,'on');
gui.handles.axesVLAmyocardiumright = plot(gui.handles.axesVLAright,VLAy,VLAx,'w-');
set(gui.handles.axesVLAmyocardiumright,'linewidth',2);
hold(gui.handles.axesVLAright,'off');


%------------------------------
function plotahasections(panel)
%------------------------------
%plot the 17 aha section model

global DATA

gui = DATA.GUI.SpectPlot2d;

if isequal(panel,'left') || isequal(panel,'both')
  if gui.showno(1) ~= 0
    switch gui.showimagetype{1}
      case 'stress'
        if gui.LVsegmentationstress
          LVcenterxapical = mean(get(gui.handles.axesSAapicalendoleft,'xdata'));
          LVcenteryapical = mean(get(gui.handles.axesSAapicalendoleft,'ydata'));
          LVcenterxmid = mean(get(gui.handles.axesSAmidendoleft,'xdata'));
          LVcenterymid = mean(get(gui.handles.axesSAmidendoleft,'ydata'));
          LVcenterxbasal = mean(get(gui.handles.axesSAbasalendoleft,'xdata'));
          LVcenterybasal = mean(get(gui.handles.axesSAbasalendoleft,'ydata'));
          b = sqrt(0.75);
          a = 0.5;
          c = 1/sqrt(2);
          %apical
          ax = gui.handles.axesSAapicalleft;
          hold(ax,'on');
          gui.handles.axesSAapicalahalineleft = plot(ax,LVcenterxapical*[1-c 1+c],LVcenteryapical*[1-c 1+c],'w', ...
            LVcenterxapical*[1-c 1+c],LVcenteryapical*[1+c 1-c],'w-');
          hold(ax,'off');
          %mid
          ax = gui.handles.axesSAmidleft;
          hold(ax,'on');
          gui.handles.axesSAmidahalineleft = plot(ax,LVcenterxmid*[0 2],LVcenterymid*[1 1] ,'w', ...
            LVcenterxmid*[1-1*a 1+1*a],LVcenterymid*[1-1*b 1+1*b],'w', ...
            LVcenterxmid*[1-1*a 1+1*a],LVcenterymid*[1+1*b 1-1*b],'w');
          hold(ax,'off');
          %basal
          ax = gui.handles.axesSAbasalleft;
          hold(ax,'on');
          gui.handles.axesSAbasalahalineleft = plot(ax,LVcenterxbasal*[0 2],LVcenterybasal*[1 1],'w', ...
            LVcenterxbasal*[1-1*a 1+1*a],LVcenterybasal*[1-1*b 1+1*b],'w', ...
            LVcenterxbasal*[1-1*a 1+1*a],LVcenterybasal*[1+1*b 1-1*b],'w-');
          hold(ax,'off');
          %HLA
          ax = gui.handles.axesHLAleft;
          xline = get(gui.handles.axesHLAleft,'xlim');
          hold(ax,'on');
          ylim = get(gui.handles.axesHLAleft,'ylim');
          yline = [(ylim(2)-gui.SAslicesstress{1}(end))*[1 1] ; (ylim(2)-gui.SAslicesstress{2}(end))*[1 1] ; ...
            (ylim(2)-gui.SAslicesstress{3}(end))*[1 1] ; (ylim(2)-gui.SAslicesstress{3}(1))*[1 1]];
          gui.handles.axesHLAahalineleft = plot(ax,xline,yline(1,:),'w',xline,yline(2,:),'w',xline,yline(3,:),'w',xline,yline(4,:),'w-');
          hold(ax,'off');
          %VLA
          ax = gui.handles.axesVLAleft;
          yline = get(gui.handles.axesVLAleft,'xlim');
          hold(ax,'on');
          xline = [gui.SAslicesstress{1}(end)*[1 1] ; gui.SAslicesstress{2}(end)*[1 1] ; ...
            gui.SAslicesstress{3}(end)*[1 1] ; gui.SAslicesstress{3}(1)*[1 1]];
          gui.handles.axesVLAahalineleft = plot(ax,xline(1,:),yline,'w',xline(2,:),yline,'w',xline(3,:),yline ,'w',xline(4,:),yline,'w-');
          hold(ax,'off');
        end
    end
  end
end
if isequal(panel,'right') || isequal(panel,'both')
  if gui.showno(2) ~= 0
    switch gui.showimagetype{2}
      case 'rest'
        if gui.LVsegmentationrest
          LVcenterxapical = mean(get(gui.handles.axesSAapicalendoright,'xdata'));
          LVcenteryapical = mean(get(gui.handles.axesSAapicalendoright,'ydata'));
          LVcenterxmid = mean(get(gui.handles.axesSAmidendoright,'xdata'));
          LVcenterymid = mean(get(gui.handles.axesSAmidendoright,'ydata'));
          LVcenterxbasal = mean(get(gui.handles.axesSAbasalendoright,'xdata'));
          LVcenterybasal = mean(get(gui.handles.axesSAbasalendoright,'ydata'));
          b = sqrt(0.75);
          a = 0.5;
          c = 1/sqrt(2);
          %apical
          ax = gui.handles.axesSAapicalright;
          hold(ax,'on');
          gui.handles.axesSAapicalahalineright = plot(ax,LVcenterxapical*[1-c 1+c],LVcenteryapical*[1-c 1+c],'w', ...
            LVcenterxapical*[1-c 1+c],LVcenteryapical*[1+c 1-c],'w-');
          hold(ax,'off');
          %mid
          ax = gui.handles.axesSAmidright;
          hold(ax,'on');
          gui.handles.axesSAmidahalineright = plot(ax,LVcenterxmid*[0 2],LVcenterymid*[1 1],'w', ...
            LVcenterxmid*[1-1*a 1+1*a],LVcenterymid*[1-1*b 1+1*b],'w', ...
            LVcenterxmid*[1-1*a 1+1*a],LVcenterymid*[1+1*b 1-1*b],'w-');
          hold(ax,'off');
          %basal
          ax = gui.handles.axesSAbasalright;
          hold(ax,'on');
          gui.handles.axesSAbasalahalineright = plot(ax,LVcenterxbasal*[0 2],LVcenterybasal*[1 1],'w', ...
            LVcenterxbasal*[1-1*a 1+1*a],LVcenterybasal*[1-1*b 1+1*b],'w', ...
            LVcenterxbasal*[1-1*a 1+1*a],LVcenterybasal*[1+1*b 1-1*b],'w-');
          hold(ax,'off');
          %HLA
          ax = gui.handles.axesHLAright;
          xline = get(gui.handles.axesHLAright,'xlim');
          hold(ax,'on');
          ylim = get(gui.handles.axesHLAright,'ylim');
          yline = [(ylim(2)-gui.SAslicesrest{1}(end))*[1 1] ; (ylim(2)-gui.SAslicesrest{2}(end))*[1 1] ; ...
            (ylim(2)-gui.SAslicesrest{3}(end))*[1 1] ; (ylim(2)-gui.SAslicesrest{3}(1))*[1 1]];
          gui.handles.axesHLAahalineright = plot(ax,xline,yline(1,:),'w',xline,yline(2,:),'w',xline,yline(3,:),'w',xline,yline(4,:),'w-');
          hold(ax,'off');
          %VLA
          ax = gui.handles.axesVLAright;
          yline = get(gui.handles.axesVLAright,'xlim');
          hold(ax,'on');
          xline = [gui.SAslicesrest{1}(end)*[1 1] ; gui.SAslicesrest{2}(end)*[1 1] ; ...
            gui.SAslicesrest{3}(end)*[1 1] ; gui.SAslicesrest{3}(1)*[1 1]];
          gui.handles.axesVLAahalineright = plot(ax,xline(1,:),yline,'w',xline(2,:),yline,'w',xline(3,:),yline,'w',xline(4,:),yline,'w-');
          hold(ax,'off');
        end
    end
  end
end
if isequal(panel,'right') || isequal(panel,'none')
  %left
  ax = gui.handles.axesSAapicalleft;
  hold(ax,'on');
  gui.handles.axesSAapicalahalineleft = plot(ax,0,0,'w',0,0,'w');
  hold(ax,'off');
  ax = gui.handles.axesSAmidleft;
  hold(ax,'on');
  gui.handles.axesSAmidahalineleft = plot(ax,0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
  ax = gui.handles.axesSAbasalleft;
  hold(ax,'on');
  gui.handles.axesSAbasalahalineleft = plot(ax,0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
  ax = gui.handles.axesHLAleft;
  hold(ax,'on');
  gui.handles.axesHLAahalineleft = plot(ax,0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
  ax = gui.handles.axesVLAleft;
  hold(ax,'on');
  gui.handles.axesVLAahalineleft = plot(ax,0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
end
if isequal(panel,'left') || isequal(panel,'none')
  %right
  ax = gui.handles.axesSAapicalright;
  hold(ax,'on');
  gui.handles.axesSAapicalahalineright = plot(ax,0,0,'w',0,0,'w-');
  hold(ax,'off');
  ax = gui.handles.axesSAmidright;
  hold(ax,'on');
  gui.handles.axesSAmidahalineright = plot(ax,0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
  ax = gui.handles.axesSAbasalright;
  hold(ax,'on');
  gui.handles.axesSAbasalahalineright = plot(ax,0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
  ax = gui.handles.axesHLAright;
  hold(ax,'on');
  gui.handles.axesHLAahalineright = plot(ax,0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
  ax = gui.handles.axesVLAright;
  hold(ax,'on');
  gui.handles.axesVLAahalineright = plot(ax,0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  hold(ax,'off');
end


%----------------------------
function plotahavalues(panel)
%----------------------------
%show the 17 aha section model values

global DATA

gui = DATA.GUI.SpectPlot2d;
SSS = 0;
SRS = 0;
SDS = 0;
pos = [17 13 16:-1:14 7 12:-1:8 1 6:-1:2];
gating = get(gui.handles.radiobuttongated,'value');

switch panel
  case 'both'
    editlabels = {gui.handles.editHLAleft, ...
      gui.handles.editapicalantleft,gui.handles.editapicallatleft, ...
      gui.handles.editapicalinfleft,gui.handles.editapicalsepleft, ...
      gui.handles.editmidantleft,gui.handles.editmidantlatleft,gui.handles.editmidinflatleft, ...
      gui.handles.editmidinfleft,gui.handles.editmidinfsepleft,gui.handles.editmidantsepleft, ...
      gui.handles.editbasalantleft,gui.handles.editbasalantlatleft,gui.handles.editbasalinflatleft, ...
      gui.handles.editbasalinfleft,gui.handles.editbasalinfsepleft,gui.handles.editbasalantsepleft, ...      
      gui.handles.editHLAright, ...
      gui.handles.editapicalantright,gui.handles.editapicallatright, ...
      gui.handles.editapicalinfright,gui.handles.editapicalsepright, ...
      gui.handles.editmidantright,gui.handles.editmidantlatright,gui.handles.editmidinflatright, ...
      gui.handles.editmidinfright,gui.handles.editmidinfsepright,gui.handles.editmidantsepright, ...
      gui.handles.editbasalantright,gui.handles.editbasalantlatright,gui.handles.editbasalinflatright, ...
      gui.handles.editbasalinfright,gui.handles.editbasalinfsepright,gui.handles.editbasalantsepright};
    pos = [pos pos];
    for editloop = 1:length(editlabels)
      set(editlabels{editloop},'visible','on');
      if editloop <= 17
        set(editlabels{editloop},'String',gui.scoringvalues{1,pos(editloop)});
      else
        set(editlabels{editloop},'String',gui.scoringvalues{2,pos(editloop)});
      end
    end
  case 'right'
    editlabels = {gui.handles.editHLAright, ...
      gui.handles.editapicalantright,gui.handles.editapicallatright, ...
      gui.handles.editapicalinfright,gui.handles.editapicalsepright, ...
      gui.handles.editmidantright,gui.handles.editmidantlatright,gui.handles.editmidinflatright, ...
      gui.handles.editmidinfright,gui.handles.editmidinfsepright,gui.handles.editmidantsepright, ...
      gui.handles.editbasalantright,gui.handles.editbasalantlatright,gui.handles.editbasalinflatright, ...
      gui.handles.editbasalinfright,gui.handles.editbasalinfsepright,gui.handles.editbasalantsepright};
    for editloop = 1:length(editlabels)
      set(editlabels{editloop},'visible','on');
      if gating == 1
        set(editlabels{editloop},'String',gui.scoringvalues{3,pos(editloop)});
      else
        set(editlabels{editloop},'String',gui.scoringvalues{2,pos(editloop)});
      end
    end
  case 'left'
    editlabels = {gui.handles.editHLAleft, ...
      gui.handles.editapicalantleft,gui.handles.editapicallatleft, ...
      gui.handles.editapicalinfleft,gui.handles.editapicalsepleft, ...
      gui.handles.editmidantleft,gui.handles.editmidantlatleft,gui.handles.editmidinflatleft, ...
      gui.handles.editmidinfleft,gui.handles.editmidinfsepleft,gui.handles.editmidantsepleft, ...
      gui.handles.editbasalantleft,gui.handles.editbasalantlatleft,gui.handles.editbasalinflatleft, ...
      gui.handles.editbasalinfleft,gui.handles.editbasalinfsepleft,gui.handles.editbasalantsepleft};
    pos = [pos pos];
    for editloop = 1:length(editlabels)
      set(editlabels{editloop},'visible','on');
      set(editlabels{editloop},'String',gui.scoringvalues{1,pos(editloop)});
    end
end

switch gating
  case 1
    for loop = 1:size(gui.scoringvalues,2)
      SRS = SRS+gui.scoringvalues{3,loop};
    end
    set(gui.handles.textScoring,'String',['SRS: ',num2str(SRS)]);
    scoringexplanation =  sprintf('%s\n%s\n%s\n%s\n%s\n%s','Scoring of presence of infarct','0: normal','1: equivocal','2: moderate','3: severe infact','4: apparent infact');
    set(gui.handles.textScoringExplanation,'String',scoringexplanation);
  case 0
    for loop = 1:size(gui.scoringvalues,2)
      SSS = SSS+gui.scoringvalues{1,loop};
      SRS = SRS+gui.scoringvalues{2,loop};
      SDS = SDS+max(0,gui.scoringvalues{1,loop}-gui.scoringvalues{2,loop});
    end
    scorevalues = sprintf('%s\n%s\n%s',['SSS: ',num2str(SSS)],['SRS: ',num2str(SRS)],['SDS: ',num2str(SDS)]);
    set(gui.handles.textScoring,'String',scorevalues);
    scoringexplanation =  sprintf('%s\n%s\n%s\n%s\n%s\n%s','Scoring of tracer uptake','0: normal','1: equivocal','2: moderate','3: severe reduction','4: apparent absence');
    set(gui.handles.textScoringExplanation,'String',scoringexplanation);
end
    


%----------------------------------
function plotSAintersections(panel)
%----------------------------------

global DATA

gui = DATA.GUI.SpectPlot2d;

if isequal(panel,'left') || isequal(panel,'both')
  if gui.showno(1) ~= 0
    switch gui.showimagetype{1}
      case 'stress'
        tf = 1;
        plotLAbasalline = mean(gui.SAslicesstress{3});
        plotLAmidline = mean(gui.SAslicesstress{2});
        plotLAapicalline = mean(gui.SAslicesstress{1});
      case 'stressgated'
        tf = gui.handles.tf;
        plotLAbasalline = mean(gui.SAslicesstressgated{3,tf});
        plotLAmidline = mean(gui.SAslicesstressgated{2,tf});
        plotLAapicalline = mean(gui.SAslicesstressgated{1,tf});
      case 'restgated'
        tf = gui.handles.tf;
        plotLAbasalline = mean(gui.SAslicesrestgated{3,tf});
        plotLAmidline = mean(gui.SAslicesrestgated{2,tf});
        plotLAapicalline = mean(gui.SAslicesrestgated{1,tf});
    end
    %HLA
    hold(gui.handles.axesHLAleft,'on');
    HLAxlim = get(gui.handles.axesHLAleft,'xlim');
    ylim = get(gui.handles.axesHLAleft,'ylim');
    gui.handles.axesHLAlinesleft = plot(gui.handles.axesHLAleft, ...
      [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))], ...
      (ylim(2)-plotLAbasalline)*[1 1],'w-', ...
      [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)], ...
      (ylim(2)-plotLAbasalline)*[1 1],'w-', ...
      [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))], ...
      (ylim(2)-plotLAmidline)*[1 1],'w-', ...
      [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)], ...
      (ylim(2)-plotLAmidline)*[1 1],'w-', ...
      [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))], ...
      (ylim(2)-plotLAapicalline)*[1 1],'w-', ...
      [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)], ...
      (ylim(2)-plotLAapicalline)*[1 1],'w-');
    set(gui.handles.axesHLAlinesleft,'linewidth',2);
    hold(gui.handles.axesHLAleft,'off');
    %VLA
    hold(gui.handles.axesVLAleft,'on');
    VLAylim = get(gui.handles.axesVLAleft,'ylim');
    gui.handles.axesVLAlinesleft = plot(gui.handles.axesVLAleft, ...
      [plotLAbasalline plotLAbasalline], ...
      [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))],'w-', ...
      [plotLAbasalline plotLAbasalline], ...
      [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)],'w-', ...
      [plotLAmidline plotLAmidline], ...
      [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))],'w-', ...
      [plotLAmidline plotLAmidline], ...
      [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)],'w-', ...
      [plotLAapicalline plotLAapicalline], ...
      [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))],'w-', ...
      [plotLAapicalline plotLAapicalline], ...
      [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)],'w-');
    set(gui.handles.axesVLAlinesleft,'linewidth',2);
    hold(gui.handles.axesVLAleft,'off');
  else
    %HLA
    hold(gui.handles.axesHLAleft,'on');
    gui.handles.axesHLAlinesleft = plot(gui.handles.axesHLAleft, ...
      [],[],'w-',[],[],'w',[],[],'w',[],[],'w',[],[],'w',[],[],'w');
    set(gui.handles.axesHLAlinesleft,'linewidth',2);
    hold(gui.handles.axesHLAleft,'off');
    %VLA
    hold(gui.handles.axesVLAleft,'on');
    gui.handles.axesVLAlinesleft = plot(gui.handles.axesVLAleft, ...
      [],[],'w-',[],[],'w',[],[],'w',[],[],'w',[],[],'w',[],[],'w');
    set(gui.handles.axesVLAlinesleft,'linewidth',2);
    hold(gui.handles.axesVLAleft,'off');
  end
end
if isequal(panel,'right') || isequal(panel,'both')
  if gui.showno(2) ~= 0
    switch gui.showimagetype{2}
      case 'rest'
        tf = 1;
        plotLAbasalline = mean(gui.SAslicesrest{3});
        plotLAmidline = mean(gui.SAslicesrest{2});
        plotLAapicalline = mean(gui.SAslicesrest{1});
      case 'restgated'
        tf = gui.handles.tf;
        plotLAbasalline = mean(gui.SAslicesrestgated{3,tf});
        plotLAmidline = mean(gui.SAslicesrestgated{2,tf});
        plotLAapicalline = mean(gui.SAslicesrestgated{1,tf});
    end
    %HLA
    hold(gui.handles.axesHLAright,'on');
    HLAxlim = get(gui.handles.axesHLAright,'xlim');
    ylim = get(gui.handles.axesHLAright,'ylim');
    gui.handles.axesHLAlinesright = plot(gui.handles.axesHLAright, ...
      [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))], ...
      (ylim(2)-plotLAbasalline)*[1 1],'w-', ...
      [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)], ...
      (ylim(2)-plotLAbasalline)*[1 1],'w-', ...
      [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))], ...
      (ylim(2)-plotLAmidline)*[1 1],'w-', ...
      [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)], ...
      (ylim(2)-plotLAmidline)*[1 1],'w-', ...
      [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))], ...
      (ylim(2)-plotLAapicalline)*[1 1],'w-', ...
      [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)], ...
      (ylim(2)-plotLAapicalline)*[1 1],'w-');
    set(gui.handles.axesHLAlinesright,'linewidth',2);
    hold(gui.handles.axesHLAright,'off');
    %VLA
    hold(gui.handles.axesVLAright,'on');
    VLAylim = get(gui.handles.axesVLAright,'ylim');
    gui.handles.axesVLAlinesright = plot(gui.handles.axesVLAright, ...
      [plotLAbasalline plotLAbasalline], ...
      [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))],'w-', ...
      [plotLAbasalline plotLAbasalline], ...
      [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)],'w-', ...
      [plotLAmidline plotLAmidline], ...
      [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))],'w-', ...
      [plotLAmidline plotLAmidline], ...
      [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)],'w-', ...
      [plotLAapicalline plotLAapicalline], ...
      [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))],'w-', ...
      [plotLAapicalline plotLAapicalline(tf)], ...
      [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)],'w-');
    set(gui.handles.axesVLAlinesright,'linewidth',2);
    hold(gui.handles.axesVLAright,'off');
  else
    hold(gui.handles.axesHLAright,'on');
    gui.handles.axesHLAlinesright = plot(gui.handles.axesHLAright,...
      [],[],'w-',[],[],'w',[],[],'w',[],[],'w',[],[],'w',[],[],'w');
    set(gui.handles.axesHLAlinesright,'linewidth',2);
    hold(gui.handles.axesHLAright,'off');
    %VLA
    hold(gui.handles.axesVLAright,'on');
    gui.handles.axesVLAlinesright = plot(gui.handles.axesVLAright,...
      [],[],'w-',[],[],'w',[],[],'w',[],[],'w',[],[],'w',[],[],'w');
    set(gui.handles.axesVLAlinesright,'linewidth',2);
    hold(gui.handles.axesVLAright,'off');
  end
end
if isequal(panel,'right') || isequal(panel,'none')
  %left
  hold(gui.handles.axesHLAleft,'on');
  gui.handles.axesHLAlinesleft = plot(gui.handles.axesHLAleft, ...
    0,0,'w-',0,0,'w',0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  set(gui.handles.axesHLAlinesleft,'linewidth',2);
  hold(gui.handles.axesHLAleft,'off');
  hold(gui.handles.axesVLAleft,'on');
  gui.handles.axesVLAlinesleft = plot(gui.handles.axesVLAleft, ...
    0,0,'w-',0,0,'w',0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  set(gui.handles.axesVLAlinesleft,'linewidth',2);
  hold(gui.handles.axesVLAleft,'off');
end
if isequal(panel,'left') || isequal(panel,'none')
  %right
  hold(gui.handles.axesHLAright,'on');
  gui.handles.axesHLAlinesright = plot(gui.handles.axesHLAright,...
    0,0,'w-',0,0,'w',0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  set(gui.handles.axesHLAlinesright,'linewidth',2);
  hold(gui.handles.axesHLAright,'off');
  hold(gui.handles.axesVLAright,'on');
  gui.handles.axesVLAlinesright = plot(gui.handles.axesVLAright,...
    0,0,'w-',0,0,'w',0,0,'w',0,0,'w',0,0,'w',0,0,'w');
  set(gui.handles.axesVLAlinesright,'linewidth',2);
  hold(gui.handles.axesVLAright,'off');
end

%--------------------
function updateimages
%--------------------
%update the images
global DATA

gui = DATA.GUI.SpectPlot2d;

if gui.showno(1) ~= 0
  switch gui.showimagetype{1}
    case 'stress'
      tf = 1;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount = gui.maxcountstress;
        normmincount = gui.mincountstress;
      else
        normmaxcount = max(gui.maxcountstress);
        normmincount = min(gui.mincountstress);
      end
      set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagestress(:,:,tf));
      set(gui.handles.axesSAbasalleft,'xlim',[1 size(gui.SAbasalimagestress(:,:,tf),2)]);
      set(gui.handles.axesSAbasalleft,'ylim',[1 size(gui.SAbasalimagestress(:,:,tf),1)]);
      set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagestress(:,:,tf));
      set(gui.handles.axesSAmidleft,'xlim',[1 size(gui.SAmidimagestress(:,:,tf),2)]);
      set(gui.handles.axesSAmidleft,'ylim',[1 size(gui.SAmidimagestress(:,:,tf),1)]);
      set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagestress(:,:,tf));
      set(gui.handles.axesSAapicalleft,'xlim',[1 size(gui.SAapicalimagestress(:,:,tf),2)]);
      set(gui.handles.axesSAapicalleft,'ylim',[1 size(gui.SAapicalimagestress(:,:,tf),1)]);
      set(gui.handles.axesHLAimageleft,'cdata',gui.HLAimagestress(:,:,tf));
      set(gui.handles.axesHLAleft,'xlim',[1 size(gui.HLAimagestress(:,:,tf),2)]);
      set(gui.handles.axesHLAleft,'ylim',[1 size(gui.HLAimagestress(:,:,tf),1)]);
      set(gui.handles.axesVLAimageleft,'cdata',gui.VLAimagestress(:,:,tf));
      set(gui.handles.axesVLAleft,'xlim',[1 size(gui.VLAimagestress(:,:,tf),2)]);
      set(gui.handles.axesVLAleft,'ylim',[1 size(gui.VLAimagestress(:,:,tf),1)]);
      textleft = 'Stress';
    case 'stressgated'
      tf = gui.handles.tf;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount = gui.maxcountstressgated(tf);
        normmincount = gui.mincountstressgated(tf);
      else
        normmaxcount = max(gui.maxcountstressgated);
        normmincount = min(gui.mincountstressgated);
      end
      set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagestressgated(:,:,tf));
      set(gui.handles.axesSAbasalleft,'xlim',[1 size(gui.SAbasalimagestressgated(:,:,tf),2)]);
      set(gui.handles.axesSAbasalleft,'ylim',[1 size(gui.SAbasalimagestressgated(:,:,tf),1)]);
      set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagestressgated(:,:,tf));
      set(gui.handles.axesSAmidleft,'xlim',[1 size(gui.SAmidimagestressgated(:,:,tf),2)]);
      set(gui.handles.axesSAmidleft,'ylim',[1 size(gui.SAmidimagestressgated(:,:,tf),1)]);
      set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagestressgated(:,:,tf));
      set(gui.handles.axesSAapicalleft,'xlim',[1 size(gui.SAapicalimagestressgated(:,:,tf),2)]);
      set(gui.handles.axesSAapicalleft,'ylim',[1 size(gui.SAapicalimagestressgated(:,:,tf),1)]);
      set(gui.handles.axesHLAimageleft,'cdata',gui.HLAimagestressgated(:,:,tf));
      set(gui.handles.axesHLAleft,'xlim',[1 size(gui.HLAimagestressgated(:,:,tf),2)]);
      set(gui.handles.axesHLAleft,'ylim',[1 size(gui.HLAimagestressgated(:,:,tf),1)]);
      set(gui.handles.axesVLAimageleft,'cdata',gui.VLAimagestressgated(:,:,tf));
      set(gui.handles.axesVLAleft,'xlim',[1 size(gui.VLAimagestressgated(:,:,tf),2)]);
      set(gui.handles.axesVLAleft,'ylim',[1 size(gui.VLAimagestressgated(:,:,tf),1)]);
      textleft = 'Stress gated';
    case 'restgated'
      tf = gui.handles.tf;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount = gui.maxcountrestgated(tf);
        normmincount = gui.mincountrestgated(tf);
      else
        normmaxcount = max(gui.maxcountrestgated);
        normmincount = min(gui.mincountrestgated);
      end
      set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagerestgated(:,:,tf));
      set(gui.handles.axesSAbasalleft,'xlim',[1 size(gui.SAbasalimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesSAbasalleft,'ylim',[1 size(gui.SAbasalimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagerestgated(:,:,tf));
      set(gui.handles.axesSAmidleft,'xlim',[1 size(gui.SAmidimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesSAmidleft,'ylim',[1 size(gui.SAmidimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagerestgated(:,:,tf));
      set(gui.handles.axesSAapicalleft,'xlim',[1 size(gui.SAapicalimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesSAapicalleft,'ylim',[1 size(gui.SAapicalimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesHLAimageleft,'cdata',gui.HLAimagerestgated(:,:,tf));
      set(gui.handles.axesHLAleft,'xlim',[1 size(gui.HLAimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesHLAleft,'ylim',[1 size(gui.HLAimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesVLAimageleft,'cdata',gui.VLAimagerestgated(:,:,tf));
      set(gui.handles.axesVLAleft,'xlim',[1 size(gui.VLAimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesVLAleft,'ylim',[1 size(gui.VLAimagerestgated(:,:,tf),1)]);
      textleft = 'Rest gated';
  end
  
  switch gui.showimagetype{1}    
    case {'stress','stressgated','restgated'}
      set(gui.handles.axesSAbasalleft,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesSAmidleft,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesSAapicalleft,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesHLAleft,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesVLAleft,'clim',[normmincount normmaxcount]);
      set(gui.handles.textleft,'String',textleft);      
      set(gui.handles.textSAbasalleft,'String','Basal SA');
      set(gui.handles.textSAmidleft,'String','Midventricular SA');
      set(gui.handles.textSAapicalleft,'String','Apical SA');
      set(gui.handles.textHLAleft,'String','HLA');
      set(gui.handles.textVLAleft,'String','VLA');
      set(gui.handles.uipanelcorrectLVleft,'visible','on');

    otherwise
      %plot the images
      set(gui.handles.axesSAbasalimageleft,'cdata',[]);
      set(gui.handles.axesSAmidimageleft,'cdata',[]);
      set(gui.handles.axesSAapicalimageleft,'cdata',[]);
      set(gui.handles.textSAbasalleft,'String',[]);
      set(gui.handles.textSAmidleft,'String',[]);
      set(gui.handles.textSAapicalleft,'String',[]);
      %plot Long axis views
      set(gui.handles.axesHLAimageleft,'cdata',[]);
      set(gui.handles.axesVLAimageleft,'cdata',[]);
      set(gui.handles.textHLAleft,'String',[]);
      set(gui.handles.textVLAleft,'String',[]);
      set(gui.handles.textleft,'String',[]);
  end
  colormap(gui.handles.axesSAapicalleft,'spect');
  colormap(gui.handles.axesSAmidleft,'spect');
  colormap(gui.handles.axesSAbasalleft,'spect');
  colormap(gui.handles.axesHLAleft,'spect');
  colormap(gui.handles.axesVLAleft,'spect');
else
  %no images to plot 
  set(gui.handles.axesSAbasalimageleft,'cdata',[]);
  set(gui.handles.axesSAmidimageleft,'cdata',[]);
  set(gui.handles.axesSAapicalimageleft,'cdata',[]);
%   axis(gui.handles.axesSAbasalleft,'off');
%   axis(gui.handles.axesSAmidleft,'off');
%   axis(gui.handles.axesSAapicalleft,'off');
  set(gui.handles.textSAbasalleft,'String',[]);
  set(gui.handles.textSAmidleft,'String',[]);
  set(gui.handles.textSAapicalleft,'String',[]);
  set(gui.handles.axesHLAimageleft,'cdata',[]);
  set(gui.handles.axesVLAimageleft,'cdata',[]);
%   axis(gui.handles.axesHLAleft,'off');
%   axis(gui.handles.axesVLAleft,'off');
  set(gui.handles.textHLAleft,'String',[]);
  set(gui.handles.textVLAleft,'String',[]);
  set(gui.handles.textleft,'String',[]);
  set(gui.handles.uipanelcorrectLVleft,'visible','off');
end

if gui.showno(2) ~= 0
  switch gui.showimagetype{2}
    case 'rest'
      tf = 1;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount =  gui.maxcountrest(tf);
        normmincount =  gui.mincountrest(tf);
      else
        normmaxcount = max(gui.maxcountrest);
        normmincount = min(gui.mincountrest);
      end
      set(gui.handles.axesSAbasalimageright,'cdata',gui.SAbasalimagerest(:,:,tf));
      set(gui.handles.axesSAbasalright,'xlim',[1 size(gui.SAbasalimagerest(:,:,tf),2)]);
      set(gui.handles.axesSAbasalright,'ylim',[1 size(gui.SAbasalimagerest(:,:,tf),1)]);
      set(gui.handles.axesSAmidimageright,'cdata',gui.SAmidimagerest(:,:,tf));
      set(gui.handles.axesSAmidright,'xlim',[1 size(gui.SAmidimagerest(:,:,tf),2)]);
      set(gui.handles.axesSAmidright,'ylim',[1 size(gui.SAmidimagerest(:,:,tf),1)]);
      set(gui.handles.axesSAapicalimageright,'cdata',gui.SAapicalimagerest(:,:,tf));
      set(gui.handles.axesSAapicalright,'xlim',[1 size(gui.SAapicalimagerest(:,:,tf),2)]);
      set(gui.handles.axesSAapicalright,'ylim',[1 size(gui.SAapicalimagerest(:,:,tf),1)]);
      set(gui.handles.axesHLAimageright,'cdata',gui.HLAimagerest(:,:,tf));
      set(gui.handles.axesHLAright,'xlim',[1 size(gui.HLAimagerest(:,:,tf),2)]);
      set(gui.handles.axesHLAright,'ylim',[1 size(gui.HLAimagerest(:,:,tf),1)]);
      set(gui.handles.axesVLAimageright,'cdata',gui.VLAimagerest(:,:,tf));
      set(gui.handles.axesVLAright,'xlim',[1 size(gui.VLAimagerest(:,:,tf),2)]);
      set(gui.handles.axesVLAright,'ylim',[1 size(gui.VLAimagerest(:,:,tf),1)]);
      textright = 'Rest';
    case 'restgated'
      tf = gui.handles.tf;
      if get(gui.handles.radiobuttonnormeachtf,'value')
        normmaxcount =  gui.maxcountrestgated(tf);
        normmincount =  gui.mincountrestgated(tf);
      else
        normmaxcount = max(gui.maxcountrestgated);
        normmincount = min(gui.mincountrestgated);
      end
      set(gui.handles.axesSAbasalimageright,'cdata',gui.SAbasalimagerestgated(:,:,tf));
      set(gui.handles.axesSAbasalright,'xlim',[1 size(gui.SAbasalimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesSAbasalright,'ylim',[1 size(gui.SAbasalimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesSAmidimageright,'cdata',gui.SAmidimagerestgated(:,:,tf));
      set(gui.handles.axesSAmidright,'xlim',[1 size(gui.SAmidimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesSAmidright,'ylim',[1 size(gui.SAmidimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesSAapicalimageright,'cdata',gui.SAapicalimagerestgated(:,:,tf));
      set(gui.handles.axesSAapicalright,'xlim',[1 size(gui.SAapicalimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesSAapicalright,'ylim',[1 size(gui.SAapicalimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesHLAimageright,'cdata',gui.HLAimagerestgated(:,:,tf));
      set(gui.handles.axesHLAright,'xlim',[1 size(gui.HLAimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesHLAright,'ylim',[1 size(gui.HLAimagerestgated(:,:,tf),1)]);
      set(gui.handles.axesVLAimageright,'cdata',gui.VLAimagerestgated(:,:,tf));
      set(gui.handles.axesVLAright,'xlim',[1 size(gui.VLAimagerestgated(:,:,tf),2)]);
      set(gui.handles.axesVLAright,'ylim',[1 size(gui.VLAimagerestgated(:,:,tf),1)]);
      textright = 'Rest gated';
  end  
  switch gui.showimagetype{2}
    case {'rest','restgated'}
      set(gui.handles.axesSAbasalright,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesSAmidright,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesSAapicalright,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesHLAright,'clim',[normmincount normmaxcount]);
      set(gui.handles.axesVLAright,'clim',[normmincount normmaxcount]);
      set(gui.handles.textright,'String',textright);   
      set(gui.handles.textSAbasalright,'String','Basal SA');
      set(gui.handles.textSAmidright,'String','Midventricular SA');
      set(gui.handles.textSAapicalright,'String','Apical SA');
      set(gui.handles.textHLAright,'String','HLA');
      set(gui.handles.textVLAright,'String','VLA');
      set(gui.handles.uipanelcorrectLVright,'visible','on');
      
    otherwise
      %plot the SA images
      set(gui.handles.axesSAbasalimageright,'cdata',[]);
      set(gui.handles.axesSAmidimageright,'cdata',[]);
      set(gui.handles.axesSAapicalimageright,'cdata',[]);
      set(gui.handles.textSAbasalright,'String',[]);
      set(gui.handles.textSAmidright,'String',[]);
      set(gui.handles.textSAapicalright,'String',[]);
      %plot Long axis views
      set(gui.handles.axesHLAimageright,'cdata',[]);
      set(gui.handles.axesVLAimageright,'cdata',[]);
      set(gui.handles.textHLAright,'String',[]);
      set(gui.handles.textVLAright,'String',[]);
      set(gui.handles.textright,'String',[]);
  end
else
  %no images to plot
  set(gui.handles.axesSAbasalimageright,'cdata',[]);
  set(gui.handles.axesSAmidimageright,'cdata',[]);
  set(gui.handles.axesSAapicalimageright,'cdata',[]);
  axis(gui.handles.axesSAbasalright,'off');
  axis(gui.handles.axesSAmidright,'off');
  axis(gui.handles.axesSAapicalright,'off');
  set(gui.handles.textSAbasalright,'String',[]);
  set(gui.handles.textSAmidright,'String',[]);
  set(gui.handles.textSAapicalright,'String',[]);
  set(gui.handles.axesHLAimageright,'cdata',[]);
  set(gui.handles.axesVLAimageright,'cdata',[]);
  axis(gui.handles.axesHLAright,'off');
  axis(gui.handles.axesVLAright,'off');
  set(gui.handles.textHLAright,'String',[]);
  set(gui.handles.textVLAright,'String',[]);
  set(gui.handles.textright,'String',[]);
  set(gui.handles.uipanelcorrectLVright,'visible','off');
end

%--------------------------------
function updateahasections(panel)
%--------------------------------
%update the 17 aha section model

global DATA

gui = DATA.GUI.SpectPlot2d;

if isequal(panel,'left') || isequal(panel,'both')
  if gui.showno(1) ~= 0
    switch gui.showimagetype{1}
      case 'stress'
        if gui.LVsegmentationstress
          LVcenterxapical = mean(get(gui.handles.axesSAapicalendoleft,'xdata'));
          LVcenteryapical = mean(get(gui.handles.axesSAapicalendoleft,'ydata'));
          LVcenterxmid = mean(get(gui.handles.axesSAmidendoleft,'xdata'));
          LVcenterymid = mean(get(gui.handles.axesSAmidendoleft,'ydata'));
          LVcenterxbasal = mean(get(gui.handles.axesSAbasalendoleft,'xdata'));
          LVcenterybasal = mean(get(gui.handles.axesSAbasalendoleft,'ydata'));
          b = sqrt(0.75);
          a = 0.5;
          c = 1/sqrt(2);
          %apical
          set(gui.handles.axesSAapicalahalineleft(1),'xdata',LVcenterxapical*[1-c 1+c]);
          set(gui.handles.axesSAapicalahalineleft(1),'ydata',LVcenteryapical*[1-c 1+c]);
          set(gui.handles.axesSAapicalahalineleft(2),'xdata',LVcenterxapical*[1-c 1+c]);
          set(gui.handles.axesSAapicalahalineleft(2),'ydata',LVcenteryapical*[1+c 1-c]);
          %mid
          set(gui.handles.axesSAmidahalineleft(1),'xdata',LVcenterxmid*[0 2]);
          set(gui.handles.axesSAmidahalineleft(1),'ydata',LVcenterymid*[1 1]);
          set(gui.handles.axesSAmidahalineleft(2),'xdata',LVcenterxmid*[1-1*a 1+1*a]);
          set(gui.handles.axesSAmidahalineleft(2),'ydata',LVcenterymid*[1-1*b 1+1*b]);
          set(gui.handles.axesSAmidahalineleft(3),'xdata',LVcenterxmid*[1-1*a 1+1*a]);
          set(gui.handles.axesSAmidahalineleft(3),'ydata',LVcenterymid*[1+1*b 1-1*b]);
          %basal
          set(gui.handles.axesSAbasalahalineleft(1),'xdata',LVcenterxbasal*[0 2]);
          set(gui.handles.axesSAbasalahalineleft(1),'ydata',LVcenterybasal*[1 1]);
          set(gui.handles.axesSAbasalahalineleft(2),'xdata',LVcenterxbasal*[1-1*a 1+1*a]);
          set(gui.handles.axesSAbasalahalineleft(2),'ydata',LVcenterybasal*[1-1*b 1+1*b]);
          set(gui.handles.axesSAbasalahalineleft(3),'xdata',LVcenterxbasal*[1-1*a 1+1*a]);
          set(gui.handles.axesSAbasalahalineleft(3),'ydata',LVcenterybasal*[1+1*b 1-1*b]);
          %HLA
          xline = get(gui.handles.axesHLAleft,'xlim');
          ylim = get(gui.handles.axesHLAleft,'ylim');
          yline = [(ylim(2)-gui.SAslicesstress{1}(end))*[1 1] ; (ylim(2)-gui.SAslicesstress{2}(end))*[1 1] ; ...
            (ylim(2)-gui.SAslicesstress{3}(end))*[1 1] ; (ylim(2)-gui.SAslicesstress{3}(1))*[1 1]];
          for k = 1:4
            set(gui.handles.axesHLAahalineleft(k),'xdata',xline);
            set(gui.handles.axesHLAahalineleft(k),'ydata',yline(k,:));
          end
          %VLA
          yline = get(gui.handles.axesVLAleft,'xlim');
          xline = [gui.SAslicesstress{1}(end)*[1 1] ; gui.SAslicesstress{2}(end)*[1 1] ; ...
            gui.SAslicesstress{3}(end)*[1 1] ; gui.SAslicesstress{3}(1)*[1 1]];
          for k = 1:4
            set(gui.handles.axesVLAahalineleft(k),'xdata',xline(k,:));
            set(gui.handles.axesVLAahalineleft(k),'ydata',yline);
          end
        end
    end
  end
end
if isequal(panel,'right') || isequal(panel,'both')
  if gui.showno(2) ~= 0
    switch gui.showimagetype{2}
      case 'rest'
        if gui.LVsegmentationrest
          LVcenterxapical = mean(get(gui.handles.axesSAapicalendoright,'xdata'));
          LVcenteryapical = mean(get(gui.handles.axesSAapicalendoright,'ydata'));
          LVcenterxmid = mean(get(gui.handles.axesSAmidendoright,'xdata'));
          LVcenterymid = mean(get(gui.handles.axesSAmidendoright,'ydata'));
          LVcenterxbasal = mean(get(gui.handles.axesSAbasalendoright,'xdata'));
          LVcenterybasal = mean(get(gui.handles.axesSAbasalendoright,'ydata'));
          b = sqrt(0.75);
          a = 0.5;
          c = 1/sqrt(2);
          %apical
          set(gui.handles.axesSAapicalahalineright(1),'xdata',LVcenterxapical*[1-c 1+c]);
          set(gui.handles.axesSAapicalahalineright(1),'ydata',LVcenteryapical*[1-c 1+c]);
          set(gui.handles.axesSAapicalahalineright(2),'xdata',LVcenterxapical*[1-c 1+c]);
          set(gui.handles.axesSAapicalahalineright(2),'ydata',LVcenteryapical*[1+c 1-c]);
          %mid
          set(gui.handles.axesSAmidahalineright(1),'xdata',LVcenterxmid*[0 2]);
          set(gui.handles.axesSAmidahalineright(1),'ydata',LVcenterymid*[1 1]);
          set(gui.handles.axesSAmidahalineright(2),'xdata',LVcenterxmid*[1-1*a 1+1*a]);
          set(gui.handles.axesSAmidahalineright(2),'ydata',LVcenterymid*[1-1*b 1+1*b]);
          set(gui.handles.axesSAmidahalineright(3),'xdata',LVcenterxmid*[1-1*a 1+1*a]);
          set(gui.handles.axesSAmidahalineright(3),'ydata',LVcenterymid*[1+1*b 1-1*b]);
          %basal
          set(gui.handles.axesSAbasalahalineright(1),'xdata',LVcenterxbasal*[0 2]);
          set(gui.handles.axesSAbasalahalineright(1),'ydata',LVcenterybasal*[1 1]);
          set(gui.handles.axesSAbasalahalineright(2),'xdata',LVcenterxbasal*[1-1*a 1+1*a]);
          set(gui.handles.axesSAbasalahalineright(2),'ydata',LVcenterybasal*[1-1*b 1+1*b]);
          set(gui.handles.axesSAbasalahalineright(3),'xdata',LVcenterxbasal*[1-1*a 1+1*a]);
          set(gui.handles.axesSAbasalahalineright(3),'ydata',LVcenterybasal*[1+1*b 1-1*b]);
          %HLA
          xline = get(gui.handles.axesHLAright,'xlim');
          ylim = get(gui.handles.axesHLAright,'ylim');
          yline = [(ylim(2)-gui.SAslicesrest{1}(end))*[1 1] ; (ylim(2)-gui.SAslicesrest{2}(end))*[1 1] ; ...
            (ylim(2)-gui.SAslicesrest{3}(end))*[1 1] ; (ylim(2)-gui.SAslicesrest{3}(1))*[1 1]];
          for k = 1:4
            set(gui.handles.axesHLAahalineright(k),'xdata',xline);
            set(gui.handles.axesHLAahalineright(k),'ydata',yline(k,:));
          end
          %VLA
          yline = get(gui.handles.axesVLAright,'xlim');
          xline = [gui.SAslicesrest{1}(end)*[1 1] ; gui.SAslicesrest{2}(end)*[1 1] ; ...
            gui.SAslicesrest{3}(end)*[1 1] ; gui.SAslicesrest{3}(1)*[1 1]];
          for k = 1:4
            set(gui.handles.axesVLAahalineright(k),'xdata',xline(k,:));
            set(gui.handles.axesVLAahalineright(k),'ydata',yline);
          end
        end
    end
  end
end


%------------------------------------
function updateSAintersections(panel)
%------------------------------------

global DATA

gui = DATA.GUI.SpectPlot2d;

if isequal(panel,'left') || isequal(panel,'both')
  if gui.showno(1) ~= 0
    switch gui.showimagetype{1}
      case 'stress'
        plotLAbasalline = mean(gui.SAslicesstress{3});
        plotLAmidline = mean(gui.SAslicesstress{2});
        plotLAapicalline = mean(gui.SAslicesstress{1});
      case 'stressgated'
        tf = gui.handles.tf;
        plotLAbasalline = mean(gui.SAslicesstressgated{3,tf});
        plotLAmidline = mean(gui.SAslicesstressgated{2,tf});
        plotLAapicalline = mean(gui.SAslicesstressgated{1,tf});
      case 'restgated'
        tf = gui.handles.tf;
        plotLAbasalline = mean(gui.SAslicesrestgated{3,tf});
        plotLAmidline = mean(gui.SAslicesrestgated{2,tf});
        plotLAapicalline = mean(gui.SAslicesrestgated{1,tf});
    end
    %HLA
    HLAxlim = get(gui.handles.axesHLAleft,'xlim');
    ylim = get(gui.handles.axesHLAleft,'ylim');
    xdataplot{1,1} = [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))];
    xdataplot{2,1} = [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)];
    xdataplot{3,1} = [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))];
    xdataplot{4,1} = [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)];
    xdataplot{5,1} = [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))];
    xdataplot{6,1} = [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)];
    ydataplot{1,1} = (ylim(2)-plotLAbasalline)*[1 1];
    ydataplot{2,1} = (ylim(2)-plotLAbasalline)*[1 1];
    ydataplot{3,1} = (ylim(2)-plotLAmidline)*[1 1];
    ydataplot{4,1} = (ylim(2)-plotLAmidline)*[1 1];
    ydataplot{5,1} = (ylim(2)-plotLAapicalline)*[1 1];
    ydataplot{6,1} = (ylim(2)-plotLAapicalline)*[1 1];
    if isempty(gui.handles.axesHLAlinesleft)  
      hold(gui.handles.axesHLAleft,'on');
      gui.handles.axesHLAlinesleft = plot(gui.handles.axesHLAleft,...
        xdataplot{1},ydataplot{1},'w',xdataplot{2},ydataplot{2},'w',xdataplot{3},ydataplot{3},'w', ...
        xdataplot{4},ydataplot{4},'w',xdataplot{5},ydataplot{5},'w',xdataplot{6},ydataplot{6},'w');
      set(gui.handles.axesHLAlinesleft,'linewidth',2);
      hold(gui.handles.axesHLAleft,'off');
    else
      for k = 1:size(xdataplot,1)
        set(gui.handles.axesHLAlinesleft(k),'xdata',xdataplot{k});
        set(gui.handles.axesHLAlinesleft(k),'ydata',ydataplot{k});
      end
    end
    %VLA
    VLAylim = get(gui.handles.axesVLAleft,'ylim');
    xdataplot{1,1} = plotLAbasalline*[1 1];
    xdataplot{2,1} = plotLAbasalline*[1 1];
    xdataplot{3,1} = plotLAmidline*[1 1];
    xdataplot{4,1} = plotLAmidline*[1 1];
    xdataplot{5,1} = plotLAapicalline*[1 1];
    xdataplot{6,1} = plotLAapicalline*[1 1];
    ydataplot{1,1} = [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))];
    ydataplot{2,1} = [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)];
    ydataplot{3,1} = [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))];
    ydataplot{4,1} = [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)];
    ydataplot{5,1} = [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))];
    ydataplot{6,1} = [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)];
    if isempty(gui.handles.axesVLAlinesleft)
      hold(gui.handles.axesVLAleft,'on');
      gui.handles.axesVLAlinesleft = plot(gui.handles.axesVLAleft,...
        xdataplot{1},ydataplot{1},'w',xdataplot{2},ydataplot{2},'w',xdataplot{3},ydataplot{3},'w', ...
        xdataplot{4},ydataplot{4},'w',xdataplot{5},ydataplot{5},'w',xdataplot{6},ydataplot{6},'w');
      set(gui.handles.axesVLAlinesleft,'linewidth',2);
      hold(gui.handles.axesVLAleft,'off');
    else
      for k = 1:size(xdataplot,1)
        set(gui.handles.axesVLAlinesleft(k),'xdata',xdataplot{k});
        set(gui.handles.axesVLAlinesleft(k),'ydata',ydataplot{k});
      end
    end
  elseif ~isempty(gui.handles.axesVLAlinesleft)
    for k = 1:size(gui.handles.axesHLAlinesleft,1)
      set(gui.handles.axesHLAlinesleft(k),'xdata',[]);
      set(gui.handles.axesHLAlinesleft(k),'ydata',[]);
      set(gui.handles.axesVLAlinesleft(k),'xdata',[]);
      set(gui.handles.axesVLAlinesleft(k),'ydata',[]);
    end
  end
end
if isequal(panel,'right') || isequal(panel,'both')
  if gui.showno(2) ~= 0
    switch gui.showimagetype{2}
      case 'rest'
        plotLAbasalline = mean(gui.SAslicesrest{3});
        plotLAmidline = mean(gui.SAslicesrest{2});
        plotLAapicalline = mean(gui.SAslicesrest{1});
      case 'restgated'
        tf = gui.handles.tf;
        plotLAbasalline = mean(gui.SAslicesrestgated{3,tf});
        plotLAmidline = mean(gui.SAslicesrestgated{2,tf});
        plotLAapicalline = mean(gui.SAslicesrestgated{1,tf});
    end
    %HLA
    HLAxlim = get(gui.handles.axesHLAright,'xlim');
    ylim = get(gui.handles.axesHLAright,'ylim');
    xdataplot{1,1} = [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))];
    xdataplot{2,1} = [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)];
    xdataplot{3,1} = [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))];
    xdataplot{4,1} = [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)];
    xdataplot{5,1} = [HLAxlim(1) round(0.05*abs(HLAxlim(2)-HLAxlim(1)))];
    xdataplot{6,1} = [HLAxlim(end)-round(0.05*abs(HLAxlim(2)-HLAxlim(1))) HLAxlim(end)];
    ydataplot{1,1} = (ylim(2)-plotLAbasalline)*[1 1];
    ydataplot{2,1} = (ylim(2)-plotLAbasalline)*[1 1];
    ydataplot{3,1} = (ylim(2)-plotLAmidline)*[1 1];
    ydataplot{4,1} = (ylim(2)-plotLAmidline)*[1 1];
    ydataplot{5,1} = (ylim(2)-plotLAapicalline)*[1 1];
    ydataplot{6,1} = (ylim(2)-plotLAapicalline)*[1 1];
    if isempty(gui.handles.axesHLAlinesright)  
      hold(gui.handles.axesHLAright,'on');
      gui.handles.axesHLAlinesright = plot(gui.handles.axesHLAright,...
        xdataplot{1},ydataplot{1},'w',xdataplot{2},ydataplot{2},'w',xdataplot{3},ydataplot{3},'w', ...
        xdataplot{4},ydataplot{4},'w',xdataplot{5},ydataplot{5},'w',xdataplot{6},ydataplot{6},'w');
      set(gui.handles.axesHLAlinesright,'linewidth',2);
      hold(gui.handles.axesHLAright,'off');
    else
      for k = 1:size(xdataplot,1)
        set(gui.handles.axesHLAlinesright(k),'xdata',xdataplot{k});
        set(gui.handles.axesHLAlinesright(k),'ydata',ydataplot{k});
      end
    end
    %VLA
    VLAylim = get(gui.handles.axesVLAright,'ylim');
    xdataplot{1,1} = [plotLAbasalline plotLAbasalline];
    xdataplot{2,1} = [plotLAbasalline plotLAbasalline];
    xdataplot{3,1} = [plotLAmidline plotLAmidline];
    xdataplot{4,1} = [plotLAmidline plotLAmidline];
    xdataplot{5,1} = [plotLAapicalline plotLAapicalline];
    xdataplot{6,1} = [plotLAapicalline plotLAapicalline];
    ydataplot{1,1} = [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))];
    ydataplot{2,1} = [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)];
    ydataplot{3,1} = [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))];
    ydataplot{4,1} = [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)];
    ydataplot{5,1} = [VLAylim(1) round(0.05*abs(VLAylim(2)-VLAylim(1)))];
    ydataplot{6,1} = [VLAylim(end)-round(0.05*abs(VLAylim(2)-VLAylim(1))) VLAylim(end)];
    if isempty(gui.handles.axesVLAlinesright)
      hold(gui.handles.axesVLAright,'on');
      gui.handles.axesVLAlinesright = plot(gui.handles.axesVLAright,...
        xdataplot{1},ydataplot{1},'w',xdataplot{2},ydataplot{2},'w',xdataplot{3},ydataplot{3},'w', ...
        xdataplot{4},ydataplot{4},'w',xdataplot{5},ydataplot{5},'w',xdataplot{6},ydataplot{6},'w');
      set(gui.handles.axesVLAlinesright,'linewidth',2);
      hold(gui.handles.axesVLAright,'off');
    else
      for k = 1:length(xdataplot)
        set(gui.handles.axesVLAlinesright(k),'xdata',xdataplot{k});
        set(gui.handles.axesVLAlinesright(k),'ydata',ydataplot{k});
      end
    end
  elseif ~isempty(gui.handles.axesVLAlinesright)
    for k = 1:size(gui.handles.axesHLAlinesright,1)
      set(gui.handles.axesHLAlinesright(k),'xdata',[]);
      set(gui.handles.axesHLAlinesright(k),'ydata',[]);
      set(gui.handles.axesVLAlinesright(k),'xdata',[]);
      set(gui.handles.axesVLAlinesright(k),'ydata',[]);
    end
  end
end


%----------------------------
function updateLVsegmentation
%----------------------------
%update the endo- and epicardial segmentation

global DATA

gui = DATA.GUI.SpectPlot2d;

if gui.showno(1) ~= 0
  switch gui.showimagetype{1}
    case 'stress'
      if gui.LVsegmentationstress
        endoX = gui.EndoXstress;
        endoY = gui.EndoYstress;
        epiX = gui.EpiXstress;
        epiY = gui.EpiYstress;
        tf = 1;
        SAslices = gui.SAslicesstress;
        HLAmyocardiumx = gui.HLAmyocardiumxstress;
        HLAmyocardiumy = gui.HLAmyocardiumystress;
        VLAmyocardiumx = gui.VLAmyocardiumxstress;
        VLAmyocardiumy = gui.VLAmyocardiumystress;
      end
    case 'stressgated'
      if gui.LVsegmentationstressgated
        endoX = gui.EndoXstressgated;
        endoY = gui.EndoYstressgated;
        epiX = gui.EpiXstressgated;
        epiY = gui.EpiYstressgated;
        tf = gui.handles.tf;
        for s = 1:3
          SAslices{s} = gui.SAslicesstressgated{s,tf};
        end
        HLAmyocardiumx = gui.HLAmyocardiumxstressgated;
        HLAmyocardiumy = gui.HLAmyocardiumystressgated;
        VLAmyocardiumx = gui.VLAmyocardiumxstressgated;
        VLAmyocardiumy = gui.VLAmyocardiumystressgated;
      end
    case 'restgated'
      if gui.LVsegmentationrestgated
        endoX = gui.EndoXrestgated;
        endoY = gui.EndoYrestgated;
        epiX = gui.EpiXrestgated;
        epiY = gui.EpiYrestgated;
        tf = gui.handles.tf;
        for s = 1:3
          SAslices{s} = gui.SAslicesrestgated{s,tf};
        end
        HLAmyocardiumx = gui.HLAmyocardiumxrestgated;
        HLAmyocardiumy = gui.HLAmyocardiumyrestgated;
        VLAmyocardiumx = gui.VLAmyocardiumxrestgated;
        VLAmyocardiumy = gui.VLAmyocardiumyrestgated;
      end
  end
  if isempty(gui.handles.axesSAbasalendoleft)
    %basal
    hold(gui.handles.axesSAbasalleft,'on');
    gui.handles.axesSAbasalendoleft = plot(gui.handles.axesSAbasalleft,endoY(:,tf,round(mean(SAslices{3}))),endoX(:,tf,round(mean(SAslices{3}))),'w-');
    set(gui.handles.axesSAbasalendoleft,'linewidth',2);
    gui.handles.axesSAbasalepileft = plot(gui.handles.axesSAbasalleft,epiY(:,tf,round(mean(SAslices{3}))),epiX(:,tf,round(mean(SAslices{3}))),'w-');
    set(gui.handles.axesSAbasalepileft,'linewidth',2);
    hold(gui.handles.axesSAbasalleft,'off');
    %mid
    hold(gui.handles.axesSAmidleft,'on');
    gui.handles.axesSAmidendoleft = plot(gui.handles.axesSAmidleft,endoY(:,tf,round(mean(SAslices{2}))),endoX(:,tf,round(mean(SAslices{2}))),'w-');
    set(gui.handles.axesSAmidendoleft,'linewidth',2);
    gui.handles.axesSAmidepileft = plot(gui.handles.axesSAmidleft,epiY(:,tf,round(mean(SAslices{2}))),epiX(:,tf,round(mean(SAslices{2}))),'w-');
    set(gui.handles.axesSAmidepileft,'linewidth',2);
    hold(gui.handles.axesSAmidleft,'off');
    %apical
    hold(gui.handles.axesSAapicalleft,'on');
    gui.handles.axesSAapicalendoleft = plot(gui.handles.axesSAapicalleft,endoY(:,tf,round(mean(SAslices{1}))),endoX(:,tf,round(mean(SAslices{1}))),'w-');
    set(gui.handles.axesSAapicalendoleft,'linewidth',2);
    gui.handles.axesSAapicalepileft = plot(gui.handles.axesSAapicalleft,epiY(:,tf,round(mean(SAslices{1}))),epiX(:,tf,round(mean(SAslices{1}))),'w-');
    set(gui.handles.axesSAapicalepileft,'linewidth',2);
    hold(gui.handles.axesSAapicalleft,'off');
    %HLA
    hold(gui.handles.axesHLAleft,'on');
    gui.handles.axesHLAmyocardiumleft = plot(gui.handles.axesHLAleft,HLAmyocardiumy{tf},HLAmyocardiumx{tf},'w-');
    set(gui.handles.axesHLAmyocardiumleft,'linewidth',2);
    hold(gui.handles.axesHLAleft,'off');
    %VLA
    hold(gui.handles.axesVLAleft,'on');
    gui.handles.axesVLAmyocardiumleft = plot(gui.handles.axesVLAleft,VLAmyocardiumy{tf},VLAmyocardiumx{tf},'w-');
    set(gui.handles.axesVLAmyocardiumleft,'linewidth',2);
    hold(gui.handles.axesVLAleft,'off');
  else
    %basal
    set(gui.handles.axesSAbasalendoleft,'xdata',endoY(:,tf,round(mean(SAslices{3}))));
    set(gui.handles.axesSAbasalendoleft,'ydata',endoX(:,tf,round(mean(SAslices{3}))));
    set(gui.handles.axesSAbasalepileft,'xdata',epiY(:,tf,round(mean(SAslices{3}))));
    set(gui.handles.axesSAbasalepileft,'ydata',epiX(:,tf,round(mean(SAslices{3}))));
    %mid
    set(gui.handles.axesSAmidendoleft,'xdata',endoY(:,tf,round(mean(SAslices{2}))));
    set(gui.handles.axesSAmidendoleft,'ydata',endoX(:,tf,round(mean(SAslices{2}))));
    set(gui.handles.axesSAmidepileft,'xdata',epiY(:,tf,round(mean(SAslices{2}))));
    set(gui.handles.axesSAmidepileft,'ydata',epiX(:,tf,round(mean(SAslices{2}))));
    %apical
    set(gui.handles.axesSAapicalendoleft,'xdata',endoY(:,tf,round(mean(SAslices{1}))));
    set(gui.handles.axesSAapicalendoleft,'ydata',endoX(:,tf,round(mean(SAslices{1}))));
    set(gui.handles.axesSAapicalepileft,'xdata',epiY(:,tf,round(mean(SAslices{1}))));
    set(gui.handles.axesSAapicalepileft,'ydata',epiX(:,tf,round(mean(SAslices{1}))));
    %HLA
    set(gui.handles.axesHLAmyocardiumleft,'xdata',HLAmyocardiumy{tf});
    set(gui.handles.axesHLAmyocardiumleft,'ydata',HLAmyocardiumx{tf});
    %VLA
    set(gui.handles.axesVLAmyocardiumleft,'xdata',VLAmyocardiumy{tf});
    set(gui.handles.axesVLAmyocardiumleft,'ydata',VLAmyocardiumx{tf});
  end
else
  %basal
  set(gui.handles.axesSAbasalendoleft,'xdata',[]);
  set(gui.handles.axesSAbasalendoleft,'ydata',[]);
  set(gui.handles.axesSAbasalepileft,'xdata',[]);
  set(gui.handles.axesSAbasalepileft,'ydata',[]);
  %mid
  set(gui.handles.axesSAmidendoleft,'xdata',[]);
  set(gui.handles.axesSAmidendoleft,'ydata',[]);
  set(gui.handles.axesSAmidepileft,'xdata',[]);
  set(gui.handles.axesSAmidepileft,'ydata',[]);
  %apical
  set(gui.handles.axesSAapicalendoleft,'xdata',[]);
  set(gui.handles.axesSAapicalendoleft,'ydata',[]);
  set(gui.handles.axesSAapicalepileft,'xdata',[]);
  set(gui.handles.axesSAapicalepileft,'ydata',[]);
  %HLA
  set(gui.handles.axesHLAmyocardiumleft,'xdata',[]);
  set(gui.handles.axesHLAmyocardiumleft,'ydata',[]);
  %VLA
  set(gui.handles.axesVLAmyocardiumleft,'xdata',[]);
  set(gui.handles.axesVLAmyocardiumleft,'ydata',[]);
end

if gui.showno(2) ~= 0
  switch gui.showimagetype{2}
    case 'rest'
      if gui.LVsegmentationrest
        endoX = gui.EndoXrest;
        endoY = gui.EndoYrest;
        epiX = gui.EpiXrest;
        epiY = gui.EpiYrest;
        tf = 1;
        SAslices = gui.SAslicesrest;
        HLAmyocardiumx = gui.HLAmyocardiumxrest;
        HLAmyocardiumy = gui.HLAmyocardiumyrest;
        VLAmyocardiumx = gui.VLAmyocardiumxrest;
        VLAmyocardiumy = gui.VLAmyocardiumyrest;
      end
    case 'restgated'
      if gui.LVsegmentationrestgated
        endoX = gui.EndoXrestgated;
        endoY = gui.EndoYrestgated;
        epiX = gui.EpiXrestgated;
        epiY = gui.EpiYrestgated;
        tf = gui.handles.tf;
        for s = 1:3
          SAslices{s} = gui.SAslicesrestgated{s,tf};
        end
        HLAmyocardiumx = gui.HLAmyocardiumxrestgated;
        HLAmyocardiumy = gui.HLAmyocardiumyrestgated;
        VLAmyocardiumx = gui.VLAmyocardiumxrestgated;
        VLAmyocardiumy = gui.VLAmyocardiumyrestgated;
      end
  end
  if isempty(gui.handles.axesSAbasalendoright)
    %basal
    hold(gui.handles.axesSAbasalright,'on');
    gui.handles.axesSAbasalendoright = plot(gui.handles.axesSAbasalright,endoY(:,tf,round(mean(SAslices{3}))),endoX(:,tf,round(mean(SAslices{3}))),'w-');
    set(gui.handles.axesSAbasalendoright,'linewidth',2);
    gui.handles.axesSAbasalepiright = plot(gui.handles.axesSAbasalright,epiY(:,tf,round(mean(SAslices{3}))),epiX(:,tf,round(mean(SAslices{3}))),'w-');
    set(gui.handles.axesSAbasalepiright,'linewidth',2);
    hold(gui.handles.axesSAbasalright,'off');
    %mid
    hold(gui.handles.axesSAmidright,'on');
    gui.handles.axesSAmidendoright = plot(gui.handles.axesSAmidright,endoY(:,tf,round(mean(SAslices{2}))),endoX(:,tf,round(mean(SAslices{2}))),'w-');
    set(gui.handles.axesSAmidendoright,'linewidth',2);
    gui.handles.axesSAmidepiright = plot(gui.handles.axesSAmidright,epiY(:,tf,round(mean(SAslices{2}))),epiX(:,tf,round(mean(SAslices{2}))),'w-');
    set(gui.handles.axesSAmidepiright,'linewidth',2);
    hold(gui.handles.axesSAmidright,'off');
    %apical
    hold(gui.handles.axesSAapicalright,'on');
    gui.handles.axesSAapicalendoright = plot(gui.handles.axesSAapicalright,endoY(:,tf,round(mean(SAslices{1}))),endoX(:,tf,round(mean(SAslices{1}))),'w-');
    set(gui.handles.axesSAapicalendoright,'linewidth',2);
    gui.handles.axesSAapicalepiright = plot(gui.handles.axesSAapicalright,epiY(:,tf,round(mean(SAslices{1}))),epiX(:,tf,round(mean(SAslices{1}))),'w-');
    set(gui.handles.axesSAapicalepiright,'linewidth',2);
    hold(gui.handles.axesSAapicalright,'off');
    %HLA
    hold(gui.handles.axesHLAright,'on');
    gui.handles.axesHLAmyocardiumright = plot(gui.handles.axesHLAright,HLAmyocardiumy{tf},HLAmyocardiumx{tf},'w-');
    set(gui.handles.axesHLAmyocardiumright,'linewidth',2);
    hold(gui.handles.axesHLAright,'off');
    %VLA
    hold(gui.handles.axesVLAright,'on');
    gui.handles.axesVLAmyocardiumright = plot(gui.handles.axesVLAright,VLAmyocardiumy{tf},VLAmyocardiumx{tf},'w-');
    set(gui.handles.axesVLAmyocardiumright,'linewidth',2);
    hold(gui.handles.axesVLAright,'off');
  else
    %basal
    set(gui.handles.axesSAbasalendoright,'xdata',endoY(:,tf,round(mean(SAslices{3}))));
    set(gui.handles.axesSAbasalendoright,'ydata',endoX(:,tf,round(mean(SAslices{3}))));
    set(gui.handles.axesSAbasalepiright,'xdata',epiY(:,tf,round(mean(SAslices{3}))));
    set(gui.handles.axesSAbasalepiright,'ydata',epiX(:,tf,round(mean(SAslices{3}))));
    %mid
    set(gui.handles.axesSAmidendoright,'xdata',endoY(:,tf,round(mean(SAslices{2}))));
    set(gui.handles.axesSAmidendoright,'ydata',endoX(:,tf,round(mean(SAslices{2}))));
    set(gui.handles.axesSAmidepiright,'xdata',epiY(:,tf,round(mean(SAslices{2}))));
    set(gui.handles.axesSAmidepiright,'ydata',epiX(:,tf,round(mean(SAslices{2}))));
    %apical
    set(gui.handles.axesSAapicalendoright,'xdata',endoY(:,tf,round(mean(SAslices{1}))));
    set(gui.handles.axesSAapicalendoright,'ydata',endoX(:,tf,round(mean(SAslices{1}))));
    set(gui.handles.axesSAapicalepiright,'xdata',epiY(:,tf,round(mean(SAslices{1}))));
    set(gui.handles.axesSAapicalepiright,'ydata',epiX(:,tf,round(mean(SAslices{1}))));
    %HLA
    set(gui.handles.axesHLAmyocardiumright,'xdata',HLAmyocardiumy{tf});
    set(gui.handles.axesHLAmyocardiumright,'ydata',HLAmyocardiumx{tf});
    %VLA
    set(gui.handles.axesVLAmyocardiumright,'xdata',VLAmyocardiumy{tf});
    set(gui.handles.axesVLAmyocardiumright,'ydata',VLAmyocardiumx{tf});
  end
else
  %basal
  set(gui.handles.axesSAbasalendoright,'xdata',[]);
  set(gui.handles.axesSAbasalendoright,'ydata',[]);
  set(gui.handles.axesSAbasalepiright,'xdata',[]);
  set(gui.handles.axesSAbasalepiright,'ydata',[]);
  %mid
  set(gui.handles.axesSAmidendoright,'xdata',[]);
  set(gui.handles.axesSAmidendoright,'ydata',[]);
  set(gui.handles.axesSAmidepiright,'xdata',[]);
  set(gui.handles.axesSAmidepiright,'ydata',[]);
  %apical
  set(gui.handles.axesSAapicalendoright,'xdata',[]);
  set(gui.handles.axesSAapicalendoright,'ydata',[]);
  set(gui.handles.axesSAapicalepiright,'xdata',[]);
  set(gui.handles.axesSAapicalepiright,'ydata',[]);
  %HLA
  set(gui.handles.axesHLAmyocardiumright,'xdata',[]);
  set(gui.handles.axesHLAmyocardiumright,'ydata',[]);
  %VLA
  set(gui.handles.axesVLAmyocardiumright,'xdata',[]);
  set(gui.handles.axesVLAmyocardiumright,'ydata',[]);
end


%--------------------------
function hideLVsegmentation
%--------------------------
%hide the endo- and epicardial segmentation

global DATA

gui = DATA.GUI.SpectPlot2d;


if gui.showno(1) ~= 0
  %basal
  set(gui.handles.axesSAbasalendoleft,'xdata',[]);
  set(gui.handles.axesSAbasalendoleft,'ydata',[]);
  set(gui.handles.axesSAbasalepileft,'xdata',[]);
  set(gui.handles.axesSAbasalepileft,'ydata',[]);
  %mid
  set(gui.handles.axesSAmidendoleft,'xdata',[]);
  set(gui.handles.axesSAmidendoleft,'ydata',[]);
  set(gui.handles.axesSAmidepileft,'xdata',[]);
  set(gui.handles.axesSAmidepileft,'ydata',[]);
  %apical
  set(gui.handles.axesSAapicalendoleft,'xdata',[]);
  set(gui.handles.axesSAapicalendoleft,'ydata',[]);
  set(gui.handles.axesSAapicalepileft,'xdata',[]);
  set(gui.handles.axesSAapicalepileft,'ydata',[]);
  %HLA
  set(gui.handles.axesHLAmyocardiumleft,'xdata',[]);
  set(gui.handles.axesHLAmyocardiumleft,'ydata',[]);
  %VLA
  set(gui.handles.axesVLAmyocardiumleft,'xdata',[]);
  set(gui.handles.axesVLAmyocardiumleft,'ydata',[]);
end

if gui.showno(2) ~= 0
  %basal
  set(gui.handles.axesSAbasalendoright,'xdata',[]);
  set(gui.handles.axesSAbasalendoright,'ydata',[]);
  set(gui.handles.axesSAbasalepiright,'xdata',[]);
  set(gui.handles.axesSAbasalepiright,'ydata',[]);
  %mid
  set(gui.handles.axesSAmidendoright,'xdata',[]);
  set(gui.handles.axesSAmidendoright,'ydata',[]);
  set(gui.handles.axesSAmidepiright,'xdata',[]);
  set(gui.handles.axesSAmidepiright,'ydata',[]);
  %apical
  set(gui.handles.axesSAapicalendoright,'xdata',[]);
  set(gui.handles.axesSAapicalendoright,'ydata',[]);
  set(gui.handles.axesSAapicalepiright,'xdata',[]);
  set(gui.handles.axesSAapicalepiright,'ydata',[]);
  %HLA
  set(gui.handles.axesHLAmyocardiumright,'xdata',[]);
  set(gui.handles.axesHLAmyocardiumright,'ydata',[]);
  %VLA
  set(gui.handles.axesVLAmyocardiumright,'xdata',[]);
  set(gui.handles.axesVLAmyocardiumright,'ydata',[]);
end


%------------------------------
function hideahasections(panel)
%------------------------------
%hide the aha section division lines

global DATA

gui = DATA.GUI.SpectPlot2d;

if isequal(panel,'left') || isequal(panel,'both')
  if gui.showno(1) ~= 0 && isfield(gui.handles,'axesSAapicalahalineleft')
    %apical
    for k = 1:2
      set(gui.handles.axesSAapicalahalineleft(k),'xdata',0);
      set(gui.handles.axesSAapicalahalineleft(k),'ydata',0);
    end
    %mid and basal
    for k = 1:3
      set(gui.handles.axesSAmidahalineleft(k),'xdata',0);
      set(gui.handles.axesSAmidahalineleft(k),'ydata',0);
      set(gui.handles.axesSAbasalahalineleft(k),'xdata',0);
      set(gui.handles.axesSAbasalahalineleft(k),'ydata',0);
    end
    %HLA and VLA
    for k = 1:4
      set(gui.handles.axesHLAahalineleft(k),'xdata',0);
      set(gui.handles.axesHLAahalineleft(k),'ydata',0);
      set(gui.handles.axesVLAahalineleft(k),'xdata',0);
      set(gui.handles.axesVLAahalineleft(k),'ydata',0);
    end
  end
end
if isequal(panel,'right') | isequal(panel,'both')
  if gui.showno(2) ~= 0 && isfield(gui.handles,'axesSAapicalahalineright')
    %apical
    for k = 1:2
      set(gui.handles.axesSAapicalahalineright(k),'xdata',0);
      set(gui.handles.axesSAapicalahalineright(k),'ydata',0);
    end
    %mid and basal
    for k = 1:3
      set(gui.handles.axesSAmidahalineright(k),'xdata',0);
      set(gui.handles.axesSAmidahalineright(k),'ydata',0);
      set(gui.handles.axesSAbasalahalineright(k),'xdata',0);
      set(gui.handles.axesSAbasalahalineright(k),'ydata',0);
    end
    %HLA and VLA
    for k = 1:4
      set(gui.handles.axesHLAahalineright(k),'xdata',0);
      set(gui.handles.axesHLAahalineright(k),'ydata',0);
      set(gui.handles.axesVLAahalineright(k),'xdata',0);
      set(gui.handles.axesVLAahalineright(k),'ydata',0);
    end
  end
end


%----------------------------
function hideahavalues(panel)
%----------------------------
%hide the aha section division values

global DATA

gui = DATA.GUI.SpectPlot2d;

switch panel
  case 'both'
    editlabels = {gui.handles.editapicalantleft,gui.handles.editapicallatleft, ...
      gui.handles.editapicalinfleft,gui.handles.editapicalsepleft, ...
      gui.handles.editmidantleft,gui.handles.editmidantlatleft,gui.handles.editmidinflatleft, ...
      gui.handles.editmidinfleft,gui.handles.editmidinfsepleft,gui.handles.editmidantsepleft, ...
      gui.handles.editbasalantleft,gui.handles.editbasalantlatleft,gui.handles.editbasalinflatleft, ...
      gui.handles.editbasalinfleft,gui.handles.editbasalinfsepleft,gui.handles.editbasalantsepleft, ...
      gui.handles.editHLAleft, ...
      gui.handles.editapicalantright,gui.handles.editapicallatright, ...
      gui.handles.editapicalinfright,gui.handles.editapicalsepright, ...
      gui.handles.editmidantright,gui.handles.editmidantlatright,gui.handles.editmidinflatright, ...
      gui.handles.editmidinfright,gui.handles.editmidinfsepright,gui.handles.editmidantsepright, ...
      gui.handles.editbasalantright,gui.handles.editbasalantlatright,gui.handles.editbasalinflatright, ...
      gui.handles.editbasalinfright,gui.handles.editbasalinfsepright,gui.handles.editbasalantsepright, ...
      gui.handles.editHLAright};
  case 'left'
    editlabels = {gui.handles.editapicalantleft,gui.handles.editapicallatleft, ...
      gui.handles.editapicalinfleft,gui.handles.editapicalsepleft, ...
      gui.handles.editmidantleft,gui.handles.editmidantlatleft,gui.handles.editmidinflatleft, ...
      gui.handles.editmidinfleft,gui.handles.editmidinfsepleft,gui.handles.editmidantsepleft, ...
      gui.handles.editbasalantleft,gui.handles.editbasalantlatleft,gui.handles.editbasalinflatleft, ...
      gui.handles.editbasalinfleft,gui.handles.editbasalinfsepleft,gui.handles.editbasalantsepleft, ...
      gui.handles.editHLAleft};
  case 'right'
    editlabels = {gui.handles.editapicalantright,gui.handles.editapicallatright, ...
      gui.handles.editapicalinfright,gui.handles.editapicalsepright, ...
      gui.handles.editmidantright,gui.handles.editmidantlatright,gui.handles.editmidinflatright, ...
      gui.handles.editmidinfright,gui.handles.editmidinfsepright,gui.handles.editmidantsepright, ...
      gui.handles.editbasalantright,gui.handles.editbasalantlatright,gui.handles.editbasalinflatright, ...
      gui.handles.editbasalinfright,gui.handles.editbasalinfsepright,gui.handles.editbasalantsepright, ...
      gui.handles.editHLAright};
  otherwise
    editlabels = {};
end
for editloop = 1:length(editlabels)
  set(editlabels{editloop},'visible','off');
end


%----------------------------------
function hideSAintersections(panel)
%----------------------------------

global DATA

gui = DATA.GUI.SpectPlot2d;

if isequal(panel,'left') || isequal(panel,'both')
  if gui.showno(1) ~= 0
    for k = 1:length(gui.handles.axesHLAlinesleft)
      set(gui.handles.axesHLAlinesleft(k),'xdata',0);
      set(gui.handles.axesHLAlinesleft(k),'ydata',0);
      set(gui.handles.axesVLAlinesleft(k),'xdata',0);
      set(gui.handles.axesVLAlinesleft(k),'ydata',0);
    end
  end
end
if isequal(panel,'right') || isequal(panel,'both')
  if gui.showno(2) ~= 0
    for k = 1:length(gui.handles.axesHLAlinesright)
      set(gui.handles.axesHLAlinesright(k),'xdata',0);
      set(gui.handles.axesHLAlinesright(k),'ydata',0);
      set(gui.handles.axesVLAlinesright(k),'xdata',0);
      set(gui.handles.axesVLAlinesright(k),'ydata',0);
    end
  end
end

%----------------------------------
function showsegmentations_Callback %#ok<DEFNU>
%----------------------------------
global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.radiobuttonshowsegmentation,'value',1);
set(gui.handles.radiobuttonhidesegmentation,'value',0);
updateLVsegmentation;
if isequal(gui.mode,'scoring')
  if ~get(gui.handles.radiobuttongated,'value')
    updateahasections('both');
  else
    updateahasections('right');
    hideahasections('left');
  end
else
  updateSAintersections('both');
end


%---------------------------------
function hidesegmentations_Callback %#ok<DEFNU>
%----------------------------------
global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.radiobuttonshowsegmentation,'value',0);
set(gui.handles.radiobuttonhidesegmentation,'value',1);
hideLVsegmentation;
hideahasections('both');


%--------------------------
function showgated_Callback %#ok<DEFNU>
%--------------------------
global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.radiobuttongated,'value',1);
set(gui.handles.radiobuttonungated,'value',0);
set(gui.handles.pushbuttonnext,'enable','on');
set(gui.handles.pushbuttonprev,'enable','on');
set(gui.handles.pushbuttonplay,'enable','on');
set(gui.handles.pushbuttonstop,'enable','on');
%normalize
set(gui.handles.radiobuttonnormeachtf,'enable','on');
set(gui.handles.radiobuttonnormall,'enable','on');
%frame rate slider
set(gui.handles.sliderframerate,'enable','on');
set(gui.handles.sliderframerate,'value',gui.framerate);

switch gui.mode
  case 'visualization'
    if ~isempty(gui.restgatedno) && ~isempty(gui.stressgatedno)
      gui.showno = [gui.stressgatedno gui.restgatedno];
      gui.showimagetype = {'stressgated','restgated'};
    elseif ~isempty(gui.stressgatedno)
      gui.showno = [gui.stressgatedno 0];
      gui.showimagetype = {'stressgated', ''};
    else
      gui.showno = [0 gui.restgatedno];
      gui.showimagetype = {'', 'restgated'};
    end
    hideahavalues('both');
  case 'scoring'
    gui.showno = [gui.restgatedno gui.restno];
    gui.showimagetype = {'restgated','rest'};
    hideahavalues('left');
    plotahavalues('right');
    if isempty(gui.showno)
      gui.showno = [gui.stressgatedno 0];
      gui.showimagetype = {'stressgated',''};
      hideahavalues('both');
    end
end
if isequal(gui.mode,'visualization')
  switch gui.showimagetype{1}
    case 'stress'
      set(gui.handles.textvolumesleft,'String',gui.handles.StringStress);
    case 'stressgated'
      set(gui.handles.textvolumesleft,'String',gui.handles.StringStressgated);
    case 'restgated'
      set(gui.handles.textvolumesleft,'String',gui.handles.StringRestgated);
    otherwise      
      set(gui.handles.textvolumesleft,'String',[]);
  end
  switch gui.showimagetype{2}
    case 'rest'
      set(gui.handles.textvolumesright,'String',gui.handles.StringRest);
    case 'restgated'
      set(gui.handles.textvolumesright,'String',gui.handles.StringRestgated);
    otherwise      
      set(gui.handles.textvolumesright,'String',[]);
  end
end
updateimages;
if get(gui.handles.radiobuttonshowsegmentation,'value')
  updateLVsegmentation;
  if isequal(gui.mode,'scoring')
    updateahasections('right');
    hideahasections('left');
    hideSAintersections('right');
  else
    hideahasections('both');
  end
end
if isequal(gui.mode,'scoring')
  updateSAintersections('left');
else
  updateSAintersections('both');
end
plotvolumecurve;


%----------------------------
function showungated_Callback %#ok<DEFNU>
%----------------------------
global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.radiobuttongated,'value',0);
set(gui.handles.radiobuttonungated,'value',1);
set(gui.handles.pushbuttonnext,'enable','off');
set(gui.handles.pushbuttonprev,'enable','off');
set(gui.handles.pushbuttonplay,'enable','off');
set(gui.handles.pushbuttonstop,'enable','off');
set(gui.handles.radiobuttonnormeachtf,'enable','off');
set(gui.handles.radiobuttonnormall,'enable','off');
set(gui.handles.sliderframerate,'enable','off');

if ~isempty(gui.restno) && ~isempty(gui.stressno)
  gui.showno = [gui.stressno gui.restno];
  gui.showimagetype = {'stress', 'rest'};
elseif ~isempty(gui.stressno)
  gui.showno = [gui.stressno 0];
  gui.showimagetype = {'stress', ''};
else
  gui.showno = [0 gui.restno];
  gui.showimagetype = {'', 'rest'};
end
if isequal(gui.mode,'visualization')
  switch gui.showimagetype{1}
    case 'stress'
      set(gui.handles.textvolumesleft,'String',gui.handles.StringStress);
    case 'stressgated'
      set(gui.handles.textvolumesleft,'String',gui.handles.StringStressgated);
    case 'restgated'
      set(gui.handles.textvolumesleft,'String',gui.handles.StringRestgated);
    otherwise
      set(gui.handles.textvolumesrleft,'String',[]);
  end
  switch gui.showimagetype{2}
    case 'rest'
      set(gui.handles.textvolumesright,'String',gui.handles.StringRest);
    case 'restgated'
      set(gui.handles.textvolumesright,'String',gui.handles.StringRestgated);
    otherwise
      set(gui.handles.textvolumesright,'String',[]);
  end
end
updateimages;
if get(gui.handles.radiobuttonshowsegmentation,'value')
  updateLVsegmentation;
  if isequal(gui.mode,'scoring')
    updateahasections('both');   
    updateSAintersections('left');
  end
end
if isequal(gui.mode,'scoring')
  if ~isempty(gui.showimagetype{1}) && ~isempty(gui.showimagetype{2})
    plotahavalues('both');
  elseif ~isempty(gui.showimagetype{1}) 
    plotahavalues('left');
  elseif ~isempty(gui.showimagetype{2})
    plotahavalues('right');
  end
  hideSAintersections('both');
else
  hideahavalues('both');
  updateSAintersections('both');
end
plotvolumecurve;


%--------------------------------------
function setnormmaxcounteachtf_Callback %#ok<DEFNU>
%--------------------------------------
global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.radiobuttonnormeachtf,'value',1);
set(gui.handles.radiobuttonnormall,'value',0);
roundtf('closest');
trackupdate;


%-----------------------------------
function setnormmaxcountall_Callback %#ok<DEFNU>
%-----------------------------------
global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.radiobuttonnormeachtf,'value',0);
set(gui.handles.radiobuttonnormall,'value',1);
roundtf('closest');
trackupdate;


%---------------------
function next_Callback
%---------------------
%one time frame forwards

global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.pushbuttonplay,'value',0);

gui.handles.tf = gui.handles.tf+1;
if gui.handles.tf > gui.nbroftimeframes 
  gui.handles.tf = 1;
end
roundtf('up');
trackupdate;


%--------------------
function prev_Callback %#ok<DEFNU>
%--------------------
%one time frame backwards

global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.pushbuttonplay,'value',0);

gui.handles.tf = gui.handles.tf-1;
if gui.handles.tf < 1
  gui.handles.tf = gui.nbroftimeframes;
end
roundtf('down');
trackupdate;


%----------------------------
function roundtf(roundmethod)
%----------------------------
%round the time frame to the closest original time frames

global DATA

gui = DATA.GUI.SpectPlot2d;

switch roundmethod
  case 'down'
    originaltf = linspace(1,gui.nbroftimeframes,gui.bothtf);
    difftf = gui.handles.tf-originaltf;
    tfindex = find(difftf >= 0,1,'last');
    gui.handles.tf = round(originaltf(tfindex));
  case 'up'
    originaltf = linspace(1,gui.nbroftimeframes,gui.bothtf);
    difftf = gui.handles.tf-originaltf;
    tfindex = find(difftf <= 0,1,'first');
    gui.handles.tf = round(originaltf(tfindex));
  case 'closest'
    gui.handles.tf = interp1(linspace(1,gui.nbroftimeframes,gui.bothtf),linspace(1,gui.nbroftimeframes,gui.bothtf),gui.handles.tf,'nearest');
end

%--------------------
function play_Callback %#ok<DEFNU>
%--------------------
%play a movie

global DATA

gui = DATA.GUI.SpectPlot2d;

DATA.StartFrame = gui.handles.tf;
DATA.StartTime = now;
try
  while get(gui.handles.pushbuttonplay,'value')
    gui.handles.tf = gui.handles.tf+1;
    if gui.handles.tf > gui.nbroftimeframes
      gui.handles.tf = 1;
    end
    trackupdate;
    pause(gui.framerate/gui.nbroftimeframes);
  end;
catch %#ok<CTCH>
end


%-----------------------
function stop_Callback %#ok<DEFNU>
%-----------------------
global DATA

gui = DATA.GUI.SpectPlot2d;

set(gui.handles.pushbuttonplay,'value',0);

roundtf('up');
trackupdate;


%-------------------
function trackupdate
%-------------------
%update the tracking

global DATA

gui = DATA.GUI.SpectPlot2d;

switch gui.showimagetype{1}
  case 'stressgated'
    %normalize maximal count
    if get(gui.handles.radiobuttonnormeachtf,'value')
      normmaxcountstressgated = gui.maxcountstressgated(gui.handles.tf);
      normmincountstressgated = gui.mincountstressgated(gui.handles.tf);
    else
      normmaxcountstressgated = max(gui.maxcountstressgated);
      normmincountstressgated = min(gui.mincountstressgated);
    end
    %plot images
    set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagestressgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagestressgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagestressgated(:,:,gui.handles.tf));
    set(gui.handles.axesHLAimageleft,'cdata',gui.HLAimagestressgated(:,:,gui.handles.tf));
    set(gui.handles.axesVLAimageleft,'cdata',gui.VLAimagestressgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAbasalleft,'clim',[normmincountstressgated normmaxcountstressgated]);
    set(gui.handles.axesSAmidleft,'clim',[normmincountstressgated normmaxcountstressgated]);
    set(gui.handles.axesSAapicalleft,'clim',[normmincountstressgated normmaxcountstressgated]);
    set(gui.handles.axesHLAleft,'clim',[normmincountstressgated normmaxcountstressgated]);
    set(gui.handles.axesVLAleft,'clim',[normmincountstressgated normmaxcountstressgated]);
  case 'restgated'
    %normalize maximal count
    if get(gui.handles.radiobuttonnormeachtf,'value')
      normmaxcountrestgated = gui.maxcountrestgated(gui.handles.tf);
      normmincountrestgated = gui.mincountrestgated(gui.handles.tf);
    else
      normmaxcountrestgated = max(gui.maxcountrestgated);
      normmincountrestgated = min(gui.mincountrestgated);
    end
    %plot images
    set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesHLAimageleft,'cdata',gui.HLAimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesVLAimageleft,'cdata',gui.VLAimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAbasalleft,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesSAmidleft,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesSAapicalleft,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesHLAleft,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesVLAleft,'clim',[normmincountrestgated normmaxcountrestgated]);
end

switch gui.showimagetype{2}
  case 'restgated'
    %normalize maximal count
    if get(gui.handles.radiobuttonnormeachtf,'value')
      normmaxcountrestgated = gui.maxcountrestgated(gui.handles.tf);
      normmincountrestgated = gui.mincountrestgated(gui.handles.tf);
    else
      normmaxcountrestgated = max(gui.maxcountrestgated);
      normmincountrestgated = min(gui.mincountrestgated);
    end
    %plot images
    set(gui.handles.axesSAbasalimageright,'cdata',gui.SAbasalimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAmidimageright,'cdata',gui.SAmidimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAapicalimageright,'cdata',gui.SAapicalimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesHLAimageright,'cdata',gui.HLAimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesVLAimageright,'cdata',gui.VLAimagerestgated(:,:,gui.handles.tf));
    set(gui.handles.axesSAbasalright,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesSAmidright,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesSAapicalright,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesHLAright,'clim',[normmincountrestgated normmaxcountrestgated]);
    set(gui.handles.axesVLAright,'clim',[normmincountrestgated normmaxcountrestgated]);
end

%set time frame
set(gui.handles.timebar,'xdata',[gui.handles.tf gui.handles.tf]);
set(gui.handles.timebar,'ydata',get(gui.handles.volumeaxes,'ylim'));
if get(gui.handles.radiobuttonshowsegmentation,'value');
  updateLVsegmentation; 
end
switch gui.mode
  case 'visualization'
    updateSAintersections('both');
  case 'scoring'
    updateSAintersections('left');
end

try
  if DATA.Record
    drawnow;
    DATA.MovieFrame = mygetframe(gui.fig);
    export('exportmovierecorder_Callback','newframe');
  end
catch %#ok<CTCH>
end


%-----------------------------
function setframerate_Callback
%-----------------------------
%update the Beat Time

global DATA

gui = DATA.GUI.SpectPlot2d;

framerate = get(gui.handles.sliderframerate,'value');
%rescale from [0 1] to [0.1 0.9] and reverse
framerate = (1-framerate*0.8+0.1);
gui.framerate = framerate;


%--------------------------------
function succeed = close_Callback
%--------------------------------
%close the gui

global DATA SET

gui = DATA.GUI.SpectPlot2d;

if get(gui.handles.radiobuttongated,'value')
  stop_Callback;
end

if isequal(gui.mode,'scoring')
  for k = 1:size(gui.scoringvalues,2)
    if isempty(gui.scoringvalues{1,k}) || ~isnumeric(gui.scoringvalues{1,k}) || gui.scoringvalues{1,k}<0 || gui.scoringvalues{1,k}>4
      myfailed('REST: Scoring value missing, or outside the allowed range [0 4]');
      succeed = false;
      return;
    end
    if isempty(gui.scoringvalues{2,k}) || ~isnumeric(gui.scoringvalues{2,k}) || gui.scoringvalues{2,k}<0 || gui.scoringvalues{2,k}>4
      myfailed('STRESS: Scoring value missing, or outside the allowed range [0 4]');
      succeed = false;
      return;
    end
    if isempty(gui.scoringvalues{3,k}) || ~isnumeric(gui.scoringvalues{3,k}) || gui.scoringvalues{3,k}<0 || gui.scoringvalues{3,k}>4
      myfailed('REST Infarct: Scoring value missing, or outside the allowed range [0 4]');
      succeed = false;
      return;
    end
  end
  if ~isempty(gui.restno) && isempty(SET(gui.restno).Perfusion)
    perfusion('initdefault',gui.restno);
    perfusion('createmyocardmask',gui.restno);
  end
  if ~isempty(gui.stressno) && isempty(SET(gui.stressno).Perfusion)
    perfusion('initdefault',gui.stressno);
    perfusion('createmyocardmask',gui.stressno);
  end
  for k = 1:length(gui.scoringvalues)
    if ~isempty(gui.stressno)
      SET(gui.stressno).Perfusion.MPS.ScoringIschManual(k) = gui.scoringvalues{1,k};
    end
    if ~isempty(gui.restno)
      SET(gui.restno).Perfusion.MPS.ScoringIschManual(k) = gui.scoringvalues{2,k};
      SET(gui.restno).Perfusion.MPS.ScoringFixManual(k) = gui.scoringvalues{3,k};
    end
  end
end
succeed = true;
mode = gui.mode;

try
  DATA.GUI.SpectPlot2d = close(DATA.GUI.SpectPlot2d);
catch %#ok<CTCH>
  close(gcf)
end

% if isequal(mode,'scoring')
DATA.Buffer.KeyStroke = {'ok'};
filemenu('saveall_Callback'); %save the file
% end


%---------------------------------------------
function newim = upsampleimagespatial(f1,f2,im)
%---------------------------------------------
%Helper function to upsample an image spatially
%input
%f1:    upsample factor in x-direction
%f2:    upsample factor in y-direction
%im:    image to upsample
%output
%newim: upsampled image

%Find new size
% newsize = size(imresize(im(:,:,1),f,'nearest'));
newsize = size(imresize(im(:,:,1),[round(size(im,1)*f1) round(size(im,2)*f2)],'nearest'));

%Reserve memory
if isa(im,'single')
  newim = repmat(single(0),[newsize(1) newsize(2) size(im,3) ]);
elseif isa(im,'int16')
    newim = repmat(int16(0),[newsize(1) newsize(2) size(im,3)]);
else
  newim = zeros(newsize(1),newsize(2),size(im,3));
end

%Loop over image volume
for tloop = 1:size(im,3)
  if isa(im,'single')
    newim(:,:,tloop) = single(imresize(im(:,:,tloop),[round(size(im,1)*f1) round(size(im,2)*f2)],'bicubic'));
  elseif isa(im,'int16')
    newim(:,:,tloop) = int16(imresize(im(:,:,tloop),[round(size(im,1)*f1) round(size(im,2)*f2)],'bicubic'));
  else
    newim(:,:,tloop) = imresize(im(:,:,tloop),[round(size(im,1)*f1) round(size(im,2)*f2)],'bicubic');
  end
end

%-------------------------------------------
function newim = upsampleimagetemporal(f,im)
%-------------------------------------------
%Helper function to upsample an image temporally.

N = size(im,3);

tnew = linspace(0,N,round(f*N)+1); 
tnew(end) = [];

newim = zeros(size(im,1),size(im,2),length(tnew));
for i=1:size(newim,3)
  g = 0;
  kfac = 6;
  for j=-kfac*size(im,3):((kfac+1)*size(im,3)+1)
    tj = j-1;
    if abs(tnew(i)-tj)<sqrt(eps)
      q = 1;
    else
      q = sin(pi*(tnew(i)-tj))/(pi*(tnew(i)-tj));
    end
    jmod = mod(j-1,size(im,3))+1;
    g = g + q*im(:,:,jmod);
  end
  newim(:,:,i) = g;
end

if not(isempty(im))
  newim(:,:,1,:) = im(:,:,1);
else
  newim = [];
end

%--------------------------------------------------------------------------
function [newline,isLVseg] = adjustsegmentation(fs,fsz,ft,line,startslices,endslices,isLVseg)
%--------------------------------------------------------------------------
%adjust segmentation after upsampling line is either EndoX, EndoY, EpiX or
%EpiY

% warnstate = warning('off'); %#ok<WNOFF>

%adjust for upsampling in X/Y direction
if fs>0
  d = (ceil(fs)-1)/2;
else
  d = 0;
end
newline = line*fs-d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust for upsmpling in Z
%Find new size
temp=squeeze(newline(1,1,:));

xi = linspace(1,length(temp),round(length(temp)*fsz));
newsize = length(xi);

% %Reserve memory
% templine1 = zeros(size(newline,1),newsize,size(newline,2));
% %Loop over image volume
% for xloop=1:size(newline,1)
%   templine1(xloop,:,:) = interp1(squeeze(newline(xloop,:,:))',xi,'linear');
% end;
% templine1 = permute(templine1,[1 3 2]);
% % newline = templine;

startslices = startslices*fsz-d;
endslices = endslices*fsz-d;

templine = zeros(size(newline,1),size(newline,2),newsize);
for xloop = 1:size(newline,1)
  for tloop = 1:size(newline,2)
    temp = squeeze(newline(xloop,tloop,:));
    if all(isnan(temp))
      templine(xloop,tloop,:) = NaN(1,newsize);
    elseif not(any(isnan(temp)))
      templine(xloop,tloop,:) = interp1(temp,xi,'linear');
    elseif sum(~isnan(temp)) == 1
      index = find(not(isnan(temp)));
      isvalid = ((1:size(templine,3))>=round(startslices(tloop)) & (1:size(templine,3))<=round(endslices(tloop)));
      templine(xloop,tloop,isvalid) = repmat(temp(index),[1 sum(isvalid)]);
      templine(xloop,tloop,not(isvalid)) = nan(1,sum(not(isvalid)));
    else
      index = find(not(isnan(temp)));
      isvalid = ((1:size(templine,3))>=round(startslices(tloop)) & (1:size(templine,3))<=round(endslices(tloop)));
      templine(xloop,tloop,isvalid) = interp1(index,temp(index),xi(isvalid),'linear');
      templine(xloop,tloop,not(isvalid)) = nan(1,sum(not(isvalid)));
    end
  end
end;
newline = templine;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust for upsampling in T
if ft > 1
  %Find new size
  temp = squeeze(newline(1,:,1));

  xi = linspace(1,length(temp),round(length(temp)*ft));
  newsize = length(xi);

  %Reserve memory
  templine = zeros(size(newline,1),newsize,size(newline,3));

  %first interpolate startslices and endslices
  if length(isLVseg) == length(startslices)
    isLVseg = 1:newsize;
    startslices = interp1(startslices,xi,'linear');
    endslices = interp1(endslices,xi,'linear');
    startslices = round(startslices);
    endslices = round(endslices);
  else
    oldisLVseg = isLVseg;
    oldstartslices = startslices;
    oldendslices = endslices;
    startslices = nan(1,newsize);
    endslices = startslices;
    isLVseg = round(oldisLVseg*(newsize+1)/length(oldstartslices)-1);
    startslices(isLVseg) = oldstartslices(oldisLVseg);
    endslices(isLVseg) = oldendslices(oldisLVseg);
  end  

  %Loop over image volume
  for xloop=1:size(newline,1)
    for zloop=1:size(newline,3)
      temp=squeeze(newline(xloop,:,zloop));
      if all(isnan(temp))
        templine(xloop,:,zloop) = NaN(1,newsize);
      elseif not(any(isnan(temp)))
        templine(xloop,:,zloop) = interp1(temp,xi,'linear');
      elseif sum(~isnan(temp)) == 1
        index=find(not(isnan(temp)));
        isvalid=(zloop>=startslices&zloop<=endslices);
        templine(xloop,isvalid,zloop) = repmat(temp(index),[1 sum(isvalid)]);
        templine(xloop,not(isvalid),zloop) = nan(1,sum(not(isvalid)));
      else
        index=find(not(isnan(temp)));
        isvalid=(zloop>=startslices&zloop<=endslices);
        templine(xloop,isvalid,zloop) = interp1(index,temp(index),xi(isvalid),'linear');
        templine(xloop,not(isvalid),zloop) = nan(1,sum(not(isvalid)));
      end
    end
  end;
  newline = templine;
end

%--------------------------------------------------------------------------
function [startslices,endslices,SAslices,HLAslice,VLAslice,cropstartslice,cropendslice] ...
  = calcslicesfromsegmentation(epix,epiy,isLVseg,mode,cropLAslices,zsize) %,SAslices,fsz)
%--------------------------------------------------------------------------

TSize = size(epix,2);
%finding slices with segmentation
startslices = nan(1,TSize);
endslices = nan(1,TSize);
for tf = isLVseg
  startslices(tf) = find(~isnan(epix(1,tf,:)),1,'first');  %the most basal slice
  endslices(tf) = find(~isnan(epix(1,tf,:)),1,'last');  %the most apical slice
end

if nargin > 4
  cropendslice = min(zsize,max(endslices)+cropLAslices);
  cropstartslice = max(1,min(startslices)-cropLAslices);
else
  cropendslice = [];
  cropstartslice = [];
end
nbrofslices = endslices-startslices+1;
switch mode
  case 'visualization'
    for tf = isLVseg
      SAslices{1,tf} = startslices(tf)+round(0.75*nbrofslices(tf));
      SAslices{2,tf} = startslices(tf)+round(0.50*nbrofslices(tf));
      SAslices{3,tf} = startslices(tf)+round(0.25*nbrofslices(tf));
    end
  case 'scoring'
    for tf = isLVseg
      SAslices{1,tf} = startslices(tf)+round(0.6333*nbrofslices(tf))+1:startslices(tf)+round(0.85*nbrofslices(tf));
      SAslices{2,tf} = startslices(tf)+round(0.4167*nbrofslices(tf))+1:startslices(tf)+round(0.6333*nbrofslices(tf));
      SAslices{3,tf} = startslices(tf)+round(0.2000*nbrofslices(tf)):startslices(tf)+round(0.4167*nbrofslices(tf));
    end
end

tempepix = epix(:,1:TSize,min(startslices):max(endslices));
tempepiy = epiy(:,1:TSize,min(startslices):max(endslices));
HLAslice = round(mynanmean(tempepix(:)));
VLAslice = round(mynanmean(tempepiy(:)));


%--------------------------------------------------------------------------
function [startslices,endslices,SAslices,HLAslice,VLAslice,startslice,endslice] = calcslicesfromstack(IM,no)  
%--------------------------------------------------------------------------
global SET

if SET(no).EndSlice > SET(no).StartSlice+3
  startslices = SET(no).StartSlice*ones(1,SET(no).TSize);
  endslices = SET(no).EndSlice*ones(1,SET(no).TSize);
else
  disp('Less than 4 slices marked, using all slices instead');
  startslices = ones(1,SET(no).TSize);
  endslices = SET(no).ZSize*ones(1,SET(no).TSize);
end
nbrofslices = endslices-startslices+1;

SAslices{1} = startslices(1)+round(0.75*nbrofslices(1));
SAslices{2} = startslices(1)+round(0.50*nbrofslices(1));
SAslices{3} = startslices(1)+round(0.25*nbrofslices(1));

HLAslice = round(size(IM,1)/2);
VLAslice = round(size(IM,2)/2);

startslice = 1;
endslice = size(IM,4);


%-----------------------
function plotvolumecurve
%-----------------------
%Calculate volume of segmentation and graphically update.
global DATA SET

gui = DATA.GUI.SpectPlot2d;
plottedright = 0;
plottedleft = 0;
yminvalue = [];
ymaxvalue = [];
cla(gui.handles.volumeaxes);

if get(gui.handles.radiobuttongated,'value')  
  switch gui.showimagetype{1}
    case 'stressgated'
      %calculate LVV for the interpolated time frames
      if length(gui.isLVsegstressgated) == gui.nbroftimeframes
        t = 1:gui.nbroftimeframes;
        lvvleft = calclvv(gui.startslicesstressgated,gui.endslicesstressgated,gui.EndoXstressgated,gui.EndoYstressgated,t,gui.resamplexystressgated,gui.resamplezstressgated);
        gui.handles.volumecurveleft = plot(gui.handles.volumeaxes,t,lvvleft,'b-');
        tpointsleft = linspace(1,gui.nbroftimeframes,SET(gui.stressgatedno).TSize);
        lvvpointsleft = interp1(1:gui.nbroftimeframes,lvvleft,tpointsleft);
        EDTleft = gui.EDTstressgated;
        ESTleft = gui.ESTstressgated;
      else
        t = gui.isLVsegstressgated;
        lvvleft = calclvv(gui.startslicesstressgated,gui.endslicesstressgated,gui.EndoXstressgated,gui.EndoYstressgated,t,gui.resamplexystressgated,gui.resamplezstressgated);
        lvvpointsleft = lvvsleft;
        tpointsleft = t;
        gui.handles.volumecurveleft = [];
        EDTleft = gui.EDTstressgated;
        ESTleft = gui.ESTstressgated;
      end
      textleft = 'stress';
    case 'restgated'
      %calculate LVV for the interpolated time frames
      if length(gui.isLVsegrestgated) == gui.nbroftimeframes
        t = 1:gui.nbroftimeframes;
        lvvleft = calclvv(gui.startslicesrestgated,gui.endslicesrestgated,gui.EndoXrestgated,gui.EndoYrestgated,t,gui.resamplexyrestgated,gui.resamplezrestgated);
        gui.handles.volumecurveleft = plot(gui.handles.volumeaxes,t,lvvleft,'b-');
        tpointsleft = linspace(1,gui.nbroftimeframes,SET(gui.restgatedno).TSize);
        lvvpointsleft= interp1(1:gui.nbroftimeframes,lvvleft,tpointsleft);
        EDTleft = gui.EDTrestgated;
        ESTleft = gui.ESTrestgated;
      else
        t = gui.isLVsegrestgated;
        lvvleft = calclvv(gui.startslicesrestgated,gui.endslicesrestgated,gui.EndoXrestgated,gui.EndoYrestgated,t,gui.resamplexyrestgated,gui.resamplezrestgated);
        lvvpointsleft = lvvsleft;
        tpointsleft = t;
        gui.handles.volumecurveleft = [];
        EDTleft = gui.EDTrestgated;
        ESTleft = gui.ESTrestgated;
      end
      textleft = 'rest';
  end
  switch gui.showimagetype{1}
    case {'stressgated','restgated'}
      hold(gui.handles.volumeaxes,'on');
      gui.handles.volumepointsleft = plot(gui.handles.volumeaxes,tpointsleft,lvvpointsleft,'b.');
      set(gui.handles.volumepointsleft,'markersize',5);
      set([gui.handles.volumecurveleft gui.handles.volumepointsleft],'ButtonDownFcn','spect.spectplot2d(''volumeaxes_Buttondown'')');
      yminvalue = [yminvalue ceil(min(lvvleft))];
      ymaxvalue = [ymaxvalue ceil(max(lvvleft))];
      plottedleft = 1;
    otherwise
      set(gui.handles.textlegendleft,'string',[]);
  end
  
  switch gui.showimagetype{2}
    case 'restgated'
      %calculate LVV for the interpolated time frames
      if length(gui.isLVsegrestgated) == gui.nbroftimeframes
        t = 1:gui.nbroftimeframes;
        lvvright = calclvv(gui.startslicesrestgated,gui.endslicesrestgated,gui.EndoXrestgated,gui.EndoYrestgated,t,gui.resamplexyrestgated,gui.resamplezrestgated);
        gui.handles.volumecurveright = plot(gui.handles.volumeaxes,t,lvvright,'r-');
        tpointsright = linspace(1,gui.nbroftimeframes,SET(gui.restgatedno).TSize);
        lvvpointsright= interp1(1:gui.nbroftimeframes,lvvright,tpointsright);
        EDTright = gui.EDTrestgated;
        ESTright = gui.ESTrestgated;
      else
        t = gui.isLVsegrestgated;
        lvvright = calclvv(gui.startslicesrestgated,gui.endslicesrestgated,gui.EndoXrestgated,gui.EndoYrestgated,t,gui.resamplexyrestgated,gui.resamplezrestgated);
        lvvpointsright = lvvsright;
        tpointsright = t;
        gui.handles.volumecurveright = [];
        EDTright = gui.EDTrestgated;
        ESTright = gui.ESTrestgated;
      end
      hold(gui.handles.volumeaxes,'on');
      gui.handles.volumepointsright = plot(gui.handles.volumeaxes,tpointsright,lvvpointsright,'r.');
      set(gui.handles.volumepointsright,'markersize',5);
      set([gui.handles.volumecurveright gui.handles.volumepointsright],'ButtonDownFcn','spect.spectplot2d(''volumeaxes_Buttondown'')');
      yminvalue = [yminvalue ceil(min(lvvright))];
      ymaxvalue = [ymaxvalue ceil(max(lvvright))];
      plottedright = 1;
      textright = 'rest';
    otherwise
      set(gui.handles.textlegendright,'string',[]);
  end
  
  %Time resolved  
  hold(gui.handles.volumeaxes,'on');  
 
  if isempty(ymaxvalue)
    %no LV segmentation
    ymaxvalue = 1;
    yminvalue = 0;
  end
  ymax = max(ymaxvalue);
  ymin = min(yminvalue);
  yaxislength = ymax-ymin;
  if yaxislength == 0
    ymax = ymax+0.1*ymax;
    ymin = ymin-0.1*ymin;
  end
  yaxislength = ymax-ymin;
  ymax = ymax+0.4*yaxislength;
  ymin = ymin-0.4*yaxislength;
  lc = [0.5 0.5 0];
  if plottedright
    gui.handles.edtextright = text(...
      'position',[EDTright ymin+0.2*yaxislength],...
      'string',sprintf('%s\n%s','ED',textright),...
      'parent',gui.handles.volumeaxes,...
      'color',lc);
    gui.handles.edlineright = plot(gui.handles.volumeaxes,...
      [EDTright EDTright],[ymin-0.2*(ymax-ymin) ymin+0.2*(ymax-ymin)],'color',lc,'linewidth',2);
    gui.handles.estextright = text(...
      'parent',gui.handles.volumeaxes,...
      'position',[ESTright ymin+0.2*yaxislength],...
      'string',sprintf('%s\n%s','ES',textright),...
      'color',lc);
    gui.handles.eslineright = plot(gui.handles.volumeaxes,...
      [ESTright ESTright],[ymin-0.2*(ymax-ymin) ymin+0.2*(ymax-ymin)],'color',lc,'linewidth',2);
    set([gui.handles.edtextright gui.handles.edlineright gui.handles.estextright gui.handles.eslineright], ...
    'ButtonDownFcn','spect.spectplot2d(''volumeaxes_Buttondown'')');
  end
  if plottedleft
    gui.handles.edtextleft = text(...
      'position',[EDTleft ymax-0.2*yaxislength],...
      'string',sprintf('%s\n%s','ED',textleft),...
      'parent',gui.handles.volumeaxes,...
      'color',lc);
    gui.handles.edlineleft = plot(gui.handles.volumeaxes,...
      [EDTleft EDTleft],[ymax+0.2*(ymax-ymin) ymax-0.2*(ymax-ymin)],'color',lc,'linewidth',2);
    gui.handles.estextleft = text(...
      'parent',gui.handles.volumeaxes,...
      'position',[ESTleft ymax-0.2*yaxislength],...
      'string',sprintf('%s\n%s','ES',textleft),...
      'color',lc);
    gui.handles.eslineleft = plot(gui.handles.volumeaxes,...
      [ESTleft ESTleft],[ymax+0.2*(ymax-ymin) ymax-0.2*(ymax-ymin)],'color',lc,'linewidth',2);
    set([gui.handles.edtextleft gui.handles.edlineleft gui.handles.estextleft gui.handles.eslineleft], ...
    'ButtonDownFcn','spect.spectplot2d(''volumeaxes_Buttondown'')');
  end
  
  %current time frame
  grid(gui.handles.volumeaxes,'on');
  gui.handles.timebar = plot(gui.handles.volumeaxes,...
    [gui.handles.tf gui.handles.tf],[ymin ymax],'k-');
  set(gui.handles.timebar,'linewidth',2);  
  set(gui.handles.volumeaxes,'ButtonDownFcn','spect.spectplot2d(''volumeaxes_Buttondown'')');  %mouse click
  set(gui.fig,'keypressfcn',@keypressed);  %arrow keys pressed
  hold(gui.handles.volumeaxes,'off');
  %X-axis properties
  set(gui.handles.volumeaxes,'xlim',[1-0.01*gui.nbroftimeframes gui.nbroftimeframes]);
  if isequal(gui.showimagetype{1},'stressgated') && isequal(gui.showimagetype{2},'stressrest')
    if SET(gui.showno(1)).TSize == SET(gui.showno(2)).TSize
      set(gui.handles.volumeaxes,'xtick',linspace(1,gui.nbroftimeframes,SET(gui.showno(2)).TSize));
      set(gui.handles.volumeaxes,'XTickLabel',1:SET(gui.showno(2)).TSize);
    else
      set(gui.handles.volumeaxes,'xtick',[]);
    end
  elseif isequal(gui.showimagetype{1},'stressgated')
    set(gui.handles.volumeaxes,'xtick',linspace(1,gui.nbroftimeframes,SET(gui.stressgatedno).TSize));
    set(gui.handles.volumeaxes,'XTickLabel',1:SET(gui.stressgatedno).TSize);
  elseif isequal(gui.showimagetype{1},'restgated') || isequal(gui.showimagetype{2},'restgated')
    set(gui.handles.volumeaxes,'xtick',linspace(1,gui.nbroftimeframes,SET(gui.restgatedno).TSize));
    set(gui.handles.volumeaxes,'XTickLabel',1:SET(gui.restgatedno).TSize);
  else
    set(gui.handles.volumeaxes,'xtick',[]);
  end
  xlabel(gui.handles.volumeaxes,'Time frame','color',DATA.GUISettings.VolumeAxesColor);
  %Y-axis properties
  if plottedleft || plottedright
    set(gui.handles.volumeaxes,'ylim',[ymin ymax]);
    ylabel(gui.handles.volumeaxes,'Volume [ml]','color',DATA.GUISettings.VolumeAxesColor);
    %Plot properties
    title(gui.handles.volumeaxes,'Left ventricular blood volume','color',DATA.GUISettings.VolumeAxesColor);
  else
    set(gui.handles.volumeaxes,'YTickLabel','');
    ylabel(gui.handles.volumeaxes,'','color',DATA.GUISettings.VolumeAxesColor);
    title(gui.handles.volumeaxes,'','color',DATA.GUISettings.VolumeAxesColor);
  end
  set(gui.handles.volumeaxes,...
    'XColor',DATA.GUISettings.VolumeAxesColor,...
    'YColor',DATA.GUISettings.VolumeAxesColor);
  set(gui.handles.volumeaxes, ...
    'Color',DATA.GUISettings.VolumeColorGraph);
  %Legend properties
  if plottedright
    legendtext = sprintf('%s%s','- ',textright);
    set(gui.handles.textlegendright,'string',legendtext,'foregroundcolor','r');
  end
  if plottedleft
    legendtext = sprintf('%s%s','- ',textleft);
    set(gui.handles.textlegendleft,'string',legendtext,'foregroundcolor','b');
  end
else
  cla(gui.handles.volumeaxes);
  gui.handles.volumecurveright = imagesc([],'parent',gui.handles.volumeaxes);
  axis(gui.handles.volumeaxes,'off');
  gui.handles.timebar = [];
  set(gui.handles.textlegendright,'string',[]);
  set(gui.handles.textlegendleft,'string',[]);
end


%-----------------------------
function volumeaxes_Buttondown %#ok<DEFNU>
%-----------------------------
%update the time frame after clicking in the volume plot
%
%written by Helen Soneson 2010-06-07

global DATA

gui = DATA.GUI.SpectPlot2d;

p = get(gca,'CurrentPoint');
t = round(p(1));
t = max([min([t gui.nbroftimeframes]) 1]);
gui.handles.tf = t;

%update image and segmentations
roundtf('closest');
trackupdate;


%----------------------------------
function moveSAline(movement,image)
%---------------------------------
%move the SA line in the LA images

global DATA

gui = DATA.GUI.SpectPlot2d;
if nargin < 2
  image = gui.LAline.imagestart;
end

if isequal(gui.mode,'scoring')
  if get(gui.handles.radiobuttonshowsegmentation,'value')
    %catch mouse click
    type = get(gui.fig,'SelectionType');
    switch type
      case 'normal'  %left button on the mouse
        switch image
          case 'HLAright'
            h = gui.handles.axesHLAright;
          case 'VLAright'
            h = gui.handles.axesVLAright;
          case 'HLAleft'
            h = gui.handles.axesHLAleft;
          case 'VLAleft'
            h = gui.handles.axesVLAleft;
        end
        switch movement
          case 'down'
            gui.LAline.imagestart = image;
            set(gui.fig,'WindowButtonMotionFcn','spect.spectplot2d(''moveSAline'',''motion'')');
            set(gui.fig,'WindowButtonUpFcn','spect.spectplot2d(''moveSAline'',''up'')');
            [x,y] = mygetcurrentpoint(gca);
            xlim = get(h,'xlim');
            ylim = get(h,'ylim');
            if y>ylim(2) || y<ylim(1) || ...
                x>xlim(2) || x<xlim(1)
              myfailed('The point is outside the image dimension',gui);
              return;
            end
            gui.LAline.ystart = round(y);
            gui.LAline.xstart = round(x);
          case 'motion'
            [x,y] = mygetcurrentpoint(gca);
          case 'up'
            [x,y] = mygetcurrentpoint(gca);
            xlim = get(h,'xlim');
            ylim = get(h,'ylim');
            if y>ylim(2) || y<ylim(1) || ...
                x>xlim(2) || x<xlim(1)
              myfailed('The point is outside the image dimension',gui);
              return;
            end
            gui.LAline.yend = round(y);
            gui.LAline.xend = round(x);
            gui.LAline.imageend = image;
            set(gui.fig,'WindowButtonMotionFcn','');
            set(gui.fig,'WindowButtonUpFcn','');
            if ~isempty(gui.LAline.xstart) && ~isempty(gui.LAline.ystart) && ...
                ~isempty(gui.LAline.xend) && ~isempty(gui.LAline.yend) %&& ...
              %isequal(gui.LAline.imagestart,gui.LAline.imageend)
              %find the corresponding line (closest to the start point)
              %update the SA intersection lines in LA images
              %recalculate and update the affected SA images and segmentations
              switch image
                case 'HLAright'
                  ahalines = get(gui.handles.axesHLAahalineright,'ydata');
                  ylim = get(gui.handles.axesHLAright,'ylim');
                  for k = 1:length(ahalines)
                    dist2line(k) = gui.LAline.ystart-ahalines{k}(1);
                  end
                  [val,mindist2lineindex] = min(abs(dist2line));
                  switch mindist2lineindex
                    case 1 %apical
                      if gui.SAslicesrest{1}(1)>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(1);
                      elseif gui.endslicesrest<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(1):gui.endslicesrest;
                      else
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(1):(ylim(2)-gui.LAline.yend);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(1),'ydata',(ylim(2)-gui.SAslicesrest{1}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(1),'xdata',gui.SAslicesrest{1}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{1})));
                      set(gui.handles.axesSAapicalimageright,'cdata',gui.SAapicalimagerest(:,:,1));
                    case 2 %apical and mid
                      if gui.SAslicesrest{2}(1)>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{1} = gui.SAslicesrest{2}(1)+1:gui.SAslicesrest{1}(end);
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(1);
                      elseif gui.SAslicesrest{1}(end)<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(end);
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(1):gui.SAslicesrest{1}(end)+1;
                      else
                        gui.SAslicesrest{1} = (ylim(2)-gui.LAline.yend)-1:gui.SAslicesrest{1}(end);
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(1):(ylim(2)-gui.LAline.yend);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(2),'ydata',(ylim(2)-gui.SAslicesrest{2}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(2),'xdata',gui.SAslicesrest{2}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{1})));
                      set(gui.handles.axesSAapicalimageright,'cdata',gui.SAapicalimagerest(:,:,1));
                      gui.SAmidimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{2})));
                      set(gui.handles.axesSAmidimageright,'cdata',gui.SAmidimagerest(:,:,1));
                    case 3 %mid and basal
                      if gui.SAslicesrest{3}(1)>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{2} = gui.SAslicesrest{3}(1)+1:gui.SAslicesrest{2}(end);
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(1);
                      elseif gui.SAslicesrest{2}(end)<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(end);
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(1):gui.SAslicesrest{2}(end)+1;
                      else
                        gui.SAslicesrest{2} = (ylim(2)-gui.LAline.yend)-1:gui.SAslicesrest{2}(end);
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(1):(ylim(2)-gui.LAline.yend);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(3),'ydata',(ylim(2)-gui.SAslicesrest{3}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(3),'xdata',gui.SAslicesrest{3}(end)*[1 1]);
                      %SA
                      gui.SAmidimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{2})));
                      set(gui.handles.axesSAmidimageright,'cdata',gui.SAmidimagerest(:,:,1));
                      gui.SAbasalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{3})));
                      set(gui.handles.axesSAbasalimageright,'cdata',gui.SAbasalimagerest(:,:,1));
                    case 4 %basal
                      if gui.startslicesrest>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{3} = gui.startslicesrest:gui.SAslicesrest{3}(end);
                      elseif gui.SAslicesrest{3}(end)<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(end);
                      else
                        gui.SAslicesrest{3} = (ylim(2)-gui.LAline.yend)-1:gui.SAslicesrest{3}(end);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(4),'ydata',(ylim(2)-gui.SAslicesrest{3}(1))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(4),'xdata',gui.SAslicesrest{3}(1)*[1 1]);
                      %SA
                      gui.SAbasalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{3})));
                      set(gui.handles.axesSAbasalimageright,'cdata',gui.SAbasalimagerest(:,:,1));
                  end
                case 'VLAright'
                  ahalines = get(gui.handles.axesVLAahalineright,'xdata');
                  ylim = get(gui.handles.axesHLAright,'ylim');
                  for k = 1:length(ahalines)
                    dist2line(k) = gui.LAline.xstart-ahalines{k}(1);
                  end
                  [val,mindist2lineindex] = min(abs(dist2line));
                  switch mindist2lineindex
                    case 1 %apical
                      if gui.SAslicesrest{1}(1)>gui.LAline.xend
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(1);
                      elseif gui.endslicesrest<gui.LAline.xend
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(1):gui.endslicesrest;
                      else
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(1):gui.LAline.xend;
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(1),'ydata',(ylim(2)-gui.SAslicesrest{1}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(1),'xdata',gui.SAslicesrest{1}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{1})));
                      set(gui.handles.axesSAapicalimageright,'cdata',gui.SAapicalimagerest(:,:,1));
                    case 2 %apical and mid
                      if gui.SAslicesrest{2}(1)>gui.LAline.xend
                        gui.SAslicesrest{1} = gui.SAslicesrest{2}(1)+1:gui.SAslicesrest{1}(end);
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(1);
                      elseif gui.SAslicesrest{1}(end)<gui.LAline.xend
                        gui.SAslicesrest{1} = gui.SAslicesrest{1}(end);
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(1):gui.SAslicesrest{1}(end)+1;
                      else
                        gui.SAslicesrest{1} = gui.LAline.xend-1:gui.SAslicesrest{1}(end);
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(1):gui.LAline.xend;
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(2),'ydata',(ylim(2)-gui.SAslicesrest{2}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(2),'xdata',gui.SAslicesrest{2}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{1})));
                      set(gui.handles.axesSAapicalimageright,'cdata',gui.SAapicalimagerest(:,:,1));
                      gui.SAmidimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{2})));
                      set(gui.handles.axesSAmidimageright,'cdata',gui.SAmidimagerest(:,:,1));
                    case 3 %mid and basal
                      if gui.SAslicesrest{3}(1)>gui.LAline.xend
                        gui.SAslicesrest{2} = gui.SAslicesrest{3}(1)+1:gui.SAslicesrest{2}(end);
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(1);
                      elseif gui.SAslicesrest{2}(end)<gui.LAline.xend
                        gui.SAslicesrest{2} = gui.SAslicesrest{2}(end);
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(1):gui.SAslicesrest{2}(end)+1;
                      else
                        gui.SAslicesrest{2} = gui.LAline.xend-1:gui.SAslicesrest{2}(end);
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(1):gui.LAline.xend;
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(3),'ydata',(ylim(2)-gui.SAslicesrest{3}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(3),'xdata',gui.SAslicesrest{3}(end)*[1 1]);
                      %SA
                      gui.SAmidimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{2})));
                      set(gui.handles.axesSAmidimageright,'cdata',gui.SAmidimagerest(:,:,1));
                      gui.SAbasalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{3})));
                      set(gui.handles.axesSAbasalimageright,'cdata',gui.SAbasalimagerest(:,:,1));
                    case 4 %basal
                      if gui.startslicesrest>gui.LAline.xend
                        gui.SAslicesrest{3} = gui.startslicesrest:gui.SAslicesrest{3}(end);
                      elseif gui.SAslicesrest{3}(end)<gui.LAline.xend
                        gui.SAslicesrest{3} = gui.SAslicesrest{3}(end);
                      else
                        gui.SAslicesrest{3} = gui.LAline.xend-1:gui.SAslicesrest{3}(end);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineright(4),'ydata',(ylim(2)-gui.SAslicesrest{3}(1))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineright(4),'xdata',gui.SAslicesrest{3}(1)*[1 1]);
                      %SA
                      gui.SAbasalimagerest = gui.restIM(:,:,1,round(mean(gui.SAslicesrest{3})));
                      set(gui.handles.axesSAbasalimageright,'cdata',gui.SAbasalimagerest(:,:,1));
                  end
                case 'HLAleft'
                  ahalines = get(gui.handles.axesHLAahalineleft,'ydata');
                  ylim = get(gui.handles.axesHLAleft,'ylim');
                  for k = 1:length(ahalines)
                    dist2line(k) = gui.LAline.ystart-ahalines{k}(1);
                  end
                  [~,mindist2lineindex] = min(abs(dist2line));
                  switch mindist2lineindex
                    case 1 %apical
                      if gui.SAslicesstress{1}(1)>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(1);
                      elseif gui.endslicesstress<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(1):gui.endslicesstress;
                      else
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(1):(ylim(2)-gui.LAline.yend);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(1),'ydata',(ylim(2)-gui.SAslicesstress{1}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(1),'xdata',gui.SAslicesstress{1}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{1})));
                      set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagestress(:,:,1));
                    case 2 %apical and mid
                      if gui.SAslicesstress{2}(1)>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{1} = gui.SAslicesstress{2}(1)+1:gui.SAslicesstress{1}(end);
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(1);
                      elseif gui.SAslicesstress{1}(end)<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(end);
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(1):gui.SAslicesstress{1}(end)+1;
                      else
                        gui.SAslicesstress{1} = (ylim(2)-gui.LAline.yend)-1:gui.SAslicesstress{1}(end);
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(1):(ylim(2)-gui.LAline.yend);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(2),'ydata',(ylim(2)-gui.SAslicesstress{2}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(2),'xdata',gui.SAslicesstress{2}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{1})));
                      set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagestress(:,:,1));
                      gui.SAmidimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{2})));
                      set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagestress(:,:,1));
                    case 3 %mid and basal
                      if gui.SAslicesstress{3}(1)>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{2} = gui.SAslicesstress{3}(1)+1:gui.SAslicesstress{2}(end);
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(1);
                      elseif gui.SAslicesstress{2}(end)<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(end);
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(1):gui.SAslicesstress{2}(end)+1;
                      else
                        gui.SAslicesstress{2} = (ylim(2)-gui.LAline.yend)-1:gui.SAslicesstress{2}(end);
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(1):(ylim(2)-gui.LAline.yend);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(3),'ydata',(ylim(2)-gui.SAslicesstress{3}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(3),'xdata',gui.SAslicesstress{3}(end)*[1 1]);
                      %SA
                      gui.SAmidimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{2})));
                      set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagestress(:,:,1));
                      gui.SAbasalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{3})));
                      set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagestress(:,:,1));
                    case 4 %basal
                      if gui.startslicesstress>(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{3} = gui.startslicesstress:gui.SAslicesstress{3}(end);
                      elseif gui.SAslicesstress{3}(end)<(ylim(2)-gui.LAline.yend)
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(end);
                      else
                        gui.SAslicesstress{3} = (ylim(2)-gui.LAline.yend)-1:gui.SAslicesstress{3}(end);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(4),'ydata',(ylim(2)-gui.SAslicesstress{3}(1))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(4),'xdata',gui.SAslicesstress{3}(1)*[1 1]);
                      %SA
                      gui.SAbasalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{3})));
                      set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagestress(:,:,1));
                  end
                case 'VLAleft'
                  ahalines = get(gui.handles.axesVLAahalineleft,'xdata');
                  ylim = get(gui.handles.axesHLAleft,'ylim');
                  for k = 1:length(ahalines)
                    dist2line(k) = gui.LAline.xstart-ahalines{k}(1);
                  end
                  [val,mindist2lineindex] = min(abs(dist2line));
                  switch mindist2lineindex
                    case 1 %apical
                      if gui.SAslicesstress{1}(1)>gui.LAline.xend
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(1);
                      elseif gui.endslicesstress<gui.LAline.xend
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(1):gui.endslicesstress;
                      else
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(1):gui.LAline.xend;
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(1),'ydata',(ylim(2)-gui.SAslicesstress{1}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(1),'xdata',gui.SAslicesstress{1}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{1})));
                      set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagestress(:,:,1));
                    case 2 %apical and mid
                      if gui.SAslicesstress{2}(1)>gui.LAline.xend
                        gui.SAslicesstress{1} = gui.SAslicesstress{2}(1)+1:gui.SAslicesstress{1}(end);
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(1);
                      elseif gui.SAslicesstress{1}(end)<gui.LAline.xend
                        gui.SAslicesstress{1} = gui.SAslicesstress{1}(end);
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(1):gui.SAslicesstress{1}(end)+1;
                      else
                        gui.SAslicesstress{1} = gui.LAline.xend-1:gui.SAslicesstress{1}(end);
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(1):gui.LAline.xend;
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(2),'ydata',(ylim(2)-gui.SAslicesstress{2}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(2),'xdata',gui.SAslicesstress{2}(end)*[1 1]);
                      %SA
                      gui.SAapicalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{1})));
                      set(gui.handles.axesSAapicalimageleft,'cdata',gui.SAapicalimagestress(:,:,1));
                      gui.SAmidimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{2})));
                      set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagestress(:,:,1));
                    case 3 %mid and basal
                      if gui.SAslicesstress{3}(1)>gui.LAline.xend
                        gui.SAslicesstress{2} = gui.SAslicesstress{3}(1)+1:gui.SAslicesstress{2}(end);
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(1);
                      elseif gui.SAslicesstress{2}(end)<gui.LAline.xend
                        gui.SAslicesstress{2} = gui.SAslicesstress{2}(end);
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(1):gui.SAslicesstress{2}(end)+1;
                      else
                        gui.SAslicesstress{2} = gui.LAline.xend-1:gui.SAslicesstress{2}(end);
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(1):gui.LAline.xend;
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(3),'ydata',(ylim(2)-gui.SAslicesstress{3}(end))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(3),'xdata',gui.SAslicesstress{3}(end)*[1 1]);
                      %SA
                      gui.SAmidimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{2})));
                      set(gui.handles.axesSAmidimageleft,'cdata',gui.SAmidimagestress(:,:,1));
                      gui.SAbasalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{3})));
                      set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagestress(:,:,1));
                    case 4 %basal
                      if gui.startslicesstress>gui.LAline.xend
                        gui.SAslicesstress{3} = gui.startslicesstress:gui.SAslicesstress{3}(end);
                      elseif gui.SAslicesstress{3}(end)<gui.LAline.xend
                        gui.SAslicesstress{3} = gui.SAslicesstress{3}(end);
                      else
                        gui.SAslicesstress{3} = gui.LAline.xend-1:gui.SAslicesstress{3}(end);
                      end
                      %HLA
                      set(gui.handles.axesHLAahalineleft(4),'ydata',(ylim(2)-gui.SAslicesstress{3}(1))*[1 1]);
                      %VLA
                      set(gui.handles.axesVLAahalineleft(4),'xdata',gui.SAslicesstress{3}(1)*[1 1]);
                      %SA
                      gui.SAbasalimagestress = gui.stressIM(:,:,1,round(mean(gui.SAslicesstress{3})));
                      set(gui.handles.axesSAbasalimageleft,'cdata',gui.SAbasalimagestress(:,:,1));
                  end
              end
              updateLVsegmentation; % plotLVsegmentation; % 
            end
        end
    end
  end
end
% 
% %------------------------------------------------
% function image = upsampleimage(image,no,SAslices)
% %------------------------------------------------
% %upsample the ungated SA image
% 
% global SET
% 
% %set upsampling factor
% resampledresolution = 1;
% %smoothing parameters
% sigma = [2 2 2];
% if mean([SET(no).ResolutionX SET(no).ResolutionY]) < resampledresolution
%   resamplexy = 1;
%   fs = mean([SET(no).ResolutionX SET(no).ResolutionY]);
% else
%   resamplexy = resampledresolution;
%   fs = mean([SET(no).ResolutionX SET(no).ResolutionY])/resampledresolution;
% end
% if SET(no).SliceThickness+SET(no).SliceGap < resampledresolution
%   resamplez = 1;
%   fsz = SET(no).SliceThickness+SET(no).SliceGap;
% else
%   resamplez = resampledresolution;
%   fsz = (SET(no).SliceThickness+SET(no).SliceGap)/resampledresolution;
% end
% %upsample images spatially
% upsampledimage = upsampleimagespatial(fs,fs,image);
% %smooth images
% upsampledimage = smoothimage(upsampledimage,sigma(1),sigma(1),sigma(1));
% upsampledimage = squeeze(mean(upsampledimage(:,:,:,SAslices),4));
                        
                        
%-------------------------------
function keypressed(fignum,evnt) %#ok<INUSL>
%-------------------------------
%move in the image slices with the keyboard arrows

global DATA

gui = DATA.GUI.SpectPlot2d;

if nargin==0
  mydisp('No inputarguments.');
  return;
end;

key = getkey(evnt);
switch key
  case 'leftarrow' %Arrow left
    prev_Callback;
  case 'rightarrow' %Arrow right
    next_Callback;
  case 'downarrow'
    framerate = get(gui.handles.sliderframerate,'value');
    framerate = max(0,min(1,framerate-0.1));
    set(gui.handles.sliderframerate,'value',framerate);
    setframerate_Callback;
  case 'uparrow'
    framerate = get(gui.handles.sliderframerate,'value');
    framerate = max(0,min(1,framerate+0.1));
    set(gui.handles.sliderframerate,'value',framerate);
    setframerate_Callback;
  case 'shift-ctrl-n'
    mode = gui.mode;
    succeed = close_Callback;
    if succeed
      succeed = filemenu('loadnext_Callback'); %close the file and open the next file
      if succeed
        spect.spectplot2d('init_Callback',mode); %open the scoring gui
      end
    end
end


%----------------------------
function editscoring_Callback %#ok<DEFNU>
%----------------------------
%update the scoring values

global DATA

gui = DATA.GUI.SpectPlot2d;

currentobj = get(gui.fig,'CurrentObject');
currenttag = get(currentobj,'Tag');
SRS = 0;
SSS = 0;
SDS = 0;
pos = [17 13 16:-1:14 7 12:-1:8 1 6:-1:2];

if isequal(currentobj.Style,'edit')
if get(gui.handles.radiobuttongated,'value')
  tags = {'editHLAright', ...
    'editapicalantright','editapicallatright', ...
    'editapicalinfright','editapicalsepright', ...
    'editmidantright','editmidantlatright','editmidinflatright', ...
    'editmidinfright','editmidinfsepright','editmidantsepright', ...
    'editbasalantright','editbasalantlatright','editbasalinflatright', ...
    'editbasalinfright','editbasalinfsepright','editbasalantsepright'};
  currentpos = find(strcmp(currenttag,tags));
  gui.scoringvalues{3,pos(currentpos)} = str2num(get(currentobj,'String'));
  for loop = 1:size(gui.scoringvalues,2)
    SRS = SRS+gui.scoringvalues{3,loop};
  end
  set(gui.handles.textScoring,'String',['SRS: ',num2str(SRS)]);
else
  tags = {'editHLAleft', ...
    'editapicalantleft','editapicallatleft', ...
    'editapicalinfleft','editapicalsepleft', ...
    'editmidantleft','editmidantlatleft','editmidinflatleft', ...
    'editmidinfleft','editmidinfsepleft','editmidantsepleft', ...
    'editbasalantleft','editbasalantlatleft','editbasalinflatleft', ...
    'editbasalinfleft','editbasalinfsepleft','editbasalantsepleft', ...
    'editHLAright', ...
    'editapicalantright','editapicallatright', ...
    'editapicalinfright','editapicalsepright', ...
    'editmidantright','editmidantlatright','editmidinflatright', ...
    'editmidinfright','editmidinfsepright','editmidantsepright', ...
    'editbasalantright','editbasalantlatright','editbasalinflatright', ...
    'editbasalinfright','editbasalinfsepright','editbasalantsepright'};
  pos = [pos pos];
  currentpos = find(strcmp(currenttag,tags));
  if currentpos <= 17
    gui.scoringvalues{1,pos(currentpos)} = str2num(get(currentobj,'String'));
  else
    gui.scoringvalues{2,pos(currentpos)} = str2num(get(currentobj,'String'));
  end
  for loop = 1:size(gui.scoringvalues,2)
    SSS = SSS+gui.scoringvalues{1,loop};
    SRS = SRS+gui.scoringvalues{2,loop};
    SDS = SDS+max(0,gui.scoringvalues{1,loop}-gui.scoringvalues{2,loop});
  end
  scoringvalues = sprintf('%s\n%s\n%s',['SSS: ',num2str(SSS)],['SRS: ',num2str(SRS)],['SDS: ',num2str(SDS)]);
  set(gui.handles.textScoring,'String',scoringvalues);
end
end

%---------------------------------------------------------------------
function lvv = calclvv(startslices,endslices,endox,endoy,t,resXY,resZ)
%---------------------------------------------------------------------
%calculate the temporal lvv based on the endo and epicardium segmentation
%
%input
%startslices:     most basal slice with LV segmentation for each time frame
%endslices:       most apical slice with LV segmentation for each time frame
%endox / y:       myocardial borders
%t:               time frames
%slicethickness:  slice thickness in the resampled image stack
%xres / yres:     resolution in the resampled image stack
%output
%lvv:             LV blood volume
%
%written by Helen Soneson 2010-06-04

lvv = zeros(1,length(t));
loop = 1;
for tloop = t
  templvv = 0;
  for sliceloop = startslices(tloop):endslices(tloop)
    if ~isnan(endox(1,tloop,sliceloop))
      endoarea = polyarea(endoy(:,tloop,sliceloop).*resXY,endox(:,tloop,sliceloop).*resXY);
      templvv = templvv+endoarea;
    end
  end
  lvv(loop) = templvv*resZ/1000;
  loop = loop+1;
end
