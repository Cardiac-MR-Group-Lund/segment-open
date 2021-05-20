function autoperfusionspect
%--------------------------
%automatic segmentation of LV and perfusion analysis in SPECT images
%perform on .mat files with at least 3 image stacks
%(rest ungated, rest gated, stress ungated)

global SET NO DATA

%close open files
segment('filecloseall_Callback',1);

%selecting the folder
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder with .mat files');
if isequal(pathname,0)
%   myfailed('Aborted.');
  return;
end
%find the .mat files
files2load = dir([pathname filesep '*.mat']);
numfiles = length(files2load);
%Error check
if numfiles==0
  myfailed('Found no files to export from.');
  return;
end
marginSA = 15;  %number of pixels around the LV seg for cropping

%output to spreadsheet
output = cell(numfiles+1,7);
output{1,1} = 'Filename';
output{1,2} = 'Rest LVM [g]';
output{1,3} = 'Stress LVM [g]';
output{1,4} = 'EDV [ml]';
output{1,5} = 'ESV [ml]';
output{1,6} = 'EF [%]';
output{1,7} = 'D-TPD';
output{1,8} = 'S-TPD';
output{1,9} = 'Manual SSS';
output{1,10} = 'Corrected Manual SSS';
output{1,11} = 'Manual SRS';
output{1,12} = 'Manual SDS';
fullanalysis = yesno('Re-analyze all data?'); % ask user

for fileloop = 1:numfiles
  %Load file
  DATA.Silent = false;
  disp(sprintf('Loading %s.',files2load(fileloop).name));
  filename = [pathname filesep files2load(fileloop).name];
  SET=[];
  load(filename,'-mat');
  SET = setstruct;
  clear setstruct;
  DATA.Preview.PreviewFile = filename;
  
  openfile('setupstacksfrommat',1);
  segment('renderstacksfrommat');
  [path,name,ext]=fileparts(filename);
  
  %find image stacks
  restno = [];
  stressno = [];
  restgatedno = [];
  stressgatedno = [];
  for loop = 1:length(SET)
    if isequal(SET(loop).ImageType,'Perfusion Rest')
      if SET(loop).TSize > 1
        restgatedno = loop;
      else
        restno = loop;
      end
    elseif isequal(SET(loop).ImageType,'Perfusion Stress')
      if SET(loop).TSize > 1
        stressgatedno = loop;
      else
        stressno = loop;
      end
    end
  end
  
  doanalysis = 1;
  if ~fullanalysis
    if isfield(SET(stressno),'Perfusion') && ~isempty(SET(stressno).Perfusion)
      if isfield(SET(stressno).Perfusion,'MPS')
        if isfield(SET(stressno).Perfusion.MPS,'TPD')
          doanalysis = 0;
        end
      end
    end
  end
  
  if doanalysis
    for loop = 1:length(SET)
      SET(loop).Linked = loop;
      viewfunctions('setview',1,1,NO,{'one'})
      tools('setcolormap_Callback','spect');
      if SET(NO).ResolutionX > 4
        tools('upsampleimage_Callback',2);
      end
      tools('setcolormap_Callback','spect');
      if SET(NO).ResolutionX > 4
        tools('upsampleimage_Callback',2);
      end
      %LV segmentation
      SET(NO).StartSlice = 1;
      SET(NO).EndSlice = 1;
      SET(NO).CurrentSlice = 1;
      spect.spectlvsegmentation('segmentlv_Callback',0,NO,NO);
      %crop in SA
      startslice = find(~isnan(SET(NO).EpiX(1,1,:)),1,'first');  %the most basal slice
      endslice = find(~isnan(SET(NO).EpiX(1,1,:)),1,'last');  %the most apical slice
      midslice = round(mean([startslice endslice]));
      diameter = max(SET(NO).EpiY(:,1,midslice))-min(SET(NO).EpiY(:,1,midslice));
      ydown = max(round(min(SET(NO).EpiX(:))-marginSA),1);
      yup = min(round(max(SET(NO).EpiX(:))+marginSA),SET(NO).XSize);
      xleft = max(round(min(SET(NO).EpiY(:))-marginSA-0.5*diameter),1);
      xright = min(round(max(SET(NO).EpiY(:))+marginSA),SET(NO).YSize);
      xind = ydown:yup;
      yind = xleft:xright;
      if length(xind) < size(SET(NO).IM,1) || length(yind) < size(SET(NO).IM,2)
        nos = tools('crophelper',NO,xind,yind);
      end
    end
    %perfusion analysis
    if ~isempty(restno)
      spect.spectperfusionsegmentation('perfusionregistration_Callback',0);
      spect.spectregistration('analysis_Callback');
    else
      spect.spectperfusionsegmentation('perfusionanalyzeStress',stressno,0);
    end
  else
    %corrected error in old images
    if isfield(SET(stressno).Perfusion.MPS,'registration')
      if ~isfield(SET(stressno).Perfusion.MPS.registration,'Linked')
        SET(stressno).Perfusion.MPS.registration.Linked = SET(stressno).Linked;
      end
    else
      SET(stressno).Perfusion.MPS.analysis.resolution = [SET(stressno).ResolutionX SET(stressno).ResolutionY];
    end
    for loop = 1:length(SET)
      SET(loop).Linked = loop;
    end
  end
  
  viewfunctions('setview')

  %saving file
  DATA.Buffer.KeyStroke = {'ok'}; %Put ok in queue for confirm box.
  filemenu('saveall_Callback')
  %output to spreadsheet
  output{1+fileloop,1} = files2load(fileloop).name;
  if ~isempty(restno)
    output{1+fileloop,2} = SET(restno).LVM*1.05;
  end
  if ~isempty(restgatedno)
    output{1+fileloop,4} = SET(restgatedno).EDV;
    output{1+fileloop,5} = SET(restgatedno).ESV;
    output{1+fileloop,6} = SET(restgatedno).EF*100;
  elseif ~isempty(stressgatedno)
    output{1+fileloop,4} = SET(stressgatedno).EDV;
    output{1+fileloop,5} = SET(stressgatedno).ESV;
    output{1+fileloop,6} = SET(stressgatedno).EF*100;
  end
  if ~isempty(restno) && ~isempty(stressno)
    output{1+fileloop,7} = SET(stressno).Perfusion.MPS.TPDchange;
  end
  if ~isempty(stressno)
      output{1+fileloop,3} = SET(stressno).LVM*1.05;
      output{1+fileloop,8} = SET(stressno).Perfusion.MPS.TPD;
      %manual scoring
      SSSisch = 0;
      SSSscores = zeros(1,17);
      SRSscores = zeros(1,17);
      SDSscores = zeros(1,17);
      if isfield(SET(stressno).Perfusion.MPS,'ScoringIschManual')
        for k = 1:17
          stressisch = SET(stressno).Perfusion.MPS.ScoringIschManual(k);
          SSSisch = SSSisch+stressisch;
          SSSscores(k) = stressisch;
        end
      end
      if isfield(SET(restno).Perfusion.MPS,'ScoringIschManual')
        for k = 1:17
          restisch = SET(restno).Perfusion.MPS.ScoringIschManual(k);
%           SRSisch = SRSisch+restisch;
          SRSscores(k) = restisch;  
        end
      end
      if isfield(SET(stressno).Perfusion.MPS,'ScoringIschManual') && isfield(SET(restno).Perfusion.MPS,'ScoringIschManual') 
        for k = 1:17
          diffisch = max(0,SET(stressno).Perfusion.MPS.ScoringIschManual(k)-SET(restno).Perfusion.MPS.ScoringIschManual(k));
%           SDSisch = SDSisch+diffisch;
          SDSscores(k) = diffisch;
        end 
      end
      output{1+fileloop,9} = sum(SSSscores);
      output{1+fileloop,11} = sum(SRSscores);
      output{1+fileloop,12} = sum(SDSscores);
      %correct for single 1-scores
      segments{1} = [2 6 7];
      segments{2} = [1 3 8];
      segments{3} = [2 4 9];
      segments{4} = [3 5 10];
      segments{5} = [4 6 11];
      segments{6} = [1 5 12];
      segments{7} = [1 8 12 13];
      segments{8} = [2 7 9 14];
      segments{9} = [3 8 10 14];
      segments{10} = [4 9 11 15];
      segments{11} = [5 10 12 16];
      segments{12} = [6 7 11 16];
      segments{13} = [7 14 16 17];
      segments{14} = [8 9 13 15 17];
      segments{15} = [10 14 16 17];
      segments{16} = [11 12 13 15 17];
      segments{17} = [13 14 15 16];
      for k = 1:17
          seg = segments{k};
          if sum(SSSscores(seg)) < 1 && SSSscores(k) == 1
              SSSscores(k) = 0;
          end
      end    
      output{1+fileloop,10} = sum(SSSscores);
  end
%   filenamesave = ['pt' sprintf('%03.0f',fileloop+284) 'te'];
%   DATA.Buffer.KeyStroke = {'ok'}; %Put ok in queue for confirm box.
%   filemenu('saveallas_helper',pathname,filenamesave);
  segment('filecloseall_Callback',1);
  SET = [];
end
segment('cell2clipboard',output);
  